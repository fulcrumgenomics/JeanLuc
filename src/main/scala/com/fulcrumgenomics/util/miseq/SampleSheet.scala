/**
  * Copyright (c) 2016, Fulcrum Genomics LLC
  * All rights reserved.
  *
  * Redistribution and use in source and binary forms, with or without
  * modification, are permitted provided that the following conditions are met:
  *
  * 1. Redistributions of source code must retain the above copyright notice,
  * this list of conditions and the following disclaimer.
  *
  * 2. Redistributions in binary form must reproduce the above copyright notice,
  * this list of conditions and the following disclaimer in the documentation
  * and/or other materials provided with the distribution.
  *
  * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
  * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  * POSSIBILITY OF SUCH DAMAGE.
  */
package com.fulcrumgenomics.util.miseq

import java.nio.file.Path

import scala.collection.mutable.ListBuffer
import scala.io.Source

/**
  * Stores information about samples from an Illumina Experiment Manager Sample Sheet (typically a MiSeq).  The
  * samples may also include derived information.  If library id is not specified, it will be set to the sample name.
  *
  * Optional fields include support for specifying the expected sample barcode for each sample.  The sample barcode
  * can be present as a sub-sequence (or sub-sequences) in the i7 or i5 read.  If additional bases are found in the i7
  * or i5 read, such as molecular barcodes, the should be included as Ns.  It is up to the developer to obtain the
  * correct read structure elsewhere to infer which bases are sample barcode and which bases are not (ex. molecular
  * identifiers).  If the sample barcode is inline in either read of a pair, the sample barcode sequence can be specified
  * here.  In this case, only include the sample barcode, and not the additional bases, such as template bases.
  *
  * NB: the lookup in the columns in the Sample Sheet is case insensitive.
  *
  * @author Nils Homer
  */
object SampleSheet {
  /** Required header names for standard fields in the Illumina Experiment Manager Sample Sheet. */
  val SampleId: String       = "Sample_ID"
  val SampleName: String     = "Sample_Name"
  val LibraryId: String      = "Library_ID"
  val SampleProject: String  = "Sample_Project"
  val Description: String    = "Description"
  /** Optional header names for standard fields in the Illumina Experiment Manager Sample Sheet. */
  val R1BarcodeBases: String = "R1_Barcode_Bases"
  val R2BarcodeBases: String = "R2_Barcode_Bases"
  val I7Bases: String        = "Index"
  val I5Bases: String        = "Index2"

  /** Attempts to clean special characters from a sample id or name as Illumina's software does. */
  def cleanMiseqSampleId(id: String): String = {
    id.map {
      case '#' => ' '
      case '_' | '+' | ' ' => '-'
      case c => c
    }
  }

  /** Gets the trimmed value of the given key (`name`) from the map.  If empty or not found, returns None */
  private def getStringField(sampleDatum: Map[String, String], name: String): Option[String] = {
    sampleDatum.get(name.toUpperCase) match {
      case None => None
      case Some(str) if str.isEmpty => None
      case Some(str) => Some(str.trim)
    }
  }
}

/** Create a sample sheet from the given file (typically 'SampleSheet.csv'). If the library id is omitted, the
  * sample name will be used in its stead. */
class SampleSheet(sampleSheet: Path) extends Iterable[Sample] {
  private val samples: ListBuffer[Sample] = new ListBuffer[Sample]()

  getSampleData(sampleSheet).zipWithIndex.foreach {
    case (sampleDatum, sampleOrdinal) => this.samples.append(makeSample(sampleOrdinal+1, sampleDatum))
  }

  def iterator: Iterator[Sample] = this.samples.iterator

  override def size: Int = this.samples.size

  def get(index: Int): Sample = this.samples(index)

  /** Reads in the the Sample data only, one per row */
  private def getSampleData(sampleSheet: Path): List[Map[String, String]] = {
    // ignore data pre-"[Data]"
    val (preData, postData) = Source.fromFile(sampleSheet.toFile).getLines().span(!_.startsWith("[Data]"))

    // get the header (keys) (skip "[Data]")
    val header = (postData.drop(1).toList.headOption match {
      case Some(line) => line.split(",", -1) // NB: include trailing empty strings
      case None => throw new IllegalArgumentException("Could not find the header for sample data.")
    }).map(_.toUpperCase.trim)

    // get the rows (values), so skip "[Data]" and the header
    val lineNumber = preData.size + 2 // 0-based
    postData.drop(2).zipWithIndex.map { case (line, rowNumber) =>
      val values = line.split(",", -1)
      // check we have the correct # of columns
      if (values.size != header.length) {
        throw new IllegalArgumentException(s"# of columns in the header and current row do not match ('${header.length}' != '${values.size}'): row #${rowNumber + 1} line number #${lineNumber + rowNumber + 1}")
      }
      header.zip(values.map(_.trim)).toMap
    }.toList
  }

  /** Creates a sample from the given row data */
  private def makeSample(sampleOrdinal: Int, sampleDatum: Map[String, String]): Sample = {
    val sampleName: String          = SampleSheet.getStringField(sampleDatum, SampleSheet.SampleName)    getOrElse (throw new IllegalArgumentException(s"Missing: ${SampleSheet.SampleName}"))
    val sampleId: String            = SampleSheet.getStringField(sampleDatum, SampleSheet.SampleId)      getOrElse (throw new IllegalArgumentException(s"Missing: ${SampleSheet.SampleId}"))
    val libraryId: Option[String]   = SampleSheet.getStringField(sampleDatum, SampleSheet.LibraryId) orElse Some(sampleName)
    val project: Option[String]     = SampleSheet.getStringField(sampleDatum, SampleSheet.SampleProject)
    val description: Option[String] = SampleSheet.getStringField(sampleDatum, SampleSheet.Description)
    new Sample(
      sampleOrdinal        = sampleOrdinal,
      sampleId             = sampleId,
      sampleName           = sampleName,
      libraryId            = libraryId,
      project              = project,
      description          = description,
      r1SampleBarcodeBases = SampleSheet.getStringField(sampleDatum, SampleSheet.R1BarcodeBases).map(_.getBytes),
      r2SampleBarcodeBases = SampleSheet.getStringField(sampleDatum, SampleSheet.R2BarcodeBases).map(_.getBytes),
      i7IndexBases         = SampleSheet.getStringField(sampleDatum, SampleSheet.I7Bases).map(_.getBytes),
      i5IndexBases         = SampleSheet.getStringField(sampleDatum, SampleSheet.I5Bases).map(_.getBytes)
    )
  }
}
