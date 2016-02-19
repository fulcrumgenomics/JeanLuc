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

package com.fulcrumgenomics.personal.nhomer

import java.nio.file.Path

import com.fulcrumgenomics.cmdline.{JeanLucTool, Personal}
import dagr.commons.CommonsDef.PathToBam
import dagr.commons.io.Io
import dagr.sopt._
import htsjdk.samtools._
import htsjdk.samtools.metrics.{MetricBase, MetricsFile}
import htsjdk.samtools.util._
import picard.PicardException

import scala.io.Source

@clp(
  description = "Searches for DNA sequences in the read pairs.",
  group = classOf[Personal]
)
class HasSequence
( @arg(doc = "Input SAM or BAM.") val input: PathToBam,
  @arg(doc = "Output SAM or BAM.") val output: PathToBam,
  @arg(doc = "File containing the DNA sequences to search for.") val sequences: Path,
  @arg(doc = "Sequence match histogram written to this file.") val metrics: Path,
  @arg(doc = "Maximum mismatches for matching a sequence.") val maxMismatches: Int = 1,
  @arg(doc = "Minimum base quality. Any bases falling below this quality will be considered a mismatch even in the bases match.", flag = "q") val minimumBaseQuality: Int = 0,
  @arg(doc = "Treat no calls in the input sequences as DNA bases.") val noCallsAreBases: Boolean = false,
  @arg(doc = "Search for mismatches at the start only") val startOnly: Boolean = true,
  @arg(doc = "The tag to store the result.") val hasSequenceTag: String = "XW"
) extends JeanLucTool {
  Io.assertReadable(input)
  Io.assertCanWriteFile(output)
  Io.assertCanWriteFile(metrics)

  // TODO: progress logging and metrics file header
  override def execute: Int = {
    val sequencesReadOne: Iterator[String] = for(line <- Source.fromFile(this.sequences.toFile).getLines()) yield line.trim
    val sequencesReadTwo: Iterator[String] = for (sequence <- sequencesReadOne) yield SequenceUtil.reverseComplement(sequence)
    val reader: SamReader = SamReaderFactory.makeDefault.open(input.toFile)
    val header: SAMFileHeader = reader.getFileHeader
    if (header.getSortOrder ne SAMFileHeader.SortOrder.queryname) {
      throw new PicardException("Expects a queryname sorted input file, was: " + header.getSortOrder.name)
    }
    val writer: SAMFileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, output.toFile)
    val histogram: Histogram[String] = new Histogram[String]
    (sequencesReadOne ++ sequencesReadTwo).foreach { seq => histogram.increment(seq, 0) }
    val iterator = reader.iterator
    while (iterator.hasNext) {
      val rec: SAMRecord = iterator.next
      if (iterator.hasNext && rec.getReadPairedFlag) {
        val rec2 = iterator.next
        annotateRecords(rec, rec2, sequencesReadOne, sequencesReadTwo, histogram, maxMismatches, minimumBaseQuality, noCallsAreBases, startOnly)
        writer.addAlignment(rec)
        writer.addAlignment(rec2)
      } else writer.addAlignment(rec)
    }
    CloserUtil.close(reader)
    writer.close()
    val outputHistogram: Histogram[String] = new Histogram[String]("sequence", "count")
    sequencesReadOne.zip(sequencesReadTwo).foreach {
      case (seqOne, seqTwo) =>
        outputHistogram.increment(seqOne, histogram.get(seqOne).getValue + histogram.get(seqTwo).getValue)
    }
    val metricsFile: MetricsFile[MetricBase, String] = new MetricsFile[MetricBase, String]()
    metricsFile.addHistogram(outputHistogram)
    metricsFile.write(metrics.toFile)
    0
  }

  private def annotateRecords(rec: SAMRecord,
                              rec2: SAMRecord,
                              readOneSequences: Iterator[String],
                              readTwoSequences: Iterator[String],
                              histogram: Histogram[String],
                              maxMismatches: Int,
                              minimumBaseQuality: Int,
                              noCallsAreBases: Boolean,
                              startOnly: Boolean) {
    val matchesReadOne: Boolean = matchesSequence(rec, readOneSequences, histogram, maxMismatches, minimumBaseQuality, noCallsAreBases, startOnly)
    val matchesReadTwo: Boolean = matchesSequence(rec2, readTwoSequences, histogram, maxMismatches, minimumBaseQuality, noCallsAreBases, startOnly)
    val value: Int = (matchesReadOne, matchesReadTwo) match {
      case (true, true)   => 3
      case (false, true)  => 2
      case (true, false)  => 1
      case (false, false) => 0
    }
    rec.setAttribute(hasSequenceTag, value)
    rec2.setAttribute(hasSequenceTag, value)
  }

  private def matchesSequence(rec: SAMRecord,
                      sequences: Iterator[String],
                      histogram: Histogram[String],
                      maxMismatches: Int,
                      minimumBaseQuality: Int,
                      noCallsAreBases: Boolean,
                      startOnly: Boolean): Boolean = {
    val countMM: (Int, Array[Byte]) => Int = countMismatches(rec.getReadBases, rec.getBaseQualities, minimumBaseQuality, noCallsAreBases)
    val sequence = if (startOnly) {
      sequences.find { sequence => countMM(0, sequence.getBytes) <= maxMismatches }
    }
    else {
      sequences.find {
        sequence => (for (i <- (0 until (rec.getReadLength - sequence.length)).view) yield countMM(i, sequence.getBytes)).exists(_ <= maxMismatches)
      }
    }
    sequence match {
      case Some(seq) => histogram.increment(seq); true
      case None => false
    }
  }

  private def countMismatches(observedBases: Array[Byte],
                      observedQualities: Array[Byte],
                      minimumBaseQuality: Int,
                      noCallsAreBases: Boolean)
                      (observedOffset: Int,
                       expectedBases: Array[Byte]): Int = {
    val minLength = Math.min(observedBases.length, expectedBases.length)
    val mismatches = for (i <- 0 until minLength) yield {
      if (noCallsAreBases || !SequenceUtil.isNoCall(expectedBases(i))) {
        if (!SequenceUtil.basesEqual(observedBases(i), expectedBases(i))) 1
        else if (observedQualities(i) < minimumBaseQuality) 1
        else 0
      }
      else 0
    }
    mismatches.sum
  }
}
