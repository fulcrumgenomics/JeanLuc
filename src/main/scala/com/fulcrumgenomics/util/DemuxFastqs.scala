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
package com.fulcrumgenomics.util

import java.io.File
import java.nio.file.Path

import com.fulcrumgenomics.cmdline.{ClpGroups, JeanLucTool}
import com.fulcrumgenomics.util.miseq.{Sample, SampleBarcode, SampleSheet}
import dagr.commons.CommonsDef.{DirPath, FilePath, PathToFastq, _}
import dagr.commons.io.Io
import dagr.commons.util.LazyLogging
import dagr.sopt.{arg, clp}
import htsjdk.samtools.SAMFileHeader.SortOrder
import htsjdk.samtools._
import htsjdk.samtools.fastq.{FastqReader, FastqRecord}
import htsjdk.samtools.metrics.MetricsFile
import htsjdk.samtools.util._
import picard.illumina.ExtractIlluminaBarcodes
import picard.illumina.parser.ReadStructure
import picard.util.IlluminaUtil

import scala.collection.JavaConversions._
import scala.collection.mutable

object DemuxFastqs {
  val UnmatchedSampleId: String = "unmatched"
  private val solexaQualityConverter: SolexaQualityConverter = SolexaQualityConverter.getSingleton

  // TODO: duplicated but private in picard.sam.FastqToSam so make it public and static there, and use it here
  /** Creates the read name to put into a SAM record from a FASTQ record. */
  private def getReadName(fastqHeader: String, paired: Boolean, stripUnpairedMateNumber: Boolean): String = {
    val idx: Int = fastqHeader.indexOf(" ")
    val readName: String = if (idx == -1) fastqHeader else fastqHeader.substring(0, idx)
    if (stripUnpairedMateNumber && !paired && (readName.endsWith("/1") || readName.endsWith("/2"))) {
      readName.substring(0, readName.length - 2)
    }
    else {
      readName
    }
  }

  /** Counts the nucleotide mismatches between two strings of the same length.  Ignores no calls in expectedBases.
    * Observed base qualities less than the minimum base quality are counted as mismatches if not a no call. Observed
    * qualities may not be null. */
  private[util] def countMismatches(observedBases: Array[Byte], observedQualities: Array[Byte], expectedBases: Array[Byte], minimumBaseQuality: Int): Int = {
    observedBases.zip(observedQualities).zip(expectedBases).count {
      case ((observedBase, observedQuality), expectedBase) =>
         (!SequenceUtil.isNoCall(expectedBase)) &&
           (!SequenceUtil.basesEqual(observedBase, expectedBase) || (observedQuality < minimumBaseQuality))
    }
  }

  /** Gets all the bases (concatenated) for the given cycle indexes (1-based), ignoring cycles not after the end
    * of `bases`.  The cycles should be in ascending order. */
  private def getBasesAtCycles(bases: Array[Byte], cycles: Array[Int]): Array[Byte] = {
    //  ignore cycles past the end of `bases`.
    cycles.filter(cycle => cycle - 1 < bases.length).map(cycle => bases(cycle-1))
  }

  /** Sets the sample barcode in the given sample sheet using the given read structures and bases from the sample sheet.
    *
    * The final sample barcode is the concatenation of the sample barcode from the i7 and i5 bases found in the sample
    * sheet in that order.  The concatenation uses `IlluminaUtil.barcodeSeqsToString`.  The specific bases extracted
    * from the i7 and i5 reads are the sample barcode bases specified in the given i7 and i5 read structures.  The
    * final sample barcode is then updated in the given sample sheet.  An exception is thrown if no sample barcode
    * bases are found in the i5 or i7 read structures.
    * */
  private def setSampleBarcode(sampleSheet: SampleSheet, r1ReadStructure: ReadStructure, r2ReadStructure: ReadStructure, i7ReadStructure: ReadStructure, i5ReadStructure: ReadStructure) {
    sampleSheet.foreach { sample =>
      val sampleBarcodes = List(
        (sample.r1SampleBarcodeBases, r1ReadStructure),
        (sample.r2SampleBarcodeBases, r2ReadStructure),
        (sample.i7IndexBases,         i7ReadStructure),
        (sample.i5IndexBases,         i5ReadStructure)
      ).filter { case (bases, readStructure) => readStructure.sampleBarcodes.getTotalCycles > 0 }
      .map { case (bases, readStructure) =>
        new String(getBasesAtCycles(bases.getOrElse(unreachable("Found no bases")), readStructure.sampleBarcodes.getCycles))
      }
      if (sampleBarcodes.isEmpty) throw new IllegalArgumentException("No sample barcodes found in the i7 or i5 read structures.")
      sample.sampleBarcode = Some(new SampleBarcode(sampleBarcodes))
    }
  }

  /** Reads multiple FASTQs at the same time, ensuring their read names match along the way.  This is
    * useful for when we don't have a single interleaved FASTQ. */
  private[util] class ParallelFastqReader(fastq: List[Path], fastqs: List[Path]*) extends Iterator[List[FastqRecord]] {

    private var stripUnpairedMateNumber: Boolean = false

    val readers: List[SequentialFastqReader] = (fastq :: fastqs.toList).map { otherFastqs =>
      if (otherFastqs.size != fastq.size) throw new IllegalArgumentException("List of list fastqs had differing lengths.")
      new SequentialFastqReader(otherFastqs)
    }

    def withStripUnpairedMateNumber(value: Boolean): Unit = stripUnpairedMateNumber = value

    def next: List[FastqRecord] = {
      val records = readers.map(_.next)
      verifySameReadNames(records)
      records
    }

    def hasNext: Boolean = {
      if (readersHaveNext(true)) true
      else if (readersHaveNext(false)) false
      else throw new IllegalArgumentException("Fastqs have differing lengths")
    }

    private def readersHaveNext(desiredValue: Boolean): Boolean = !readers.exists(_.hasNext != desiredValue)

    private def verifySameReadNames(records: List[FastqRecord]): Unit = {
      val readName: String = getReadName(records.head.getReadHeader, paired=true, stripUnpairedMateNumber=stripUnpairedMateNumber)
      records.drop(0).foreach { record =>
        val otherReadName: String = getReadName(record.getReadHeader, paired=true, stripUnpairedMateNumber=stripUnpairedMateNumber)
        if (readName.compareTo(otherReadName) != 0) throw new IllegalArgumentException(String.format("Mismatching read names in FASTQS:\n%s\n%s\n", readName, otherReadName))
      }
    }

    def close() = readers.foreach(_.close())
  }

  private class ReadStructureInfo(readStructureString: String) extends ReadStructure(readStructureString) {
    val templateCycles: Array[Int] = this.templates.getCycles
    val sampleBarcodeCycles: Array[Int] = this.sampleBarcodes.getCycles
    val molecularBarcodeCycles: Array[Int] = this.molecularBarcode.getCycles
  }

  /** Helper class to convert fastq records into SAM records and write them out.
    *
    * @param sampleSheet the sample sheet describing all the samples
    * @param readStructureInfos the read structures info, one for each expected fastq record that will be passed
    *                           into `convertFastqsAndWriteRecords`.
    */
  private class FastqConverter(val sampleSheet: SampleSheet,
                               val umiTag: String,
                               val qualityFormat: FastqQualityFormat,
                               val stripUnpairedMateNumber: Boolean,
                               val minQ: Int,
                               val maxQ: Int,
                               val maxMismatches: Int,
                               val minMismatchDelta: Int,
                               val maxNoCalls: Int,
                               val minBaseQuality: Int,
                               val readStructureInfos: ReadStructureInfo*) {
    /** Sample barcode to its associate metric */
    private val barcodeToMetrics: mutable.Map[String, ExtractIlluminaBarcodes.BarcodeMetric] = new mutable.LinkedHashMap[String, ExtractIlluminaBarcodes.BarcodeMetric]()

    sampleSheet.foreach { sample =>
      val barcode: String = sample.sampleBarcode match {
        case Some(b) => b.concatenatedBarcode
        case None => throw new IllegalArgumentException(s"Sample with id '${sample.sampleId}' did not have a sampleBarcode")
      }
      val libraryId: String = sample.libraryId match {
        case Some(id) => id
        case None => throw new IllegalArgumentException(s"Sample with id '${sample.sampleId}' did not have a library id")
      }
      val metric: ExtractIlluminaBarcodes.BarcodeMetric = new ExtractIlluminaBarcodes.BarcodeMetric(sample.sampleName, libraryId, barcode, Array[String](barcode))
      barcodeToMetrics.put(barcode, metric)
    }

    private val noMatchBarcodeMetric: ExtractIlluminaBarcodes.BarcodeMetric = {
      val noMatchBarcodes = readStructureInfos.filter {
        readStructureInfo => readStructureInfo.sampleBarcodes.getTotalCycles > 0
      }.map {
        readStructureInfo => StringUtil.repeatCharNTimes('N', readStructureInfo.sampleBarcodes.getTotalCycles)
      }.toList
      val noMatchBarcode = IlluminaUtil.barcodeSeqsToString(noMatchBarcodes)
      val metric = new ExtractIlluminaBarcodes.BarcodeMetric(UnmatchedSampleId, UnmatchedSampleId, noMatchBarcode, Array[String](noMatchBarcode))
      barcodeToMetrics.put(noMatchBarcode, noMatchBarcodeMetric)
      metric
    }

    /** Converts the given FASTQ records to SAM records and writes them out to the output file based on the
      * sample barcode. */
    def convertFastqsAndWriteRecords(fastqRecords: List[FastqRecord], writers: List[SAMFileWriter]): List[SAMRecord] = {
      val sampleBarcodeReadBases: Array[Byte] = getSampleBarcode(fastqRecords)
      val sampleBarcodeReadQualities: Array[Byte] = getSampleBarcodeQualities(fastqRecords)
      var sampleOrdinal: Int = getSampleOrdinalFromSampleBarcode(sampleBarcodeReadBases, sampleBarcodeReadQualities, sampleSheet)
      if (sampleOrdinal == -1) sampleOrdinal = writers.length
      val writer: SAMFileWriter = writers(sampleOrdinal - 1)
      val header: SAMFileHeader = writer.getFileHeader
      val sampleId: String = writers(sampleOrdinal - 1).getFileHeader.getReadGroups.get(0).getId
      val molecularBarcodeTag: String = getMolecularBarcode(fastqRecords)

      fastqRecords.zipWithIndex.flatMap {
        case (record, recordIndex) =>
          val readStructureInfo = readStructureInfos(recordIndex)
          if (readStructureInfo.templates.isEmpty) {
            None
          }
          else {
            if (1 < recordIndex) throw new IllegalArgumentException("Found more than two ends with template bases")
            val samRecord: SAMRecord = makeSAMRecord(header, record, readStructureInfo, sampleId, recordIndex == 0)
            samRecord.setAttribute("BC", new String(sampleBarcodeReadBases))
            samRecord.setAttribute(umiTag, molecularBarcodeTag)
            writer.addAlignment(samRecord)
            Some(samRecord)
          }
        }
      }

    def getBarcodeMetrics: Iterable[ExtractIlluminaBarcodes.BarcodeMetric] = {
      ExtractIlluminaBarcodes.finalizeMetrics(barcodeToMetrics, noMatchBarcodeMetric)
      this.barcodeToMetrics.values
    }

    /** Searches for a matching sample given the observed barcode.  Returns -1 if no match is found, given the
      * sample sheet and matching parameters */
    private def getSampleOrdinalFromSampleBarcode(sampleBarcodeReadBases: Array[Byte], sampleBarcodeReadQualities: Array[Byte], sampleSheet: SampleSheet): Int = {
      val numNoCalls: Int = sampleBarcodeReadBases.count(base => SequenceUtil.isNoCall(base))
      var bestSampleOrdinal: Int = -1
      var bestMismatches: Int = Integer.MAX_VALUE
      var secondBestMismatches: Int = Integer.MAX_VALUE

      if (numNoCalls <= maxNoCalls) {
        sampleSheet.foreach {
          sample =>
            val barcodeBytes = sample.sampleBarcode match {
              case Some(barcode) => barcode.barcodeBytes
              case None => throw new IllegalArgumentException(s"Sample barcode required for sample with id '${sample.sampleId}'")
            }
            val numMismatches: Int = countMismatches(sampleBarcodeReadBases, sampleBarcodeReadQualities, barcodeBytes, minBaseQuality)
            if (numMismatches < bestMismatches) {
              bestSampleOrdinal = sample.sampleOrdinal
              secondBestMismatches = bestMismatches
              bestMismatches = numMismatches
            }
            else if (numMismatches < secondBestMismatches) {
              secondBestMismatches = numMismatches
            }
        }
      }
      // Make sure we are within the parameter limits and update barcode metrics if necessary
      if (maxMismatches < bestMismatches || maxNoCalls < numNoCalls || (secondBestMismatches - bestMismatches) < minMismatchDelta) {
        bestSampleOrdinal = -1
        noMatchBarcodeMetric.READS += 1
        noMatchBarcodeMetric.PF_READS += 1
      }
      else {
        val sampleBarcode = sampleSheet.get(bestSampleOrdinal - 1).sampleBarcode.getOrElse(unreachable(s"Sample barcode required for sample"))
        val bestBarcodeMetric: ExtractIlluminaBarcodes.BarcodeMetric = this.barcodeToMetrics.getOrElse(sampleBarcode.concatenatedBarcode, unreachable("No metric for sample"))
        bestBarcodeMetric.READS += 1
        bestBarcodeMetric.PF_READS += 1
        if (bestMismatches == 0) {
          bestBarcodeMetric.PERFECT_MATCHES += 1
          bestBarcodeMetric.PF_PERFECT_MATCHES += 1
        }
        else if (bestMismatches == 1) {
          bestBarcodeMetric.ONE_MISMATCH_MATCHES += 1
          bestBarcodeMetric.PF_ONE_MISMATCH_MATCHES += 1
        }
      }
      bestSampleOrdinal
    }

    /** Creates a SAM record (always paired) */
    private def makeSAMRecord(header: SAMFileHeader, fastqRecord: FastqRecord, readStructureInfo: ReadStructureInfo, readGroupId: String, firstOfPair: Boolean): SAMRecord = {
      val readBases: Array[Byte] = getBasesAtCycles(fastqRecord.getReadString.getBytes, readStructureInfo.templateCycles)
      val baseQualities: Array[Byte] = getBasesAtCycles(fastqRecord.getBaseQualityString.getBytes, readStructureInfo.templateCycles)
      val record: SAMRecord = new SAMRecord(header)
      record.setReadName(getReadName(fastqRecord.getReadHeader, paired=true, stripUnpairedMateNumber=stripUnpairedMateNumber))
      record.setReadBases(readBases)
      record.setReadUnmappedFlag(true)
      record.setAttribute(ReservedTagConstants.READ_GROUP_ID, readGroupId)
      convertQuality(baseQualities, qualityFormat)
      baseQualities.foreach {
        qual =>
          val uQual: Int = qual & 0xff
          if (uQual < minQ || uQual > maxQ) {
            throw new IllegalStateException(s"Base quality $uQual is not in the range $minQ ... $maxQ for read ${fastqRecord.getReadHeader}")
          }
      }
      record.setBaseQualities(baseQualities)
      record.setReadPairedFlag(true)
      record.setMateUnmappedFlag(true)
      if (firstOfPair) record.setFirstOfPairFlag(true)
      else record.setSecondOfPairFlag(true)
      record
    }

    /** Based on the type of quality scores coming in, converts them to a numeric byte[] in phred scale. */
    private[util] def convertQuality(quals: Array[Byte], version: FastqQualityFormat) {
      version match {
        case FastqQualityFormat.Standard =>
          SAMUtils.fastqToPhred(quals)
        case FastqQualityFormat.Solexa =>
          solexaQualityConverter.convertSolexaQualityCharsToPhredBinary(quals)
        case FastqQualityFormat.Illumina =>
          solexaQualityConverter.convertSolexa_1_3_QualityCharsToPhredBinary(quals)
      }
    }

    /** Gets an array of byte arrays representing the sample barcode found in each record, according to their
      * individual read structures.  This will include empty byte[]s for any reads that don't contain a sample
      * barcode read.
      */
    private def getSampleBarcode(records: List[FastqRecord]): Array[Byte] = {
      records.zip(readStructureInfos).flatMap {
        case (record, readStructureInfo) =>
            val readBarcode: Array[Byte] = getBasesAtCycles(record.getReadString.getBytes, readStructureInfo.sampleBarcodeCycles)
          readBarcode
      }.toArray
    }

    private def getSampleBarcodeQualities(records: List[FastqRecord]): Array[Byte] = {
      records.zip(readStructureInfos).flatMap {
        case (record, readStructureInfo) =>
          val readBarcode: Array[Byte] = getBasesAtCycles(record.getBaseQualityString.getBytes, readStructureInfo.sampleBarcodeCycles)
          readBarcode
      }.toArray
    }

    private def getMolecularBarcode(records: List[FastqRecord]): String = {
      val molecularBarcodes = records.zip(readStructureInfos).map {
        case (record, readStructureInfo) =>
          getBasesAtCycles(record.getReadString.getBytes, readStructureInfo.molecularBarcodeCycles)
      }.toArray
      IlluminaUtil.barcodeSeqsToString(molecularBarcodes)
    }
  }
}

/**
  * This tool de-multiplexes a set of FASTQs based on the given Illumina Experiment Manager Sample Sheet for dual-indexed
  * sequencing runs.
  *
  * See the USAGE for a detailed description.
  *
  * Possible Future Improvements:
  * - adapter trimming (see [[IlluminaUtil.IlluminaAdapterPair]]
  * - more metadata stored in the SampleSheet.csv
  */
@clp(
  description = """|Demultiplexes FASTQs based on sample barcodes and annotates molecular barcode information.
                   |
                   |This tool demultiplexes a set of FASTQs based on the given Illumina Experiment Manager Sample Sheet for dual-indexed
                   |sequencing runs.
                   |
                   |Fastqs and read structures for read one, read two, i7 read, and i5 read should be given.  The read structures
                   |may contain sample barcodes bases ('B'), molecular identifier bases ('M'), and template bases ('T'), with the latter
                   |only for template reads (read one or read two).  Any molecular identifiers will be concatenated using the '-'
                   |delimiter and placed in the given SAM record tag ("RX" by default).  The order of concatenation will be read one,
                   |read two, i7 read, and i5 read, only considering reads with molecular identifiers.  Similarly, the sample
                   |barcode bases from the given read will be placed in the "BC" tag, using the same rules as molecular identifiers,.
                   |but applied to sample barcodes.
                   |
                   |The output directory will contain one BAM file per sample in the sample sheet, plus a BAM for reads that could
                   |not be assigned to a sample given the criteria.  The output files will be the concatenation of sample id, sample
                   |name, and sample barcode bases (expected not observed), delimited by ".".  A metrics file will also be output
                   |for sample barcodes.  More information about these metrics can be found here:
                   |  https://broadinstitute.github.io/picard/picard-metric-definitions.html#ExtractIlluminaBarcodes.BarcodeMetric
                   |
                   |The read group's sample id, sample name, and library id all correspond to the similarly named values in the sample sheet.
                   |Library id will be the sample id if not found, and the platform unit will be the sample name concatenated with the sample
                   |barcode bases delimited by a ".".
                   |
                   |The sample section of the sample sheet should contain information related to each sample with the following keys and values:
                   |  - Sample Identifier:  Sample_ID
                   |  - Sample Name:        Sample_Name
                   |  - Library Identifier: Library_ID
                   |  - Sample Project:     Sample_Project
                   |  - Description:        Description
                   |
                   |The following are optional values may included in the sample sheet to specify sample barcodes in any read:
                   |  - Read One Inline Sample Barcode Bases: R1_Barcode_Bases
                   |  - Read One Inline Sample Barcode Bases: R2_Barcode_Bases
                   |  - i7 Sample Barcode Bases:              Index
                   |  - i5 Sample Barcode Bases:              Index2
                   |
                   |The i7 and i5 Sample Barcode Bases may include extra bases before or after the sample barcode (if present
                   |and specified by the read structure).
                   |
                   |The read structures will be used to extract the observed sample barcode and molecular identifiers from each
                   |read.  The observed sample barcode will be matched to the sample barcodes extracted from the bases in the sample sheet
                   |and associated read structures.  Please see the following link for the read structure:
                   |  https://broadinstitute.github.io/picard/javadoc/picard/picard/illumina/parser/ReadStructure.html
                   |
                   |""",
  group = ClpGroups.Utilities)
class DemuxFastqs
(
  @arg(flag="f1", doc="Input fastq file (optionally gzipped) for the first read of paired end data.")
  val fastq1: List[PathToFastq],
  @arg(flag="f2", doc="Input fastq file (optionally gzipped) for the second read of paired end data.")
  val fastq2: List[PathToFastq],
  @arg(flag="i7", doc="Input fastq file (optionally gzipped) for the index read of the Illumina i7 sequencing primer.  This is typically the I1 FASTQ (index read one for dual-indexed data).")
  val fastq7: List[PathToFastq],
  @arg(flag="i5", doc="Input fastq file (optionally gzipped) for the index read of the Illumina i5 sequencing primer.  This is typically the I2 FASTQ (index read two for dual-indexed data).")
  val fastq5: List[PathToFastq],
  @arg(doc="Read structure for the first read of paired end data.  Set to \"1000T\" or some large value to have all bases found be template bases. See the DemuxFastqs help message for more details.")
  val rs1: String = "1000T",
  @arg(doc="Read structure for the second read of paired end data.  Set to \"1000T\" or some large value to have all bases found be template bases. See the DemuxFastqs help message for more details.")
  val rs2: String = "1000T",
  @arg(doc="Read structure for the index read of the Illumina i7 sequencing primer.  The total bases must match the length of the index read. See the DemuxFastqs help message for more details.")
  val rs7: String,
  @arg(doc="Read structure for the index read of the Illumina i5 sequencing primer.  The total bases must match the length of the index read.See the DemuxFastqs help message for more details.")
  val rs5: String,
  @arg(flag="ss", doc="The Sample Sheet (SampleSheet.csv).  This must include information about the i7/i5 index reads' ids and bases.")
  val sampleSheet: FilePath,
  @arg(flag="o", doc="The output directory, with a BAM per sample.")
  val output: DirPath,
  @arg(flag="m", doc="Per-barcode and per-lane metrics written to this file.")
  val metrics: FilePath,
  @arg(flag="u", doc="Output BAM file name for the unmatched records.")
  val unmatched: String,
  @arg(flag="t", doc="The SAM tag for the molecular barcode.")
  val umiTag: String = "RX",
  @arg(flag="q", doc="A value describing how the quality values are encoded in the fastq.  Either Solexa for pre-pipeline 1.3 style scores (solexa scaling + 66), Illumina for pipeline 1.3 and above (phred scaling + 64) or Standard for phred scaled scores with a character shift of 33.  If this value is not specified, the quality format will be detected automatically.")
  var qualityFormat: Option[FastqQualityFormat] = None,
  @arg(flag="pl", doc="The platform type (e.g. illumina, solid) to insert into the read group header") val platform: Option[String] = Some("ILLUMINA"),
  @arg(flag="cn", doc="The sequencing center from which the data originated") val sequencingCenter: Option[String] = None,
  @arg(flag="pi", doc="Predicted median insert size, to insert into the read group header") val predictedInsertSize: Option[Integer] = None,
  @arg(flag="pg", doc="Program group to insert into the read group header.") val programGroup: Option[String] = None,
  @arg(flag="pm", doc="Platform model to insert into the group header (free-form text providing further details of the platform/technology used)") val platformModel: Option[String] = None,
  @arg(flag="co", doc="Comment(s) to include in the merged output file's header.", minElements = 0) val comments: List[String] = Nil,
  @arg(flag="ds", doc="Inserted into the read group header") val description: Option[String] = None,
  @arg(flag="dt", doc="Date the run was produced, to insert into the read group header") val runDate: Option[Iso8601Date] = None,
  @arg(flag="so", doc="The sort order for the output sam/bam file.") val sortOrder: SortOrder = SortOrder.queryname,
  @arg(doc="Minimum quality allowed in the input fastq.  An exception will be thrown if a quality is less than this value.") val minQ: Int = 0,
  @arg(doc="Maximum quality allowed in the input fastq.  An exception will be thrown if a quality is greater than this value.") val maxQ: Int = SAMUtils.MAX_PHRED_SCORE,
  @arg(doc="If true and this is an unpaired fastq any occurance of '/1' will be removed from the end of a read name.") val stripUnpairedMateName: Boolean = false,
  @arg(doc="Allow (and ignore) empty lines") val allowAndIgnoreEmptyLines: Boolean = false,
  @arg(doc="Maximum mismatches for a barcode to be considered a match.") val maxMismatches: Int = 1,
  @arg(doc="Minimum difference between number of mismatches in the best and second best barcodes for a barcode to be considered a match.") val minMismatchDelta: Int = 2,
  @arg(doc="Maximum allowable number of no-calls in a barcode read before it is considered unmatchable.") val maxNoCalls: Int = 2,
  @arg(doc="Minimum base quality. Any barcode bases falling below this quality will be considered a mismatch even in the bases match.") val minBaseQuality: Int = 0
) extends JeanLucTool with LazyLogging {

  import DemuxFastqs._

  override def execute(): Unit = {
    // Check all input and output files are readable or writable respectively.
    Io.assertReadable(fastq1 ++ fastq2 ++ fastq7 ++ fastq5)
    Io.assertReadable(this.sampleSheet)
    Io.assertWritableDirectory(this.output)
    Io.assertCanWriteFile(this.metrics)

    // Make some data structures to help us on our way
    val sampleSheet: SampleSheet = new SampleSheet(this.sampleSheet)
    val readOneReadStructureInfo = new ReadStructureInfo(rs1)
    val readTwoReadStructureInfo = new ReadStructureInfo(rs2)
    val i7ReadStructureInfo      = new ReadStructureInfo(rs7)
    val i5ReadStructureInfo      = new ReadStructureInfo(rs5)

    if (readOneReadStructureInfo.templateCycles.length == 0) throw new IllegalArgumentException("No template bases found in READ_ONE_READ_STRUCTURE.")
    if (readTwoReadStructureInfo.templateCycles.length == 0) throw new IllegalArgumentException("No template bases found in READ_TWO_READ_STRUCTURE.")
    if (i7ReadStructureInfo.templateCycles.length > 0)       throw new IllegalArgumentException("Template bases not allowed in the i7 read.")
    if (i5ReadStructureInfo.templateCycles.length > 0)       throw new IllegalArgumentException("Template bases not allowed in the i5 read.")

    // Set the barcode bases for each sample
    setSampleBarcode(sampleSheet, readOneReadStructureInfo, readTwoReadStructureInfo, i7ReadStructureInfo, i5ReadStructureInfo)

    // Determine the quality format
    determineQualityFormat()

    // Create the files for writing
    val writers = createSamFileWriters(sampleSheet)

    // Open the input FASTQs for reading
    val reader = new ParallelFastqReader(fastq1, fastq2, fastq7, fastq5)

    // Read in the quadruples of FASTQ records
    val fastqConverter = new FastqConverter(sampleSheet,
      umiTag,
      qualityFormat.getOrElse(unreachable("Quality format not set")),
      stripUnpairedMateName,
      minQ,
      maxQ,
      maxMismatches,
      minMismatchDelta,
      maxNoCalls,
      minBaseQuality,
      readOneReadStructureInfo,
      readTwoReadStructureInfo,
      i7ReadStructureInfo,
      i5ReadStructureInfo)

    reader.foreach { records => fastqConverter.convertFastqsAndWriteRecords(records, writers) }

    // close them all
    reader.close()
    writers.foreach(_.close())

    // write the metrics
    val metrics: MetricsFile[ExtractIlluminaBarcodes.BarcodeMetric, Integer] = new MetricsFile[ExtractIlluminaBarcodes.BarcodeMetric, Integer]()
    fastqConverter.getBarcodeMetrics.foreach(metrics.addMetric)
    metrics.write(this.metrics.toFile)
  }

  private def determineQualityFormat(): Unit = {
    val readers = fastq1.map { fastq => new FastqReader(fastq.toFile, allowAndIgnoreEmptyLines) }
    val detector: QualityEncodingDetector = new QualityEncodingDetector
    detector.add(QualityEncodingDetector.DEFAULT_MAX_RECORDS_TO_ITERATE, readers:_*)
    readers.foreach(_.close())
    val format = detector.generateBestGuess(QualityEncodingDetector.FileContext.FASTQ, qualityFormat.orNull)
    if (detector.isDeterminationAmbiguous) {
      logger.warning("Making ambiguous determination about fastq's quality encoding; more than one format possible based on observed qualities.")
    }
    logger.info(String.format("Auto-detected quality format as: %s.", format))
    this.qualityFormat = Some(format)
  }

  /** Creates SAM file writers for each sample in the Sample Sheet, returning them in the same order as specified
    * by each sample's ordinal.  An extra writer is appended for records that do not match a sample barcode. */
  private def createSamFileWriters(sampleSheet: SampleSheet): List[SAMFileWriter] = {
    val writers: List[SAMFileWriter] = sampleSheet.map { sample =>
      val outputFile: File = new File(this.output.toFile, IOUtil.makeFileNameSafe(String.format("%s.%s.%s.bam", sample.sampleId, sample.sampleName, sample.sampleBarcode)))
      val header: SAMFileHeader = createSamFileHeader(sample)
      val writer: SAMFileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, outputFile)
      writer
    }.toList
    val unmatchedWriter = {
      val outputFile: File = new File(this.output.toFile, this.unmatched)
      val header: SAMFileHeader = createSamFileHeader(new Sample(writers.length - 1, DemuxFastqs.UnmatchedSampleId, DemuxFastqs.UnmatchedSampleId, Some(DemuxFastqs.UnmatchedSampleId), None, None))
      val writer: SAMFileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, outputFile)
       writer
    }
    writers ++ List(unmatchedWriter)
  }

  /** Creates a read group for the SAM header with the values provided on the command line.. */
  private def createSamReadGroupRecord(readGroupId: String, sampleName: String, libraryName: String, platformUnit: String): SAMReadGroupRecord = {
    val rgroup: SAMReadGroupRecord = new SAMReadGroupRecord(readGroupId)
    rgroup.setSample(sampleName)
    rgroup.setLibrary(libraryName)
    rgroup.setPlatformUnit(platformUnit)
    platform.foreach(rgroup.setPlatform)
    sequencingCenter.foreach(rgroup.setSequencingCenter)
    predictedInsertSize.foreach(rgroup.setPredictedMedianInsertSize)
    description.foreach(rgroup.setDescription)
    runDate.foreach(rgroup.setRunDate)
    platformModel.foreach(rgroup.setPlatformModel)
    programGroup.foreach(rgroup.setProgramGroup)
    rgroup
  }

  /** Creates a simple header with the values provided on the command line. */
  private def createSamFileHeader(sample: Sample): SAMFileHeader = {
    val header: SAMFileHeader = new SAMFileHeader
    val libraryId: String = sample.libraryId match {
      case Some(id) => id
      case None => sample.sampleName
    }
    val platformUnit: String = sample.sampleBarcode match {
      case None => sample.sampleName
      case Some(sampleBarcode) => s"${sample.sampleName}.$sampleBarcode"
    }
    header.addReadGroup(createSamReadGroupRecord(sample.sampleId, sample.sampleName, libraryId, platformUnit))
    this.comments.foreach(header.addComment)
    header.setSortOrder(this.sortOrder)
    header
  }

}