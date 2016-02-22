/*
 * The MIT License
 *
 * Copyright (c) 2016 Fulcrum Genomics LLC
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package com.fulcrumgenomics.bam

import java.nio.file.Path
import java.text.DecimalFormat

import com.fulcrumgenomics.cmdline.{ClpGroups, JeanLucTool}
import dagr.commons.CommonsDef.{PathToIntervals, PathToBam}
import dagr.commons.io.Io
import dagr.commons.util.LazyLogging
import dagr.sopt._
import htsjdk.samtools._
import htsjdk.samtools.util._

import scala.collection.JavaConversions._

/**
  * Program which takes in a BAM file and filters out all reads for templates that match one or more
  * criteria.  Designed to be used to filter out reads that might confuse variant callers and lead
  * to false positive variant calls.
  *
  * @author Tim Fennell
  */
@clp(description = "Filters reads out of a BAM file. Remove reads that may not be useful in downstream processing, in order\n" +
  "to reduce the size of the file. By default will remove unmapped reads, read with MAPQ=0, records\n" +
  "marked as secondary alignments, records marked as duplicates, and if a set of Intervals are provided\n" +
  "records that do not overlap any of the intervals.\n\n" +
  "NOTE: this will usually produce a BAM file in which some mate-pairs are orphaned (i.e. read 1 or\n" +
  "read 2 is included, but not both), but does not update any flag fields.",
  group = ClpGroups.SamOrBam)
class FilterBam
( @arg(doc = "If supplied, remove all reads that do not overlap the provided intervals.") var intervals: Option[PathToIntervals] = None,
  @arg(doc = "Input BAM file.")                                         var input: PathToBam,
  @arg(doc = "Output BAM file.")                                        var output: PathToBam,
  @arg(doc = "If true remove all reads that are marked as duplicates.") var removeDuplicates: Boolean = true,
  @arg(doc = "Remove all unmapped reads.")                              var removeUnmappedReads: Boolean = true,
  @arg(doc = "Remove all reads with MAPQ lower than this number.")      var minMapQ: Int = 1,
  @arg(doc = "Remove all reads marked as secondary alignments.")        var removeSecondaryAlignments: Boolean = true
) extends JeanLucTool with LazyLogging {

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)
  intervals.foreach(Io.assertReadable)

  override def execute(): Unit = {
    //val progress: ProgressLogger = new ProgressLogger(log)
    val in: SamReader = SamReaderFactory.make.open(input.toFile)
    val iterator: SAMRecordIterator = buildInputIterator(in, intervals)
    val out: SAMFileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(in.getFileHeader, true, output.toFile)
    val kept: Long = iterator.filterNot { rec =>
      (removeDuplicates && rec.getDuplicateReadFlag) ||
        (removeUnmappedReads && rec.getReadUnmappedFlag) ||
        (rec.getMappingQuality < minMapQ) ||
        (removeSecondaryAlignments && !rec.getReadUnmappedFlag && rec.getNotPrimaryAlignmentFlag)
    }.map { rec =>
      out.addAlignment(rec)
      1.toLong
    }.sum[Long]
    logger.info("Kept " + new DecimalFormat("#,##0").format(kept) + " records.")
    out.close()
    CloserUtil.close(iterator)
  }
  /**
    * If intervalListFile is null return an interator over all the input, otherwise returns an
    * iterator over only those reads that overlap one or more of the intervals in the file.
    */
  protected def buildInputIterator(in: SamReader, intervalListFile: Option[Path]): SAMRecordIterator = {
    intervalListFile match {
      case None => in.iterator()
      case Some(file) =>
        val intervals: IntervalList = IntervalList.fromFile(file.toFile).uniqued
        val dict: SAMSequenceDictionary = intervals.getHeader.getSequenceDictionary
        val qs: Array[QueryInterval] = intervals.getIntervals.map(interval =>
          new QueryInterval(dict.getSequenceIndex(interval.getContig), interval.getStart, interval.getEnd)).toArray
        in.queryOverlapping(qs)
    }
  }
}
