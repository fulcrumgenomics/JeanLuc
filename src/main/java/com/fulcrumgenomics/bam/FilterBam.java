/*
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
package com.fulcrumgenomics.bam;

import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.SamOrBam;

import java.io.File;
import java.text.DecimalFormat;

/**
 * Program which takes in a BAM file and filters out all reads for templates that match one or more
 * criteria.  Designed to be used to filter out reads that might confuse variant callers and lead
 * to false positive variant calls.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        usage = "Filters a BAM file to remove reads that may not be useful in downstream processing, in order\n" +
                "to reduce the size of the file. By default will remove unmapped reads, read with MAPQ=0, records\n" +
                "marked as secondary alignments, records marked as duplicates, and if a set of Intervals are provided\n" +
                "records that do not overlap any of the intervals.\n\n" +
                "NOTE: this will usually produce a BAM file in which some mate-pairs are orphaned (i.e. read 1 or\n" +
                "read 2 is included, but not both), but does not update any flag fields.",
        usageShort = "Filters reads out of a BAM file.",
        programGroup = SamOrBam.class
)
public class FilterBam extends CommandLineProgram {
    @Option(doc="If true remove all reads that are marked as duplicates.")
    public boolean REMOVE_DUPLICATES = true;

    @Option(doc="Remove all unmapped reads.")
    public boolean REMOVE_UNMAPPED_READS = true;

    @Option(doc="Remove all reads with MAPQ lower than this number.")
    public int MIN_MAPQ = 1;

    @Option(doc="Remove all reads marked as secondary alignments.")
    public boolean REMOVE_SECONDARY_ALIGNMENTS = true;

    @Option(doc="If supplied, remove all reads that do not overlap the provided intervals.", optional=true)
    public File INTERVALS;

    @Option(doc="Input BAM file.", shortName= StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Option(doc="Output BAM file.", shortName= StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;

    private final Log log = Log.getInstance(FilterBam.class);

    // Default index creation to on
    public FilterBam() { this.CREATE_INDEX = true; }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        if (INTERVALS != null) IOUtil.assertFileIsReadable(INTERVALS);

        final ProgressLogger progress = new ProgressLogger(log);

        final SamReader in = SamReaderFactory.make().open(INPUT);
        final SAMRecordIterator iterator = buildInputIterator(in, INTERVALS);
        final SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(in.getFileHeader(), true, OUTPUT);

        long kept = 0;
        while (iterator.hasNext()) {
            final SAMRecord rec = iterator.next();
            progress.record(rec);

            if (REMOVE_DUPLICATES && rec.getDuplicateReadFlag()) continue;
            if (REMOVE_UNMAPPED_READS && rec.getReadUnmappedFlag()) continue;
            if (rec.getMappingQuality() < MIN_MAPQ) continue;
            if (REMOVE_SECONDARY_ALIGNMENTS && !rec.getReadUnmappedFlag() && rec.getNotPrimaryAlignmentFlag()) continue;

            ++kept;
            out.addAlignment(rec);
        }

        log.info("Kept " + new DecimalFormat("#,##0").format(kept) + " records.");
        out.close();
        CloserUtil.close(iterator);
        return 0;
    }

    /**
     * If intervalListFile is null return an interator over all the input, otherwise returns an
     * iterator over only those reads that overlap one or more of the intervals in the file.
     */
    protected SAMRecordIterator buildInputIterator(final SamReader in, final File intervalListFile) {
        if (intervalListFile == null) {
            return in.iterator();
        }
        else {
            IntervalList intervals = IntervalList.fromFile(intervalListFile).uniqued();
            final SAMSequenceDictionary dict = intervals.getHeader().getSequenceDictionary();
            final QueryInterval[] qs = intervals.getIntervals().stream()
                    .map(i -> new QueryInterval(dict.getSequenceIndex(i.getContig()), i.getStart(), i.getEnd()))
                    .toArray(QueryInterval[]::new);
            return in.queryOverlapping(qs);
        }
    }
}
