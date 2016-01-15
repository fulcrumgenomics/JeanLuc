/**
 * Copyright (c) 2015, Fulcrum Genomics LLC
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

package com.fulcrumgenomics.metrics;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.filter.SecondaryAlignmentFilter;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.StringUtil;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Metrics;
import picard.util.MathUtil;
import picard.util.RExecutor;

import java.io.File;
import java.util.*;

/**
 * Computes a number of metrics that are useful for evaluating coverage and performance of whole genome sequencing experiments.
 *
 * The main difference between this tool and Picard's CollectWgsMetrics, is this one produces a histogram to view, as well as adds a few metrics.
 *
 * @author nhomer
 */
@CommandLineProgramProperties(
        usage = "Computes a number of metrics that are useful for evaluating coverage and performance of " +
                "whole genome sequencing experiments.",
        usageShort = "Writes whole genome sequencing-related metrics for a SAM or BAM file",
        programGroup = Metrics.class
)
public class CollectWgsMetricsWithHistogram extends CommandLineProgram {

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input SAM or BAM file.")
    public File INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output metrics file.")
    public File OUTPUT;

    @Option(shortName = "CHART", doc = "A file (with .pdf extension) to write the chart to.")
    public File CHART_OUTPUT;

    @Option(shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc="The reference sequence fasta aligned to.")
    public File REFERENCE_SEQUENCE;

    @Option(shortName = "MQ", doc = "Minimum mapping quality for a read to contribute coverage.", overridable = true)
    public int MINIMUM_MAPPING_QUALITY = 20;

    @Option(shortName = "Q", doc = "Minimum base quality for a base to contribute coverage.", overridable = true)
    public int MINIMUM_BASE_QUALITY = 20;

    @Option(shortName = "CAP", doc = "Treat bases with coverage exceeding this value as if they had coverage at this value.", overridable = true)
    public int COVERAGE_CAP = 250;

    @Option(doc = "For debugging purposes, stop after processing this many genomic bases.")
    public long STOP_AFTER = -1;

    @Option(doc = "Determines whether to include the base quality histogram in the metrics file.")
    public boolean INCLUDE_BQ_HISTOGRAM = false;

    @Option(doc="If true, count unpaired reads, and paired reads with one end unmapped")
    public boolean COUNT_UNPAIRED = false;

    private final Log log = Log.getInstance(CollectWgsMetricsWithHistogram.class);

    /** Metrics for evaluating the performance of whole genome sequencing experiments. */
    public static class WgsMetrics extends MetricBase {

        public enum Category { WHOLE_GENOME, NON_ZERO_REGIONS }

        /** One of either WHOLE_GENOME or NON_ZERO_REGIONS */
        public Category CATEGORY;

        /** The number of non-N bases in the genome reference over which coverage will be evaluated. */
        public long GENOME_TERRITORY;
        /** The mean coverage in bases of the genome territory, after all filters are applied. */
        public double MEAN_COVERAGE;
        /** The standard deviation of coverage of the genome after all filters are applied. */
        public double SD_COVERAGE;
        /** The median coverage in bases of the genome territory, after all filters are applied. */
        public double MEDIAN_COVERAGE;
        /** The median absolute deviation of coverage of the genome after all filters are applied. */
        public double MAD_COVERAGE;

        /** The fraction of aligned bases that were filtered out because they were in reads with low mapping quality (default is < 20). */
        public double PCT_EXC_MAPQ;
        /** The fraction of aligned bases that were filtered out because they were in reads marked as duplicates. */
        public double PCT_EXC_DUPE;
        /** The fraction of aligned bases that were filtered out because they were in reads without a mapped mate pair. */
        public double PCT_EXC_UNPAIRED;
        /** The fraction of aligned bases that were filtered out because they were of low base quality (default is < 20). */
        public double PCT_EXC_BASEQ;
        /** The fraction of aligned bases that were filtered out because they were the second observation from an insert with overlapping reads. */
        public double PCT_EXC_OVERLAP;
        /** The fraction of aligned bases that were filtered out because they would have raised coverage above the capped value (default cap = 250x). */
        public double PCT_EXC_CAPPED;
        /** The total fraction of aligned bases excluded due to all filters. */
        public double PCT_EXC_TOTAL;

        /** The fraction of bases that attained at least 1X sequence coverage in post-filtering bases. */
        public double PCT_1X;
        /** The fraction of bases that attained at least 5X sequence coverage in post-filtering bases. */
        public double PCT_5X;
        /** The fraction of bases that attained at least 10X sequence coverage in post-filtering bases. */
        public double PCT_10X;
        /** The fraction of bases that attained at least 15X sequence coverage in post-filtering bases. */
        public double PCT_15X;
        /** The fraction of bases that attained at least 20X sequence coverage in post-filtering bases. */
        public double PCT_20X;
        /** The fraction of bases that attained at least 25X sequence coverage in post-filtering bases. */
        public double PCT_25X;
        /** The fraction of bases that attained at least 30X sequence coverage in post-filtering bases. */
        public double PCT_30X;
        /** The fraction of bases that attained at least 40X sequence coverage in post-filtering bases. */
        public double PCT_40X;
        /** The fraction of bases that attained at least 50X sequence coverage in post-filtering bases. */
        public double PCT_50X;
        /** The fraction of bases that attained at least 60X sequence coverage in post-filtering bases. */
        public double PCT_60X;
        /** The fraction of bases that attained at least 70X sequence coverage in post-filtering bases. */
        public double PCT_70X;
        /** The fraction of bases that attained at least 80X sequence coverage in post-filtering bases. */
        public double PCT_80X;
        /** The fraction of bases that attained at least 90X sequence coverage in post-filtering bases. */
        public double PCT_90X;
        /** The fraction of bases that attained at least 100X sequence coverage in post-filtering bases. */
        public double PCT_100X;
    }

    private class WgsMetricsCollector {

        private final long[] histogramArray;
        private final long[] baseQHistogramArray;
        private long basesExcludedByBaseq = 0;
        private long basesExcludedByOverlap = 0;
        private long basesExcludedByCapping = 0;
        private final int coverageCap;
        private final Map<WgsMetrics.Category, Histogram<Integer>> histograms = new HashMap<>();

        public WgsMetricsCollector(final int coverageCap) {
            histogramArray = new long[coverageCap + 1];
            baseQHistogramArray = new long[Byte.MAX_VALUE];
            this.coverageCap = coverageCap;
        }

        public void addInfo(final SamLocusIterator.LocusInfo info, final ReferenceSequence ref) {
            // Check that the reference is not N
            final byte base = ref.getBases()[info.getPosition() - 1];
            if (base == 'N') return;

            // Figure out the coverage while not counting overlapping reads twice, and excluding various things
            final HashSet<String> readNames = new HashSet<>(info.getRecordAndPositions().size());
            int pileupSize = 0;
            for (final SamLocusIterator.RecordAndOffset recs : info.getRecordAndPositions()) {

                if (recs.getBaseQuality() < MINIMUM_BASE_QUALITY)                   { ++basesExcludedByBaseq;   continue; }
                if (!readNames.add(recs.getRecord().getReadName()))                 { ++basesExcludedByOverlap; continue; }
                pileupSize++;
                if (pileupSize <= coverageCap) {
                    baseQHistogramArray[recs.getRecord().getBaseQualities()[recs.getOffset()]]++;
                }
            }

            final int depth = Math.min(readNames.size(), coverageCap);
            if (depth < readNames.size()) basesExcludedByCapping += readNames.size() - coverageCap;
            histogramArray[depth]++;
        }

        public void addMetricsToFile(final MetricsFile<WgsMetrics, Integer> file,
                                     final boolean includeBQHistogram,
                                     final CountingFilter dupeFilter,
                                     final CountingFilter mapqFilter,
                                     final CountingPairedFilter pairFilter) {
            addToMetricsFile(file, WgsMetrics.Category.WHOLE_GENOME, dupeFilter, mapqFilter, pairFilter);
            addToMetricsFile(file, WgsMetrics.Category.NON_ZERO_REGIONS, dupeFilter, mapqFilter, pairFilter);

            if (includeBQHistogram) {
                addBaseQHistogram(file);
            }
        }

        private void addBaseQHistogram(final MetricsFile<WgsMetrics, Integer> file) {
            // Construct and write the outputs
            final Histogram<Integer> baseQHistogram = new Histogram<>("value", "baseq_count");

            for (int i = 0; i < baseQHistogramArray.length; ++i) {
                baseQHistogram.increment(i, baseQHistogramArray[i]);
            }

            file.addHistogram(baseQHistogram);
        }

        private void addToMetricsFile(final MetricsFile<WgsMetrics, Integer> file,
                                      final WgsMetrics.Category category,
                                      final CountingFilter dupeFilter,
                                      final CountingFilter mapqFilter,
                                      final CountingPairedFilter pairFilter) {
            // Construct and write the outputs
            final Histogram<Integer> histogram = new Histogram<>("coverage", "count_" + category.name());

            for (int i=0; i<histogramArray.length; ++i) {
                if (category == WgsMetrics.Category.WHOLE_GENOME ||
                        (category == WgsMetrics.Category.NON_ZERO_REGIONS && 0 < i)) {
                    histogram.increment(i, histogramArray[i]);
                }
            }

            histograms.put(category, histogram);

            final WgsMetrics metrics = new WgsMetrics();
            metrics.CATEGORY = category;
            metrics.GENOME_TERRITORY = (long) histogram.getSumOfValues();
            metrics.MEAN_COVERAGE = histogram.getMean();
            metrics.SD_COVERAGE = histogram.getStandardDeviation();
            metrics.MEDIAN_COVERAGE = histogram.getMedian();
            metrics.MAD_COVERAGE = histogram.getMedianAbsoluteDeviation();

            final long basesExcludedByDupes   = dupeFilter.getFilteredBases();
            final long basesExcludedByMapq    = mapqFilter.getFilteredBases();
            final long basesExcludedByPairing = pairFilter.getFilteredBases();
            final double total             = histogram.getSum();
            final double totalWithExcludes = total + basesExcludedByDupes + basesExcludedByMapq + basesExcludedByPairing + basesExcludedByBaseq + basesExcludedByOverlap + basesExcludedByCapping;
            metrics.PCT_EXC_DUPE = basesExcludedByDupes / totalWithExcludes;
            metrics.PCT_EXC_MAPQ = basesExcludedByMapq / totalWithExcludes;
            metrics.PCT_EXC_UNPAIRED = basesExcludedByPairing / totalWithExcludes;

            metrics.PCT_EXC_BASEQ    = basesExcludedByBaseq   / totalWithExcludes;
            metrics.PCT_EXC_OVERLAP  = basesExcludedByOverlap / totalWithExcludes;
            metrics.PCT_EXC_CAPPED   = basesExcludedByCapping / totalWithExcludes;
            metrics.PCT_EXC_TOTAL    = (totalWithExcludes - total) / totalWithExcludes;

            metrics.PCT_1X     = MathUtil.sum(histogramArray, 1, histogramArray.length)   / (double) metrics.GENOME_TERRITORY;
            metrics.PCT_5X     = MathUtil.sum(histogramArray, 5, histogramArray.length)   / (double) metrics.GENOME_TERRITORY;
            metrics.PCT_10X    = MathUtil.sum(histogramArray, 10, histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
            metrics.PCT_15X    = MathUtil.sum(histogramArray, 15, histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
            metrics.PCT_20X    = MathUtil.sum(histogramArray, 20, histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
            metrics.PCT_25X    = MathUtil.sum(histogramArray, 25, histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
            metrics.PCT_30X    = MathUtil.sum(histogramArray, 30, histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
            metrics.PCT_40X    = MathUtil.sum(histogramArray, 40, histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
            metrics.PCT_50X    = MathUtil.sum(histogramArray, 50, histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
            metrics.PCT_60X    = MathUtil.sum(histogramArray, 60, histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
            metrics.PCT_70X    = MathUtil.sum(histogramArray, 70, histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
            metrics.PCT_80X    = MathUtil.sum(histogramArray, 80, histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
            metrics.PCT_90X    = MathUtil.sum(histogramArray, 90, histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
            metrics.PCT_100X   = MathUtil.sum(histogramArray, 100, histogramArray.length) / (double) metrics.GENOME_TERRITORY;

            file.addMetric(metrics);
            file.addHistogram(histogram);
        }

        public boolean areHistogramsEmpty() {
            return this.histograms.get(WgsMetrics.Category.WHOLE_GENOME).isEmpty();
        }
    }

    public static void main(final String[] args) {
        new CollectWgsMetricsWithHistogram().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
        IOUtil.assertFileIsWritable(CHART_OUTPUT);

        // Setup all the inputs
        final ProgressLogger progress = new ProgressLogger(log, 10000000, "Processed", "loci");
        final ReferenceSequenceFileWalker refWalker = new ReferenceSequenceFileWalker(REFERENCE_SEQUENCE);
        final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);
        final SamLocusIterator iterator = getLocusIterator(in);


        final List<SAMReadGroupRecord> readGroups = in.getFileHeader().getReadGroups();
        String plotSubtitle = "";
        if (readGroups.size() == 1) {
            plotSubtitle = StringUtil.asEmptyIfNull(readGroups.get(0).getLibrary());
        }


        final List<SamRecordFilter> filters   = new ArrayList<>();
        final CountingFilter dupeFilter       = new CountingDuplicateFilter();
        final CountingFilter mapqFilter       = new CountingMapQFilter(MINIMUM_MAPPING_QUALITY);
        final CountingPairedFilter pairFilter = new CountingPairedFilter();
        filters.add(mapqFilter);
        filters.add(dupeFilter);
        if (!COUNT_UNPAIRED) {
            filters.add(pairFilter);
        }
        filters.add(new SecondaryAlignmentFilter()); // Not a counting filter because we never want to count reads twice
        iterator.setSamFilters(filters);
        iterator.setEmitUncoveredLoci(true);
        iterator.setMappingQualityScoreCutoff(0); // Handled separately because we want to count bases
        iterator.setQualityScoreCutoff(0);        // Handled separately because we want to count bases
        iterator.setIncludeNonPfReads(false);

        final WgsMetricsCollector collector = new WgsMetricsCollector(COVERAGE_CAP);

        final boolean usingStopAfter = STOP_AFTER > 0;
        final long stopAfter = STOP_AFTER - 1;
        long counter = 0;

        // Loop through all the loci
        while (iterator.hasNext()) {
            final SamLocusIterator.LocusInfo info = iterator.next();
            final ReferenceSequence ref = refWalker.get(info.getSequenceIndex());
            // add to the collector
            collector.addInfo(info, ref);
            // Record progress and perhaps stop
            progress.record(info.getSequenceName(), info.getPosition());
            if (usingStopAfter && ++counter > stopAfter) break;
        }

        final MetricsFile<WgsMetrics, Integer> out = getMetricsFile();
        collector.addMetricsToFile(out, INCLUDE_BQ_HISTOGRAM, dupeFilter, mapqFilter, pairFilter);
        out.write(OUTPUT);

        if (collector.areHistogramsEmpty()) {
            log.warn("No valid bases found in input file. No plot will be produced.");
        } else {
            final int rResult = RExecutor.executeFromClasspath("scripts/com/fulcrumgenomics/metrics/wgsHistogram.R",
                    OUTPUT.getAbsolutePath(),
                    CHART_OUTPUT.getAbsolutePath(),
                    INPUT.getName(),
                    plotSubtitle);
            if (rResult != 0) {
                throw new PicardException("R script wgsHistogram.R failed with return code " + rResult);
            }
        }

        return 0;
    }
    protected WgsMetrics generateWgsMetrics() {
        return new WgsMetrics();
    }

    protected long getBasesExcludedBy(final CountingFilter filter) {
        return filter.getFilteredBases();
    }

    protected SamLocusIterator getLocusIterator(final SamReader in) {
        return new SamLocusIterator(in);
    }
}

/**
 * A SamRecordFilter that counts the number of aligned bases in the reads which it filters out. Abstract and designed
 * to be subclassed to implement the desired filter.
 */
abstract class CountingFilter implements SamRecordFilter {
    private long filteredRecords = 0;
    private long filteredBases = 0;

    /** Gets the number of records that have been filtered out thus far. */
    public long getFilteredRecords() { return this.filteredRecords; }

    /** Gets the number of bases that have been filtered out thus far. */
    public long getFilteredBases() { return this.filteredBases; }

    @Override
    public final boolean filterOut(final SAMRecord record) {
        final boolean filteredOut = reallyFilterOut(record);
        if (filteredOut) {
            ++filteredRecords;
            for (final AlignmentBlock block : record.getAlignmentBlocks()) {
                this.filteredBases += block.getLength();
            }
        }
        return filteredOut;
    }

    abstract public boolean reallyFilterOut(final SAMRecord record);

    @Override
    public boolean filterOut(final SAMRecord first, final SAMRecord second) {
        throw new UnsupportedOperationException();
    }
}

/** Counting filter that discards reads that have been marked as duplicates. */
class CountingDuplicateFilter extends CountingFilter {
    @Override
    public boolean reallyFilterOut(final SAMRecord record) { return record.getDuplicateReadFlag(); }
}

/** Counting filter that discards reads below a configurable mapping quality threshold. */
class CountingMapQFilter extends CountingFilter {
    private final int minMapq;

    CountingMapQFilter(final int minMapq) { this.minMapq = minMapq; }

    @Override
    public boolean reallyFilterOut(final SAMRecord record) { return record.getMappingQuality() < minMapq; }
}

/** Counting filter that discards reads that are unpaired in sequencing and paired reads who's mates are not mapped. */
class CountingPairedFilter extends CountingFilter {
    @Override
    public boolean reallyFilterOut(final SAMRecord record) { return !record.getReadPairedFlag() || record.getMateUnmappedFlag(); }
}

