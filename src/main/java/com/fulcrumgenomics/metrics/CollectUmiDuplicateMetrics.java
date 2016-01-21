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

package com.fulcrumgenomics.metrics;

import htsjdk.samtools.*;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.*;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Metrics;
import picard.sam.DuplicationMetrics;
import picard.util.IlluminaUtil;

import java.io.File;
import java.util.*;

// Developer note: it would be fairly trivial to output a BAM file (or BAM files) with duplicates marked
// for each mismatch threshold.
@CommandLineProgramProperties(
        usage = CollectUmiDuplicateMetrics.USAGE,
        usageShort = CollectUmiDuplicateMetrics.USAGE_SHORT,
        programGroup = Metrics.class
)
public class CollectUmiDuplicateMetrics extends CommandLineProgram {

    public final static String USAGE_SHORT = "Calculates duplication metrics for unique molecular identifiers with mismatch tolerances.";

    public final static String USAGE = USAGE_SHORT;

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "A SAM or BAM file to process.")
    public File INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The metrics file to write.")
    public File OUTPUT;

    @Option(doc = "Unique molecular identifier (UMI) SAM tag (ex. RX).")
    public String BARCODE_TAG = "RX";

    @Option(doc = "If multiple sequences are stored in the UMI SAM tag, use the ith one (0-based).", optional = true)
    public int BARCODE_TAG_INDEX = -1;

    @Option(doc = "Maximum mismatches for the UMIs to be considered equivalent.")
    public int MAX_MISMATCHES = 1;

    private final Log log = Log.getInstance(CollectUmiDuplicateMetrics.class);

    /** An enum to provide type-safe keys for transient attributes this tool puts on SAMRecords. */
    private enum Attr {
        BarcodeTag
    }

    private static char BARCODE_DELIMITER = IlluminaUtil.BARCODE_DELIMITER.charAt(0);

    @Override
    public int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        // the main loop
        final UmiDuplicateMetricsCollector metrics = new UmiDuplicateMetricsCollector(BARCODE_TAG, BARCODE_TAG_INDEX, MAX_MISMATCHES);
        final SamReader reader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);
        final DuplicateSetIterator iterator = new DuplicateSetIterator(reader.iterator(), reader.getFileHeader());
        final ProgressLogger progress = new ProgressLogger(log, 10000);
        for (final DuplicateSet duplicateSet : new IterableAdapter<>(iterator)) {
            metrics.add(duplicateSet);
            for (final SAMRecord record : duplicateSet.getRecords()) {
                progress.record(record);
            }
        }
        CloserUtil.close(reader);

        // finalize the metrics
        metrics.finalizeMetrics();

        // write the metrics file
        final MetricsFile<DuplicationMetrics, ?> metricsFile = getMetricsFile();
        metrics.addToMetricsFile(metricsFile);
        metricsFile.write(OUTPUT);

        return 0;
    }

    /** Collects duplication metrics for various scenarios: not considering the barcode tag, and considering the barcode
     * tag up to a certain mismatch threshold.
     */
    class UmiDuplicateMetricsCollector {

        private final String barcodeTag;
        private final int barcodeTagIndex;
        private final int maxMismatches;
        private final DuplicationMetrics noBarcodeTagMetrics = new DuplicationMetrics();
        private final Map<Integer, DuplicationMetrics> withMismatchesMetrics = new HashMap<>();

        public UmiDuplicateMetricsCollector(final String barcodeTag, final int barcodeTagIndex, final int maxMismatches) {
            this.barcodeTag      = barcodeTag;
            this.barcodeTagIndex = barcodeTagIndex;
            this.maxMismatches   = maxMismatches;
            noBarcodeTagMetrics.LIBRARY = "NO_BARCODE_TAG";
            for (int i = 0; i <= maxMismatches; i++) {
                final DuplicationMetrics metrics = new DuplicationMetrics();
                metrics.LIBRARY = String.format("%d_MISMATCHES", i);
                withMismatchesMetrics.put(i, metrics);
            }
        }

        public void add(final DuplicateSet duplicateSet) {
            // get the list of records on which to do analysis
            final List<SAMRecord> records = new ArrayList<>(duplicateSet.size());
            for (final SAMRecord record : duplicateSet.getRecords()) {
                if (!record.isSecondaryOrSupplementary()) {
                    // First bring the simple metrics up to date
                    if (record.getReadUnmappedFlag()) {
                        ++noBarcodeTagMetrics.UNMAPPED_READS;
                        withMismatchesMetrics.values().stream().forEach(metric -> metric.UNMAPPED_READS++);

                    } else if (!record.getReadPairedFlag() || record.getMateUnmappedFlag()) {
                        ++noBarcodeTagMetrics.UNPAIRED_READS_EXAMINED;
                        withMismatchesMetrics.values().stream().forEach(metric -> metric.UNPAIRED_READS_EXAMINED++);
                    } else {
                        // will need to be divided by 2 at the end
                        ++noBarcodeTagMetrics.READ_PAIRS_EXAMINED;
                        withMismatchesMetrics.values().stream().forEach(metric -> metric.READ_PAIRS_EXAMINED++);
                    }

                    // record the duplicates in noBarcodeTagMetrics
                    if (record.getDuplicateReadFlag()) {
                        // Update the duplication metrics
                        if (!record.getReadPairedFlag() || record.getMateUnmappedFlag()) {
                            ++noBarcodeTagMetrics.UNPAIRED_READ_DUPLICATES;
                        } else {
                            ++noBarcodeTagMetrics.READ_PAIR_DUPLICATES;// will need to be divided by 2 at the end
                        }
                    }

                    if (!record.getReadUnmappedFlag()) records.add(record);
                }
            }

            // update and count the duplicates
            for (int numMismatchesAllowed = 0; numMismatchesAllowed <= maxMismatches; numMismatchesAllowed++) {
                final Set<Set<SAMRecord>> sets = markDuplicatesWithBarcodeTag(records, barcodeTag, barcodeTagIndex, numMismatchesAllowed);
                final DuplicationMetrics metrics = withMismatchesMetrics.get(numMismatchesAllowed);

                // Update the duplication metrics
                sets.stream().forEach(set ->
                        set.stream().filter(SAMRecord::getDuplicateReadFlag).forEach(record -> {
                            // Update the duplication metrics
                            if (!record.getReadPairedFlag() || record.getMateUnmappedFlag()) {
                                ++metrics.UNPAIRED_READ_DUPLICATES;
                            } else {
                                ++metrics.READ_PAIR_DUPLICATES;// will need to be divided by 2 at the end
                            }
                        })
                );
            }
        }

        public void finalizeMetrics() {
            // We double counted read pairs
            noBarcodeTagMetrics.READ_PAIRS_EXAMINED  = noBarcodeTagMetrics.READ_PAIRS_EXAMINED / 2;
            noBarcodeTagMetrics.READ_PAIR_DUPLICATES = noBarcodeTagMetrics.READ_PAIR_DUPLICATES / 2;
            noBarcodeTagMetrics.calculateDerivedMetrics();
            for (final DuplicationMetrics metrics : withMismatchesMetrics.values()) {
                metrics.READ_PAIRS_EXAMINED  = metrics.READ_PAIRS_EXAMINED / 2;
                metrics.READ_PAIR_DUPLICATES = metrics.READ_PAIR_DUPLICATES / 2;
                metrics.calculateDerivedMetrics();
            }
        }

        public void addToMetricsFile(final MetricsFile<DuplicationMetrics, ?> metricsFile) {
            metricsFile.addMetric(this.noBarcodeTagMetrics);
            for (int i = 0; i <= maxMismatches; i++) {
                metricsFile.addMetric(withMismatchesMetrics.get(i));
            }
        }
    }

    /**
     * Partitions the duplicate set by the barcode tag.
     *
     * For each record within a set, we guarantee that there is at least one other record in that set that
     * has the same barcode given the mismatch tolerance.
     *
     * For two records in different sets, we guarantee that the number of mismatches is above the tolerance.
     *
     * Please note that there may be (and likely will be) two records in the same set whose barcodes exceed the
     * mismatch tolerance.  This is because for either record, there is another record in that set whose barcode
     * does not exceed the mismatch tolerance when compared.
     *
     * The duplicate flags will be updated.
     */
    private Set<Set<SAMRecord>> markDuplicatesWithBarcodeTag(final List<SAMRecord> records, final String barcodeTag, final int barcodeTagIndex, final int maxMismatches) {
        final Set<Set<SAMRecord>> sets = new HashSet<>();

        for (final SAMRecord record : records) {
            record.setDuplicateReadFlag(false);

            // find the sets to which this belongs.  If it can belong to more than one set, merge them.
            final Iterator<Set<SAMRecord>> iterator = sets.iterator();
            Set<SAMRecord> prevSet = null;
            while (iterator.hasNext()) {
                // determine if any record within this set can match the given record
                final Set<SAMRecord> set = iterator.next();
                for (final SAMRecord setRecord : set) {
                    if (compareBarcodes(record, setRecord, barcodeTag, barcodeTagIndex, maxMismatches) == 0) {
                        if (prevSet == null) {
                            // we found our first matching set
                            set.add(record);
                            prevSet = set;
                        }
                        else {
                            // merge the previous set and this set, and remove this set from the list of sets
                            prevSet.addAll(set);
                            iterator.remove();
                        }
                        // no need to continue, since we found a comparable record
                        break;
                    }
                }
            }
            // if we did not find a set, create a new one
            if (prevSet == null) {
                final Set<SAMRecord> set = new HashSet<>();
                set.add(record);
                sets.add(set);
            }
        }

        // set the duplicate flag for all but the first record in the set
        for (final Set<SAMRecord> set : sets) {
            boolean firstRecord = true;
            for (final SAMRecord record : set) {
                if (firstRecord) record.setDuplicateReadFlag(false);
                else record.setDuplicateReadFlag(true);
                firstRecord = false;
            }
        }

        return sets;
    }

    /**
     * Gets the barcode tag.  If `barcodeTagIndex`, the barcode tag will be split by '-' delimiter and the ith (0-based)
     * token will be returned.
     */
    private static String getBarcodeTag(final SAMRecord record, final String barcodeTag, final int barcodeTagIndex) {
        final String tagValue = record.getStringAttribute(barcodeTag);
        if (barcodeTagIndex < 0) return tagValue;
        final String[] tokens = tagValue.split(IlluminaUtil.BARCODE_DELIMITER);
        if (tokens.length <= barcodeTagIndex) {
            throw new PicardException(String.format("Not enough tokens (%d) for SAM tag (%s) in record: %s", tokens.length, barcodeTag, record.getSAMString()));
        }
        return tokens[barcodeTagIndex];
    }

    /**
     * Populates the set of transient attributes on SAMRecords if they are not already there.
     */
    private void populateTransientAttributes(final String barcodeTag, final int barcodeTagIndex, final SAMRecord... recs) {
        for (final SAMRecord rec : recs) {
            if (rec.getTransientAttribute(Attr.BarcodeTag) == null) {
                rec.setTransientAttribute(Attr.BarcodeTag, getBarcodeTag(rec, barcodeTag, barcodeTagIndex));
            }
        }
    }

    /**
     * Compares the values in the barcode tags.  This will also populate transient tags.
     */
    private int compareBarcodes(final SAMRecord samRecord1, final SAMRecord samRecord2, final String barcodeTag, final int barcodeTagIndex, final int maxMismatches) {
        populateTransientAttributes(barcodeTag, barcodeTagIndex, samRecord1, samRecord2);
        final String barcodeTagString1 = (String) samRecord1.getTransientAttribute(Attr.BarcodeTag);
        final String barcodeTagString2 = (String) samRecord2.getTransientAttribute(Attr.BarcodeTag);
        if (barcodeTagString1.length() != barcodeTagString2.length()) {
            throw new PicardException(String.format("Barcode lengths did not match for records '%s' and '%s'",
                    samRecord1.getReadName(), samRecord2.getReadName()));
        }
        return compareBarcodes(barcodeTagString1, barcodeTagString2, maxMismatches);
    }

    /**
     * Compares the values in the barcode strings.  Does not verify the lengths are the same, but they should be.
     *
     * Ignores any character that matches the delimiter in either string.
     */
    static int compareBarcodes(final String barcodeTagString1,
                               final String barcodeTagString2,
                               final int maxMismatches) {
        final byte[] barcodeTag1 = barcodeTagString1.getBytes();
        final byte[] barcodeTag2 = barcodeTagString2.getBytes();
        for (int i = 0, numMismatches = 0; i < barcodeTag1.length; i++) {
            if (!SequenceUtil.basesEqual(barcodeTag1[i], barcodeTag2[i]) &&
                    BARCODE_DELIMITER != barcodeTag1[i] &&
                    BARCODE_DELIMITER != barcodeTag2[i]) {
                numMismatches++;
                if (numMismatches > maxMismatches) {
                    final int cmp = barcodeTagString1.compareTo(barcodeTagString2);
                    if (cmp < 0) return -1;
                    else if (cmp == 0) return 0;
                    else return 1;
                }
            }
        }
        return 0;
    }
}
