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

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.SequenceUtil;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Metrics;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

/**
 * Computes metrics and histograms for assessing the accuracy of quality scores.
 *
 * The histograms show for each predicted quality score, the (computed) observed quality score.  We
 * also have histograms for each nucleotide as follows:
 * 1. a histogram for quality score assessment when the base call was X
 * 2. a histogram for quality score assessment when the reference nucleotide was X
 * for X in {A,C,G,T}.
 */
@CommandLineProgramProperties(
        usage = "Calculates the accuracy of quality scores.",
        usageShort = "Calculates the accuracy of quality scores.",
        programGroup = Metrics.class
)
public class CalculateQualityScoreAccuracy extends CommandLineProgram {
    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "A SAM or BAM file to process.")
    public File INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc = "The metrics file to write.")
    public File OUTPUT;

    @Option(doc="Consider only this cycle (zero-based). Set to -1 to examine all cycles.", optional = true)
    public int CYCLE = -1;

    @Option(doc="The number of bases at the start of the read not in the SAM/BAM file due to barcodes or clipping (ex. when running with bowtie).", optional=true)
    public int CYCLE_OFFSET = 0;

    /** Stock main method for a command line program. */
    public static void main(final String[] argv) {
        new CalculateQualityScoreAccuracy().instanceMainWithExit(argv);
    }

    /**
     * Main method for the program.  Checks that all input files are present and
     * readable and that the output file can be written to.  Then iterates through
     * all the records accumulating metrics.  Finally writes metrics file
     */
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        final SamReader reader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);

        final Histogram<Integer> mismatchesHist = new Histogram<Integer>("Predicted", "Mismatches");
        final Histogram<Integer> totalHist = new Histogram<Integer>("Predicted", "Total_Bases");
        final Map<String,Histogram> mismatchesByTypeHist = new HashMap<String, Histogram>();
        final Map<String,Histogram> totalByTypeHist = new HashMap<String, Histogram>();

        // Set up the histograms
        byte[] bases = {'A', 'C', 'G', 'T'};
        for (final byte base: bases) {
            final Histogram<Integer> h = new Histogram<Integer>("Predicted", (char)base + ">");
            mismatchesByTypeHist.put((char)base + ">", h);
            final Histogram<Integer> h2 = new Histogram<Integer>("Predicted", ">" + (char)base);
            mismatchesByTypeHist.put(">" + (char)base, h2);
        }
        for (final byte base: bases) {
            final Histogram<Integer> h = new Histogram<Integer>("Predicted", (char)base + ">");
            totalByTypeHist.put((char)base + ">", h);
            final Histogram<Integer> h2 = new Histogram<Integer>("Predicted", ">" + (char)base);
            totalByTypeHist.put(">" + (char)base, h2);
        }

        for (final SAMRecord record : reader) {
            // Ignore these as we don't know the truth
            if (record.getReadUnmappedFlag() || record.isSecondaryOrSupplementary()) {
                continue;
            }
            final byte[] readBases = record.getReadBases();
            final byte[] readQualities = record.getBaseQualities();
            final byte[] refBases = SequenceUtil.makeReferenceFromAlignment(record, false);

            // We've seen stranger things
            if (readQualities.length != readBases.length) {
                throw new PicardException("Missing Qualities ("
                        + readQualities.length
                        + "," + readBases.length
                        + ") : "
                        + record.getSAMString());
            }

            if (refBases.length != readBases.length) {
                throw new PicardException("The read length did not match the inferred reference length, please check your MD and CIGAR.");
            }

            int cycleIndex; // zero-based
            if (record.getReadNegativeStrandFlag()) {
                cycleIndex = readBases.length - 1 + CYCLE_OFFSET;
            }
            else {
                cycleIndex = CYCLE_OFFSET;
            }

            for (int i = 0; i < readBases.length; i++) {
                if (-1 == CYCLE || cycleIndex == CYCLE) {
                    if ('-' != refBases[i] && '0' != refBases[i]) { // not insertion and not soft-clipped
                        if (!SequenceUtil.basesEqual(readBases[i], refBases[i])) { // mismatch
                            mismatchesHist.increment((int) readQualities[i]);
                            if (SequenceUtil.isValidBase(refBases[i])) {
                                mismatchesByTypeHist.get((char) refBases[i] + ">").increment((int)readQualities[i]);
                            }
                            if (SequenceUtil.isValidBase(readBases[i])) {
                                mismatchesByTypeHist.get(">" + (char) readBases[i]).increment((int)readQualities[i]);
                            }
                        } else {
                            mismatchesHist.increment((int) readQualities[i], 0); // to make sure the bin will exist
                        }
                        totalHist.increment((int) readQualities[i]);
                        if (SequenceUtil.isValidBase(readBases[i])) {
                            totalByTypeHist.get(">" + (char) readBases[i]).increment((int)readQualities[i]);
                        }
                        if (SequenceUtil.isValidBase(refBases[i])) {
                            totalByTypeHist.get((char) refBases[i] + ">").increment((int)readQualities[i]);
                        }
                    }
                }
                cycleIndex += record.getReadNegativeStrandFlag() ? -1 : 1;
            }
        }
        CloserUtil.close(reader);

        final Histogram<Integer> hist = new Histogram<Integer>("Predicted", "Observed");

        double sumOfSquaresError = 0.0;

        // compute the aggregate phred values
        for (final Integer key : mismatchesHist.keySet()) {
            final double numMismatches = mismatchesHist.get(key).getValue();
            final double numBases = totalHist.get(key).getValue();
            final double phredErr = Math.log10(numMismatches / numBases) * -10.0;
            sumOfSquaresError += (0 == numMismatches) ? 0.0 : (key - phredErr) * (key - phredErr);
            hist.increment(key, phredErr);

            // make sure the bin will exist
            for (final byte base : bases) {
                mismatchesByTypeHist.get(">" + (char)base).increment(key, 0.0);
                mismatchesByTypeHist.get((char)base + ">").increment(key, 0.0);
                totalByTypeHist.get(">" + (char)base).increment(key, 0.0);
                totalByTypeHist.get((char)base + ">").increment(key, 0.0);
            }
        }

        final QualityScoreAccuracyMetrics metrics = new QualityScoreAccuracyMetrics();
        metrics.SUM_OF_SQUARE_ERROR = sumOfSquaresError;

        final MetricsFile<QualityScoreAccuracyMetrics, Integer> out = getMetricsFile();
        out.addMetric(metrics);
        out.addHistogram(hist);
        for (final byte base : bases) {
            // >base : histograms for mismatches *to* the given base
            Histogram<Integer> m = mismatchesByTypeHist.get(">" + (char)base);
            Histogram<Integer> t = totalByTypeHist.get(">" + (char)base);
            Histogram<Integer> h = new Histogram<Integer>(m.getBinLabel(), m.getValueLabel());
            for (final Integer key : m.keySet()) {
                final double numMismatches = m.get(key).getValue();
                final double numBases = t.get(key).getValue();
                final double phredErr = Math.log10(numMismatches / numBases) * -10.0;
                h.increment(key, phredErr);
            }
            out.addHistogram(h);

            // base> : histograms for mismatches *from* the given base
            m = mismatchesByTypeHist.get((char)base + ">");
            t = totalByTypeHist.get(">" + (char)base);
            h = new Histogram<Integer>(m.getBinLabel(), m.getValueLabel());
            for (final Integer key : m.keySet()) {
                final double numMismatches = m.get(key).getValue();
                final double numBases = t.get(key).getValue();
                final double phredErr = Math.log10(numMismatches / numBases) * -10.0;
                h.increment(key, phredErr);
            }
            out.addHistogram(h);
        }

        out.addHistogram(mismatchesHist);
        out.addHistogram(totalHist);
        out.write(OUTPUT);

        return 0;
    }

    public static class QualityScoreAccuracyMetrics extends MetricBase {
        /** The sum of square error for the base qualities */
        public double SUM_OF_SQUARE_ERROR;
    }
}
