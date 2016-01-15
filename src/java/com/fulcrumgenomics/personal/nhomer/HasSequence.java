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

package com.fulcrumgenomics.personal.nhomer;

import com.fulcrumgenomics.cmdline.Personal;
import htsjdk.samtools.*;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.*;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

@CommandLineProgramProperties(
        usage = "Searches for DNA sequences in the read pairs.",
        usageShort = "Searches for DNA sequences in the read pairs.",
        programGroup = Personal.class
)
public class HasSequence extends CommandLineProgram {

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "Input SAM or BAM.")
    public File INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output SAM or BAM.")
    public File OUTPUT;

    @Option(shortName = "S", doc = "File containing the DNA sequences to search for.")
    public File SEQUENCES;

    @Option(doc = "Sequence match histogram written to this file.", shortName = StandardOptionDefinitions.METRICS_FILE_SHORT_NAME)
    public File METRICS_FILE;

    @Option(doc = "Maximum mismatches for matching a sequence.")
    public int MAX_MISMATCHES = 1;

    @Option(shortName = "Q", doc = "Minimum base quality. Any bases falling below this quality will be considered a mismatch even in the bases match.")
    public int MINIMUM_BASE_QUALITY = 0;

    @Option(doc = "Treat no calls in the input sequences as DNA bases.")
    public boolean NO_CALLS_ARE_BASES = false;

    @Option(doc = "Search for mismatches at the start only")
    public boolean START_ONLY = true;

    final Log log = Log.getInstance(HasSequence.class);

    public static final String HAS_SEQUENCE_TAG = "XW";

    @Override
    public int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        IOUtil.assertFileIsReadable(SEQUENCES);

        final List<String> sequencesReadOne = readSequences();
        final List<String> sequencesReadTwo = new ArrayList<>();
        for (final String sequence : sequencesReadOne) {
            sequencesReadTwo.add(SequenceUtil.reverseComplement(sequence));
        }

        final SamReader reader = SamReaderFactory.makeDefault().open(INPUT);
        final SAMFileHeader header = reader.getFileHeader();
        if (header.getSortOrder() != SAMFileHeader.SortOrder.queryname) {
            throw new PicardException("Expects a queryname sorted input file, was: " + header.getSortOrder().name());
        }
        final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, OUTPUT);

        final Histogram<String> histogram = new Histogram<>();
        for (int i = 0; i < sequencesReadOne.size(); i++) {
            histogram.increment(sequencesReadOne.get(i), 0);
            histogram.increment(sequencesReadTwo.get(i), 0);
        }

        final SAMRecordIterator iterator = reader.iterator();
        final ProgressLogger progress = new ProgressLogger(log, 1000000);
        while (iterator.hasNext()) {
            final SAMRecord rec = iterator.next();
            final SAMRecord rec2 = (iterator.hasNext() && rec.getReadPairedFlag()) ? iterator.next() : null;

            if (null == rec2) { // ignore
                writer.addAlignment(rec);
                progress.record(rec);
            } else {
                annotateRecords(rec, rec2,
                        sequencesReadOne, sequencesReadTwo,
                        histogram,
                        MAX_MISMATCHES, MINIMUM_BASE_QUALITY,
                        NO_CALLS_ARE_BASES, START_ONLY);
                writer.addAlignment(rec);
                writer.addAlignment(rec2);
                progress.record(rec, rec2);
            }
        }

        CloserUtil.close(reader);
        writer.close();

        // finalize the histogram
        final Histogram<String> outputHistogram = new Histogram<>("sequence", "count");
        for (int i = 0; i < sequencesReadOne.size(); i++) {
            final String sequentReadOne = sequencesReadOne.get(i);
            outputHistogram.increment(sequentReadOne, histogram.get(sequentReadOne).getValue());
        }
        for (int i = 0; i < sequencesReadTwo.size(); i++) {
            final String sequentReadOne = sequencesReadOne.get(i);
            final String sequenceReadTwo = sequencesReadTwo.get(i);
            outputHistogram.increment(sequentReadOne, histogram.get(sequenceReadTwo).getValue());
        }

        // write metrics
        final MetricsFile<?, String> metricsFile = getMetricsFile();
        metricsFile.addHistogram(outputHistogram);
        metricsFile.write(METRICS_FILE);

        return 0;
    }

    private List<String> readSequences() {
        final List<String> sequences = new ArrayList<>();

        // read in the sequences
        try {
            final BufferedReader reader = new BufferedReader(new FileReader(SEQUENCES));
            String line;
            while ((line = reader.readLine()) != null) {
                sequences.add(line.trim());
            }
            CloserUtil.close(reader);

        } catch (final Exception ex) {
            throw new PicardException(ex.getMessage());
        }

        return sequences;
    }

    private static void annotateRecords(final SAMRecord rec,
                                        final SAMRecord rec2,
                                        final List<String> readOneSequences,
                                        final List<String> readTwoSequences,
                                        final Histogram<String> histogram,
                                        final int maxMismatches,
                                        final int minimumBaseQuality,
                                        final boolean noCallsAreBases,
                                        final boolean startOnly) {
        final boolean matchesReadOne = matchesSequence(rec, readOneSequences, histogram, maxMismatches, minimumBaseQuality, noCallsAreBases, startOnly);
        final boolean matchesReadTwo = matchesSequence(rec2, readTwoSequences, histogram, maxMismatches, minimumBaseQuality, noCallsAreBases, startOnly);

        final int value;
        if (matchesReadOne) {
            if (matchesReadTwo) value = 3;
            else value = 1;
        } else if (matchesReadTwo) {
            value = 2;
        } else {
            value = 0;
        }

        rec.setAttribute(HAS_SEQUENCE_TAG, value);
        rec2.setAttribute(HAS_SEQUENCE_TAG, value);
    }

    private static boolean matchesSequence(final SAMRecord rec,
                                           final List<String> sequences,
                                           final Histogram<String> histogram,
                                           final int maxMismatches,
                                           final int minimumBaseQuality,
                                           final boolean noCallsAreBases,
                                           final boolean startOnly) {
        if (startOnly) {
            for (final String sequence : sequences) {
                final int numMismatches = countMismatches(
                        rec.getReadBases(),
                        rec.getBaseQualities(),
                        sequence.getBytes(),
                        minimumBaseQuality,
                        noCallsAreBases,
                        0
                );
                if (numMismatches <= maxMismatches) {
                    histogram.increment(sequence);
                    return true;
                }
            }
        }
        else {
            for (final String sequence : sequences) {
                for (int i = 0; i < rec.getReadLength() - sequence.length(); i++) {
                    final int numMismatches = countMismatches(
                            rec.getReadBases(),
                            rec.getBaseQualities(),
                            sequence.getBytes(),
                            minimumBaseQuality,
                            noCallsAreBases,
                            i
                    );
                    if (numMismatches <= maxMismatches) {
                        histogram.increment(sequence);
                        return true;
                    }
                }
            }
        }
        return false;
    }

    public static int countMismatches(final byte[] observedBases,
                                      final byte[] observedQualities,
                                      final byte[] expectedBases,
                                      final int minimumBaseQuality,
                                      final boolean noCallsAreBases,
                                      final int observedOffset) {
        int mismatches = 0;
        int observedBaseNumber = observedOffset;
        int expectedBaseNumber = 0;
        while (observedBaseNumber < observedBases.length && expectedBaseNumber < expectedBases.length) {
            if (noCallsAreBases || !SequenceUtil.isNoCall(expectedBases[expectedBaseNumber])) {
                if (!SequenceUtil.basesEqual(observedBases[observedBaseNumber], expectedBases[expectedBaseNumber])) ++mismatches;
                else if (observedQualities[observedBaseNumber] < minimumBaseQuality) ++mismatches;
            }
            observedBaseNumber++;
            expectedBaseNumber++;
        }
        return mismatches;
    }
}

