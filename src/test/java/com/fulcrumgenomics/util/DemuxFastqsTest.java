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


package com.fulcrumgenomics.util;

import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;
import picard.illumina.ExtractIlluminaBarcodes;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;

/** Tests DemuxFastqs */
public class DemuxFastqsTest {

    private static File createTestDirectory() {
        try {
            final File directory = Files.createTempDirectory("DemuxFastqsTest").toFile();
            directory.deleteOnExit();
            return directory;
        } catch (final Exception ex) {
            throw new PicardException(ex.getMessage());
        }
    }

    private final static File TEST_DIR = new File("testdata/com/fulcrumgenomics/util/DemuxFastqsTest/");

    @Test(dataProvider = "testDemuxFastqsTestDataProvider")
    public void testDemuxFastqsTest(final File     readOneFastq,
                                    final File     readTwoFastq,
                                    final File     i7Fastq,
                                    final File     i5Fastq,
                                    final String   readOneReadStructure,
                                    final String   readTwoReadStructure,
                                    final String   i7ReadStructure,
                                    final String   i5ReadStructure,
                                    final File     sampleSheet,
                                    final int      maxMismatches,
                                    final String[] barcodes,
                                    final int[]    numReadsPerBarcode) throws IOException {
        final File output = createTestDirectory();
        final File metricsFile = new File(output, "DemuxFastqsTest.metrics.txt");

        final DemuxFastqs program = new DemuxFastqs();
        program.READ_ONE_FASTQ          = Arrays.asList(readOneFastq);
        program.READ_TWO_FASTQ          = Arrays.asList(readTwoFastq);
        program.I7_FASTQ                = Arrays.asList(i7Fastq);
        program.I5_FASTQ                = Arrays.asList(i5Fastq);
        program.READ_ONE_READ_STRUCTURE = readOneReadStructure;
        program.READ_TWO_READ_STRUCTURE = readTwoReadStructure;
        program.I7_READ_STRUCTURE       = i7ReadStructure;
        program.I5_READ_STRUCTURE       = i5ReadStructure;
        program.SAMPLE_SHEET            = sampleSheet;
        program.MAX_MISMATCHES          = maxMismatches;
        program.OUTPUT                  = output;
        program.METRICS_FILE            = metricsFile;
        program.MIN_MISMATCH_DELTA      = 1; // for testing

        Assert.assertEquals(program.doWork(), 0);

        // check that all samples are represented in the metrics file
        final List<ExtractIlluminaBarcodes.BarcodeMetric> metrics = MetricsFile.readBeans(metricsFile);
        Assert.assertEquals(metrics.size(), barcodes.length);
        for (int i = 0; i < barcodes.length; i++) {
            final ExtractIlluminaBarcodes.BarcodeMetric metric = metrics.get(i);
            final String barcode = barcodes[i];
            final int numReads = numReadsPerBarcode[i];
            Assert.assertEquals(metric.BARCODE, barcode);
            final String sampleName = (i == (barcodes.length-1)) ? DemuxFastqs.UNMATCHED_SAMPLE_ID : ("Sample_Name_"+(i+1));
            Assert.assertEquals(metric.BARCODE_NAME, sampleName);
            Assert.assertEquals(metric.LIBRARY_NAME, sampleName);
            if (metric.READS != numReads) System.err.println(String.format("ERROR: %s %d != %d\n", sampleName, metric.READS, numReads));
            Assert.assertEquals(metric.READS, numReads);
        }

        IOUtil.deleteDirectoryTree(output);
    }

    @DataProvider(name = "testDemuxFastqsTestDataProvider")
    public Object[][] testDemuxFastqsTestDataProvider() {
        final File sampleSheet = new File(TEST_DIR + "/SampleSheet.csv");
        return new Object[][]{
                // molecular id on one index read, sample barcode on the other
                {getFastq(0), getFastq(1), getFastq(2), getFastq(3), "151T", "151T", "7B", "7M", sampleSheet, 1,
                        new String[] {"GATTACA", "GATTACC", "NNNNNNN"}, new int[] {1, 1, 1}},
                // both molecular id and sample barcode present on both index reads
                {getFastq(0), getFastq(1), getFastq(2), getFastq(3), "151T", "151T", "6M1B", "1B6M", sampleSheet, 1,
                        new String[] {"A-A", "C-C", "N-N"}, new int[] {1, 1, 1}},
                // molecular id on one index read, sample barcode on the other, tests matching with at most one mismatch (third read is two mismatches away from sample two, so no match)
                {getFastq(0), getFastq(1), getFastq(2), getFastq(3), "151T", "151T", "7M", "7B", sampleSheet, 1,
                        new String[] {"ACATTAG", "CGGGGGG", "NNNNNNN"}, new int[] {2, 0, 1}},
                // molecular id on one index read, sample barcode on the other, tests matching with two mismatches (now the third read matches sample two)
                {getFastq(0), getFastq(1), getFastq(2), getFastq(3), "151T", "151T", "7M", "7B", sampleSheet, 2,
                        new String[] {"ACATTAG", "CGGGGGG", "NNNNNNN"}, new int[] {2, 1, 0}}
        };
    }

    private File getFastq(final int index) {
        if (index == 0)      return new File(TEST_DIR, "/whole_run_S0_L001_R1_001.fastq");
        else if (index == 1) return new File(TEST_DIR, "/whole_run_S0_L001_R2_001.fastq");
        else if (index == 2) return new File(TEST_DIR, "/whole_run_S0_L001_I1_001.fastq");
        else if (index == 3) return new File(TEST_DIR, "/whole_run_S0_L001_I2_001.fastq");
        else throw new PicardException("Index out of range: " + index);
    }

    @Test(dataProvider = "testCountMismatchesDataProvider")
    public void testCountMismatches(final byte[] observedBases,
                                    byte[] observedQualities,
                                    final byte[] expectedBases,
                                    final int minimumBaseQuality,
                                    final int numExpectedMismatches) {
        final int numActualMismatches = DemuxFastqs.countMismatches(observedBases,
                observedQualities,
                expectedBases,
                minimumBaseQuality);
        Assert.assertEquals(numActualMismatches, numExpectedMismatches);
    }

    @DataProvider(name = "testCountMismatchesDataProvider")
    public Object[][] testCountMismatchesDataProvider() {
        return new Object[][] {
                // no mismatches
                {"GATTACA".getBytes(), SAMUtils.fastqToPhred("@@@@@@@"), "GATTACA".getBytes(), 0, 0},
                // two mismatches
                {"GATTACA".getBytes(), SAMUtils.fastqToPhred("@@@@@@@"), "GACCACA".getBytes(), 0, 2},
                // no calls
                {"GATTACA".getBytes(), SAMUtils.fastqToPhred("@@@@@@@"), "GANNACA".getBytes(), 0, 0},
                // all mismatches
                {"GATTACA".getBytes(), SAMUtils.fastqToPhred("@@@@@@@"), "CTAATGT".getBytes(), 0, 7},
                // base qualities are too low for comparison
                {"GATTACA".getBytes(), SAMUtils.fastqToPhred("@@@@@@@"), "GATTACA".getBytes(), 32, 7},
        };
    }
}
