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

package com.fulcrumgenomics.util.miseq;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public class SampleSheetTest {

    private final static File TEST_DIR = new File("testdata/com/fulcrumgenomics/util/miseq/");

    @Test
    public void testSimpleSampleSheet() {
        final SampleSheet sampleSheet = new SampleSheet(new File(TEST_DIR, "SampleSheet.csv"));
        Assert.assertEquals(sampleSheet.size(), 12);

        for (final Sample sample : sampleSheet) {
            Assert.assertEquals(sample.sampleId,             "20000101-EXPID-" + sample.sampleOrdinal);
            Assert.assertEquals(sample.sampleName,           "Sample_Name_"    + sample.sampleOrdinal);
            Assert.assertEquals(sample.libraryId,            "Sample_Name_"    + sample.sampleOrdinal);
            Assert.assertEquals(sample.project,              "Sample_Project_" + sample.sampleOrdinal);
            Assert.assertEquals(sample.description,          "Description_"    + sample.sampleOrdinal);
            Assert.assertEquals(sample.r1SampleBarcodeBases, "GATTACAG".getBytes());
            Assert.assertEquals(sample.r2SampleBarcodeBases, "GATTACAGA".getBytes());
            Assert.assertEquals(sample.i7IndexBases,         "GATTACAACGT".getBytes());
            Assert.assertEquals(sample.i5IndexBases,         "GATTACA".getBytes());
        }
    }

    @Test
    public void testOnlyRequiredSampleSheet() {
        final SampleSheet sampleSheet = new SampleSheet(new File(TEST_DIR, "SampleSheetOnlyRequired.csv"));
        Assert.assertEquals(sampleSheet.size(), 12);

        for (final Sample sample : sampleSheet) {
            Assert.assertEquals(sample.sampleId,             "20000101-EXPID-" + sample.sampleOrdinal);
            Assert.assertEquals(sample.sampleName,           "Sample_Name_"    + sample.sampleOrdinal);
            Assert.assertEquals(sample.libraryId,            "Sample_Name_"    + sample.sampleOrdinal);
            Assert.assertEquals(sample.project,              "Sample_Project_" + sample.sampleOrdinal);
            Assert.assertEquals(sample.description,          "Description_"    + sample.sampleOrdinal);
            Assert.assertEquals(sample.r1SampleBarcodeBases, null);
            Assert.assertEquals(sample.r2SampleBarcodeBases, null);
            Assert.assertEquals(sample.i7IndexBases,         null);
            Assert.assertEquals(sample.i5IndexBases,         null);
        }
    }
}
