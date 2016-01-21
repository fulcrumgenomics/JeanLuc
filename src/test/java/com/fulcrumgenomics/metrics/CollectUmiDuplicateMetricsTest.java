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
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class CollectUmiDuplicateMetricsTest {

    @Test(dataProvider = "testCompareBarcodesDataProvider")
    public void testCompareBarcodes(final String barcodeTagString1,
                                    final String barcodeTagString2,
                                    final int maxMismatches,
                                    final int expectedValue) {
        final int actualValue = CollectUmiDuplicateMetrics.compareBarcodes(
                barcodeTagString1,
                barcodeTagString2,
                maxMismatches);
        Assert.assertEquals(actualValue, expectedValue);
    }

    @DataProvider(name = "testCompareBarcodesDataProvider")
    public Object[][] testCompareBarcodesDataProvider() {
        return new Object[][] {
                // perfect match
                {"GATTACA", "GATTACA", 0, 0},
                // one mismatch, but max mismatches is 1
                {"AATTACA", "GATTACA", 0, -1},
                // one mismatch, but max mismatches is 1
                {"GATTACA", "AATTACA", 0, 1},
                // one mismatch, with max mismatches set to 1
                {"AATTACA", "GATTACA", 1, 0},
                // one mismatch, but max mismatches is 1
                {"ATTTACA", "GATTACA", 1, -1},
                // one mismatch, with max mismatches set to 2
                {"ATTTACA", "GATTACA", 2, 0},
                // perfect match with a delimeter
                {"GATTACA-GATTACA", "GATTACA-GATTACA", 0, 0},
                // perfect match with a delimeter and N-base
                {"GATTACA-GATTACA", "GATTACANGATTACA", 0, 0},
        };
    }
}
