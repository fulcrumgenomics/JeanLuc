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

import picard.util.IlluminaUtil;

import java.util.List;

/**
 * A little class to store the sample barcode as bytes and concatenated with `IlluminaUtil.barcodeSeqsToString`.
 */
public class SampleBarcode {
    public final byte[] barcodeBytes;
    public final String concatenatedBarcode;

    public SampleBarcode(final List<String> sampleBarcodes) {
        final int numBases = sampleBarcodes.stream().mapToInt(String::length).sum();
        this.barcodeBytes = new byte[numBases];
        int barcodeBytesSrcOffset = 0;
        for (final String sampleBarcode : sampleBarcodes) {
            System.arraycopy(sampleBarcode.getBytes(), 0, barcodeBytes, barcodeBytesSrcOffset, sampleBarcode.length());
            barcodeBytesSrcOffset += sampleBarcode.length();
        }
        this.concatenatedBarcode = IlluminaUtil.barcodeSeqsToString(sampleBarcodes);
    }

    @Override
    public boolean equals(final Object o) {
        return o instanceof SampleBarcode && this.hashCode() == o.hashCode();
    }

    @Override
    public int hashCode() { return concatenatedBarcode.hashCode(); }

    @Override
    public String toString() { return concatenatedBarcode; }
}