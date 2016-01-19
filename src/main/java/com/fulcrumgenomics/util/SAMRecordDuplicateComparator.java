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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import java.util.List;

/**
 * Compares records based on if they should be considered PCR Duplicates (see MarkDuplicates), being barcode aware.
 */
public class SAMRecordDuplicateComparator extends htsjdk.samtools.SAMRecordDuplicateComparator {

    /** An enum to provide type-safe keys for transient attributes the comparator puts on SAMRecords. */
    private enum Attr {
        BarcodeTag, ReadOneBarcodeTag, ReadTwoBarcodeTag
    }

    private final String barcodeTag;
    private final String readOneBarcodeTag;
    private final String readTwoBarcodeTag;
    private final String[] barcodeTags;

    public SAMRecordDuplicateComparator(final List<SAMFileHeader> headers,
                                        final String barcodeTag,
                                        final String readOneBarcodeTag,
                                        final String readTwoBarcodeTag) {
        super(headers);
        this.barcodeTag        = barcodeTag;
        this.readOneBarcodeTag = readOneBarcodeTag;
        this.readTwoBarcodeTag = readTwoBarcodeTag;
        this.barcodeTags       = new String[] { barcodeTag, readOneBarcodeTag, readTwoBarcodeTag };
    }

    /**
     * See {@link  htsjdk.samtools.SAMRecordDuplicateComparator#compare(Object, Object)}.
     */
    @Override
    public int compare(final SAMRecord samRecord1, final SAMRecord samRecord2) {
        int cmp = compareBarcodes(samRecord1, samRecord2);
        if (cmp == 0) cmp = super.compare(samRecord1, samRecord2);
        return cmp;
    }

    /**
     * See {@link  htsjdk.samtools.SAMRecordDuplicateComparator#duplicateSetCompare(SAMRecord, SAMRecord)}.
     */
    @Override
    public int duplicateSetCompare(final SAMRecord samRecord1, final SAMRecord samRecord2) {
        int cmp = compareBarcodes(samRecord1, samRecord2);
        if (cmp == 0) cmp = super.duplicateSetCompare(samRecord1, samRecord2);
        return cmp;    }

    /**
     * See {@link  htsjdk.samtools.SAMRecordDuplicateComparator#fileOrderCompare(SAMRecord, SAMRecord)}.
     */
    @Override
    public int fileOrderCompare(final SAMRecord samRecord1, final SAMRecord samRecord2) {
        int cmp = compareBarcodes(samRecord1, samRecord2);
        if (cmp == 0) cmp = super.fileOrderCompare(samRecord1, samRecord2);
        return cmp;
    }

    /**
     * Populates the set of transient attributes on SAMRecords if they are not already there.
     */
    private void populateTransientAttributes(final SAMRecord... recs) {
        for (final SAMRecord rec : recs) {
            if (barcodeTag != null && rec.getTransientAttribute(Attr.BarcodeTag) != null) continue;
            if (this.barcodeTag != null)        rec.setTransientAttribute(Attr.BarcodeTag,        rec.getStringAttribute(this.barcodeTag));
            if (this.readOneBarcodeTag != null) rec.setTransientAttribute(Attr.ReadOneBarcodeTag, rec.getStringAttribute(this.readOneBarcodeTag));
            if (this.readTwoBarcodeTag != null) rec.setTransientAttribute(Attr.ReadTwoBarcodeTag, rec.getStringAttribute(this.readTwoBarcodeTag));
        }
    }

    /**
     * Compares the values in the barcode tags.  This will also populate transient tags.
     */
    private int compareBarcodes(final SAMRecord samRecord1, final SAMRecord samRecord2) {
        populateTransientAttributes(samRecord1, samRecord2);
        for (final String barcodeTag : barcodeTags) {
            if (barcodeTag != null) {
                final int cmp = samRecord1.getStringAttribute(barcodeTag).compareTo(samRecord2.getStringAttribute(barcodeTag));
                if (cmp != 0) return cmp;
            }
        }
        return 0;
    }
}
