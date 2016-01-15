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

/**
 * Represents information about a sample within an Illumina Experiment Manager sample sheet.
 *
 * Optionally contains information about dual indexes: i7 and i5.
 *
 * @author Nils Homer
 */
public class Sample {
    public final int sampleOrdinal;

    /**
     * Required fields
     */
    public final String sampleId;
    public final String sampleName;
    public final String libraryId;
    public final String project;
    public final String description;

    /**
     * Optional fields
     */
    public final byte[] r1SampleBarcodeBases; /** the inline sample barcode bases in the first read of a pair */
    public final byte[] r2SampleBarcodeBases; /** the inline sample barcode bases in the second read of a pair */
    public final byte[] i7IndexBases; /** the sample barcode bases in the i7 read */
    public final byte[] i5IndexBases; /** the sample barcode bases in the i7 read */

    /**
     * Non-standard fields
     */
    public SampleBarcode sampleBarcode;

    public Sample(final int sampleOrdinal,
                  final String sampleId,
                  final String sampleName,
                  final String libraryId,
                  final String project,
                  final String description) {
        this(sampleOrdinal, sampleId, sampleName, libraryId, project, description, null, null, null, null);
    }

    public Sample(final int sampleOrdinal,
                  final String sampleId,
                  final String sampleName,
                  final String libraryId,
                  final String project,
                  final String description,
                  final String r1SampleBarcodeBases,
                  final String r2SampleBarcodeBases,
                  final String i7IndexBases,
                  final String i5IndexBases) {
        this.sampleOrdinal = sampleOrdinal;
        this.sampleId = sampleId;
        this.sampleName = sampleName;
        this.libraryId = libraryId;
        this.project = project;
        this.description = description;
        this.r1SampleBarcodeBases = (r1SampleBarcodeBases == null) ? null : r1SampleBarcodeBases.getBytes();
        this.r2SampleBarcodeBases = (r2SampleBarcodeBases == null) ? null : r2SampleBarcodeBases.getBytes();
        this.i7IndexBases = (i7IndexBases == null) ? null : i7IndexBases.getBytes();
        this.i5IndexBases = (i5IndexBases == null) ? null : i5IndexBases.getBytes();
    }
}

