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
package com.fulcrumgenomics.util.miseq

/**
  * Represents information about a sample within an Illumina Experiment Manager sample sheet.
  *
  * Optionally contains information about dual indexes: i7 and i5.
  *
  * @author Nils Homer
  *
  * @param sampleOrdinal the sample ordinal if this sample belongs to a sample set.
  * @param sampleId the unique sample identifier.
  * @param sampleName the sample name.
  * @param libraryId the library identifier.
  * @param project the project identifier.
  * @param description the sample description.
  * @param r1SampleBarcodeBases the inline sample barcode bases in the first read of a pair.
  * @param r2SampleBarcodeBases the inline sample barcode bases in the second read of a pair.
  * @param i7IndexBases the sample barcode bases in the i7 read.
  * @param i5IndexBases the sample barcode bases in the i5 read.
  * @param sampleBarcode the full sample barcode for optimized access.
  */
class Sample(val sampleOrdinal: Int,
             val sampleId: String,
             val sampleName: String,
             val libraryId: Option[String],
             val project: Option[String],
             val description: Option[String],
             val r1SampleBarcodeBases: Option[Array[Byte]] = None,
             val r2SampleBarcodeBases: Option[Array[Byte]] = None,
             val i7IndexBases: Option[Array[Byte]] = None,
             val i5IndexBases: Option[Array[Byte]] = None,
             var sampleBarcode: Option[SampleBarcode] = None) {
}

