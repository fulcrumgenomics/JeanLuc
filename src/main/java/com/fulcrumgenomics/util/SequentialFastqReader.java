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


import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;

import java.io.File;
import java.util.*;

/** Reads multiple FASTQs one after the other. This is useful for when we have more than one FASTQ for a set of reads. */
public class SequentialFastqReader extends FastqReader {
    private final List<FastqReader> readers = new ArrayList<>();
    private int readerIndex = 0;

    public SequentialFastqReader(final List<File> inputs, final boolean skipBlankLines) {
        super(inputs.get(0), skipBlankLines); // ignore anything in FastqReader
        for (final File input : inputs) {
            readers.add(new FastqReader(input, skipBlankLines));
        }
    }

    public SequentialFastqReader(final List<File> inputs) {
        this(inputs, false);

    }

    @Override
    public boolean hasNext() {
        while (readerIndex < readers.size()) {
            if (readers.get(readerIndex).hasNext()) return true;
            readerIndex++;
        }
        return false;
    }

    @Override
    public FastqRecord next() {
        if (readers.size() <= readerIndex) {
            throw new NoSuchElementException();
        }
        return readers.get(readerIndex).next();
    }

    @Override
    public void close() {
        readers.stream().forEach(FastqReader::close);
    }

    @Override
    public int getLineNumber() { return readers.get(readerIndex).getLineNumber(); }

    @Override
    public File getFile() { return readers.get(readerIndex).getFile(); }

    @Override
    public String toString() {
        return "FastqReader["+(this.getFile() == null?"":this.getFile())+ " Line:"+getLineNumber()+"]";
    }
}

