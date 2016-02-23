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
package com.fulcrumgenomics.util

import java.io.File
import java.nio.file.Path

import htsjdk.samtools.fastq.FastqReader
import htsjdk.samtools.fastq.FastqRecord

/** Reads multiple FASTQs one after the other. This is useful for when we have more than one FASTQ for a set of reads. */
class SequentialFastqReader(inputs: List[Path], skipBlankLines: Boolean) extends FastqReader(inputs.head.toFile, skipBlankLines) {
  private val readers: List[FastqReader] = inputs.map(input => new FastqReader(input.toFile, skipBlankLines))
  private val readersIterator: Iterator[FastqReader] = readers.toIterator
  private var reader: Option[FastqReader] = Some(readersIterator.next()) // should always have at least one reader (the passing of inputs.head.toFile to FastqReader above should have failed otherwise).

  def this(inputs: List[Path]) = this(inputs, false)

  override def hasNext: Boolean = {
    reader match {
      case Some(r) if r.hasNext => true
      case Some(r) if readersIterator.hasNext =>
        reader = Some(readersIterator.next())
        hasNext
      case Some(r) => false
      case None if readersIterator.hasNext =>
        reader = Some(readersIterator.next())
        hasNext
      case None => false
    }
  }

  override def next: FastqRecord = {
    if (!hasNext) throw new NoSuchElementException
    reader match {
      case Some(r) => r.next()
      case None => throw new NoSuchElementException
    }
  }

  override def close() = readers.foreach(_.close())

  override def getLineNumber: Int = reader match {
      case Some(r) => r.getLineNumber
      case None => throw new NoSuchElementException
    }
  override def getFile: File = reader match {
    case Some(r) => r.getFile
    case None => throw new NoSuchElementException
  }

  override def toString: String = "FastqReader[" + (if (this.getFile == null) "" else this.getFile) + " Line:" + getLineNumber + "]"
}

