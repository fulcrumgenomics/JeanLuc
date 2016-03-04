/*
 * The MIT License
 *
 * Copyright (c) 2016 Fulcrum Genomics LLC
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package com.fulcrumgenomics.cmdline

import com.fulcrumgenomics.cmdline.JeanLucMain.FailureException
import dagr.commons.util.LazyLogging
import dagr.sopt.cmdline.CommandLineParser

/**
  * Main program for JeanLuc that loads everything up and runs the appropriate sub-command
  */
object JeanLucMain {
  /** The main method */
  def main(args: Array[String]): Unit = new JeanLucMain().makeItSo(args)

  /**
    * Exception class intended to be used by [[JeanLucMain]] and [[JeanLucTool]] to communicate
    * non-exceptional(!) failures when running a tool.
    */
  private[cmdline] case class FailureException(exit:Int = 1, message:Option[String] = None) extends RuntimeException
}

class JeanLucMain extends LazyLogging {
  /** The main method */
  def makeItSo(args: Array[String]): Unit = {
    val parser = new CommandLineParser[JeanLucTool]("JeanLuc")

    parser.parseSubCommand(args=args, packageList=packageList) match {
      case None => System.exit(1)
      case Some(tool) =>
        try {
          tool.execute
          System.exit(0)
        }
        catch {
          case ex: FailureException =>
            ex.message.foreach(logger.fatal)
            System.exit(ex.exit)
        }
    }
  }

  /** The packages we wish to include in our command line **/
  protected def packageList: List[String] = List[String]("com.fulcrumgenomics")
}
