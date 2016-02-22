package com.fulcrumgenomics.util

import scala.collection.JavaConversions._

import com.fulcrumgenomics.cmdline.{ClpGroups, JeanLucTool}
import dagr.commons.CommonsDef.{DirPath, FilePath}
import dagr.sopt.{arg, clp}
import picard.util.IlluminaUtil.IlluminaAdapterPair

/**
  * Program for picking sets of indices of arbitrary length that meet certain constraints
  * and attempt to maximize the edit distance between all members of the set picked.
  *
  * @author Tim Fennell
  */
@clp(
 description = "Picks a set of molecular indices that should work well together.",
 group = ClpGroups.Utilities)
class PickIlluminaIndices
(
  @arg(flag="l", doc="The length of each barcode sequence.")                           val length: Int = 8,
  @arg(flag="n", doc="The number of indices desired.")                                 val indices: Int,
  @arg(flag="e", doc="The minimum edit distance between two indices in the set.")      val editDistance: Int = 3,
  @arg(flag="o", doc="File to write indices to.")                                      val output: FilePath,
  @arg(          doc="Allow indices that are lexical reverses of one another")         val allowReverses: Boolean = false,
  @arg(          doc="Allow indices that are reverse complements of one another")      val allowReverseComplements: Boolean = false,
  @arg(          doc="Allow indices that are palindromic (bases == rev(bases)).")      val allowPalindromes: Boolean = false,
  @arg(          doc="Reject indices with a homopolymer of greater than this length.") val maxHomopolymer: Int = 2,
  @arg(          doc="The minimum GC fraction for a barcode to be accepted.")          val minGc: Double = 0,
  @arg(          doc="The maximum GC fraction for a barcode to be accepted.")          val maxGc: Double = 0.7,
  @arg(flag="t", doc="Number of threads to use.")                                      val threads: Int = 4,
  @arg(          doc="The installation directory for ViennaRNA.")                      val viennaRnaDir: DirPath,
  @arg(          doc="The lowest acceptable secondary structure deltaG.")              val minDeltaG: Double = -10,
  @arg(          doc="The indexed adapter sequence into which the indices will be integrated.")
  val adapters: Seq[String] = Seq(IlluminaAdapterPair.DUAL_INDEXED.get5PrimeAdapter, IlluminaAdapterPair.DUAL_INDEXED.get3PrimeAdapter),
  @arg(          doc="Sequences that should be avoided.  Any kmer of 'length' that appears in these sequences and their " + "reverse complements will be thrown out.")
  val avoidSequence: Seq[String] = IlluminaAdapterPair.values().flatMap(p => Seq(p.get5PrimeAdapter(), p.get3PrimeAdapter()))
) extends JeanLucTool{

  override def execute(): Unit = {
    val cmd = new PickIlluminaIndicesCommand()
    cmd.INDEX_LENGTH              = length
    cmd.N_INDICES                 = indices
    cmd.EDIT_DISTANCE             = editDistance
    cmd.ALLOW_REVERSES            = allowReverses
    cmd.ALLOW_REVERSE_COMPLEMENTS = allowReverseComplements
    cmd.ALLOW_PALINDROMES         = allowPalindromes
    cmd.MAX_HOMOPOLYMER           = maxHomopolymer
    cmd.MIN_GC                    = minGc
    cmd.MAX_GC                    = maxGc
    cmd.OUTPUT                    = output.toFile
    cmd.NUM_THREADS               = threads
    cmd.VIENNA_RNA_DIR            = viennaRnaDir.toFile
    cmd.MIN_DELTAG                = minDeltaG
    cmd.INDEX_ADAPTER             = adapters
    cmd.AVOID_SEQUENCE            = avoidSequence
    cmd.execute()
  }
}
