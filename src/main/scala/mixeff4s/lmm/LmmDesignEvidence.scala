package mixeff4s.lmm

import mixeff4s.design.CompiledDesign

/** Immutable row-aligned design evidence; no dense n by random-coefficient matrix is retained. */
final class LmmDesignEvidence private (
    val fixedColumns: Vector[Vector[Double]],
    val randomBlocks: Vector[RandomDesignEvidence],
    private val response: Vector[Double],
    private val criterion: Criterion,
    private val theta: Vector[Double],
    private val objective: Double,
    private val parmap: Vector[(Int, Int, Int)]
):
  private[mixeff4s] def matches(fit: LmmFit): Boolean =
    fit.n == response.length && fit.y == response && fit.options.criterion == criterion &&
      fit.theta == theta && fit.objective == objective && fit.parmap == parmap &&
      fit.beta.length == fixedColumns.length && fit.validatedContext.isRight

/** A grouping partition and its within-group basis, with either diagonal or full covariance. */
final class RandomDesignEvidence private[lmm] (
    val partition: Vector[Int],
    val columns: Vector[Vector[Double]],
    val diagonal: Boolean
)

object LmmDesignEvidence:
  private[lmm] def capture(
      design: CompiledDesign,
      criterion: Criterion,
      theta: Vector[Double],
      objective: Double,
      response: Vector[Double]
  ): LmmDesignEvidence =
    val blocks = design.reterms.map: term =>
      val labels = scala.collection.mutable.Map.empty[Int, Int]
      val partition = term.refs.map(ref => labels.getOrElseUpdate(ref, labels.size))
      new RandomDesignEvidence(
        partition,
        Vector.tabulate(term.vsize)(j => Vector.tabulate(design.n)(i => term.z(j, i))),
        term.nTheta == term.vsize
      )
    new LmmDesignEvidence(
      Vector.tabulate(design.p)(j => Vector.tabulate(design.n)(i => design.xy.xy(i, j))),
      blocks,
      response,
      criterion,
      theta,
      objective,
      design.parmap
    )
