package mixeff4s.pathology

import gale.linalg.Matrix
import gale.spectral.{Eigen, EigenSelection, EigenVectors}
import mixeff4s.error.{FitResult, MixedModelError}

/** Symmetric eigendecomposition and PSD square root for truth covariances of any q. */
object SymmetricPsd:
  private val PsdTol = 1e-10

  def eigvals(cov: Vector[Vector[Double]]): FitResult[Vector[Double]] =
    decompose(cov, wantVectors = false).map(_._1)

  def sqrt(cov: Vector[Vector[Double]]): FitResult[Vector[Vector[Double]]] =
    decompose(cov, wantVectors = true).flatMap: (vals, vecs) =>
      if vals.exists(_ < -PsdTol) then Left(MixedModelError.InvalidArgument("symmetric matrix is not PSD"))
      else
        val scales = vals.map(v => math.sqrt(math.max(v, 0.0)))
        val q = vals.length
        Right(
          Vector.tabulate(q, q): (i, j) =>
            var acc = 0.0
            var k = 0
            while k < q do
              acc += vecs(i)(k) * scales(k) * vecs(j)(k)
              k += 1
            acc
        )

  private def decompose(
      cov: Vector[Vector[Double]],
      wantVectors: Boolean
  ): FitResult[(Vector[Double], Vector[Vector[Double]])] =
    val q = cov.length
    if q == 0 then Right((Vector.empty, Vector.empty))
    else if cov.exists(_.length != q) then Left(MixedModelError.InvalidArgument("symmetric matrix is not square"))
    else
      val a = Matrix.tabulate(q, q)((i, j) => 0.5 * (cov(i)(j) + cov(j)(i)))
      val vectors = if wantVectors then EigenVectors.Right else EigenVectors.ValuesOnly
      Eigen.eigSymmetric(a, EigenSelection.All, vectors) match
        case Left(err) =>
          Left(MixedModelError.InvalidArgument(s"symmetric eigen failed: $err"))
        case Right(eig) =>
          val vals = Vector.tabulate(eig.size)(i => eig.eigenvalues(i))
          val vecs =
            if !wantVectors || eig.eigenvectors.cols == 0 then Vector.empty
            else Vector.tabulate(q, q)((i, j) => eig.eigenvectors(i, j))
          Right((vals, vecs))
