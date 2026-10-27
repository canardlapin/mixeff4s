package mixeff4s.design

/** Column identity for rank-truncated fixed effects.
  *
  * Names refer to expanded design columns, not formula terms. Retained indices are zero-based positions in the
  * requested design, in the order used by beta and both covariance axes. Dropped columns are reported in requested
  * order. Their coefficients are not separately identified by the retained solution.
  */
final case class FixedEffectLayout(requestedNames: Vector[String], retainedIndices: Vector[Int]):
  require(retainedIndices.distinct.length == retainedIndices.length, "retained column indices must be distinct")
  require(retainedIndices.forall(i => i >= 0 && i < requestedNames.length), "retained column index out of bounds")

  val retainedNames: Vector[String] = retainedIndices.map(requestedNames)
  val droppedIndices: Vector[Int] =
    val retained = retainedIndices.toSet
    requestedNames.indices.filterNot(retained.contains).toVector
  val droppedNames: Vector[String] = droppedIndices.map(requestedNames)
  def rank: Int = retainedIndices.length
