package mixeff4s.stats

import mixeff4s.design.CompiledDesign

/** Exact row-aligned witnesses. Unknown equivalent reparameterizations are declined. */
private[stats] object LikelihoodNesting:
  final case class Block(partition: Vector[Int], columns: Vector[Vector[Double]], diagonal: Boolean)
  final case class Evidence(fixed: Vector[Vector[Double]], random: Vector[Block])
  def capture(design: CompiledDesign): Evidence =
    Evidence(
      Vector.tabulate(design.p)(j => Vector.tabulate(design.n)(i => design.xy.xy(i, j))),
      design.reterms.map: term =>
        val labels = scala.collection.mutable.Map.empty[Int, Int]
        Block(
          term.refs.map(r => labels.getOrElseUpdate(r, labels.size)),
          Vector.tabulate(term.vsize)(j => Vector.tabulate(design.n)(i => term.z(j, i))),
          term.nTheta == term.vsize
        )
    )
  def nested(a: Evidence, b: Evidence): Boolean =
    includes(a.fixed, b.fixed)(_ == _) && includes(a.random, b.random)(covered)
  def sameRandom(a: Evidence, b: Evidence): Boolean =
    a.random.length == b.random.length && includes(a.random, b.random): (x, y) =>
      x.diagonal == y.diagonal && x.columns.length == y.columns.length && covered(x, y)
  private def covered(a: Block, b: Block): Boolean =
    a.partition == b.partition && (a.diagonal || !b.diagonal) && includes(a.columns, b.columns)(_ == _)
  private def includes[A](a: Vector[A], b: Vector[A])(covers: (A, A) => Boolean): Boolean =
    if a.length > b.length then false
    else
      val owner = Array.fill(b.length)(-1)
      def assign(i: Int, seen: Array[Boolean]): Boolean = b.indices.exists: j =>
        if seen(j) || !covers(a(i), b(j)) then false
        else
          seen(j) = true
          if owner(j) < 0 || assign(owner(j), seen) then { owner(j) = i; true }
          else false
      a.indices.forall(i => assign(i, Array.fill(b.length)(false)))
