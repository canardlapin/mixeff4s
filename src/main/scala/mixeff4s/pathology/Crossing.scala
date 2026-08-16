package mixeff4s.pathology

/** Secondary intercept-only grouping factor. Observations come from `cells`, or the full Cartesian product. */
final case class CrossedSpec(
    name: String,
    nLevels: Int,
    reVar: Double,
    cells: Option[Vector[(Int, Int)]] = None
)

/** Bipartite cell-graph summary. Orphans are levels that appear in no cell. */
final case class CrossedSummary(
    nPrimary: Int,
    nSecondary: Int,
    nCells: Int,
    nComponents: Int,
    primaryOrphans: Vector[Int],
    secondaryOrphans: Vector[Int]
)

object Crossing:
  def summarise(nPrimary: Int, nSecondary: Int, cells: Vector[(Int, Int)]): CrossedSummary =
    val primaryPresent = Array.fill(nPrimary)(false)
    val secondaryPresent = Array.fill(nSecondary)(false)
    cells.foreach: (i, j) =>
      if i >= 0 && i < nPrimary then primaryPresent(i) = true
      if j >= 0 && j < nSecondary then secondaryPresent(j) = true
    val total = nPrimary + nSecondary
    val parent = Array.tabulate(total)(identity)
    def find(x: Int): Int =
      var root = x
      while parent(root) != root do root = parent(root)
      var cur = x
      while parent(cur) != root do
        val next = parent(cur)
        parent(cur) = root
        cur = next
      root
    cells.foreach: (i, j) =>
      if i >= 0 && i < nPrimary && j >= 0 && j < nSecondary then
        val a = find(i)
        val b = find(nPrimary + j)
        if a != b then parent(a) = b
    val roots = scala.collection.mutable.TreeSet.empty[Int]
    primaryPresent.zipWithIndex.foreach: (present, i) =>
      if present then roots += find(i)
    secondaryPresent.zipWithIndex.foreach: (present, j) =>
      if present then roots += find(nPrimary + j)
    CrossedSummary(
      nPrimary,
      nSecondary,
      cells.length,
      roots.size,
      (0 until nPrimary).filterNot(primaryPresent).toVector,
      (0 until nSecondary).filterNot(secondaryPresent).toVector
    )

  def fullCross(name: String, nLevels: Int, reVar: Double): CrossedSpec =
    CrossedSpec(name, nLevels, reVar, cells = None)

  def blockDiagonal(name: String, blockSize: Int, nBlocks: Int, reVar: Double): (Vector[Int], CrossedSpec) =
    val nLevels = blockSize * nBlocks
    val cells = (0 until nBlocks).flatMap: b =>
      val start = b * blockSize
      (0 until blockSize).flatMap: i =>
        (0 until blockSize).map(j => (start + i, start + j))
    (Vector.fill(nLevels)(1), CrossedSpec(name, nLevels, reVar, Some(cells.toVector)))
