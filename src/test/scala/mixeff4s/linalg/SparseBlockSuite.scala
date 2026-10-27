package mixeff4s.linalg

class SparseBlockSuite extends munit.FunSuite:
  private def sparse(rows: Int, cols: Int, entries: (Int, Int)*): SparseWorkMat =
    SparseWorkMat.fromCoordinates(rows, cols, entries)

  test("structural zeros persist and clones isolate mutable values"):
    val mat = sparse(2, 2, (0, 1), (1, 0))
    assert(mat.hasEntry(0, 1))
    assertEquals(mat(0, 1), 0.0)
    intercept[IllegalArgumentException](mat.update(0, 0, 0.0))
    mat.update(0, 1, 2.0)
    val cloned = mat.cloneMat
    cloned.update(0, 1, 7.0)
    assertEquals(mat(0, 1), 2.0)
    assertEquals(cloned(0, 1), 7.0)

  test("sparse subtraction preflight prevents late partial mutation"):
    val target = sparse(2, 2, (0, 0))
    target.update(0, 0, 9.0)
    val left = sparse(2, 1, (0, 0), (1, 0))
    val right = sparse(2, 1, (0, 0), (1, 0))
    left.update(0, 0, 2.0); left.update(1, 0, 3.0)
    right.update(0, 0, 5.0); right.update(1, 0, 7.0)
    val promoted = BlockedCholesky.prepareSubtractTarget(
      MatrixBlock.Sparse(target),
      MatrixBlock.Sparse(left),
      MatrixBlock.Sparse(right)
    )
    assert(promoted.isInstanceOf[MatrixBlock.Dense])
    BlockedCholesky.subtractProduct(promoted, MatrixBlock.Sparse(left), MatrixBlock.Sparse(right))
    val dense = promoted.asDense
    assertEquals(dense(0, 0), -1.0)
    assertEquals(dense(0, 1), -14.0)
    assertEquals(dense(1, 0), -15.0)
    assertEquals(dense(1, 1), -21.0)
    intercept[IllegalArgumentException](
      BlockedCholesky.subtractProduct(MatrixBlock.Sparse(target), MatrixBlock.Sparse(left), MatrixBlock.Sparse(right))
    )
    assertEquals(target(0, 0), 9.0)

  test("CSC subtraction updates a dense target without sparse densification"):
    val left = sparse(2, 1, (0, 0), (1, 0))
    val right = sparse(2, 1, (0, 0), (1, 0))
    left.update(0, 0, 2.0); left.update(1, 0, 3.0)
    right.update(0, 0, 5.0); right.update(1, 0, 7.0)
    val target = WorkMat.zeros(2, 2)
    BlockedCholesky.subtractProduct(MatrixBlock.Dense(target), MatrixBlock.Sparse(left), MatrixBlock.Sparse(right))
    assertEquals(target(0, 0), -10.0)
    assertEquals(target(0, 1), -14.0)
    assertEquals(target(1, 0), -15.0)
    assertEquals(target(1, 1), -21.0)

  test("CSC rank-k uses structural support for a diagonal target"):
    val source = sparse(2, 2, (0, 0), (1, 1))
    source.update(0, 0, 2.0); source.update(1, 1, 3.0)
    val target = Array(10.0, 20.0)
    BlockedCholesky.rankKDowndate(MatrixBlock.Diagonal(target), MatrixBlock.Sparse(source))
    assertEquals(target.toSeq, Seq(6.0, 11.0))

  test("dense copy from sparse clears stale fill"):
    val src = sparse(2, 2, (0, 1))
    src.update(0, 1, 4.0)
    val dst = WorkMat.tabulate(2, 2)((_, _) => 8.0)
    MatrixBlock.copyBlock(MatrixBlock.Dense(dst), MatrixBlock.Sparse(src))
    assertEquals(dst(0, 0), 0.0)
    assertEquals(dst(0, 1), 4.0)
    assertEquals(dst(1, 0), 0.0)
    assertEquals(dst(1, 1), 0.0)

  test("sparse off-diagonal scaling visits stored support, not structural gaps"):
    val source = sparse(2, 2, (0, 1))
    source.update(0, 1, 6.0)
    val target = sparse(2, 2, (0, 1))
    val one = WorkMat.tabulate(1, 1)((_, _) => 2.0)
    BlockedCholesky.copyAndScaleOffdiag(MatrixBlock.Sparse(target), MatrixBlock.Sparse(source), one, one)
    assertEquals(target(0, 1), 24.0)
    assertEquals(target(1, 0), 0.0)

  test("sparse vector off-diagonal scaling transforms only occupied rectangles"):
    val source = sparse(4, 4, (0, 2))
    source.update(0, 2, 5.0)
    val rectangle = Seq((0, 2), (0, 3), (1, 2), (1, 3))
    val target = SparseWorkMat.fromCoordinates(4, 4, rectangle)
    val identity = WorkMat.zeros(2, 2)
    identity(0, 0) = 1.0; identity(1, 1) = 1.0
    BlockedCholesky.copyAndScaleOffdiag(MatrixBlock.Sparse(target), MatrixBlock.Sparse(source), identity, identity)
    assertEquals(target(0, 2), 5.0)
    assertEquals(target(1, 3), 0.0)

  test("structured matvec avoids dense block traversal semantics"):
    val b0 = WorkMat.tabulate(2, 2)((r, c) => if r == c then 2.0 else 1.0)
    val b1 = WorkMat.tabulate(2, 2)((r, c) => if r == c then 3.0 else 0.0)
    val dst = Array.fill(4)(0.0)
    MatrixBlock.BlockDiagonal(Array(b0, b1)).addMatVec(Array(1.0, 2.0, 3.0, 4.0), dst)
    assertEquals(dst.toSeq, Seq(4.0, 5.0, 9.0, 12.0))

  test("sparse right solve preserves diagonal and block-diagonal closures"):
    val diagonal = sparse(1, 2, (0, 0), (0, 1))
    diagonal.update(0, 0, 4.0); diagonal.update(0, 1, 9.0)
    BlockedCholesky.rdivLowerTranspose(MatrixBlock.Sparse(diagonal), MatrixBlock.Diagonal(Array(2.0, 3.0)))
    assertEquals(diagonal(0, 0), 2.0)
    assertEquals(diagonal(0, 1), 3.0)
    val block = WorkMat.zeros(2, 2)
    block(0, 0) = 2.0; block(1, 0) = 1.0; block(1, 1) = 3.0
    val solved = sparse(1, 2, (0, 0), (0, 1))
    solved.update(0, 0, 4.0); solved.update(0, 1, 11.0)
    BlockedCholesky.rdivLowerTranspose(MatrixBlock.Sparse(solved), MatrixBlock.BlockDiagonal(Array(block)))
    assertEquals(solved(0, 0), 2.0)
    assertEquals(solved(0, 1), 3.0)

  test("sparse block solve skips entirely absent row-level rectangles"):
    val source = sparse(3, 4, (0, 0), (0, 1), (1, 2), (1, 3))
    source.update(0, 0, 4.0); source.update(0, 1, 11.0)
    source.update(1, 2, 8.0); source.update(1, 3, 17.0)
    val dense = source.asDense
    val lower = WorkMat.zeros(2, 2)
    lower(0, 0) = 2.0; lower(1, 0) = 1.0; lower(1, 1) = 3.0
    val divisor = MatrixBlock.BlockDiagonal(Array(lower, lower.cloneMat))
    BlockedCholesky.rdivLowerTranspose(MatrixBlock.Sparse(source), divisor)
    BlockedCholesky.rdivLowerTranspose(MatrixBlock.Dense(dense), divisor)
    var row = 0
    while row < 3 do
      var col = 0
      while col < 4 do
        assertEquals(source(row, col), dense(row, col))
        col += 1
      row += 1
    assert(!source.hasEntry(2, 0))
    assert(!source.hasEntry(2, 2))

  test("dense right divisor promotes sparse target before mutation"):
    val target = sparse(1, 2, (0, 0))
    target.update(0, 0, 2.0)
    val divisor = WorkMat.zeros(2, 2)
    divisor(0, 0) = 1.0; divisor(1, 0) = 1.0; divisor(1, 1) = 1.0
    val prepared = BlockedCholesky.prepareRightSolveTarget(MatrixBlock.Sparse(target), MatrixBlock.Dense(divisor))
    assert(prepared.isInstanceOf[MatrixBlock.Dense])
    BlockedCholesky.rdivLowerTranspose(prepared, MatrixBlock.Dense(divisor))
    assertEquals(prepared(0, 0), 2.0)
    assertEquals(prepared(0, 1), -2.0)
    assertEquals(target(0, 0), 2.0)
