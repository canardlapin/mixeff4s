package mixeff4s.glmm

import mixeff4s.data.ModelFrame
import mixeff4s.formula.Formula
import mixeff4s.model.Family

/** Independent dense Newton/Cholesky oracle, assembled from raw rows. */
class JointDenseSuite extends munit.FunSuite:
  test("crossed vector terms and reduced-rank covariance match independent dense conditional solve"):
    val x = Vector.tabulate(18)(i => (i % 3 - 1).toDouble)
    val ys = Vector.tabulate(18)(i => if Set(1, 4, 5, 7, 9, 10, 12, 15, 16).contains(i) then 1.0 else 0.0)
    val groups = Vector.tabulate(18)(i => s"g${i / 6}")
    val crossed = Vector.tabulate(18)(i => s"h${i % 2}")
    val frame = ModelFrame
      .of(
        "y" -> ModelFrame.numeric(ys),
        "x" -> ModelFrame.numeric(x),
        "g" -> ModelFrame.factor(groups),
        "h" -> ModelFrame.factor(crossed)
      )
      .toOption
      .get
    val formula = Formula.parse("y ~ 1 + x + (1 + x | g) + (1 | h)").toOption.get
    Vector(0.5, 0.0).foreach: slopeDiagonal =>
      val theta = Vector(0.8, -0.25, slopeDiagonal, 0.6)
      val actual = Glmm
        .jointDeviance(formula, frame, Family.Bernoulli, Vector(-0.3, 0.4), theta)
        .fold(e => fail(e.message), identity)
      val a = Array.tabulate(18, 8): (i, j) =>
        if j < 6 && j / 2 == i / 6 then if j % 2 == 0 then 0.8 - 0.25 * x(i) else slopeDiagonal * x(i)
        else if j >= 6 && j - 6 == i % 2 then 0.6
        else 0.0
      val u = Array.fill(8)(0.0)
      def state(): (Array[Double], Array[Array[Double]], Array[Double]) =
        val eta = Array.tabulate(18)(i => -0.3 + 0.4 * x(i) + (0 until 8).map(j => a(i)(j) * u(j)).sum)
        val mu = eta.map(e => 1.0 / (1.0 + math.exp(-e)))
        val h = Array.tabulate(8, 8): (j, k) =>
          (if j == k then 1.0 else 0.0) + (0 until 18).map(i => a(i)(j) * a(i)(k) * mu(i) * (1.0 - mu(i))).sum
        val score = Array.tabulate(8)(j => u(j) + (0 until 18).map(i => a(i)(j) * (mu(i) - ys(i))).sum)
        (eta, h, score)
      def cholesky(h: Array[Array[Double]]): Array[Array[Double]] =
        val l = Array.fill(8, 8)(0.0)
        (0 until 8).foreach: i =>
          (0 to i).foreach: j =>
            val residual = h(i)(j) - (0 until j).map(k => l(i)(k) * l(j)(k)).sum
            l(i)(j) = if i == j then math.sqrt(residual) else residual / l(j)(j)
        l
      (0 until 12).foreach: _ =>
        val (_, h, score) = state()
        val l = cholesky(h)
        val step = score.clone()
        (0 until 8).foreach(i => step(i) = (step(i) - (0 until i).map(j => l(i)(j) * step(j)).sum) / l(i)(i))
        (0 until 8).reverse.foreach(i =>
          step(i) = (step(i) - (i + 1 until 8).map(j => l(j)(i) * step(j)).sum) / l(i)(i)
        )
        (0 until 8).foreach(i => u(i) -= step(i))
      val (eta, h, score) = state()
      assert(score.map(math.abs).max < 1e-12)
      val logdet = 2.0 * cholesky(h).indices.map(i => math.log(cholesky(h)(i)(i))).sum
      val objective = eta.indices.map(i => 2.0 * (math.log1p(math.exp(eta(i))) - ys(i) * eta(i))).sum +
        u.map(v => v * v).sum + logdet
      assertEqualsDouble(actual.objective, objective, 1e-9)
      assertEqualsDouble(actual.logdet, logdet, 1e-9)
      actual.sphericalModes.flatten.zip(u).foreach((observed, expected) => assertEqualsDouble(observed, expected, 1e-8))
