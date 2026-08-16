package mixeff4s.fixtures

import mixeff4s.data.ModelFrame

/** MixedModels.jl / mixeff-rs contra fixture (1934 binary observations). */
object Contra:
  def frame: ModelFrame =
    val source = Option(getClass.getResourceAsStream("/mixeff4s/contra.csv")).getOrElse:
      throw IllegalStateException("missing test resource /mixeff4s/contra.csv")
    val lines =
      try scala.io.Source.fromInputStream(source).getLines().toVector
      finally source.close()
    val useNum = Vector.newBuilder[Double]
    val age = Vector.newBuilder[Double]
    val age2 = Vector.newBuilder[Double]
    val urban = Vector.newBuilder[String]
    val livch = Vector.newBuilder[String]
    val urbanDist = Vector.newBuilder[String]
    lines.foreach: line =>
      val parts = line.split(",", -1)
      useNum += parts(0).toDouble
      age += parts(1).toDouble
      age2 += parts(2).toDouble
      urban += parts(3)
      livch += parts(4)
      urbanDist += parts(5)
    ModelFrame
      .of(
        "use_num" -> ModelFrame.numeric(useNum.result()),
        "age" -> ModelFrame.numeric(age.result()),
        "age2" -> ModelFrame.numeric(age2.result()),
        "urban" -> ModelFrame.factor(urban.result()),
        "livch" -> ModelFrame.factor(livch.result()),
        "urban_dist" -> ModelFrame.factor(urbanDist.result())
      )
      .fold(err => throw IllegalStateException(err.message), identity)
