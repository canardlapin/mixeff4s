package mixeff4s.fixtures

import mixeff4s.data.ModelFrame

object Pastes:
  val strength: Vector[Double] = Vector(
    62.8, 62.6, 60.1, 62.3, 62.7, 63.1, 60.0, 61.4, 57.5, 56.9, 61.1, 58.9, 58.7, 57.5, 63.9, 63.1,
    65.4, 63.7, 57.1, 56.4, 56.9, 58.6, 64.7, 64.5, 55.1, 55.1, 54.7, 54.2, 58.8, 57.5, 63.4, 64.9,
    59.3, 58.1, 60.5, 60.0, 62.5, 62.6, 61.0, 58.7, 56.9, 57.7, 59.2, 59.4, 65.2, 66.0, 64.8, 64.1,
    54.8, 54.8, 64.0, 64.0, 57.7, 56.8, 58.3, 59.3, 59.2, 59.2, 58.9, 56.6
  )

  val batch: Vector[String] =
    "ABCDEFGHIJ".map(_.toString).toVector.flatMap(b => Vector.fill(6)(b))

  val cask: Vector[String] =
    Vector.fill(10)(Vector("a", "a", "b", "b", "c", "c")).flatten

  def frame: ModelFrame =
    ModelFrame
      .of(
        "strength" -> ModelFrame.numeric(strength),
        "batch" -> ModelFrame.factor(batch),
        "cask" -> ModelFrame.factor(cask)
      )
      .fold(err => throw IllegalStateException(err.message), identity)
