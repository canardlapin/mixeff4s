package mixeff4s.fixtures

import mixeff4s.data.ModelFrame

object Penicillin:
  private val plates: Vector[String] =
    ('a' to 'x').map(_.toString).toVector

  private val samples: Vector[String] =
    Vector("A", "B", "C", "D", "E", "F")

  val diameter: Vector[Double] = Vector(
    27, 23, 26, 23, 23, 21, 27, 23, 26, 23, 23, 21, 25, 21, 25, 24, 24, 20, 26, 23, 25, 23, 23, 20,
    25, 22, 26, 22, 23, 20, 24, 22, 25, 23, 22, 19, 24, 20, 23, 21, 22, 19, 26, 22, 26, 24, 24, 21,
    24, 21, 24, 22, 22, 20, 24, 21, 24, 23, 22, 19, 26, 23, 26, 24, 24, 21, 25, 22, 26, 24, 24, 20,
    26, 24, 26, 24, 25, 22, 26, 23, 26, 23, 23, 20, 26, 23, 25, 24, 24, 22, 25, 22, 25, 23, 23, 20,
    25, 21, 24, 23, 23, 20, 25, 22, 24, 23, 23, 19, 24, 21, 23, 21, 21, 19, 26, 23, 26, 24, 24, 21,
    25, 21, 24, 22, 22, 18, 25, 22, 25, 22, 22, 20, 24, 21, 24, 22, 24, 19, 24, 21, 24, 22, 21, 18
  ).map(_.toDouble)

  val plate: Vector[String] =
    plates.flatMap(p => Vector.fill(6)(p))

  val sample: Vector[String] =
    Vector.fill(24)(samples).flatten

  def frame: ModelFrame =
    ModelFrame
      .of(
        "diameter" -> ModelFrame.numeric(diameter),
        "plate" -> ModelFrame.factor(plate),
        "sample" -> ModelFrame.factor(sample)
      )
      .fold(err => throw IllegalStateException(err.message), identity)
