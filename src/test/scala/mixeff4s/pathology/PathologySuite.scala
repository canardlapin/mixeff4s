package mixeff4s.pathology

import mixeff4s.data.ModelFrame
import mixeff4s.design.Design
import mixeff4s.error.MixedModelError
import mixeff4s.fixtures.Sleepstudy
import mixeff4s.formula.Formula

class PathologySuite extends munit.FunSuite:
  test("sleepstudy is an easy design-time certificate"):
    val cert = certify("reaction ~ 1 + days + (1 + days | subj)", Sleepstudy.frame)
    assertEquals(cert.stratum, Stratum.Easy)
    assertEquals(cert.fitStatus, FitStatus.NotAssessed)
    assertEquals(cert.expectedStatuses, Vector(FitStatus.ConvergedInterior))
    assertEquals(cert.n, 180)
    assertEquals(cert.nParams, 6)
    assertEquals(cert.minGroupSize, 10)
    assertEquals(cert.maxGroupSize, 10)
    assertEquals(cert.structuralIssue, None)

  test("no random effects maps to not_identifiable"):
    val frame = ModelFrame
      .of(
        "y" -> ModelFrame.numeric(Vector(1.0, 2.0, 3.0)),
        "x" -> ModelFrame.numeric(Vector(0.0, 1.0, 0.0))
      )
      .getOrElse(fail("frame"))
    compile("y ~ 1 + x", frame) match
      case Left(err) =>
        val cert = Pathology.fromError(err)
        assertEquals(err, MixedModelError.NoRandomEffects)
        assertEquals(cert.fitStatus, FitStatus.NotIdentifiable)
        assertEquals(cert.stratum, Stratum.Refusal)
      case Right(other) =>
        fail(s"expected no-RE refusal, got $other")

  test("unsupported covariance maps to not_optimized"):
    compile("reaction ~ 1 + days + cs(1 + days | subj)", Sleepstudy.frame) match
      case Left(err: MixedModelError.Unsupported) =>
        assertEquals(Pathology.mapError(err), FitStatus.NotOptimized)
      case other =>
        fail(s"expected unsupported covariance, got $other")

  test("singletons with a random slope are a refusal certificate"):
    val frame = ModelFrame
      .of(
        "y" -> ModelFrame.numeric(Vector(1.0, 2.0, 3.0, 4.0)),
        "x" -> ModelFrame.numeric(Vector(0.0, 1.0, 0.0, 1.0)),
        "g" -> ModelFrame.factor(Vector("a", "b", "c", "d"))
      )
      .getOrElse(fail("frame"))
    val cert = certify("y ~ 1 + x + (1 + x | g)", frame)
    val issue = cert.structuralIssue.getOrElse(fail("issue"))
    assertEquals(cert.stratum, Stratum.Refusal)
    assertEquals(issue.code, "singletons_with_slope")
    assertEquals(cert.expectedStatuses, Vector(FitStatus.NotIdentifiable, FitStatus.NotOptimized))

  test("a single grouping level is a refusal certificate"):
    val frame = ModelFrame
      .of(
        "y" -> ModelFrame.numeric(Vector(1.0, 2.0, 3.0)),
        "g" -> ModelFrame.factor(Vector("a", "a", "a"))
      )
      .getOrElse(fail("frame"))
    val cert = certify("y ~ 1 + (1 | g)", frame)
    val issue = cert.structuralIssue.getOrElse(fail("issue"))
    assertEquals(issue.code, "few_random_effect_levels")

  test("interior theta is converged_interior"):
    val cert = certify("reaction ~ 1 + days + (1 + days | subj)", Sleepstudy.frame)
    val assessed =
      Pathology.assessFit(cert, Vector(0.9292, 0.0182, 0.2226), Vector((0, 0, 0), (0, 1, 0), (0, 1, 1)))
    assertEquals(assessed.fitStatus, FitStatus.ConvergedInterior)
    assertEquals(assessed.expectedStatuses, Vector(FitStatus.ConvergedInterior))

  test("a diagonal theta at the lower bound is converged_boundary"):
    val cert = certify("reaction ~ 1 + days + (1 + days | subj)", Sleepstudy.frame)
    val assessed =
      Pathology.assessFit(cert, Vector(0.0, 0.1, 0.2), Vector((0, 0, 0), (0, 1, 0), (0, 1, 1)))
    assertEquals(assessed.fitStatus, FitStatus.ConvergedBoundary)

  test("an off-diagonal near zero is not a variance boundary"):
    val cert = certify("reaction ~ 1 + days + (1 + days | subj)", Sleepstudy.frame)
    val assessed =
      Pathology.assessFit(cert, Vector(0.9, 0.0, 0.2), Vector((0, 0, 0), (0, 1, 0), (0, 1, 1)))
    assertEquals(assessed.fitStatus, FitStatus.ConvergedInterior)

  test("refusal designs do not claim a converged fit status"):
    val frame = ModelFrame
      .of(
        "y" -> ModelFrame.numeric(Vector(1.0, 2.0, 3.0)),
        "g" -> ModelFrame.factor(Vector("a", "a", "a"))
      )
      .getOrElse(fail("frame"))
    val cert = certify("y ~ 1 + (1 | g)", frame)
    val assessed = Pathology.assessFit(cert, Vector(0.5), Vector((0, 0, 0)))
    assertEquals(assessed.fitStatus, FitStatus.NotAssessed)
    assert(assessed.notes.exists(_.contains("refusal")), clues(assessed.notes))

  private def compile(source: String, frame: ModelFrame) =
    Formula.parse(source) match
      case Left(err)      => Left(MixedModelError.Formula(err))
      case Right(formula) => Design.compile(formula, frame)

  private def certify(source: String, frame: ModelFrame): Certificate =
    compile(source, frame).map(Pathology.certify).getOrElse(fail(s"compile $source"))
