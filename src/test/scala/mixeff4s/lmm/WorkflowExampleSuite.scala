package mixeff4s.lmm
import mixeff4s.examples.LmmWorkflow
class WorkflowExampleSuite extends munit.FunSuite:
  test("documented public workflow fits, assesses, predicts, simulates and produces qualified intervals"):
    val result = LmmWorkflow.run().fold(e => fail(s"${e.code}: ${e.message}"), identity)
    assert(result.assessment.eligibleForBootstrap)
    assertEquals(result.predictions.length, 3)
    assertEquals(result.residuals.length, 180)
    assertEquals(result.simulated.length, 180)
    assertEquals(result.bootstrap.qualifiedCount, 39)
    val interval = result.interval.fold(e => fail(e.message), identity)
    assert(interval.lower < interval.upper)
