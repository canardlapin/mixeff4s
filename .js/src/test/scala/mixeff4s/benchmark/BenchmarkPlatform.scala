package mixeff4s.benchmark

import scala.scalajs.js

object BenchmarkPlatform:
  private val process = js.Dynamic.global.process
  val enabled: Boolean = process.env.MIXEFF4S_BENCHMARK.asInstanceOf[js.UndefOr[String]].contains("1")
  val glmmResponseCount: Int =
    process.env.MIXEFF4S_GLMM_BENCH_RESPONSES.asInstanceOf[js.UndefOr[String]].toOption.fold(12)(_.toInt)
  val name = "js"
  val runtime: String = process.version.asInstanceOf[String]
  def allocatedBytes: Long = -1L
  def heapBytes: Long = process.memoryUsage().heapUsed.asInstanceOf[Double].toLong
