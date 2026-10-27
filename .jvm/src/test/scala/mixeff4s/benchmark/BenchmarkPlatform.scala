package mixeff4s.benchmark

object BenchmarkPlatform:
  val enabled: Boolean = sys.env.get("MIXEFF4S_BENCHMARK").contains("1")
  val name = "jvm"
  val runtime: String = System.getProperty("java.vm.name") + " " + System.getProperty("java.version")
  private val bean = java.lang.management.ManagementFactory.getThreadMXBean match
    case value: com.sun.management.ThreadMXBean if value.isThreadAllocatedMemorySupported =>
      if !value.isThreadAllocatedMemoryEnabled then value.setThreadAllocatedMemoryEnabled(true)
      Some(value)
    case _ => None
  def allocatedBytes: Long = bean.fold(-1L)(_.getThreadAllocatedBytes(Thread.currentThread().threadId()))
  def heapBytes: Long = Runtime.getRuntime.totalMemory() - Runtime.getRuntime.freeMemory()
