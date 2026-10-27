package mixeff4s.lmm

import mixeff4s.design.CompiledDesign

/** Private construction seals the complete fitted snapshot. No mutable workspace is retained. */
final class LmmFitContext private[lmm] (
    private[mixeff4s] val design: CompiledDesign,
    private val expected: LmmFit,
    private[mixeff4s] val modes: Vector[Vector[Double]]
):
  private[mixeff4s] def matches(fit: LmmFit): Boolean = fit.copy(workflow = None) == expected
  private[mixeff4s] def compiled: CompiledLmm = new CompiledLmm(design)
