package mixeff4s.glmm

import mixeff4s.design.CompiledDesign
import mixeff4s.model.{Family, Link}

/** Private immutable training contract. No mutable numerical workspace is retained. */
final class GlmmFitContext private[glmm] (
    private[mixeff4s] val design: CompiledDesign,
    private[mixeff4s] val response: GlmmResponse,
    private[mixeff4s] val options: GlmmOptions,
    private[mixeff4s] val family: Family,
    private[mixeff4s] val link: Link
):
  private[mixeff4s] def compiled: CompiledGlmm = new CompiledGlmm(design, family, link)
