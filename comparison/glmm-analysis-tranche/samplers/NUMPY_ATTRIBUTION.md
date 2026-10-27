# Count-sampler source attribution

`src/main/scala/mixeff4s/glmm/CountSamplers.scala` adapts the binomial BTPE
and Poisson PTRS/inversion decision paths from NumPy's
`numpy/random/src/distributions/distributions.c`, pinned at commit
`b833c1b9bf6712cc8a03bb301b77c48984a2373d` (retrieved 2026-10-27).

Source SHA-256: `8526b63d5dc3541fad81505e4ad6cf59dd78722aa35abf256337f6447e46be22`.
The relevant upstream functions are `random_binomial_btpe`,
`random_binomial_inversion`, `random_binomial`, `random_poisson_mult`, and
`random_poisson_ptrs`. In particular, BTPE's third and fourth Stirling terms
are subtracted, matching NumPy's correction note.

NumPy random is dual-licensed under the University of Illinois/NCSA Open
Source License and the 3-Clause BSD License. This adapted source retains the
required attribution and disclaimer through this notice; see
`NUMPY_RANDOM_LICENSE.md`.

