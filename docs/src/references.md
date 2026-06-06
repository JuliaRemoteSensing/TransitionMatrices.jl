# Methods & references

Each solver and numerical trick follows the published literature below; the
in-source docstrings and comments point to the specific equations used. (This
list mirrors the `Methods & references` section of the repository `README.md`,
which is the authoritative copy.)

**T-Matrix framework & far-field conventions**

- P. C. Waterman, *Symmetry, unitarity, and geometry in electromagnetic
  scattering*, [Phys. Rev. D **3**, 825–839 (1971)](https://doi.org/10.1103/PhysRevD.3.825)
  — the null-field / EBCM origin of the T-Matrix method.
- M. I. Mishchenko, L. D. Travis & A. A. Lacis, *Scattering, Absorption, and
  Emission of Light by Small Particles*, Cambridge University Press (2002) — the
  amplitude / phase / scattering matrices, cross sections, asymmetry parameter,
  T-Matrix rotation, analytic random-orientation average, and the Wigner-d
  recursions (the pervasive `Eq. (x.y)` references throughout the source).

**Mie & coated spheres** (`bhmie`, `bhcoat`)

- C. F. Bohren & D. R. Huffman, *Absorption and Scattering of Light by Small
  Particles*, Wiley (1983) — the `bhmie` and `bhcoat` algorithms. The Mie
  T-Matrix uses Mishchenko et al. (2002), Eqs. (5.42)–(5.44).

**EBCM for axisymmetric shapes** (spheroids, cylinders, Chebyshev particles)

- P. C. Waterman (1971), above.
- M. I. Mishchenko & L. D. Travis, *Capabilities and limitations of a current
  FORTRAN implementation of the T-matrix method for randomly oriented,
  rotationally symmetric scatterers*, JQSRT **60**, 309–324 (1998) — the automatic
  convergence procedure (`routine_mishchenko`: the choice of `nₘₐₓ` and the
  Gauss division `Ng`), which the axisymmetric assembly is translated from.

**Numerically stable EBCM — the F⁺ formulation** (`stable = true`)

- W. R. C. Somerville, B. Auguié & E. C. Le Ru, *Severe loss of precision in
  calculations of T-matrix integrals*, JQSRT **113**, 524 (2012) — the
  catastrophic-cancellation diagnosis.
- W. R. C. Somerville, B. Auguié & E. C. Le Ru, JQSRT **123**, 153 (2013),
  [doi:10.1016/j.jqsrt.2012.07.017](https://doi.org/10.1016/j.jqsrt.2012.07.017)
  — the cancellation-free `F⁺_{nk}(s,x)` projection (their Eq. 45–62, Table 2)
  used for high-aspect-ratio spheroids.

**Sh-matrix moment separation** (`prepare_sh`, fast `(λ, mᵣ)` sweeps)

- The separation of the size / material parameters from the particle geometry is
  the *parameter separation* of V. G. Farafonov, V. B. Il'in & M. S. Prokopjeva,
  *Light scattering by multilayered nonspherical particles: a set of methods*,
  JQSRT **79–80**, 599–626 (2003),
  [doi:10.1016/S0022-4073(02)00310-2](https://doi.org/10.1016/S0022-4073(02)00310-2);
  the "Sh-matrix" name and formalism are due to D. Petrov, Yu. Shkuratov, E. Zubko
  & G. Videen, *Sh-matrices method as applied to scattering by particles with
  layered structure*, JQSRT **106**, 437–454 (2007),
  [doi:10.1016/j.jqsrt.2007.01.027](https://doi.org/10.1016/j.jqsrt.2007.01.027)
  (extended in Petrov, Shkuratov & Videen, *The Sh-matrices method applied to
  light scattering by small lenses*, JQSRT **110**, 1448–1459 (2009),
  [doi:10.1016/j.jqsrt.2009.01.016](https://doi.org/10.1016/j.jqsrt.2009.01.016)).
- The term-by-term radial integration reuses the `F⁺` projection of Somerville
  et al. (2013) and the Riccati–Bessel power series of
  [DLMF §10.53](https://dlmf.nist.gov/10.53).

**Invariant Imbedding T-Matrix (IITM)** (axisymmetric, N-fold, and arbitrary shapes)

- B. R. Johnson, *Invariant imbedding T matrix approach to electromagnetic
  scattering*, [Appl. Opt. **27**, 4861–4873 (1988)](https://opg.optica.org/ao/abstract.cfm?uri=ao-27-23-4861)
  — Eq. (97).
- L. Bi, P. Yang, G. W. Kattawar & M. I. Mishchenko, *Efficient implementation of
  the invariant imbedding T-matrix method and the separation of variables method
  applied to large nonspherical inhomogeneous particles*, JQSRT **116**, 169–183
  (2013) — Eq. (38).
- A. Doicu & T. Wriedt, *The Invariant Imbedding T Matrix Approach*, in *The
  Generalized Multipole Technique for Light Scattering* (Springer, 2018),
  [doi:10.1007/978-3-319-74890-0_2](https://doi.org/10.1007/978-3-319-74890-0_2)
  — Eq. (2.40).
- B. Sun, L. Bi, P. Yang, M. Kahnert & G. Kattawar, *Invariant Imbedding T-matrix
  Method for Light Scattering by Nonspherical and Inhomogeneous Particles*,
  Elsevier (2019) — Eq. (4.2.36).
- S. Hu (胡帅), *Research on the Numerical Computational Models and Application of
  the Scattering Properties of Nonspherical Atmospheric Particles*, PhD
  dissertation, National University of Defense Technology (2018) — Eq. (5.71).

**Far-field observables & orientation averaging**

- Mishchenko et al. (2002) — cross sections (Eqs. (5.102), (5.107), (5.140), (5.141)),
  phase matrix (Eqs. (2.106)–(2.121)), asymmetry parameter (Eq. (4.92)), and the
  analytic random-orientation T-Matrix (Eq. (5.96)).
- L. Bi & P. Yang, *Accurate simulation of the optical properties of atmospheric
  ice crystals with the invariant imbedding T-matrix method*, JQSRT **138**, 17–35
  (2014), [doi:10.1016/j.jqsrt.2014.01.013](https://doi.org/10.1016/j.jqsrt.2014.01.013)
  — the scattering-matrix expansion coefficients for a general T-Matrix (Eqs. (24)–(74)).

**Shapes**

- A. Mugnai & W. J. Wiscombe, *Scattering from nonspherical Chebyshev particles. 1*,
  [Appl. Opt. **25**, 1235 (1986)](https://opg.optica.org/ao/abstract.cfm?uri=ao-25-7-1235)
  — the Chebyshev-particle definition `r(ϑ) = r₀(1 + ε·Tₙ(cos ϑ))`.

**Linearization & Jacobians**

- R. Spurr, J. Wang, J. Zeng & M. I. Mishchenko, *Linearized T-matrix and Mie
  scattering computations*, JQSRT (2012).
- F. Xu & A. B. Davis, *Derivatives of light scattering properties of a
  nonspherical particle computed with the T-matrix method*, Opt. Lett. **36**(22),
  4464–4466 (2011), [doi:10.1364/OL.36.004464](https://doi.org/10.1364/OL.36.004464).
- B. Sun, M. Gao, L. Bi & R. Spurr, *Analytical Jacobians of single scattering
  optical properties using the invariant imbedding T-matrix method*, Opt. Express
  **29**(6), 9635–9669 (2021), [doi:10.1364/OE.421886](https://doi.org/10.1364/OE.421886).
- M. Gao & B. Sun, *Improvement and application of linearized invariant imbedding
  T-matrix scattering method*, JQSRT **290**, 108322 (2022),
  [doi:10.1016/j.jqsrt.2022.108322](https://doi.org/10.1016/j.jqsrt.2022.108322).

**Numerical kernels (via dependencies)**

- Wigner 3-j symbols via `Wigxjpf.jl`: H. T. Johansson & C. Forssén, *Fast and
  accurate evaluation of Wigner 3j, 6j, and 9j symbols using prime factorization
  and multiword integer arithmetic*, SIAM J. Sci. Comput. **38**(1), A376–A384 (2016).
- Gauss–Legendre quadrature via
  [`FastGaussQuadrature.jl`](https://github.com/JuliaApproximation/FastGaussQuadrature.jl).
