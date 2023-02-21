var documenterSearchIndex = {"docs":
[{"location":"api/","page":"API","title":"API","text":"CurrentModule = TransitionMatrices","category":"page"},{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"","category":"page"},{"location":"api/","page":"API","title":"API","text":"Modules = [TransitionMatrices]","category":"page"},{"location":"api/#TransitionMatrices.AbstractTransitionMatrix","page":"API","title":"TransitionMatrices.AbstractTransitionMatrix","text":"A general T-Matrix T_m n m^prime n^prime^k l stored in a 6-dimensional array, in the order (m n m^prime n^prime k l).\n\n\n\n\n\n","category":"type"},{"location":"api/#TransitionMatrices.Chebyshev","page":"API","title":"TransitionMatrices.Chebyshev","text":"A Chebyshev scatterer defined by\n\nr(theta phi)=r_0(1+varepsilon T_n(costheta))\n\nwhere T_n(costheta)=cos ntheta.\n\nAttributes:\n\nr₀: radius of the base sphere.\nε: deformation parameter, which satisfies -1levarepsilon1.\nn: degree of the Chebyshev polynomial.\nm: relative complex refractive index.\n\n\n\n\n\n","category":"type"},{"location":"api/#TransitionMatrices.Cylinder","page":"API","title":"TransitionMatrices.Cylinder","text":"A cylindrical scatterer.\n\nAttributes:\n\nr: radius of the cylinder base\nh: height of the cylinder\nm: relative complex refractive index\n\n\n\n\n\n","category":"type"},{"location":"api/#TransitionMatrices.MieTransitionMatrix","page":"API","title":"TransitionMatrices.MieTransitionMatrix","text":"According to Eq. (5.42) – Eq. (5.44) in Mishchenko et al. (2002), the T-Matrix for a Mie particle can be written as:\n\nbeginarrayl\nT_m n m^prime n^prime^12(P) equiv 0 quad T_m n m^prime n^prime^21(P) equiv 0 \nT_m n m^prime n^prime^11(P)=-delta_m m^prime delta_n n^prime b_n \nT_m n m^prime n^prime^22(P)=-delta_m m^prime delta_n n^prime a_n \nendarray\n\nMieTransitionMatrix{CT, N}(x::Real, m::Number)\n\nGenerate the T-Matrix from the Mie coefficients of a homogeneous sphere.\n\nMieTransitionMatrix{CT, N}(x_core::Real, x_mantle::Real, m_core::Number, m_mantle::Number)\n\nGenerate the T-Matrix from the Mie coefficients of a coated sphere.\n\nThis struct provides the T-Matrix API for a Mie particle.\n\n\n\n\n\n","category":"type"},{"location":"api/#TransitionMatrices.OrderDegreeIterator","page":"API","title":"TransitionMatrices.OrderDegreeIterator","text":"Iterator for the order-degree pairs of the given maximum order nₘₐₓ.\n\nExample of nₘₐₓ=2:\n\njulia> collect(OrderDegreeIterator(2))\n8-element Vector{Tuple{Int64, Int64}}:\n (1, -1)\n (1, 0)\n (1, 1)\n (2, -2)\n (2, -1)\n (2, 0)\n (2, 1)\n (2, 2)\n\n\n\n\n\n","category":"type"},{"location":"api/#TransitionMatrices.RandomOrientationTransitionMatrix-Union{Tuple{AbstractTransitionMatrix{CT, N}}, Tuple{N}, Tuple{CT}} where {CT, N}","page":"API","title":"TransitionMatrices.RandomOrientationTransitionMatrix","text":"RandomOrientationTransitionMatrix(𝐓::AbstractTransitionMatrix{CT, N}) where {CT, N}\n\nCalculate the random-orientation T-Matrix according to Eq. (5.96) in Mishchenko et al. (2002).\n\nleftlangle T_m n m^prime n^prime^k l(L)rightrangle=fracdelta_n n^prime delta_m m^prime2 n+1 sum_m_1=-n^n T_m_1 n m_1 n^k l(P) quad k l=12\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionMatrices.Spheroid","page":"API","title":"TransitionMatrices.Spheroid","text":"A spheroidal scatterer.\n\nAttributes:\n\na: length of the semi-major axis\nc: length of the semi-minor axis\nm: relative complex refractive index\n\n\n\n\n\n","category":"type"},{"location":"api/#TransitionMatrices.TransitionMatrix","page":"API","title":"TransitionMatrices.TransitionMatrix","text":"Concrete type for a general T-Matrix.\n\n\n\n\n\n","category":"type"},{"location":"api/#TransitionMatrices.amplitude_matrix-Union{Tuple{N}, Tuple{CT}, Tuple{AbstractTransitionMatrix{CT, N}, Any, Any, Any, Any}, Tuple{AbstractTransitionMatrix{CT, N}, Any, Any, Any, Any, Any}} where {CT, N}","page":"API","title":"TransitionMatrices.amplitude_matrix","text":"Calculate the amplitude matrix of the given T-Matrix 𝐓 at the given incidence and scattering angles. k₁ is the wavenumber of the incident wave in the host medium, which should be calculated by k₁ = 2π * mₕ / λ, where mₕ is the refractive index of the host medium and λ is the wavelength of the incident wave. The default value is k₁ = 1.0.\n\nGeneral T-Matrix\n\namplitude_matrix(𝐓::AbstractTransitionMatrix{CT, N}, ϑᵢ, φᵢ, ϑₛ, φₛ, k₁=1.0)\n\nFor a general T-Matrix, Eq. (5.11) – Eq. (5.17) in Mishchenko et al. (2002) is used as a fallback.\n\nbeginarrayl\nS_11left(hatmathbfn^text sca  hatmathbfn^text inc right)=frac1k_1 sum_n=1^infty sum_n^prime=1^infty sum_m=-n^n sum_m^prime=-n^prime^n^prime alpha_m n m^prime n^primeleftT_m n m^prime n^prime^11 pi_m nleft(vartheta^text sca right) pi_m^prime n^primeleft(vartheta^text inc right)right \n+T_m n m^prime n^prime^21 tau_m nleft(vartheta^text sca right) pi_m^prime n^primeleft(vartheta^mathrmincright)+T_m n m^prime n^prime^12 pi_m nleft(vartheta^text sca right) tau_m^prime n^primeleft(vartheta^mathrmincright) \nleft+T_m n m^prime n^2^22 tau_m nleft(vartheta^text sca right) tau_m^prime n^primeleft(vartheta^text inc right)right exp leftmathrmileft(m varphi^text sca -m^prime varphi^text inc right)right text   \nS_12left(hatmathbfn^mathrmsca hatmathbfn^mathrmincright)=frac1mathrmi k_1 sum_n=1^infty sum_n^prime=1^infty sum_m=-n^n sum_m^prime=-n^prime^n^prime alpha_m n m^prime n^primeleftT_m n m^prime n^prime^11 pi_m nleft(vartheta^mathrmscaright) tau_m^prime n^primeleft(vartheta^mathrmincright)right \n+T_m n m^prime n^prime^21 tau_m nleft(vartheta^text sca right) tau_m^prime n^primeleft(vartheta^mathrmincright)+T_m n m^prime n^prime^12 pi_m nleft(vartheta^text sca right) pi_m^prime n^primeleft(vartheta^text inc right) \nleft+T_m n m^prime n^prime^22 tau_m nleft(vartheta^text sca right) pi_m^prime n^primeleft(vartheta^text inc right)right exp leftmathrmileft(m varphi^text sca -m^prime varphi^text inc right)right text   \nS_21left(hatmathbfn^mathrmsca hatmathbfn^mathrmincright)=fracmathrmik_1 sum_n=1^infty sum_n^prime=1^infty sum_m=-n^n sum_m^prime=-n^prime^n^prime alpha_m n m^prime n^primeleftT_m n m^prime n^prime^11 tau_m nleft(vartheta^mathrmscaright) pi_m^prime n^primeleft(vartheta^mathrmincright)right \n+T_m n m^prime n^prime^21 pi_m nleft(vartheta^text sca right) pi_m^prime n^primeleft(vartheta^mathrmincright)+T_m n m^prime n^prime^12 tau_m nleft(vartheta^text sca right) tau_m^prime n^primeleft(vartheta^text inc right) \nleft+T_m n m^prime n^prime^22 pi_m nleft(vartheta^text sca right) tau_m^prime n^primeleft(vartheta^text inc right)right exp leftmathrmileft(m varphi^text sca -m^prime varphi^text inc right)right text   \nS_22left(hatmathbfn^text sca  hatmathbfn^mathrmincright)=frac1k_1 sum_n=1^infty sum_n^prime=1^infty sum_m=-n^n sum_m^prime=-n^prime^n^prime alpha_m n m^prime n^primeleftT_m n n^prime n^prime^11 tau_m nleft(vartheta^text sca right) tau_m^prime n^primeleft(vartheta^text inc right)right \n+T_m n m^prime n^prime^21 pi_m nleft(vartheta^text sca right) tau_m^prime n^primeleft(vartheta^text inc right)+T_m n m^prime n^12^12 tau_m nleft(vartheta^text sca right) pi_m^prime n^primeleft(vartheta^text inc right) \nleft+T_m n m^prime n^prime^22 pi_m nleft(vartheta^text sca right) pi_m^prime n^primeleft(vartheta^text inc right)right exp leftmathrmileft(m varphi^text sca -m^prime varphi^text inc right)right \nendarray\n\nWhere\n\nbeginarrayl\nalpha_m n m^prime n^prime=mathrmi^n^prime-n-1(-1)^m+m^primeleftfrac(2 n+1)left(2 n^prime+1right)n(n+1) n^primeleft(n^prime+1right)right^1  2 \npi_m n(vartheta)=fracm d_0 m^n(vartheta)sin vartheta quad pi_-m n(vartheta)=(-1)^m+1 pi_m n(vartheta) \ntau_m n(vartheta)=fracmathrmd d_0 m^n(vartheta)mathrmd vartheta quad tau_-m n(vartheta)=(-1)^m tau_m n(vartheta)\nendarray\n\nAxisymmetric T-Matrix\n\nMie T-Matrix\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionMatrices.bhcoat-NTuple{5, Any}","page":"API","title":"TransitionMatrices.bhcoat","text":"bhcoat([T=Float64,], xᵢₙ, xₒᵤₜ, mᵢₙ, mₒᵤₜ; nₘₐₓ, tolerance = 1e-8)\n\nInputs:\n\nT: Type used for calculation. All real numbers will be stored as T, while complex numbers will be stored as C = complex(T).\nxᵢₙ: Size parameter of the inner sphere. Defined as frac2pi rlambda\nxₒᵤₜ: Size parameter of the coated sphere. xₒᵤₜ >= xᵢₙ should hold.\nmᵢₙ: Refractive index of the inner sphere, relative to the host medium.\nmₒᵤₜ: Refractive index of the mantle, relative to the host medium.\n\nKeyword arguments:\n\nnₘₐₓ: Maximum order of the Mie coefficients. Default to max(x_m + 4sqrt3x_m + 2 max(m_c m_m)x_m).\ntolerance: Error tolerance. Default is 1e-8.\n\nOutputs:\n\na, b: Mie coefficients. Both are of type Vector{C}.\n\nReferences:\n\nBohren, C.F., Huffman, D.R., 1983. Absorption and scattering of light by small particles. John Wiley & Sons.\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionMatrices.bhmie-Tuple{Any, Any, Any}","page":"API","title":"TransitionMatrices.bhmie","text":"bhmie([T=Float64,], x, m; nₘₐₓ)\n\nThis is a simplified version of Bohren (1983)'s bhmie routine. Only the Mie coefficients are evaluated and returned.\n\nInputs:\n\nT: Type used for calculation. Real numbers will be stored as T, while complex numbers will be stored as C = complex(T).\nx: Size parameter of the sphere scatterer. Defined as frac2pi rlambda\nm: Relative refractive index of the scatterer.\n\nKeyword arguments:\n\nnₘₐₓ: Maximum order of the Mie coefficients. Default to max(x + 4sqrt3x + 2 mx).\n\nOutputs:\n\na, b: Mie coefficients. Both are Vector{C} with nₘₐₓ elements.\n\nReferences:\n\nBohren, C.F., Huffman, D.R., 1983. Absorption and scattering of light by small particles. John Wiley & Sons.\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionMatrices.clear_factorial_table!-Union{Tuple{Type{T}}, Tuple{T}} where T","page":"API","title":"TransitionMatrices.clear_factorial_table!","text":"clear_factorial_table!(T)\n\nClear the cached factorial table for a specific type. This is necessary after you  change the precision of a dynamic-precision type, e.g., Arb or BigFloat.\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionMatrices.factorial-Union{Tuple{T}, Tuple{Type{T}, Any}} where T","page":"API","title":"TransitionMatrices.factorial","text":"factorial([T=Float64,], n)\n\nCalculate factorials with a global cache for every type.\n\nnote: Note\nThis function is not thread-safe. To use it in a multi-threaded environment, you should pre-compute the maximum factorial you need for a specific type.For dynamic-precision types, e.g., BigFloat or Arb, after you change the  precision, you should call clear_factorial_table!(T) to clear the cache, and then re-compute the maximum factorial you need for the new precision.By default, the maximum factorial for the following types are pre-computed:Float64: 150\nDouble64: 150\nFloat128: 300\nArb: 500\nBigFloat: 500\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionMatrices.gausslegendre-Tuple{Type{Float64}, Integer}","page":"API","title":"TransitionMatrices.gausslegendre","text":"gausslegendre([T=Float64,], n)\n\nCalculate the n-point Gauss-Legendre quadrature nodes and weights.\n\nFor Float64, we use the mathcalO(n) implementation from FastGaussQuadrature.jl.\nFor other types, we use the arbitrary precision implementation from Arblib.jl, then convert the results to the desired type.\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionMatrices.order-Union{Tuple{AbstractTransitionMatrix{CT, N}}, Tuple{N}, Tuple{CT}} where {CT, N}","page":"API","title":"TransitionMatrices.order","text":"Get the maximum order of a T-Matrix.\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionMatrices.orientation_average-Union{Tuple{N}, Tuple{CT}, Tuple{AbstractTransitionMatrix{CT, N}, Any}} where {CT, N}","page":"API","title":"TransitionMatrices.orientation_average","text":"Calculate the orientation average of a transition matrix, given the orientation distribution function p_o(alphabetagamma). \n\nGeneral T-Matrix\n\norientation_average(𝐓::AbstractTransitionMatrix{CT, N}, pₒ; Nα = 10, Nβ = 10, Nγ = 10) where {CT, N}\n\nFor a general T-Matrix and a general orientation distribution function, we use numerical integration to calculate the orientation average. The orientation average is given by\n\nlangle T_m n m^prime n^prime^p p^prime(L)rangle = int_0^2pimathrmdalphaint_0^pimathrmdbetasinbetaint_0^2pimathrmdgamma p_o(alphabetagamma) T_m n m^prime n^prime^p p^prime(L alphabetagamma)\n\nParameters:\n\n𝐓: The T-Matrix to be orientation averaged.\npₒ: The orientation distribution function. Note that the sinbeta part is already included.\nNα: The number of points used in the numerical integration of alpha. Default to 10.\nNβ: The number of points used in the numerical integration of beta. Default to 10.\nNγ: The number of points used in the numerical integration of gamma. Default to 10.\n\nnote: Note\nThis is the fallback method and does not utilize any symmetry, so it is expected to be slow. You should use specified versions of this function, or implement your own if there is no suited version for your combination of T-Matrix and orientation distribution function.You may also need to test the convergence of Nα, Nβ and Nγ manually. If any one is too small, there will be large errors in the results.\n\nRandom orientation T-Matrix and Mie T-Matrix\n\norientation_average(mie::MieTransitionMatrix, _pₒ; _kwargs...)\norientation_average(mie::MieTransitionMatrix, _pₒ; _kwargs...)\n\nBoth types are invariant under rotation. Therefore, the original T-Matrix will be returned.\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionMatrices.phase_matrix-Tuple{AbstractMatrix}","page":"API","title":"TransitionMatrices.phase_matrix","text":"phase_matrix(𝐒::AbstractMatrix)\n\nCalculate the phase matrix 𝐙 from the amplitude matrix 𝐒, according to Eq. (2.106) – Eq. (2.121) in Mishchenko et al. (2002).\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionMatrices.pi_func-Union{Tuple{T}, Tuple{Type{T}, Integer, Integer, Number}} where T","page":"API","title":"TransitionMatrices.pi_func","text":"pi_func([T=Float64,], m::Integer, n::Integer, ϑ::Number; d=nothing)\n\nCalculate\n\npi_m n(vartheta)=fracm d_0 m^n(vartheta)sin vartheta\n\nIf d is given, it is used as the value of d_0 m^n(vartheta).\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionMatrices.ricattibesselj!-NTuple{6, Any}","page":"API","title":"TransitionMatrices.ricattibesselj!","text":"ricattibesselj!(ψₙ, ψₙ′, z, nₘₐₓ, nₑₓₜᵣₐ, x)\n\nParameters:\n\nψ: Output array for psi_n(x).\nψ′: Output array for psi^prime_n(x).\nz: Auxiliary array for downward recursion.\nnₘₐₓ: Maximum order of psi_n(x).\nnₑₓₜᵣₐ: Extra terms for downward recursion of z.\nx: Argument.\n\nIn-place version of ricattibesselj.\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionMatrices.ricattibesselj-Tuple{Integer, Integer, Any}","page":"API","title":"TransitionMatrices.ricattibesselj","text":"ricattibesselj(nₘₐₓ::Integer, nₑₓₜᵣₐ::Integer, x)\n\nCalculate Ricatti-Bessel function of the first kind, psi_n(x), and its derivative psi^prime_n(x) for 1le nle n_max.\n\npsi_n(x) is defined as\n\npsi_n(x) = xj_n(x)\n\nwhere j_n(x) is the Bessel function of the first kind.\n\nParameters:\n\nnₘₐₓ: Maximum order of psi_n(x).\nnₑₓₜᵣₐ: Extra terms for downward recursion of z.\nx: Argument.\n\nReturns: (ψ, ψ′)\n\nψ: Array of psi_n(x).\nψ′: Array of psi^prime_n(x).\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionMatrices.rotate-Union{Tuple{N}, Tuple{CT}, Tuple{AbstractTransitionMatrix{CT, N}, Rotations.Rotation{3}}} where {CT, N}","page":"API","title":"TransitionMatrices.rotate","text":"Rotate the given T-Matrix 𝐓 by the Euler angle rot and generate a new T-Matrix.\n\nGeneral T-Matrix\n\nrotate(𝐓::AbstractTransitionMatrix{CT, N}, rot::Rotation{3})\n\nFor a general T-Matrix, Eq. (5.29) in Mishchenko et al. (2002) is used as a fallback. A TransitionMatrix will be returned, which is the most general yet concrete type.\n\nT_m n m^prime n^prime^p p(L  alpha beta gamma)=sum_m_1=-n^n sum_m_2=-n^prime^n^prime D_m m_1^n(alpha beta gamma) T_m_1 n m_2 n^prime^p p(P) D_m_2 m^prime^n^prime(-gamma-beta-alpha)quad pp=12\n\nAxisymmetric T-Matrix\n\nMie T-Matrix\n\nrotate(𝐓::MieTransitionMatrix{CT, N}, rot::Rotation{3})\n\nFor a MieTransitionMatrix, the underlying Mie coefficients are copied and a new MieTransitionMatrix will be returned.\n\nOrientation-averaged T-Matrix\n\nrotate(oa::OrientationAveragedTransitionMatrix{CT, N, V}, ::Rotation{3}) where {CT, N, V} =\n    typeof(oa)(copy(oa.𝐓))\n\nFor an OrientationAveragedTransitionMatrix, the underlying T-Matrix is copied and a new OrientationAveragedTransitionMatrix will be returned.\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionMatrices.tau_func-Union{Tuple{T}, Tuple{Type{T}, Integer, Integer, Number}} where T","page":"API","title":"TransitionMatrices.tau_func","text":"tau_func([T=Float64,], m::Integer, n::Integer, ϑ::Number)\n\nCalculate\n\ntau_m n(vartheta)=fracmathrmd d_0 m^n(vartheta)mathrmd vartheta\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionMatrices.transition_matrix-Union{Tuple{CT}, Tuple{T}, Tuple{AbstractAxisymmetricShape{T, CT}, Any, Any, Any}} where {T, CT}","page":"API","title":"TransitionMatrices.transition_matrix","text":"transition_matrix(s::AbstractAxisymmetricShape{T, CT}, λ, nₘₐₓ, Ng) where {T, CT}\n\nCalculate the T-Matrix for a given scatterer and wavelength.\n\nParameters:\n\ns: An axisymmetricsScatterer\nλ: Wavelength\nnₘₐₓ: Maximum order of the T-Matrix\nNg: Number of Gauss-Legendre quadrature points\n\nReturns:\n\n𝐓: An AxisymmetricTransitionMatrix struct containing the T-Matrix\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionMatrices.wigner_D-Union{Tuple{T}, Tuple{Type{T}, Integer, Integer, Integer, Number, Number, Number}} where T","page":"API","title":"TransitionMatrices.wigner_D","text":"wigner_D(::Type{T}, m::Integer, m′::Integer, n::Integer, α::Number, β::Number, γ::Number) where {T}\n\nCalculate the Wigner D-function D^j_mn(theta), which is defined as:\n\nD_m m^prime^n(alpha beta gamma)=mathrme^-mathrmi m alpha d_m m^prime^n(beta) mathrme^-mathrmi m^prime gamma\n\nwhere\n\n0 leq alpha2 pi quad 0 leq beta leq pi quad 0 leq gamma2 pi\n\nwarning: Warning\nThis function easily overflows for large values of s, and it is no faster than the recursive method. It is provided here only for checking the correctness of the recursive method. Users are recommended to use wigner_D_recursion instead.\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionMatrices.wigner_D_recursion!-Union{Tuple{CT}, Tuple{AbstractVector{CT}, Integer, Integer, Integer, Number, Number, Number}} where CT","page":"API","title":"TransitionMatrices.wigner_D_recursion!","text":"wigner_D_recursion!(d::AbstractVector{CT}, m::Integer, m′::Integer, nmax::Integer, α::Number, β::Number, γ::Number)\n\nCalculate the Wigner D-function recursively, in place.\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionMatrices.wigner_D_recursion-Union{Tuple{T}, Tuple{Type{T}, Integer, Integer, Integer, Number, Number, Number}} where T","page":"API","title":"TransitionMatrices.wigner_D_recursion","text":"wigner_D_recursion([T=Float64,], m::Integer, m′::Integer, nmax::Integer, α::Number, β::Number, γ::Number)\n\nCalculate the Wigner D-function recursively (use wigner_d_recursion).\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionMatrices.wigner_d-Union{Tuple{T}, Tuple{Type{T}, Integer, Integer, Integer, Number}} where T","page":"API","title":"TransitionMatrices.wigner_d","text":"wigner_d([T=Float64,], m::Integer, n::Integer, s::Integer, ϑ::Number) where {T}\n\nCalculate Wigner (small) d-function d_mn^s(theta) for a single (m n s) combination, using Eq. (B.1) of Mishchenko et al. (2002).\n\nbeginaligned\nd_m n^s(vartheta)= sqrt(s+m) (s-m) (s+n) (s-n)  \n times sum_k=max(0m-n)^min(s + m s - n)(-1)^k fracleft(cos frac12 varthetaright)^2 s-2 k+m-nleft(sin frac12 varthetaright)^2 k-m+nk (s+m-k) (s-n-k) (n-m+k) \nendaligned\n\nwarning: Warning\nThis function easily overflows for large values of s, and it is no faster than the recursive method. It is provided here only for checking the correctness of the recursive method. Users are recommended to use wigner_d_recursion` instead.\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionMatrices.wigner_d_recursion!-Union{Tuple{T}, Tuple{AbstractVector{T}, Integer, Integer, Integer, Number}} where T","page":"API","title":"TransitionMatrices.wigner_d_recursion!","text":"wigner_d_recursion!(d::AbstractVector{T}, m::Integer, n::Integer, smax::Integer, ϑ::Number; deriv=nothing)\n\nCalculate the Wigner d-function recursively, in place.\n\n\n\n\n\n","category":"method"},{"location":"api/#TransitionMatrices.wigner_d_recursion-Union{Tuple{T}, Tuple{Type{T}, Integer, Integer, Integer, Number}} where T","page":"API","title":"TransitionMatrices.wigner_d_recursion","text":"wigner_d_recursion([T=Float64,], m::Integer, n::Integer, smax::Integer, ϑ::Number; deriv::Bool = false) where {T}\n\nCalculate Wigner (small) d-function d_mn^s(theta) for sins_min=max(m n)s_max (alternatively, its derivative as well) via upward recursion, using Eq. (B.22) of Mishchenko et al. (2002).\n\nbeginaligned\nd_m n^s+1(vartheta)= frac1s sqrt(s+1)^2-m^2 sqrt(s+1)^2-n^2left(2 s+1)s(s+1) x-m n d_m n^s(vartheta)right\nleft-(s+1) sqrts^2-m^2 sqrts^2-n^2 d_m n^s-1(vartheta)right quad s geq s_min \nendaligned\n\nThe initial terms are given by Eq. (B.23) and Eq. (B.24).\n\nbeginarrayl\nd_m n^s_min -1(vartheta)=0 \nd_m n^s_min (vartheta)=xi_m n 2^-s_min leftfracleft(2 s_min right) (m-n) (m+n) right^1  2(1-x)^m-n  2(1+x)^m+n  2\nendarray\n\nwhere\n\nxi_m n=leftbeginarrayll\n1  text  for  n geq m \n(-1)^m-n  text  for  nm\nendarrayright\n\n\n\n\n\n","category":"method"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = TransitionMatrices","category":"page"},{"location":"#TransitionMatrices.jl","page":"Home","title":"TransitionMatrices.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The transition matrix method, or T-Matrix method, is one of the most powerful and widely used tools for rigorously computing electromagnetic scattering by single and compounded particles. As a package focusing on this method, TransitionMatrices.jl provides the following features:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Calculate the T-Matrix of various types of scatterers\nHomogeneous sphere (via bhmie)\nCoated sphere (via bhcoat)\nCalculate far-field scattering properties using the T-Matrix","category":"page"}]
}
