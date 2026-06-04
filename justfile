# Development tasks for TransitionMatrices.jl.

default:
    @just --list

[group("development")]
dev:
    julia --project=. -e 'using Pkg; Pkg.instantiate()'
    julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'

[group("testing")]
test:
    julia --project=. -e 'using Pkg; Pkg.test()'
alias t := test

[group("testing")]
test-coverage:
    julia --project=. -e 'using Pkg; Pkg.test(; coverage=true)'

[group("docs")]
docs:
    julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate(); ENV["DOCUMENTER_SKIP_DEPLOY"]="true"; include("docs/make.jl")'
alias d := docs

[group("docs")]
doctest:
    julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate(); using Documenter: DocMeta, doctest; using TransitionMatrices; DocMeta.setdocmeta!(TransitionMatrices, :DocTestSetup, :(using TransitionMatrices); recursive=true); doctest(TransitionMatrices)'

[group("checks")]
ci:
    @just test
    @just doctest
    @just docs

[group("quality")]
format-julia:
    julia -e 'using Pkg; Pkg.add("JuliaFormatter"); using JuliaFormatter; format(".")'
alias fj := format-julia

[group("quality")]
format-julia-check:
    julia -e 'using Pkg; Pkg.add("JuliaFormatter"); using JuliaFormatter; format(".", overwrite=false) || exit(1)'

[group("development")]
clean:
    julia -e 'for path in ("docs/build", "lcov", "coverage-lcov.info", "lcov.info"); if isdir(path); rm(path; recursive=true, force=true); elseif isfile(path); rm(path; force=true); end; end; for (root, dirs, files) in walkdir("."); for file in files; if endswith(file, ".jl.cov") || endswith(file, ".jl.mem") || (occursin(".jl.", file) && endswith(file, ".cov")); rm(joinpath(root, file); force=true); end; end; end'
