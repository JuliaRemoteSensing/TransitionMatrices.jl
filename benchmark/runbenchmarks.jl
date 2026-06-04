# Local convenience runner for the benchmark suite.
#
#   julia --project=benchmark benchmark/runbenchmarks.jl [output.json]
#
# On first run it tunes the suite and caches the parameters in params.json so
# subsequent runs (and CI) are reproducible. Results are saved as JSON for
# later comparison with `BenchmarkTools.judge` / `ratio`.

using BenchmarkTools

include(joinpath(@__DIR__, "benchmarks.jl"))

const PARAMS_PATH = joinpath(@__DIR__, "params.json")

if isfile(PARAMS_PATH)
    @info "Loaded tuned parameters" path=PARAMS_PATH
else
    @info "Tuning benchmark parameters (first run; this is slow)…"
    tune!(SUITE)
    BenchmarkTools.save(PARAMS_PATH, params(SUITE))
    @info "Saved tuned parameters" path=PARAMS_PATH
end

results = run(SUITE; verbose = true)

output = isempty(ARGS) ? joinpath(@__DIR__, "results.json") : ARGS[1]
BenchmarkTools.save(output, results)
@info "Saved benchmark results" path=output

println("\n==== median timings ====")
display(median(results))
println()
