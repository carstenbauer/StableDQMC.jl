# Runs the StableDQMC.jl benchmarks.
# A branch name or commit id can be supplied as an argument.

using Pkg
Pkg.activate(".")

using PkgBenchmark

branch_or_commit = nothing
if length(ARGS) > 0
    branch_or_commit = ARGS[1]
end

out = isnothing(branch_or_commit) ? "master" : string(branch_or_commit)

config = BenchmarkConfig(id = out,
                         juliacmd = `julia -O3 --project=.`,
                         env = Dict("JULIA_NUM_THREADS" => 1,
                                    "OMP_NUM_THREADS" => Sys.CPU_THREADS,
                                    "MKL_NUM_THREADS" => Sys.CPU_THREADS))


r = benchmarkpkg("StableDQMC", config;
             script = "benchmark/benchmarks.jl",
             resultfile = string(out, ".bench"))

export_markdown(string(out, ".md"), r)