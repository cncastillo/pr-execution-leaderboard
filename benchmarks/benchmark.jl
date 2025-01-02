using BenchmarkTools
using MyPkg

# Get GitHub username from environment variable
github_username = get(ENV, "GITHUB_ACTOR", "UNKNOWN")

suite = BenchmarkGroup()
suite["function1"][github_username] = @benchmarkable expensive_computation(300000)
tune!(suite)
results = run(suite, verbose=true)

# Extract median and standard deviation
results_summary = Dict()
for (name, result) in results
    results_summary[name] = Dict(
        "median" => median(result),
        "std" => std(result)
    )
end

# Save the benchmark results to a JSON file
using JSON
open("benchmarks/benchmark_results.json", "w") do io
    JSON.print(io, results_summary)
end