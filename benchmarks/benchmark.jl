using BenchmarkTools
using MyPkg
using JSON
using Pkg

# Get GitHub username from environment variable
github_username = get(ENV, "GITHUB_ACTOR", "UNKNOWN")

# Define the benchmark
suite = BenchmarkGroup()
suite["function1"][github_username] = @benchmarkable expensive_computation(30000)
suite["function2"][github_username] = @benchmarkable fast_function()

# Tune and run the benchmark
tune!(suite)
results = run(suite, verbose=true)
median_results = median(results)
median_results["function1"]["UNKNOWN"].gctime = std(results["function1"][github_username]).time
median_results["function2"]["UNKNOWN"].gctime = std(results["function2"][github_username]).time
BenchmarkTools.save("benchmarks/benchmark_results.json", median_results)
