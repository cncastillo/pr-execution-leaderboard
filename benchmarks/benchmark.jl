using BenchmarkTools
push!(LOAD_PATH, "src")
using MyModule

# Get GitHub username from environment variable
github_username = get(ENV, "GITHUB_ACTOR", "UNKNOWN")

suite = BenchmarkGroup()
suite["function"] = BenchmarkGroup(["USERNAME"])  
suite["function1"][github_username] = @benchmarkable expensive_computation()
tune!(suite)
results = run(suite, verbose=true)

# Save the benchmark results to a JSON file
BenchmarkTools.save("benchmarks/benchmark_results.json", median(results))
