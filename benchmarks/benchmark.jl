using BenchmarkTools
push!(LOAD_PATH, "src")
using MyModule

# Get GitHub username from environment variable
github_username = get(ENV, "GITHUB_ACTOR", "UNKNOWN")

suite = BenchmarkGroup()
suite[github_username] = BenchmarkGroup(["username"])  
suite[github_username][10] = @benchmarkable expensive_computation(10)
tune!(suite)
results = run(suite, verbose=true)

# Save the benchmark results to a JSON file
BenchmarkTools.save("benchmarks/benchmark_results.json", median(results))
