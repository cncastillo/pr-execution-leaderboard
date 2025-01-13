using BenchmarkTools
using MyPkg

# Get GitHub username from environment variable
github_username = get(ENV, "GITHUB_ACTOR", "UNKNOWN")

# Inputs
m0   = [1.0, 0.0, 0.0]
Δt   = 0.001
tmax = 3.0

# Define the benchmark
suite = BenchmarkGroup()
# suite["function1"][github_username] = @benchmarkable solve(m0, Δt, tmax, ForwardEuler())
suite["function2"][github_username] = @benchmarkable solve(m0, Δt, tmax, params, ForwardEuler())

# Tune and run the benchmark
tune!(suite)
results = run(suite, verbose = true)
median_results = median(results)
# median_results["function1"][github_username].gctime = std(results["function1"][github_username]).time
median_results["function2"][github_username].gctime = std(results["function2"][github_username]).time
BenchmarkTools.save("benchmarks/benchmark_results.json", median_results)
