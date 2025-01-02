using BenchmarkTools
using MyPkg
using JSON
using Pkg

# Get GitHub username from environment variable
github_username = get(ENV, "GITHUB_ACTOR", "UNKNOWN")

# Define the benchmark
suite = BenchmarkGroup()
suite[github_username] = @benchmarkable expensive_computation(30000)

# Tune and run the benchmark
tune!(suite)
results = run(suite, verbose=true)

# Prepare the benchmark data
benchmark_data = [
    Dict("Julia" => string(VERSION), "BenchmarkTools" => string(Pkg.installed()["BenchmarkTools"])),
    [
        [
            "BenchmarkGroup",
            Dict(
                "data" => Dict(
                    github_username => [
                        "TrialEstimate",
                        Dict(
                            "allocs" => median(results[github_username]).allocs,
                            "time" => median(results[github_username]).time,
                            "std" => std(results[github_username]).time,
                            "memory" => median(results[github_username]).memory,
                            "params" => median(results[github_username]).params,
                            "gctime" => median(results[github_username]).gctime
                        )
                    ]
                ),
                "tags" => []
            )
        ]
    ]
]

# Save the benchmark results to a JSON file
open("benchmarks/benchmark_results.json", "w") do io
    write(io, JSON.json(benchmark_data))
end

