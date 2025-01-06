using BenchmarkTools
using MyPkg
using JSON
using Pkg

# Get GitHub username from environment variable
github_username = get(ENV, "GITHUB_ACTOR", "UNKNOWN")
#Initial values
m0 = [1.0, 0.0, 0.0]            #Initial magnetization
b = [0.0, 0.0, 1e-7]            #Magnetic field B=(0,0,B_z)
Δt = 0.001                      #Simulation step size for time
t_max = 5.0                    #Simulation maximum time

# Define the benchmark
suite = BenchmarkGroup()
suite["function1"][github_username] = @benchmarkable solver(m0, b, Δt, t_max, FEuler())
suite["function2"][github_username] = @benchmarkable solver(m0, b, Δt, t_max, Reference())

# Tune and run the benchmark
tune!(suite)
results = run(suite, verbose=true)
median_results = median(results)
#median_results["function1"][github_username].gctime = std(results["function1"][github_username]).time
median_results["function2"][github_username].gctime = std(results["function2"][github_username]).time
BenchmarkTools.save("benchmarks/benchmark_results.json", median_results)
