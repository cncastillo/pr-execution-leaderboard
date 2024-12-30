using BenchmarkTools
push!(LOAD_PATH, "src")
using MyModule  # Import your module

# Benchmark the expensive_computation function
results = @benchmark expensive_computation(20)  # Benchmark for n = 30

# Save the benchmark results to a JSON file
BenchmarkTools.save("benchmarks/benchmark_results.json", results)