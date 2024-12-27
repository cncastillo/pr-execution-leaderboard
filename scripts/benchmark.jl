using BenchmarkTools
using JSON

function my_function_to_benchmark()
    for i in 1:15
        println("hola")
    end
end

# Benchmark the function
b = @benchmark my_function_to_benchmark()

# Extract the minimum execution time (or any other relevant metric)
execution_time = minimum(b).time / 1e9  # Convert to seconds

# Save results to JSON
open("results.json", "w") do f
    JSON.print(f, Dict("benchmark_results" => [Dict("time" => execution_time)]))
end
