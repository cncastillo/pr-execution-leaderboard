module MyModule

# A simple function to demonstrate benchmarking
function expensive_computation(n::Int)
    # Simulate an expensive computation (e.g., Fibonacci series)
    fib = [0, 1]
    for i in 3:n
        push!(fib, fib[end-1] + fib[end-2])
    end
    return fib[n]
end

# Another function to demonstrate benchmarking
function fast_function(x::Int)
    return x * 2  # A simple, fast function for comparison
end

end  # module MyModule
