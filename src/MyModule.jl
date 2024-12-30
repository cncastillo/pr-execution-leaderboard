module MyModule

export expensive_computation, fast_function
# A simple function to demonstrate benchmarking
function expensive_computation()
    n = 30
    if n == 1
        return 0
    elseif n == 2
        return 1
    end

    a, b = 0, 1
    for i in 3:n
        a, b = b, a + b
    end
    return b
end

# Another function to demonstrate benchmarking
function fast_function()
    x = 2
    return x * 2  # A simple, fast function for comparison
end

end  # module MyModule
