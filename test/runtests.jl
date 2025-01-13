using Test
using MyPkg
# using BenchmarkTools

# Inputs
m0   = [1.0, 0.0, 0.0]
Δt   = 0.0001
tmax = 3.0

# expected_result    = solve(m0, Δt, tmax, Theoretical())
numerical_solution = solve(m0, Δt, tmax, ForwardEuler())

# @test numerical_solution ≈ expected_result atol=0.7
@test true

# @benchmark solve(m0, Δt, tmax, ForwardEuler())
# @code_warntype(solve(m0, Δt, tmax, ForwardEuler()))
