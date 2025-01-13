using Test
using MyPkg
using Cthulhu

# Inputs
M0 = 1.0
tmax = 3.0
m0 = [M0, 0.0, 0.0]
Δt = 0.001

expected_result    = solve(m0, Δt, tmax, Theoretical())
@time numerical_solution = solve(m0, Δt, tmax, ForwardEuler());

# @test numerical_solution ≈ expected_result atol=1
