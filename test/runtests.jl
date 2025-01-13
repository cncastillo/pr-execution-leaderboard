using Test
using MyPkg

# Inputs
M0 = 1
tmax = 3
m0 = [M0, 0, 0]
Δt = 1/1000000

expected_result    = solve(m0, Δt, tmax, Theoretical())
numerical_solution = solve(m0, Δt, tmax, ForwardEuler())

@test numerical_solution ≈ expected_result atol=0.1
