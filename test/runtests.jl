using Test
using MyPkg

# Inputs
m0   = [1.0, 0.0, 0.0]
Δt   = 0.000001
tmax = 3.0


expected_result    = solve(m0, Δt, tmax, Theoretical())
numerical_solution = solve(m0, Δt, tmax, ForwardEuler())

@test numerical_solution ≈ expected_result atol=1e-1