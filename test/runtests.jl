using Test
using MyPkg

# Inputs
m0   = [1.0, 0.0, 0.0]
Δt   = 1e-5
tmax = 3.0

expected_result    = solve(m0, Δt, tmax, Theoretical())
numerical_solution = solve(m0, Δt, tmax, ForwardEuler())
maxError           = maximum(abs.(numerical_solution .- expected_result))

@test maxError < 1e-3
