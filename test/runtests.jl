using Test
using MyPkg

# Inputs
m0   = [1.0, 0.0, 0.0]
Δt   = 0.000001
tmax = 3.0

expected_result    = solve(m0, Δt, tmax, Theoretical())
numerical_solution = solve(m0, Δt, tmax, ForwardEuler())

@test sum((expected_result-numerical_solution).^2)/length(expected_result)<0.5