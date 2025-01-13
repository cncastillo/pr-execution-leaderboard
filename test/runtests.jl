using Test
using MyPkg

# Inputs
m0   = [1.0, 0.0, 0.0]
dt = 1e-3
tmax = 3.0

expected_result    = solve(m0, Δt, tmax, Theoretical())
numerical_solution = solve(m0, Δt, tmax, ForwardEuler())

@test abs(solve(m0, dt, 3.0, ForwardEuler())[1, end] - solve(M0, dt, 3.0, Theoretical())[end, 1]) <= 1e-3
