using Test
using MyPkg
#Initial values
m0 = [1.0, 0.0, 0.0]            #Initial magnetization
b = [0.0, 0.0, 1e-7]            #Magnetic field B=(0,0,B_z)
Δt = 0.001                      #Simulation step size for time
t_max = 5.0                    #Simulation maximum Time

expected_result1    = solver(m0, b, Δt, t_max, Reference()) #solve(..., Theoretical())
expected_result2    = solver(m0, b, Δt, t_max, FEuler()) #solve(..., Theoretical())
numerical_solution1 = solver(m0, b, Δt, t_max, FEuler()) #solve(..., ForwardEuler())
numerical_solution2 = solver(m0, b, Δt, t_max, Reference()) #solve(..., ForwardEuler())

@test numerical_solution1 ≈ expected_result2

@test numerical_solution2 ≈ expected_result1