using Test
using MyPkg


#Initialization of variables and constants
γ = 2.675e8         #Frequency of hydrogen (rad/(s*T))
T1 = 1          #Longitudinal relaxation time (s)
T2 = 1             #Transverse relaxation time (s)
M_0 = 1             #Equilibrium magnetization in z axis
#Initial values
m0 = [1.0, 0.0, 0.0]            #Initial magnetization
b = [0.0, 0.0, 1e-7]            #Magnetic field B=(0,0,B_z)
Δt = 0.001                      #Simulation step size for time
t_max = 5.0                    #Simulation maximum Time

#Expected results
t = 0:Δt:t_max
Mx = @. m0[1] * exp(-t / T2) * cos(γ * b[3] * t) + m0[2] * exp(-t / T2) * sin(γ * b[3] * t)
My = @. m0[2] * exp(-t / T2) * cos(γ * b[3] * t) - m0[1] * exp(-t / T2) * sin(γ * b[3] * t)
Mz = @. M_0 + (m0[3] - M_0) * exp(-t / T1)
M = hcat(Mx, My, Mz)'

expected_result = M

numerical_solution1 = solver(m0, b, Δt, t_max, FEuler())
numerical_solution2 = solver(m0, b, Δt, t_max, RK2())

@test numerical_solution1 ≈ expected_result atol=8
@test numerical_solution2 ≈ expected_result atol=1e-1