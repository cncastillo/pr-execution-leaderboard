module MyPkg

export FEuler, RK2, solver, SimulationMethod

# A simple function to demonstrate benchmarking
#Definition of solvers
abstract type SimulationMethod end
struct FEuler <: SimulationMethod end
struct RK2 <: SimulationMethod end

#Initialization of variables and constants
γ = 2.675e8         #Frequency of hydrogen (rad/(s*T))
T1 = 1          #Longitudinal relaxation time (s)
T2 = 1             #Transverse relaxation time (s)
M_0 = 1             #Equilibrium magnetization in z axis

#Definition of Bloch model
function bloch_model!(dm, m, b)
    dm[1] = γ * (m[2] * b[3] - m[3] * b[2]) - (m[1] / T2)
    dm[2] = γ * (m[3] * b[1] - m[1] * b[3]) - (m[2] / T2)
    dm[3] = γ * (m[1] * b[2] - m[2] * b[1]) - (m[3] - M_0) / T1
end

#Definition of forward euler method
function step!(Δt, m_new, m_old, b, method::FEuler)
    dm = [0.0, 0.0, 0.0]
    bloch_model!(dm, m_old, b)
    m_new .= m_old .+ Δt .* dm
    return nothing
end

#Definition of runge kutta 2 method (trapezoidal)
function step!(Δt, m_new, m_old, b, method::RK2)
    #Prediction step
    dm = [0.0, 0.0, 0.0]
    bloch_model!(dm, m_old, b)
    m_pred = m_old .+ Δt .* dm
    #Correction step
    dm_pred = [0.0, 0.0, 0.0]
    bloch_model!(dm_pred, m_pred, b)
    m_new .= m_old .+ (Δt / 2) .* (dm + dm_pred)
    return nothing
end

function solver(m0, b, Δt, t_max, method::SimulationMethod)
    #Calculate the number of steps in the simulation using t_max and Δt
    N_steps = Int(t_max / Δt)           #Set the number of simulation steps of time
    m = zeros(3, N_steps + 1)           #Instanciate result array
    m[:, 1] .= m0                       #Set the initial value
    #Iterate to generate results
    @assert N_steps > 1
    for i in 1:N_steps
        @views step!(Δt, m[:, i+1], m[:, i], b, method)
    end
    return m
end

end  # module MyModule
