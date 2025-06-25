module MyPkg
using LinearAlgebra: cross

γ = 2π * 42.58e6
m0 = 1
T1 = 1.0
T2 = 100
B = 1e-6

function step(dt, m, method) 
    # This is a generic function that will be dispatched based on method type
    error("Method not implemented")
end

struct ForwardEuler
end

function bloch(m)
    cross_term = γ * cross(m, [0, 0, B])
    relaxation_term = [m[1] / T2, m[2] / T2, (m[3] - m0) / T1]
    dM = cross_term - relaxation_term
    return dM
end

function step(dt, m, method::ForwardEuler)
    dM = bloch(m)
    m_next = m + dM * dt
    return m_next
end


# Task 4
"Solves the Bloch equation using the specified method."
function solve(m0_vec::Vector, dt::Real, tmax::Real, method)
    Nsteps = Int(round(tmax / dt))
    m = copy(m0_vec)
    mt = zeros(3, Nsteps + 1)
    mt[:, 1] = m
    for i in 1:Nsteps
        m = step(dt, m, method)
        mt[:, i + 1] = m
    end
    return mt
end

struct Theoretical
end
function solve(m0_vec::Vector, dt::Real, tmax::Real, ::Theoretical)
    ts = 0:dt:tmax
    
    # Extract initial magnetization components
    Mx0, My0, Mz0 = m0_vec[1], m0_vec[2], m0_vec[3]
    
    Mx = Mx0 .* cos.(γ * B .* ts) .* exp.(-ts ./ T2) - My0 .* sin.(γ * B .* ts) .* exp.(-ts ./ T2)
    My = Mx0 .* sin.(γ * B .* ts) .* exp.(-ts ./ T2) + My0 .* cos.(γ * B .* ts) .* exp.(-ts ./ T2)
    Mz = Mz0 .* exp.(-ts ./ T1) + m0 .* (1 .- exp.(-ts ./ T1))  # m0 is the global equilibrium magnetization = 1

    mt = [Mx'; My'; Mz']
    return mt
end

export step, ForwardEuler, RungeKutta2, solve, bloch, Theoretical
end # module MyPkg
