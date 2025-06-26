module MyPkg

using LinearAlgebra: cross
using Cthulhu
using StaticArrays

const γ = 2π * 42.58e6
const M0 = 1.0
const T1 = 1.0
const T2 = 100.0
const B = 1e-6

function step(dt, m, method) 
    # This is a generic function that will be dispatched based on method type
    error("Method not implemented")
end

struct ForwardEuler
end

function bloch(m)
    cross_term = γ * cross(m, SA[0, 0, B])
    relaxation_term = [m[1] / T2, m[2] / T2, (m[3] - M0) / T1]
    dM = cross_term - relaxation_term
    return dM
end

function step(dt, m, ::ForwardEuler)
    dM = bloch(m)
    m_next = m + dM .* dt
    return m_next
end


# Task 4
"Solves the Bloch equation using the specified method."
function solve(m0, dt, tmax, method)
    Nsteps = Int(round(tmax / dt))
    mt = Vector{SVector{3, Float64}}(undef, Nsteps + 1)
    mt[1] = SVector(m0[1], m0[2], m0[3])
    for i in 1:Nsteps
        mt[i + 1] = step(dt, mt[i], method)
    end
    return mt
end

struct Theoretical
end
function solve(m0, dt, tmax, method::Theoretical)
    ts = 0:dt:tmax
    
    # Extract initial magnetization components
    Mx0, My0, Mz0 = m0[1], m0[2], m0[3]
    
    γB = γ * B
    tst2 = @. ts / T2
    tst1 = @. ts / T1
    Mx = @. Mx0 * cos(γB * ts) * exp(-tst2) - My0 * sin(γB * ts) * exp(-tst2)
    My = @. Mx0 * sin(γB * ts) * exp(-tst2) + My0 * cos(γB * ts) * exp(-tst2)
    Mz = @. Mz0 * exp(-tst1) + M0 * (1 - exp(-tst1)) 

    return [Mx, My, Mz]
end

export step, ForwardEuler, RungeKutta2, solve, bloch, bloch2, Theoretical
end # module MyPkg
