module MyPkg
using LinearAlgebra: cross
using StaticArrays

export solve, step, Theoretical, ForwardEuler

const gamma = 2pi * 42.58 * 1e6
const M0 = 1.0
const T1 = 1.0
const T2 = 0.5
const Bz = 1e-7

const Ti = SA[1/T2, 1/T2, 1/T1]
const m0t1 = SA[0.0,0.0, M0/T1]
const B = SA[0.0, 0.0, Bz]
const tmax = 3.0
const m0 = SA[M0, 0.0, 0.0]

struct ForwardEuler

end

struct Theoretical

end

function step(dt, m, method::ForwardEuler)
    return m .+ dt .* bloch(m)
end

function bloch(m)
    crossproduct = cross(m, B)
    return gamma .* crossproduct .- (Ti .* m) .- m0t1
end

"""
    solve(m0::AbstractVector, dt::Real, tmax::Real, method::Function) -> AbstractMatrix

Solves the time evolution of a system using a specified method.

# Arguments
- `m0::AbstractVector`: Initial state vector of the system.
- `dt::Real`: Time step for the simulation.
- `tmax::Real`: Total simulation time.
- `method::Singleton`: Struct that determines the time evolution step.

# Returns
- `AbstractMatrix`: Matrix where each column represents the state vector at each time step.
"""
function solve(m0, dt, tmax, method)
    Nsteps = Int(ceil(tmax/dt))
    m = SVector{3}(m0)
    mt = zeros(3, Int(ceil(tmax/dt)) + 1)
    for i in 1:Nsteps
        m = step(dt, m, method)
        mt[:, i] = m
    end
    return mt
end

function solve(M0, dt, tmax, method::Theoretical)
    mt = zeros(3, Int(ceil(tmax/dt)) + 1)
    t = 0:dt:tmax
    mt[1, :] .= M0 .* cos.(gamma .* Bz .* t) .* exp.(-t ./ T2)
    mt[2, :] .= -M0 .* sin.(gamma .* Bz .* t) .* exp.(-t ./ T2)
    mt[3, :] .= M0 .* (1 .- exp.(-t ./ T1))
    return mt
end

end
