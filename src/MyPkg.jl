module MyPkg

import LinearAlgebra

export step, solve, Theoretical, ForwardEuler

M0 = 1
γ = 2π * 42.58 * 10^6
T1 = 1
T2 = 0.5

Bz = 10^(-7)
B = [0, 0, Bz]

tmax = 3
m0 = [M0, 0, 0]

struct Theoretical
end

struct ForwardEuler
end

struct RungeKutta2
end

function bloch(M)
    DM = γ * LinearAlgebra.cross(M, B) - [M[1], M[2], 0] / T2 - [0, 0, M[3] - M0] / T1
    return DM
end

function step(dt, m, method::ForwardEuler)
    Mi = m
    Mf = Mi + dt * bloch(Mi)
    return Mf
end

"""
solve(m0, dt, tmax, method)

Solve the Bloch equations for a given initial magnetization `m0`, time step `dt`, and maximum time `tmax` using the specified method.

# Arguments
- `m0::Vector{Float64}`: Initial magnetization.
- `dt::Float64`: Time step.
- `tmax::Float64`: Maximum time.
- `method::Union{Theoretical, ForwardEuler}`: Method to use.

# Returns
- `mt::Matrix{Float64}`: Magnetization at each time step.
"""
function solve(m0, dt, tmax, method)
    Nsteps = Int(floor(tmax / dt))
    m = m0
    mt = zeros(3, Nsteps)
    for i in 1:Nsteps
        m = step(dt, m, method)
        mt[:, i] = m
    end
    return mt
end

function solve(m0, dt, tmax, method::Theoretical)
    Nsteps = Int(floor(tmax / dt))

    t = (0:Nsteps-1) * dt

    M = zeros(3, Nsteps)

    @. M[1, :] = M0 * cos(γ * Bz * t) * exp(-t / T2)
    @. M[2, :] = -M0 * sin(γ * Bz * t) * exp(-t / T2)
    @. M[3, :] = M0 * (1 - exp(-t / T1))

    return M
end

end
