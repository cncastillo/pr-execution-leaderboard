module MyPkg
using LinearAlgebra
using StaticArrays
struct Theoretical
end
struct ForwardEuler
end
const GAMMA = 2 * π * 42.58 * 1e6 # rad / sT
const M0 = 1.0
const T1 = 1.0 # s
const T2 = 0.5 # s
const Bz = 1e-7 # T (tesla)
const GAMMABz = GAMMA * Bz # T (tesla)
const NEG1_OVER_T1 = -1.0 / T1
const NEG1_OVER_T2 = -1.0 / T2
const M0_OVER_T1 = M0 / T1
# const LinearBlochPart = GAMMABz * Bcross + M2
# LinearBlochPart = [[-1.0 / T2, -GAMMABz, 0] [GAMMABz, -1.0 / T2, 0] [0, 0, -1.0 / T1]]
# NEG1_OVER_T2,           GAMMABz,                0,
#       -GAMMABz,    NEG1_OVER_T2,                0,
#            0                  0      NEG1_OVER_T1

"""Return the magnetization vector for each timestep as a matrix.

Progress the state of the magnetization m from its initial state m0
according to the type of the input argument `method`.

Included method types are Theoretical and ForwardEuler.
"""
function solve(m0, dt, tmax, method)
    Nsteps = floor(Int64, tmax / dt)
    mt = zeros(3, Nsteps + 1)
    mt[:, 1] .= m0
    for j in axes(mt, 2)
        j > 1 && @inbounds begin
            mbefore = view(mt, :, j-1)
            mt[:, j] .= step(dt, mbefore, method)
        end
    end
    return mt
end

function theoretical_m(t)
    mx = M0 * cos(GAMMA * Bz * t) * exp(-t / T2)
    my = -M0 * sin(GAMMA * Bz * t) * exp(-t / T2)
    mz = M0 * (1 - exp(-t / T1))
    return [mx, my, mz]
end

"""Return the magnetization vector for each timestep as a matrix.

Progress the state of the magnetization m from its initial state m0
according to the Theoretical solution using this Package's initial
conditions.
"""
function solve(m0, dt, tmax, method::Theoretical)
    Nsteps = floor(Int64, tmax / dt)
    times = [dt * i for i in 1:Nsteps]'
    broadcasted = theoretical_m.(times)
    # The broadcasted result is a vector of vectors
    # We have to horizontally concatenate these to make a matrix 3xNsteps matrix
    mt = hcat(broadcasted...)
    # Note: missing m0 in the mt matrix!! This is consistent with the other one
    return mt
end

function step(dt, m, method::ForwardEuler)
    @inbounds begin
        mx, my, mz = m
        return (
            mx + dt * (NEG1_OVER_T2 * mx + GAMMABz * my),
            my + dt * (NEG1_OVER_T2 * my - GAMMABz * mx),
            mz + dt * (NEG1_OVER_T1 * mz + M0_OVER_T1)
        )
    end
end

export solve, step, Theoretical, ForwardEuler
end
