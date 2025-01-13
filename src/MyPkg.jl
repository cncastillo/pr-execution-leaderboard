module MyPkg
using LinearAlgebra

export solve, step, Theoretical, ForwardEuler

const gamma = 2 * π * 42.58 * 1e6 # rad/(s T)
const M_0 = 1
const T_1 = 1
const T_2 = 0.5
const B_z = 1e-7
const B = [0.0, 0.0, B_z]

struct RungeKutta2 end
struct ForwardEuler end
struct Theoretical end

function solve(m0, dt, tmax, method)
    Nsteps = Int(ceil(tmax/dt) + 1)
    mt = zeros(3, Nsteps)
    m = m0
    for i in 1:Nsteps
        m = step(dt, m, method)
        mt[:, i] = m
    end 
    return mt
end

function bloch(m)
    mxb = LinearAlgebra.cross(m, B)
    mx, my, mz = m
    dmx = gamma * mxb[1] - mx/T_2
    dmy = gamma * mxb[2] - my/T_2
    dmz = gamma * mxb[3] - (mz - M_0)/T_1
    return [dmx, dmy, dmz]
end

function step(dt, m, ::ForwardEuler)
    return m .+ dt .* bloch(m)
end

# function step(dt, m, ::RungeKutta2)
#     k1 = bloch(m)
#     k2 = bloch(m .+ dt .* k1)
#     return m .+ dt .* (k1 .+ k2) ./ 2
# end

function solve(m0, dt, tmax, ::Theoretical)
    t = 0:dt:tmax
    mt = zeros(3, Int(ceil(tmax/dt) + 1))
    mt[1,:] .= m0[1] .* cos.(gamma .* B_z .* t) .* exp.(-t ./ T_2)
    mt[2,:] .= m0[1] .* sin.(gamma .* B_z .* t) .* exp.(-t ./ T_2)
    mt[3,:] .= m0[1] .* (1 .- exp.(-t ./ T_1))
    return mt
end

end # module MyPkg