module MyPkg
using LinearAlgebra

struct Theoretical
end
struct ForwardEuler
end
struct RungeKutta2
end


function step(dt, m, method::ForwardEuler)
    return m + bloch(m)*dt
end


function step(dt, m, method::RungeKutta2)
    m1 = m + bloch(m)*dt
    return m + dt*(bloch(m1)+bloch(m))/2
end

function solve(m0, dt, tmax, method)
    Nsteps = Int(round(tmax/dt))
    m = m0
    mt = zeros(3, Nsteps)
    for i in 1:Nsteps
        m = step(dt, m, method)
        mt[:, i] = m
    end
    return mt
end

function solve(m0, dt, tmax, method::Theoretical)
    Bz = [0.0, 0.0, 1e-7]       # T
    γ = 2π*42.58e6              # rad/(sT)
    T1 = 1                      # s
    T2 = 0.5                    # s
    t = dt:dt:tmax
    m0 = m0[1]
    mx = @. m0*cos(γ*Bz[1].*t)*exp.(-t./T2)
    my = @. -m0*sin(γ*Bz[2]*t)*exp(-t/T2)
    mz = @. m0*(1-exp(-t/T1))
    return [mx my mz]'
end

function bloch(m)
    Bz = [0.0, 0.0, 1e-7]       # T
    γ = 2π*42.58e6              # rad/(sT)
    T1 = 1                      # s
    T2 = 0.5                    # s
    M0 = [1.0, 0.0, 0.0]

    MB = cross(γ*m,Bz)
    Mx = MB[1] - m[1]/T2
    My = MB[2] - m[2]/T2
    Mz = MB[3] - (m[3] - M0[1])/T1

    return [Mx, My, Mz]
end

export solve, step, Theoretical, ForwardEuler, RungeKutta2
end  # module MyPkg
