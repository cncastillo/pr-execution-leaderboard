module MyPkg

export solve, step, Theoretical, ForwardEuler


struct Theoretical
end
struct ForwardEuler
end


function solve(m0, dt, tmax, method)
    Nsteps = Int(tmax/dt)
    m = m0
    mt = zeros(Float64,3,Nsteps)
    for i in 1:Nsteps
        m = step(dt, m, method)
        mt[:,i] = m
    end
    return mt
end

function solve(m0, dt, tmax, method::Theoretical)
    t = (dt:dt:tmax)
    g::Float64 = 2 * π * 42.58 * 10^6
    Bz::Float64 = 10^(-7)
    T1::Float64 = 1
    T2::Float64 = 0.5
    M0::Float64 = 1
    mt = hcat((x -> [m0[1]*cos(g*Bz*x)*exp(-x/T2); -m0[1]*sin(g*Bz*x)*exp(-x/T2); m0[1]*(1-exp(-x/T1))]).(t))
    result = zeros(Float64, length(mt), 3)
    for i in 1:length(mt)
        result[i,:] = mt[i, 1]
    end
    return result'
end

function bloch(M)
    g::Float64 = 2 * π * 42.58 * 10^6
    Bz::Float64 = 10^(-7)
    T1::Float64 = 1
    T2::Float64 = 0.5
    M0::Float64 = 1
    return [-1/T2 g*Bz 0; -g*Bz -1/T2 0; 0 0 -1/T1] * M + [0 ; 0 ; M0/T1]
end

step(dt, m, method::ForwardEuler) = m + dt*bloch(m)

end # module MyPkg