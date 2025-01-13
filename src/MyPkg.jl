module MyPkg
using LinearAlgebra: cross
using StaticArrays

export solve, step, Theoretical, ForwardEuler

const gamma = 2pi * 42.58 * 1e6
const M0 = 1.0
const T1 = 1.0
const T2 = 0.5
const Bz = 1e-7
const gammaBz = gamma * Bz

const M0T1 = -M0/T1

const Ti = SA[1/T2, 1/T2, 1/T1]
const m0t1 = SA[0.0,0.0, M0/T1]
const B = SA[0.0, 0.0, Bz]
const tmax = 3.0
const m0 = SA[M0, 0.0, 0.0]

struct ForwardEuler

end

struct Theoretical

end

function step(dt, m, gammaBzdt, M0T1dt, Tidt, method::ForwardEuler)
    return m .+ bloch(m, gammaBzdt, M0T1dt, Tidt)
end

function aux(a, gammaBzdt, M0T1dt)
    return @SVector [a[2] * gammaBzdt, -a[1] * gammaBzdt, M0T1dt]
  end

function bloch(m, gammaBzdt, M0T1dt, Tidt)
    return aux(m, gammaBzdt, M0T1dt) .- (Tidt .* m) 
end

function zeros_via_calloc(::Type{T}, dims::Integer...) where T
    ptr = Ptr{T}(Libc.calloc(prod(dims), sizeof(T)))
    return unsafe_wrap(Array{T}, ptr, dims; own=true)
 end


function solve(m0, dt, tmax, method)
    Nsteps = Int(ceil(tmax/dt))
    m = SVector{3}(m0)
    mt = zeros_via_calloc(Float64, ceil(Int64,tmax/dt) + 1, 3)#zeros(Float64, (ceil(Int64,tmax/dt) + 1, 3))
    gammaBzdt, M0T1dt, Tidt = gammaBz * dt, M0T1 * dt, Ti .* dt
    @inbounds for i in 1:Nsteps
        m = step(dt, m,gammaBzdt, M0T1dt, Tidt, method)
        mt[i, :] .= m
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
