module MyPkg
using StaticArrays

struct Theoretical
end
struct ForwardEuler
end
struct RungeKutta2
end

const Bz = SVector{3, Float64}(0.0, 0.0, 1e-7)
const γ = 2π*42.58e6::Float64
const Bzz_γ = Bz[3]*γ::Float64
const T1 = 1.00::Float64 
const T1_i = 1/T1::Float64
const T2 =  0.5::Float64
const T2_i = 1/T2::Float64
const M0 = 1.0::Float64
const bloch_matrix = @SMatrix[  -T1_i      Bzz_γ    0;
                                -Bzz_γ    -T2_i     0;
                                0           0       -T1_i]
const independent = @SVector[0, 0, M0*T1_i]

function step!(m, Δt::Float64,  method::ForwardEuler)
    m + bloch(m)*Δt
end


function step!(m, Δt::Float64, method::RungeKutta2)
    m1 = m + bloch(m)*Δt
    m + Δt*(bloch(m1)+bloch(m))/2
end

function solve(m0, Δt, tmax, method)
    Nsteps = Int(round(tmax/Δt))
    m = m0
    mt = zeros(Nsteps, 3)
    mt[1, :] = m0
    for i in axes(mt, 1)
        m = step!(m, Δt, method)
        mt[i, :] = m
    end
    return mt
end

function solve(m0, Δt, tmax, method::Theoretical)
    t = Δt:Δt:tmax
    mx = @. m0[1].*cos.(γ*Bz[1].*t)*exp.(-t./T2)
    my = @. -m0[1].*sin.(γ*Bz[2]*t)*exp.(-t./T2)
    mz = @. m0[1]*(1-exp(-t/T1))
    return [mx my mz]
end

function bloch(m)
    M_new = bloch_matrix*m+independent
    return M_new
end

export solve, step, Theoretical, ForwardEuler, RungeKutta2
end  # module MyPkg
