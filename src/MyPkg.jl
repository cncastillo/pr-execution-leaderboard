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


function step!(m::SVector{3, Float64}, Δt::Float64,  method::ForwardEuler)
    SVector{3, Float64}(m) + SVector{3, Float64}(bloch(m)*Δt)
end


function step!(m::SVector{3, Float64}, Δt::Float64, method::RungeKutta2)
    m + bloch(m)*Δt
    m + Δt*(bloch(m1)+bloch(m))*0.5
end

function solve(m0::Vector{Float64}, Δt::Float64, tmax::Float64, method)
    Nsteps = Int(round(tmax/Δt))
    m = SVector{3, Float64}(m0)
    mt = zeros(Float64, Nsteps, 3)
    mt[1, :] = m0
    for i in axes(mt, 1)
        m = step!(m, Δt, method)
        mt[i, :] = m
    end
    return mt
end

function solve(m0::Vector{Float64}, Δt::Float64, tmax::Float64, method::Theoretical)
    t = Δt:Δt:tmax
    mx = @. m0[1].*cos.(γ*Bz[1].*t)*exp.(-t./T2)
    my = @. -m0[1].*sin.(γ*Bz[2]*t)*exp.(-t./T2)
    mz = @. m0[1]*(1-exp(-t/T1))
    return [mx my mz]
end

function bloch(m)
    #const bloch_matrix = @SMatrix[  -T1_i      Bzz_γ    0;
    #                            -Bzz_γ    -T2_i     0;
    #                            0           0       -T1_i]
    M_new = SVector{3, Float64}(-T1_i*m[1]+Bzz_γ*m[2], -Bzz_γ*m[1]-T2_i*m[2], -T1_i*m[3]) + independent
    return M_new
end

export solve, step, Theoretical, ForwardEuler, RungeKutta2
end  # module MyPkg
