module MyPkg

using StaticArrays
# using LinearAlgebra

const γ = 2 * π * 42.58 * 1e6 # rad s^-1 T^-1
const T₁ = 1.0  # s
const T₂ = 0.5  # s
const B₃ = 1e-7  # T
const B = SA[0.0; 0.0; B₃]  # T
const M₀ = 1.0
const m0 = SA[M₀; 0.0; 0.0]


const Ts = SA[1/T₂, 1/T₂, 1/T₁]
const Mconst = SA[0.0, 0.0, M₀/T₁]

struct Theoretical
end

struct ForwardEuler
end

function solve(m0, dt, tmax, method::Theoretical)
	t = collect(Float32, dt:dt:tmax)
	return SA[
		@. M₀ * cos(γ * B₃ * t) * exp(-t / T₂);
		@. -M₀ * sin(γ * B₃ * t) * exp(-t / T₂);
		@. M₀ * (1 - exp(-t / T₁))
	]'
end

function crossMB(m)
	M₁, M₂, M₃ = m
	return SA[M₂*B₃, -M₁*B₃, 0.0]
end

function bloch(m)::SVector{3, Float32}
	return γ .* crossMB(m) .- m .* Ts + Mconst
end

function step(dt, m, method::ForwardEuler)
	return m .+ dt * bloch(m)
end

function solve(m0, dt, tmax, method)
	Nsteps = Int(tmax / dt)
	m = SVector{3}(m0)
	mt = zeros(Float64, (Nsteps, 3))
	for i in 1:Nsteps
		m = step(dt, m, method)
		mt[i, :] = m
	end
	return mt
end

export solve, ForwardEuler, Theoretical

end # module MyPkg
