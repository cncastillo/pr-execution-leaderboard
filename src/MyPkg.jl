module MyPkg

using StaticArrays
using LinearAlgebra

const γ = 2 * π * 42.58 * 1e6 # rad s^-1 T^-1
const T₁ = 1.0  # s
const T₂ = 0.5  # s
const B₃ = 1e-7  # T
const B = SA[0.0; 0.0; B₃]  # T
const M₀ = 1.0
const m0 = SA[M₀; 0.0; 0.0]

struct Theoretical
end

function solve(m0, dt, tmax, method::Theoretical)
	t = collect(Float32, dt:dt:tmax)
	return SA[
		@. M₀ * cos(γ * B₃ * t) * exp(-t / T₂);
		@. -M₀ * sin(γ * B₃ * t) * exp(-t / T₂);
		@. M₀ * (1 - exp(-t / T₁))
	]'
end

function step(dt, m)
	return m .+ dt * bloch(m)
end

function solve(m0, dt, tmax, method)
	Nsteps = Int(tmax / dt)
	m = SVector{3}(m0)
	mt = zeros(Float64, (Nsteps, 3))
	for i in 1:Nsteps
		M₁, M₂, M₃ = m
		m = m .+ (dt * γ) .* SA[M₂*B₃-M₁/T₂, -M₁*B₃+M₂/T₂, (M₃-M₀)/T₁]
		mt[i, :] = m
	end
	return mt
end

struct ForwardEuler
end

"Returns the magnetization vector at the next time step using the Forward Euler method."
function step(dt, m, method::ForwardEuler)
	return m .+ dt * bloch(m)
end

"Returns the time derivative of the magnetization vector m."
function bloch(m)
	M₁, M₂, M₃ = m
	B₁, B₂, B₃ = B
	M₀ = m0[1]
	crossmB = [M₂ * B₃ - M₃ * B₂; M₃ * B₁ - M₁ * B₃; M₁ * B₂ - M₂ * B₁]
	return γ * crossmB - [M₁ / T₂; M₂ / T₂; (M₃ - M₀) / T₁]
end

export solve, step, ForwardEuler, Theoretical

end # module MyPkg
