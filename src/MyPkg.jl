module MyPkg

γ = 2 * π * 42.58 * 1e6  # rad s^-1 T^-1
T₁ = 1.  # s
T₂ = 0.5  # s
B = [0.; 0.; 1e-7]  # T
m0 = [1.; 0.; 0.]

struct Theoretical
end

"Returns the theoretical solution of the Bloch equations for a given initial magnetization m0, time step dt, maximum time tmax."
function solve(m0, dt, tmax, method::Theoretical)
    t = dt:dt:tmax
    M₀ = m0[1]
    B₃ = B[3]
    M₁ = M₀ * cos.(γ * B₃ * t) .* exp.(-t / T₂)
    M₂ = -M₀ * sin.(γ * B₃ * t) .* exp.(-t / T₂)
    M₃ = M₀ * (1 .- exp.(-t / T₁))
    return [M₁'; M₂'; M₃']
end

"Simulates the Bloch equations for a given initial magnetization m0, time step dt, maximum time tmax and a simulation method."
function solve(m0, dt, tmax, method)
    Nsteps = Int(tmax / dt)
    m = m0
    mt = zeros(eltype(m), (length(m), Nsteps))
    for i in 1:Nsteps
        m = step(dt, m, method)
        mt[:, i] = m
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
