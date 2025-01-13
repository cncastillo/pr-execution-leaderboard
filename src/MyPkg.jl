module MyPkg

struct Theoretical end
struct ForwardEuler end
struct ForwardEulerOptim end

gamma = 2 * π * 42.58 * 10^6
T1 = 1
T2 = 0.5
m0 = 1
Bz = 10^-7

function crossBz(a)
  return [a[2] * Bz, -a[1] * Bz, 0]
end

function bloch(m)
  v1 = gamma .* crossBz(m)
  v2 = [m[1] / T2, m[2] / T2, (m[3] - m0) / T1]
  return v1 .- v2
end

function step(dt, m, method::ForwardEulerOptim)
  return m .+ dt .* bloch(m)
end


"""Solve with Theoretical solution."""
function solve(m0, dt, tmax, method::Theoretical)
  t = 0:dt:tmax
  x = cos.(gamma .* Bz .* t) .* exp.(-t ./ T2)
  y = -sin.(gamma .* Bz .* t) .* exp.(-t ./ T2)
  z = 1 .- exp.(-t ./ T1)
  return m0 .* [x'; y'; z']
end

"Solve with a numerical method."
function solve(m0, dt, tmax, method)
  n_steps = length(0:dt:tmax)
  m = m0
  mt = zeros(3, n_steps)
  mt[:, 1] .= m
  for i in 2:n_steps
    m = step(dt, m, method)
    mt[:, i] .= m
  end
  return mt
end

export solve, step, Theoretical, ForwardEuler, ForwardEulerOptim, m0

end  # module MyPkg
