module MyPkg

using StaticArrays

struct Theoretical end
struct ForwardEuler end
struct ForwardEulerOptim end

const gamma = 2 * π * 42.58 * 10^6
const T1 = 1
const T2 = 0.5
const m0 = 1
const Bz = 10^-7
const gammaBz = gamma * Bz

const V1 = SVector(1 / T2, 1 / T2, 1 / T1)
const V2 = SVector(0, 0, m0 / T1)

function crossBz(a)
  return @SVector [a[2] * gammaBz, -a[1] * gammaBz, 0.0]
end

function bloch(m)
  return crossBz(m) .- m .* V1 .+ V2
end

function step(dt, m, method::ForwardEulerOptim)
  return m .+ dt .* bloch(m)
end


"""Solve with Theoretical solution."""
function solve(m0, dt, tmax, method::Theoretical)
  t = 0:dt:tmax
  x = cos.(gammaBz .* t) .* exp.(-t ./ T2)
  y = -sin.(gammaBz .* t) .* exp.(-t ./ T2)
  z = 1 .- exp.(-t ./ T1)
  return m0 .* [x'; y'; z']
end

"Solve with a numerical method."
function solve(m, dt, tmax, method)
  n_steps = Int64(ceil(tmax / dt) + 1)
  m = SVector{3,Float64}(m)
  mt = zeros(3, n_steps)
  for i in axes(mt, 2)
    m = step(dt, m, method)
    mt[:, i] .= m
  end
  return mt
end

export solve, step, Theoretical, ForwardEuler, ForwardEulerOptim, m0

end  # module MyPkg
