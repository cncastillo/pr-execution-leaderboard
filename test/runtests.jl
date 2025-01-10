using Test
using MyPkg

m0 = 1
dt = 0.0000001
tmax = 3

euler_soln = solve(m0, dt, tmax, ForwardEuler())
theo_soln = solve(m0, dt, tmax, Theoretical())
diff = abs.(euler_soln .- theo_soln)

@test all(diff .< 1)

