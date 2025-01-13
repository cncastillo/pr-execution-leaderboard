using Test
using MyPkg

m0 = [1.0, 0, 0]
dt = 0.0000001
tmax = 3

euler_soln = solve(m0, dt, tmax, ForwardEulerOptim())
theo_soln = solve(m0, dt, tmax, Theoretical())
diff = abs.(euler_soln .- theo_soln)

@test all(diff .< 10000)
