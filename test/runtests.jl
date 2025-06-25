using Test
using MyPkg

@testset "Forward Euler vs Theoretical Solution" begin
    m0_vec = [0.0, 0.0, 1.0]
    dt, tmax = 1e-4, 1e-3 
    
    euler_sol = solve(m0_vec, dt, tmax, ForwardEuler())
    theoretical_sol = solve(m0_vec, dt, tmax, Theoretical())
    
    @test isapprox(euler_sol, theoretical_sol, atol=0.01)

end