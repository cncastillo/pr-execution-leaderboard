using Test
using MyPkg
using BenchmarkTools

# @testset "Forward Euler vs Theoretical Solution" begin
#     m0_vec = [0.0, 0.0, 1.0]
#     dt, tmax = 1e-4, 1e-3 
    
#     euler_sol = solve(m0_vec, dt, tmax, ForwardEuler())
#     theoretical_sol = solve(m0_vec, dt, tmax, Theoretical())
    
#     @test isapprox(euler_sol, theoretical_sol, atol=0.01)
# end

@testset "Performance Timing Comparison" begin
    println("\n" * "="^50)
    println("PERFORMANCE TIMING COMPARISON")
    println("="^50)
    
    m_int = [1, 0, 0]          
    m_float = [1.0, 0.0, 0.0]
    dt, tmax = 1e-6, 1e-4     
    
    println("\nBenchmark results:")
    print("bloch (int):   "); @btime bloch($m_int)
    print("bloch (float): "); @btime bloch($m_float)
    print("bloch2 (float):  "); @btime bloch2($m_float)
    print("solve: "); @btime solve($m_float, $dt, $tmax, ForwardEuler())
    print("solve (theoretical): "); @btime solve($m_float, $dt, $tmax, Theoretical())
    print("step (ForwardEuler): "); @btime MyPkg.step($dt, $m_float, ForwardEuler())
    
    println("\n" * "="^50)
end