using Test
using MyPkg

expected_result1    = 34 #solve(..., Theoretical())
numerical_solution1 = expensive_computation(10) #solve(..., ForwardEuler())
expected_result2    = 4 #solve(..., Theoretical())
numerical_solution2 = fast_function() #solve(..., ForwardEuler())

@test numerical_solution1 ≈ expected_result1


@test numerical_solution2 ≈ expected_result2