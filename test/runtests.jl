using Test
using MyPkg
using StaticArrays

const gamma = 2pi * 42.58 * 1e6
const M0 = 1.0
const T1 = 1.0
const T2 = 0.5
const Bz = 1e-7

const Ti = SA[1/T2, 1/T2, 1/T1]
const m0t1 = SA[0.0,0.0, M0/T1]
const B = SA[0.0, 0.0, Bz]
const tmax = 3.0
const m0 = SA[M0, 0.0, 0.0]

dt = 1e-1

@test abs(solve(m0, dt, 3.0, ForwardEuler())[1, end] - solve(M0, dt, 3.0, Theoretical())[1, end]) <= 1e-3