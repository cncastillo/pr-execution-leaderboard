module MyPkg
    using LinearAlgebra
    using StaticArrays
    struct Theoretical
    end
    struct ForwardEuler
    end
    const GAMMA = 2 * π * 42.58 * 1e6 # rad / sT
    const M0 = 1.0
    const T1 = 1.0 # s
    const T2 = 0.5 # s
    const Bz = 1e-7 # T (tesla)
    const B = SVector{3}([0., 0., Bz])
    const xhat = SVector{3}([1., 0., 0.])
    const yhat = SVector{3}([0., 1., 0.])
    const zhat = SVector{3}([0., 0., 1.])

    """Return the magnetization vector for each timestep as a matrix.

    Progress the state of the magnetization m from its initial state m0
    according to the type of the input argument `method`.

    Included method types are Theoretical and ForwardEuler.
    """
    function solve(m0, dt, tmax, method)
        Nsteps = floor(Int64, tmax / dt)
        m = m0
        mt = zeros(3, Nsteps)
        # Note: missing m0 in the mt matrix!!
        range = 1:Nsteps
        for i in range
            m = step(dt, m, method)
            mt[:, i] = m
        end
        return mt
    end

    function theoretical_m(t)
        mx = M0 * cos(GAMMA * Bz * t) * exp(-t / T2)
        my = -M0 * sin(GAMMA * Bz * t) * exp(-t / T2)
        mz = M0 * (1 - exp(-t / T1))
        return [mx, my, mz]
    end

    """Return the magnetization vector for each timestep as a matrix.

    Progress the state of the magnetization m from its initial state m0
    according to the Theoretical solution using this Package's initial
    conditions.
    """
    function solve(m0, dt, tmax, method::Theoretical)
        Nsteps = floor(Int64, tmax / dt)
        times = [dt * i for i in 1:Nsteps]'
        broadcasted = theoretical_m.(times)
        # The broadcasted result is a vector of vectors
        # We have to horizontally concatenate these to make a matrix 3xNsteps matrix
        mt = hcat(broadcasted...)
        # Note: missing m0 in the mt matrix!! This is consistent with the other one
        return mt
    end

    function bloch(m)
        mx, my, mz = m

        first = GAMMA * cross(m, B)
        second = (mx * xhat + my * yhat) / T2
        third = (mz - M0) / T1 * zhat

        return first - second - third
    end

    function step(dt, m, method::ForwardEuler)
        nextm = m + dt * bloch(m)
        return nextm
    end

    export solve, step, Theoretical, ForwardEuler
end
