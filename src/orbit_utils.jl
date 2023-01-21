""" `true_anomaly_from_eci(x, e, GM)`
Find the true anomaly given an ECI state vector `x`, the eccentricity `e`, and the gravitational constant `GM`.
"""
function true_anomaly_from_eci(x, e, GM)
    r = x[1:3]
    v = x[4:6]

    a = 1.0 / (2.0 / norm(r) - norm(v)^2 / GM)
    n = sqrt(GM / (a^3)) # mean motion

    E = atan(dot(r, v) / (n * a^2), (1 - norm(r) / a)) # eccentric anomaly
    theta = atan(sqrt(1 - e^2) * sin(E), cos(E) - e) # true anomaly

    return theta
end

function E_from_M(M, e; tol=1e-10)
    f(E) = E - e * sin(E) - M # = 0
    f_prime(E) = 1 - e * cos(E)

    # solve for E using Newton's method
    E_n = M
    while abs(f(E_n)) > tol
        E_n = E_n - f(E_n) / f_prime(E_n)
    end

    return E_n
end

function theta_from_E(E, e)
    return atan(sqrt(1 - e^2) * sin(E), cos(E) - e)
end

function E_from_theta(theta, e)
    return atan(sqrt(1 - e^2) * sin(theta), cos(theta) + e)
end

function M_from_E(E, e)
    return E - e * sin(E)
end

"""
Source: Followed https://github.com/sisl/SatelliteDynamics.jl/blob/master/src/astrodynamics.jl#L248
"""
function state_CART_to_OSC(x_st, GM)

    r = x_st[1:3]
    v = x_st[4:6]

    h = cross(r, v)
    W = h / norm(h)

    i = atan(sqrt(W[1]^2 + W[2]^2), W[3])
    Omega = atan(W[1], -W[2])

    p = norm(h)^2 / GM
    a = 1.0 / (2.0 / norm(r) - norm(v)^2 / GM)

    # numerical stability hack for circular/near circular orbits
    # ensures that (1-p/a) is always positive
    if p > a
        if !isapprox(a, p, atol=1e-9, rtol=1e-8)
            @warn "Approximating $p = $a"
        end
        p = a
    end
    e = sqrt(1 - p / a) # eccentricity

    n = sqrt(GM / (a^3)) # mean motion
    E = atan(dot(r, v) / (n * a^2), (1 - norm(r) / a)) # eccentric anomaly
    M = M_from_E(E, e)
    u = atan(r[3], -r[1] * W[2] + r[2] * W[1]) # mean longitude
    theta = atan(sqrt(1 - e^2) * sin(E), cos(E) - e) # true anomaly
    omega = u - theta # argument of perigee

    omega = rem2pi(omega, RoundDown)
    Omega = rem2pi(Omega, RoundDown)
    theta = rem2pi(theta, RoundDown)

    return [a, e, i, Omega, omega, M]
end

