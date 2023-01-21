
""" `position_inertial_to_ks(x)`
Given a 3-vector `x`, correxponding to an inertial position,
return the ks-transformed 4-vector `p`.
"""
function position_inertial_to_ks(x)
    x1, _, _ = x
    r = sqrt(x'x)

    if x1 >= 0
        return position_cart_to_ks_plus(x)
    else
        return position_cart_to_ks_minus(x)
    end
end


""" position_cart_to_ks_plus
Convert a vector to ks, uses the transform best suited for vectors with x1 > 0
"""
function position_cart_to_ks_plus(x)
    x1, x2, x3 = x
    r = sqrt(x'x)

    p4 = 0.0
    p1 = sqrt(0.5 * (r + x1) - p4^2)
    p2 = (x2 * p1 + x3 * p4) / (x1 + r)
    p3 = (x3 * p1 - x2 * p4) / (x1 + r)
    return [p1, p2, p3, p4]
end

""" position_cart_to_ks_minus
Convert a vector to ks, uses the transform best suited for vectors with x1 < 0
"""
function position_cart_to_ks_minus(x)
    x1, x2, x3 = x
    r = sqrt(x'x)

    p3 = 0.0
    p2 = sqrt(0.5 * (r - x1) - p3^2)
    p1 = (x2 * p2 + x3 * p3) / (r - x1)
    p4 = (x3 * p2 - x2 * p3) / (r - x1)
    return [p1, p2, p3, p4]
end

""" position_cart_to_ks_comparision
Convert two cartesian vectors to KS vectors, making sure to use the same method to convert both so the difference between them is small.
"""
function position_cart_to_ks_comparison(x, y)
    x1, _, _ = x
    y1, _, _ = y

    if x1 >= 0
        if y1 >= 0
            return position_cart_to_ks_plus(x), position_cart_to_ks_plus(y)
        else
            if x1 + y1 >= 0
                return position_cart_to_ks_plus(x), position_cart_to_ks_plus(y)
            else
                return position_cart_to_ks_minus(x), position_cart_to_ks_minus(y)
            end
        end
    else
        if y1 < 0
            return position_cart_to_ks_minus(x), position_cart_to_ks_minus(y)
        else
            if x1 + y1 >= 0
                return position_cart_to_ks_plus(x), position_cart_to_ks_plus(y)
            else
                return position_cart_to_ks_minus(x), position_cart_to_ks_minus(y)
            end
        end
    end
end

function position_cart_to_ks_plus_near(x, p_prev)
    x1, x2, x3 = x
    r = sqrt(x'x)

    p1, _, _, p4 = p_prev

    # find a single solution for q1, q4
    q4 = 0.0
    q1 = sqrt(0.5 * (r + x1) - q4^2)

    qbar = [q1; q4]
    pbar = [p1; p4]

    # solution nearest p is pbar scaled to have norm = qbar
    alpha = sqrt((qbar'qbar) / (pbar'pbar))

    qnear = alpha * pbar

    q1, q4 = qnear

    q2 = (x2 * q1 + x3 * q4) / (x1 + r)
    q3 = (x3 * q1 - x2 * q4) / (x1 + r)
    return [q1, q2, q3, q4]
end

function position_cart_to_ks_minus_near(x, p_prev)
    x1, x2, x3 = x
    r = sqrt(x'x)

    _, p2, p3, _ = p_prev

    # find a single solution for q2, q3
    q3 = 0.0
    q2 = sqrt(0.5 * (r - x1) - q3^2)

    qbar = [q2; q3]
    pbar = [p2; p3]

    alpha = sqrt((qbar'qbar) / (pbar'pbar))

    qnear = alpha * pbar

    q2, q3 = qnear

    q1 = (x2 * q2 + x3 * q3) / (r - x1)
    q4 = (x3 * q2 - x2 * q3) / (r - x1)
    return [q1, q2, q3, q4]
end

function position_cart_to_ks_comparison_near(x, px_prev, y, py_prev)
    x1, _, _ = x
    y1, _, _ = y

    if x1 >= 0
        if y1 >= 0
            return position_cart_to_ks_plus_near(x, px_prev), position_cart_to_ks_plus_near(y, py_prev)
        else
            if x1 + y1 >= 0
                return position_cart_to_ks_plus_near(x, px_prev), position_cart_to_ks_plus_near(y, py_prev)
            else
                return position_cart_to_ks_minus_near(x, px_prev), position_cart_to_ks_minus_near(y, py_prev)
            end
        end
    else
        if y1 < 0
            return position_cart_to_ks_minus_near(x, px_prev), position_cart_to_ks_minus_near(y, py_prev)
        else
            if x1 + y1 >= 0
                return position_cart_to_ks_plus_near(x, px_prev), position_cart_to_ks_plus_near(y, py_prev)
            else
                return position_cart_to_ks_minus_near(x, px_prev), position_cart_to_ks_minus_near(y, py_prev)
            end
        end
    end
end

""" `position_ks_to_inertial()`
Given a 4-vector `p`, correxponding to a KS position,
return the inertial 3-vector `x`.
"""
function position_ks_to_inertial(p)
    p1 = p[1]
    p2 = p[2]
    p3 = p[3]
    p4 = p[4]

    x1 = p1^2 - p2^2 - p3^2 + p4^2
    x2 = 2.0 * (p1 * p2 - p3 * p4)
    x3 = 2.0 * (p1 * p3 + p2 * p4)

    return [x1, x2, x3]
end

""" `velocity_inertial_to_ks(p, x_dot)`
Given a 4-vector `p`, corresponding to the ks position, and a
3-vector `x_dot`, corresponding to inertial velocities, 
return the 4-vector of ks-transformed velocities `p_prime`.
The inertial velocities are wrt real time (t), whereas the ks-transformed velocities
are wrt fictitious time (s), where dt = ||x||ds.
"""
function velocity_inertial_to_ks(p, x_dot)

    p1 = p[1]
    p2 = p[2]
    p3 = p[3]
    p4 = p[4]

    x1_dot = x_dot[1]
    x2_dot = x_dot[2]
    x3_dot = x_dot[3]

    p1_prime = 0.5 * (p1 * x1_dot + p2 * x2_dot + p3 * x3_dot)
    p2_prime = 0.5 * (-p2 * x1_dot + p1 * x2_dot + p4 * x3_dot)
    p3_prime = 0.5 * (-p3 * x1_dot - p4 * x2_dot + p1 * x3_dot)
    p4_prime = 0.5 * (p4 * x1_dot - p3 * x2_dot + p2 * x3_dot)

    return [p1_prime, p2_prime, p3_prime, p4_prime]
end

""" `velocity_ks_to_inertial(p, p_prime)`
Given 4-vectors `p` and `p_prime`, corresponding to the ks transformed
position and velocity, respectively, return the 3-vector inertial velocity `x`.
"""
function velocity_ks_to_inertial(p, p_prime)

    p1 = p[1]
    p2 = p[2]
    p3 = p[3]
    p4 = p[4]

    p1_prime = p_prime[1]
    p2_prime = p_prime[2]
    p3_prime = p_prime[3]
    p4_prime = p_prime[4]

    pTp = p'p

    x1_dot = (2.0 / pTp) * (p1 * p1_prime - p2 * p2_prime - p3 * p3_prime + p4 * p4_prime)
    x2_dot = (2.0 / pTp) * (p2 * p1_prime + p1 * p2_prime - p4 * p3_prime - p3 * p4_prime)
    x3_dot = (2.0 / pTp) * (p3 * p1_prime + p4 * p2_prime + p1 * p3_prime + p2 * p4_prime)

    return [x1_dot, x2_dot, x3_dot]
end


acceleration_ks_to_inertial(p, u) = (2.0 / (p'p)) * inv(ks_L(p)') * u

""" `state_ks_to_inertial(p_state)`
Given an 8-vector `p_state = [p, p_prime]` compute the 6-vector cartesian state `x_state = [x, x_dot]`.
"""
function state_ks_to_inertial(p_state)
    p = p_state[1:4]
    p_prime = p_state[5:8]

    x = position_ks_to_inertial(p)
    x_dot = velocity_ks_to_inertial(p, p_prime)

    return [x; x_dot]
end

state_ks_to_cart(p_state) = state_ks_to_inertial(p_state)

""" `state_inertial_to_ks(x_state)`
Given a 6-vector `x_state = [x, x_dot]` compute the 8-vector KS state `p_state = [p, p_prime]`
"""
function state_inertial_to_ks(x_state)
    x = x_state[1:3]
    x_dot = x_state[4:6]

    p = position_inertial_to_ks(x)
    p_prime = velocity_inertial_to_ks(p, x_dot)

    return [p; p_prime]
end

""" `state_inertial_to_ks_comparison(x_state, y_state)`
Given a 6-vector `x_state = [x, x_dot]` compute the 8-vector KS state `p_state = [p, p_prime]`
"""
function state_inertial_to_ks_comparison(x_state, y_state)
    x = x_state[1:3]
    x_dot = x_state[4:6]

    y = y_state[1:3]
    y_dot = y_state[4:6]

    px, py = position_cart_to_ks_comparison(x, y)
    px_prime = velocity_inertial_to_ks(px, x_dot)
    py_prime = velocity_inertial_to_ks(py, y_dot)


    return [px; px_prime], [py; py_prime]
end


""" `ks_L(p)`
Given a ks quaternion `p`, compute the Left multiply Matrix L(p) where p*r = L(p)r
"""
function ks_L(p)

    p1 = p[1]
    p2 = p[2]
    p3 = p[3]
    p4 = p[4]

    return [p1 -p2 -p3 p4
        p2 p1 -p4 -p3
        p3 p4 p1 p2
        p4 -p3 p2 -p1]
end

""" `ks_R(p)`
Given a ks quaternion `p`, compute the Right multiply matrix R(p) where r*p = R(p)r
"""
function ks_R(p)

    p1 = p[1]
    p2 = p[2]
    p3 = p[3]
    p4 = p[4]

    return [p1 p2 p3 -p4
        p2 -p1 p4 p3
        p3 -p4 -p1 -p2
        p4 p3 -p2 p1]
end

function state_cart_to_ks_newton(x_state; p_near=[1.0; 0.0; 0.0; 0.0], tol=1e-12, max_iter=100)
    x = x_state[1:3]
    x_dot = x_state[4:6]

    p = position_cart_to_ks_newton(x; p_near=p_near, tol=tol, max_iter=max_iter)
    p_prime = velocity_inertial_to_ks(p, x_dot)

    return [p; p_prime]
end

function position_cart_to_ks_newton(x; p_near=[1.0; 0.0; 0.0; 0.0], tol=1e-12, max_iter=100, return_verbose=false)

    """
    kkt matrix
    """
    function kkt(p, λ)
        return [(I(4).+2ks_R(λ)') 2ks_L(p)'; 2ks_L(p) zeros(4, 4)]
    end

    """
    rhs vector
    """
    function rhs(x, p_near, p, λ)
        return [-p .+ p_near .- (2ks_L(p)'λ); -ks_L(p) * p .+ [x; 0]]
    end

    y = [p_near; zeros(4)]
    i = 0
    while i < max_iter
        i += 1
        F = kkt(y[1:4], y[5:8])
        b = rhs(x, p_near, y[1:4], y[5:8])
        Δy = F \ b
        y += Δy
        if norm(Δy) < tol
            break
        end
    end

    if return_verbose
        return y, i
    end

    return y[1:4]
end

function ks_h_energy(p, p_prime, GM)
    return (GM - 2 * (p_prime'p_prime)) / (p'p)
end

function ks_A_matrix(h)
    return [zeros(4, 4) I(4); -0.5*h*I(4) zeros(4, 4)]
end

function ks_B_matrix(p)
    B1 = zeros(4, 3)
    B2 = ((p'p/2)*ks_L(p)')[:, 1:3]
    return [B1; B2]
end
