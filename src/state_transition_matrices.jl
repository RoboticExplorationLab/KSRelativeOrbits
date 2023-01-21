using LinearAlgebra
using ForwardDiff
using SatelliteDynamics
using OrdinaryDiffEq

include("ks_transform.jl")
include("sim_dynamics.jl")
include("ks_dynamics.jl")
include("kgd_stm.jl")
include("orbit_utils.jl")

function cw_stm(n, T)
    return exp(T * cw_A_matrix(n))
end

function cw_A_matrix(n)
    return [0 0 0 1 0 0
        0 0 0 0 1 0
        0 0 0 0 0 1
        3n^2 0 0 0 2n 0
        0 0 0 -2n 0 0
        0 0 -n^2 0 0 0.0
    ]
end

""" `ya_stm(e, f, μ, h, T)`
Compute the Yamanaka-Ankerson state transition matrix.
See Alfriend (Spacecraft Formation Flying) pg 109

  * `e` is the orbit eccentricity
  * `f` is the true anomaly
  * `μ` is the gravitational constant GM
  * `h` is the magnitude of the angular momentum vector, h = ||r × ṙ ||

The expected state vector ordering is `[x, y, z, x', y', z']`
"""
function ya_stm(e, f0, fT, μ, h, T)
    # phi * inv(phi0) = (phi0' \ phi)'
    return ya_permute()' * (ya_phi(e, fT, μ, h, T) * ya_phi_inv(e, f0, μ, h, T)) * ya_permute()
end

""" `ya_permute`
convert from `[x, y, z; x', y', z']` to `[x, x', y, y', z, z']`
"""
function ya_permute()
    return P = [
        1 0 0 0 0 0
        0 0 0 1 0 0
        0 1 0 0 0 0
        0 0 0 0 1 0
        0 0 1 0 0 0
        0 0 0 0 0 1
    ]
end

""" `ya_phi(e, f, μ, h, T)`
Compute the Yamanaka-Ankerson phi matrix.
See Alfriend (Spacecraft Formation Flying) pg 109

  * `e` is the orbit eccentricity
  * `f` is the true anomaly
  * `μ` is the gravitational constant GM
  * `h` is the magnitude of the angular momentum vector, h = ||r × ṙ ||
  * `T` is the time over which the stm applies
"""
function ya_phi(e, f, μ, h, T)
    k = 1 + e * cos(f)
    s = k * sin(f)
    c = k * cos(f)
    # sp = cos(f) + e * cos(2 * f)
    # cp = -(sin(f) + e * sin(2 * f))
    sp = cos(f) - e + 2 * e * cos(f)^2
    cp = -sin(f) * (2 * e * cos(f) + 1)
    J = (μ^2 / h^3) * T

    P11 = s
    P12 = c
    P13 = 2 - 3 * e * s * J
    P14 = P15 = P16 = 0.0

    P21 = sp
    P22 = cp
    P23 = -3 * e * (sp * J + (s / k^2))
    P24 = P25 = P26 = 0.0

    P31 = c * (1 + (1 / k))
    P32 = -s * (1 + 1 / k)
    P33 = -3 * k^2 * J
    P34 = 1
    P35 = P36 = 0.0

    P41 = -2 * s
    P42 = e - 2 * c
    P43 = -3 * (1 - 2 * e * s * J)
    P44 = P45 = P46 = 0.0

    P51 = P52 = P53 = P54 = 0.0
    P55 = cos(f)
    P56 = sin(f)

    P61 = P62 = P63 = P64 = 0.0
    P65 = -sin(f)
    P66 = cos(f)

    P1 = [P11 P12 P13 P14 P15 P16]
    P2 = [P21 P22 P23 P24 P25 P26]
    P3 = [P31 P32 P33 P34 P35 P36]
    P4 = [P41 P42 P43 P44 P45 P46]
    P5 = [P51 P52 P53 P54 P55 P56]
    P6 = [P61 P62 P63 P64 P65 P66]

    ϕ = [P1; P2; P3; P4; P5; P6]

    return ϕ
end

""" `ya_phi_inv(e, f, μ, h, T)`
Compute the Yamanaka-Ankerson phi inverse matrix.
See Alfriend (Spacecraft Formation Flying) pg 109

  * `e` is the orbit eccentricity
  * `f` is the true anomaly
  * `μ` is the gravitational constant GM
  * `h` is the magnitude of the angular momentum vector, h = ||r × ṙ ||
"""
function ya_phi_inv0(e, f)
    k = 1 + e * cos(f)
    s = k * sin(f)
    c = k * cos(f)
    # sp = cos(f) + e * cos(2 * f)
    # cp = -(sin(f) + e * sin(2 * f))
    sp = cos(f) - e + 2 * e * cos(f)^2
    cp = -sin(f) * (2 * e * cos(f) + 1)
    η = sqrt(1 - e^2)

    P11 = -3 * s * (k + e^2) / k^2
    P12 = c - 2 * e
    P13 = 0
    P14 = -s * (k + 1) / k
    P15 = 0
    P16 = 0

    P21 = -3 * (e + c / k)
    P22 = -s
    P23 = 0
    P24 = -(c * (k + 1) / k + e)
    P25 = 0
    P26 = 0

    P31 = (3 * k - η^2)
    P32 = e * s
    P33 = 0
    P34 = k^2
    P35 = 0
    P36 = 0

    P41 = -3 * e * s * (k + 1) / k^2
    P42 = -2 + e * c
    P43 = η^2
    P44 = -e * s * (k + 1) / k
    P45 = 0
    P46 = 0

    P51 = P52 = P53 = P54 = 0
    P55 = η^2 * cos(f)
    P56 = -η^2 * sin(f)

    P61 = P62 = P63 = P64 = 0
    P65 = η^2 * sin(f)
    P66 = η^2 * cos(f)

    P1 = [P11 P12 P13 P14 P15 P16]
    P2 = [P21 P22 P23 P24 P25 P26]
    P3 = [P31 P32 P33 P34 P35 P36]
    P4 = [P41 P42 P43 P44 P45 P46]
    P5 = [P51 P52 P53 P54 P55 P56]
    P6 = [P61 P62 P63 P64 P65 P66]

    ϕ_inv = (1 / η^2) * [P1; P2; P3; P4; P5; P6]

    return ϕ_inv
end

function ks_stm(h, S)
    return exp(S * Matrix(ks_A_matrix(h)))
end

""" `transition_relative_state_via_cw(x_bar, x_tilde)`
Apply the state transition to the Cartesian relative state x_tilde.
"""
function transition_relative_state_via_cw(x_tilde_rtn, GM, a_bar, T)
    n = sqrt(GM / a_bar^3)
    return cw_stm(n, T) * x_tilde_rtn
end

""" `ya_stm(e, f, μ, h, T)`
Compute the Yamanaka-Ankerson state transition matrix.
See Alfriend (Spacecraft Formation Flying) pg 109

  * `x_tilde` cartesian relative state
  * `x_bar` cartesian ECI state
  * `fT` is the final true anomaly
  * `μ` is the gravitational constant GM
  * `h` is the magnitude of the angular momentum vector, h = ||r × ṙ ||

The expected state vector ordering is `[x, y, z, x', y', z']`
"""
function transition_relative_state_via_ya(x_tilde_rtn, x_bar, f0, fT, GM, T)
    @warn "Not properly implemented"
    a0, e0, i0, Ω0, ω0, M0 = SatelliteDynamics.sCARTtoOSC(x_bar, GM=GM)

    h0 = norm(cross(x_bar[1:3], x_bar[4:6]))

    x_tilde_T = ya_stm(e0, f0, fT, GM, h0, T) * x_tilde_rtn

    return x_tilde_T
end

""" `transition_relative_state_via_ks(x_tilde, x_bar, GM, T)`
Transition a Cartesian relative state using the KS dynamics defined at `x_bar`.
  * `x_tilde` cartesian relative state
  * `x_bar` cartesian ECI state
  * `GM` gravitational constant
  * `T` orbit time
"""
function transition_relative_state_via_ks(x_tilde, x_bar, GM, T)

    q_bar = state_cart_to_ks_newton(x_bar)
    q_hat = state_cart_to_ks_newton(x_bar .+ x_tilde; p_near=q_bar[1:4])

    h = ks_h_energy(q_bar[1:4], q_bar[5:8], GM)
    ΦT = ks_stm(h, T)
    q_bar_T = ΦT * q_bar
    q_hat_T = ΦT * q_hat

    x_bar_T = state_ks_to_cart(q_bar_T)
    x_hat_T = state_ks_to_cart(q_hat_T)

    x_tilde_T = x_hat_T - x_bar_T

    return x_tilde_T
end

function transition_relative_state_via_ode(x_tilde, x_bar, GM, T)

    function scaled_keplarian_dynamics!(ẋ, x, p, t)
        ẋ .= keplarian_dynamics(x[1:6], zeros(3), t; GM=GM) # no control inputs
    end

    x_hat = x_bar + x_tilde

    prob_hat = ODEProblem(scaled_keplarian_dynamics!, x_hat, (0.0, T))
    prob_bar = ODEProblem(scaled_keplarian_dynamics!, x_bar, (0.0, T))


    sol_hat = solve(prob_hat, Tsit5(); abstol=1e-12, reltol=1e-13)
    sol_bar = solve(prob_bar, Tsit5(); abstol=1e-12, reltol=1e-13)

    x_hat_T = sol_hat.u[end][1:6]
    x_bar_T = sol_bar.u[end][1:6]

    x_tilde_T = x_hat_T - x_bar_T

    return x_tilde_T, x_bar_T
end

function transition_relative_state_via_rk4(x_tilde, x_bar, GM, T, dt)

    function scaled_keplarian_dynamics(x_k, u_k, t)
        keplarian_dynamics(x_k, u_k, t; GM=GM)
    end

    x_hat = x_bar + x_tilde

    x_hat_T = x_hat
    x_bar_T = x_bar
    Nk = Int(floor(T / dt))

    for k = 1:Nk-1
        x_hat_T = rk4(scaled_keplarian_dynamics, x_hat_T, zeros(3), 0.0, dt)
        x_bar_T = rk4(scaled_keplarian_dynamics, x_bar_T, zeros(3), 0.0, dt)
    end

    x_tilde_T = x_hat_T - x_bar_T

    return x_tilde_T
end

function transition_relative_state_via_ode_perturbed(x_tilde, x_bar, sat_model::SatelliteModel, T)

    function perturbed_dynamics!(ẋ, x, p, t)
        ẋ .= cartesian_J2_drag_perturbed_dynamics(x, zeros(3), t, sat_model)
    end

    x_hat = x_bar + x_tilde

    prob_hat = ODEProblem(perturbed_dynamics!, x_hat, (0.0, T))
    prob_bar = ODEProblem(perturbed_dynamics!, x_bar, (0.0, T))


    sol_hat = solve(prob_hat, Tsit5(); abstol=1e-12, reltol=1e-13)
    sol_bar = solve(prob_bar, Tsit5(); abstol=1e-12, reltol=1e-13)

    x_hat_T = sol_hat.u[end]
    x_bar_T = sol_bar.u[end]

    x_tilde_T = x_hat_T - x_bar_T

    return x_tilde_T, x_bar_T
end

function transition_relative_state_via_rk4_perturbed(x_tilde, x_bar, sat_model::SatelliteModel, T, dt)

    function perturbed_dynamics(x_k, u_k, t)
        cartesian_J2_drag_perturbed_dynamics(x_k, u_k, t, sat_model)
    end

    x_hat = x_bar + x_tilde

    x_hat_T = x_hat
    x_bar_T = x_bar
    Nk = Int(floor(T / dt))

    for k = 1:Nk-1
        x_hat_T = rk4(perturbed_dynamics, x_hat_T, zeros(3), 0.0, dt)
        x_bar_T = rk4(perturbed_dynamics, x_bar_T, zeros(3), 0.0, dt)
    end

    x_tilde_T = x_hat_T - x_bar_T

    return x_tilde_T, x_bar_T
end

function ks_J2_perturbed_stm_in_ks(x_bar_ks, T, sat_model)
    # scaling
    d_scale = sat_model.distance_scale
    t_scale = sat_model.time_scale

    GM_scaled = SD.GM_EARTH / (d_scale^3 / t_scale^2)
    h = ks_h_energy(x_bar_ks[1:4], x_bar_ks[5:8], GM_scaled)

    # vector including control inputs
    Nu = 3
    z0_ks = [x_bar_ks; h; zeros(Nu)]
    Nz = length(z0_ks)
    xh0_ks = z0_ks[1:9]

    # dynamics with single vector input for ForwardDiff
    ks_full_dynamics_vec(z) = [ks_full_dynamics(sat_model, z[1:9], z[10:12]); zeros(Nu)]

    # state + stm dynamics
    function ks_stm_dynamics(xh_Φ, p, s)
        xh = xh_Φ[1:9]
        u = zeros(3)
        Φ = reshape(xh_Φ[10:end], Nz, Nz)
        z = [xh; u]
        z_dot = ks_full_dynamics_vec(z)

        ∂f∂z = ForwardDiff.jacobian(ks_full_dynamics_vec, z)

        Φ_dot = ∂f∂z * Φ

        return [z_dot[1:9]; vec(Φ_dot)]
    end

    function ks_stm_dynamics!(xh_Φ_dot, xh_Φ, p, s)
        xh_Φ_dot .= ks_stm_dynamics(xh_Φ, p, s)
    end

    xh_Φ_0 = [xh0_ks; vec(I(Nz))]

    prob = ODEProblem(ks_stm_dynamics!, xh_Φ_0, (0.0, T))
    sol = solve(prob, Tsit5(); abstol=1e-12, reltol=1e-13)

    xh_Φ = sol.u[end]

    ks_ode_stm = reshape(xh_Φ[10:end], Nz, Nz)
    xhT_ks = xh_Φ[1:9]

    return ks_ode_stm, xh0_ks, xhT_ks
end

function ks_J2_perturbed_stm(x_bar_scaled, T; with_control=false,
    sat_model=SatelliteModel( # No Drag, With J2
        1.0, # satellite_mass::Float64, -- Unused, but should be nonzero
        0.0, # drag_area::Float64, -- No Drag
        d_scale, # distance_scale::Float64, 
        t_scale, # time_scale::Float64
        1.0, # u_scale::Float64
        1.0, # rel_scale
        # GM::Float64=SatelliteDynamics.GM_EARTH,
        # J2::Float64=SatelliteDynamics.J2_EARTH,
        # R_earth::Float64=SatelliteDynamics.R_EARTH,
        # rho::Float64=1e-12, # TODO: need better drag model
        # CD::Float64=1.0,
        add_perturbations=true
    ))

    x_bar_ks = state_cart_to_ks_newton(x_bar_scaled)
    ks_ode_stm, xh0_ks, xhT_ks = ks_J2_perturbed_stm_in_ks(x_bar_ks, T, sat_model)

    if with_control
        return ks_ode_stm, xh0_ks, xhT_ks
    else
        return ks_ode_stm[1:9, 1:9], xh0_ks, xhT_ks
    end
end

function transition_relative_state_via_ks_h(q_tilde, q_bar, q_bar_T, ΦT, GM)
    q_hat = q_tilde .+ q_bar
    h_bar = ks_h_energy(q_bar[1:4], q_bar[5:8], GM)
    h_tilde = ks_h_energy(q_hat[1:4], q_hat[5:8], GM) - h_bar

    qh_tilde_T = ΦT * [q_tilde; h_tilde]
    q_tilde_T = qh_tilde_T[1:8]

    q_hat_T = q_tilde_T + q_bar_T

    x_bar_T = state_ks_to_cart(q_bar_T)
    x_hat_T = state_ks_to_cart(q_hat_T)

    x_tilde_T = x_hat_T - x_bar_T

    return x_tilde_T
end

function transition_relative_state_via_ks_h_from_cart(x_tilde, x_bar, q_bar_T, ΦT, GM)
    q_bar = state_cart_to_ks_newton(x_bar)
    q_hat = state_cart_to_ks_newton(x_bar .+ x_tilde; p_near=q_bar[1:4])
    q_tilde = q_hat .- q_bar
    return transition_relative_state_via_ks_h(q_tilde, q_bar, q_bar_T, ΦT, GM)
end

"""

"""
function transition_relative_state_via_kgd(x_tilde, x_bar, x_bar_T, T, sat_model::SatelliteModel)
    d_scale = sat_model.distance_scale
    t_scale = sat_model.time_scale
    GM = sat_model.GM
    RE = sat_model.R_earth
    J2 = sat_model.J2 # nondimensional

    # unscale states
    unscale = [d_scale * ones(3); d_scale * ones(3) / t_scale]
    x_tilde_unscaled = x_tilde .* unscale
    x_bar_unscaled = x_bar .* unscale
    x_bar_T_unscaled = x_bar_T .* unscale

    # transform states to kgd orbital elements
    α_d = kgd_cart_to_nonsingular_orbital_elements(x_bar_unscaled .+ x_tilde_unscaled, GM)
    α_c = kgd_cart_to_nonsingular_orbital_elements(x_bar_unscaled, GM)

    α_c_osc = kgd_osc_orbital_elements_from_nonsingular(α_c)

    # find kgd relative orbital elements (ROEs)
    δα = kgd_nonsingular_relative_orbital_elements(α_d, α_c)

    # get kgd STM
    Φ = kgd_nonsingular_stm(α_c_osc, T, J2, RE, GM)

    # transform kgd ROEs
    δα_T = Φ * δα

    # convert from ROEs to OEs
    α_c_T = kgd_cart_to_nonsingular_orbital_elements(x_bar_T_unscaled, GM) # this is introducing the error
    α_d_T = kgd_d_osc_from_nonsingular_relative(δα_T, α_c_T)

    # convert from OEs to cartesian
    x_hat_T_unscaled = kgd_nonsingular_orbital_elements_to_cart(α_d_T, GM)

    # relative cartesian state
    x_tilde_T_unscaled = x_hat_T_unscaled - x_bar_T_unscaled

    x_tilde_T = x_tilde_T_unscaled ./ unscale

    return x_tilde_T

end
