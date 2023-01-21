using SatelliteDynamics
using LinearAlgebra

""" inertial_accel_perturbations(x, u, model, epc)
    Compute the inertial frame perturbation (non keplarian) accelerations at state `x = [position, velocity]`,
    `u = [u_drag]`, and epoch `epc`
    for the satellite defined by `model`.

    Based off: https://github.com/RoboticExplorationLab/GravNav/blob/main/jacobian_demo/scaled_dynamics.jl
"""
# function inertial_accel_perturbations(x::Array{<:Real}, u::Array{<:Real}, model::KS_Satellite, epc::Epoch)

#     r = x[1:3]
#     v = x[4:6]

#     u_drag = u[1] # between [0,1] where 0 corresponds to minimum drag and 1 corresponds to maximum drag

#     # Compute ECI to ECEF Transformation -> IAU2010 Theory
#     PN = bias_precession_nutation(epc)
#     E  = earth_rotation(epc)
#     W  = polar_motion(epc)
#     R  = W * E * PN

#     # Compute sun and moon position
#     r_sun  = sun_position(epc)
#     r_moon = moon_position(epc)

#     # Compute acceleration (eltype(x) makes this forward diff friendly)
#     a = zeros(eltype(x), 3)

#     # spherical harmonic gravity
#     a += accel_gravity(x, R, model.n_gravity, model.m_gravity)

#     # atmospheric drag
#     ρ = density_harris_priester(epc, r)
#     controlled_drag_area = (model.min_drag_area + (u * (model.max_drag_area - model.min_drag_area)))
#     a += accel_drag(
#         [r;v], 
#         ρ, 
#         model.satellite_mass, 
#         controlled_drag_area, 
#         model.drag_coefficient,
#         Array{Real, 2}(PN))

#     # SRP
#     nu = eclipse_conical(x, r_sun)
#     a += (nu * accel_srp(
#         x, 
#         r_sun, 
#         model.satellite_mass, 
#         model.srp_area,
#         model.srp_coefficient))

#     # third body sun
#     a += accel_thirdbody_sun(x, r_sun)

#     # third body moon
#     a += accel_thirdbody_moon(x, r_moon)

#     return a
# end

struct SatelliteModel
    # fixed parameters
    satellite_mass::Float64
    drag_area::Float64
    distance_scale::Float64
    time_scale::Float64
    u_scale::Float64
    rel_scale::Float64
    # constants
    GM::Float64
    J2::Float64
    R_earth::Float64
    rho::Float64
    CD::Float64
    add_perturbations::Bool
end

# default constructor to set fixed parameters
function SatelliteModel(
    satellite_mass::Float64,
    drag_area::Float64,
    distance_scale::Float64,
    time_scale::Float64,
    u_scale::Float64,
    rel_scale::Float64;
    GM::Float64=SatelliteDynamics.GM_EARTH,
    J2::Float64=SatelliteDynamics.J2_EARTH,
    R_earth::Float64=SatelliteDynamics.R_EARTH,
    rho::Float64=1e-12, # TODO: need better drag model
    CD::Float64=1.0,
    add_perturbations=true
)
    return SatelliteModel(
        satellite_mass,
        drag_area,
        distance_scale,
        time_scale,
        u_scale,
        rel_scale,
        GM,
        J2,
        R_earth,
        rho,
        CD,
        add_perturbations,
    )
end

function inertial_J2_perturbation(model::SatelliteModel, x)
    r = x[1:3]
    r_mag = sqrt(r' * r)

    J2 = model.J2
    μ = model.GM
    R_earth = model.R_earth

    a_J2 = -(3.0 / 2.0) * J2 * (μ / r_mag^2) * (R_earth / r_mag)^2 * [
               (1 - 5 * (r[3] / r_mag)^2) * r[1] / r_mag,
               (1 - 5 * (r[3] / r_mag)^2) * r[2] / r_mag,
               (3 - 5 * (r[3] / r_mag)^2) * r[3] / r_mag]

    return a_J2
end

function inertial_drag_perturbation(model::SatelliteModel, x)

    v = x[4:6]
    v_mag = sqrt(v' * v)

    ρ = model.rho
    C_D = model.CD
    mass = model.satellite_mass
    A = model.drag_area

    F_drag = -0.5 * ρ * C_D * A * v_mag .* v
    a_drag = (1.0 / mass) * F_drag

    return a_drag
end

function keplarian_dynamics(x_k, u_k, t; GM=SD.GM_EARTH)
    r = x_k[1:3]
    rmag = norm(r)
    v = x_k[4:6]

    rdot = v
    vdot = u_k - ((GM_scaled / rmag^3) * r)

    return [rdot; vdot]
end

function cartesian_J2_drag_perturbed_dynamics(x, u_scaled, t, model)
    r = x[1:3]
    rmag = norm(r)
    v = x[4:6]

    u_k = u_scaled / model.u_scale

    rdot = v
    vdot = u_k - ((GM_scaled / rmag^3) * r)


    # perturbed dynamics
    if model.add_perturbations
        # unscale cartesian state
        x_state = [r * model.distance_scale
            v * (model.distance_scale / model.time_scale)]

        # compute perturbation accelerations
        a_j2 = inertial_J2_perturbation(model, x_state)
        a_drag = inertial_drag_perturbation(model, x_state)

        a_unscaled = a_j2 + a_drag

        # scale acceleration
        vdot .+= a_unscaled * (model.distance_scale / (model.time_scale^2))
    end

    return [rdot; vdot]
end


function rk4(dynamics, x_k, u_k, t, dt)

    k1 = dynamics(x_k, u_k, t)
    k2 = dynamics(x_k .+ 0.5 .* dt .* k1, u_k, t + 0.5 * dt)
    k3 = dynamics(x_k .+ 0.5 .* dt .* k2, u_k, t + 0.5 * dt)
    k4 = dynamics(x_k .+ dt .* k3, u_k, t + dt)

    x_kp1 = x_k .+ (1.0 / 6.0) * dt .* (k1 .+ (2 .* k2) .+ (2 .* k3) .+ k4)

    return x_kp1
end