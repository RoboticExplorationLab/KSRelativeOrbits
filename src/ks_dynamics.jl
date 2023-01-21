using SatelliteDynamics
using LinearAlgebra

using ForwardDiff

include("ks_transform.jl")
include("sim_dynamics.jl")

function ks_full_dynamics(model::SatelliteModel, p_state, u_scaled)
    p = p_state[1:4]
    p_prime = p_state[5:8]
    h = p_state[9]

    # unperturbed KS dynamics
    p_pprime = ((-h / 2) * p)
    h_prime = 0.0

    # control inputs are cartesian accelerations
    u = u_scaled / model.u_scale
    a = u

    # perturbed KS dynamics
    if model.add_perturbations
        # compute cartesian state
        x_scaled = state_ks_to_inertial(p_state[1:8])

        # unscale cartesian state
        x_state = [x_scaled[1:3] * model.distance_scale
            x_scaled[4:6] * (model.distance_scale / model.time_scale)]

        # compute perturbation accelerations
        a_j2 = inertial_J2_perturbation(model, x_state)
        a_drag = inertial_drag_perturbation(model, x_state)

        a_unscaled = a_j2 + a_drag

        # scale acceleration
        a = a + a_unscaled * (model.distance_scale / (model.time_scale^2))
    end

    L = ks_L(p)
    p_pprime += ((p'p) / 2) * L' * [a; 0]
    h_prime += -2 * p_prime' * L' * [a; 0]
    return [p_prime; p_pprime; h_prime]
end
