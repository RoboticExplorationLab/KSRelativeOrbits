using LinearAlgebra
using ForwardDiff
using SatelliteDynamics
using OSQP
using SparseArrays
using StaticArrays

include("trajopt_qp.jl")
include("ks_transform.jl")
include("state_transition_matrices.jl")
include("sim_dynamics.jl")

""" `discrete_relative_dynamics(q_hat_1, q_bar_1, ds)`
"""
function discrete_relative_dynamics_perturbed(qh_bar_0, Nk, ds, model::SatelliteModel)
    Nx = length(qh_bar_0)
    Nu = 3

    A_bar_list = [spzeros(Nx, Nx) for k = 1:Nk-1]
    B_bar_list = [spzeros(Nx, Nu) for k = 1:Nk-1]
    qh_bar_list = [zeros(Nx) for k = 1:Nk]
    qh_bar_list[1] .= qh_bar_0

    for k = 1:Nk-1

        q_bar_k = qh_bar_list[k][1:8]

        A_bar, B_bar, _, qh_bar_kp1 = discretize_ks_dynamics_perturbed(q_bar_k, ds, model)

        # get ECI to RTN transformation
        x_bar_k = state_ks_to_inertial(q_bar_k)
        R = SatelliteDynamics.rRTNtoECI(x_bar_k)

        A_bar_list[k] .= A_bar
        B_bar_list[k] .= B_bar * R
        qh_bar_list[k+1] .= qh_bar_kp1

    end

    return A_bar_list, B_bar_list, qh_bar_list
end

""" `discretize_ks_dynamics(q, ds, GM)`
"""
function discretize_ks_dynamics_perturbed(q, ds, model::SatelliteModel)

    ks_rk_stm, xh0_ks, xhT_ks = ks_J2_perturbed_stm_in_ks(q, ds, sat_model)

    Ad = ks_rk_stm[1:9, 1:9]
    Bd = ks_rk_stm[1:9, 10:12]

    return Ad, Bd, xh0_ks, xhT_ks
end

function rollout_discrete_trajectory_perturbed(q, ds, GM, Nk)
    Nx = length(q)
    qtraj = [zeros(Nx) for k = 1:Nk]
    qtraj[1] .= q

    for k = 1:Nk-1
        Ad, _ = discretize_ks_dynamics(qtraj[k], ds, GM)
        qtraj[k+1] = Ad * qtraj[k]
    end

    return qtraj
end

function solve_relative_maneuver_perturbed(q, q_bar, ds, GM, Nk, Q, R, Qf, sat_model; num_constraints=0, constraints_setup_fn=nothing)
    Nx = length(q) + 1 # include h_tilde in solve
    Nu = 3
    qp = TrajOptQP_OSQP(Nx, Nu, Nk; Nd_additional=num_constraints)

    if num_constraints > 0 && !isnothing(constraints_setup_fn)
        constraints_setup_fn(qp, sat_model)
    end

    h_bar = ks_h_energy(q_bar[1:4], q_bar[5:8], GM)
    h = ks_h_energy(q[1:4], q[5:8], GM)

    qh_0 = [q; h]
    qh_bar_0 = [q_bar; h_bar]

    qh_tilde_0 = qh_0 .- qh_bar_0

    A_bar_list, B_bar_list, qh_bar_traj = discrete_relative_dynamics_perturbed(qh_bar_0, Nk, ds, sat_model)
    D_tilde_list = [zeros(Nx) for k = 1:Nk-1]

    # scale everything going into the QP
    @show qh_tilde_0
    D_tilde_list ./= sat_model.rel_scale
    qh_tilde_0 ./= sat_model.rel_scale
    @show qh_tilde_0

    buildQP!(qp, A_bar_list, B_bar_list, D_tilde_list, qh_tilde_0, Q, R, Qf)

    qh_sol_traj, u_sol_traj = solve_trajopt(qp)

    @show qh_sol_traj[1]
    qh_sol_traj .*= sat_model.rel_scale
    u_sol_traj .*= sat_model.rel_scale
    @show qh_sol_traj[1]

    q_hat_traj = [qh_sol_traj[k][1:8] .+ qh_bar_traj[k][1:8] for k = 1:Nk]

    tvec = rollout_timesteps_perturbed(q_hat_traj, ds)
    return qh_sol_traj, u_sol_traj, qh_bar_traj, tvec
end

function rollout_timesteps_perturbed(qtraj, ds)
    Nk = length(qtraj)
    tvec = zeros(Nk)
    for k = 1:Nk-1
        r = norm(qtraj[k][1:4])^2
        tvec[k+1] = tvec[k] + r * ds
    end
    return tvec
end