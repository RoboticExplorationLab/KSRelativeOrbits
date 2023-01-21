using LinearAlgebra
using ForwardDiff
using SatelliteDynamics
using OSQP
using SparseArrays
using StaticArrays

struct TrajOptQP{S}
  Pcost::SparseMatrixCSC{Float64,Int}
  qcost::Vector{Float64}
  Adyn::SparseMatrixCSC{Float64,Int}
  Acon::SparseMatrixCSC{Float64,Int}
  lb_dyn::Vector{Float64}
  ub_dyn::Vector{Float64}
  lb_con::Vector{Float64}
  ub_con::Vector{Float64}
  Nsteps::Int
  solver::S
  Xref::Vector{Vector{Float64}}
  Uref::Vector{Vector{Float64}}
  times::Vector{Float64}
  Nx::Int
  Nu::Int
  xinds::Vector{Vector{Int}}
  uinds::Vector{Vector{Int}}
  x1::Vector{Float64}
end

function TrajOptQP_OSQP(Nx::Integer, Nu::Integer, Nsteps::Integer; Nd_additional::Integer=0)
  Np = (Nsteps - 1) * (Nx + Nu)   # number of primals
  Nd_dynamics = (Nsteps - 1) * Nx # duals for dynamics
  Nd = Nd_dynamics + Nd_additional # total number of dual variables
  Pcost = spzeros(Np, Np)
  qcost = zeros(Np)
  Adyn = spzeros(Nd_dynamics, Np)
  Acon = spzeros(Nd_additional, Np)
  lb_dyn = zeros(Nd_dynamics)
  ub_dyn = zeros(Nd_dynamics)
  lb_con = zeros(Nd_additional)
  ub_con = zeros(Nd_additional)
  Xref = [zeros(Nx) for k = 1:Nsteps]
  Uref = [zeros(Nu) for k = 1:Nsteps]
  tref = zeros(Nsteps)
  xinds = [((Nu+1):(Nu+Nx)) .+ (k - 1) * (Nx + Nu) for k = 1:Nsteps-1]
  uinds = [(1:Nu) .+ (k - 1) * (Nx + Nu) for k = 1:Nsteps-1]
  x1 = zeros(Nx)
  solver = OSQP.Model()
  TrajOptQP{OSQP.Model}(Pcost, qcost, Adyn, Acon, lb_dyn, ub_dyn, lb_con, ub_con, Nsteps, solver, Xref, Uref, tref, Nx, Nu, xinds, uinds, x1)
end

""" `buildQP_cost!(qp::TrajOptQP, Q, R, Qf)`

  * `qp::TrajOptQP`: The TrajOptQP struct to store the optimization problem in
  * `Q`: An `(Nx, Nx)` state cost matrix, must be positive definite
  * `R`: An `(Nu, Nu` control cost matrix, must be positive definite
  * `Qf`: The `(Nx, Nx)` terminal cost matrix, must be positive definite
"""
function buildQP_cost!(qp::TrajOptQP, Q, R, Qf)
  Nx = qp.Nx
  Nu = qp.Nu
  Ns = qp.Nsteps

  # set up P cost matrix
  for k = 1:Ns-1
    xi = qp.xinds[k]
    ui = qp.uinds[k]
    qp.Pcost[xi, xi] .= Q
    qp.Pcost[ui, ui] .= R
  end
  xi_end = qp.xinds[end]
  qp.Pcost[xi_end, xi_end] .= Qf

end

""" `buildQP_dynamics_constraints!(qp::TrajOptQP, Adyn_list, Bdyn_list, Ddyn_list; kwargs...)`

  * `qp::TrajOptQP`: The TrajOptQP struct to store the optimization problem in
  * `Adyn_list`: An `Nsteps` long list of `(Nx, Nx)` dimension dynamics matrices (discrete time)
  * `Bdyn_list`: An `Nsteps` long list of `(Nx, Nu)` dimension control matrices (discrete time)
  * `Ddyn_list`: An `Nsteps` long list of `(Nx,)` dimension constant disturbance vectors
  * `x1`: The initial state
"""
function buildQP_dynamics!(qp::TrajOptQP, Adyn_list, Bdyn_list, Ddyn_list, x1)
  Nx = qp.Nx
  Nu = qp.Nu
  Ns = qp.Nsteps

  row_inds = [(1:Nx) .+ (k - 1) * Nx for k = 1:Ns-1]

  # first row is [B -I ...] * [u1; x2] = [-Ax1 - D1]
  ui_1 = qp.uinds[1]
  rr_1 = row_inds[1]
  xi_1 = qp.xinds[1]
  qp.Adyn[rr_1, ui_1] .= Bdyn_list[1]
  qp.Adyn[rr_1, xi_1] .= -I(Nx)

  qp.lb_dyn[rr_1] .= -Adyn_list[1] * x1 - Ddyn_list[1]
  qp.ub_dyn[rr_1] .= -Adyn_list[1] * x1 - Ddyn_list[1]

  # all other dynamics constraints are of the form 
  # [... Ak Bk -I ...] * [xk; uk; xk+1] = [-Dk]
  for k = 2:Ns-1
    rri = row_inds[k]
    ui = qp.uinds[k]
    xi = qp.xinds[k-1]
    xi_n = qp.xinds[k]

    # A constraint matrix
    qp.Adyn[rri, ui] .= Bdyn_list[k]
    qp.Adyn[rri, xi] .= Adyn_list[k]
    qp.Adyn[rri, xi_n] .= -I(Nx)

    # upper/lower bounds
    qp.lb_dyn[rri] .= -Ddyn_list[k]
    qp.ub_dyn[rri] .= -Ddyn_list[k]
  end

  qp.x1 .= x1

end


""" `buildQP!(qp::TrajOptQP, Adyn_list, Bdyn_list, Ddyn_list, x1, Q, R, Qf; kwargs...)`

  * `qp::TrajOptQP`: The TrajOptQP struct to store the optimization problem in
  * `Adyn_list`: An `Nsteps` long list of `(Nx, Nx)` dimension dynamics matrices (discrete time)
  * `Bdyn_list`: An `Nsteps` long list of `(Nx, Nu)` dimension control matrices (discrete time)
  * `Ddyn_list`: An `Nsteps` long list of `(Nx,)` dimension constant disturbance vectors
  * `x1`: The initial state
  * `Q`: An `(Nx, Nx)` state cost matrix, must be positive definite
  * `R`: An `(Nu, Nu` control cost matrix, must be positive definite
  * `Qf`: The `(Nx, Nx)` terminal cost matrix, must be positive definite
  * `kwargs`: Arguments to pass to the OSQP solver
"""
function buildQP!(qp::TrajOptQP, Adyn_list, Bdyn_list, Ddyn_list, x1, Q, R, Qf; kwargs...)
  buildQP_cost!(qp, Q, R, Qf)
  buildQP_dynamics!(qp, Adyn_list, Bdyn_list, Ddyn_list, x1)
  initialize_solver!(qp; kwargs...)
  return nothing
end

""" `initialize_solver!(qp::TrajOptQP; tol=1e-6, verbose=false)` 

Initialize the internal solver once the QP matrices are initialized in 
the `qp` problem struct.
"""
function initialize_solver!(qp::TrajOptQP; tol=1e-5, verbose=false)
  A = [qp.Adyn; qp.Acon]
  lb = [qp.lb_dyn; qp.lb_con]
  ub = [qp.ub_dyn; qp.ub_con]
  OSQP.setup!(qp.solver, P=qp.Pcost, q=qp.qcost, A=A, l=lb, u=ub,
    verbose=verbose, eps_rel=tol, eps_abs=tol, eps_prim_inf=1e-5, eps_dual_inf=1e-5, polish=1, max_iter=10000)
end

""" `solve_trajopt(qp::TrajOptQP)`
  
  * `qp::TrajOptQP`: The problem definition struct
  
Returns:

  * `xtraj`: A `Nsteps` long vector of `(Nx,)` state vectors
  * `utraj`: A `Nsteps` long vector of `(Nu,)` control vectors
"""
function solve_trajopt(qp::TrajOptQP)
  results = OSQP.solve!(qp.solver)

  if results.info.status == :Solved
    @info "Solve completed successfully, $(results.info)"
  else
    @warn "Solve Failed with status $(results.info.status), $(results.info)"
  end

  utraj = [results.x[qp.uinds[k]] for k = 1:qp.Nsteps-1]
  xtraj = [qp.x1, [results.x[qp.xinds[k]] for k = 1:qp.Nsteps-1]...]

  return xtraj, utraj
end

function setup_u_bound_constraints!(qp::TrajOptQP, umin, umax)
  Nk = qp.Nsteps
  Nu = qp.Nu
  # Constraints on u
  for k = 1:Nk-1
    ui = qp.uinds[k]
    for l = 1:Nu
      rr = (k - 1) * Nu + l
      qp.Acon[rr, ui[l]] = 1.0
      qp.lb_con[rr] = umin
      qp.ub_con[rr] = umax
    end
  end
end