# Implementation of the nonsingular state transition matrix from 
# New State Transition Matrices for Spacecraft Relative Motion in Perturbed Orbits
# by
# Koenig, Guffanti, and D'Amico (KGD)

include("orbit_utils.jl")

""" `kgd_nonsingular_relative_orbital_elements(α_d, α_c)`
Compute the relative orbital elements from the nonsingular orbital elements
  * `α_d` - deputy nonsingular orbital elements
  * `α_c` - chief nonsingular orbital elements

nonsingular orbital elements: `a, l, ex, ey, ix, iy`
"""
function kgd_nonsingular_relative_orbital_elements(α_d, α_c)
    δα = α_d .- α_c
    δα[1] = δα[1] / α_c[1] # scale 

    return δα
end

function kgd_d_osc_from_nonsingular_relative(δα, α_c)
    # unscale a
    δα_unscaled = δα .* [α_c[1], 1, 1, 1, 1, 1]

    α_d = δα_unscaled .+ α_c

    return α_d
end


""" `kgd_nonsingular_orbital_elements_from_osc(α_osc)`
"""
function kgd_nonsingular_orbital_elements_from_osc(α_osc)
    replace!(α_osc, -0.0 => 0.0)
    a, e, i, Ω, ω, M = α_osc

    l = M + ω + Ω
    l = rem2pi(l, RoundDown)
    ex = e * cos(ω + Ω)
    ey = e * sin(ω + Ω)
    ix = tan(i / 2) * cos(Ω)
    iy = tan(i / 2) * sin(Ω)

    α_ns = [a, l, ex, ey, ix, iy]
    replace!(α_ns, -0.0 => 0.0)
    return α_ns
end

""" `kgd_osc_orbital_elements_from_nonsingular(  `
"""
function kgd_osc_orbital_elements_from_nonsingular(α_ns)
    replace!(α_ns, -0.0 => 0.0)
    a, l, ex, ey, ix, iy = α_ns

    e = sqrt(ex^2 + ey^2)
    i = atan(2 * sqrt(ix^2 + iy^2), 1 - ix^2 - iy^2)
    Ω = atan(iy, ix) #acos(ix / tan(i / 2))
    ω = atan((ey * ix - ex * iy), (ex * ix + ey * iy)) #acos(ex / e) - Ω
    M = l - Ω - ω

    α_osc = [a, e, i, Ω, ω, M]
    replace!(α_osc, -0.0 => 0.0)
    return α_osc
end

function kgd_cart_to_nonsingular_orbital_elements(x_cart, GM)
    replace!(x_cart, -0.0 => 0.0)
    α_osc = state_CART_to_OSC(x_cart, GM)
    return kgd_nonsingular_orbital_elements_from_osc(α_osc)
end

function kgd_nonsingular_orbital_elements_to_cart(α_ns, GM)
    α_osc = kgd_osc_orbital_elements_from_nonsingular(α_ns)

    replace!(α_osc, -0.0 => 0.0)
    x_cart = SD.sOSCtoCART(α_osc, use_degrees=false, GM=GM)

    replace!(x_cart, -0.0 => 0.0)
    return x_cart
end


function kgd_nonsingular_stm(α_chief_osc, τ, J2, Re, GM)
    replace!(α_chief_osc, -0.0 => 0.0)
    a, e, i, Ω_i, ω_i, M_i = α_chief_osc

    # equation 9
    n = sqrt(GM) / a^(3 / 2)

    # Equations 14 - 16
    η = sqrt(1 - e^2)
    κ = (0.75) * J2 * (Re^2) * sqrt(GM) / ((a^(7 / 2)) * η^4)
    E = 1 + η
    F = 4 + 3η
    G = 1 / (η^2)
    P = 3 * (cos(i)^2) - 1
    Q = 5 * (cos(i)^2 - 1)
    R = cos(i)
    S = sin(2 * i)
    T = sin(i)^2
    U = sin(i)
    V = tan(i / 2)
    W = cos(i / 2)^2

    # Appendix A1
    ω_dot = κ * Q
    Ω_dot = -2 * κ * R
    Ω_f = Ω_i + Ω_dot * τ
    ω_f = ω_i + ω_dot * τ

    # really e_^star
    exi = e * cos(ω_i + Ω_i)
    eyi = e * sin(ω_i + Ω_i)
    exf = e * cos(ω_f + Ω_f)
    eyf = e * sin(ω_f + Ω_f)

    # really i_^star
    ixi = tan(i / 2) * cos(Ω_i)
    iyi = tan(i / 2) * sin(Ω_i)
    ixf = tan(i / 2) * cos(Ω_f)
    iyf = tan(i / 2) * sin(Ω_f)

    # Appendix A4

    # row 1
    P11 = 1
    P12 = P13 = P14 = P15 = P16 = 0.0

    P1 = [P11 P12 P13 P14 P15 P16]

    # row 2
    P21 = -((3 * n / 2) + (7 * κ / 2) * (η * P + Q - 2 * R)) * τ
    P22 = 1.0
    P23 = κ * exi * G * (3 * η * P + 4 * Q - 8 * R) * τ
    P24 = κ * eyi * G * (3 * η * P + 4 * Q - 8 * R) * τ
    P25 = 2 * κ * W * (-(3 * η + 5) * S + 2 * U) * cos(Ω_i) * τ
    P26 = 2 * κ * W * (-(3 * η + 5) * S + 2 * U) * sin(Ω_i) * τ
    P2 = [P21 P22 P23 P24 P25 P26]

    # row 3
    P31 = (7 / 2) * κ * eyf * (Q - 2 * R) * τ
    P32 = 0
    P33 = cos((ω_dot + Ω_dot) * τ) - 4 * κ * eyf * exi * G * (Q - 2 * R) * τ
    P34 = -sin((ω_dot + Ω_dot) * τ) - 4 * κ * eyf * eyi * G * (Q - 2 * R) * τ
    P35 = -2 * κ * eyf * W * (-5 * S + 2 * U) * cos(Ω_i) * τ
    P36 = -2 * κ * eyf * W * (-5 * S + 2 * U) * sin(Ω_i) * τ
    P3 = [P31 P32 P33 P34 P35 P36]

    # row 4
    P41 = -(7 / 2) * κ * exf * (Q - 2 * R) * τ
    P42 = 0
    P43 = sin((ω_dot + Ω_dot) * τ) + 4 * κ * exf * exi * G * (Q - 2 * R) * τ
    P44 = cos((ω_dot + Ω_dot) * τ) + 4 * κ * exf * eyi * G * (Q - 2 * R) * τ
    P45 = 2 * κ * exf * W * (-5 * S + 2 * U) * cos(Ω_i) * τ
    P46 = 2 * κ * exf * W * (-5 * S + 2 * U) * sin(Ω_i) * τ
    P4 = [P41 P42 P43 P44 P45 P46]

    # row 5
    P51 = -7 * κ * iyf * R * τ
    P52 = 0
    P53 = 8 * κ * exi * iyf * G * R * τ
    P54 = 8 * κ * eyi * iyf * G * R * τ
    P55 = cos(Ω_dot * τ) - 4 * κ * iyf * U * W * cos(Ω_i) * τ
    P56 = -sin(Ω_dot * τ) - 4 * κ * iyf * U * W * sin(Ω_i) * τ
    P5 = [P51 P52 P53 P54 P55 P56]

    # row 6
    P61 = 7 * κ * ixf * R * τ
    P62 = 0
    P63 = -8 * κ * exi * ixf * G * R * τ
    P64 = -8 * κ * eyi * ixf * G * R * τ
    P65 = sin(Ω_dot * τ) + 4 * κ * ixf * U * W * cos(Ω_i) * τ
    P66 = cos(Ω_dot * τ) + 4 * κ * ixf * U * W * sin(Ω_i) * τ
    P6 = [P61 P62 P63 P64 P65 P66]

    Φ = [P1; P2; P3; P4; P5; P6]

    return Φ
end
