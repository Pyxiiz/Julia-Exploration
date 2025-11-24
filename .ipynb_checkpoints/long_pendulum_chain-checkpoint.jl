<<<<<<< HEAD:infinite_pendulum.jl
using LinearAlgebra
=======
# long_pendulum_chain.jl
# Simulate a long chain of coupled pendula and animate.
#
# Model:
#   I * θ̈_n = - m g L sin(θ_n) + K (θ_{n+1} - 2 θ_n + θ_{n-1})
# Convert to first-order system for use with DifferentialEquations.jl
#
# Requires: DifferentialEquations, Plots
using Pkg
Pkg.add("DifferentialEquations")
Pkg.add("Plots")


using DifferentialEquations
>>>>>>> d3a161496e11eedcbc64a4fc9a9518821a3b177c:long_pendulum_chain.jl
using Plots
#To make a pendulum, use make_pendulum(n= number of links, fixed_pivot = true/false).
#Recomend lower n and true pivot
#False pivot allows for the pivot to move around. 
#After rendering, will save a gif file which you can either save or display in the notebook. 
# -------------------------
# User adjustable options
# -------------------------
const FPS = 60
const DT = 0.01
const T_FINAL = 30.0

# pivot behavior (true means fixed at origin)
const fixed_pivot = true

# pivot motion if not fixed
pivot_motion(t) = (0.5 * sin(0.6 * t), 0.0)

# number of links
const n = 3

# lengths and masses
L = ones(n) .* 1.0
m = ones(n) .* 1.0
g = 9.81

# -------------------------
# Forward kinematics
# -------------------------
function forward_kinematics(theta, L, px, py)
    N = length(theta)
    xs = zeros(N+1)
    ys = zeros(N+1)
    xs[1] = px
    ys[1] = py

    for k in 1:N
        xs[k+1] = xs[k] + L[k] * sin(theta[k])
        ys[k+1] = ys[k] - L[k] * cos(theta[k])
    end

    return xs, ys
end

# -------------------------
# Mass matrix
# -------------------------
function mass_matrix(theta, L, m)
    N = length(theta)
    M = zeros(N, N)

    for k in 1:N
        for i in 1:k
            for j in 1:k
                M[i,j] += m[k] * (
                    L[i]*cos(theta[i]) * L[j]*cos(theta[j]) +
                    L[i]*sin(theta[i]) * L[j]*sin(theta[j])
                )
            end
        end
    end

    return M
end

# -------------------------
# Gravity vector
# -------------------------
function gravity_vector(theta, L, m, g)
    N = length(theta)
    G = zeros(N)

    downstream = zeros(N)
    downstream[N] = m[N]
    for i in N-1:-1:1
        downstream[i] = downstream[i+1] + m[i]
    end

    for i in 1:N
        G[i] = L[i] * sin(theta[i]) * g * downstream[i]
    end

    return G
end

# -------------------------
# Coriolis vector using numeric Christoffel symbols
# -------------------------
function coriolis_vector(theta, thetadot, L, m)
    N = length(theta)
    eps = 1e-7

    M0 = mass_matrix(theta, L, m)
    dM = Array{Float64,3}(undef, N, N, N)

    for k in 1:N
        thp = copy(theta)
        thm = copy(theta)
        thp[k] += eps
        thm[k] -= eps
        Mp = mass_matrix(thp, L, m)
        Mm = mass_matrix(thm, L, m)
        dM[:,:,k] = (Mp .- Mm) ./ (2eps)
    end

    C = zeros(N)

    for i in 1:N
        for j in 1:N
            for k in 1:N
                Γ = 0.5 * (dM[i,j,k] + dM[i,k,j] - dM[j,k,i])
                C[i] += Γ * thetadot[j] * thetadot[k]
            end
        end
    end

    return C
end

# -------------------------
# Full dynamics
# -------------------------
function acceleration(theta, thetadot, L, m, g)
    M = mass_matrix(theta, L, m)
    G = gravity_vector(theta, L, m, g)
    C = coriolis_vector(theta, thetadot, L, m)
    return M \ -(G .+ C)
end

# -------------------------
# RK4
# -------------------------
function rk4_step(theta, thetadot, dt, L, m, g)
    f = (th, thd) -> acceleration(th, thd, L, m, g)

    k1_thd = thetadot
    k1_thdd = f(theta, thetadot)

    k2_thd = thetadot .+ 0.5*dt .* k1_thdd
    k2_thdd = f(theta .+ 0.5*dt .* k1_thd, k2_thd)

    k3_thd = thetadot .+ 0.5*dt .* k2_thdd
    k3_thdd = f(theta .+ 0.5*dt .* k2_thd, k3_thd)

    k4_thd = thetadot .+ dt .* k3_thdd
    k4_thdd = f(theta .+ dt .* k3_thd, k4_thd)

    theta_new = theta .+ dt .* (k1_thd .+ 2k2_thd .+ 2k3_thd .+ k4_thd) ./ 6
    thetadot_new = thetadot .+ dt .* (k1_thdd .+ 2k2_thdd .+ 2k3_thdd .+ k4_thdd) ./ 6

    return theta_new, thetadot_new
end

# -------------------------
# Simulation and rendering
# -------------------------
function make_pendulum(n; T=T_FINAL, dt=DT, fps=FPS, fixed_pivot=true)

    # local filename helper to avoid global name conflicts
    filename(n) = "pendulum_$(n).gif"

    L = fill(1.0, n)
    m = fill(1.0, n)

    theta = zeros(n)
    for i in 1:n
        theta[i] = pi/2 + 0.2*(-1)^i*(i/n)
    end

    thetadot = zeros(n)

    anim = @animate for step in 1:Int(round(T/dt))
        t = step * dt
        theta, thetadot = rk4_step(theta, thetadot, dt, L, m, g)

        px, py = fixed_pivot ? (0.0, 0.0) : pivot_motion(t)

        xs, ys = forward_kinematics(theta, L, px, py)

        plt = plot(xs, ys,
            seriestype = :path,
            linewidth = 3,
            xlim = (-sum(L)-0.5, sum(L)+0.5),
            ylim = (-sum(L)-0.5, sum(L)+0.5),
            aspect_ratio = 1,
            legend = false,
            title = "n link pendulum (n = $n)"
        )

        scatter!(plt, xs, ys, markersize = 4)
        scatter!(plt, [px], [py], marker = :star5, markersize = 6)
        plt
    end

    gif(anim, filename(n), fps = fps)
    println("Saved GIF as: $(filename(n))")
end

