# long_pendulum_chain.jl
# Simulate a long chain of coupled pendula and animate.
#
# Model:
#   I * θ̈_n = - m g L sin(θ_n) + K (θ_{n+1} - 2 θ_n + θ_{n-1})
# Convert to first-order system for use with DifferentialEquations.jl
#
# Requires: DifferentialEquations, Plots

using DifferentialEquations
using Plots

# Uncomment and run once if you need to install
# import Pkg; Pkg.add(["DifferentialEquations","Plots"])

gr()  # GR backend

# -------------------------
# Physical / numerical params
# -------------------------
const N = 120               # number of pendula (try 120; increase/decrease for performance)
const L = 1.0               # length (m)
const m = 1.0               # mass (kg)
const I = m * L^2           # moment of inertia for point mass at distance L
const g = 9.81              # gravity m/s^2
const K = 50.0              # torsional coupling stiffness (N*m/rad)
const spacing = 1.0         # pure label: physical spacing a can be used to get continuum limits
const tmax = 12.0           # total simulation time in seconds

# Derived
const ω0sq = g / L
const coupling_coef = K / I   # multiplies discrete Laplacian

# -------------------------
# Initial condition: localized kick
# -------------------------
# Example: give a gaussian-shaped initial angle profile
using Random
rng = MersenneTwister(1234)

x = range(0, stop=1, length=N)
x0 = 0.25
sigma = 0.06
θ0 = 0.5 .* exp.(.-(x .- x0).^2 ./ (2*sigma^2))      # initial angular displacement (radians)
# Option: add a little noise if desired
θ0 .+= 0.0 .* (randn(rng, N) .* 1e-3)

ω0 = zeros(N)   # initial angular velocities

# state vector u = [θ1...θN, ω1...ωN]
u0 = vcat(θ0, ω0)

# -------------------------
# Right-hand side
# -------------------------
function chain_dynamics!(du, u, p, t)
    θ = @view u[1:N]
    ω = @view u[N+1:2N]
    dθ = @view du[1:N]
    dω = @view du[N+1:2N]

    # θ_dot = ω
    @inbounds for i in 1:N
        dθ[i] = ω[i]
    end

    # Compute discrete Laplacian with fixed-end boundary conditions:
    # θ_0 = 0 and θ_{N+1} = 0 (you can pick free-end by copying neighbor values)
    @inbounds for i in 1:N
        left  = (i == 1) ? 0.0 : θ[i-1]
        right = (i == N) ? 0.0 : θ[i+1]
        lap = right - 2.0*θ[i] + left
        dω[i] = - ω0sq * sin(θ[i]) + coupling_coef * lap
    end
end

# -------------------------
# Solve
# -------------------------
tspan = (0.0, tmax)
prob = ODEProblem(chain_dynamics!, u0, tspan)

# For oscillatory mechanical systems, Tsit5 or Vern9 works fine.
# Increase abstol/reltol for accuracy if desired.
sol = solve(prob, Tsit5(); abstol=1e-8, reltol=1e-8)

# -------------------------
# Animation: show angles as vertical displacements for clarity
# -------------------------
# We'll plot a 2D "string" where horizontal axis is pendulum index and vertical axis is angle or y-pos.
# To make it visually like pendula, compute bob positions (x_i, y_i) with small amplitude mapping.

# Mapping params for visuals
x_positions = collect(1:N) .- (N+1)/2        # center horizontally
scale_angle_to_y = 0.9 * (1.0 / maximum(abs.(θ0) .+ 1e-6))  # scale so initial shape is visible
# better mapping: convert angle θ -> bob x = pivot_x + L sin θ, y = -L cos θ
function bob_coords(θ_vector)
    xs = similar(θ_vector)
    ys = similar(θ_vector)
    for i in 1:length(θ_vector)
        xs[i] = x_positions[i] + L * sin(θ_vector[i]) * 0.6   # horizontal wiggle scaled
        ys[i] = -L * cos(θ_vector[i])                         # vertical
    end
    return xs, ys
end

# Build animation frames
anim = @animate for (ti, ui) in zip(sol.t, sol.u)
    θ = ui[1:N]
    xs, ys = bob_coords(θ)

    plt = scatter(xs, ys; markersize=4, ylim=(-L-0.2, 0.2), xlim=(minimum(x_positions)-1, maximum(x_positions)+1),
        markerstrokecolor = :black, legend=false, aspect_ratio = 0.35,
        title = @sprintf("Coupled pendula chain, t = %.2f s, N = %d", ti, N),
        xlabel = "index (spatial)", ylabel = "vertical position (m)")

    # draw strings from pivot line y=0 to bob
    for i in 1:N
        plot!([x_positions[i], xs[i]], [0.0, ys[i]]; linewidth=0.8)
    end

    # draw bobs on top
    scatter!(xs, ys; markersize=8, markercolor = :blue)

    plt
end fps=20

# Save to gif
outfile = "pendula_chain_N$(N).gif"
gif(anim, outfile, fps = 20)
println("Saved animation to: ", outfile)
