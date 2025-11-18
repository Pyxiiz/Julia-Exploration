# newtons_cradle.jl
# Newton's cradle simulation using pendulums + instantaneous elastic collisions.
# Requires: DifferentialEquations, Plots

using DifferentialEquations
using Plots

# Uncomment and run once to install required packages:
# import Pkg; Pkg.add(["DifferentialEquations","Plots"])

# Use GR backend for Plots (default usually GR)
gr()

# Parameters
const N = 5                     # number of bobs
const L = 1.0                   # pendulum length (m)
const m = 1.0                   # mass of each bob (kg), equal masses
const r = 0.06                  # radius of each bob (m)
const g = 9.81                  # gravity (m/s^2)
const spacing = 2r + 0.01       # horizontal spacing between pivots (small gap)
const pivot_xs = collect(0:spacing:(N-1)*spacing) .- ((N-1)*spacing)/2
# initial condition: give the leftmost bob an initial angle
thetas0 = zeros(N)
thetas0[1] = 0.6                # initial angle for bob 1 in radians (small-ish)
omegas0 = zeros(N)

# State vector u = [theta_1, ..., theta_N, omega_1, ..., omega_N]
u0 = vcat(thetas0, omegas0)

# ODE: for each pendulum (uncoupled except collisions)
function cradle_dynamics!(du, u, p, t)
    θ = @view u[1:N]
    ω = @view u[N+1:2N]
    dθ = @view du[1:N]
    dω = @view du[N+1:2N]

    # θ_dot = ω
    @inbounds for i in 1:N
        dθ[i] = ω[i]
    end

    # ω_dot = - (g/L) * sin(θ)    (simple pendulum, no damping)
    @inbounds for i in 1:N
        dω[i] = - (g/L) * sin(θ[i])
    end
end

# Helper: compute cartesian coordinates of bob i (x,y) given θ_i
function bob_pos(θi, pivotx)
    x = pivotx + L * sin(θi)
    y = -L * cos(θi)
    return x, y
end

# Helper: horizontal linear velocity vx of bob i from angular velocity ω_i
# vx ≈ L * ω * cos(θ) (exact for derivative of x = L sin θ)
function horiz_vx(ωi, θi)
    return L * ωi * cos(θi)
end

# Callback root functions: for each adjacent pair i and i+1, we detect when distance <= 2r
g_funs = []
for i in 1:(N-1)
    function root_pair(u, t, integrator)  # capture i
        θ = @view u[1:N]
        xi, _ = bob_pos(θ[i], pivot_xs[i])
        xj, _ = bob_pos(θ[i+1], pivot_xs[i+1])
        return (xi - xj) + 2r  # detect when (xj - xi) - 2r == 0 => xi - xj + 2r == 0
        # Note: roots are triggered when this crosses zero. We'll use condition sign to ensure approach.
    end
    push!(g_funs, root_pair)
end

# Affect: called when a root triggers. We'll swap horizontal velocities (equal-mass elastic collision).
function affect_pair!(integrator, i)
    u = integrator.u
    θ = @view u[1:N]
    ω = @view u[N+1:2N]

    # compute horizontal velocities
    vx_i = horiz_vx(ω[i], θ[i])
    vx_j = horiz_vx(ω[i+1], θ[i+1])

    # Only handle collision if approaching: relative velocity along x negative (i moving right relative to j)
    # Here we check approach: (vx_i - vx_j) > 0 means i is moving right faster than j (closing)
    # But sign depends on geometry; we take absolute check to avoid double triggers
    if (vx_i - vx_j) > 0
        # For equal masses, elastic collision swaps the velocities along the collision axis.
        new_vx_i = vx_j
        new_vx_j = vx_i

        # Convert back to angular velocities: ω = vx / (L * cos(θ))
        # Guard against cos(θ) ≈ 0 (we assume small angles; clamp if necessary)
        c_i = cos(θ[i]); c_j = cos(θ[i+1])
        if abs(c_i) < 1e-6 || abs(c_j) < 1e-6
            return  # skip if degenerate; hope this does not happen for typical initial angles
        end

        ω[i]   = new_vx_i / (L * c_i)
        ω[i+1] = new_vx_j / (L * c_j)

        # write back
        for k in 1:N
            integrator.u[N+k] = ω[k]
        end
    end
    return
end

# Build callbacks for each adjacent pair
callbacks = CallbackSet()
for i in 1:(N-1)
    condition(u,t,integrator) = begin
        θ = @view u[1:N]
        xi, _ = bob_pos(θ[i], pivot_xs[i])
        xj, _ = bob_pos(θ[i+1], pivot_xs[i+1])
        return (xj - xi) - 2r  # positive when separated, zero when touching
    end
    # only fire when approaching: directional = -1 means trigger when function crosses zero going negative? 
    # We'll check approach in affect to be safe. Use ContinuousCallback.
    cb = ContinuousCallback(condition,
        (integrator)->affect_pair!(integrator, i);
        rootfind = true, interp_neg = false, interp_pos = false)
    push!(callbacks.callbacks, cb)
end

# Solve
tspan = (0.0, 6.0)  # seconds
prob = ODEProblem(cradle_dynamics!, u0, tspan)
sol = solve(prob, Tsit5(), callback=callbacks, abstol=1e-8, reltol=1e-8)

# Prepare animation frames
anim = @animate for (ti, ui) in zip(sol.t, sol.u)
    θ = ui[1:N]
    xs = Float64[]
    ys = Float64[]
    for i in 1:N
        x, y = bob_pos(θ[i], pivot_xs[i])
        push!(xs, x)
        push!(ys, y)
    end

    # draw strings
    plt = scatter(xs, ys; xlims = (-1.2, 1.2), ylims = (-1.2, 0.2),
        markersize = 0, leg=false, aspect_ratio = :equal, xlabel = "", ylabel = "",
        title = @sprintf("Newton's cradle t = %.2f s", ti))

    # Draw strings and bobs
    for i in 1:N
        # line from pivot to bob
        plot!([pivot_xs[i], xs[i]], [0.0, ys[i]]; linewidth = 1.5)
        # bob circle: plot with scatter marker
        scatter!([xs[i]], [ys[i]]; markersize = 18, markerstrokecolor = :black)
    end

    # optional: baseline
    hline!( [-L-0.1], alpha=0)  # no-op to keep axes consistent
    plt
end fps=30

# Save gif
gif(anim, "newtons_cradle.gif", fps = 30)
println("Saved animation to newtons_cradle.gif")
