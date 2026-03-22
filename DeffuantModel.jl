module DeffuantModel

using Plots
using Graphs
using GraphRecipes
using NetworkLayout
using Random
using Statistics

export make_lattice, simulate, scan, plot_evolution, plot_snapshot, plot_scan, run_animation



"""
    make_lattice(type, L; periodic=true) → (g, x_coords, y_coords)

Build an L×L lattice of the requested type.
"""
function make_lattice(type::Symbol, L::Int; periodic::Bool=true)
    N = L * L
    g = SimpleGraph(N)
    xs = zeros(Float64, N)
    ys = zeros(Float64, N)

    idx(r, c) = (r - 1) * L + c

    offsets_fn = _offsets(type)
    coord_fn = _coords(type)

    for r in 1:L, c in 1:L
        u = idx(r, c)
        xs[u], ys[u] = coord_fn(r, c, L)

        for (dr, dc) in offsets_fn(r, c)
            nr, nc = r + dr, c + dc

            if periodic
                nr = mod1(nr, L)
                nc = mod1(nc, L)
            else
                (1 ≤ nr ≤ L && 1 ≤ nc ≤ L) || continue
            end

            nr == r && nc == c && continue
            add_edge!(g, u, idx(nr, nc))
        end
    end

    return g, xs, ys
end

# ── Neighbour offset rules ──────────────────────────────────

function _offsets(type::Symbol)
    if type == :triangular      # degree 3 — honeycomb bipartite
        return (r, c) -> iseven(r + c) ?
                         [(1, 0), (0, 1), (0, -1)] :
                         [(-1, 0), (0, 1), (0, -1)]

    elseif type == :square      # degree 4 — Von Neumann
        return (r, c) -> [(1, 0), (-1, 0), (0, 1), (0, -1)]

    elseif type == :hexagonal   # degree 6 — triangular lattice
        return (r, c) -> [(1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)]

    elseif type == :octagonal   # degree 8 — Moore
        return (r, c) -> [(1, 0), (-1, 0), (0, 1), (0, -1),
            (1, 1), (1, -1), (-1, 1), (-1, -1)]

    elseif type == :decagonal   # degree 10 — Moore + 2 axis
        return (r, c) -> [(1, 0), (-1, 0), (0, 1), (0, -1),
            (1, 1), (1, -1), (-1, 1), (-1, -1),
            (2, 0), (-2, 0)]

    elseif type == :dodecagonal # degree 12 — Moore + 4 axis
        return (r, c) -> [(1, 0), (-1, 0), (0, 1), (0, -1),
            (1, 1), (1, -1), (-1, 1), (-1, -1),
            (2, 0), (-2, 0), (0, 2), (0, -2)]

    else
        error("Unknown lattice: $type. Choose: " *
              ":triangular, :square, :hexagonal, :octagonal, :decagonal, :dodecagonal")
    end
end

# ── Visual coordinate helpers ───────────────────────────────

function _coords(type::Symbol)
    if type in (:triangular, :hexagonal)
        # Staggered rows for visual clarity
        return (r, c, L) -> begin
            cx = Float64(c) + (iseven(r) ? 0.5 : 0.0)
            cy = Float64(r) * (sqrt(3) / 2)
            (cx, cy)
        end
    else
        return (r, c, L) -> (Float64(c), Float64(r))
    end
end

# ============================================================
#  SIMULATION  (static & adaptive unified)
# ============================================================

"""
    simulate(g, MCS, ϵ, μ; kwargs...) → (history, g_history)

Run the Deffuant model on graph `g` for `MCS` Monte-Carlo steps.

Keyword arguments:
- `adaptive`  : enable co-evolutionary link rewiring (default `false`)
- `w`         : rewiring probability (used when adaptive=true, default 0.5)
- `random_μ`  : draw convergence param uniformly from [0.0, 0.5] (default `true`)
- `extremists`: pin nodes 1 and N as fixed extremists (default `false`)
- `save_interval`: graph snapshot frequency for adaptive runs (default 10)

Returns:
- `history`   : (MCS+1) × N matrix of opinions at each step
- `g_history` : vector of graph snapshots (empty for static runs)
"""
function simulate(g::AbstractGraph, MCS::Int, ϵ::Float64, μ::Float64;
    adaptive::Bool=false,
    w::Float64=0.5,
    random_μ::Bool=true,
    extremists::Bool=false,
    save_interval::Int=10)

    N = nv(g)
    g_cur = adaptive ? copy(g) : g

    opinions = rand(N) .* 2.0 .- 1.0
    if extremists
        opinions[1] = 0.95
        opinions[N] = -0.95
    end

    history = zeros(Float64, MCS + 1, N)
    history[1, :] .= opinions

    g_history = adaptive ? [copy(g_cur)] : typeof(g_cur)[]

    for t in 1:MCS
        for _ in 1:N
            i = rand(1:N)
            nbrs = neighbors(g_cur, i)
            isempty(nbrs) && continue
            j = rand(nbrs)
            updated = _step!(opinions, N, i, j, ϵ, μ, random_μ, extremists)

            if adaptive && !updated && rand() < w
                for _ in 1:50
                    k = rand(1:N)
                    if k ≠ i && k ≠ j && !has_edge(g_cur, i, k) &&
                       abs(opinions[i] - opinions[k]) ≤ ϵ
                        rem_edge!(g_cur, i, j)
                        add_edge!(g_cur, i, k)
                        break
                    end
                end
            end
        end

        history[t+1, :] .= opinions
        if adaptive && t % save_interval == 0
            push!(g_history, copy(g_cur))
        end
    end

    return history, g_history
end

function _step!(opinions, N, i, j, ϵ, μ, random_μ, extremists)
    abs(opinions[i] - opinions[j]) ≥ ϵ && return false
    diff = opinions[j] - opinions[i]
    μ_eff = random_μ ? (rand() * 0.5) : μ

    if extremists && (i == 1 || i == N)
        opinions[j] -= μ_eff * diff   # extremist pulls neighbour
    elseif extremists && (j == 1 || j == N)
        opinions[i] += μ_eff * diff
    else
        opinions[i] += μ_eff * diff
        opinions[j] -= μ_eff * diff
    end
    return true
end

# ============================================================
#  THERMODYNAMIC SCAN
# ============================================================

"""
    scan(g, ϵ_range, MCS, μ, num_trials; kwargs...) → NamedTuple

Sweep over tolerance values and compute ensemble averages.
Returns `(ϵ, O, χ, sum_abs, n_clusters)`.
"""
function scan(g::AbstractGraph, ϵ_range, MCS::Int, μ::Float64, num_trials::Int;
    adaptive::Bool=false,
    w::Float64=0.5,
    random_μ::Bool=true,
    extremists::Bool=false,
    cluster_threshold::Float64=1e-2)

    N = nv(g)
    ne = length(ϵ_range)

    O_arr = zeros(ne)
    χ_arr = zeros(ne)
    sa_arr = zeros(ne)
    nc_arr = zeros(ne)

    for (i, ϵ) in enumerate(ϵ_range)
        O_list = zeros(num_trials)
        sa_list = zeros(num_trials)
        nc_list = zeros(num_trials)

        for trial in 1:num_trials
            hist, _ = simulate(g, MCS, ϵ, μ;
                adaptive, w, random_μ, extremists, save_interval=MCS)
            final_ops = hist[end, :]

            # --- cluster analysis on sorted opinions ---
            # Anchor-based clustering: fix the first opinion in each cluster
            # as the reference. If a subsequent opinion deviates from the
            # anchor by more than cluster_threshold, start a new cluster.
            sorted = sort(final_ops)
            sizes = Int[]
            anchor = sorted[1]
            cur = 1
            for k in 2:N
                if abs(sorted[k] - anchor) < cluster_threshold
                    cur += 1
                else
                    push!(sizes, cur)
                    anchor = sorted[k]
                    cur = 1
                end
            end
            push!(sizes, cur)

            O_list[trial] = maximum(sizes) / N
            sa_list[trial] = sum(abs, final_ops) / N
            p = sizes ./ N
            nc_list[trial] = 1.0 / sum(p .^ 2)
        end

        O_mean = mean(O_list)
        O2_mean = mean(O_list .^ 2)
        O_arr[i] = O_mean
        χ_arr[i] = N * (O2_mean - O_mean^2)
        sa_arr[i] = mean(sa_list)
        nc_arr[i] = mean(nc_list)
    end

    return (ϵ=collect(ϵ_range), O=O_arr, χ=χ_arr,
        sum_abs=sa_arr, n_clusters=nc_arr)
end

# ============================================================
#  PLOTTING
# ============================================================

"""
    plot_evolution(history; title, color_grad) → Plot

Line plot of all agent opinions over time, coloured by final state.
"""
function plot_evolution(history; title="Opinion Evolution", color_grad=:bwr, filename=nothing)
    MCS = size(history, 1) - 1
    N = size(history, 2)
    final_ops = history[end, :]

    p = plot(title=title, xlabel="Time (MCS)", ylabel="Opinion",
        ylims=(-1.05, 1.05), legend=false, grid=true,
        background_color=:black, foreground_color=:white)
    for i in 1:N
        cv = (final_ops[i] + 1.0) / 2.0
        plot!(p, 0:MCS, history[:, i];
            linecolor=cgrad(color_grad)[cv], linealpha=0.3, linewidth=1.2)
    end
    if !isnothing(filename)
        savefig(p, filename)
        println("✓ Plot saved → $filename")
    end
    return p
end

"""
    plot_snapshot(g, history; x_coords, y_coords, t, title, color_grad) → Plot

Network plot coloured by opinions at time step `t` (default: last step).
Pass `x_coords` / `y_coords` from `make_lattice` for proper spatial layout.
"""
function plot_snapshot(g::AbstractGraph, history;
    x_coords=nothing,
    y_coords=nothing,
    t=nothing,
    title="Network Snapshot",
    color_grad=:bwr,
    filename=nothing)

    ts = isnothing(t) ? size(history, 1) : t
    opinions = history[ts, :]

    has_coords = !isnothing(x_coords) && !isnothing(y_coords)

    if has_coords
        p = graphplot(g; x=x_coords, y=y_coords, size=(550, 550),
            markersize=0.01, markeralpha=0.0,
            linecolor=:gray, linealpha=0.4,
            title=title, colorbar=false, axis=false, grid=false,
            background_color=:black, foreground_color=:white)
        msizes = fill(5.0, nv(g))
    else
        degs = Float64.(degree(g))
        msizes = 4.0 .+ 16.0 .* (degs ./ (maximum(degs) + 1e-6))
        p = graphplot(g; size=(550, 550),
            markersize=0.01, markeralpha=0.0,
            linecolor=:gray, linealpha=0.1, method=:sfdp,
            title=title, colorbar=false, axis=false, grid=false,
            background_color=:black, foreground_color=:white)
    end

    ns = p.series_list[end]
    scatter!(p, ns[:x], ns[:y];
        markersize=msizes, marker_z=opinions, markercolor=color_grad,
        clims=(-1.0, 1.0), markershape=:circle,
        markerstrokewidth=0.4, markerstrokecolor=:gray40, label="")
    if !isnothing(filename)
        savefig(p, filename)
        println("✓ Plot saved → $filename")
    end
    return p
end

"""
    plot_scan(result; N, title_prefix) → Plot

Four-panel diagnostic: order parameter, susceptibility, Σ|opinions|, effective clusters.
`result` is the NamedTuple returned by `scan()`.
"""
function plot_scan(result; N=nothing, title_prefix="", filename=nothing)
    ϵ = result.ϵ
    ns = isnothing(N) ? "" : " (N=$N)"

    p1 = plot(ϵ, result.O;
        marker=:circle, color=:steelblue,
        title="$(title_prefix)Phase Diagram$ns",
        ylabel="Order Parameter", legend=false, lw=2,
        background_color=:black, foreground_color=:white)
    p2 = plot(ϵ, result.χ;
        marker=:square, color=:firebrick,
        ylabel="Susceptibility χ", legend=false, lw=2,
        background_color=:black, foreground_color=:white)
    p3 = plot(ϵ, result.sum_abs;
        marker=:hexagon, color=:darkorange,
        ylabel="Σ|oᵢ|/N", legend=false, lw=2,
        background_color=:black, foreground_color=:white)
    p4 = plot(ϵ, result.n_clusters;
        marker=:star5, color=:mediumpurple,
        xlabel="Tolerance ϵ", ylabel="Effective Clusters",
        legend=false, lw=2,
        series_annotations=text.(round.(result.n_clusters, digits=2), 9, :mediumpurple, :bottom),
        background_color=:black, foreground_color=:white)

    p_final = plot(p1, p2, p3, p4; layout=(4, 1), size=(1000, 1600),
        left_margin=8Plots.mm, bottom_margin=5Plots.mm,
        background_color=:black, foreground_color=:white)
    if !isnothing(filename)
        savefig(p_final, filename)
        println("✓ Plot saved → $filename")
    end
    return p_final
end

# ============================================================
#  ANIMATION
# ============================================================

"""
    run_animation(g, history, g_history; x_coords, y_coords, fps, filename, step, color_grad)

Produce an animated GIF of opinion dynamics.
Works for **both** static (fixed network) and adaptive (co-evolving network) runs.

Arguments:
- `g`         : original graph (used when `g_history` is empty)
- `history`   : (MCS+1) × N opinion matrix from `simulate()`
- `g_history` : vector of graph snapshots (empty vector → static run)
- `x_coords`  : x positions from `make_lattice()` (optional but recommended)
- `y_coords`  : y positions from `make_lattice()`
- `fps`       : frames per second (default 15)
- `filename`  : output file path (default "opinion_anim.gif")
- `step`      : only animate every `step`-th time-step (default 5)
- `color_grad`: opinion colour gradient (default :bwr)
"""
function run_animation(g::AbstractGraph, history, g_history=nothing;
    x_coords=nothing,
    y_coords=nothing,
    fps::Int=7,
    filename::String="opinion_anim.gif",
    step::Int=5,
    color_grad=:bwr)

    MCS = size(history, 1) - 1
    is_adaptive = !isnothing(g_history) && length(g_history) > 1
    has_coords = !isnothing(x_coords) && !isnothing(y_coords)

    # ── STATIC: fixed lattice coords, iterate over every `step` time-steps ──
    if !is_adaptive
        # Use degree-based sizing only for random/irregular graphs (no provided coords)
        # Regular lattices (with provided coords) keep a uniform size (e.g. 7.0)
        use_degree_sizing = !has_coords

        # If coordinates are missing, compute them once to ensure a fixed orientation
        if !has_coords
            layout = NetworkLayout.sfdp(g)
            x_coords = [p[1] for p in layout]
            y_coords = [p[2] for p in layout]
            has_coords = true
        end

        if use_degree_sizing
            node_degrees = Float64.(degree(g))
            msizes = 2.0 .+ 8.0 .* (node_degrees ./ (maximum(node_degrees) + 1e-6))
            actual_markersize = 2.0 .* msizes
        else
            actual_markersize = 7.0
        end

        frames = 1:step:(MCS+1)
        anim = @animate for t in frames
            opinions = history[t, :]
            time_mcs = t - 1

            # We now always have coords (either passed in or computed above)
            p = graphplot(g; x=x_coords, y=y_coords, size=(600, 570),
                markersize=0.01, markeralpha=0.0,
                linecolor=:gray70, linealpha=0.4,
                title="t = $time_mcs MCS",
                colorbar=false, axis=false, grid=false,
                background_color=:white, foreground_color=:black)

            ns = p.series_list[end]
            scatter!(p, ns[:x], ns[:y];
                markersize=actual_markersize, marker_z=opinions, markercolor=color_grad,
                clims=(-1.0, 1.0), markershape=:circle,
                markerstrokewidth=0.4, markerstrokecolor=:black, label="")
        end

        # ── ADAPTIVE: continuous spring layout to prevent inter-frame rotation ──
    else
        num_snaps = length(g_history)

        # Initialise layout from the very first snapshot
        current_layout = NetworkLayout.spring(g_history[1])

        anim = @animate for snap_idx in 1:num_snaps
            g_snap = g_history[snap_idx]
            t_step = (snap_idx - 1) * max(1, MCS ÷ max(1, num_snaps - 1)) + 1
            opinions = history[min(t_step, size(history, 1)), :]
            time_mcs = t_step - 1

            # Update layout continuously — small local adjustments only,
            # no global re-orientation thanks to initialpos
            current_layout = NetworkLayout.spring(g_snap; initialpos=current_layout)

            x_cur = [p[1] for p in current_layout]
            y_cur = [p[2] for p in current_layout]

            node_degrees = Float64.(degree(g_snap))
            msizes = 2.0 .+ 8.0 .* (node_degrees ./ (maximum(node_degrees) + 1e-6))

            p = graphplot(g_snap; x=x_cur, y=y_cur, size=(600, 600),
                markersize=0.01, markeralpha=0.0,
                linecolor=:gray, linealpha=0.15,
                title="t = $time_mcs MCS",
                colorbar=false, axis=false, grid=false,
                background_color=:white, foreground_color=:black)

            scatter!(p, x_cur, y_cur;
                markersize=2.0 .* msizes, marker_z=opinions, markercolor=color_grad,
                clims=(-1.0, 1.0), markershape=:circle,
                markerstrokewidth=0.5, markerstrokecolor=:black,
                markeralpha=1.0, label="")
        end
    end

    gif(anim, filename; fps=fps)
    println("✓ Animation saved → $filename")
    return filename
end


end  # module DeffuantModel
