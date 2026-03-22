include("DeffuantModel.jl")
using .DeffuantModel
using Graphs, Plots, Statistics, Random
gr()  # fast raster backend

# ==============================================================================
# 0. GLOBAL SETUP & BASE PARAMETERS
# ==============================================================================
# ── Global parameters ──────────────────────────────────────────
L = 40          # grid side → N = L² = 256 agents
MCS = 6000         # Monte-Carlo steps per simulation
μ = 0.5         # convergence parameter
w = 0.5         # rewiring probability (adaptive)
save_int = 1          # graph snapshot interval (adaptive)

# Scan parameters
ϵ_range = 0.15:0.025:0.65
MCS_scan = 1000
num_trials = 250

# Barabási–Albert parameters
N_ba = L * L       # same number of agents as lattice runs
k_ba = 4           # edges to attach per new node

# Epsilon values for evolution & animation
eps_vals = [0.3, 0.4, 0.5]

# Base seed for deterministic reproducibility across plot blocks and animation blocks
base_seed = 42

# Regular lattice types
lattice_types = [:triangular, :square, :hexagonal, :octagonal, :decagonal, :dodecagonal]
# Two lattice types for evolution/animation sections
evo_lattices = [:triangular, :hexagonal]

# Colour palette for lattice comparison plots
lattice_colors = [:steelblue, :firebrick, :darkgreen, :darkorange, :mediumpurple, :teal]

println("Setup complete. N = $N_ba agents per graph.")


# # ==============================================================================
# # 1. REGULAR LATTICES: THERMODYNAMIC SCANS
# # ==============================================================================
# L = 40
# μ = 0.5
# w = 0.5
# ϵ_range = 0.15:0.025:0.65
# MCS_scan = 10000
# num_trials = 200


# # ── Run thermodynamic scan for every regular lattice ──────────
# scan_regular = Dict{Symbol,NamedTuple}()

# for lt in lattice_types
#     println("Scanning $lt ...")
#     g, _, _ = make_lattice(lt, L)
#     scan_regular[lt] = scan(g, ϵ_range, MCS_scan, μ, num_trials; adaptive=false, w=w, random_μ=true, extremists=false, cluster_threshold=1e-2)
# end
# println("All regular lattice scans complete.")


# # ── Plot: all regular lattices on the same four-panel figure ──
# p_O = plot(title="Order Parameter", xlabel="Tolerance ϵ", ylabel="Order param.", legend=:bottomleft, lw=2)
# p_χ = plot(title="Susceptibility χ", xlabel="Tolerance ϵ", ylabel="χ", legend=:topright, lw=2)
# p_sa = plot(title="Σ|opinions|/N", xlabel="Tolerance ϵ", ylabel="Σ|oᵢ|/N", legend=:bottomleft, lw=2)
# p_nc = plot(title="Effective Clusters Neff", xlabel="Tolerance ϵ", ylabel="Neff", legend=:topright, lw=2)

# for (lt, col) in zip(lattice_types, lattice_colors)
#     r = scan_regular[lt]
#     g_tmp, _, _ = make_lattice(lt, L)
#     N_tmp = nv(g_tmp)
#     lbl = string(lt)
#     plot!(p_O, r.ϵ, r.O; label=lbl, color=col, marker=:circle, markersize=3)
#     plot!(p_χ, r.ϵ, r.χ; label=lbl, color=col, marker=:square, markersize=3)
#     plot!(p_sa, r.ϵ, r.sum_abs; label=lbl, color=col, marker=:diamond, markersize=3)
#     plot!(p_nc, r.ϵ, r.n_clusters; label=lbl, color=col, marker=:utriangle, markersize=3, series_annotations=text.(round.(r.n_clusters, digits=2), 6, :col, :bottom))
# end

# plot(p_O, p_χ;
#     layout=(2, 1), size=(1000, 1000), dpi=300,
#     plot_title="Regular Lattices — Thermodynamic Scan (no extremists)",
#     left_margin=6Plots.mm, bottom_margin=4Plots.mm)
# savefig("Thermodynamic_Scan_Regular_O_Chi_No_Extremists_size_$(L)_MCS_$(MCS_scan)_trials_$(num_trials).png")


# plot(p_sa, p_nc;
#     layout=(2, 1), size=(1000, 1000), dpi=300,
#     plot_title="Regular Lattices — Thermodynamic Scan (no extremists)",
#     left_margin=6Plots.mm, bottom_margin=4Plots.mm)
# savefig("Thermodynamic_Scan_Regular_SA_Neff_No_Extremists_size_$(L)_MCS_$(MCS_scan)_trials_$(num_trials).png")


# # ── Scan all regular lattices with extremists ──────────────────
# scan_regular_ext = Dict{Symbol,NamedTuple}()

# for lt in lattice_types
#     println("Scanning $lt (extremists) ...")
#     g, _, _ = make_lattice(lt, L)
#     scan_regular_ext[lt] = scan(g, ϵ_range, MCS_scan, μ, num_trials; adaptive=false, w=w, random_μ=true, extremists=true, cluster_threshold=1e-2)
# end
# println("All regular lattice scans (extremists) complete.")


# # ── Plot: all regular lattices on the same four-panel figure ──
# p_O = plot(title="Order Parameter", xlabel="Tolerance ϵ", ylabel="Order param.", legend=:bottomleft, lw=2)
# p_χ = plot(title="Susceptibility χ", xlabel="Tolerance ϵ", ylabel="χ", legend=:topright, lw=2)
# p_sa = plot(title="Σ|opinions|/N", xlabel="Tolerance ϵ", ylabel="Σ|oᵢ|/N", legend=:bottomleft, lw=2)
# p_nc = plot(title="Effective Clusters Neff", xlabel="Tolerance ϵ", ylabel="Neff", legend=:topright, lw=2)

# for (lt, col) in zip(lattice_types, lattice_colors)
#     r = scan_regular_ext[lt]
#     g_tmp, _, _ = make_lattice(lt, L)
#     N_tmp = nv(g_tmp)
#     lbl = string(lt)
#     plot!(p_O, r.ϵ, r.O; label=lbl, color=col, marker=:circle, markersize=3)
#     plot!(p_χ, r.ϵ, r.χ; label=lbl, color=col, marker=:square, markersize=3)
#     plot!(p_sa, r.ϵ, r.sum_abs; label=lbl, color=col, marker=:diamond, markersize=3)
#     plot!(p_nc, r.ϵ, r.n_clusters; label=lbl, color=col, marker=:utriangle, markersize=3, series_annotations=text.(round.(r.n_clusters, digits=1), 6, :col, :bottom))
# end

# plot(p_O, p_χ;
#     layout=(2, 1), size=(1000, 1000), dpi=300,
#     plot_title="Regular Lattices — Thermodynamic Scan (extremists=true)",
#     left_margin=6Plots.mm, bottom_margin=4Plots.mm)
# savefig("Thermodynamic_Scan_Regular_O_Chi_Extremists_size_$(L)_MCS_$(MCS_scan)_trials_$(num_trials).png")



# plot(p_sa, p_nc;
#     layout=(2, 1), size=(1000, 1000), dpi=300,
#     plot_title="Regular Lattices — Thermodynamic Scan (extremists=true)",
#     left_margin=6Plots.mm, bottom_margin=4Plots.mm)
# savefig("Thermodynamic_Scan_Regular_SA_Neff_Extremists_size_$(L)_MCS_$(MCS_scan)_trials_$(num_trials).png")



# # ==============================================================================
# # 2. REGULAR LATTICES: OPINION EVOLUTION PLOTS
# # ==============================================================================
# L           = 40          # grid side → N = L² 
# MCS         = 4000         # Monte-Carlo steps per simulation
# μ           = 0.5         # convergence parameter
# w           = 0.5         # rewiring probability (adaptive)
# save_int    = 1          # graph snapshot interval (adaptive)


# # ── Opinion evolution: regular lattices (square & hexagonal) ──
# # One panel per (lattice, ϵ) — arranged as rows=lattice, cols=ϵ
# evo_plots_regular_ext = []

# for lt in evo_lattices
#     g_lt, _, _ = make_lattice(lt, L)
#     for ϵ in eps_vals
#         hist, _ = simulate(g_lt, MCS, ϵ, μ; extremists=false)
#         p = plot_evolution(hist; title="$(lt) | ϵ=$(ϵ) | no ext")
#         push!(evo_plots_regular_ext, p)
#     end
# end

# plot(evo_plots_regular_ext...;
#      layout=(length(evo_lattices), length(eps_vals)),
#      size=(1400, 900), dpi=300, titlefontsize=9, tickfontsize=7,
#      guidefontsize=8, legendfontsize=6,
#      plot_title="Opinion Evolution — Regular Lattices (no extremists)",
#      left_margin=5Plots.mm, bottom_margin=3Plots.mm)
# savefig("Opinion_Evolution_Regular_Lattices_No_Extremists_size_$(L)_MCS_$(MCS).png")


# # ── Opinion evolution: regular lattices — extremists ──────────
# evo_plots_regular_ext = []

# for lt in evo_lattices
#     g_lt, _, _ = make_lattice(lt, L)
#     for ϵ in eps_vals
#         hist, _ = simulate(g_lt, MCS, ϵ, μ; extremists=true)
#         p = plot_evolution(hist; title="$(lt) | ϵ=$(ϵ) | ext")
#         push!(evo_plots_regular_ext, p)
#     end
# end

# plot(evo_plots_regular_ext...;
#      layout=(length(evo_lattices), length(eps_vals)),
#      size=(1400, 900), dpi=300, titlefontsize=9, tickfontsize=7,
#      guidefontsize=8, legendfontsize=6,
#      plot_title="Opinion Evolution — Regular Lattices (extremists=true)",
#      left_margin=5Plots.mm, bottom_margin=3Plots.mm)
# savefig("Opinion_Evolution_Regular_Lattices_Extremists_size_$(L)_MCS_$(MCS).png")



# # ==============================================================================
# # 3. REGULAR LATTICES: ANIMATIONS
# # ==============================================================================
# L           = 40          # grid side → N = L² = 256 agents
# MCS         = 500         # Monte-Carlo steps per simulation
# μ           = 0.5         # convergence parameter
# w           = 0.5         # rewiring probability (adaptive)
# save_int    = 1          # graph snapshot interval (adaptive)


# anim_step = 1

# # ── Animations: regular lattices ─────────────────────────────


# for lt in evo_lattices
#     g_lt, x_lt, y_lt = make_lattice(lt, L)
#     for ϵ in eps_vals
#         hist, _ = simulate(g_lt, MCS, ϵ, μ)
#         fname = "anim_$(lt)_eps$(ϵ).gif"
#         run_animation(g_lt, hist;
#             x_coords=x_lt, y_coords=y_lt,
#             fps=7, step=anim_step,
#             filename=fname)
#     end
# end

# # ── Animations: regular lattices — extremists ─────────────────
# for lt in evo_lattices
#     g_lt, x_lt, y_lt = make_lattice(lt, L)
#     for ϵ in eps_vals
#         hist, _ = simulate(g_lt, MCS, ϵ, μ; extremists=true)
#         fname = "anim_$(lt)_eps$(ϵ)_ext.gif"
#         run_animation(g_lt, hist;
#             x_coords=x_lt, y_coords=y_lt,
#             fps=7, step=anim_step,
#             filename=fname)
#     end
# end



# # ==============================================================================
# # 4. BARABÁSI-ALBERT: THERMODYNAMIC SCANS (STATIC & ADAPTIVE)
# # ==============================================================================
# L = 40
# μ = 0.5
# w = 0.5
# ϵ_range = 0.15:0.025:0.65
# MCS_scan = 1000
# num_trials = 200


# #── Build a Barabási–Albert graph ─────────────────────────────
Random.seed!(base_seed)
g_ba = barabasi_albert(N_ba, k_ba)
println("BA graph: N=$(nv(g_ba)), edges=$(ne(g_ba)), mean degree=$(round(mean(degree(g_ba)), digits=2))")

# # Static scan
# println("Scanning BA static ...")
# scan_ba_static = scan(g_ba, ϵ_range, MCS_scan, μ, num_trials; adaptive=false, w=w, random_μ=true, extremists=false, cluster_threshold=1e-2)

# # Adaptive scan
# println("Scanning BA adaptive ...")
# scan_ba_adaptive = scan(g_ba, ϵ_range, MCS_scan, μ, num_trials; adaptive=true, w=w, random_μ=true, extremists=false, cluster_threshold=1e-2)

# println("BA scans complete.")


# # ── Plot: BA static vs adaptive ───────────────────────────────
# ba_colors = [:steelblue, :firebrick]
# ba_labels = ["BA static", "BA adaptive"]
# ba_results = [scan_ba_static, scan_ba_adaptive]

# q_O = plot(title="Order Parameter", xlabel="Tolerance ϵ", ylabel="Order param.", legend=:bottomleft, lw=2)
# q_χ = plot(title="Susceptibility χ", xlabel="Tolerance ϵ", ylabel="χ", legend=:topright, lw=2)
# q_sa = plot(title="Σ|opinions|/N", xlabel="Tolerance ϵ", ylabel="Σ|oᵢ|/N", legend=:bottomleft, lw=2)
# q_nc = plot(title="Effective Clusters Neff", xlabel="Tolerance ϵ", ylabel="Neff", legend=:topright, lw=2)

# N_ba_agents = nv(g_ba)

# for (r, lbl, col) in zip(ba_results, ba_labels, ba_colors)
#     plot!(q_O, r.ϵ, r.O; label=lbl, color=col, marker=:circle, markersize=3)
#     plot!(q_χ, r.ϵ, r.χ; label=lbl, color=col, marker=:square, markersize=3)
#     plot!(q_sa, r.ϵ, r.sum_abs; label=lbl, color=col, marker=:diamond, markersize=3)
#     plot!(q_nc, r.ϵ, r.n_clusters; label=lbl, color=col, marker=:utriangle, markersize=3, series_annotations=text.(round.(r.n_clusters, digits=1), 6, :col, :bottom))
# end

# plot(q_O, q_χ;
#     layout=(2, 1), size=(1000, 1000), dpi=300,
#     plot_title="Barabási–Albert — Static vs Adaptive (no extremists)",
#     left_margin=6Plots.mm, bottom_margin=4Plots.mm)
# savefig("Thermodynamic_Scan_BA_O_Chi_No_Extremists_size_$(L)_MCS_$(MCS_scan)_trials_$(num_trials).png")


# plot(q_sa, q_nc;
#     layout=(2, 1), size=(1000, 1000), dpi=300,
#     plot_title="Barabasi-Albert — Thermodynamic Scan (no extremists)",
#     left_margin=6Plots.mm, bottom_margin=4Plots.mm)
# savefig("Thermodynamic_Scan_BA_SA_Neff_No_Extremists_size_$(L)_MCS_$(MCS_scan)_trials_$(num_trials).png")



# MCS_scan = 20000
# println("Scanning BA static (extremists) ...")
# scan_ba_static_ext = scan(g_ba, ϵ_range, MCS_scan, μ, num_trials; adaptive=false, w=w, random_μ=true, extremists=true, cluster_threshold=1e-2)

# println("Scanning BA adaptive (extremists) ...")
# scan_ba_adaptive_ext = scan(g_ba, ϵ_range, MCS_scan, μ, num_trials; adaptive=true, w=w, random_μ=true, extremists=true, cluster_threshold=1e-2)

# println("BA scans (extremists) complete.")


# # ── Plot: BA static vs adaptive ───────────────────────────────
# ba_colors = [:steelblue, :firebrick]
# ba_labels = ["BA static", "BA adaptive"]
# ba_results = [scan_ba_static_ext, scan_ba_adaptive_ext]

# q_O = plot(title="Order Parameter", xlabel="Tolerance ϵ", ylabel="Order param.", legend=:bottomleft, lw=2)
# q_χ = plot(title="Susceptibility χ", xlabel="Tolerance ϵ", ylabel="χ", legend=:topright, lw=2)
# q_sa = plot(title="Σ|opinions|/N", xlabel="Tolerance ϵ", ylabel="Σ|oᵢ|/N", legend=:bottomleft, lw=2)
# q_nc = plot(title="Effective Clusters Neff", xlabel="Tolerance ϵ", ylabel="Neff", legend=:topright, lw=2)

# N_ba_agents = nv(g_ba)

# for (r, lbl, col) in zip(ba_results, ba_labels, ba_colors)
#     plot!(q_O, r.ϵ, r.O; label=lbl, color=col, marker=:circle, markersize=3)
#     plot!(q_χ, r.ϵ, r.χ; label=lbl, color=col, marker=:square, markersize=3)
#     plot!(q_sa, r.ϵ, r.sum_abs; label=lbl, color=col, marker=:diamond, markersize=3)
#     plot!(q_nc, r.ϵ, r.n_clusters; label=lbl, color=col, marker=:utriangle, markersize=3, series_annotations=text.(round.(r.n_clusters, digits=1), 6, :col, :bottom))
# end

# plot(q_O, q_χ;
#     layout=(2, 1), size=(1000, 1000), dpi=300,
#     plot_title="Barabási–Albert — Static vs Adaptive (extremists=true)",
#     left_margin=6Plots.mm, bottom_margin=4Plots.mm)
# savefig("Thermodynamic_Scan_BA_O_Chi_Extremists_size_$(L)_MCS_$(MCS_scan)_k_$(k_ba)_trials_$(num_trials).png")

# plot(q_sa, q_nc;
#     layout=(2, 1), size=(1000, 1000), dpi=300,
#     plot_title="Barabasi-Albert — Thermodynamic Scan (extremists=true)",
#     left_margin=6Plots.mm, bottom_margin=4Plots.mm)
# savefig("Thermodynamic_Scan_BA_SA_Neff_Extremists_size_$(L)_MCS_$(MCS_scan)_k_$(k_ba)_trials_$(num_trials).png")


# # ==============================================================================
# # 5. BARABÁSI-ALBERT: OPINION EVOLUTION PLOTS
# # ==============================================================================
# L = 40          # grid side → N = L² = 256 agents
# MCS = 2000         # Monte-Carlo steps per simulation
# μ = 0.5         # convergence parameter
# w = 0.5         # rewiring probability (adaptive)
# save_int = 1          # graph snapshot interval (adaptive)



# # ── Opinion evolution: BA static ──────────────────────────────
# evo_plots_ba_static = []

# for ϵ in eps_vals
#     Random.seed!(base_seed + round(Int, ϵ * 100))
#     hist_bas, _ = simulate(g_ba, MCS, ϵ, μ; adaptive=false)
#     p = plot_evolution(hist_bas; title="BA static | ϵ=$(ϵ)")
#     push!(evo_plots_ba_static, p)
# end

# plot(evo_plots_ba_static...;
#     layout=(1, length(eps_vals)),
#     size=(1100, 350), dpi=300, titlefontsize=9, tickfontsize=7,
#     guidefontsize=8, legendfontsize=6,
#     plot_title="Opinion Evolution — BA Static (no extremists)",
#     left_margin=5Plots.mm, bottom_margin=3Plots.mm)
# savefig("Opinion_Evolution_BA_Static_No_Extremists_size_$(L)_MCS_$(MCS)_k_$(k_ba).png")


# # ── Opinion evolution: BA adaptive ────────────────────────────
# evo_plots_ba_adaptive = []

# for ϵ in eps_vals
#     Random.seed!(base_seed + round(Int, ϵ * 100))
#     hist_baa, g_snaps_baa = simulate(g_ba, MCS, ϵ, μ; adaptive=true, w=w, save_interval=save_int)
#     p = plot_evolution(hist_baa; title="BA adaptive | ϵ=$(ϵ)")
#     push!(evo_plots_ba_adaptive, p)
# end

# plot(evo_plots_ba_adaptive...;
#     layout=(1, length(eps_vals)),
#     size=(1100, 350), dpi=300,
#     plot_title="Opinion Evolution — BA Adaptive (no extremists)",
#     left_margin=5Plots.mm, bottom_margin=3Plots.mm)
# savefig("Opinion_Evolution_BA_Adaptive_No_Extremists_size_$(L)_MCS_$(MCS)_k_$(k_ba).png")

# # ── Opinion evolution: BA static — extremists ─────────────────
# evo_plots_ba_static_ext = []

# for ϵ in eps_vals
#     Random.seed!(base_seed + round(Int, ϵ * 100))
#     hist_bas, _ = simulate(g_ba, MCS, ϵ, μ; adaptive=false, extremists=true)
#     p = plot_evolution(hist_bas; title="BA static | ϵ=$(ϵ) | ext")
#     push!(evo_plots_ba_static_ext, p)
# end

# plot(evo_plots_ba_static_ext...;
#     layout=(1, length(eps_vals)),
#     size=(1100, 350), dpi=300,
#     plot_title="Opinion Evolution — BA Static (extremists=true)",
#     left_margin=5Plots.mm, bottom_margin=3Plots.mm)
# savefig("Opinion_Evolution_BA_Static_Extremists_size_$(L)_MCS_$(MCS)_k_$(k_ba).png")

# # ── Opinion evolution: BA adaptive — extremists ───────────────
# evo_plots_ba_adaptive_ext = []

# for ϵ in eps_vals
#     Random.seed!(base_seed + round(Int, ϵ * 100))
#     hist_baa, _ = simulate(g_ba, MCS, ϵ, μ; adaptive=true, w=w, extremists=true, save_interval=save_int)
#     p = plot_evolution(hist_baa; title="BA adaptive | ϵ=$(ϵ) | ext")
#     push!(evo_plots_ba_adaptive_ext, p)
# end

# plot(evo_plots_ba_adaptive_ext...;
#     layout=(1, length(eps_vals)),
#     size=(1100, 350), dpi=300,
#     plot_title="Opinion Evolution — BA Adaptive (extremists=true)",
#     left_margin=5Plots.mm, bottom_margin=3Plots.mm)
# savefig("Opinion_Evolution_BA_Adaptive_Extremists_size_$(L)_MCS_$(MCS)_k_$(k_ba).png")



# # ==============================================================================
# # 6. BARABÁSI-ALBERT: ANIMATIONS
# # ==============================================================================
# L = 40          # grid side → N = L² = 256 agents
# MCS = 500         # Monte-Carlo steps per simulation
# μ = 0.5         # convergence parameter
# w = 0.5         # rewiring probability (adaptive)
# save_int = 1          # graph snapshot interval (adaptive)

# anim_step = 1



# # ── Animations: BA static ─────────────────────────────────────
# for ϵ in eps_vals
#     Random.seed!(base_seed + round(Int, ϵ*100))
#     hist_bas, _ = simulate(g_ba, MCS, ϵ, μ; adaptive=false)
#     fname = "anim_BA_static_eps$(ϵ).gif"
#     run_animation(g_ba, hist_bas;
#         fps=7, step=anim_step,
#         filename=fname)
# end

# # ── Animations: BA adaptive ───────────────────────────────────
# for ϵ in eps_vals
#     Random.seed!(base_seed + round(Int, ϵ*100))
#     hist_baa, g_snaps_baa = simulate(g_ba, MCS, ϵ, μ;
#         adaptive=true, w=w, save_interval=save_int)
#     fname = "anim_BA_adaptive_eps$(ϵ).gif"
#     run_animation(g_ba, hist_baa, g_snaps_baa;
#         fps=10, step=anim_step,
#         filename=fname)
# end

# # ── Animations: BA static — extremists ───────────────────────
# for ϵ in eps_vals
#     Random.seed!(base_seed + round(Int, ϵ*100))
#     hist_bas, _ = simulate(g_ba, MCS, ϵ, μ; adaptive=false, extremists=true)
#     fname = "anim_BA_static_eps$(ϵ)_ext.gif"
#     run_animation(g_ba, hist_bas;
#         fps=7, step=anim_step,
#         filename=fname)
# end

# # ── Animations: BA adaptive — extremists ─────────────────────
# for ϵ in eps_vals
#     Random.seed!(base_seed + round(Int, ϵ*100))
#     hist_baa, g_snaps_baa = simulate(g_ba, MCS, ϵ, μ;
#         adaptive=true, w=w, extremists=true, save_interval=save_int)
#     fname = "anim_BA_adaptive_eps$(ϵ)_ext.gif"
#     run_animation(g_ba, hist_baa, g_snaps_baa;
#         fps=10, step=anim_step,
#         filename=fname)
# end