using Random
using Statistics
using Plots

"""
Simulates the baseline 1D Deffuant-Weisbuch opinion dynamics model 
and tracks the full time evolution for plotting.
"""
function simulate_and_track_deffuant(N::Int, MCS::Int, ϵ::Float64, μ::Float64)
    # Initialize opinions uniformly between -1.0 and 1.0
    opinions = rand(N) .* 2.0 .- 1.0

    # Pre-allocate a matrix to store the history: (MCS+1) rows, N columns
    history = zeros(Float64, MCS + 1, N)
    history[1, :] .= opinions

    # Random Sequential Updating
    for t in 1:MCS
        for _ in 1:N
            i = rand(1:N)
            j = rand(1:N)
            while i == j
                j = rand(1:N)
            end

            # Evaluate Bounded Confidence
            if abs(opinions[i] - opinions[j]) < ϵ
                diff = opinions[j] - opinions[i]
                opinions[i] += μ * diff
                opinions[j] -= μ * diff
            end
        end
        # Save the macroscopic state after exactly 1 full MCS
        history[t+1, :] .= opinions
    end

    return history
end

# 1. Run the Simulation
N = 500
MCS = 500 # 50 MCS is usually enough to see 1D convergence
ϵ = 0.2
μ = 0.5
history = simulate_and_track_deffuant(N, MCS, ϵ, μ)

# 2. Plotting the Time Evolution
# We extract the final opinions to color-code the lines based on their absorbing cluster
final_opinions = history[end, :]

# Create a base plot
plt = plot(
    title="Time Evolution of Opinions (1D Deffuant-Weisbuch)",
    xlabel="Time (Monte Carlo Steps)",
    ylabel="Opinion (x)",
    ylims=(-1.05, 1.05),
    legend=false,
    grid=true,
    size=(800, 500)
)

# Plot each agent's trajectory over time
# We map the final opinion to a color gradient (e.g., :viridis) for visual clarity
for i in 1:N
    # Normalize final opinion to a [0, 1] range for the color palette
    color_val = (final_opinions[i] + 1.0) / 2.0
    plot!(
        plt,
        0:MCS,
        history[:, i],
        linecolor=cgrad(:viridis)[color_val],
        linealpha=0.3, # Transparency to see density
        linewidth=1.5
    )
end

# Display the plot
display(plt)

# Optional: Save the plot to your working directory
# savefig(plt, "deffuant_evolution.png")

using Random
using Statistics
using Plots

"""
Simulates the Deffuant model on a 2D square lattice with a Moore neighborhood,
tracking the opinions of all agents over Monte Carlo steps.
"""
function simulate_and_track_2d_deffuant(L::Int, MCS::Int, ϵ::Float64, μ::Float64)
    N = L * L
    # Initialize L x L grid with uniform random opinions between -1.0 and 1.0
    grid = rand(L, L) .* 2.0 .- 1.0

    # History: (MCS+1) rows, N columns
    # We flatten the L x L grid into a vector of size N for each step
    history = zeros(Float64, MCS + 1, N)
    history[1, :] .= reshape(grid, :)

    # Random Sequential Updating
    for t in 1:MCS
        for _ in 1:N
            # 1. Select a random agent (i, j)
            i = rand(1:L)
            j = rand(1:L)

            # 2. Select a random neighbor from the Moore neighborhood (8 connections)
            dx, dy = 0, 0
            while dx == 0 && dy == 0
                dx = rand(-1:1)
                dy = rand(-1:1)
            end

            # Apply periodic boundary conditions using mod1
            ni = mod1(i + dx, L)
            nj = mod1(j + dy, L)

            # 3. Evaluate Bounded Confidence
            if abs(grid[i, j] - grid[ni, nj]) < ϵ
                diff = grid[ni, nj] - grid[i, j]
                grid[i, j] += μ * diff
                grid[ni, nj] -= μ * diff
            end
        end
        # Save the flattened state after each full MCS
        history[t+1, :] .= reshape(grid, :)
    end

    return history
end

# 1. Run the Simulation
L = 30          # 30x30 grid = 900 agents (smaller size for clearer plot)
MCS = 5000
ϵ = 0.2
μ = 0.5

println("Starting simulation...")
history = simulate_and_track_2d_deffuant(L, MCS, ϵ, μ)
println("Simulation completed.")

# 2. Plotting the Time Evolution
final_opinions = history[end, :]
N = L * L

plt = plot(
    title="Time Evolution of Opinions (2D Grid Deffuant Model)",
    xlabel="Time (Monte Carlo Steps)",
    ylabel="Opinion (x)",
    ylims=(-1.05, 1.05),
    legend=false,
    grid=true,
    size=(900, 600)
)

println("Generating plot...")
# Plot each agent's trajectory over time
for i in 1:N
    # Normalize final opinion to a [0, 1] range for the color palette
    color_val = (final_opinions[i] + 1.0) / 2.0
    plot!(
        plt,
        0:MCS,
        history[:, i],
        linecolor=cgrad(:viridis)[color_val],
        linealpha=0.2, # Transparency to see density
        linewidth=1.2
    )
end

display(plt)
println("Plot displayed.")

# Optional: Save the plot
# savefig(plt, "2d_evolution_plot.png")

