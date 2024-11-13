# Import necessary packages
using CSV, DataFrames
using DifferentialEquations
using Turing, Distributions
using Random
using Plots

# Load the data
data = CSV.read("osci4_sundials_0_250_250_simulation.csv", DataFrame)
time = data.Time
species_data = Matrix(data[:, Not(:Time)])  # Extract species data without the Time column

# Define the Generalized Lotka-Volterra (GLV) model
function glv!(du, u, p, t)
    α = p[1:4]  # growth rates
    β = reshape(p[5:end], 4, 4)  # interaction matrix reshaped for 4 species

    for i in 1:4
        du[i] = u[i] * (α[i] + sum(β[i, j] * u[j] for j in 1:4))
    end
end

# Set up initial conditions and time span
u0 = species_data[1, :]  # Initial conditions for the species
tspan = (time[1], time[end])

# Bayesian Model Definition with Turing
@model function bayesian_glv_model(data, time, u0)
    # Priors for growth rates (α) - Normal with mean 0 and std 1
    α ~ filldist(Normal(0, 1), 4)

    # Priors for interaction coefficients (β) - Normal with mean 0 and std 0.1
    β ~ filldist(Normal(0, 0.1), 16)  # 4x4 matrix flattened to a vector of length 16

    # Prior for noise term
    σ ~ InverseGamma(2, 3)

    # Combine α and β into a single parameter vector
    p = vcat(α, β)

    # Define the ODE problem using the GLV model with the current parameters
    prob = ODEProblem(glv!, u0, tspan, p)

    # Solve the ODE with a solver; use Tsit5 for stability
    sol = solve(prob, Tsit5(), saveat=time)

    # Convert the solution to a matrix
    sol_mat = Array(sol)

    # Likelihood: Gaussian with observed data and noise σ
    for t in 1:length(time)
        for i in 1:4
            data[t, i] ~ Normal(sol_mat[i, t], σ)
        end
    end
end

# Instantiate the model with the data
bayesian_model = bayesian_glv_model(species_data, time, u0)

# Set random seed for reproducibility
Random.seed!(1234)

# Sample from the posterior using NUTS (No-U-Turn Sampler)
chain = sample(bayesian_model, NUTS(0.65), 1000)

# Summary and Diagnostics
println(chain)

# Plot Posterior Distributions
plot(chain)

# Posterior Predictive Checks: Simulate with mean parameters
# Extract posterior means for α and β
α_est = Array(mean(chain[:α], dims=1))'
β_est = Array(mean(chain[:β], dims=1))'

# Simulate model with posterior mean parameters
p_estimated = vcat(α_est, β_est)
prob_estimated = ODEProblem(glv!, u0, tspan, p_estimated)
sol_estimated = solve(prob_estimated, Tsit5(), saveat=time)

# Plot results
plot(time, species_data, label=["sp1 data" "sp2 data" "sp3 data" "sp4 data"])
plot!(time, hcat(sol_estimated.u...)', linestyle=:dash, label=["sp1 est" "sp2 est" "sp3 est" "sp4 est"])