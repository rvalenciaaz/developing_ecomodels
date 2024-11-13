# Import necessary packages
using CSV, DataFrames
using DifferentialEquations
using Turing, Distributions
using Random
using Zygote  # For automatic differentiation
using Plots

# Load the csv data as a dataframe
data = CSV.read("osci4_sundials_0_250_250_simulation.csv", DataFrame)
time = data.Time
species_data = Matrix(data[:, Not(:Time)])  # Extract species data without the Time column

# Normalize species data, depends of what are you looking for?
#max_values = maximum(species_data, dims=1)
#species_data_scaled = species_data ./ max_values
#u0_scaled = species_data_scaled[1, :]

species_data_scaled = species_data
u0_scaled = species_data_scaled[1, :]

# Define the Generalized Lotka-Volterra (GLV) model
function glv!(du, u, p, t)
    α = p[1:4]  # Growth rates
    β = reshape(p[5:end], 4, 4)  # Interaction matrix reshaped for 4 species

    for i in 1:4
        du[i] = u[i] * (α[i] + sum(β[i, j] * u[j] for j in 1:4))
    end
end

# Set up time span
tspan = (time[1], time[end])

# Bayesian Model Definition with Turing
@model function bayesian_glv_model(data, time, u0)
    # Priors for growth rates (α) - LogNormal for positive rates
    α ~ filldist(LogNormal(0, 0.25), 4)
    
    # Priors for interaction coefficients (β)
    β ~ filldist(Normal(-1, 2), 16)
    
    # Prior for noise term
    σ ~ InverseGamma(2, 3)
    
    # Combine α and β into a single parameter vector
    p = vcat(α, β)
    
    # Define the ODE problem
    prob = ODEProblem(glv!, u0, tspan, p)
    
    # Solve the ODE with forced stepping at data times
    sol = solve(prob, Tsit5(), saveat=time, tstops=time)
    
    # Convert the solution to a matrix
    sol_mat = Array(sol)
    
    # Check for non-finite or negative values
    if any(isnan, sol_mat) || any(isinf, sol_mat) || any(sol_mat .< 0)
        Turing.@addlogprob!(-Inf)
        return
    end
    
    # Likelihood: Gaussian with observed data and noise σ
    for t in 1:length(time)
        for i in 1:4
            data[t, i] ~ Normal(sol_mat[i, t], σ)
        end
    end
end

# Instantiate the model with the scaled data
bayesian_model = bayesian_glv_model(species_data_scaled, time, u0_scaled)

# Set random seed for reproducibility
Random.seed!(1234)

# Adjust sampler settings
Turing.setadbackend(:zygote)
chain = sample(bayesian_model, NUTS(), 1000; num_warmup=1000, max_depth=15, init_ϵ=0.01)

# Summary and Diagnostics
println(chain)

# Plot Posterior Distributions
plot(chain)

# Posterior Predictive Checks: Simulate with mean parameters
α_est = Array(mean(chain[:α], dims=1))[:]
β_est = Array(mean(chain[:β], dims=1))[:]

# Simulate model with posterior mean parameters
p_estimated = vcat(α_est, β_est)
prob_estimated = ODEProblem(glv!, u0_scaled, tspan, p_estimated)
sol_estimated = solve(prob_estimated, Tsit5(), saveat=time, tstops=time)

# Rescale the solution
sol_estimated_mat = Array(sol_estimated) .* max_values'

# Plot results
plot(time, species_data, label=["sp1 data", "sp2 data", "sp3 data", "sp4 data"])
plot!(time, sol_estimated_mat', linestyle=:dash, label=["sp1 est", "sp2 est", "sp3 est", "sp4 est"])