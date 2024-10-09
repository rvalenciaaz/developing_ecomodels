import numpy as np
import pandas as pd
import os
import sys
from scipy.integrate import odeint
from botorch.optim import optimize_acqf
from botorch.models import SingleTaskGP
from botorch.acquisition import qUpperConfidenceBound
from botorch import fit_gpytorch_model
from gpytorch.mlls import ExactMarginalLogLikelihood
from torch.quasirandom import SobolEngine
import torch
import matplotlib.pyplot as plt

# Ensure proper usage
if len(sys.argv) != 2:
    print("Usage: python glv_model_fit.py <input_csv_file>")
    sys.exit(1)

# Get the input file from command-line arguments
input_filename = sys.argv[1]
input_file_path = os.path.join('simulation_results', input_filename)

# Ensure input file exists
if not os.path.exists(input_file_path):
    print(f"Error: File '{input_file_path}' not found.")
    sys.exit(1)

# Extract simulation_id from the input filename
simulation_id = input_filename.split('_')[0]

# Load the data
data = pd.read_csv(input_file_path)

# Prepare time points and species data
time_points = data['Time'].values
species_data = data.iloc[:, 1:].values  # Species data (sp1 to sp9)

#this looks ok

# Define the generalized Lotka-Volterra model
def glv_model(X, t, alpha, beta):
    N = len(X)
    dXdt = np.zeros(N)
    for i in range(N):
        interaction = np.sum(beta[i] * X)
        dXdt[i] = X[i] * (alpha[i] + interaction)
    return dXdt

# Define a function to compute the sum of squared residuals
def compute_residuals(params):
    N = species_data.shape[1]
    alpha = params[:N]
    beta = params[N:].reshape((N, N))
    
    # Integrate the GLV system
    X0 = species_data[0, :]  # Initial conditions from the data
    sol = odeint(glv_model, X0, time_points, args=(alpha, beta))
    
    # Compute the difference between the model and the observed data
    residuals = (sol - species_data).ravel()
    return np.sum(residuals**2)

# Set up the Bayesian optimization problem
num_species = species_data.shape[1]
initial_alpha = np.zeros(num_species, dtype=np.float64)  # Assume no initial growth rate bias
initial_beta = np.random.uniform(-0.01, 0.01, (num_species, num_species)).astype(np.float64)  # Small interaction terms

# Define the parameter space (alpha and beta) using float64 precision
bounds = torch.stack([torch.full((num_species + num_species**2,), -0.1, dtype=torch.float64),
                      torch.full((num_species + num_species**2,), 0.1, dtype=torch.float64)])

# Convert initial parameters to torch tensors with float64 precision
initial_params = torch.tensor(np.concatenate([initial_alpha, initial_beta.ravel()]), dtype=torch.float64)

# Sobol sequence for initial sampling (also in float64 precision)
sobol = SobolEngine(dimension=num_species + num_species**2, scramble=True)
initial_samples = sobol.draw(10, dtype=torch.float64) * (bounds[1] - bounds[0]) + bounds[0]

# Evaluate initial samples with float64 precision
Y = torch.tensor([compute_residuals(x.numpy()) for x in initial_samples], dtype=torch.float64).unsqueeze(-1)

# Fit a GP model using float64 precision
gp = SingleTaskGP(initial_samples, Y)
mll = ExactMarginalLogLikelihood(gp.likelihood, gp)
fit_gpytorch_model(mll)

# Define batch acquisition function for batch Bayesian optimization (using float64 precision)
batch_size = 3  # Number of candidates to evaluate in each iteration
q_ucb = qUpperConfidenceBound(gp, beta=0.1)

# Optimization loop
best_values = []
possible_minimums = []
possible_minimum_values = []

for i in range(30):  # 30 iterations of Bayesian optimization
    # Optimize the batch acquisition function
    candidates, acq_value = optimize_acqf(q_ucb, bounds=bounds, q=batch_size, num_restarts=5, raw_samples=20, dtype=torch.float64)
    
    # Evaluate the candidates with float64 precision
    new_y = torch.tensor([compute_residuals(candidate.detach().numpy()) for candidate in candidates], dtype=torch.float64).unsqueeze(-1)

    # Update the model with new observations
    gp = gp.condition_on_observations(candidates, new_y)
    
    # Refit the GP model
    mll = ExactMarginalLogLikelihood(gp.likelihood, gp)
    fit_gpytorch_model(mll)
    
    # Record the best value from the current batch
    best_value_in_batch = torch.min(new_y).item()
    best_values.append(best_value_in_batch)
    
    # Store the possible minimum parameters and their values
    for j, candidate in enumerate(candidates):
        possible_minimums.append(candidate.detach().numpy())
        possible_minimum_values.append(new_y[j].item())

# Plot the best value over iterations
plt.figure()
plt.plot(best_values, label='Best value')
plt.xlabel('Iteration')
plt.ylabel('Objective value (sum of squared residuals)')
plt.title('Batch Bayesian Optimization of GLV model (float64)')
plt.legend()
plt.show()

# Save the list of possible minimums and their values (using float64 precision)
possible_minimums_df = pd.DataFrame(possible_minimums, columns=[f'param_{i+1}' for i in range(num_species + num_species**2)])
possible_minimums_df['Objective Value'] = possible_minimum_values
possible_minimums_df.to_csv(f"{simulation_id}_possible_minimums_and_values_float64.csv", index=False)

print(f"Possible minimums and their values saved to: {simulation_id}_possible_minimums_and_values_float64.csv")