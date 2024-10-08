import numpy as np
import pandas as pd
import os
import sys
from scipy.integrate import odeint
from scipy.optimize import least_squares

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

# Define the generalized Lotka-Volterra model
def glv_model(X, t, alpha, beta):
    N = len(X)
    dXdt = np.zeros(N)
    for i in range(N):
        interaction = np.sum(beta[i] * X)
        dXdt[i] = X[i] * (alpha[i] + interaction)
    return dXdt

# Define a function to compute residuals for least squares fitting
def residuals(params, t, species_data):
    N = species_data.shape[1]
    alpha = params[:N]
    beta = params[N:].reshape((N, N))
    
    # Integrate the GLV system
    X0 = species_data[0, :]  # Initial conditions from the data
    sol = odeint(glv_model, X0, t, args=(alpha, beta))
    
    # Compute the difference between the model and the observed data
    return (sol - species_data).ravel()

# Initialize parameters (alpha and beta) for fitting
num_species = species_data.shape[1]
initial_alpha = np.zeros(num_species)  # Assume no initial growth rate bias
initial_beta = np.random.uniform(-0.01, 0.01, (num_species, num_species))  # Small interaction terms
initial_params = np.concatenate([initial_alpha, initial_beta.ravel()])

# Perform least squares fitting
result = least_squares(residuals, initial_params, args=(time_points, species_data), method='trf')

# Extract the fitted parameters
fitted_alpha = result.x[:num_species]
fitted_beta = result.x[num_species:].reshape((num_species, num_species))

# Prepare the least_squares folder for saving results
output_dir = 'least_squares'
os.makedirs(output_dir, exist_ok=True)

# Prepare the simulation_params folder for saving initial conditions and parameters
params_dir = 'simulation_params'
os.makedirs(params_dir, exist_ok=True)

# Define filenames using the pattern "simulation_id + '_LS1_' + original_name"
initial_abundances_filename = f"{simulation_id}LS1_initial_abundances.csv"
growth_rate_filename = f"{simulation_id}LS1_growth_rates.csv"
interaction_matrix_filename = f"{simulation_id}LS1_interaction_matrix.csv"
alpha_filename = f"{simulation_id}LS1_fitted_alpha.txt"
beta_filename = f"{simulation_id}LS1_fitted_beta.txt"
simulated_filename = f"{simulation_id}LS1_{('_'.join(input_filename.split('_')[1:]))[0:-4]}_least_squares.csv"

# Save the initial abundances as CSV with the format "Species", "Initial Abundance"
initial_abundances_path = os.path.join(params_dir, initial_abundances_filename)
initial_abundances_df = pd.DataFrame({
    'Species': [f'sp{i+1}' for i in range(num_species)],
    'Initial Abundance': species_data[0, :]
})
initial_abundances_df.to_csv(initial_abundances_path, index=False)

# Save the fitted alpha (growth rate) as CSV with the format "Species", "Growth rate"
growth_rate_path = os.path.join(params_dir, growth_rate_filename)
growth_rate_df = pd.DataFrame({
    'Species': [f'sp{i+1}' for i in range(num_species)],
    'Growth Rate': fitted_alpha
})
growth_rate_df.to_csv(growth_rate_path, index=False)

# Save the fitted beta (interaction matrix) as CSV with species labels as headers and the first cell empty
interaction_matrix_path = os.path.join(params_dir, interaction_matrix_filename)
interaction_matrix_df = pd.DataFrame(fitted_beta, columns=[f'sp{i+1}' for i in range(num_species)])
interaction_matrix_df.insert(0, '', [f'sp{i+1}' for i in range(num_species)])
interaction_matrix_df.to_csv(interaction_matrix_path, index=False)

# Save the fitted alpha and beta parameters in text files for the least squares folder
alpha_output_path = os.path.join(output_dir, alpha_filename)
beta_output_path = os.path.join(output_dir, beta_filename)

np.savetxt(alpha_output_path, fitted_alpha, fmt='%.6f', header='Fitted alpha (growth rates)')
np.savetxt(beta_output_path, fitted_beta, fmt='%.6f', header='Fitted beta (interaction coefficients)')

# Now simulate the system using the fitted parameters
def simulate_glv(t, initial_conditions, alpha, beta):
    return odeint(glv_model, initial_conditions, t, args=(alpha, beta))

# Use the initial conditions from the dataset and simulate
initial_conditions = species_data[0, :]
simulated_data = simulate_glv(time_points, initial_conditions, fitted_alpha, fitted_beta)

# Save the simulated data to the least_squares folder
simulated_output_path = os.path.join(output_dir, simulated_filename)
simulated_df = pd.DataFrame(simulated_data, columns=[f'sp{i+1}' for i in range(num_species)])
simulated_df.insert(0, 'Time', time_points)  # Insert time as the first column
simulated_df.to_csv(simulated_output_path, index=False)

# Notify the user that the results have been saved
print(f"Initial abundances saved to: {initial_abundances_path}")
print(f"Growth rates (alpha) saved to: {growth_rate_path}")
print(f"Interaction matrix (beta) saved to: {interaction_matrix_path}")
print(f"Fitted alpha saved to: {alpha_output_path}")
print(f"Fitted beta saved to: {beta_output_path}")
print(f"Simulated species data saved to: {simulated_output_path}")