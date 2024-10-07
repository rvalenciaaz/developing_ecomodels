import pandas as pd
import os

# Directory where parameter files will be saved
PARAMS_DIR = 'simulation_params'
os.makedirs(PARAMS_DIR, exist_ok=True)  # Create the directory if it doesn't exist

# Example data for the CSV files

# 1. Growth rates
growth_rates_data = {
    "Species": ["sp1", "sp2", "sp3"],
    "Growth Rate": [0.5, 0.8, 0.6]
}
growth_rates_df = pd.DataFrame(growth_rates_data)
growth_rates_df.to_csv(os.path.join(PARAMS_DIR, "simulation1_growth_rates.csv"), index=False)

# 2. Interaction matrix (with row and column names)
interaction_matrix_data = {
    "sp1": [0.0, 0.1, -0.2],
    "sp2": [0.1, 0.0, 0.3],
    "sp3": [-0.2, 0.3, 0.0]
}
interaction_matrix_df = pd.DataFrame(interaction_matrix_data, index=["sp1", "sp2", "sp3"])
interaction_matrix_df.to_csv(os.path.join(PARAMS_DIR, "simulation1_interaction_matrix.csv"))

# 3. Initial abundances
initial_abundances_data = {
    "Species": ["sp1", "sp2", "sp3"],
    "Initial Abundance": [1000, 1500, 1200]
}
initial_abundances_df = pd.DataFrame(initial_abundances_data)
initial_abundances_df.to_csv(os.path.join(PARAMS_DIR, "simulation1_initial_abundances.csv"), index=False)

print(f"Files saved in directory: {PARAMS_DIR}")
