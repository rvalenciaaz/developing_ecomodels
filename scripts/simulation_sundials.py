import sqlite3
import pandas as pd
import logging
import sys
import os

from lib import glv_functions

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(message)s')

# Define the directories to store simulation results and parameter files
RESULTS_DIR = "simulation_results"
PARAMS_DIR = "simulation_params"

def simulate_metagenomic_data(simulation_id):
    # Ensure the results directory exists
    if not os.path.exists(RESULTS_DIR):
        os.makedirs(RESULTS_DIR)
        logging.info(f"Created directory: {RESULTS_DIR}")

    # Connect to SQLite database
    try:
        conn = sqlite3.connect('simulations.db')

        # Use Pandas to validate that the simulation_id exists
        df_simulations = pd.read_sql_query("SELECT * FROM simulations WHERE simulation_id = ?", conn, params=(simulation_id,))
        if df_simulations.empty:
            logging.error(f"Simulation ID {simulation_id} does not exist.")
            return

        # Retrieve parameter filenames for the given simulation_id using Pandas
        query = """
            SELECT growth_rates_filename, interaction_matrix_filename, initial_abundances_filename 
            FROM simulation_parameters 
            WHERE simulation_id = ?
        """
        df_params = pd.read_sql_query(query, conn, params=(simulation_id,))
        if df_params.empty:
            logging.error(f"No parameter files found for simulation ID {simulation_id}.")
            return

        growth_rates_file = os.path.join(PARAMS_DIR, df_params['growth_rates_filename'].values[0])
        interaction_matrix_file = os.path.join(PARAMS_DIR, df_params['interaction_matrix_filename'].values[0])
        initial_abundances_file = os.path.join(PARAMS_DIR, df_params['initial_abundances_filename'].values[0])

        # Check if the files exist
        if not os.path.exists(growth_rates_file):
            logging.error(f"Growth rates file not found: {growth_rates_file}")
            return
        if not os.path.exists(interaction_matrix_file):
            logging.error(f"Interaction matrix file not found: {interaction_matrix_file}")
            return
        if not os.path.exists(initial_abundances_file):
            logging.error(f"Initial abundances file not found: {initial_abundances_file}")
            return

        # Read the parameter files using pandas
        try:
            growth_rates = pd.read_csv(growth_rates_file)
            interaction_matrix = pd.read_csv(interaction_matrix_file)
            initial_abundances = pd.read_csv(initial_abundances_file)
            logging.info(f"Successfully read parameter files for simulation ID {simulation_id}.")
        except Exception as e:
            logging.error(f"Error reading parameter files: {e}")
            return

        growth_rates= growth_rates["Growth Rate"].values.flatten()
        interaction_matrix= interaction_matrix.iloc[:,1:].values
        initial_abundances= initial_abundances["Initial Abundance"].values.flatten()

        # Define simulation parameters
        n_steps = 500  # Number of time steps
        dt = 0.1  # Time step for integration

        logging.info(f"Growth rates: {growth_rates}")
        logging.info(f"Interaction matrix: {interaction_matrix}")
        logging.info(f"Initial abundances: {initial_abundances}")

        # Run the GLV simulation
        try:
            simulation_results = glv_functions.run_glv_simulation(growth_rates, interaction_matrix, initial_abundances, n_steps, dt)
        except Exception as e:
            logging.error(f"Error during GLV simulation: {e}")
            return

        # Create a DataFrame for simulation results
        time_points = pd.Series([dt * step for step in range(n_steps)], name="Time")
        df_results = pd.DataFrame(simulation_results)
        df_results.insert(0, 'Time', time_points)


        # Save results metadata into the simulation_results table
        method_label = "sundials"
        plot_filename = None  # Placeholder if plot is generated in the future

        # Define the path for saving the simulation results
        results_filename = os.path.join(RESULTS_DIR, f'{simulation_id}_{method_label}_simulation.csv')

        # Save simulation results to CSV
        try:
            df_results.to_csv(results_filename, index=False)
            logging.info(f"Simulation data has been saved to {results_filename}.")
        except Exception as e:
            logging.error(f"Failed to save CSV: {e}")
            return

        

        try:
            # Prepare the metadata as a DataFrame for insertion
            df_metadata = pd.DataFrame({
                'simulation_id': [simulation_id],
                'method_label': [method_label],
                'results_filename': [results_filename],
                'plot_filename': [plot_filename],
                'dt': [dt],
                'n_steps': [n_steps]
            })

            # Append metadata to the simulation_results table using pandas
            df_metadata.to_sql('simulation_results', conn, if_exists='append', index=False)
            logging.info(f"Inserted results metadata for simulation ID {simulation_id}.")
        except sqlite3.Error as e:
            logging.error(f"Error inserting simulation results metadata: {e}")
            return

        # Commit changes to the database
        conn.commit()

    except sqlite3.Error as e:
        logging.error(f"SQLite error: {e}")
    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}")
    finally:
        conn.close()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python simulate_metagenomic_data.py <simulation_id>")
        sys.exit(1)

    # Get the simulation_id from the command-line argument
    simulation_id = sys.argv[1]

    # Run the simulation
    simulate_metagenomic_data(simulation_id)