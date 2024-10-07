import sqlite3
import pandas as pd
import logging
import sys
import os
from glv_functions import run_glv_simulation

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(message)s')

# Define the directory to store simulation results
RESULTS_DIR = "simulation_results"

def simulate_metagenomic_data(simulation_id):
    # Ensure the results directory exists
    if not os.path.exists(RESULTS_DIR):
        os.makedirs(RESULTS_DIR)
        logging.info(f"Created directory: {RESULTS_DIR}")

    # Connect to SQLite database
    try:
        conn = sqlite3.connect('simulations.db')
        cursor = conn.cursor()

        # Validate that the simulation_id exists
        cursor.execute("SELECT COUNT(1) FROM simulations WHERE simulation_id = ?", (simulation_id,))
        if cursor.fetchone()[0] == 0:
            logging.error(f"Simulation ID {simulation_id} does not exist.")
            return

        # Retrieve parameter filenames for the given simulation_id
        cursor.execute("""
            SELECT growth_rates_filename, interaction_matrix_filename, initial_abundances_filename 
            FROM simulation_parameters 
            WHERE simulation_id = ?
        """, (simulation_id,))
        result = cursor.fetchone()
        if not result:
            logging.error(f"No parameter files found for simulation ID {simulation_id}.")
            return

        growth_rates_file, interaction_matrix_file, initial_abundances_file = result

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
            growth_rates = pd.read_csv(growth_rates_file).values.flatten()
            interaction_matrix = pd.read_csv(interaction_matrix_file).values
            initial_abundances = pd.read_csv(initial_abundances_file).values.flatten()

            logging.info(f"Successfully read parameter files for simulation ID {simulation_id}.")
        except Exception as e:
            logging.error(f"Error reading parameter files: {e}")
            return

        # Define simulation parameters
        n_steps = 500  # Number of time steps
        dt = 0.1  # Time step for integration

        logging.info(f"Growth rates: {growth_rates}")
        logging.info(f"Interaction matrix: {interaction_matrix}")
        logging.info(f"Initial abundances: {initial_abundances}")

        # Run the GLV simulation
        try:
            simulation_results = run_glv_simulation(growth_rates, interaction_matrix, initial_abundances, n_steps, dt)
        except Exception as e:
            logging.error(f"Error during GLV simulation: {e}")
            return

        # Create a DataFrame for simulation results
        time_points = pd.Series([dt * step for step in range(n_steps)], name="Time")
        df_results = pd.DataFrame(simulation_results)
        df_results.insert(0, 'Time', time_points)

        # Define the path for saving the simulation results
        results_filename = os.path.join(RESULTS_DIR, f'simulated_metagenomic_data_{simulation_id}.csv')

        # Save simulation results to CSV
        try:
            df_results.to_csv(results_filename, index=False)
            logging.info(f"Simulation data has been saved to {results_filename}.")
        except Exception as e:
            logging.error(f"Failed to save CSV: {e}")
            return

        # Save results metadata into the simulation_results table
        method_label = "sundials"
        plot_filename = None  # Placeholder if plot is generated in the future

        try:
            cursor.execute("""
            INSERT INTO simulation_results (simulation_id, method_label, results_filename, plot_filename, dt, n_steps)
            VALUES (?, ?, ?, ?, ?, ?)
            """, (simulation_id, method_label, results_filename, plot_filename, dt, n_steps))
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
        cursor.close()
        conn.close()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python simulate_metagenomic_data.py <simulation_id>")
        sys.exit(1)

    # Get the simulation_id from the command-line argument
    simulation_id = sys.argv[1]

    # Run the simulation
    simulate_metagenomic_data(simulation_id)