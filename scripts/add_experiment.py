import sqlite3
import pandas as pd
import os
import sys

# Directories where parameter files are saved
PARAMS_DIR = 'simulation_params'

def get_available_bacteria(cursor, species_ids):
    """Fetch bacterial species by their IDs from the database."""
    cursor.execute(f"""
    SELECT species_id, genus, species, strain 
    FROM bacterial_species 
    WHERE species_id IN ({','.join(['?']*len(species_ids))});
    """, species_ids)
    return cursor.fetchall()

def read_parameter_file(filepath):
    """Helper function to read a CSV file using pandas."""
    return pd.read_csv(filepath)

def verify_species_ids(growth_rates, interaction_matrix, initial_abundances):
    """Verify that the species IDs match across the different parameter files."""
    growth_species_ids = growth_rates.iloc[:, 0].tolist()
    interaction_species_ids = interaction_matrix.columns.tolist()[1:]
    abundance_species_ids = initial_abundances.iloc[:, 0].tolist()

    if not (growth_species_ids == interaction_species_ids == abundance_species_ids):
        raise ValueError("Species IDs mismatch between parameter files.")
    
    return abundance_species_ids

def create_simulation(cursor, simulation_id, simulation_description, number_species, species_ids, growth_rates_file, interaction_matrix_file, initial_abundances_file):
    """Create a new simulation experiment and save it to the database, including parameter file names."""
    
    # Check if species_ids exist in the database
    available_bacteria = get_available_bacteria(cursor, species_ids)
    available_species_ids = [str(row[0]) for row in available_bacteria]

    missing_species = set(species_ids) - set(available_species_ids)
    if missing_species:
        raise ValueError(f"The following species IDs are missing in the database: {missing_species}")

    # Insert the simulation data into the simulations table
    cursor.execute("""
    INSERT INTO simulations (simulation_id, simulation_description, number_species)
    VALUES (?, ?, ?);
    """, (simulation_id, simulation_description, number_species))

    # Link simulation with the species using the species_ids list
    for species_id in species_ids:
        cursor.execute("""
        INSERT INTO simulation_species (simulation_id, species_id)
        VALUES (?, ?);
        """, (simulation_id, species_id))

    # Insert the parameter file names into the simulation_parameters table
    cursor.execute("""
    INSERT INTO simulation_parameters (simulation_id, growth_rates_filename, interaction_matrix_filename, initial_abundances_filename)
    VALUES (?, ?, ?, ?);
    """, (simulation_id, growth_rates_file, interaction_matrix_file, initial_abundances_file))

    return simulation_id

def main():
    if len(sys.argv) != 2:
        print("Usage: python add_experiment.py <simulation_id>")
        sys.exit(1)

    # Get the simulation_id from the command-line argument
    simulation_id = sys.argv[1]

    # Connect to SQLite database
    conn = sqlite3.connect('simulations.db')
    cursor = conn.cursor()

    # Read simulation parameters from files
    try:
        # Define file paths
        growth_rates_file = f"{simulation_id}_growth_rates.csv"
        interaction_matrix_file = f"{simulation_id}_interaction_matrix.csv"
        initial_abundances_file = f"{simulation_id}_initial_abundances.csv"
        
        growth_rates_path = os.path.join(PARAMS_DIR, growth_rates_file)
        interaction_matrix_path = os.path.join(PARAMS_DIR, interaction_matrix_file)
        initial_abundances_path = os.path.join(PARAMS_DIR, initial_abundances_file)

        # Read files
        growth_rates = read_parameter_file(growth_rates_path)
        interaction_matrix = read_parameter_file(interaction_matrix_path)
        initial_abundances = read_parameter_file(initial_abundances_path)

        # Verify species IDs consistency between files
        extracted_species_ids = verify_species_ids(growth_rates, interaction_matrix, initial_abundances)
        number_species = len(extracted_species_ids)

        print(f"Verified species IDs (from parameter files): {extracted_species_ids}")

        # Verify species IDs exist in the database and then create a new simulation experiment
        create_simulation(cursor,
                          simulation_id=simulation_id,
                          simulation_description="Custom simulation using input files",
                          number_species=number_species,
                          species_ids=extracted_species_ids,
                          growth_rates_file=growth_rates_file,
                          interaction_matrix_file=interaction_matrix_file,
                          initial_abundances_file=initial_abundances_file)

        # Commit changes and close connection
        conn.commit()
        print(f"Simulation {simulation_id} successfully created!")
    
    except FileNotFoundError as e:
        print(f"Error: {e}. Ensure the parameter files are present in the {PARAMS_DIR} directory.")
    
    except ValueError as ve:
        print(f"Validation error: {ve}")
    
    finally:
        cursor.close()
        conn.close()

if __name__ == "__main__":
    main()