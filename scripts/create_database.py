import sqlite3
import os

# Create folders for simulation parameters and results if they don't exist
def create_folders():
    os.makedirs('simulation_params', exist_ok=True)
    os.makedirs('simulation_results', exist_ok=True)

def create_bacterial_species_table(cursor):
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS bacterial_species (
        species_id TEXT PRIMARY KEY,
        genus TEXT,
        species TEXT,
        strain TEXT,
        taxonomy TEXT,
        genome_filename TEXT
    );
    """)

def create_simulations_table(cursor):
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS simulations (
        simulation_id TEXT PRIMARY KEY,
        simulation_description TEXT,
        number_species INTEGER
    );
    """)

def create_simulation_species_link_table(cursor):
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS simulation_species (
        link_id INTEGER PRIMARY KEY AUTOINCREMENT,
        simulation_id TEXT,
        species_id TEXT,
        FOREIGN KEY(simulation_id) REFERENCES simulations(simulation_id),
        FOREIGN KEY(species_id) REFERENCES bacterial_species(species_id)
    );
    """)

def create_simulation_parameters_table(cursor):
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS simulation_parameters (
        param_id INTEGER PRIMARY KEY AUTOINCREMENT,
        simulation_id TEXT,
        growth_rates_filename TEXT,
        interaction_matrix_filename TEXT,
        initial_abundances_filename TEXT,
        FOREIGN KEY(simulation_id) REFERENCES simulations(simulation_id)
    );
    """)

def create_simulation_results_table(cursor):
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS simulation_results (
        result_id INTEGER PRIMARY KEY AUTOINCREMENT,
        simulation_id TEXT,
        method_label TEXT,
        results_filename TEXT,
        plot_filename TEXT,
        dt REAL,  -- Time step (dt)
        n_steps INTEGER,  -- Number of steps (n_steps)
        FOREIGN KEY(simulation_id) REFERENCES simulations(simulation_id)
    );
    """)

def create_database():
    # Create necessary folders
    create_folders()
    # Connect to SQLite database (creates file if it doesn't exist)
    conn = sqlite3.connect('simulations.db')
    cursor = conn.cursor()

    # Create tables
    create_bacterial_species_table(cursor)
    create_simulations_table(cursor)
    create_simulation_species_link_table(cursor)
    create_simulation_parameters_table(cursor)
    create_simulation_results_table(cursor)

    # Commit changes and close connection
    conn.commit()
    cursor.close()
    conn.close()

if __name__ == "__main__":
    create_database()