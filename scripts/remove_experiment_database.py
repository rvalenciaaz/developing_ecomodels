import sqlite3
import argparse

# Remove the simulation experiment and all its related data from the database
def remove_simulation(simulation_id):
    # Connect to the database
    conn = sqlite3.connect('simulations.db')
    cursor = conn.cursor()

    try:
        # Delete entries from the simulation_species table
        cursor.execute("DELETE FROM simulation_species WHERE simulation_id = ?", (simulation_id,))

        # Delete entries from the simulation_parameters table
        cursor.execute("DELETE FROM simulation_parameters WHERE simulation_id = ?", (simulation_id,))

        # Delete entries from the simulation_results table
        cursor.execute("DELETE FROM simulation_results WHERE simulation_id = ?", (simulation_id,))

        # Finally, delete the simulation entry itself
        cursor.execute("DELETE FROM simulations WHERE simulation_id = ?", (simulation_id,))

        # Commit the changes
        conn.commit()
        print(f"Simulation {simulation_id} and all associated database entries have been removed.")
    except Exception as e:
        conn.rollback()
        print(f"An error occurred: {e}")
    finally:
        # Close the connection
        cursor.close()
        conn.close()

if __name__ == "__main__":
    # Set up argument parser to accept simulation_id from the command line
    parser = argparse.ArgumentParser(description='Remove a simulation and related database entries.')
    parser.add_argument('simulation_id', type=str, help='The ID of the simulation to remove.')

    # Parse the arguments
    args = parser.parse_args()

    # Call the function to remove the simulation using the provided simulation_id
    remove_simulation(args.simulation_id)