import pandas as pd
import matplotlib.pyplot as plt
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(message)s')

def plot_simulation_results(simulation_id):
    """
    Function to plot the results from a simulation CSV file.
    
    Args:
        simulation_filename (str): Path to the CSV file containing the simulation data.
    """
    # Load the data
    try:
        df = pd.read_csv(simulation_filename)
        logging.info(f"Loaded data from {simulation_filename}")
    except FileNotFoundError as e:
        logging.error(f"File not found: {simulation_filename}")
        return
    
    # Extract the species columns (everything except 'Time')
    time = df['Time']
    species_columns = [col for col in df.columns if col != 'Time']

    # Plot the species abundances
    plt.figure(figsize=(10, 6))
    
    for specie in species_columns:
        plt.plot(time, df[specie], label=specie)
    
    # Add title and labels
    plt.title('Species Abundances Over Time')
    plt.xlabel('Time')
    plt.ylabel('Abundance')
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.tight_layout()
    # Save the plot to a file
    plt.savefig(f'{simulation_id}_abundances.png', dpi=800, transparent=True,bbox_inches="tight")
    logging.info(f'Species abundances plot saved to {simulation_id}_abundances.png')

def main():
    """
    Main function to execute the plotting process.
    """
    simulation_filename = 'S92_sundials_simulation.csv'  # Adjust this if necessary
    plot_simulation_results(simulation_filename)

if __name__ == "__main__":
    main()