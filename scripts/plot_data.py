import pandas as pd
import matplotlib.pyplot as plt
import logging
import argparse
import os

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(message)s')

def plot_simulation_results(filename):
    """
    Function to plot the results from a simulation CSV file.

    Args:
        filename (str): Filename of the CSV containing simulation results.
    """
    # Define the directory containing the simulation results
    input_dir = 'simulation_results'
    
    # Construct the full path to the file
    filepath = os.path.join(input_dir, filename)
    
    # Load the data
    try:
        df = pd.read_csv(filepath)
        logging.info(f"Loaded data from {filepath}")
    except FileNotFoundError as e:
        logging.error(f"File not found: {filepath}")
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

    # Extract the base filename (without the directory) to generate the output filename
    base_filename = os.path.basename(filename)
    output_filename = os.path.splitext(base_filename)[0]

    # Ensure the directory 'simulation_plots' exists
    output_dir = 'simulation_plots'
    os.makedirs(output_dir, exist_ok=True)
    
    # Save the plot to the 'simulation_plots' directory
    plot_filename = os.path.join(output_dir, f'{output_filename}_abundances.png')
    plt.savefig(plot_filename, dpi=800, transparent=True, bbox_inches="tight")
    logging.info(f'Species abundances plot saved to {plot_filename}')

def main():
    """
    Main function to execute the plotting process.
    """
    # Use argparse to get command-line arguments
    parser = argparse.ArgumentParser(description='Plot simulation results from a CSV file.')
    
    # Define expected argument: the CSV filename
    parser.add_argument('filename', type=str, help='The CSV file (located in the simulation_results folder) containing the simulation results')
    
    # Parse the arguments
    args = parser.parse_args()

    # Call the plotting function with the parsed filename
    plot_simulation_results(args.filename)

if __name__ == "__main__":
    main()