import pandas as pd
import matplotlib.pyplot as plt
import os
import sys

def plot_comparison(simulation_filename, least_squares_filename):
    # Define the file paths
    simulation_file = os.path.join('simulation_results', simulation_filename)
    least_squares_file = os.path.join('least_squares', least_squares_filename)

    # Output folder and filename
    output_folder = 'least_squares'
    output_filename = simulation_filename.replace('simulation.csv', 'comparison.png')
    output_filepath = os.path.join(output_folder, output_filename)

    # Load the CSV files
    simulation_df = pd.read_csv(simulation_file)
    least_squares_df = pd.read_csv(least_squares_file)

    # Define the time column
    time_column = 'Time'

    # Create a scatter plot for the simulation data
    plt.figure(figsize=(10, 6))
    for column in simulation_df.columns[1:]:  # Exclude the 'Time' column
        plt.scatter(simulation_df[time_column], simulation_df[column], label=f'Simulation {column}', s=10)

    # Create a line plot for the least squares data
    for column in least_squares_df.columns[1:]:  # Exclude the 'Time' column
        plt.plot(least_squares_df[time_column], least_squares_df[column], label=f'Least Squares {column}', linestyle='-')

    # Add labels and title
    plt.xlabel('Time')
    plt.ylabel('Values')
    plt.title('Comparison of Simulation vs Least Squares Inferred Parameters')

    # Add legend
    plt.legend(loc='upper right', fontsize='small', ncol=2)

    # Ensure the output folder exists
    os.makedirs(output_folder, exist_ok=True)

    # Save the plot
    plt.savefig(output_filepath)

    # Show the plot
    plt.show()

    print(f"Plot saved as {output_filepath}")

if __name__ == "__main__":
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 3:
        print("Usage: python script.py <simulation_filename> <least_squares_filename>")
        sys.exit(1)

    # Get the filenames from the command line
    simulation_filename = sys.argv[1]
    least_squares_filename = sys.argv[2]

    # Run the plot comparison
    plot_comparison(simulation_filename, least_squares_filename)
