import ctypes
import numpy as np
import logging
import os

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(message)s')

# Get the directory where glv_functions.py is located
lib_dir = os.path.dirname(os.path.abspath(__file__))

# Construct the full path to libglv.so
so_file = os.path.join(lib_dir, 'libglv.so')  # Change to 'glv.dll' on Windows

# Load the shared library using the absolute path
lib = ctypes.CDLL(so_file)

# Define argument and return types for the C functions
lib.set_parameters.argtypes = (
    ctypes.c_int,  # Number of species
    ctypes.POINTER(ctypes.c_double),  # Growth rates array
    ctypes.POINTER(ctypes.c_double)   # Interaction matrix array
)

lib.solve_gLV.argtypes = (
    ctypes.POINTER(ctypes.c_double),  # Results array (output)
    ctypes.c_int,  # Number of steps
    ctypes.c_double,  # Time step
    ctypes.POINTER(ctypes.c_double)   # Initial abundances array
)

lib.free_parameters.restype = None  # No return value for free_parameters

# Function to run the generalized Lotka-Volterra (gLV) simulation using the C library
def run_glv_simulation(growth_rates, interaction_matrix, initial_abundances, n_steps, dt):
    """
    Run the gLV simulation using a C library.

    Args:
        growth_rates (np.ndarray): Array of growth rates for each species.
        interaction_matrix (np.ndarray): Interaction matrix of species.
        initial_abundances (np.ndarray): Initial abundances of species.
        n_steps (int): Number of time steps.
        dt (float): Time step size.

    Returns:
        np.ndarray: Simulation results, a 2D array with shape (n_steps, num_species).
    """
    num_species = len(initial_abundances)  # Number of species
    print(num_species)

    # Convert input arrays to ctypes for the C library
    growth_rates_c = (ctypes.c_double * num_species)(*growth_rates)
    interaction_matrix_c = (ctypes.c_double * (num_species * num_species))(*interaction_matrix.flatten())
    initial_abundances_c = (ctypes.c_double * num_species)(*initial_abundances)
    results = (ctypes.c_double * (n_steps * num_species))()  # Array to store results

    # Set parameters in the C library
    logging.info("Setting parameters in the C library.")
    lib.set_parameters(num_species, growth_rates_c, interaction_matrix_c)

    # Solve the system using the C library
    logging.info("Solving the gLV system.")
    lib.solve_gLV(results, n_steps, dt, initial_abundances_c)

    # Free allocated memory in the C library
    logging.info("Freeing parameters in the C library.")
    lib.free_parameters()

    # Convert results back to a NumPy array
    results_np = np.array(results).reshape((n_steps, num_species))

    return results_np
