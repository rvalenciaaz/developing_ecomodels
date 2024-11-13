# developing_ecomodels
Simulate time series data using a generalised Lotka-Volterra (gLV) model or consumer-resource models in Python/sundials (C)

Commands:

1. `python -m scripts.create_database`

It runs a Python script (create_database) within the scripts module to create a new database. The database serves as the storage system for subsequent bacterial species, experiment definition (generalised Lotka-Volterra paramters), and simulation data.

2. `python -m scripts.add_bacteria_from_table generic_species_list.csv`
   
This command adds bacterial species to the database using data from a CSV file (eg. generic_species_list.csv), which is located in the species_list_tables folder. The script (add_bacteria_from_table) reads the file, extracts species information, and inserts the relevant data into the database.

3. `python -m scripts.remove_experiment_database osci4`

(Optional) This command removes an existing experiment (example id: osci4) from the database. It uses the script (remove_experiment_database) to locate and delete the data associated with the experiment in the database.

4. `python -m scripts.add_experiment osci4`

It adds a new experiment (id: osci4) to the database. The script (add_experiment) sets up the structure for storing the experimental data in the database, including any parameters or configurations needed for the experiment. Here, the parameters should be saved as csv files in the simulation_params folder and the first string before a `_` in the filename should match the expected ID string.

6. `python -m scripts.simulation_sundials osci4 0 500 500`

This command runs a generalised Lotka Volterra (gLV) simulation (using the C/sundials solver) for the experiment osci4. The arguments (0 500 500) represent the start time, end time, and the number of steps for the simulation, which will generate data over the specified range.

8. `python -m scripts.plot_data osci4_sundials_0_500_500_simulation.csv`
