# developing_ecomodels
Simulate time series data using gLV and consumer-resource models in Python/sundials (C)

Commands:

1. `python -m scripts.create_database`

It runs a Python script (create_database) within the scripts module to create a new database. The database serves as the storage system for subsequent bacterial species, experiment definition (generalised Lotka-Volterra paramters), and simulation data.

2. `python -m scripts.add_bacteria_from_table generic_species_list.csv`
   
This command adds bacterial species to the database using data from a CSV file (eg. generic_species_list.csv), which is located in the species_list_tables folder. The script (add_bacteria_from_table) reads the file, extracts species information, and inserts the relevant data into the database.

3. `python -m scripts.remove_experiment_database osci4`

(Optional) This command removes an existing experiment (example id: osci4) from the database. It uses the script (remove_experiment_database) to locate and delete the data associated with the experiment in the database.

4. `python -m scripts.add_experiment osci4`


5. `python -m scripts.simulation_sundials osci4 0 500 500`

6. `python -m scripts.plot_data osci4_sundials_0_500_500_simulation.csv`
