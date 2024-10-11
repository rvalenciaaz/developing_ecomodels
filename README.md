# developing_ecomodels
Simulate time series data using gLV and consumer-resource models in Python/sundials (C)

Commands:

1. `python -m scripts.create_database`
2. `python -m scripts.add_bacteria_from_table generic_species_list.csv`
This command adds bacterial species to the database using data from a CSV file (eg. generic_species_list.csv), which is located in the species_list_tables folder. The script (add_bacteria_from_table) reads the file, extracts species information, and inserts the relevant data into the database.
4. `python -m scripts.remove_experiment_database osci4`
5. `python -m scripts.add_experiment osci4`
6. `python -m scripts.simulation_sundials osci4 0 500 500`
7. `python -m scripts.plot_data osci4_sundials_0_500_500_simulation.csv`
