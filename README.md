# developing_ecomodels
Simulate time series data using gLV and consumer-resource models in Python/sundials (C)

Commands:

1. `python -m scripts.create_database`
2. `python -m scripts.add_bacteria_from_table generic_species_list.csv`
3. `python -m scripts.remove_experiment_database osci4`
4. `python -m scripts.add_experiment osci4`
5. `python -m scripts.simulation_sundials osci4 0 500 500`
6. `python -m scripts.plot_data osci4_sundials_0_500_500_simulation.csv`
