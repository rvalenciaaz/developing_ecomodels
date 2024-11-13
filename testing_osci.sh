#python -m scripts.create_database
#python -m scripts.add_bacteria_from_table generic_species_list.csv
python -m scripts.remove_experiment_database osci4
python -m scripts.add_experiment osci4
python -m scripts.simulation_sundials osci4 0 500 500
python -m scripts.plot_data osci4_sundials_0_500_500_simulation.csv
#python -m scripts.glv_least_squares osci4_sundials_0_500_500_simulation.csv
#python -m scripts.plot_comparison osci4_sundials_0_500_500_simulation.csv osci4LS1_sundials_0_500_500_simulation_least_squares.csv
