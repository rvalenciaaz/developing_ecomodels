python -m scripts.remove_experiment_database osci4
python -m scripts.add_experiment osci4
python -m scripts.simulation_sundials osci4 0 500 500
python -m scripts.plot_data osci4_sundials_0_500_500_simulation.csv