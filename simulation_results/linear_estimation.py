import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

# Load the data
file_path = 'osci4_sundials_0_250_250_simulation.csv'
data = pd.read_csv(file_path)

# Calculate the time intervals (assumes evenly spaced data)
time_intervals = data['Time'].diff().dropna().values
mean_time_interval = np.mean(time_intervals)  # Assume a constant time step for simplicity

# Compute the derivative of the logarithm of the abundances
log_abundances = np.log(data[['sp1', 'sp2', 'sp3', 'sp4']])
log_abundance_derivative = log_abundances.diff() / mean_time_interval

# Drop the first row due to differencing
log_abundance_derivative = log_abundance_derivative.iloc[1:].reset_index(drop=True)
abundances = data[['sp1', 'sp2', 'sp3', 'sp4']].iloc[1:].reset_index(drop=True)

# Concatenate to form the transformed dataset for linear regression
transformed_data = pd.concat([log_abundance_derivative, abundances], axis=1)
transformed_data.columns = [
    'd_log_sp1', 'd_log_sp2', 'd_log_sp3', 'd_log_sp4', 'sp1', 'sp2', 'sp3', 'sp4'
]

# Perform linear regression for each species derivative
coefficients = {}
intercepts = {}

for species in ['sp1', 'sp2', 'sp3', 'sp4']:
    # Response variable and predictors
    response = transformed_data[f'd_log_{species}']
    predictors = transformed_data[['sp1', 'sp2', 'sp3', 'sp4']]
    
    # Fit the linear regression model
    model = LinearRegression()
    model.fit(predictors, response)
    
    # Store the coefficients and intercept for each species
    coefficients[species] = model.coef_
    intercepts[species] = model.intercept_

# Organize the coefficients and intercepts into a DataFrame for easy interpretation
regression_results = pd.DataFrame(coefficients, index=['coef_sp1', 'coef_sp2', 'coef_sp3', 'coef_sp4']).T
regression_results['intercept'] = pd.Series(intercepts)

# Initialize dictionary to store predicted abundances with initial values adjusted
predicted_abundances = {species: [abundances.iloc[0][species]] for species in ['sp1', 'sp2', 'sp3', 'sp4']}

# Predict log-derivatives for each time step and integrate to get abundances
for t in range(1, len(abundances)):
    # Current abundances as predictors
    current_abundances = {species: predicted_abundances[species][t - 1] for species in ['sp1', 'sp2', 'sp3', 'sp4']}
    
    for species in ['sp1', 'sp2', 'sp3', 'sp4']:
        # Calculate predicted log-derivative using regression coefficients
        d_log_abundance = (
            sum(regression_results.loc[species, f'coef_{sp}'] * current_abundances[sp] for sp in ['sp1', 'sp2', 'sp3', 'sp4'])
            + regression_results.loc[species, 'intercept']
        )
        
        # Update abundance based on the previous abundance and predicted log-derivative
        previous_abundance = predicted_abundances[species][t - 1]
        predicted_abundance = previous_abundance * np.exp(d_log_abundance * mean_time_interval)
        
        # Append the new predicted abundance
        predicted_abundances[species].append(predicted_abundance)

# Convert the predictions to a DataFrame for plotting
predicted_abundances_df = pd.DataFrame(predicted_abundances)
predicted_abundances_df['Time'] = data['Time'][1:].values

# Adjust the time array for plotting to match the derivatives (starting from second time point)
time_adjusted = data['Time'][1:].values  # Ignore the first time point

# Plot actual vs predicted abundances for each species
for species in ['sp1', 'sp2', 'sp3', 'sp4']:
    plt.figure()
    plt.plot(time_adjusted, abundances[species], label='Actual', marker='o', linestyle='-')
    plt.plot(predicted_abundances_df['Time'], predicted_abundances_df[species], label='Predicted', linestyle='--')
    plt.title(f'Actual vs Predicted Abundance for {species} (Adjusted Initial)')
    plt.xlabel('Time')
    plt.ylabel('Abundance')
    plt.legend()
    plt.show()
