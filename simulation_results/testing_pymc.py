import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import pymc as pm
from pymc.ode import DifferentialEquation
import arviz as az
import matplotlib.pyplot as plt

# Load the CSV data into a pandas dataframe
data = pd.read_csv("osci4_sundials_0_250_250_simulation.csv")
time = data['Time'].values
species_data = data.drop(columns=['Time']).values  # Exclude the 'Time' column

# Define the GLV model ODE function compatible with PyMC's DifferentialEquation
def glv_ode(y, t, theta):
    α = theta[:4]
    β = theta[4:].reshape((4,4))
    interaction = pm.math.dot(β, y)
    du_dt = y * (α + interaction)
    return du_dt

# Set up the DifferentialEquation object
glv_ode_model = DifferentialEquation(
    func=glv_ode,
    times=time,
    n_states=4,
    n_theta=20,  # 4 α and 16 β
    t0=time[0]
)

# Bayesian Model Definition with PyMC
with pm.Model() as model:
    # Priors for growth rates α
    α = pm.LogNormal('α', mu=0, sigma=0.25, shape=4)
    
    # Priors for interaction coefficients β
    β = pm.Normal('β', mu=-1, sigma=2, shape=16)
    
    # Prior for noise term σ
    σ = pm.InverseGamma('σ', alpha=2, beta=3)
    
    # Combine α and β into a single parameter vector
    theta = pm.math.concatenate([α, β])
    
    # Solve the ODE
    glv_solution = glv_ode_model(y0=u0, theta=theta)
    
    # Likelihood: Gaussian with observed data and noise σ
    Y_obs = pm.Normal('Y_obs', mu=glv_solution, sigma=σ, observed=species_data)
    
    # Sample from the posterior
    trace = pm.sample(1000, tune=1000, chains=4, cores=4, init='adapt_diag', target_accept=0.9)
    
# Summary and Diagnostics
print(az.summary(trace))

# Plot Posterior Distributions
az.plot_trace(trace)
plt.show()

# Extract posterior means
α_est = trace.posterior['α'].mean(dim=('chain', 'draw')).values
β_est = trace.posterior['β'].mean(dim=('chain', 'draw')).values

# Simulate the model with posterior mean parameters
theta_est = np.concatenate([α_est, β_est])

# Define the ODE function for scipy
def glv_ode_for_scipy(t, y, p):
    α = p[:4]
    β = p[4:].reshape((4,4))
    interaction = β @ y
    du_dt = y * (α + interaction)
    return du_dt

# Solve the ODE with estimated parameters
sol_estimated = solve_ivp(
    fun=lambda t, y: glv_ode_for_scipy(t, y, theta_est),
    t_span=(time[0], time[-1]),
    y0=u0,
    t_eval=time
)

sol_estimated_mat = sol_estimated.y.T  # Transpose to match data dimensions

# Plot results
plt.figure(figsize=(10,6))
for i in range(4):
    plt.plot(time, species_data[:,i], label=f"Species {i+1} Data")
    plt.plot(time, sol_estimated_mat[:,i], linestyle='--', label=f"Species {i+1} Estimate")
plt.legend()
plt.xlabel('Time')
plt.ylabel('Species Concentration')
plt.title('GLV Model Fit to Data')
plt.show()
