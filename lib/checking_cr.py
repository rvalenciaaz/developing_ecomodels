import numpy as np
import ctypes, pathlib, os

# ------------------------------------------------------------------ load lib
lib = np.ctypeslib.load_library("libcrsim", str(pathlib.Path(__file__).parent))

# ------------------------------------------------------------------ signatures
lib.set_parameters.argtypes = [
    ctypes.c_int, ctypes.c_int,                       # n_res, n_con
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1), # S
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1), # d
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1), # m
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1), # C  (row-major j×i)
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1)  # E  (row-major j×i)
]
lib.set_parameters.restype = None

lib.solve_consumer_resource.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1), # results
    ctypes.c_int, ctypes.c_double,                    # n_steps, dt
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1), # init_R
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1)  # init_N
]
lib.solve_consumer_resource.restype = None

lib.free_parameters.argtypes = []
lib.free_parameters.restype  = None

# ------------------------------------------------------------------ sample run
n_res, n_con   = 2, 3
n_steps, dt    = 500, 0.05             # 25 time-units total
num_state      = n_res + n_con

S = np.array([1.0, 0.5],           np.float64)        # supply
d = np.array([0.1, 0.1],           np.float64)        # decay
m = np.array([0.2, 0.3, 0.15],     np.float64)        # mortality

# consumption & efficiency (row-major: consumer j by resource i)
C = np.asarray([
    1.0, 0.2,
    0.3, 0.9,
    0.5, 0.4
], np.float64)

E = np.full_like(C, 0.6)                            # 60 % efficient everywhere

init_R = np.array([1.0, 0.8],        np.float64)
init_N = np.array([0.1, 0.05, 0.07], np.float64)

lib.set_parameters(n_res, n_con, S, d, m, C, E)

results = np.empty(n_steps * num_state, dtype=np.float64)
lib.solve_consumer_resource(results, n_steps, dt, init_R, init_N)
lib.free_parameters()

traj = results.reshape((n_steps, num_state))
t     = np.linspace(0, dt*(n_steps-1), n_steps)

print(traj[:5])          # first five rows
