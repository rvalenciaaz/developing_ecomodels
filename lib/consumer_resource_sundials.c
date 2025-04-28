// consumer_resource_sim.c
// -----------------------------------------------------------------------------
// Minimal change version of the original gLV example that now simulates a
// MacArthur consumer‑resource model.  It keeps the same CVODE boiler‑plate and
// general structure, so you can compile/run it exactly the same way you did
// before.
//
//   State vector ordering  : [ R_0 .. R_{n_res-1}, N_0 .. N_{n_con-1} ]
//   Equations (for i = resources, j = consumers):
//       dR_i/dt =  S_i − d_i R_i − R_i * Σ_j C_{j,i} N_j
//       dN_j/dt =  N_j * ( Σ_i E_{j,i} C_{j,i} R_i − m_j )
//
// -----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <cvode/cvode.h>             // CVODE functions & constants
#include <nvector/nvector_serial.h>  // serial N_Vector interface
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <cvode/cvode_direct.h>
#include <sundials/sundials_types.h>

// -----------------------------------------------------------------------------
// Global sizes
// -----------------------------------------------------------------------------
int num_resources = 0;   // # resource species
int num_consumers = 0;   // # consumer species
int num_species   = 0;   // total state length ( = num_resources + num_consumers )

// -----------------------------------------------------------------------------
// Model parameters (set via set_parameters)
// -----------------------------------------------------------------------------
static realtype *S_rates = NULL;   // external supply rates  S_i      [n_res]
static realtype *d_rates = NULL;   // resource decay rates   d_i      [n_res]
static realtype *m_rates = NULL;   // consumer mortalities   m_j      [n_con]
static realtype **C      = NULL;   // consumption matrix   C_{j,i} [n_con][n_res]
static realtype **E      = NULL;   // conversion efficiency E_{j,i} [n_con][n_res]

// Helpers to access entries in y / ydot
#define R_IDX(i)          (i)                       // resource i
#define N_IDX(j)          (num_resources + (j))     // consumer j

// -----------------------------------------------------------------------------
// Derivative function: consumer‑resource ODE system
// -----------------------------------------------------------------------------
static int cr_ode(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    (void)t;           // unused
    (void)user_data;   // signature compliance

    // 1) Resources ---------------------------------------------------------
    for (int i = 0; i < num_resources; ++i) {
        realtype Ri = NV_Ith_S(y, R_IDX(i));
        realtype cons = 0.0;
        for (int j = 0; j < num_consumers; ++j) {
            realtype Nj = NV_Ith_S(y, N_IDX(j));
            cons += C[j][i] * Nj;
        }
        NV_Ith_S(ydot, R_IDX(i)) = S_rates[i] - d_rates[i] * Ri - Ri * cons;
    }

    // 2) Consumers ---------------------------------------------------------
    for (int j = 0; j < num_consumers; ++j) {
        realtype Nj = NV_Ith_S(y, N_IDX(j));
        realtype intake = 0.0;
        for (int i = 0; i < num_resources; ++i) {
            realtype Ri = NV_Ith_S(y, R_IDX(i));
            intake += E[j][i] * C[j][i] * Ri;
        }
        NV_Ith_S(ydot, N_IDX(j)) = Nj * (intake - m_rates[j]);
    }

    return 0;  // success
}

// -----------------------------------------------------------------------------
// Solver wrapper – almost identical to original solve_gLV() --------------------
// -----------------------------------------------------------------------------
void solve_consumer_resource(double *results, int n_steps, double dt,
                             double *initial_resources,
                             double *initial_consumers)
{
    num_species = num_resources + num_consumers;

    realtype t0 = 0.0, t;

    // Build initial state vector
    N_Vector y = N_VNew_Serial(num_species);
    for (int i = 0; i < num_resources; ++i)
        NV_Ith_S(y, R_IDX(i)) = initial_resources[i];
    for (int j = 0; j < num_consumers; ++j)
        NV_Ith_S(y, N_IDX(j)) = initial_consumers[j];

    // Create CVODE memory block (unchanged API usage)
    void *cvode_mem = CVodeCreate(CV_ADAMS, CV_NEWTON);
    if (cvode_mem == NULL) { fprintf(stderr, "CVodeCreate error\n"); return; }

    if (CVodeInit(cvode_mem, cr_ode, t0, y) != CV_SUCCESS) {
        fprintf(stderr, "CVodeInit error\n"); return; }

    CVodeSStolerances(cvode_mem, 1e-4, 1e-8);

    // Dense linear solver, exactly as before
    SUNMatrix A = SUNDenseMatrix(num_species, num_species);
    SUNLinearSolver LS = SUNDenseLinearSolver(y, A);
    CVodeSetLinearSolver(cvode_mem, LS, A);

    // Integration loop
    t = t0;
    for (int step = 0; step < n_steps; ++step) {
        realtype tout = t0 + step * dt;
        if (CVode(cvode_mem, tout, y, &t, CV_NORMAL) < 0) {
            fprintf(stderr, "CVode failure at step %d\n", step);
            break;
        }
        // copy to results (row‑major)
        for (int k = 0; k < num_species; ++k)
            results[step * num_species + k] = NV_Ith_S(y, k);
    }

    // Housekeeping
    N_VDestroy(y);
    CVodeFree(&cvode_mem);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
}

// -----------------------------------------------------------------------------
// Parameter setup – simple C arrays, no SUNDIALS types -------------------------
// -----------------------------------------------------------------------------
void set_parameters(int n_res, int n_con,
                    double *S_in, double *d_in, double *m_in,
                    double *C_in, double *E_in)
{
    num_resources = n_res;
    num_consumers = n_con;
    num_species   = n_res + n_con;

    // 1‑D arrays -----------------------------------------------------------
    S_rates = (realtype*)malloc(n_res * sizeof(realtype));
    d_rates = (realtype*)malloc(n_res * sizeof(realtype));
    m_rates = (realtype*)malloc(n_con * sizeof(realtype));

    // 2‑D matrices ---------------------------------------------------------
    C = (realtype**)malloc(n_con * sizeof(realtype*));
    E = (realtype**)malloc(n_con * sizeof(realtype*));
    for (int j = 0; j < n_con; ++j) {
        C[j] = (realtype*)malloc(n_res * sizeof(realtype));
        E[j] = (realtype*)malloc(n_res * sizeof(realtype));
    }

    // Copy scalars ---------------------------------------------------------
    for (int i = 0; i < n_res; ++i) {
        S_rates[i] = S_in[i];
        d_rates[i] = d_in[i];
    }
    for (int j = 0; j < n_con; ++j) {
        m_rates[j] = m_in[j];
        for (int i = 0; i < n_res; ++i) {
            C[j][i] = C_in[j * n_res + i];
            E[j][i] = E_in[j * n_res + i];
        }
    }
}

// -----------------------------------------------------------------------------
// Clean‑up ---------------------------------------------------------------------
// -----------------------------------------------------------------------------
void free_parameters(void)
{
    free(S_rates);  free(d_rates);  free(m_rates);
    for (int j = 0; j < num_consumers; ++j) {
        free(C[j]); free(E[j]);
    }
    free(C);  free(E);

    S_rates = d_rates = m_rates = NULL;
    C = E = NULL;
}
