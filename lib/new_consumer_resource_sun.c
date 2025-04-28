// consumer_resource_sundials.c --------------------------------------------------
// Consumer–resource ODE simulator (MacArthur form) that:
//   * Compiles on SUNDIALS ≤5 and ≥6.
//   * Avoids the "tout too close to t0" error by recording the initial
//     condition directly into the results array (step 0) and starting CVODE
//     integration at step 1.
// -----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_config.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <cvode/cvode_direct.h>
#include <sundials/sundials_types.h>

// -----------------------------------------------------------------------------
// Version‑compatibility shims --------------------------------------------------
#if SUNDIALS_VERSION_MAJOR >= 6
static SUNContext sunctx = NULL;
#define NV_NEW_SERIAL(len)              N_VNew_Serial((len), sunctx)
#define CVODE_CREATE(lmm)               CVodeCreate((lmm), sunctx)
#define SUN_DENSE_MATRIX(m,n)           SUNDenseMatrix((m), (n), sunctx)
#define SUN_DENSE_LSOL(y,A)             SUNLinSol_Dense((y), (A), sunctx)
#define CVODE_SET_LINSOL(cv,LS,A)       CVodeSetLinearSolver((cv), (LS), (A))
#define FREE_SUNCTX()                   do { if (sunctx) SUNContext_Free(&sunctx); } while (0)
#else
#define NV_NEW_SERIAL(len)              N_VNew_Serial((len))
#define CVODE_CREATE(lmm)               CVodeCreate((lmm), CV_NEWTON)
#define SUN_DENSE_MATRIX(m,n)           SUNDenseMatrix((m), (n))
#define SUN_DENSE_LSOL(y,A)             SUNDenseLinearSolver((y), (A))
#define CVODE_SET_LINSOL(cv,LS,A)       CVDlsSetLinearSolver((cv), (LS), (A))
#define FREE_SUNCTX()                   /* nothing */
#endif

// -----------------------------------------------------------------------------
// Model sizes & parameters -----------------------------------------------------
int num_resources = 0;   // # R
int num_consumers = 0;   // # N
int num_species   = 0;   // total state length

static realtype *S_rates = NULL;   // supply S_i
static realtype *d_rates = NULL;   // decay  d_i
static realtype *m_rates = NULL;   // mortality m_j
static realtype **C = NULL;        // consumption
static realtype **E = NULL;        // efficiency

#define R_IDX(i)  (i)
#define N_IDX(j)  (num_resources + (j))

// -----------------------------------------------------------------------------
// ODE right‑hand side ----------------------------------------------------------
static int cr_ode(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    (void)t; (void)user_data;

    // Resources
    for (int i = 0; i < num_resources; ++i) {
        realtype Ri = NV_Ith_S(y, R_IDX(i));
        realtype cons = 0.0;
        for (int j = 0; j < num_consumers; ++j)
            cons += C[j][i] * NV_Ith_S(y, N_IDX(j));
        NV_Ith_S(ydot, R_IDX(i)) = S_rates[i] - d_rates[i] * Ri - Ri * cons;
    }

    // Consumers
    for (int j = 0; j < num_consumers; ++j) {
        realtype Nj = NV_Ith_S(y, N_IDX(j));
        realtype intake = 0.0;
        for (int i = 0; i < num_resources; ++i)
            intake += E[j][i] * C[j][i] * NV_Ith_S(y, R_IDX(i));
        NV_Ith_S(ydot, N_IDX(j)) = Nj * (intake - m_rates[j]);
    }

    return 0;
}

// -----------------------------------------------------------------------------
// Solver wrapper --------------------------------------------------------------
void solve_consumer_resource(double *results, int n_steps, double dt,
                             double *init_R, double *init_N)
{
    num_species = num_resources + num_consumers;

#if SUNDIALS_VERSION_MAJOR >= 6
    if (!sunctx && SUNContext_Create(NULL, &sunctx)) {
        fprintf(stderr, "SUNContext_Create failed\n");
        return;
    }
#endif

    if (n_steps < 1) return;

    // Build initial vector ------------------------------------------------
    N_Vector y = NV_NEW_SERIAL(num_species);
    for (int i = 0; i < num_resources; ++i) NV_Ith_S(y, R_IDX(i)) = init_R[i];
    for (int j = 0; j < num_consumers; ++j) NV_Ith_S(y, N_IDX(j)) = init_N[j];

    // Step 0: save initial state -----------------------------------------
    for (int k = 0; k < num_species; ++k)
        results[0 * num_species + k] = NV_Ith_S(y, k);

    if (n_steps == 1) { N_VDestroy(y); FREE_SUNCTX(); return; }

    // CVODE setup ---------------------------------------------------------
    void *cvode_mem = CVODE_CREATE(CV_ADAMS);
    if (!cvode_mem) { fprintf(stderr, "CVodeCreate error\n"); return; }

    realtype t0 = 0.0, t;
    if (CVodeInit(cvode_mem, cr_ode, t0, y) != CV_SUCCESS) {
        fprintf(stderr, "CVodeInit error\n"); return; }

    CVodeSStolerances(cvode_mem, 1e-4, 1e-8);

    SUNMatrix A = SUN_DENSE_MATRIX(num_species, num_species);
    SUNLinearSolver LS = SUN_DENSE_LSOL(y, A);
    CVODE_SET_LINSOL(cvode_mem, LS, A);

    // Integration loop (start at step 1) ----------------------------------
    for (int step = 1; step < n_steps; ++step) {
        realtype tout = t0 + step * dt;  // strictly > t0
        if (CVode(cvode_mem, tout, y, &t, CV_NORMAL) < 0) {
            fprintf(stderr, "CVode error at step %d\n", step);
            break;
        }
        for (int k = 0; k < num_species; ++k)
            results[step * num_species + k] = NV_Ith_S(y, k);
    }

    // Cleanup -------------------------------------------------------------
    N_VDestroy(y);
    CVodeFree(&cvode_mem);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
    FREE_SUNCTX();
}

// -----------------------------------------------------------------------------
// Parameter I/O ---------------------------------------------------------------
void set_parameters(int n_res, int n_con,
                    double *S_in, double *d_in, double *m_in,
                    double *C_in, double *E_in)
{
    num_resources = n_res;
    num_consumers = n_con;
    num_species   = n_res + n_con;

    S_rates = (realtype*)malloc(n_res * sizeof(realtype));
    d_rates = (realtype*)malloc(n_res * sizeof(realtype));
    m_rates = (realtype*)malloc(n_con * sizeof(realtype));

    C = (realtype**)malloc(n_con * sizeof(realtype*));
    E = (realtype**)malloc(n_con * sizeof(realtype*));
    for (int j = 0; j < n_con; ++j) {
        C[j] = (realtype*)malloc(n_res * sizeof(realtype));
        E[j] = (realtype*)malloc(n_res * sizeof(realtype));
    }

    for (int i = 0; i < n_res; ++i) { S_rates[i] = S_in[i]; d_rates[i] = d_in[i]; }
    for (int j = 0; j < n_con; ++j) {
        m_rates[j] = m_in[j];
        for (int i = 0; i < n_res; ++i) {
            C[j][i] = C_in[j * n_res + i];
            E[j][i] = E_in[j * n_res + i];
        }
    }
}

void free_parameters(void)
{
    free(S_rates); free(d_rates); free(m_rates);
    for (int j = 0; j < num_consumers; ++j) { free(C[j]); free(E[j]); }
    free(C); free(E);
    S_rates = d_rates = m_rates = NULL; C = E = NULL;
#if SUNDIALS_VERSION_MAJOR >= 6
    FREE_SUNCTX();
#endif
}
