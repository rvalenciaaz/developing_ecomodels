// glv_sundials_compat.c -------------------------------------------------------
// (patched) Avoids the “tout too close to t0” error by **not** asking CVODE to
// advance at step 0; instead we copy the initial condition directly into the
// results array and start integration at step 1.
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

int num_species = 0;
static realtype *r      = NULL;
static realtype **alpha = NULL;

static int gLV(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    (void)t; (void)user_data;
    for (int i = 0; i < num_species; ++i) {
        realtype inter = 0.0;
        for (int j = 0; j < num_species; ++j)
            inter += alpha[i][j] * NV_Ith_S(y, j);
        NV_Ith_S(ydot, i) = NV_Ith_S(y, i) * (r[i] + inter);
    }
    return 0;
}

void solve_gLV(double *results, int n_steps, double dt, double *init_abund)
{
#if SUNDIALS_VERSION_MAJOR >= 6
    if (!sunctx && SUNContext_Create(NULL, &sunctx)) {
        fprintf(stderr, "SUNContext_Create failed\n");
        return;
    }
#endif

    if (n_steps < 1) return; // nothing to do

    N_Vector y = NV_NEW_SERIAL(num_species);
    for (int i = 0; i < num_species; ++i) NV_Ith_S(y, i) = init_abund[i];

    // Store step 0 (initial condition) ------------------------------------
    for (int k = 0; k < num_species; ++k)
        results[0 * num_species + k] = NV_Ith_S(y, k);

    if (n_steps == 1) { N_VDestroy(y); FREE_SUNCTX(); return; }

    // Set up CVODE ---------------------------------------------------------
    void *cvode_mem = CVODE_CREATE(CV_ADAMS);
    if (!cvode_mem) { fprintf(stderr, "CVodeCreate error\n"); return; }

    realtype t0 = 0.0, t;
    if (CVodeInit(cvode_mem, gLV, t0, y) != CV_SUCCESS) {
        fprintf(stderr, "CVodeInit error\n"); return; }

    CVodeSStolerances(cvode_mem, 1e-4, 1e-8);

    SUNMatrix A = SUN_DENSE_MATRIX(num_species, num_species);
    SUNLinearSolver LS = SUN_DENSE_LSOL(y, A);
    CVODE_SET_LINSOL(cvode_mem, LS, A);

    // Advance starting from step 1 ---------------------------------------
    for (int step = 1; step < n_steps; ++step) {
        realtype tout = t0 + step * dt;   // tout > t0 for step ≥1
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

void set_parameters(int num, double *growth, double *inter_mat)
{
    num_species = num;
    r = (realtype*)malloc(num * sizeof(realtype));
    alpha = (realtype**)malloc(num * sizeof(realtype*));
    for (int i = 0; i < num; ++i) alpha[i] = (realtype*)malloc(num * sizeof(realtype));

    for (int i = 0; i < num; ++i) r[i] = growth[i];
    for (int i = 0; i < num; ++i)
        for (int j = 0; j < num; ++j)
            alpha[i][j] = inter_mat[i * num + j];
}

void free_parameters(void)
{
    free(r);
    for (int i = 0; i < num_species; ++i) free(alpha[i]);
    free(alpha);
#if SUNDIALS_VERSION_MAJOR >= 6
    FREE_SUNCTX();
#endif
}
