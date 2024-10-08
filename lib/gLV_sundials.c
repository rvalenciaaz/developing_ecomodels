#include <stdio.h>
#include <stdlib.h>
#include <cvode/cvode.h>             // prototypes for CVODE functions and constants
#include <nvector/nvector_serial.h>  // serial N_Vector types, functions, and macros
#include <sunmatrix/sunmatrix_dense.h> // access to dense SUNMatrix
#include <sunlinsol/sunlinsol_dense.h> // access to dense SUNLinearSolver
#include <cvode/cvode_direct.h> // access to CVDls interface
#include <sundials/sundials_types.h>  // defs. of realtype, sunindextype

// Global parameters for the model
int num_species;
realtype *r;  // growth rates
realtype **alpha;  // interaction coefficients

// Function to compute the derivatives (gLV system)
int gLV(realtype t, N_Vector y, N_Vector ydot, void *user_data) {
    for (int i = 0; i < num_species; i++) {
        realtype interaction = 0.0;
        for (int j = 0; j < num_species; j++) {
            interaction += alpha[i][j] * NV_Ith_S(y, j);
        }
        NV_Ith_S(ydot, i) = NV_Ith_S(y, i) * (r[i] + interaction);
    }
    return 0;
}

// Function to solve the gLV system and store results in an array
void solve_gLV(double *results, int n_steps, double dt, double *initial_abundances) {
    realtype t0 = 0.0;
    realtype t;
    realtype T = dt * (n_steps - 1);

    // Create a serial vector for the initial abundances
    N_Vector y = N_VNew_Serial(num_species);
    for (int i = 0; i < num_species; i++) {
        NV_Ith_S(y, i) = initial_abundances[i];
    }

    // Create the CVODE memory block
    void *cvode_mem = CVodeCreate(CV_ADAMS, CV_NEWTON);
    if (cvode_mem == NULL) {
        fprintf(stderr, "Error in CVodeCreate\n");
        return;
    }

    // Initialize CVODE
    int flag = CVodeInit(cvode_mem, gLV, t0, y);
    if (flag != CV_SUCCESS) {
        fprintf(stderr, "Error in CVodeInit\n");
        return;
    }

    // Specify the relative and absolute tolerances
    flag = CVodeSStolerances(cvode_mem, 1e-4, 1e-8);
    if (flag != CV_SUCCESS) {
        fprintf(stderr, "Error in CVodeSStolerances\n");
        return;
    }

    // Create the dense SUNMatrix
    SUNMatrix A = SUNDenseMatrix(num_species, num_species);
    if (A == NULL) {
        fprintf(stderr, "Error in SUNDenseMatrix\n");
        return;
    }

    // Create the dense SUNLinearSolver
    SUNLinearSolver LS = SUNDenseLinearSolver(y, A);
    if (LS == NULL) {
        fprintf(stderr, "Error in SUNLinSol_Dense\n");
        return;
    }

    // Attach the linear solver to CVODE
    flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
    if (flag != CV_SUCCESS) {
        fprintf(stderr, "Error in CVodeSetLinearSolver\n");
        return;
    }

    // Time-stepping loop
    t = t0;
    for (int i = 0; i < n_steps; i++) {
        flag = CVode(cvode_mem, t + dt, y, &t, CV_NORMAL);
        if (flag != CV_SUCCESS) {
            fprintf(stderr, "Error in CVode\n");
            return;
        }
        for (int j = 0; j < num_species; j++) {
            results[i * num_species + j] = NV_Ith_S(y, j);
        }
    }

    // Free memory
    N_VDestroy(y);
    CVodeFree(&cvode_mem);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
}

// Function to set model parameters from Python
void set_parameters(int num, double *growth_rates, double *interaction_matrix) {
    num_species = num;
    
    // Allocate memory for growth rates and interaction matrix
    r = (realtype *)malloc(num_species * sizeof(realtype));
    alpha = (realtype **)malloc(num_species * sizeof(realtype *));
    for (int i = 0; i < num_species; i++) {
        alpha[i] = (realtype *)malloc(num_species * sizeof(realtype));
    }

    // Set growth rates
    for (int i = 0; i < num_species; i++) {
        r[i] = growth_rates[i];
    }

    // Set interaction coefficients
    for (int i = 0; i < num_species; i++) {
        for (int j = 0; j < num_species; j++) {
            alpha[i][j] = interaction_matrix[i * num_species + j];
        }
    }
}

// Function to free the memory for parameters
void free_parameters() {
    free(r);
    for (int i = 0; i < num_species; i++) {
        free(alpha[i]);
    }
    free(alpha);
}
