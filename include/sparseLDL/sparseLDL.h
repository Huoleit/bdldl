#pragma once

#include <stddef.h>

#include "sparseLDL/sparseLDL_types.h"

#ifdef __cplusplus
extern "C" {
#endif

struct qp_size {
  int *nx;  // number of states
  int *nu;  // number of inputs
  int *ng;  // number of constraints
  int N;    // horizon length
  size_t memsize;
};

struct qp_data {
  // Dynamics
  ldl_matrix *A;
  ldl_matrix *B;
  ldl_matrix *C;
  ldl_matrix *D;

  // Cost
  ldl_matrix *Q;
  ldl_matrix *R;
  ldl_matrix *H;
  ldl_vector *q;
  ldl_vector *r;

  // Constraints
  ldl_matrix *E;
  ldl_matrix *F;
  ldl_vector *g;

  size_t memsize;
};

struct ldl_workspace {
  ldl_matrix *L;
  ldl_matrix *D;
  ldl_matrix *DInv;

  size_t memsize;
};

/////////////////////////////////////////////
// Problem size
/////////////////////////////////////////////

/**
 * Get the size of the memory required to store the problem size
 *
 * @param N horizon length
 * @return size_t size of memory required
 */
size_t qp_size_getRequiredMemorySize(const int N);

/**
 * Map the problem size struct to pre-allocated memory
 *
 * @param N horizon length
 * @param dim problem size struct storing the size information of every time step
 * @param memory pointer to pre-allocated memory
 */
void qp_size_mapToMemory(const int N, struct qp_size *dim, void *memory);

/**
 * Initialize the problem size struct
 *
 * @param nx number of states
 * @param nu number of inputs
 * @param ng number of constraints
 * @param dim problem size struct storing the size information of every time step
 */
void qp_size_setAll(const int *nx, const int *nu, const int *ng, struct qp_size *dim);

/////////////////////////////////////////////
// Problem data
/////////////////////////////////////////////

/**
 * Get the size of the memory required to store the problem data
 *
 * @param dim problem size struct storing the size information of every time step
 * @return size_t
 */
size_t qp_data_getRequiredMemorySize(const struct qp_size *dim);

/**
 * Map the problem data struct to pre-allocated memory
 *
 * @param data problem data struct storing the data of every time step
 * @param memory pointer to pre-allocated memory
 */
void qp_data_mapToMemory(struct qp_data *data, void *memory);

/**
 * Initialize the problem data struct
 *
 * @param A dynamics matrix A
 * @param B dynamics matrix B
 * @param C dynamics matrix C
 * @param D dynamics matrix D
 * @param Q cost matrix Q
 * @param R cost matrix R
 * @param H cost matrix H
 * @param q cost vector q
 * @param r cost vector r
 * @param E constraints matrix E
 * @param F constraints matrix F
 * @param g constraints vector g
 * @param dim problem size struct storing the size information of every time step
 * @param data problem data struct storing the data of every time step
 */
void qp_data_setAll(const ldl_float *A, const ldl_float *B, const ldl_float *C,
                    const ldl_float *D, const ldl_float *Q, const ldl_float *R,
                    const ldl_float *H, const ldl_float *q, const ldl_float *r,
                    const ldl_float *E, const ldl_float *F, const ldl_float *g,
                    const struct qp_size *dim, struct qp_data *data);

/////////////////////////////////////////////
// factorization
/////////////////////////////////////////////

/**
 * Get the size of the memory required to store the factorization workspace
 *
 * @param dim
 * @return size_t
 */
size_t sparseLDL_getRequiredMemorySize(const struct qp_size *dim);

/**
 * Map the factorization workspace to pre-allocated memory
 *
 * @param workspace
 * @param memory
 */
void sparseLDL_mapToMemory(struct ldl_workspace *ws, void *memory);

/**
 * Factorize the KKT matrix
 *
 * @param dim problem size struct storing the size information of every time step
 * @param data problem data struct storing the data of every time step
 * @param workspace factorization workspace
 */
void sparseLDL_factor(const struct qp_data *data, struct ldl_workspace *ws);

/**
 * Solve the KKT system in place
 *
 * @param workspace factorization workspace containing the factorization
 * @param b solution vector
 */
void sparseLDL_solveInPlace(const struct ldl_workspace *ws, ldl_vector *b);

/**
 * Solve the KKT system (non destructive)
 *
 * @param ws factorization workspace containing the factorization
 * @param b RHS vector
 * @param result solution vector
 */
void sparseLDL_solve(const struct ldl_workspace *ws, const ldl_vector *b,
                     ldl_vector *result);

#ifdef __cplusplus
}  // #extern "C"
#endif