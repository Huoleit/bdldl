/**
 * @file bdldl.h
 * @author Fu Zhengyu (zhengfuaj@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-08-14
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

#include <stddef.h>

#include "bdldl/bdldl_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/////////////////////////////////////////////
// Forward declaration
/////////////////////////////////////////////

struct ldl_Solver;

/////////////////////////////////////////////
// Solver status
/////////////////////////////////////////////

typedef enum {
  bdldl_OK = 0,                      //!<  OK: No Error
  bdldl_ERROR_SINGULAR,              //!<  ERROR: KKT matrix is singular
  bdldl_ERROR_OUT_OF_MEMORY,         //!<  ERROR: Memory allocated is less than required
  bdldl_ERROR_OUT_OF_BOUNDARY,       //!<  ERROR: Index out of boundary
  bdldl_ERROR_DIMENSION_MISMATCH,    //!<  ERROR: Dimension of augments mismatch
  bdldl_ERROR_FREE_EXTERNAL_MEMORY,  //!<  ERROR: Cannot free memory allocated externally

} bdldl_status;

/////////////////////////////////////////////
// Factorization & Solve
/////////////////////////////////////////////

/**
 * @brief Factorize the KKT matrix
 *
 * @param[in, out] solver Pointer to the solver
 * @param[in] primal_regularization Regularization for the primal variables
 * @param[in] dual_regularization Regularization for the dual variables
 * @return Status of the factorization. See \ref bdldl_status for details.
 */

bdldl_status bdldl_factor(struct ldl_Solver *solver, const ldl_float primal_regularization,
                          const ldl_float dual_regularization);

/**
 * @brief Solve the KKT system (destructive)
 *
 * @param[in] solver Pointer to the solver
 * @param[in, out] b Vector to solve for. The solution is stored in this vector.
 * @return Status of the solve. See \ref bdldl_status for details.
 */
bdldl_status bdldl_solveInPlace(const struct ldl_Solver *solver, ldl_vector *b);

/**
 * Solve the KKT system (non destructive)
 *
 * @param[in] solver Pointer to the solver
 * @param[in] b Vector to solve fo
 * @param[out] result solution vector
 * @return Status of the solve. See \ref bdldl_status for details.
 */
bdldl_status bdldl_solve(const struct ldl_Solver *solver, const ldl_vector *b,
                         ldl_vector *result);

/////////////////////////////////////////////
// bdldl solver
/////////////////////////////////////////////

/**
 * Get required memory size for the solver
 *
 * @param[in] num_states State dimension of every stage. Size `N + 1`
 * @param[in] num_inputs State dimension of every stage. Size `N + 1`
 * @param[in] num_constraints Constraint dimension of every stage. Size `N + 1`
 * @param[in] num_horizon Number of intermediate states `N`. For example, if N = 2, then
 * x0(initial state) and x1 are intermediate states, x2 is the final state.
 * @param[in] use_diagonal_costs Use diagonal costs. If true, only memory for diagonal
 * entries is allocated; otherwise, memory for dense Q and R matrices is allocated.
 * @param[out] status Status of the calculation. See \ref bdldl_status for details.
 * @return size_t Required memory size
 */
size_t bdldl_getRequiredMemorySize(const int *num_states, const int *num_inputs,
                                   const int *num_constraints, int num_horizon,
                                   ldl_bool use_diagonal_costs, bdldl_status *status);

/**
 * (Allocate required memory &) map the solver to the given memory(or the internally
 * allocated memory).
 *
 * @param[in] num_states State dimension of every stage. Size `N + 1`
 * @param[in] num_inputs State dimension of every stage. Size `N + 1`
 * @param[in] num_constraints Constraint dimension of every stage. Size `N + 1`
 * @param[in] num_horizon Number of intermediate states `N`. See \ref
 * bdldl_getRequiredMemorySize for details.
 * @param[in] use_diagonal_costs Use diagonal costs. If true, only memory for diagonal
 * entries is allocated; otherwise, memory for dense Q and R matrices is allocated.
 * @param[in] memory If NULL, required memory is allocated internally. Users should call
 * \ref bdldl_freeSolver to free the allocated memory explicitly. If a valid memory
 * address, solver is mapped to this pre-allocated memory space. Users should manage the
 * memory on their own.
 * @param[out] status Status of the allocation. See \ref bdldl_status for details.
 * @return struct ldl_Solver* Pointer to the solver. Return NULL if errors occur.
 */
struct ldl_Solver *bdldl_newSolver(const int *num_states, const int *num_inputs,
                                   const int *num_constraints, int num_horizon,
                                   ldl_bool use_diagonal_costs, void *memory,
                                   bdldl_status *status);

/**
 * Release the memory allocated by the solver and NULL all internal pointers.
 *
 * @param[in] solver Pointer to the solver
 * @return Status of the free operation. See \ref bdldl_status for details.
 */
bdldl_status bdldl_freeSolver(struct ldl_Solver *solver);

/////////////////////////////////////////////
// Costs setters
/////////////////////////////////////////////

/**
 * @brief Set the quadratic and affine cost terms on the state for a single time step.
 *
 * @param[in] solver Pointer to the solver
 * @param[in] Q Quadratic state cost. Size is determined by whether the costs are dense or
 * diagonal. If NULL, no quadratic cost is set and the corresponding memory remains
 * unchanged.
 * @param[in] q Linear state cost. Size `num_states[k]`. If NULL, no
 * linear state cost is set and the corresponding memory remains unchanged.
 * @param[in] k Time step at which to set the cost.
 * @return Status of the set operation. See \ref bdldl_status for details.
 */
bdldl_status bdldl_setStateCosts(struct ldl_Solver *solver, const ldl_float *Q,
                                 const ldl_float *q, int k);

/**
 * @brief Set the quadratic and affine cost terms on the input for a single time step.
 *
 * @param[in] solver Pointer to the solver
 * @param[in] R Quadratic input cost. Size is determined by whether the costs are dense or
 * diagonal. If NULL, no quadratic cost is set and the corresponding memory remains
 * unchanged.
 * @param[in] r Linear input cost. Size `num_inputs[k]`. If NULL, no linear input cost is
 * set and the corresponding memory remains unchanged.
 * @param[in] k Time step at which to set the cost.
 * @return Status of the set operation. See \ref bdldl_status for details.
 */
bdldl_status bdldl_setInputCosts(struct ldl_Solver *solver, const ldl_float *R,
                                 const ldl_float *r, int k);

/**
 * @brief Set the cross-term cost between the states and inputs.
 *
 * @param[in] solver Pointer to the solver
 * @param[in] H Cross-term cost. Size `(num_inputs[k], num_states[k])`.
 * @param[in] k Time step at which to set the cost.
 * @return Status of the set operation. See \ref bdldl_status for details.
 */
bdldl_status bdldl_setStateInputCosts(struct ldl_Solver *solver, const ldl_float *Hux,
                                      int k);

/////////////////////////////////////////////
// Dynamics setters
/////////////////////////////////////////////

/**
 * @brief Set the dynamics at a given time step
 *
 * @param[in]  solver Pointer to the solver
 * @param[in]  A State transition matrix. Size `(num_states[k+1], num_states[k])`. If NULL,
 * A is not set and the corresponding memory remains unchanged.
 * @param[in]  B Input matrix. Size `(num_states[k+1], num_inputs[k])`. If NULL, B is not
 * set and the corresponding memory remains unchanged.
 * @param[in]  C Next state matrix. Size (`num_states[k+1], num_states[k+1])`. If NULL, C is
 * not set and the corresponding memory remains unchanged.
 * @param[in]  D Next input matrix. Size (`num_states[k+1], num_inputs[k+1])`. If NULL, D is
 * not set and the corresponding memory remains unchanged.
 * @param[in]  k Time step at which to set the dynamics.
 * @return Status of the set operation. See \ref bdldl_status for details.
 */
bdldl_status bdldl_setDynamics(struct ldl_Solver *solver, const ldl_float *A,
                               const ldl_float *B, const ldl_float *C, const ldl_float *D,
                               int k);

/////////////////////////////////////////////
// Constraints setters
/////////////////////////////////////////////

/**
 * @brief
 *
 * @param solver
 * @param E
 * @param F
 * @param g
 * @param k Time step at which to set the constraint.
 * @return Status of the set operation. See \ref bdldl_status for details.
 */
bdldl_status bdldl_setConstraint(struct ldl_Solver *solver, const ldl_float *E,
                                 const ldl_float *F, const ldl_float *g, int k);

#ifdef __cplusplus
}  // #extern "C"
#endif