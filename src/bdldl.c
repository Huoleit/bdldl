#include "bdldl/bdldl.h"

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

struct ldl_Solver {
  struct qp_size *size;
  struct qp_size *data;
  struct ldl_workspace *ws;
};
