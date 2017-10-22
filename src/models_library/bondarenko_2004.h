#include "model_common.h"
#include <assert.h>
#include <stdlib.h>

#define NEQ 41
#define INITIAL_V (-82.4202f)

void solve_model_ode_cpu(real dt, real *sv, real stim_current);
void RHS_cpu(const real *sv, real *rDY_, real stim_current);