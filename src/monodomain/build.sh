MONODOMAIN_SOURCE_FILES="ode_solver.c monodomain_solver.c"
MONODOMAIN_HEADER_FILES="ode_solver.h monodomain_solver.h"

COMPILE_STATIC_LIB "solvers" "$MONODOMAIN_SOURCE_FILES" "$MONODOMAIN_HEADER_FILES"
