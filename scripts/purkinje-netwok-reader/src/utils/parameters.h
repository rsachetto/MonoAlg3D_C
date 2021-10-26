#ifndef MONOALG3D_PARAMETERS_H
#define MONOALG3D_PARAMETERS_H

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

static const uint32_t NUM_THREADS = 6;
static const double DT_PDE = 0.02;
static const double SIMULATION_TIME = 400.0;
static const bool ABORT_ON_NO_ACTIVITY = false;
static const bool USE_ADAPTIVITY = false;
static const bool CALC_ACTIVATION_TIME = true;

static const char* UPDATE_MONODOMAIN_FUNCTION = "update_monodomain_default";

static const uint32_t PRINT_RATE = 100;
static const char* OUTPUT_DIR = "outputs/lucas_cheat_purkinje_coupled";
static const char* SAVE_MESH_FUNCTION = "save_as_vtu";
static const bool SAVE_PVD = true;
static const char* FILE_PREFIX = "V";
static const bool BINARY = false;
static const bool COMPRESS = false;

static const double SIGMA_X = 0.00001334;
static const double SIGMA_Y = 0.0000176;
static const double SIGMA_Z = 0.0000176;
//static const double SIGMA_X = 0.00001334;
//static const double SIGMA_Y = 0.00001334;
//static const double SIGMA_Z = 0.00001334;
static const char* ASSEMBLY_MATRIX_LIBRARY = "shared_libs/libdefault_matrix_assembly.so";
static const char* ASSEMBLY_MATRIX_FUNCTION = "homogeneous_sigma_assembly_matrix";

static const double TOLERANCE = 1e-16;
static const bool USE_PRECONDITIONER = true;
static const uint32_t MAX_ITERATIONS = 200;
static const char* LINEAR_SYSTEM_LIBRARY = "shared_libs/libdefault_linear_system_solver.so";
static const char* LINEAR_SYSTEM_FUNCTION = "conjugate_gradient";

static const double REFINEMENT_BOUND = 0.11;
static const double DEREFINEMENT_BOUND = 0.10;
static const uint32_t REFINE_EACH = 1;
static const uint32_t DEREFINE_EACH = 1;

static const char* DOMAIN_NAME = "Plain Mesh";
static const uint32_t NUM_LAYER = 1;
static const double START_DX = 200.0;
static const double START_DY = 200.0;
static const double START_DZ = 200.0;
static const double SIDE_LENGTH = 20000.0;
static const char* DOMAIN_FUNCTION = "initialize_grid_with_square_mesh";

static const double DT_ODE = 0.02;
static const bool USE_GPU = true;
static const uint32_t GPU_ID = 0;
static const char* ODE_LIBRARY = "shared_libs/libten_tusscher_2006.so";

static const double STIM_DURATION = 0.5;
static const double STIM_CURRENT = -50.0;
static const double STIM_PERIOD = 500.0;
static const char* STIM_FUNCTION = "stim_if_inside_circle_than";

static const double PMJ_RADIUS = 500.0;

#endif