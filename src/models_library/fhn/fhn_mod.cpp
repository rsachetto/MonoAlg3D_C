#include "fhn_mod.h"
#include "../../ode_solver/ode_solver.h"
#include <cstdio>
#include <dpct/dpct.hpp>
#include <sycl/sycl.hpp>

void RHS_sycl(real *Y, real stim_current, real *dY, int sv_id, int num_cells);

extern "C" GET_CELL_MODEL_DATA(init_cell_model_data) {

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;
}

extern "C" SET_ODE_INITIAL_CONDITIONS_GPU(set_model_initial_conditions_gpu) {
    printf("Using FitzHugh-Nagumo 1961 SYCL model\n");

    dpct::device_ext &dev_ct1 = dpct::get_current_device();
    sycl::queue &q_ct1 = dev_ct1.default_queue();
    printf("Running on '%s'\n", q_ct1.get_device().get_info<sycl::info::device::name>().c_str());

    uint32_t num_cells = solver->original_num_cells;

    solver->sv = sycl::malloc_device<real>(num_cells * NEQ * sizeof(real), q_ct1);

    if(solver->sv) {
        const int BLOCK_SIZE = 32;
        const int GRID = (num_cells + BLOCK_SIZE - 1) / BLOCK_SIZE;
        const int TOTAL_THREADS = BLOCK_SIZE * GRID;

        try {
            q_ct1
                .submit([&](sycl::handler &h) {
                    real *sv = solver->sv;
                    h.parallel_for(sycl::nd_range<1>(sycl::range<1>(TOTAL_THREADS), sycl::range<1>(BLOCK_SIZE)), [=](sycl::nd_item<1> item) {
                        int i = item.get_global_id(0);
                        if(i < num_cells) {
                            sv[0 * num_cells + i] = 0.0; // Vm millivolt
                            sv[1 * num_cells + i] = 0.0; // v dimensionless
                        }
                    });
                })
                .wait();
            return 0;
        } catch(sycl::exception &e) {
            printf("SYCL exception: %s\n", e.what());
            return -1;
        }
    }
}

extern "C" SET_ODE_INITIAL_CONDITIONS_GPU(set_model_initial_conditions_cpu) {
    set_model_initial_conditions_gpu(solver, ode_extra_config);
}

extern "C" SOLVE_MODEL_ODES(solve_model_odes_gpu) {

    size_t num_cells_to_solve = ode_solver->num_cells_to_solve;
    uint32_t *cells_to_solve = ode_solver->cells_to_solve;
    real *sv = ode_solver->sv;
    real dt = ode_solver->min_dt;
    uint32_t num_steps = ode_solver->num_steps;

    dpct::device_ext &dev_ct1 = dpct::get_current_device();
    sycl::queue &q_ct1 = dev_ct1.default_queue();

    // Using Unifed Shared Memory (USM) instead of access buffer to improve performance
    real *d_stim = sycl::malloc_device<real>(num_cells_to_solve, q_ct1);
    uint32_t *d_cells_to_solve = NULL;
    if(cells_to_solve) {
        d_cells_to_solve = sycl::malloc_device<uint32_t>(num_cells_to_solve, q_ct1);
        q_ct1.memcpy(d_cells_to_solve, cells_to_solve, num_cells_to_solve * sizeof(uint32_t));
    }

    // Copy initial data to device
    q_ct1.memcpy(d_stim, stim_currents, num_cells_to_solve * sizeof(real)).wait();

    // Define block and grid sizes
    const int BLOCK_SIZE = 32;
    const int GRID = (num_cells_to_solve + BLOCK_SIZE - 1) / BLOCK_SIZE;

    try {
        q_ct1
            .submit([&](sycl::handler &hand) {
                // num_cells_to_solve -> 'i'
                hand.parallel_for(sycl::nd_range<1>(sycl::range<1>(GRID * BLOCK_SIZE), sycl::range<1>(BLOCK_SIZE)),
                                  [=](sycl::nd_item<1> item) [[sycl::reqd_work_group_size(32)]] {
                                      size_t i = item.get_global_id(0);
                                      if(i >= num_cells_to_solve)
                                          return;

                                      int sv_id = (d_cells_to_solve) ? d_cells_to_solve[i] : i;

                                      // Private memory for fast access
                                      real rDY[NEQ];

                                      for(int n = 0; n < num_steps; ++n) {
                                          RHS_sycl(sv, d_stim[i], rDY, sv_id, num_cells_to_solve);

                                          for(int j = 0; j < NEQ; j++) {
                                              sv[j * num_cells_to_solve + i] += dt * rDY[j];
                                          }
                                      }
                                  });
            })
            .wait();

        // Free device memory
        sycl::free(d_stim, q_ct1);
        if(d_cells_to_solve)
            sycl::free(d_cells_to_solve, q_ct1);

    } catch(sycl::exception &e) {
        printf("SYCL exception: %s\n", e.what());
    }
}

extern "C" SOLVE_MODEL_ODES(solve_model_odes_cpu) {
    solve_model_odes_gpu(ode_solver, ode_extra_config, current_t, stim_currents);
}

inline void RHS_sycl(real *Y, real stim_current, real *dY, int sv_id, int num_cells) {
    // State variables
    const real u = Y[0 * num_cells + sv_id];
    const real v = Y[1 * num_cells + sv_id];

    const real a = 0.2f;
    const real b = 0.5f;
    const real k = 36.0;
    const real epsilon = 0.00040;

    dY[0] = k * (u * (1.0f - u) * (u - a) - u * v) + stim_current;
    dY[1] = k * epsilon * (b * u - v);
}
