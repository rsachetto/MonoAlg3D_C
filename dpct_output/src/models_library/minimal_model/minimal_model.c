#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "minimal_model.h"
#include <stdio.h>

GET_CELL_MODEL_DATA(init_cell_model_data) {

   assert(cell_model);

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;
}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    log_info("Using Minimal Model Mixed CPU model\n");

    uint32_t num_cells = solver->original_num_cells;
    solver->sv = (real*)malloc(NEQ*num_cells*sizeof(real));

    OMP(parallel for)
        for(uint32_t i = 0; i < num_cells; i++) {
            real *sv = &solver->sv[i * NEQ];
            sv[0] = 0.0f; //u
            sv[1] = 1.0f; //v
            sv[2] = 1.0f; //w
            sv[3] = 0.0f; //s
        }
}

SOLVE_MODEL_ODES(solve_model_odes_cpu) {

    uint32_t sv_id;

    uint32_t num_cells_to_solve = ode_solver->num_cells_to_solve;
    size_t num_cells_to_solve_size = num_cells_to_solve*sizeof(uint32_t);
    uint32_t * cells_to_solve = ode_solver->cells_to_solve;
    real *sv = ode_solver->sv;
    real dt = ode_solver->min_dt;
    uint32_t num_steps = ode_solver->num_steps;

    int i;

    // Mapping array
    uint32_t *mapping = NULL;
    
    // Extra data function will tag the cells in the grid
    if (ode_solver->ode_extra_data) {
        mapping = (uint32_t*)ode_solver->ode_extra_data; 
    
        OMP(parallel for private(sv_id))
        for (i = 0; i < num_cells_to_solve; i++) {
            if(cells_to_solve)
                sv_id = cells_to_solve[i];
            else
                sv_id = i;

            for (int j = 0; j < num_steps; ++j) {
                solve_model_ode_cpu(dt, sv + (sv_id * NEQ), stim_currents[i], mapping[i]);
            }
        }
    }
    // Default: All cells are ENDO type
    else {
        OMP(parallel for private(sv_id))
        for (i = 0; i < num_cells_to_solve; i++) {
            if(cells_to_solve)
                sv_id = cells_to_solve[i];
            else
                sv_id = i;

            for (int j = 0; j < num_steps; ++j) {
                solve_model_ode_cpu(dt, sv + (sv_id * NEQ), stim_currents[i], 0);
            }
        }
    }
}


void solve_model_ode_cpu(real dt, real *sv, real stim_current, int type_cell)  {

    assert(sv);

    real rY[NEQ], rDY[NEQ];

    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    RHS_cpu(rY, rDY, stim_current, dt, type_cell);

    // Transmembrane potencial is solved using Explicit Euler
    sv[0] = dt*rDY[0] + sv[0];

    // The other gate variables are solved using Rush-Larsen
    for(int i = 1; i < NEQ; i++) {
        sv[i] = rDY[i];
    }
}

void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt, int type_cell) {
    const real u   = sv[0];
    const real v   = sv[1];
    const real w   = sv[2];
    const real s   = sv[3];

    const real u_o = 0.0;
    const real theta_v = 0.3;
    const real theta_w = 0.13;
    const real tau_vplus = 1.4506;
    const real tau_s1 = 2.7342;
    const real k_s = 2.0994;
    const real u_s = 0.9087;

    real u_u, theta_vminus, theta_o, tau_v1minus, tau_v2minus, tau_w1minus, tau_w2minus;
    real k_wminus, u_wminus, tau_wplus, tau_fi, tau_o1, tau_o2, tau_so1, tau_so2;
    real k_so, u_so, tau_s2, tau_si, tau_winf, w_infstar;

    if (type_cell == 0) {        // ENDO
        u_u = 1.56; theta_vminus = 0.2; theta_o = 0.006; tau_v1minus = 75.0; tau_v2minus = 10.0;
        tau_w1minus = 6.0; tau_w2minus = 140.0; k_wminus = 200.0; u_wminus = 0.016;
        tau_wplus = 280.0; tau_fi = 1.5*0.1; tau_o1 = 470.0; tau_o2 = 6.0; tau_so1 = 40.0;
        tau_so2 = 1.2; k_so = 2.0; u_so = 0.65; tau_s2 = 2.0; tau_si = 2.9013;
        tau_winf = 0.0273; w_infstar = 0.78;
    } else if (type_cell == 1) { // MYO
        u_u = 1.61; theta_vminus = 0.1; theta_o = 0.005; tau_v1minus = 80.0; tau_v2minus = 1.4506;
        tau_w1minus = 70.0; tau_w2minus = 8.0; k_wminus = 200.0; u_wminus = 0.016;
        tau_wplus = 280.0; tau_fi = 1.5*0.078; tau_o1 = 410.0; tau_o2 = 7.0; tau_so1 = 91.0;
        tau_so2 = 0.8; k_so = 2.1; u_so = 0.6; tau_s2 = 4.0; tau_si = 3.3849;
        tau_winf = 0.01; w_infstar = 0.5;
    } else {                    // EPI
        u_u = 1.55; theta_vminus = 0.006; theta_o = 0.006; tau_v1minus = 60.0; tau_v2minus = 1150.0;
        tau_w1minus = 60.0; tau_w2minus = 15.0; k_wminus = 65.0; u_wminus = 0.03;
        tau_wplus = 200.0; tau_fi = 1.5*0.11; tau_o1 = 400.0; tau_o2 = 6.0; tau_so1 = 30.0181;
        tau_so2 = 0.9957; k_so = 2.0458; u_so = 0.65; tau_s2 = 16.0; tau_si = 1.8875;
        tau_winf = 0.07; w_infstar = 0.94;
    }

    real H = (u - theta_v > 0) ? 1.0 : 0.0;
    real h_o = (u - theta_o > 0) ? 1.0 : 0.0;
    real h_w = (u - theta_w > 0) ? 1.0 : 0.0;
    real h_v_minus = (u - theta_vminus > 0) ? 1.0 : 0.0;

    real tau_o = (1.0 - h_o) * tau_o1 + h_o * tau_o2;
    real tau_so = tau_so1 + (tau_so2 - tau_so1) * (1.0 + tanh(k_so * (u - u_so))) * 0.5;
    real tau_vminus = (1.0 - h_v_minus) * tau_v1minus + h_v_minus * tau_v2minus;

    real J_fi = -v * H * (u - theta_v) * (u_u - u) / tau_fi;
    real J_so = (u - u_o) * (1.0 - h_w) / tau_o + h_w / tau_so;
    real J_si = -h_w * w * s / tau_si;

    rDY_[0] = -(J_fi + J_so + J_si) + stim_current;

    real v_inf = (u < theta_vminus) ? 1.0 : 0.0;
    real tau_v_rl = (tau_vplus * tau_vminus) / (tau_vplus - tau_vplus * H + tau_vminus * H);
    real v_inf_rl = (tau_vplus * v_inf * (1 - H)) / (tau_vplus - tau_vplus * H + tau_vminus * H);

    if (tau_v_rl > 1e-10) {
        rDY_[1] = v_inf_rl - (v_inf_rl - v) * exp(-dt / tau_v_rl);
    } else {
        rDY_[1] = dt * ((1.0 - H) * (v_inf - v) / tau_vminus - H * v / tau_vplus) + v;
    }

    real w_inf = (1.0 - h_o) * (1.0 - (u / tau_winf)) + h_o * w_infstar;
    real tau_wminus = tau_w1minus + (tau_w2minus - tau_w1minus) * (1.0 + tanh(k_wminus * (u - u_wminus))) * 0.5;
    real tau_w_rl = (tau_wplus * tau_wminus) / (tau_wplus - tau_wplus * h_w + tau_wminus * h_w);
    real w_inf_rl = (tau_wplus * w_inf * (1 - h_w)) / (tau_wplus - tau_wplus * h_w + tau_wminus * h_w);

    if (tau_w_rl > 1e-10) {
        rDY_[2] = w_inf_rl - (w_inf_rl - w) * exp(-dt / tau_w_rl);
    } else {
        rDY_[2] = dt * ((1.0 - h_w) * (w_inf - w) / tau_wminus - h_w * w / tau_wplus) + w;
    }

    real tau_s = (1.0 - h_w) * tau_s1 + h_w * tau_s2;
    real s_inf_rl = (1.0 + tanh(k_s * (u - u_s))) / 2;

    if (tau_s > 1e-10) {
        rDY_[3] = s_inf_rl - (s_inf_rl - s) * exp(-dt / tau_s);
    } else {
        rDY_[3] = dt * (((1.0 + tanh(k_s * (u - u_s))) * 0.5 - s) / tau_s) + s;
    }
}
