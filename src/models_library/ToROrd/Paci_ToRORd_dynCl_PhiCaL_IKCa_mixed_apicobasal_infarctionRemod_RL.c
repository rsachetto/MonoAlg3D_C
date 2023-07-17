#include "Paci_ToRORd_dynCl_PhiCaL_IKCa_mixed_apicobasal_infarctionRemod_RL.h"
#include <stdlib.h>
#include <stdio.h>

GET_CELL_MODEL_DATA(init_cell_model_data) {
    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;
}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    log_info("Using Paci_ToRORd_dynCl_PhiCaL_IKCa_mixed_apicobasal_infarctionRemod CPU model\n");

    uint32_t num_cells = solver->original_num_cells;
    solver->sv = (real*)malloc(NEQ*num_cells*sizeof(real));
    bool adpt = solver->adaptive;

    if(adpt) {
        solver->ode_dt = (real*)malloc(num_cells*sizeof(real));

        OMP(parallel for)
        for(int i = 0; i < num_cells; i++) {
            solver->ode_dt[i] = solver->min_dt;
        }

        solver->ode_previous_dt = (real*)calloc(num_cells, sizeof(real));
        solver->ode_time_new    = (real*)calloc(num_cells, sizeof(real));
        log_info("Using Adaptive Euler model to solve the ODEs\n");
    } else {
        log_info("Using Euler model to solve the ODEs\n");
    }

    real *extra_data = NULL;
    if(solver->ode_extra_data) {
        extra_data = (real*)solver->ode_extra_data;
    } else {
        log_error_and_exit("You need to specify a mask function when using a mixed model!\n");
    }

    int infarct_stage = (int) extra_data[5 * num_cells];

    OMP(parallel for)
    for(uint32_t i = 0; i < num_cells; i++) {
        int is_paci = (int) extra_data[i];
        int layer = (int) extra_data[i + num_cells];
        int infarct_zone = (int) extra_data[i + (2*num_cells)];
        real *sv = &solver->sv[i * NEQ];

        #include "Paci_ToRORd_dynCl_PhiCaL_IKCa_mixed_apicobasal_infarctionRemod_SS.common.c"
        if(is_paci > 0) {
            sv[0] = v;             // V;         millivolt
            sv[1] = caSR;          // Ca_SR;     millimolar
            sv[2] = cai;        // Ca_i;      millimolar
            sv[3] = d;       // d;         dimensionless
            sv[4] = f1;          // f1;        dimensionless
            sv[5] = f2;          // f2;        dimensionless
            sv[6] = fCa;          // fCa;       dimensionless
            sv[7] = Xr1;        // Xr1;       dimensionless
            sv[8] = Xr2;          // Xr2;       dimensionless
            sv[9] = Xs;         // Xs;        dimensionless
            sv[10] = h;         // h;         dimensionless
            sv[11] = j;         // j;         dimensionless
            sv[12] = m;        // m;         dimensionless
            sv[13] = Xf;         // Xf;        dimensionless
            sv[14] = q;         // q;         dimensionless
            sv[15] = r;       // r;         dimensionless
            sv[16] = nai;          // Na_i;      millimolar
            sv[17] = mL;       // m_L;       dimensionless
            sv[18] = hL;         // h_L;       dimensionless
            sv[19] = RyRa;        // RyRa;      dimensionless
            sv[20] = RyRo;      // RyRo;      dimensionless
            sv[21] = RyRc;         // RyRc;      dimensionless
            sv[22] = 0.0;
            sv[23] = 0.0;
            sv[24] = 0.0;
            sv[25] = 0.0;
            sv[26] = 0.0;
            sv[27] = 0.0;
            sv[28] = 0.0;
            sv[29] = 0.0;
            sv[30] = 0.0;
            sv[31] = 0.0;
            sv[32] = 0.0;
            sv[33] = 0.0;
            sv[34] = 0.0;
            sv[35] = 0.0;
            sv[36] = 0.0;
            sv[37] = 0.0;
            sv[38] = 0.0;
            sv[39] = 0.0;
            sv[40] = 0.0;
            sv[41] = 0.0;
            sv[42] = 0.0;
            sv[43] = 0.0;
            sv[44] = 0.0;
        } else if (is_paci == 0) {
            sv[ 0]= v;
            sv[ 1]= nai;
            sv[ 2]= nass;
            sv[ 3]= ki;
            sv[ 4]= kss;
            sv[ 5]= cai;
            sv[ 6]= cass;
            sv[ 7]= cansr;
            sv[ 8]= cajsr;
            sv[ 9]= m;
            sv[10]= hp;
            sv[11]= h;
            sv[12]= j;
            sv[13]= jp;
            sv[14]= mL;
            sv[15]= hL;
            sv[16]= hLp;
            sv[17]= a;
            sv[18]= iF;
            sv[19]= iS;
            sv[20]= ap;
            sv[21]= iFp;
            sv[22]= iSp;
            sv[23]= d;
            sv[24]= ff;
            sv[25]= fs;
            sv[26]= fcaf;
            sv[27]= fcas;
            sv[28]= jca;
            sv[29]= nca;
            sv[30]= nca_i;
            sv[31]= ffp;
            sv[32]= fcafp;
            sv[33]= xs1;
            sv[34]= xs2;
            sv[35]= Jrel_np;
            sv[36]= CaMKt;
            sv[37]= ikr_c0;
            sv[38]= ikr_c1;
            sv[39]= ikr_c2;
            sv[40]= ikr_o;
            sv[41]= ikr_i;
            sv[42]= Jrel_p;
            sv[43]= cli;
            sv[44]= clss;
        }
    }
}

SOLVE_MODEL_ODES(solve_model_odes_cpu) {

    uint32_t sv_id;

    size_t num_cells_to_solve = ode_solver->num_cells_to_solve;
    uint32_t * cells_to_solve = ode_solver->cells_to_solve;
    real *sv = ode_solver->sv;
    real dt = ode_solver->min_dt;
    uint32_t num_steps = ode_solver->num_steps;
    bool adpt = ode_solver->adaptive;

    real *extra_data = NULL;
    if(ode_solver->ode_extra_data) {
        extra_data = (real*)ode_solver->ode_extra_data;
    } else {
        log_error_and_exit("You need to specify a mask function when using this mixed model!\n");
    }

    int infarct_stage = (int) extra_data[5 * num_cells_to_solve];
    real current_scaling[45];
    for (int i = 1; i < 45; i++) {
        current_scaling[i] = (real) extra_data[5 * num_cells_to_solve + i];
    }

    #pragma omp parallel for private(sv_id)
    for (int i = 0; i < num_cells_to_solve; i++) {
        int is_paci = (int) extra_data[i];
        int layer = (int) extra_data[i + num_cells_to_solve];
        int infarct_zone = (int) extra_data[i + (2*num_cells_to_solve)];
        real apicobasal = extra_data[i + (3*num_cells_to_solve)];

        if(cells_to_solve)
            sv_id = cells_to_solve[i];
        else
            sv_id = i;

        if(adpt) {
            solve_RL_cpu_adpt(sv + (sv_id * NEQ), stim_currents[i], current_t + dt, sv_id, ode_solver, is_paci, layer, infarct_zone, infarct_stage, apicobasal, current_scaling);
            //solve_forward_euler_cpu_adpt(sv + (sv_id * NEQ), stim_currents[i], current_t + dt, sv_id, ode_solver, is_paci, layer, infarct_zone, infarct_stage, apicobasal, current_scaling);
        } else {
            for (int j = 0; j < num_steps; ++j) {
                solve_model_ode_cpu(dt, sv + (sv_id * NEQ), stim_currents[i], is_paci, layer, infarct_zone, infarct_stage, apicobasal, current_scaling);
            }
        }
    }
}

void solve_model_ode_cpu(real dt, real *sv, real stim_current, int is_paci, int layer, int infarct_zone, int infarct_stage, real apicobasal, real *current_scaling)  {
    if((infarct_zone != 1) || (infarct_stage != 3)) {
        const real TOLERANCE = 1e-8;
        real rY[NEQ], rDY[NEQ];

        for(int i = 0; i < NEQ; i++)
            rY[i] = sv[i];

        real a[NEQ], b[NEQ];
        RHS_cpu(a, b, sv, rDY, stim_current, dt, is_paci, layer, infarct_zone, infarct_stage, apicobasal, current_scaling);


        bool is_rush_larsen[NEQ];
        //for (int i = 0; i < NEQ; i++)
        //       is_rush_larsen[i] = false;
        if((is_paci > 0) && (is_paci < 4)) {
           //for (int i = 0; i < NEQ; i++)
            //    is_rush_larsen[i] = false;
           for (int i = 0; i < NEQ_PACI; i++)
                is_rush_larsen[i] = true;
           for (int i = NEQ_PACI; i < NEQ; i++)
                is_rush_larsen[i] = false;
           is_rush_larsen[0] = false;
           is_rush_larsen[1] = false;
           is_rush_larsen[2] = false;
           is_rush_larsen[16] = false;
        } else {
           for (int i = 0; i < 9; i++)
               is_rush_larsen[i] = false;
           for (int i = 9; i < 29; i++)
               is_rush_larsen[i] = true;
           is_rush_larsen[29] = false;
           is_rush_larsen[30] = false;
           for (int i = 31; i < 36; i++)
               is_rush_larsen[i] = true;
           for (int i = 36; i < 42; i++)
               is_rush_larsen[i] = false;
           is_rush_larsen[42] = true;
           is_rush_larsen[43] = false;
           is_rush_larsen[44] = false;
        }

        for (int i = 0; i < NEQ; i++) {
            if (is_rush_larsen[i]) {
                SOLVE_EQUATION_RUSH_LARSEN_CPU(i);
            }
            else {
                SOLVE_EQUATION_EULER_CPU(i);
            }
        }
    }
    // if(is_paci == 1) {
    //     for(int i = 0; i < NEQ_PACI; i++)
    //         sv[i] = dt*rDY[i] + rY[i];
    // } else if(is_paci == 0) {
    //     for(int i = 0; i < NEQ; i++)
    //         sv[i] = dt*rDY[i] + rY[i];
    // }
}

void solve_forward_euler_cpu_adpt(real *sv, real stim_curr, real final_time, int sv_id, struct ode_solver *solver, int is_paci, int layer, int infarct_zone, int infarct_stage, real apicobasal, real *current_scaling) {
    const real _beta_safety_ = 0.8;

    int numEDO;
    if(is_paci > 0) {
        numEDO = NEQ_PACI;
    } else {
        numEDO = NEQ;
    }

    real rDY[numEDO];

    real _tolerances_[numEDO];
    real _aux_tol = 0.0;
    // initializes the variables
    solver->ode_previous_dt[sv_id] = solver->ode_dt[sv_id];

    real edos_old_aux_[numEDO];
    real edos_new_euler_[numEDO];
    real *_k1__ = (real *)malloc(sizeof(real) * numEDO);
    real *_k2__ = (real *)malloc(sizeof(real) * numEDO);
    real *a_ = (real *)malloc(sizeof(real) * numEDO);
    real *b_ = (real *)malloc(sizeof(real) * numEDO);
    real *a_new = (real *)malloc(sizeof(real) * numEDO);
    real *b_new = (real *)malloc(sizeof(real) * numEDO);
    real *_k_aux__;

    real *dt = &solver->ode_dt[sv_id];
    real *time_new = &solver->ode_time_new[sv_id];
    real *previous_dt = &solver->ode_previous_dt[sv_id];

    // Keep 'dt' inside the adaptive interval
    if(*time_new + *dt > final_time) {
        *dt = final_time - *time_new;
    }

    RHS_cpu(a_, b_, sv, rDY, stim_curr, *dt, is_paci, layer, infarct_zone, infarct_stage, apicobasal, current_scaling);
    *time_new += *dt;

    for(int i = 0; i < numEDO; i++) {
        _k1__[i] = rDY[i];
    }

    const real rel_tol = solver->rel_tol;
    const real abs_tol = solver->abs_tol;

    const real __tiny_ = powf(abs_tol, 2.0);

    real min_dt = solver->min_dt;
    real max_dt = solver->max_dt;

    while(1) {

        for(int i = 0; i < numEDO; i++) {
            // stores the old variables in a vector
            edos_old_aux_[i] = sv[i];
            // computes euler method
            edos_new_euler_[i] = _k1__[i] * *dt + edos_old_aux_[i];
            // steps ahead to compute the rk2 method
            sv[i] = edos_new_euler_[i];
        }

        *time_new += *dt;
        RHS_cpu(a_new, b_new, sv, rDY, stim_curr, *dt, is_paci, layer, infarct_zone, infarct_stage, apicobasal, current_scaling);
        *time_new -= *dt; // step back

        double greatestError = 0.0, auxError = 0.0;
        for(int i = 0; i < numEDO; i++) {
            // stores the new evaluation
            _k2__[i] = rDY[i];
            _aux_tol = fabs(edos_new_euler_[i]) * rel_tol;
            _tolerances_[i] = (abs_tol > _aux_tol) ? abs_tol : _aux_tol;
            // finds the greatest error between  the steps
            auxError = fabs(((*dt / 2.0) * (_k1__[i] - _k2__[i])) / _tolerances_[i]);
            greatestError = (auxError > greatestError) ? auxError : greatestError;
        }
        /// adapt the time step
        greatestError += __tiny_;
        *previous_dt = *dt;
        /// adapt the time step
        *dt = _beta_safety_ * (*dt) * sqrtf(1.0f / greatestError);

        if(*dt < min_dt) {
            *dt = min_dt;
        } else if(*dt > max_dt) {
            *dt = max_dt;
        }

        if(*time_new + *dt > final_time) {
            *dt = final_time - *time_new;
        }

        // it doesn't accept the solution
        if(greatestError >= 1.0f && *dt > min_dt) {
            // restore the old values to do it again
            for(int i = 0; i < numEDO; i++) {
                sv[i] = edos_old_aux_[i];
            }
            // throw the results away and compute again
        } else {
            // it accepts the solutions
            if(greatestError >= 1.0) {
                printf("Accepting solution with error > %lf \n", greatestError);
            }

            _k_aux__ = _k2__;
            _k2__ = _k1__;
            _k1__ = _k_aux__;

            // it steps the method ahead, with euler solution
            for(int i = 0; i < numEDO; i++) {
                sv[i] = edos_new_euler_[i];
            }

            if(*time_new + *previous_dt >= final_time) {
                if(final_time == *time_new) {
                    break;
                } else if(*time_new < final_time) {
                    *dt = *previous_dt = final_time - *time_new;
                    *time_new += *previous_dt;
                    break;
                }
            } else {
                *time_new += *previous_dt;
            }
        }
    }

    free(_k1__);
    free(_k2__);
}

void solve_RL_cpu_adpt(real *sv, real stim_curr, real final_time, int sv_id, struct ode_solver *solver, int is_paci, int layer, int infarct_zone, int infarct_stage, real apicobasal, real *current_scaling) {
    int numEDO = NEQ;
    real rDY[numEDO];
    solver->ode_previous_dt[sv_id] = solver->ode_dt[sv_id];

    real edos_old_aux_[numEDO];
    real edos_new_euler_[numEDO];
    real *_k1__ = (real *)malloc(sizeof(real) * numEDO);
    real *_k2__ = (real *)malloc(sizeof(real) * numEDO);
    real *a_ = (real *)malloc(sizeof(real) * numEDO);
    real *b_ = (real *)malloc(sizeof(real) * numEDO);
    real *a_new = (real *)malloc(sizeof(real) * numEDO);
    real *b_new = (real *)malloc(sizeof(real) * numEDO);
    real *_k_aux__, *_a_aux__, *_b_aux__;

    real *dt = &solver->ode_dt[sv_id];
    real *time_new = &solver->ode_time_new[sv_id];
    real *previous_dt = &solver->ode_previous_dt[sv_id];

    // Keep 'dt' inside the adaptive interval
    if(*time_new + *dt > final_time) {
        *dt = final_time - *time_new;
    }

    RHS_cpu(a_, b_, sv, rDY, stim_curr, *dt, is_paci, layer, infarct_zone, infarct_stage, apicobasal, current_scaling);
    *time_new += *dt;

    for(int i = 0; i < numEDO; i++) {
        _k1__[i] = rDY[i];
    }

    const real rel_tol = solver->rel_tol;
    const real abs_tol = solver->abs_tol;
    const real __tiny_ = powf(abs_tol, 2.0);

    real min_dt = solver->min_dt;
    real max_dt = solver->max_dt;

    bool is_rush_larsen[NEQ];
    if((is_paci > 0) && (is_paci < 4)) {
        //for (int i = 0; i < NEQ; i++)
          //   is_rush_larsen[i] = false;
         for (int i = 0; i < NEQ_PACI; i++)
              is_rush_larsen[i] = true;
         for (int i = NEQ_PACI; i < NEQ; i++)
              is_rush_larsen[i] = false;
         is_rush_larsen[0] = false;
         is_rush_larsen[1] = false;
         is_rush_larsen[2] = false;
         is_rush_larsen[16] = false;
    } else {
        for (int i = 0; i < 9; i++)
            is_rush_larsen[i] = false;
        for (int i = 9; i < 29; i++)
            is_rush_larsen[i] = true;
        is_rush_larsen[29] = false;
        is_rush_larsen[30] = false;
        for (int i = 31; i < 36; i++)
            is_rush_larsen[i] = true;
        for (int i = 36; i < 42; i++)
            is_rush_larsen[i] = false;
        is_rush_larsen[42] = true;
        is_rush_larsen[43] = false;
        is_rush_larsen[44] = false;
    }

    while(1) {

        for (int i = 0; i < NEQ; i++) {
            if (is_rush_larsen[i]) {
                SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(i);
            }
            else {
                SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(i);
            }
        }

        *time_new += *dt;
        RHS_cpu(a_new, b_new, sv, rDY, stim_curr, *dt, is_paci, layer, infarct_zone, infarct_stage, apicobasal, current_scaling);
        *time_new -= *dt; // step back

        double greatestError = 0.0, auxError = 0.0;
        real as, bs, f, y_2nd_order;
        for (int i = 0; i < NEQ; i++) {
            if (is_rush_larsen[i]) {
                SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(i);
            }
            else {
                SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(i);
            }
        }

        /// adapt the time step
        greatestError += __tiny_;
        *previous_dt = *dt;
        /// adapt the time step
        *dt = (*dt) * sqrt(0.5 * rel_tol / greatestError);            // Jhonny`s formula

        if(*dt < min_dt) {
            *dt = min_dt;
        } else if(*dt > max_dt) {
            *dt = max_dt;
        }

        if(*time_new + *dt > final_time) {
            *dt = final_time - *time_new;
        }

        // it doesn't accept the solution
        if(greatestError >= 1.0f && *dt > min_dt) {
            // restore the old values to do it again
            for(int i = 0; i < numEDO; i++) {
                sv[i] = edos_old_aux_[i];
            }
            // throw the results away and compute again
        } else {
            // it accepts the solutions
            if(greatestError >= 1.0) {
                printf("Accepting solution with error > %lf \n", greatestError);
            }

            _k_aux__ = _k2__;
            _k2__ = _k1__;
            _k1__ = _k_aux__;

            _a_aux__ = a_;
            a_ = a_new;
            a_new = _a_aux__;

            _b_aux__ = b_;
            b_ = b_new;
            b_new = _b_aux__;

            // it steps the method ahead, with euler solution
            for(int i = 0; i < numEDO; i++) {
                sv[i] = edos_new_euler_[i];
            }

            if(*time_new + *previous_dt >= final_time) {
                if(final_time == *time_new) {
                    break;
                } else if(*time_new < final_time) {
                    *dt = *previous_dt = final_time - *time_new;
                    *time_new += *previous_dt;
                    break;
                }
            } else {
                *time_new += *previous_dt;
            }
        }
    }

    free(_k1__);
    free(_k2__);
    free(a_);
    free(b_);
    free(a_new);
    free(b_new);
}

void RHS_cpu(real *a_, real *b_, const real *sv, real *rDY, real stim_current, real dt, int is_paci, int layer, int infarct_zone, int infarct_stage, real apicobasal, real *current_scaling) {
    // State variables
    real v;
    real nai;
    real nass;
    real ki;
    real kss;
    real cai;
    real cass;
    real cansr;
    real cajsr;
    real m;
    real hp;
    real h;
    real j;
    real jp;
    real mL;
    real hL;
    real hLp;
    real a;
    real iF;
    real iS;
    real ap;
    real iFp;
    real iSp;
    real d;
    real ff;
    real fs;
    real fcaf;
    real fcas;
    real jca;
    real nca;
    real nca_i;
    real ffp;
    real fcafp;
    real xs1;
    real xs2;
    real Jrel_np;
    real CaMKt;
    real ikr_c0;
    real ikr_c1;
    real ikr_c2;
    real ikr_o;
    real ikr_i;
    real Jrel_p;
    real caSR;
    real f1;
    real f2;
    real fCa;
    real Xr1;
    real Xr2;
    real Xs;
    real Xf;
    real q;
    real r;
    real RyRa;
    real RyRo;
    real RyRc;
    real cli;
    real clss;

    real INa_mult = current_scaling[1];
    real IKr_mult = current_scaling[2];
    real ICaL_mult = current_scaling[3];
    real INaL_mult = current_scaling[4];
    real IKs_mult = current_scaling[5];
    real Ito_mult = current_scaling[6];
    real IK1_mult = current_scaling[7];

    real INa_mult_RZ = current_scaling[8];
    real IKr_mult_RZ = current_scaling[9];
    real ICaL_mult_RZ = current_scaling[10];
    real INaL_mult_RZ = current_scaling[11];
    real IKs_mult_RZ = current_scaling[12];
    real Ito_mult_RZ = current_scaling[13];
    real IK1_mult_RZ = current_scaling[14];
    real aCaMK_mult_RZ = current_scaling[15];
    real tau_relp_mult_RZ = current_scaling[16];
    real ICab_mult_RZ = current_scaling[17];
    real Iup_mult_RZ = current_scaling[18];
    real IKCa_mult_RZ = current_scaling[19];
    real IClCa_mult_RZ = current_scaling[20];

    real INa_mult_IZ = current_scaling[21];
    real IKr_mult_IZ = current_scaling[22];
    real ICaL_mult_IZ = current_scaling[23];
    real INaL_mult_IZ = current_scaling[24];
    real IKs_mult_IZ = current_scaling[25];
    real Ito_mult_IZ = current_scaling[26];
    real IK1_mult_IZ = current_scaling[27];
    real aCaMK_mult_IZ = current_scaling[28];
    real tau_relp_mult_IZ = current_scaling[29];
    real ICab_mult_IZ = current_scaling[30];
    real Ko_mult_IZ = current_scaling[31];

    real INa_mult_BZ = current_scaling[32];
    real IKr_mult_BZ = current_scaling[33];
    real ICaL_mult_BZ = current_scaling[34];
    real INaL_mult_BZ = current_scaling[35];
    real IKs_mult_BZ = current_scaling[36];
    real Ito_mult_BZ = current_scaling[37];
    real IK1_mult_BZ = current_scaling[38];
    real aCaMK_mult_BZ = current_scaling[39];
    real tau_relp_mult_BZ = current_scaling[40];
    real ICab_mult_BZ = current_scaling[41];
    real Iup_mult_BZ = current_scaling[42];
    real IKCa_mult_BZ = current_scaling[43];
    real IClCa_mult_BZ = current_scaling[44];

    real atrial_INa = 0.5698;
    real atrial_INaL = 3.4781;
    real atrial_Ito = 3.1498;
    real atrial_IKr = 4.2080;
    real atrial_IKs = 0.6331;
    real atrial_IK1 = 0.3464;
    real atrial_ICaL = 1.7471;
    real atrial_INaCa = 1.5376;
    real atrial_INaK = 1.1698;
    real atrial_Iup = 2.9601;
    real atrial_If = 1.1746;
    real atrial_Irel = 0.0943;
    real atrial_Ileak = 0.1118;

    real nodal_INa = 1.1739;
    real nodal_INaL = 1.8433;
    real nodal_Ito = 2.2774;
    real nodal_IKr = 3.15;
    real nodal_IKs = 0.57;
    real nodal_IK1 = 0.01;
    real nodal_ICaL = 0.95;
    real nodal_INaCa = 1.1768;
    real nodal_INaK = 0.95;
    real nodal_Iup = 0.05;
    real nodal_If = 1.25;
    real nodal_Irel = 1.2028;
    real nodal_Ileak = 2.9065;

    if(is_paci > 0) {
        v = sv[0];
        caSR = sv[1];
        cai = sv[2];
        d = sv[3];
        f1 = sv[4];
        f2 = sv[5];
        fCa = sv[6];
        Xr1 = sv[7];
        Xr2 = sv[8];
        Xs = sv[9];
        h = sv[10];
        j = sv[11];
        m = sv[12];
        Xf = sv[13];
        q = sv[14];
        r = sv[15];
        nai = sv[16];
        mL = sv[17];
        hL = sv[18];
        RyRa = sv[19];
        RyRo = sv[20];
        RyRc = sv[21];
    } else if (is_paci == 0) {
        v = sv[0];
        nai = sv[1];
        nass = sv[2];
        ki = sv[3];
        kss = sv[4];
        cai = sv[5];
        cass = sv[6];
        cansr = sv[7];
        cajsr = sv[8];
        m = sv[9];
        hp = sv[10];
        h = sv[11];
        j = sv[12];
        jp = sv[13];
        mL = sv[14];
        hL = sv[15];
        hLp = sv[16];
        a = sv[17];
        iF = sv[18];
        iS = sv[19];
        ap = sv[20];
        iFp = sv[21];
        iSp = sv[22];
        d = sv[23];
        ff = sv[24];
        fs = sv[25];
        fcaf = sv[26];
        fcas = sv[27];
        jca = sv[28];
        nca = sv[29];
        nca_i = sv[30];
        ffp = sv[31];
        fcafp = sv[32];
        xs1 = sv[33];
        xs2 = sv[34];
        Jrel_np = sv[35];
        CaMKt = sv[36];
        ikr_c0 = sv[37];
        ikr_c1 = sv[38];
        ikr_c2 = sv[39];
        ikr_o = sv[40];
        ikr_i = sv[41];
        Jrel_p = sv[42];
        cli = sv[43];
        clss = sv[44];
    }

    #include "Paci_ToRORd_dynCl_PhiCaL_IKCa_mixed_apicobasal_infarctionRemod_RL.common.c"
}
