#include "trovato_2020.h"
#include <stdlib.h>

GET_CELL_MODEL_DATA(init_cell_model_data) {

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ;
}

SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    log_info("Using Trovato_2020 CPU model\n");

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

    OMP(parallel for)
    for(uint32_t i = 0; i < num_cells; i++) {
        
        real *sv = &solver->sv[i * NEQ];

        // Default initial conditions (1000 beats at BCL=1000ms) 
        sv[0]       = -86.550102957989600;
        sv[1]       = 0.005060490773142;
        sv[2]       = 1.017658486729359e-04;
        sv[3]       = 8.231857731510893;     
        sv[4]       = 8.231553454361393;     
        sv[5]       = 8.231561539013950;
        sv[6]       = 1.437673447706863e+02;  
        sv[7]       = 1.437677533510394e+02;
        sv[8]       = 1.437677771226899e+02;  
        sv[9]       = 4.360080908582633e-05; 
        sv[10]      = 1.020101597047671e-04; 
        sv[11]      = 1.263525645841406;   
        sv[12]      = 1.248146625349512;    
        sv[13]      = 1.265185617721750;
        sv[14]      = 0;                    
        sv[15]      = 0;
        sv[16]      = 0.006341207769833;       
        sv[17]      = 0.788541761218318;       
        sv[18]      = 0.788474863764949;
        sv[19]      = 0.790412100577539;       
        sv[20]      = 0.579594138981772;      
        sv[21]      = 0.790885181877794;
        sv[22]      = 0.0;                    
        sv[23]      = 0.463478975643765;       
        sv[24]      = 0.240123148015689;
        sv[25]      = 0.0;                     
        sv[26]      = 0.649386868151536;        
        sv[27]      = 0.989963717273401;
        sv[28]      = 0.0;                     
        sv[29]      = 0.999999963501871;       
        sv[30]      = 0.926598485399264;
        sv[31]      = 0.999999963493016;    
        sv[32]      = 0.999834100785052;     
        sv[33]      = 0.999979200703676;
        sv[34]      = 0.999999963494571;      
        sv[35]      = 0.999999963612862;
        sv[36]      = 0.005470852996192;     
        sv[37]      = 0.0;                     
        sv[38]      = 0.994211562437775;        
        sv[39]      = 0.0;
        sv[40]      = 0.568856244015729;     
        sv[41]      = 0.191294664752654;      
        sv[42]      = 0.0;
        sv[43]      = 0.233014639857230;        
        sv[44]      = 0.997085416662044;     
        sv[45]      = 0.466232550741101;       
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

    OMP(parallel for private(sv_id))
    for (u_int32_t i = 0; i < num_cells_to_solve; i++) {

        if(cells_to_solve)
            sv_id = cells_to_solve[i];
        else
            sv_id = i;

        if(adpt) {
            //solve_forward_euler_cpu_adpt(sv + (sv_id * NEQ), stim_currents[i], current_t + dt, sv_id, ode_solver);
            solve_rush_larsen_cpu_adpt(sv + (sv_id * NEQ), stim_currents[i], current_t + dt, sv_id, ode_solver);
        }
        else {
            for (int j = 0; j < num_steps; ++j) {
                solve_model_ode_cpu(dt, sv + (sv_id * NEQ), stim_currents[i]);
            }
        }
    }
}

void solve_model_ode_cpu(real dt, real *sv, real stim_current) {

    const real TOLERANCE = 1e-8;
    real rY[NEQ], rDY[NEQ];

    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    // Compute 'a', 'b' coefficients alongside 'rhs'
    real a[NEQ], b[NEQ];
    RHS_RL_cpu(a, b, sv, rDY, stim_current, dt);

    // Solve variables based on its type:
    //  Non-linear = Euler
    //  Hodkin-Huxley = Rush-Larsen || Euler (if 'a' coefficient is too small)
    SOLVE_EQUATION_EULER_CPU(0);        // v        
    SOLVE_EQUATION_EULER_CPU(1);        // CaMKt    
    SOLVE_EQUATION_EULER_CPU(2);        // cass 
    SOLVE_EQUATION_EULER_CPU(3);        // nai  
    SOLVE_EQUATION_EULER_CPU(4);        // nasl 
    SOLVE_EQUATION_EULER_CPU(5);        // nass 
    SOLVE_EQUATION_EULER_CPU(6);        // ki 
    SOLVE_EQUATION_EULER_CPU(7);        // kss
    SOLVE_EQUATION_EULER_CPU(8);        // ksl
    SOLVE_EQUATION_EULER_CPU(9);        // cai
    SOLVE_EQUATION_EULER_CPU(10);       // casl
    SOLVE_EQUATION_EULER_CPU(11);       // cansr
    SOLVE_EQUATION_EULER_CPU(12);       // cajsr
    SOLVE_EQUATION_EULER_CPU(13);       // cacsr
    SOLVE_EQUATION_EULER_CPU(14);       // Jrel1
    SOLVE_EQUATION_EULER_CPU(15);       // Jrel2
    SOLVE_EQUATION_RUSH_LARSEN_CPU(16); // m
    SOLVE_EQUATION_RUSH_LARSEN_CPU(17); // hf
    SOLVE_EQUATION_RUSH_LARSEN_CPU(18); // hs
    SOLVE_EQUATION_RUSH_LARSEN_CPU(19); // j
    SOLVE_EQUATION_RUSH_LARSEN_CPU(20); // hsp
    SOLVE_EQUATION_RUSH_LARSEN_CPU(21); // jp
    SOLVE_EQUATION_RUSH_LARSEN_CPU(22); // mL
    SOLVE_EQUATION_RUSH_LARSEN_CPU(23); // hL
    SOLVE_EQUATION_RUSH_LARSEN_CPU(24); // hLp
    SOLVE_EQUATION_RUSH_LARSEN_CPU(25); // a
    SOLVE_EQUATION_RUSH_LARSEN_CPU(26); // i1
    SOLVE_EQUATION_RUSH_LARSEN_CPU(27); // i2
    SOLVE_EQUATION_RUSH_LARSEN_CPU(28); // d
    SOLVE_EQUATION_RUSH_LARSEN_CPU(29); // ff
    SOLVE_EQUATION_RUSH_LARSEN_CPU(30); // fs
    SOLVE_EQUATION_RUSH_LARSEN_CPU(31); // fcaf
    SOLVE_EQUATION_RUSH_LARSEN_CPU(32); // fcas
    SOLVE_EQUATION_RUSH_LARSEN_CPU(33); // jca
    SOLVE_EQUATION_RUSH_LARSEN_CPU(34); // ffp
    SOLVE_EQUATION_RUSH_LARSEN_CPU(35); // fcafp
    SOLVE_EQUATION_EULER_CPU(36);       // nca
    SOLVE_EQUATION_RUSH_LARSEN_CPU(37); // b
    SOLVE_EQUATION_RUSH_LARSEN_CPU(38); // g
    SOLVE_EQUATION_RUSH_LARSEN_CPU(39); // xrf
    SOLVE_EQUATION_RUSH_LARSEN_CPU(40); // xrs
    SOLVE_EQUATION_RUSH_LARSEN_CPU(41); // xs1
    SOLVE_EQUATION_RUSH_LARSEN_CPU(42); // xs2
    SOLVE_EQUATION_RUSH_LARSEN_CPU(43); // y
    SOLVE_EQUATION_RUSH_LARSEN_CPU(44); // xk1
    SOLVE_EQUATION_EULER_CPU(45);       // u
}

void solve_forward_euler_cpu_adpt(real *sv, real stim_curr, real final_time, int sv_id, struct ode_solver *solver) {

    const real _beta_safety_ = 0.8;
    int numEDO = NEQ;

    real rDY[numEDO];

    real _tolerances_[numEDO];
    real _aux_tol = 0.0;
    // initializes the variables
    solver->ode_previous_dt[sv_id] = solver->ode_dt[sv_id];

    real edos_old_aux_[numEDO];
    real edos_new_euler_[numEDO];
    real *_k1__ = (real *)malloc(sizeof(real) * numEDO);
    real *_k2__ = (real *)malloc(sizeof(real) * numEDO);
    real *_k_aux__;

    real *dt = &solver->ode_dt[sv_id];
    real *time_new = &solver->ode_time_new[sv_id];
    real *previous_dt = &solver->ode_previous_dt[sv_id];

    // Keep 'dt' inside the adaptive interval
    if(*time_new + *dt > final_time) {
        *dt = final_time - *time_new;
    }

    RHS_cpu(sv, rDY, stim_curr, *dt);
    *time_new += *dt;

    for(int i = 0; i < numEDO; i++) {
        _k1__[i] = rDY[i];
    }

    const real rel_tol = solver->rel_tol;
    const real abs_tol = solver->abs_tol;

    const real __tiny_ = pow(abs_tol, 2.0);

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
        RHS_cpu(sv, rDY, stim_curr, *dt);
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
        *dt = _beta_safety_ * (*dt) * sqrt(1.0f / greatestError);

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

void solve_rush_larsen_cpu_adpt(real *sv, real stim_curr, real final_time, int sv_id, struct ode_solver *solver) {
    
    int numEDO = NEQ;
    real rDY[numEDO];

    // initializes the variables
    solver->ode_previous_dt[sv_id] = solver->ode_dt[sv_id];

    real edos_old_aux_[numEDO];
    real edos_new_euler_[numEDO];
    real *_k1__ = (real *)malloc(sizeof(real) * numEDO);
    real *_k2__ = (real *)malloc(sizeof(real) * numEDO);
    real *a_ = (real *)malloc(sizeof(real) * numEDO);
    real *a_new = (real *)malloc(sizeof(real) * numEDO);
    real *b_ = (real *)malloc(sizeof(real) * numEDO);
    real *b_new = (real *)malloc(sizeof(real) * numEDO);
    real *_k_aux__, *_a_aux__, *_b_aux__;

    real *dt = &solver->ode_dt[sv_id];
    real *time_new = &solver->ode_time_new[sv_id];
    real *previous_dt = &solver->ode_previous_dt[sv_id];

    // Keep 'dt' inside the adaptive interval
    if(*time_new + *dt > final_time) {
        *dt = final_time - *time_new;
    }

    RHS_RL_cpu(a_, b_, sv, rDY, stim_curr, *dt);
    *time_new += *dt;

    for(int i = 0; i < numEDO; i++) {
        _k1__[i] = rDY[i];
    }

    const real rel_tol = solver->rel_tol;
    const real abs_tol = solver->abs_tol;

    const real __tiny_ = pow(abs_tol, 2.0);

    real min_dt = solver->min_dt;
    real max_dt = solver->max_dt;

    while(1) {

        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(0);        // v        
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(1);        // CaMKt    
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(2);        // cass 
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(3);        // nai  
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(4);        // nasl 
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(5);        // nass 
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(6);        // ki 
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(7);        // kss
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(8);        // ksl
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(9);        // cai
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(10);       // casl
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(11);       // cansr
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(12);       // cajsr
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(13);       // cacsr
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(14);       // Jrel1
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(15);       // Jrel2
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(16);          // m
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(17);          // hf
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(18);          // hs
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(19);          // j
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(20);          // hsp
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(21);          // jp
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(22);          // mL
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(23);          // hL
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(24);          // hLp
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(25);          // a
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(26);          // i1
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(27);          // i2
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(28);          // d
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(29);          // ff
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(30);          // fs
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(31);          // fcaf
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(32);          // fcas
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(33);          // jca
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(34);          // ffp
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(35);          // fcafp
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(36);       // nca
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(37);          // b
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(38);          // g
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(39);          // xrf
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(40);          // xrs
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(41);          // xs1
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(42);          // xs2
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(43);          // y
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_RL_CPU(44);          // xk1
        SOLVE_EQUATION_ADAPT_RUSH_LARSEN_EULER_CPU(45);       // u

        *time_new += *dt;
        RHS_RL_cpu(a_new, b_new, sv, rDY, stim_curr, *dt);
        *time_new -= *dt; // step back

        double greatestError = 0.0, auxError = 0.0;
        real as, bs, f, y_2nd_order;
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(0);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(1);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(2);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(3);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(4);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(5);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(6);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(7);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(8);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(9);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(10);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(11);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(12);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(13);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(14);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(15);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(16);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(17);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(18);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(19);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(20);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(21);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(22);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(23);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(24);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(25);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(26);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(27);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(28);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(29);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(30);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(31);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(32);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(33);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(34);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(35);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(36);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(37);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(38);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(39);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(40);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(41);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(42);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(43);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_RL_CPU(44);
        SOLVE_ERROR_ADAPT_RUSH_LARSEN_EULER_CPU(45);

        /// adapt the time step
        greatestError += __tiny_;
        *previous_dt = *dt;
        /// adapt the time step
        *dt = (*dt) * sqrt(0.5 * rel_tol / greatestError);

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

            // Swap pointers
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

void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt) {

    // Get the stimulus current from the current cell
    real calc_I_stim = stim_current;

    //State variables
    const real v = sv[0];
    const real CaMKt = sv[1];
    const real cass = sv[2];
    const real nai = sv[3];
    const real nasl = sv[4];
    const real nass = sv[5];
    const real ki = sv[6];
    const real kss = sv[7];
    const real ksl = sv[8];
    const real cai = sv[9];
    const real casl = sv[10];
    const real cansr = sv[11];
    const real cajsr = sv[12];
    const real cacsr = sv[13];
    const real Jrel1 = sv[14];
    const real Jrel2 = sv[15];
    const real m = sv[16];
    const real hf = sv[17];
    const real hs = sv[18];
    const real j = sv[19];
    const real hsp = sv[20];
    const real jp = sv[21];
    const real mL = sv[22];
    const real hL = sv[23];
    const real hLp = sv[24];
    const real a = sv[25];
    const real i1 = sv[26];
    const real i2 = sv[27];
    const real d = sv[28];
    const real ff = sv[29];
    const real fs = sv[30];
    const real fcaf = sv[31];
    const real fcas = sv[32];
    const real jca = sv[33];
    const real ffp = sv[34];
    const real fcafp = sv[35];
    const real nca = sv[36];
    const real b = sv[37];
    const real g = sv[38];
    const real xrf = sv[39];
    const real xrs = sv[40];
    const real xs1 = sv[41];
    const real xs2 = sv[42];
    const real y = sv[43];
    const real xk1 = sv[44];
    const real u = sv[45];

    #include "trovato_2020_common.inc.c"
}

void RHS_RL_cpu (real *a_, real *b_, real *sv, real *rDY_, real stim_current, real dt) {
    
    // Get the stimulus current from the current cell
    real calc_I_stim = stim_current;

    //State variables
    const real v = sv[0];
    const real CaMKt = sv[1];
    const real cass = sv[2];
    const real nai = sv[3];
    const real nasl = sv[4];
    const real nass = sv[5];
    const real ki = sv[6];
    const real kss = sv[7];
    const real ksl = sv[8];
    const real cai = sv[9];
    const real casl = sv[10];
    const real cansr = sv[11];
    const real cajsr = sv[12];
    const real cacsr = sv[13];
    const real Jrel1 = sv[14];
    const real Jrel2 = sv[15];
    const real m = sv[16];
    const real hf = sv[17];
    const real hs = sv[18];
    const real j = sv[19];
    const real hsp = sv[20];
    const real jp = sv[21];
    const real mL = sv[22];
    const real hL = sv[23];
    const real hLp = sv[24];
    const real a = sv[25];
    const real i1 = sv[26];
    const real i2 = sv[27];
    const real d = sv[28];
    const real ff = sv[29];
    const real fs = sv[30];
    const real fcaf = sv[31];
    const real fcas = sv[32];
    const real jca = sv[33];
    const real ffp = sv[34];
    const real fcafp = sv[35];
    const real nca = sv[36];
    const real b = sv[37];
    const real g = sv[38];
    const real xrf = sv[39];
    const real xrs = sv[40];
    const real xs1 = sv[41];
    const real xs2 = sv[42];
    const real y = sv[43];
    const real xk1 = sv[44];
    const real u = sv[45];

    #include "trovato_2020_RL_common.inc.c"
}