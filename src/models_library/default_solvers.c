void solve_model_ode_cpu(real dt, real *sv, real stim_current)  {

    real rY[NEQ], rDY[NEQ];

    for(int i = 0; i < NEQ; i++)
        rY[i] = sv[i];

    RHS_cpu(rY, rDY, stim_current, dt);

    for(int i = 0; i < NEQ; i++)
        sv[i] = dt*rDY[i] + rY[i];
}

void solve_forward_euler_cpu_adpt(real *sv, real stim_curr, real final_time, int sv_id) {

    const real _beta_safety_ = 0.8;
    int numEDO = NEQ;

    real rDY[numEDO];

    real _tolerances_[numEDO];
    real _aux_tol = 0.0;
    //initializes the variables
    ode_previous_dt[sv_id] = ode_dt[sv_id];

    real edos_old_aux_[numEDO];
    real edos_new_euler_[numEDO];
    real *_k1__ = (real*) malloc(sizeof(real)*numEDO);
    real *_k2__ = (real*) malloc(sizeof(real)*numEDO);
    real *_k_aux__;

    real *dt = &ode_dt[sv_id];
    real *time_new = &ode_time_new[sv_id];
    real *previous_dt = &ode_previous_dt[sv_id];

    if(*time_new + *dt > final_time) {
       *dt = final_time - *time_new;
    }

    RHS_cpu(sv, rDY, stim_curr, *dt);
    *time_new += *dt;

    for(int i = 0; i < numEDO; i++){
        _k1__[i] = rDY[i];
    }

    const double __tiny_ = pow(abstol, 2.0);

    int count = 0;

    int count_limit = (final_time - *time_new)/min_step;

    int aux_count_limit = count_limit+2000000;

    if(aux_count_limit > 0) {
        count_limit = aux_count_limit;
    }

    while(1) {

        for(int i = 0; i < numEDO; i++) {
            //stores the old variables in a vector
            edos_old_aux_[i] = sv[i];
            //computes euler method
            edos_new_euler_[i] = _k1__[i] * *dt + edos_old_aux_[i];
            //steps ahead to compute the rk2 method
            sv[i] = edos_new_euler_[i];
        }

        *time_new += *dt;
        RHS_cpu(sv, rDY, stim_curr, *dt);
        *time_new -= *dt;//step back

        double greatestError = 0.0, auxError = 0.0;
        for(int i = 0; i < numEDO; i++) {
            //stores the new evaluation
            _k2__[i] = rDY[i];
            _aux_tol = fabs(edos_new_euler_[i])*reltol;
            _tolerances_[i] = (abstol > _aux_tol )?abstol:_aux_tol;
            //finds the greatest error between  the steps
            auxError = fabs(( (*dt/2.0)*(_k1__[i] - _k2__[i])) / _tolerances_[i]);

            greatestError = (auxError > greatestError) ? auxError : greatestError;
        }
        ///adapt the time step
        greatestError += __tiny_;
        *previous_dt = *dt;
        ///adapt the time step
        *dt = _beta_safety_ * (*dt) * sqrt(1.0f/greatestError);

        if (*time_new + *dt > final_time) {
            *dt = final_time - *time_new;
        }

        //it doesn't accept the solution
        if ( count < count_limit  && (greatestError >= 1.0f)) {
            //restore the old values to do it again
            for(int i = 0;  i < numEDO; i++) {
                sv[i] = edos_old_aux_[i];
            }

            count++;
            //throw the results away and compute again
        } else{//it accepts the solutions


            if(greatestError >=1.0) {
                printf("Accepting solution with error > %lf \n", greatestError);
            }

            //printf("%e %e\n", _ode->time_new, edos_new_euler_[0]);
            if (*dt < min_step) {
                *dt = min_step;
            }

            else if (*dt > max_step && max_step != 0) {
                *dt = max_step;
            }

            if (*time_new + *dt > final_time) {
                *dt = final_time - *time_new;
            }

            _k_aux__ = _k2__;
            _k2__	= _k1__;
            _k1__	= _k_aux__;

            //it steps the method ahead, with euler solution
            for(int i = 0; i < numEDO; i++){
                sv[i] = edos_new_euler_[i];
            }

            if(*time_new + *previous_dt >= final_time){
                if((fabs(final_time - *time_new) < 1.0e-5) ){
                    break;
                }else if(*time_new < final_time){
                    *dt = *previous_dt = final_time - *time_new;
                    *time_new += *previous_dt;
                    break;

                }else{
                    printf("Error: time_new %.20lf final_time %.20lf diff %e \n", *time_new , final_time, fabs(final_time - *time_new) );
                    break;
                }
            }else{
                *time_new += *previous_dt;
            }

        }
    }

    free(_k1__);
    free(_k2__);
}

