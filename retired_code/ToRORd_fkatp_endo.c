#include "ToRORd_fkatp_endo.h"
#include <stdlib.h>
GET_CELL_MODEL_DATA(init_cell_model_data) {

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = NEQ; //for count and m
}



SET_ODE_INITIAL_CONDITIONS_CPU(set_model_initial_conditions_cpu) {

    sv[0] = -88.763800f; //v millivolt 
    sv[1] = 0.011100f; //CaMKt millimolar 
    sv[2] = 12.102500f; //nai millimolar 
    sv[3] = 12.102900f; //nass millimolar 
    sv[4] = 142.300200f; //ki millimolar 
    sv[5] = 142.300200f; //kss millimolar 
    sv[6] = 0.000082f; //cai millimolar 
    sv[7] = 0.000070f; //cass millimolar 
    sv[8] = 1.521100f; //cansr millimolar 
    sv[9] = 1.521400f; //cajsr millimolar 
    sv[10] = 0.000806f; //m dimensionless 
    sv[11] = 0.828600f; //h dimensionless 
    sv[12] = 0.828400f; //j dimensionless
    sv[13] = 0.670700f; //hp dimensionless 
    sv[14] = 0.828100f; //jp dimensionless 
    sv[15] = 0.000163f; //mL dimensionless 
    sv[16] = 0.525500f; //hL dimensionless 
    sv[17] = 0.287200f; //hLp dimensionless 
    sv[18] = 0.000951f; //a dimensionless 
    sv[19] = 0.999600f; //iF dimensionless 
    sv[20] = 0.593600f; //iS dimensionless 
    sv[21] = 0.000485f; //ap dimensionless 
    sv[22] = 0.999600f; //iFp dimensionless 
    sv[23] = 0.653800f; //iSp dimensionless 
    sv[24] = 0.000000f; //d dimensionless 
    sv[25] = 1.000000f; //ff dimensionless 
    sv[26] = 0.939000f; //fs dimensionless 
    sv[27] = 1.000000f; //fcaf dimensionless 
    sv[28] = 0.999900f; //fcas dimensionless 
    sv[29] = 1.000000f; //jca dimensionless 
    sv[30] = 1.000000f; //ffp dimensionless 
    sv[31] = 1.000000f; //fcafp dimensionless 
    sv[32] = 0.000665f; //nca_ss dimensionless 
    sv[33] = 0.001200f; //nca_i dimensionless 
    sv[34] = 0.998100f; //C3 dimensionless 
    sv[35] = 0.000851f; //C2 dimensionless 
    sv[36] = 0.000703f; //C1 dimensionless 
    sv[37] = 0.000376f; //O dimensionless 
    sv[38] = 0.000013f; //I dimensionless 
    sv[39] = 0.248000f; //xs1 dimensionless 
    sv[40] = 0.000177f; //xs2 dimensionless 
    sv[41] = 0.000000f; //Jrel_np millimolar_per_millisecond 
    sv[42] = 0.000000f; //Jrel_p millimolar_per_millisecond 
}

real *ode_dt, *ode_previous_dt, *ode_time_new;
static bool first_call = true;
void addt2(real *sv, real final_time, real *dt, real *previous_dt, real *time_new, real Istim, real abstol, real reltol, real min_step, real max_step);

SOLVE_MODEL_ODES_CPU(solve_model_odes_cpu) {

    uint32_t sv_id;

    if(first_call) {
        ode_dt          = (real*)malloc(num_cells_to_solve*sizeof(real));

        for(int i = 0; i < num_cells_to_solve; i++) {
            ode_dt[i] = dt;
        }

        ode_previous_dt = (real*)calloc(num_cells_to_solve, sizeof(real));
        ode_time_new    = (real*)calloc(num_cells_to_solve, sizeof(real));
        first_call = false;
    }

    #pragma omp parallel for private(sv_id)
    for (u_int32_t i = 0; i < num_cells_to_solve; i++) {

        if(cells_to_solve)
            sv_id = cells_to_solve[i];
        else
            sv_id = i;

        addt2(sv + (sv_id * NEQ), current_t + dt, &ode_dt[sv_id], &ode_previous_dt[sv_id], &ode_time_new[sv_id], stim_currents[i], 1e-10, 1e-10, 0.0001, 0.01);

    }
}

void addt2(real *sv, real final_time, real *dt, real *previous_dt, real *time_new, real Istim, real abstol, real reltol, real min_step, real max_step) {

    const real _beta_safety_ = 0.8;
    int numEDO = NEQ;

    real rDY[numEDO];

    real _tolerances_[numEDO];
    real _aux_tol = 0.0;
    //initializes the variables
    *previous_dt = *dt;

    real edos_old_aux_[numEDO];
    real edos_new_euler_[numEDO];
    real *_k1__ = (real*) malloc(sizeof(real)*numEDO);
    real *_k2__ = (real*) malloc(sizeof(real)*numEDO);
    real *_k_aux__;

    if(*time_new + *dt > final_time) {
        *dt = final_time - *time_new;
    }

    //lado direito da primeira vez
    RHS_cpu(sv, rDY, Istim, *dt);
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
        RHS_cpu(sv, rDY, Istim, *dt);
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


        //TALVEZ LIMITAR O DT AQUI MESMO PQ NAO FAZ SENTIDO IR MENOR DO QUE NOS SABEMOS QUE DA CERTO
        if (*time_new + *dt > final_time) {
            *dt = final_time - *time_new;
        }

        //it doesn't accept the solution
        if ( count < count_limit  && (greatestError >= 1.0f)) {
            //restore the old values to do it again
            for(int i = 0;  i < numEDO; i++) {
                sv[i] = edos_old_aux_[i];
            }

            //printf("COUNT: %d of %d using dt = %lf\n", count, count_limit, _ode->__dt);
            count++;
            //throw the results away and compute again
        } else{//it accepts the solutions


            // printf("Time %lf of %lf \n", _ode->time_new, final_time);
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
            //change vectors k1 e k2 , para que k2 seja aproveitado como k1 na proxima iteração
            _k_aux__ = _k2__;
            _k2__	= _k1__;
            _k1__	= _k_aux__;

            //it steps the method ahead, with euler solution
            for(int i = 0; i < numEDO; i++){
                sv[i] = edos_new_euler_[i];
            }

            //verifica se o incremento para a próxima iteração ultrapassa o tempo de salvar, q neste caso é o tempo final
            if(*time_new + *previous_dt >= final_time){
                //se são iguais, ja foi calculada a iteração no ultimo passo de tempo e deve-se para o laço
                //nao usar igualdade - usar esta conta, pode-se mudar a tolerância
                if((fabs(final_time - *time_new) < 1.0e-5) ){
                    break;
                    //se o passo de tempo + time_new ultrapassam o momento final, ajusta-se dt para que a próxima iteração do método ocorra exatamente no passo final.
                }else if(*time_new < final_time){
                    *dt = *previous_dt = final_time - *time_new;
                    *time_new += *previous_dt;
                    break;

                }else{
                    printf("PAU - nao eh bom chegar aqui time_new %.20lf final_time %.20lf diff %e \n", *time_new , final_time, fabs(final_time - *time_new) );
                    break;
                }
            }else{//incremento normal - nao esta na "borda" do tempo
                *time_new += *previous_dt;
            }

        }
    }//fim-while

    free(_k1__);
    free(_k2__);
}

#define IFNUMBER_1(name)if((celltype==1.000000000000000e+00f)) { (name) = (cmdnmax_b*1.300000000000000e+00f);    }  else{ (name) = cmdnmax_b;    }
#define IFNUMBER_2(name)if((v_old_>=(-4.000000000000000e+01f))) { (name) = 0.000000000000000e+00f;    }  else{ (name) = (5.700000000000000e-02f*expf(((-(v_old_+8.000000000000000e+01f))/6.800000000000000e+00f)));    }
#define IFNUMBER_3(name)if((v_old_>=(-4.000000000000000e+01f))) { (name) = (7.700000000000000e-01f/(1.300000000000000e-01f*(1.000000000000000e+00f+expf(((-(v_old_+1.066000000000000e+01f))/1.110000000000000e+01f)))));    }  else{ (name) = ((2.700000000000000e+00f*expf((7.900000000000000e-02f*v_old_)))+(3.150000000000000e+00f*expf((3.485000000000000e-01f*v_old_))));    }
#define IFNUMBER_4(name)if((v_old_>=(-4.000000000000000e+01f))) { (name) = 0.000000000000000e+00f;    }  else{ (name) = (((((-2.542840000000000e+00f)*expf((2.444000000000000e-01f*v_old_)))-(6.948000000000000e+00f*expf(((-4.391000000000000e-02f)*v_old_))))*(v_old_+3.778000000000000e+01f))/(1.000000000000000e+00f+expf((3.110000000000000e-01f*(v_old_+7.923000000000000e+01f)))));    }
#define IFNUMBER_5(name)if((v_old_>=(-4.000000000000000e+01f))) { (name) = ((6.000000000000000e-01f*expf((5.700000000000000e-02f*v_old_)))/(1.000000000000000e+00f+expf(((-1.000000000000000e-01f)*(v_old_+3.200000000000000e+01f)))));    }  else{ (name) = ((2.424000000000000e-02f*expf(((-1.052000000000000e-02f)*v_old_)))/(1.000000000000000e+00f+expf(((-1.378000000000000e-01f)*(v_old_+4.014000000000000e+01f)))));    }
#define IFNUMBER_6(name)if((celltype==1.000000000000000e+00f)) { (name) = (GNaL_b*6.000000000000000e-01f);    }  else{ (name) = GNaL_b;    }
#define IFNUMBER_7(name)if((celltype==1.000000000000000e+00f)) { (name) = (1.000000000000000e+00f-(9.500000000000000e-01f/(1.000000000000000e+00f+expf(((v_old_+EKshift+7.000000000000000e+01f)/5.000000000000000e+00f)))));    }  else{ (name) = 1.000000000000000e+00f;    }
#define IFNUMBER_8(name)if((celltype==1.000000000000000e+00f)) { (name) = (Gto_b*2.000000000000000e+00f);    }  else if((celltype==2.000000000000000e+00f)){ (name) = (Gto_b*2.000000000000000e+00f);    } else{ (name) = Gto_b;    }
#define IFNUMBER_9(name)if((v_old_>=3.149780000000000e+01f)) { (name) = 1.000000000000000e+00f;    }  else{ (name) = (1.076300000000000e+00f*expf(((-1.007000000000000e+00f)*expf(((-8.290000000000000e-02f)*v_old_)))));    }
#define IFNUMBER_10(name)if((celltype==1.000000000000000e+00f)) { (name) = (PCa_b*1.200000000000000e+00f);    }  else if((celltype==2.000000000000000e+00f)){ (name) = (PCa_b*2.000000000000000e+00f);    } else{ (name) = PCa_b;    }
#define IFNUMBER_11(name)if((celltype==1.000000000000000e+00f)) { (name) = (GKr_b*1.300000000000000e+00f);    }  else if((celltype==2.000000000000000e+00f)){ (name) = (GKr_b*8.000000000000000e-01f);    } else{ (name) = GKr_b;    }
#define IFNUMBER_12(name)if((celltype==1.000000000000000e+00f)) { (name) = (GKs_b*1.400000000000000e+00f);    }  else{ (name) = GKs_b;    }
#define IFNUMBER_13(name)if((celltype==1.000000000000000e+00f)) { (name) = (GK1_b*1.200000000000000e+00f);    }  else if((celltype==2.000000000000000e+00f)){ (name) = (GK1_b*1.300000000000000e+00f);    } else{ (name) = GK1_b;    }
#define IFNUMBER_14(name)if((celltype==1.000000000000000e+00f)) { (name) = (Gncx_b*1.100000000000000e+00f);    }  else if((celltype==2.000000000000000e+00f)){ (name) = (Gncx_b*1.400000000000000e+00f);    } else{ (name) = Gncx_b;    }
#define IFNUMBER_15(name)if((celltype==1.000000000000000e+00f)) { (name) = (Pnak_b*9.000000000000000e-01f);    }  else if((celltype==2.000000000000000e+00f)){ (name) = (Pnak_b*7.000000000000000e-01f);    } else{ (name) = Pnak_b;    }
#define IFNUMBER_16(name)if((celltype==1.000000000000000e+00f)) { (name) = (GKb_b*6.000000000000000e-01f);    }  else{ (name) = GKb_b;    }
#define IFNUMBER_17(name)if((celltype==2.000000000000000e+00f)) { (name) = (calc_Jrel_inf_b*1.700000000000000e+00f);    }  else{ (name) = calc_Jrel_inf_b;    }
#define IFNUMBER_18(name)if((calc_tau_rel_b<1.000000000000000e-03f)) { (name) = 1.000000000000000e-03f;    }  else{ (name) = calc_tau_rel_b;    }
#define IFNUMBER_19(name)if((celltype==2.000000000000000e+00f)) { (name) = (calc_Jrel_infp_b*1.700000000000000e+00f);    }  else{ (name) = calc_Jrel_infp_b;    }
#define IFNUMBER_20(name)if((calc_tau_relp_b<1.000000000000000e-03f)) { (name) = 1.000000000000000e-03f;    }  else{ (name) = calc_tau_relp_b;    }
#define IFNUMBER_21(name)if((celltype==1.000000000000000e+00f)) { (name) = 1.300000000000000e+00f;    }  else{ (name) = 1.000000000000000e+00f;    }

void RHS_cpu(const real *sv, real *rDY_, real stim_current, real dt) {

    //State variables
    const real v_old_ = sv[0];
    const real CaMKt_old_ = sv[1];
    const real nai_old_ = sv[2];
    const real nass_old_ = sv[3];
    const real ki_old_ = sv[4];
    const real kss_old_ = sv[5];
    const real cai_old_ = sv[6];
    const real cass_old_ = sv[7];
    const real cansr_old_ = sv[8];
    const real cajsr_old_ = sv[9];
    const real m_old_ = sv[10];
    const real h_old_ = sv[11];
    const real j_old_ = sv[12];
    const real hp_old_ = sv[13];
    const real jp_old_ = sv[14];
    const real mL_old_ = sv[15];
    const real hL_old_ = sv[16];
    const real hLp_old_ = sv[17];
    const real a_old_ = sv[18];
    const real iF_old_ = sv[19];
    const real iS_old_ = sv[20];
    const real ap_old_ = sv[21];
    const real iFp_old_ = sv[22];
    const real iSp_old_ = sv[23];
    const real d_old_ = sv[24];
    const real ff_old_ = sv[25];
    const real fs_old_ = sv[26];
    const real fcaf_old_ = sv[27];
    const real fcas_old_ = sv[28];
    const real jca_old_ = sv[29];
    const real ffp_old_ = sv[30];
    const real fcafp_old_ = sv[31];
    const real nca_ss_old_ = sv[32];
    const real nca_i_old_ = sv[33];
    const real C3_old_ = sv[34];
    const real C2_old_ = sv[35];
    const real C1_old_ = sv[36];
    const real O_old_ = sv[37];
    const real I_old_ = sv[38];
    const real xs1_old_ = sv[39];
    const real xs2_old_ = sv[40];
    const real Jrel_np_old_ = sv[41];
    const real Jrel_p_old_ = sv[42];

    #include "ToROrd_common.inc.c"
}