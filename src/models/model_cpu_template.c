#include "model_common.h"
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>


void RHS_cpu(const Real *sv, Real *rDY_, Real stim_current, Real time, Real dt);


void init_cell_model_data(struct cell_model_data* cell_model, bool get_initial_v, bool get_neq) {

    assert(cell_model);

    if(get_initial_v)
        cell_model->initial_v = INITIAL_V;
    if(get_neq)
        cell_model->number_of_ode_equations = 19;

}

void set_model_initial_conditions_cpu(Real *sv) {

    sv[0] = -85.23f;   // V;       millivolt
    //Set the rest of initial condions here
}

void solve_model_ode_cpu(Real dt, Real *sv, Real stim_current, Real time, int neq, void *extra_data)  {

    assert(sv);
    extra_data = NULL;

    Real rY[neq], rDY[neq];

    for(int i = 0; i < neq; i++)
        rY[i] = sv[i];

    RHS_cpu(rY, rDY, stim_current, time, dt);

    //THIS MODEL USES THE Rush Larsen Method TO SOLVE THE EDOS
    sv[0] = dt*rDY[0] + rY[0];

    sv[1]  = rDY[1];
    sv[2]  = rDY[2];
    sv[3]  = rDY[3];
    sv[4]  = rDY[4];
    sv[5]  = rDY[5];
    sv[6]  = rDY[6];
    sv[7]  = rDY[7];
    sv[8]  = rDY[8];
    sv[9]  = rDY[9];
    sv[10] = rDY[10];
    sv[11] = rDY[11];
    sv[12] = rDY[12];


    sv[13] = dt*rDY[13] + rY[13];
    sv[14] = dt*rDY[14] + rY[14];
    sv[15] = dt*rDY[15] + rY[15];
    sv[16] = dt*rDY[16] + rY[16];
    sv[17] = dt*rDY[17] + rY[17];
    sv[18] = dt*rDY[18] + rY[18];
}


void RHS_cpu(const Real *sv, Real *rDY_, Real stim_current, Real time, Real dt) {



}