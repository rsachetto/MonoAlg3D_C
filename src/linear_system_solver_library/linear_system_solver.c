//
// Created by sachetto on 04/10/17.
//

#include "../config/linear_system_solver_config.h"
#include "../libraries_common/config_helpers.h"
#include "../libraries_common/common_data_structures.h"

bool initialized = false;
bool use_jacobi;
int max_its = 50;
double tol = 1e-16;

SOLVE_LINEAR_SYSTEM(conjugate_gradient) {

    if(!initialized) {
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(double, tol, config->config_data.config, "tolerance");
        GET_PARAMETER_BINARY_VALUE_OR_USE_DEFAULT(use_jacobi, config->config_data.config, "use_preconditioner");
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, max_its, config->config_data.config, "max_iterations");
        initialized = true;
    }


    double  rTr,
            r1Tr1,
            pTAp,
            alpha,
            beta,
            precision = tol,
            rTz,
            r1Tz1;


    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node** ac = the_grid->active_cells;

    *error = 1.0;
    *number_of_iterations = 1;

    //__________________________________________________________________________
    //Computes int_vector A*x, residue r = b - Ax, scalar rTr = r^T * r and
    //sets initial search direction p.

    rTr = 0.0;
    rTz = 0.0;

    struct element element;
    int i;

    #pragma omp parallel for private (element) reduction(+:rTr,rTz)
    for (i = 0; i < num_active_cells; i++) {

        if(CG_INFO(ac[i]) == NULL) {
            INITIALIZE_CONJUGATE_GRADIENT_INFO(ac[i]);
        }

        struct element *cell_elements = ac[i]->elements;
        ac[i]->Ax = 0.0;

        size_t max_el = sb_count(cell_elements);

        for(int el = 0; el < max_el; el++) {
            element = cell_elements[el];
            ac[i]->Ax += element.value * element.cell->v;
        }

        CG_R(ac[i]) = ac[i]->b - ac[i]->Ax;
        if(use_jacobi) {
            double value = cell_elements[0].value;
            if(value == 0.0) value = 1.0;
            CG_Z(ac[i]) = (1.0/value) * CG_R(ac[i]); // preconditioner
            rTz += CG_R(ac[i]) * CG_Z(ac[i]);
            CG_P(ac[i]) = CG_Z(ac[i]);
        }
        else {
            CG_P(ac[i]) = CG_R(ac[i]);
        }

        rTr += CG_R(ac[i]) * CG_R(ac[i]);
    }

    *error = rTr;
    //__________________________________________________________________________
    //Conjugate gradient iterations.
    if( *error >= precision ) {
        while( *number_of_iterations < max_its ) {
            //__________________________________________________________________
            // Computes Ap and pTAp. Uses Ax to store Ap.
            pTAp = 0.0;

#pragma omp parallel for private(element) reduction(+ : pTAp)
            for (i = 0; i < num_active_cells; i++) {

                ac[i]->Ax = 0.0;
                struct element *cell_elements = ac[i]->elements;

                size_t max_el = sb_count(cell_elements);
                for(int el = 0; el < max_el; el++) {
                    element = cell_elements[el];
                    ac[i]->Ax += element.value * CG_P(element.cell);
                }

                pTAp += CG_P(ac[i]) * ac[i]->Ax;
            }

            //__________________________________________________________________
            // Computes alpha.
            if(use_jacobi) {
                alpha = rTz/pTAp;
            }
            else {
                alpha = rTr/pTAp;
            }
            //__________________________________________________________________


            r1Tr1 = 0.0;
            r1Tz1 = 0.0;

            // Computes new value of solution: u = u + alpha*p.
#pragma omp parallel for reduction (+:r1Tr1,r1Tz1)
            for (i = 0; i < num_active_cells; i++) {
                ac[i]->v += alpha * CG_P(ac[i]);

                CG_R(ac[i]) -= alpha * ac[i]->Ax;

                if(use_jacobi) {
                    double value = ac[i]->elements[0].value;
                    if(value == 0.0) value = 1.0;
                    CG_Z(ac[i]) = (1.0/value) * CG_R(ac[i]);
                    r1Tz1 += CG_Z(ac[i]) * CG_R(ac[i]);
                }

                r1Tr1 += CG_R(ac[i]) * CG_R(ac[i]);
            }
            //__________________________________________________________________
            //Computes beta.
            if(use_jacobi) {
                beta = r1Tz1/rTz;
            }
            else {
                beta = r1Tr1/rTr;
            }

            *error = r1Tr1;
            *number_of_iterations = *number_of_iterations + 1;
            if( *error <= precision ) {
                break;
            }
            //__________________________________________________________________
            //Computes int_vector p1 = r1 + beta*p and uses it to upgrade p.
#pragma omp parallel for
            for (i = 0; i < num_active_cells; i++) {
                if(use_jacobi) {
                    CG_P1(ac[i]) = CG_Z(ac[i]) + beta * CG_P(ac[i]);
                }
                else {
                    CG_P1(ac[i]) = CG_R(ac[i]) + beta * CG_P(ac[i]);
                }
                CG_P(ac[i]) = CG_P1(ac[i]);
            }

            rTz = r1Tz1;
            rTr = r1Tr1;

        }

    }//end of conjugate gradient iterations.

}//end conjugateGradient() function.

// Berg's code
SOLVE_LINEAR_SYSTEM(jacobi) {


    if(!initialized) {
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(double, tol, config->config_data.config, "tolerance");
        max_its = 500;
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, max_its, config->config_data.config, "max_iterations");
        initialized = true;
    }


    double  sigma,
            precision = tol;

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node** ac = the_grid->active_cells;

    *error = 1.0;
    *number_of_iterations = 1;

    struct element element;
    int i;

    if (*error >= precision)
    {
        //__________________________________________________________________________
        //Jacobi iterations.
        while (*number_of_iterations < max_its)
        {
            #pragma omp parallel for private (element,sigma)
            for (i = 0; i < num_active_cells; i++)
            {
                if(JACOBI_INFO(ac[i]) == NULL) {
                    INITIALIZE_JACOBI_INFO(ac[i]);
                }

                struct element *cell_elements = ac[i]->elements;
                sigma = 0.0;

                size_t max_el = sb_count(cell_elements);

                // Do not take the diagonal element
                for(int el = 1; el < max_el; el++)
                {
                    element = cell_elements[el];
                    sigma += element.value * element.cell->v;
                }

                double value = cell_elements[0].value;
                JACOBI_X_AUX(ac[i]) = (1.0/value)*(ac[i]->b - sigma);
            }
            double residue = 0.0;
            double sum;
            #pragma omp parallel for private (element,sum) reduction (+:residue)
            for (i = 0; i < num_active_cells; i++)
            {
                struct element *cell_elements = ac[i]->elements;

                size_t max_el = sb_count(cell_elements);

                // Do not take the diagonal element
                sum = 0.0;
                for(int el = 0; el < max_el; el++)
                {
                    element = cell_elements[el];
                    sum += element.value * JACOBI_X_AUX(element.cell);
                }

                ac[i]->v = JACOBI_X_AUX(ac[i]);
                residue += pow(ac[i]->b - sum,2);
            }
            // The error is norm of the residue
            residue = sqrt(residue);
            *error = residue;

            *number_of_iterations = *number_of_iterations + 1;
            if( *error <= precision )
                break;
        }
    }
}

//// Berg's code
SOLVE_LINEAR_SYSTEM(biconjugate_gradient)
{


    if(!initialized) {
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(double, tol, config->config_data.config, "tolerance");
        char *preconditioner_char;
        GET_PARAMETER_VALUE_CHAR (preconditioner_char, config->config_data.config, "use_preconditioner");
        if (preconditioner_char != NULL)
        {
            use_jacobi = ((strcmp (preconditioner_char, "yes") == 0) || (strcmp (preconditioner_char, "true") == 0));
        }

        max_its = 100;
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, max_its, config->config_data.config, "max_iterations");
        initialized = true;
    }


    double  rTr,
            r1Tr1,
            pTAp,
            alpha,
            beta,
            precision = tol,
            rTz,
            r1Tz1;


    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node** ac = the_grid->active_cells;

    *error = 1.0;
    *number_of_iterations = 1;

    struct element element;
    int i;

    //__________________________________________________________________________
    // Zero all entries on the int_vector x*A
    // And initialize the second guess vector x_aux
    #pragma omp parallel for
    for (i = 0; i < num_active_cells; i++)
    {

        if(BCG_INFO(ac[i]) == NULL) {
            INITIALIZE_BICONJUGATE_GRADIENT_INFO(ac[i]);
        }

        BCG_XA(ac[i]) = 0.0;
        BCG_X_AUX(ac[i]) = ac[i]->v;
    }


    //__________________________________________________________________________
    //Computes int_vector A*x, x*A
    //xA must be fully calculated to start doing anything over the r_aux vector
    #pragma omp parallel for private (element)
    for (i = 0; i < num_active_cells; i++)
    {
        struct element *cell_elements = ac[i]->elements;
        ac[i]->Ax = 0.0;

        size_t max_el = sb_count(cell_elements);

        for(int el = 0; el < max_el; el++)
        {
            element = cell_elements[el];
            uint32_t col = element.column;
            ac[i]->Ax += element.value * element.cell->v;

            #pragma omp critical
            BCG_XA(ac[col]) += element.value * BCG_X_AUX(ac[i]);
        }
    }

    rTr = 0.0;
    rTz = 0.0;

    //__________________________________________________________________________
    //Computes residues r, r_aux
    //scalar rTr = r^T * r_aux and
    //sets initial search directions p and p_aux.
    #pragma omp parallel for private (element) reduction(+:rTr,rTz)
    for (i = 0; i < num_active_cells; i++)
    {
        struct element *cell_elements = ac[i]->elements;

        BCG_R(ac[i]) = ac[i]->b - ac[i]->Ax;
        BCG_R_AUX(ac[i]) = ac[i]->b - BCG_XA(ac[i]);

        if(use_jacobi)
        {
            double value = cell_elements[0].value;
            if(value == 0.0) value = 1.0;
            BCG_Z(ac[i]) = (1.0/value) * BCG_R(ac[i]); // preconditioner
            BCG_Z_AUX(ac[i]) = (1.0/value) * BCG_R_AUX(ac[i]);
            rTz += BCG_R_AUX(ac[i]) * BCG_Z(ac[i]);
            BCG_P(ac[i]) = BCG_Z(ac[i]);
            BCG_P_AUX(ac[i]) = BCG_Z_AUX(ac[i]);
        }
        else
        {
            BCG_P(ac[i]) = BCG_R(ac[i]);
            BCG_P_AUX(ac[i])= BCG_R_AUX(ac[i]);
        }
        rTr += BCG_R_AUX(ac[i]) * BCG_R(ac[i]);
    }

    *error = rTr;

    //__________________________________________________________________________
    //Biconjugate gradient iterations.
    if( *error >= precision )
    {
        while( *number_of_iterations < max_its )
        {
            //__________________________________________________________________
            // Computes Ap, pA and pTAp. Uses Ax to store Ap and xA to store pA
            pTAp = 0.0;

            #pragma omp parallel for
            for (i = 0; i < num_active_cells; i++)
                BCG_XA(ac[i]) = 0.0;

            #pragma omp parallel for private(element) reduction(+ : pTAp)
            for (i = 0; i < num_active_cells; i++)
            {
                ac[i]->Ax = 0.0;
                struct element *cell_elements = ac[i]->elements;

                size_t max_el = sb_count(cell_elements);
                for(int el = 0; el < max_el; el++)
                {
                    element = cell_elements[el];
                    uint32_t col = element.column;
                    ac[i]->Ax += element.value * BCG_P(element.cell);

                    #pragma omp critical
                    BCG_XA(ac[col]) += element.value * BCG_P_AUX(ac[i]);
                }

                pTAp += BCG_P_AUX(ac[i]) * ac[i]->Ax;
            }

            //__________________________________________________________________
            // Computes alpha.
            if(use_jacobi)
            {
                alpha = rTz/pTAp;
            }
            else
            {
                alpha = rTr/pTAp;
            }
            //__________________________________________________________________

            r1Tr1 = 0.0;
            r1Tz1 = 0.0;

            // Computes new value of solution: u = u + alpha*p.
            //                                 u_aux = u_aux + alpha*p_aux
            #pragma omp parallel for reduction (+:r1Tr1,r1Tz1)
            for (i = 0; i < num_active_cells; i++)
            {
                ac[i]->v += alpha * BCG_P(ac[i]);
                BCG_X_AUX(ac[i]) += alpha * BCG_P_AUX(ac[i]);

                BCG_R(ac[i]) -= alpha * ac[i]->Ax;
                BCG_R_AUX(ac[i]) -= alpha * BCG_XA(ac[i]);

                if(use_jacobi)
                {
                    double value = ac[i]->elements[0].value;
                    if(value == 0.0) value = 1.0;
                    BCG_Z(ac[i]) = (1.0/value) * BCG_R(ac[i]);
                    BCG_Z_AUX(ac[i]) = (1.0/value) * BCG_R_AUX(ac[i]);
                    r1Tz1 += BCG_Z(ac[i]) * BCG_R_AUX(ac[i]);
                }

                r1Tr1 += BCG_R(ac[i]) * BCG_R_AUX(ac[i]);
            }
            //__________________________________________________________________
            //Computes beta.
            if(use_jacobi)
            {
                beta = r1Tz1/rTz;
            }
            else
            {
                beta = r1Tr1/rTr;
            }

            *error = r1Tr1;
            *number_of_iterations = *number_of_iterations + 1;
            if( *error <= precision )
            {
                break;
            }

            //__________________________________________________________________
            //Computes int_vector p1 = r1 + beta*p and uses it to upgrade p.
            #pragma omp parallel for
            for (i = 0; i < num_active_cells; i++)
            {
                if(use_jacobi)
                {
                    BCG_P1(ac[i]) = BCG_Z(ac[i]) + beta * BCG_P(ac[i]);
                    BCG_P1_AUX(ac[i]) = BCG_Z_AUX(ac[i]) + beta * BCG_P_AUX(ac[i]);
                }
                else
                {
                    BCG_P1(ac[i]) = BCG_R(ac[i]) + beta * BCG_P(ac[i]);
                    BCG_P1_AUX(ac[i]) = BCG_R_AUX(ac[i]) + beta * BCG_P_AUX(ac[i]);
                }
                BCG_P(ac[i]) = BCG_P1(ac[i]);
                BCG_P_AUX(ac[i]) = BCG_P1_AUX(ac[i]);
            }

            rTz = r1Tz1;
            rTr = r1Tr1;

        }

    }//end of biconjugate gradient iterations.

}//end biconjugateGradient() function.