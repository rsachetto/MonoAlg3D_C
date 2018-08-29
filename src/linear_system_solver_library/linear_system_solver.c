//
// Created by sachetto on 04/10/17.
//

#include "../config/linear_system_solver_config.h"
#include "../libraries_common/config_helpers.h"

SOLVE_LINEAR_SYSTEM(conjugate_gradient) {

    double tol = 1e-16;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(double, tol, config->config_data.config, "tolerance");


    bool use_jacobi = true;
    char *preconditioner_char;
    GET_PARAMETER_VALUE_CHAR (preconditioner_char, config->config_data.config, "use_preconditioner");
    if (preconditioner_char != NULL) {
        use_jacobi = ((strcmp (preconditioner_char, "yes") == 0) || (strcmp (preconditioner_char, "true") == 0));
    }

    int max_its = 50;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, max_its, config->config_data.config, "max_iterations");


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
        struct element *cell_elements = ac[i]->elements;
        ac[i]->Ax = 0.0;

        size_t max_el = sb_count(cell_elements);

        for(int el = 0; el < max_el; el++) {
            element = cell_elements[el];
            ac[i]->Ax += element.value * element.cell->v;
        }

        ac[i]->r = ac[i]->b - ac[i]->Ax;
        if(use_jacobi) {
            double value = cell_elements[0].value;
            if(value == 0.0) value = 1.0;
            ac[i]->z = (1.0/value) * ac[i]->r; // preconditioner
            rTz += ac[i]->r * ac[i]->z;
            ac[i]->p = ac[i]->z;
        }
        else {
            ac[i]->p = ac[i]->r;
        }

        rTr += ac[i]->r * ac[i]->r;
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
                    ac[i]->Ax += element.value * element.cell->p;
                }

                pTAp += ac[i]->p * ac[i]->Ax;
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
                ac[i]->v += alpha * ac[i]->p;

                ac[i]->r -= alpha * ac[i]->Ax;

                if(use_jacobi) {
                    double value = ac[i]->elements[0].value;
                    if(value == 0.0) value = 1.0;
                    ac[i]->z = (1.0/value) * ac[i]->r;
                    r1Tz1 += ac[i]->z * ac[i]->r;
                }

                r1Tr1 += ac[i]->r * ac[i]->r;
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
                    ac[i]->p1 = ac[i]->z + beta * ac[i]->p;
                }
                else {
                    ac[i]->p1 = ac[i]->r + beta * ac[i]->p;
                }
                ac[i]->p = ac[i]->p1;
            }

            rTz = r1Tz1;
            rTr = r1Tr1;

        }

    }//end of conjugate gradient iterations.

}//end conjugateGradient() function.

// Berg's code
SOLVE_LINEAR_SYSTEM(jacobi) {

    double tol = 1e-08;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(double, tol, config->config_data.config, "tolerance");

    int max_its = 500;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, max_its, config->config_data.config, "max_iterations");


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
                ac[i]->x_aux = (1/value)*(ac[i]->b - sigma);
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
                    sum += element.value * element.cell->x_aux;
                }

                ac[i]->v = ac[i]->x_aux;
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

// Berg's code
SOLVE_LINEAR_SYSTEM(biconjugate_gradient)
{

    double tol = 1e-16;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(double, tol, config->config_data.config, "tolerance");


    bool use_jacobi = true;
    char *preconditioner_char;
    GET_PARAMETER_VALUE_CHAR (preconditioner_char, config->config_data.config, "use_preconditioner");
    if (preconditioner_char != NULL)
    {
        use_jacobi = ((strcmp (preconditioner_char, "yes") == 0) || (strcmp (preconditioner_char, "true") == 0));
    }

    int max_its = 100;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, max_its, config->config_data.config, "max_iterations");


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
        ac[i]->xA = 0.0;
        ac[i]->x_aux = ac[i]->v;
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
            ac[col]->xA += element.value * ac[i]->x_aux;
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

        ac[i]->r = ac[i]->b - ac[i]->Ax;
        ac[i]->r_aux = ac[i]->b - ac[i]->xA;

        if(use_jacobi)
        {
            double value = cell_elements[0].value;
            if(value == 0.0) value = 1.0;
            ac[i]->z = (1.0/value) * ac[i]->r; // preconditioner
            ac[i]->z_aux = (1.0/value) * ac[i]->r_aux;
            rTz += ac[i]->r_aux * ac[i]->z;
            ac[i]->p = ac[i]->z;
            ac[i]->p_aux = ac[i]->z_aux;
        }
        else
        {
            ac[i]->p = ac[i]->r;
            ac[i]->p_aux = ac[i]->r_aux;
        }
        rTr += ac[i]->r_aux * ac[i]->r;
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
                ac[i]->xA = 0.0;

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
                    ac[i]->Ax += element.value * element.cell->p;

                    #pragma omp critical
                    ac[col]->xA += element.value * ac[i]->p_aux;
                }

                pTAp += ac[i]->p_aux * ac[i]->Ax;
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
                ac[i]->v += alpha * ac[i]->p;
                ac[i]->x_aux += alpha * ac[i]->p_aux;

                ac[i]->r -= alpha * ac[i]->Ax;
                ac[i]->r_aux -= alpha * ac[i]->xA;

                if(use_jacobi)
                {
                    double value = ac[i]->elements[0].value;
                    if(value == 0.0) value = 1.0;
                    ac[i]->z = (1.0/value) * ac[i]->r;
                    ac[i]->z_aux = (1.0/value) * ac[i]->r_aux;
                    r1Tz1 += ac[i]->z * ac[i]->r_aux;
                }

                r1Tr1 += ac[i]->r * ac[i]->r_aux;
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
                    ac[i]->p1 = ac[i]->z + beta * ac[i]->p;
                    ac[i]->p1_aux = ac[i]->z_aux + beta * ac[i]->p_aux;
                }
                else
                {
                    ac[i]->p1 = ac[i]->r + beta * ac[i]->p;
                    ac[i]->p1_aux = ac[i]->r_aux + beta * ac[i]->p_aux;
                }
                ac[i]->p = ac[i]->p1;
                ac[i]->p_aux = ac[i]->p1_aux;
            }

            rTz = r1Tz1;
            rTr = r1Tr1;

        }

    }//end of biconjugate gradient iterations.

}//end biconjugateGradient() function.
