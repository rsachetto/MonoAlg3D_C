//
// Created by sachetto on 04/10/17.
//

#include "linear_system_solver.h"

uint32_t conjugate_gradient(struct grid *the_grid, int max_its, double tol, bool use_jacobi, double *error) {

    double    rTr,
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
    uint32_t number_of_iterations = 1;

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
        while( number_of_iterations < max_its ) {
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
            number_of_iterations++;
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

    return number_of_iterations;

}//end conjugateGradient() function.
