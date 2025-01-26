//
// Created by sachetto on 29/09/17.
//

#ifndef MONOALG3D_CELL_H
#define MONOALG3D_CELL_H

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../common_types/common_types.h"
#include "../../monodomain/constants.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

enum cell_type { CELL_NODE, TRANSITION_NODE };

enum transition_direction {
    FRONT, // z 0
    BACK,  // z 1
    TOP,   // y 2
    DOWN,  // y 3
    RIGHT, // x 4
    LEFT,  // x 5
    FRONT_TOP,
    FRONT_DOWN,
    FRONT_RIGHT,
    FRONT_LEFT,
    FRONT_TOP_RIGHT,
    FRONT_TOP_LEFT,
    FRONT_DOWN_RIGHT,
    FRONT_DOWN_LEFT,
    BACK_TOP,
    BACK_DOWN,
    BACK_RIGHT,
    BACK_LEFT,
    BACK_TOP_RIGHT,
    BACK_TOP_LEFT,
    BACK_DOWN_RIGHT,
    BACK_DOWN_LEFT,
    TOP_RIGHT,
    TOP_LEFT,
    DOWN_RIGHT,
    DOWN_LEFT,
    NUM_DIRECTIONS
};

#define CENTER_EQUALS(c1, c2) (c1.x == c2.x && c1.y == c2.y && c1.z == c2.z)
#define CENTER_X_EQUALS(c1, x) (c1.x == x)
#define CENTER_Y_EQUALS(c1, x) (c1.x == y)
#define CENTER_Z_EQUALS(c1, x) (c1.x == z)
#define NUM_NEIGHBOURS 6

#define VALID_SIMPLE_DIRECTION(direction) (direction >= FRONT && direction <= LEFT)
#define VALID_DIAGONAL_POINT_DIRECTION(direction) (direction >= FRONT_TOP && direction <= DOWN_LEFT)
#define STARTS_WITH_FRONT(direction) (direction >= FRONT_TOP && direction <= FRONT_DOWN_LEFT)
#define STARTS_WITH_BACK(direction) (direction >= BACK_TOP && direction <= BACK_DOWN_LEFT)
#define STARTS_WITH_TOP(direction) (direction >= TOP_RIGHT && direction <= TOP_LEFT)
#define STARTS_WITH_DOWN(direction) (direction >= DOWN_RIGHT && direction <= DOWN_LEFT)

#define TOP_IS_VISIBLE    1
#define RIGHT_IS_VISIBLE  2
#define DOWN_IS_VISIBLE   4
#define LEFT_IS_VISIBLE   8
#define BACK_IS_VISIBLE  16
#define FRONT_IS_VISIBLE 32

struct element {
#ifdef ENABLE_DDM
    enum transition_direction direction;
#endif
    real_cpu value;
    real_cpu value_ecg; //used for ecg computation
    uint32_t column; // Column of the matrix to which this element belongs.
    struct cell_node *cell;
};

struct basic_cell_data {
    enum cell_type type;
    uint8_t level; // This should be enough for the refinement levels
};

struct cell_node {

    struct basic_cell_data cell_data; // DO NOT CHANGE THIS MEMBER POSITION

    uint64_t bunch_number; // Bunch identifier

    struct point_3d center;

    void **neighbours;

    struct cell_node *previous; // Previous cell in the Hilbert curve ordering.
    struct cell_node *next;     // Next cell of in the Hilbert curve ordering.

    // Indicates position of cell on grid according to ordering provided by
    // the modified Hilbert curve.
    uint32_t grid_position;

    // Original position in the mesh file. Used for fibers.
    uint32_t original_position_in_file;

    // Variable used to storage the form of the Hilbert curve in the  bunch
    // created when this cell is refined.
    uint8_t hilbert_shape_number;

    // Fluxes used to decide if a cell should be refined or if a bunch
    // should be derefined.
    real_cpu front_flux, // Flux coming from north direction.
        back_flux,       // Flux coming from south direction.
        top_flux,        // Flux coming from east direction.
        down_flux,       // Flux coming from west direction.
        right_flux,      // Flux coming from front direction.
        left_flux;       // Flux coming from back direction.

    /* The matrix row. The elements[0] corresponds to the diagonal element of the row. */
    element_array elements;

    struct point_3d discretization;

    bool active;
    bool can_change;
    bool visited;
    uint8_t visible;

    //______________________________________________________________________________
    /* Variables used in solving the discretized system Ax = b through the conjugate gradient
   method.
   The grid discretization matrix and its resolution are directly implemented on the grid,
   which improves performance. There is no independent linear algebra package. */
    real_cpu Ax; /* Element of vector Ax = b associated to this cell. Also plays the role of Ap.*/
    real_cpu b;  /* In Ax = b, corresponds to the element in vector b associated to this cell. */

    void *linear_system_solver_extra_info;
    size_t linear_system_solver_extra_info_size;

    uint32_t sv_position;

    void *mesh_extra_info;
    size_t mesh_extra_info_size;

    // Variables used by some applications of partial differential equations.
    real_cpu v;

    struct condutivity sigma;

#ifdef ENABLE_DDM
    struct point_3d kappa;
#endif

#if defined(_OPENMP)
    omp_lock_t updating;
#endif
};

struct terminal {
    bool active;

    struct node *purkinje_cell;

    struct cell_node **tissue_cells;
};
// ----------------------------------------------------

struct transition_node {
    struct basic_cell_data cell_data; // DO NOT CHANGE THIS STRUCT POSITION

    void *single_connector;
    void *quadruple_connector1;
    void *quadruple_connector2;
    void *quadruple_connector3;
    void *quadruple_connector4;

    enum transition_direction direction;
};

struct cell_node *new_cell_node();

void init_cell_node(struct cell_node *cell_node);

void free_cell_node(struct cell_node *cell_node);

void lock_cell_node(struct cell_node *cell_node);

void unlock_cell_node(struct cell_node *cell_node);

struct transition_node *new_transition_node();

void init_transition_node(struct transition_node *transition_node);

void set_transition_node_data(struct transition_node *the_transition_node, uint16_t level, enum transition_direction direction, void *single_connector,
                              void *quadruple_connector1, void *quadruple_connector2, void *quadruple_connector3, void *quadruple_connector4);

void set_cell_node_data(struct cell_node *the_cell, struct point_3d discretization, uint64_t bunch_number, void **neighbours, void *previous, void *next,
                        uint32_t grid_position, uint8_t hilbert_shape_number, struct point_3d center);

void set_cell_flux(struct cell_node *the_cell, enum transition_direction direction);
real_cpu get_cell_maximum_flux(struct cell_node *the_cell);

void set_refined_cell_data(struct cell_node *the_cell, struct cell_node *other_cell, struct point_3d discretization, struct point_3d center,
                           uint64_t bunch_number, ui32_array free_sv_positions, ui32_array *refined_this_step);

void set_refined_transition_node_data(struct transition_node *the_node, struct cell_node *other_node, enum transition_direction direction);

void simplify_refinement(struct transition_node *transition_node);
void refine_cell(struct cell_node *cell, ui32_array free_sv_positions, ui32_array *refined_this_step);

bool cell_needs_derefinement(struct cell_node *grid_cell, real_cpu derefinement_bound);
struct cell_node *get_front_northeast_cell(struct cell_node *first_bunch_cell);
uint8_t get_father_bunch_number(struct cell_node *first_bunch_cell);
void simplify_derefinement(struct transition_node *transition_node);

void derefine_cell_bunch(struct cell_node *first_bunch_cell, ui32_array *free_sv_positions);

struct cell_node *get_cell_neighbour(struct cell_node *grid_cell, void *neighbour_grid_cell);
bool cell_has_neighbour(struct cell_node *grid_cell, void *neighbour_grid_cell);

struct cell_node *get_cell_neighbour_with_same_refinement_level(struct cell_node *grid_cell, enum transition_direction direction);

enum transition_direction get_inverse_direction(enum transition_direction direction);

int find_neighbour_index(struct cell_node *grid_cell, struct cell_node *neighbour);

uint8_t get_visibility_mask(struct cell_node *grid_cell);

#endif // MONOALG3D_CELL_H
