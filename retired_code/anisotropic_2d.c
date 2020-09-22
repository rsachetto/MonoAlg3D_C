//
// Created by sachetto on 10/06/2020.
//

static void fill_elements_aniso(struct cell_node *grid_cell, struct cell_node *neighbours[26]) {

    struct element *elements = grid_cell->elements;

    real_cpu dx = grid_cell->discretization.x;
    real_cpu dy = grid_cell->discretization.y;

    real_cpu dx_div_dy = dx / dy;
    real_cpu dy_div_dx = dy / dx;

    bool has_all_n_for_ip1_j_flux1 = neighbours[TOP_RIGHT] && neighbours[RIGHT] && neighbours[TOP];
    bool has_all_n_for_ip1_j_flux2 = neighbours[RIGHT] && neighbours[DOWN_RIGHT] && neighbours[DOWN];

    bool has_ip1_j_flux1_sx = neighbours[RIGHT] || (neighbours[TOP] && neighbours[TOP_RIGHT]);
    bool has_ip1_j_flux1_sxy = neighbours[TOP] || (neighbours[RIGHT] && neighbours[TOP_RIGHT]);

    bool has_ip1_j_flux1 = has_ip1_j_flux1_sx || has_ip1_j_flux1_sxy;

    bool has_ip1_j_flux2_sx = neighbours[RIGHT] || (neighbours[DOWN] && neighbours[DOWN_RIGHT]);
    bool has_ip1_j_flux2_sxy = neighbours[DOWN] || (neighbours[RIGHT] && neighbours[DOWN_RIGHT]);
    bool has_ip1_j_flux2 = has_ip1_j_flux2_sx || has_ip1_j_flux2_sxy;

    bool has_all_n_for_im1_j_flux1 = neighbours[TOP] && neighbours[TOP_LEFT] && neighbours[LEFT];
    bool has_all_n_for_im1_j_flux2 = neighbours[DOWN] && neighbours[LEFT] && neighbours[DOWN_LEFT];

    bool has_im1_j_flux1_sx = neighbours[LEFT] || (neighbours[TOP] && neighbours[TOP_LEFT]);
    bool has_im1_j_flux1_sxy = neighbours[TOP] || (neighbours[LEFT] && neighbours[TOP_LEFT]);
    bool has_im1_j_flux1 = has_im1_j_flux1_sx || has_im1_j_flux1_sxy;

    bool has_im1_j_flux2_sx = neighbours[LEFT] || (neighbours[DOWN] && neighbours[DOWN_LEFT]);
    bool has_im1_j_flux2_sxy = neighbours[DOWN] || (neighbours[LEFT] && neighbours[DOWN_LEFT]);
    bool has_im1_j_flux2 = has_im1_j_flux2_sx || has_im1_j_flux2_sxy;

    bool has_i_jp1_flux1 = has_ip1_j_flux1;
    bool has_i_jp1_flux2 = has_im1_j_flux1;

    bool has_i_jm1_flux1 = has_im1_j_flux2;
    bool has_i_jm1_flux2 = has_ip1_j_flux2;

    // Jx_i+1/2,j
    real_cpu sigma_x_tr = 0.0;
    real_cpu sigma_xy_tr = 0.0;

    real_cpu sigma_x_dr = 0.0;
    real_cpu sigma_xy_dr = 0.0;

    // Jx_i-1/2,j
    real_cpu sigma_x_tl = 0.0;
    real_cpu sigma_xy_tl = 0.0;

    real_cpu sigma_x_dl = 0.0;
    real_cpu sigma_xy_dl = 0.0;

    // Jy_i, j+1/2
    real_cpu sigma_y_tr = 0.0;
    real_cpu sigma_y_tl = 0.0;

    // Jy_i, j-1/2
    real_cpu sigma_y_dl = 0.0;
    real_cpu sigma_y_dr = 0.0;

    //Main diagonal
    // Jx_i+1/2,j
    if(has_all_n_for_ip1_j_flux1) {
        if(neighbours[TOP_RIGHT]->sigma.x > 0 && neighbours[RIGHT]->sigma.x > 0 && neighbours[TOP]->sigma.x > 0 && grid_cell->sigma.x > 0) {
            sigma_x_tr = 4.0 / (1.0 / neighbours[TOP_RIGHT]->sigma.x + 1.0 / neighbours[RIGHT]->sigma.x +
                                1.0 / neighbours[TOP]->sigma.x + 1.0 / grid_cell->sigma.x);
        }

        if(neighbours[TOP_RIGHT]->sigma.xy > 0 && neighbours[RIGHT]->sigma.xy > 0 && neighbours[TOP]->sigma.xy > 0 && grid_cell->sigma.xy > 0) {
            sigma_xy_tr = 4.0 / (1.0 / neighbours[TOP_RIGHT]->sigma.xy + 1.0 / neighbours[RIGHT]->sigma.xy +
                                 1.0 / neighbours[TOP]->sigma.xy + 1.0 / grid_cell->sigma.xy);
        }

        if(neighbours[TOP_RIGHT]->sigma.y > 0 && neighbours[RIGHT]->sigma.y > 0 && neighbours[TOP]->sigma.y > 0 && grid_cell->sigma.y > 0) {
            sigma_y_tr = 4.0 / (1.0 / neighbours[TOP_RIGHT]->sigma.y + 1.0 / neighbours[RIGHT]->sigma.y +
                                1.0 / neighbours[TOP]->sigma.y + 1.0 / grid_cell->sigma.y);
        }

    }

    if(has_all_n_for_ip1_j_flux2) {
        if(neighbours[DOWN_RIGHT]->sigma.x > 0 && neighbours[RIGHT]->sigma.x > 0 && neighbours[DOWN]->sigma.x > 0 && grid_cell->sigma.x > 0) {
            sigma_x_dr = 4.0 / (1.0 / neighbours[DOWN_RIGHT]->sigma.x + 1.0 / neighbours[RIGHT]->sigma.x +
                                1.0 / neighbours[DOWN]->sigma.x + 1.0 / grid_cell->sigma.x);
        }

        if(neighbours[DOWN_RIGHT]->sigma.xy > 0 && neighbours[RIGHT]->sigma.xy > 0 && neighbours[DOWN]->sigma.xy > 0 &&
           grid_cell->sigma.xy > 0) {
            sigma_xy_dr = 4.0 / (1.0 / neighbours[DOWN_RIGHT]->sigma.xy + 1.0 / neighbours[RIGHT]->sigma.xy +
                                 1.0 / neighbours[DOWN]->sigma.xy + 1.0 / grid_cell->sigma.xy);
        }

        if(neighbours[DOWN_RIGHT]->sigma.y > 0 && neighbours[RIGHT]->sigma.y > 0 && neighbours[DOWN]->sigma.y > 0 && grid_cell->sigma.y > 0) {
            sigma_y_dr = 4.0 / (1.0 / neighbours[DOWN_RIGHT]->sigma.y + 1.0 / neighbours[RIGHT]->sigma.y +
                                1.0 / neighbours[DOWN]->sigma.y + 1.0 / grid_cell->sigma.y);
        }

    }

    // Jx_i-1/2,j
    if(has_all_n_for_im1_j_flux1) {

        if(neighbours[TOP]->sigma.x > 0 && neighbours[TOP_LEFT]->sigma.x > 0 && neighbours[LEFT]->sigma.x > 0 && grid_cell->sigma.x > 0) {
            sigma_x_tl = 4.0 / (1.0 / neighbours[TOP]->sigma.x + 1.0 / neighbours[TOP_LEFT]->sigma.x +
                                1.0 / neighbours[LEFT]->sigma.x + 1.0 / grid_cell->sigma.x);
        }

        if(neighbours[TOP]->sigma.xy > 0 && neighbours[TOP_LEFT]->sigma.xy > 0 && neighbours[LEFT]->sigma.xy > 0 && grid_cell->sigma.xy > 0) {
            sigma_xy_tl = 4.0 / (1.0 / neighbours[TOP]->sigma.xy + 1.0 / neighbours[TOP_LEFT]->sigma.xy +
                                 1.0 / neighbours[LEFT]->sigma.xy + 1.0 / grid_cell->sigma.xy);
        }
        if(neighbours[TOP]->sigma.y > 0 && neighbours[TOP_LEFT]->sigma.y > 0 && neighbours[LEFT]->sigma.y > 0 && grid_cell->sigma.y > 0) {
            sigma_y_tl = 4.0 / (1.0 / neighbours[TOP]->sigma.y + 1.0 / neighbours[TOP_LEFT]->sigma.y +
                                1.0 / neighbours[LEFT]->sigma.y + 1.0 / grid_cell->sigma.y);
        }
    }

    if(has_all_n_for_im1_j_flux2) {

        if(neighbours[LEFT]->sigma.x > 0 && neighbours[DOWN]->sigma.x > 0 && neighbours[DOWN_LEFT]->sigma.x > 0 && grid_cell->sigma.x > 0) {
            sigma_x_dl = 4.0 / (1.0 / neighbours[LEFT]->sigma.x + 1.0 / neighbours[DOWN]->sigma.x +
                                1.0 / neighbours[DOWN_LEFT]->sigma.x + 1.0 / grid_cell->sigma.x);
        }
        if(neighbours[LEFT]->sigma.xy > 0 && neighbours[DOWN]->sigma.xy > 0 && neighbours[DOWN_LEFT]->sigma.xy > 0 && grid_cell->sigma.xy > 0) {
            sigma_xy_dl = 4.0 / (1.0 / neighbours[LEFT]->sigma.xy + 1.0 / neighbours[DOWN]->sigma.xy +
                                 1.0 / neighbours[DOWN_LEFT]->sigma.xy + 1.0 / grid_cell->sigma.xy);
        }

        if(neighbours[LEFT]->sigma.y > 0 && neighbours[DOWN]->sigma.y > 0 && neighbours[DOWN_LEFT]->sigma.y > 0 && grid_cell->sigma.y > 0) {

            sigma_y_dl = 4.0 / (1.0 / neighbours[LEFT]->sigma.y + 1.0 / neighbours[DOWN]->sigma.y +
                                1.0 / neighbours[DOWN_LEFT]->sigma.y + 1.0 / grid_cell->sigma.y);
        }

    }

    real_cpu partial_x;
    real_cpu partial_y1;
    real_cpu partial_y2;

    // Main diagonal

    partial_x = 0.0;
    partial_y1 = 0.0;
    partial_y2 = 0.0;

    if(has_all_n_for_ip1_j_flux1) {
        // J_x_i+1/2,j f1
        partial_x -= sigma_x_tr / 2.0;
        partial_x -= dx_div_dy * sigma_xy_tr / 2.0;

        // J_y_i,j+1/2 f1
        partial_y1 -= dy_div_dx * sigma_xy_tr / 2.0;
        partial_y1 -= sigma_y_tr / 2.0;

    } else {
        if(has_ip1_j_flux1_sx) {
            if(neighbours[RIGHT]) {
                // J_x_i+1/2,j f1
                real_cpu local_sigma_x_tr = 2.0 / (1.0 / neighbours[RIGHT]->sigma.x + 1.0 / grid_cell->sigma.x);
                partial_x -= local_sigma_x_tr;

                // J_y_i,j+1/2 f1
                if(neighbours[RIGHT]->sigma.xy > 0 && grid_cell->sigma.xy > 0) {
                    real_cpu local_sigma_xy_tr = 2.0 / (1.0 / neighbours[RIGHT]->sigma.xy + 1.0 / grid_cell->sigma.xy);
                    partial_y1 -= dy_div_dx * local_sigma_xy_tr;
                }

            }
        }
        if(has_ip1_j_flux1_sxy) {

            if(neighbours[TOP]) {
                // J_x_i+1/2,j f1
                if(neighbours[TOP]->sigma.xy > 0 && grid_cell->sigma.xy > 0) {
                    real_cpu local_sigma_xy_tr = 2.0 / (1.0 / neighbours[TOP]->sigma.xy + 1.0 / grid_cell->sigma.xy);
                    partial_x -= dx_div_dy * local_sigma_xy_tr;
                }

                // J_y_i,j+1/2 f1
                real_cpu local_sigma_y_tr = 2.0 / (1.0 / neighbours[TOP]->sigma.y + 1.0 / grid_cell->sigma.y);
                partial_y1 -= local_sigma_y_tr;
            }
        }
    }

    if(has_all_n_for_ip1_j_flux2) {
        // J_x_i+1/2,j f2
        partial_x -= sigma_x_dr / 2.0;
        partial_x += dx_div_dy * sigma_xy_dr / 2.0;

        // J_y_i,j-1/2 f2
        partial_y2 += dy_div_dx * sigma_xy_dr / 2.0;
        partial_y2 -= sigma_y_dr / 2.0;

    } else {
        if(has_ip1_j_flux2_sx) {

            if(neighbours[RIGHT]) {
                // J_x_i+1/2,j f2
                real_cpu local_sigma_x_dr = 2.0 / (1.0 / neighbours[RIGHT]->sigma.x + 1.0 / grid_cell->sigma.x);
                partial_x -= local_sigma_x_dr;

                if(neighbours[RIGHT]->sigma.xy > 0 && grid_cell->sigma.xy > 0) {
                    // J_y_i,j-1/2 f2
                    real_cpu local_sigma_xy_dr = 2.0 / (1.0 / neighbours[RIGHT]->sigma.xy + 1.0 / grid_cell->sigma.xy);
                    partial_y2 += dy_div_dx * local_sigma_xy_dr;
                }
            }
        }

        if(has_ip1_j_flux2_sxy) {
            // J_x_i+1/2,j f2
            if(neighbours[DOWN]) {
                if(neighbours[DOWN]->sigma.xy > 0 && grid_cell->sigma.xy > 0) {
                    real_cpu local_sigma_xy_dr = 2.0 / (1.0 / neighbours[DOWN]->sigma.xy + 1.0 / grid_cell->sigma.xy);
                    partial_x += dx_div_dy * local_sigma_xy_dr;
                }

                // J_y_i,j-1/2 f2
                real_cpu local_sigma_y_dr = 2.0 / (1.0 / neighbours[DOWN]->sigma.y + 1.0 / grid_cell->sigma.y);
                partial_y2 -= local_sigma_y_dr;
            }
        }
    }

    if(has_ip1_j_flux1 && has_ip1_j_flux2) {
        partial_x /= 2.0;
    }

    elements[0].value += -partial_x;

    partial_x = 0;
    if(has_all_n_for_im1_j_flux1) {
        // J_x_i-1/2,j f1
        partial_x -= sigma_x_tl / 2.0;
        partial_x += dx_div_dy * sigma_xy_tl / 2.0;

        // J_y_i,j+1/2 f2
        partial_y1 += dy_div_dx * sigma_xy_tl / 2.0;
        partial_y1 -= sigma_y_tl / 2.0;

    } else {
        if(has_im1_j_flux1_sx) {
            if(neighbours[LEFT]) {
                // J_x_i-1/2,j f1
                real_cpu local_sigma_x_tl = 2.0 / (1.0 / neighbours[LEFT]->sigma.x + 1.0 / grid_cell->sigma.x);
                partial_x -= local_sigma_x_tl;

                if(neighbours[LEFT]->sigma.xy > 0 && grid_cell->sigma.xy > 0) {
                    // J_y_i,j+1/2 f2
                    real_cpu local_sigma_xy_tl = 2.0 / (1.0 / neighbours[LEFT]->sigma.xy + 1.0 / grid_cell->sigma.xy);
                    partial_y1 += dy_div_dx * local_sigma_xy_tl;
                }
            }
        }

        if(has_im1_j_flux1_sxy) {

            if(neighbours[TOP]) {
                if(neighbours[TOP]->sigma.xy > 0 && grid_cell->sigma.xy > 0) {
                    // J_x_i-1/2,j f1
                    real_cpu local_sigma_xy_tl = 2.0 / (1.0 / neighbours[TOP]->sigma.xy + 1.0 / grid_cell->sigma.xy);
                    partial_x += dx_div_dy * local_sigma_xy_tl;
                }

                // J_y_i,j+1/2 f2
                real_cpu local_sigma_y_tl = 2.0 / (1.0 / neighbours[TOP]->sigma.y + 1.0 / grid_cell->sigma.y);
                partial_y1 -= local_sigma_y_tl;
            }
        }
    }

    if(has_i_jp1_flux1 && has_i_jp1_flux2) {
        partial_y1 /= 2.0;
    }

    elements[0].value += -partial_y1;

    if(has_all_n_for_im1_j_flux2) {
        // J_x_i-1/2,j f2
        partial_x -= sigma_x_dl / 2.0;
        partial_x -= dx_div_dy * sigma_xy_dl / 2.0;

        // J_y_i,j-1/2 f1
        partial_y2 -= dy_div_dx * sigma_xy_dl / 2.0;
        partial_y2 -= sigma_y_dl / 2.0;

    } else {
        if(has_im1_j_flux2_sx) {
            if(neighbours[LEFT]) {
                if(neighbours[LEFT]->sigma.x > 0 && grid_cell->sigma.x > 0) {
                    // J_x_i-1/2,j f2
                    real_cpu local_sigma_x_dl = 2.0 / (1.0 / neighbours[LEFT]->sigma.x + 1.0 / grid_cell->sigma.x);
                    partial_x -= local_sigma_x_dl;
                }

                if(neighbours[LEFT]->sigma.xy > 0 && grid_cell->sigma.xy > 0) {
                    // J_y_i,j-1/2 f1
                    real_cpu local_sigma_xy_dl = 2.0 / (1.0 / neighbours[LEFT]->sigma.xy + 1.0 / grid_cell->sigma.xy);
                    partial_y2 -= dy_div_dx * local_sigma_xy_dl;
                }
            }
        }
        if(has_im1_j_flux2_sxy) {
            if(neighbours[DOWN]) {
                if(neighbours[DOWN]->sigma.xy > 0 && grid_cell->sigma.xy > 0) {
                    // J_x_i-1/2,j f2
                    real_cpu local_sigma_xy_dl = 2.0 / (1.0 / neighbours[DOWN]->sigma.xy + 1.0 / grid_cell->sigma.xy);
                    partial_x -= dx_div_dy * local_sigma_xy_dl;
                }

                // J_y_i,j-1/2 f1
                real_cpu local_sigma_y_dl = 2.0 / (1.0 / neighbours[DOWN]->sigma.y + 1.0 / grid_cell->sigma.y);
                partial_y2 -= local_sigma_y_dl;
            }
        }
    }

    if(has_im1_j_flux1 && has_im1_j_flux2) {
        partial_x /= 2.0;
    }

    if(has_i_jm1_flux1 && has_i_jm1_flux2) {
        partial_y2 /= 2.0;
    }

    elements[0].value += -(partial_x + partial_y2);

    // All neighbours
    for(int direction = 0; direction < NUM_DIRECTIONS; direction++) {
        if(neighbours[direction]) {

            struct element new_element;
            new_element.value = 0.0;

            switch(direction) {
            case FRONT:
                break;
            case BACK:
                break;
            case TOP: // i, j+1
                partial_x = 0.0;
                partial_y1 = 0.0;
                if(has_all_n_for_ip1_j_flux1) {
                    //[TOP_RIGHT];
                    partial_x -= sigma_x_tr / 2.0;
                    partial_x += dx_div_dy * sigma_xy_tr / 2.0;

                    partial_y1 -= dy_div_dx * sigma_xy_tr / 2.0;
                    partial_y1 += sigma_y_tr / 2.0;

                } else {
                    if(has_ip1_j_flux1_sx) {
                        if(!neighbours[RIGHT]) {
                            real_cpu local_sigma_x_tr =
                                2.0 / (1.0 / neighbours[TOP]->sigma.x + 1.0 / neighbours[TOP_RIGHT]->sigma.x);
                            partial_x -= local_sigma_x_tr;

                            if(neighbours[TOP]->sigma.xy > 0 && neighbours[TOP_RIGHT]->sigma.xy > 0) {
                                real_cpu local_sigma_xy_tr =
                                    2.0 / (1.0 / neighbours[TOP]->sigma.xy + 1.0 / neighbours[TOP_RIGHT]->sigma.xy);
                                partial_y1 -= dy_div_dx * local_sigma_xy_tr;
                            }
                        }
                    }

                    if(has_ip1_j_flux1_sxy) {
                        if(neighbours[TOP]->sigma.xy > 0 && grid_cell->sigma.xy > 0) {
                            real_cpu local_sigma_xy_tr =
                                2.0 / (1.0 / neighbours[TOP]->sigma.xy + 1.0 / grid_cell->sigma.xy);
                            partial_x += dx_div_dy * local_sigma_xy_tr;
                        }

                        real_cpu local_sigma_y_tr = 2.0 / (1.0 / neighbours[TOP]->sigma.y + 1.0 / grid_cell->sigma.y);
                        partial_y1 += local_sigma_y_tr;
                    }
                }

                if(has_ip1_j_flux1 && has_ip1_j_flux2) {
                    partial_x /= 2.0;
                }

                new_element.value += -partial_x;

                partial_x = 0.0;
                if(has_all_n_for_im1_j_flux1) {
                    //[TOP_LEFT];
                    partial_x -= sigma_x_tl / 2.0;
                    partial_x -= dx_div_dy * sigma_xy_tl / 2.0;

                    partial_y1 += dy_div_dx * sigma_xy_tl / 2.0;
                    partial_y1 += sigma_y_tl / 2.0;

                } else {

                    if(has_im1_j_flux1_sx) {
                        if(!neighbours[LEFT]) {
                            real_cpu local_sigma_x_tl =
                                2.0 / (1.0 / neighbours[TOP]->sigma.x + 1.0 / neighbours[TOP_LEFT]->sigma.x);
                            partial_x -= local_sigma_x_tl;

                            if(neighbours[TOP]->sigma.xy > 0 && neighbours[TOP_LEFT]->sigma.xy > 0) {
                                real_cpu local_sigma_xy_tl =
                                    2.0 / (1.0 / neighbours[TOP]->sigma.xy + 1.0 / neighbours[TOP_LEFT]->sigma.xy);
                                partial_y1 += dy_div_dx * local_sigma_xy_tl;
                            }
                        }
                    }

                    if(has_im1_j_flux1_sxy) {
                        if(neighbours[TOP]->sigma.xy > 0 && grid_cell->sigma.xy > 0) {
                            real_cpu local_sigma_xy_tl =
                                2.0 / (1.0 / neighbours[TOP]->sigma.xy + 1.0 / grid_cell->sigma.xy);
                            partial_x -= dx_div_dy * local_sigma_xy_tl;
                        }

                        real_cpu local_sigma_y_tl = 2.0 / (1.0 / neighbours[TOP]->sigma.y + 1.0 / grid_cell->sigma.y);
                        partial_y1 += local_sigma_y_tl;
                    }
                }

                if(has_im1_j_flux1 && has_im1_j_flux2) {
                    partial_x /= 2.0;
                }

                if(has_i_jp1_flux1 && has_i_jp1_flux2) {
                    partial_y1 /= 2.0;
                }

                new_element.value += -(partial_x + partial_y1);

                break;
            case DOWN: // i, j-1

                partial_x = 0.0;
                partial_y1 = 0.0;

                //[DOWN_RIGHT];
                if(has_all_n_for_ip1_j_flux2) {
                    partial_x -= sigma_x_dr / 2.0;
                    partial_x -= dx_div_dy * sigma_xy_dr / 2.0;

                    partial_y1 += dy_div_dx * sigma_xy_dr / 2.0;
                    partial_y1 += sigma_y_dr / 2.0;

                } else {
                    if(has_ip1_j_flux2_sx) {
                        if(!neighbours[RIGHT]) {
                            real_cpu local_sigma_x_dr =
                                2.0 / (1.0 / neighbours[DOWN]->sigma.x + 1.0 / neighbours[DOWN_RIGHT]->sigma.x);
                            partial_x -= local_sigma_x_dr;

                            if(neighbours[DOWN]->sigma.xy > 0 && neighbours[DOWN_RIGHT]->sigma.xy > 0) {
                                real_cpu local_sigma_xy_dr =
                                    2.0 / (1.0 / neighbours[DOWN]->sigma.xy + 1.0 / neighbours[DOWN_RIGHT]->sigma.xy);
                                partial_y1 += dy_div_dx * local_sigma_xy_dr;
                            }
                        }
                    }

                    if(has_ip1_j_flux2_sxy) {

                        if(neighbours[DOWN]->sigma.xy > 0 && grid_cell->sigma.xy > 0) {
                            real_cpu local_sigma_xy_dr =
                                2.0 / (1.0 / neighbours[DOWN]->sigma.xy + 1.0 / grid_cell->sigma.xy);
                            partial_x -= dx_div_dy * local_sigma_xy_dr;
                        }

                        real_cpu local_sigma_y_dr = 2.0 / (1.0 / neighbours[DOWN]->sigma.y + 1.0 / grid_cell->sigma.y);
                        partial_y1 += local_sigma_y_dr;
                    }
                }

                if(has_ip1_j_flux1 && has_ip1_j_flux2) {
                    partial_x /= 2.0;
                }

                new_element.value += -partial_x;

                partial_x = 0.0;
                //[DOWN_LEFT];
                if(has_all_n_for_im1_j_flux2) {
                    partial_x -= sigma_x_dl / 2.0;
                    partial_x += dx_div_dy * sigma_xy_dl / 2.0;

                    partial_y1 -= dy_div_dx * sigma_xy_dl / 2.0;
                    partial_y1 += sigma_y_dl / 2.0;

                } else {
                    if(has_im1_j_flux2_sx) {
                        if(!neighbours[LEFT]) {
                            real_cpu local_sigma_x_dl =
                                2.0 / (1.0 / neighbours[DOWN]->sigma.x + 1.0 / neighbours[DOWN_LEFT]->sigma.x);
                            partial_x -= local_sigma_x_dl;

                            if(neighbours[DOWN]->sigma.xy > 0 && neighbours[DOWN_LEFT]->sigma.xy > 0) {
                                real_cpu local_sigma_xy_dl =
                                    2.0 / (1.0 / neighbours[DOWN]->sigma.xy + 1.0 / neighbours[DOWN_LEFT]->sigma.xy);
                                partial_y1 -= dy_div_dx * local_sigma_xy_dl;
                            }
                        }
                    }

                    if(has_im1_j_flux2_sxy) {
                        if(neighbours[DOWN]->sigma.xy > 0 && grid_cell->sigma.xy > 0) {
                            real_cpu local_sigma_xy_dl =
                                2.0 / (1.0 / neighbours[DOWN]->sigma.xy + 1.0 / grid_cell->sigma.xy);
                            partial_x += dx_div_dy * local_sigma_xy_dl;
                        }

                        real_cpu local_sigma_y_dl = 2.0 / (1.0 / neighbours[DOWN]->sigma.y + 1.0 / grid_cell->sigma.y);
                        partial_y1 += local_sigma_y_dl;
                    }
                }

                if(has_i_jm1_flux1 && has_i_jm1_flux2) {
                    partial_y1 /= 2.0;
                }

                if(has_im1_j_flux1 && has_im1_j_flux2) {
                    partial_x /= 2.0;
                }

                new_element.value += -(partial_x + partial_y1);
                break;

            case RIGHT: // i+1, j
                partial_x = 0.0;
                partial_y1 = 0.0;
                partial_y2 = 0.0;

                if(has_all_n_for_ip1_j_flux1) {
                    //[TOP_RIGHT];
                    partial_x += sigma_x_tr / 2.0;
                    partial_x -= dx_div_dy * sigma_xy_tr / 2.0;

                    partial_y1 += dy_div_dx * sigma_xy_tr / 2.0;
                    partial_y1 -= sigma_y_tr / 2.0;

                } else {
                    if(has_ip1_j_flux1_sx) {
                        real_cpu local_sigma_x_tr = 2.0 / (1.0 / neighbours[RIGHT]->sigma.x + 1.0 / grid_cell->sigma.x);
                        partial_x += local_sigma_x_tr;

                        if(neighbours[RIGHT]->sigma.xy > 0 && grid_cell->sigma.xy > 0) {
                            real_cpu local_sigma_xy_tr =
                                2.0 / (1.0 / neighbours[RIGHT]->sigma.xy + 1.0 / grid_cell->sigma.xy);
                            partial_y1 += dy_div_dx * local_sigma_xy_tr;
                        }
                    }

                    if(has_ip1_j_flux1_sxy) {
                        if(!neighbours[TOP]) {
                            if(neighbours[RIGHT]->sigma.xy > 0 && neighbours[TOP_RIGHT]->sigma.xy > 0) {
                                real_cpu local_sigma_xy_tr =
                                    2.0 / (1.0 / neighbours[RIGHT]->sigma.xy + 1.0 / neighbours[TOP_RIGHT]->sigma.xy);
                                partial_x -= dx_div_dy * local_sigma_xy_tr;
                            }

                            real_cpu local_sigma_y_tr =
                                2.0 / (1.0 / neighbours[RIGHT]->sigma.y + 1.0 / neighbours[TOP_RIGHT]->sigma.y);
                            partial_y1 -= local_sigma_y_tr;
                        }
                    }
                }

                if(has_all_n_for_ip1_j_flux2) {
                    partial_x += sigma_x_dr / 2.0;
                    partial_x += dx_div_dy * sigma_xy_dr / 2.0;

                    partial_y2 -= dy_div_dx * sigma_xy_dr / 2.0;
                    partial_y2 -= sigma_y_dr / 2.0;

                } else {
                    if(has_ip1_j_flux2_sx) {
                        real_cpu local_sigma_x_dr = 2.0 / (1.0 / neighbours[RIGHT]->sigma.x + 1.0 / grid_cell->sigma.x);
                        partial_x += local_sigma_x_dr;

                        if(neighbours[RIGHT]->sigma.xy > 0 && grid_cell->sigma.xy > 0) {
                            real_cpu local_sigma_xy_dr =
                                2.0 / (1.0 / neighbours[RIGHT]->sigma.xy + 1.0 / grid_cell->sigma.xy);
                            partial_y2 -= dy_div_dx * local_sigma_xy_dr;
                        }
                    }

                    if(has_ip1_j_flux2_sxy) {
                        if(!neighbours[DOWN]) {
                            if(neighbours[RIGHT]->sigma.xy > 0 && neighbours[DOWN_RIGHT]->sigma.xy > 0) {
                                real_cpu local_sigma_xy_dr =
                                    2.0 / (1.0 / neighbours[RIGHT]->sigma.xy + 1.0 / neighbours[DOWN_RIGHT]->sigma.xy);
                                partial_x += dx_div_dy * local_sigma_xy_dr;
                            }

                            real_cpu local_sigma_y_dr =
                                2.0 / (1.0 / neighbours[RIGHT]->sigma.y + 1.0 / neighbours[DOWN_RIGHT]->sigma.y);
                            partial_y2 -= local_sigma_y_dr;
                        }
                    }
                }

                if(has_ip1_j_flux1 && has_ip1_j_flux2) {
                    partial_x /= 2.0;
                }

                if(has_i_jp1_flux1 && has_i_jp1_flux2) {
                    partial_y1 /= 2.0;
                }

                if(has_i_jm1_flux1 && has_i_jm1_flux2) {
                    partial_y2 /= 2.0;
                }

                new_element.value += -(partial_x + partial_y1 + partial_y2);

                break;

            case LEFT: // i-1, j

                partial_x = 0;
                partial_y1 = 0.0;
                partial_y2 = 0.0;

                if(has_all_n_for_im1_j_flux1) {
                    partial_x += sigma_x_tl / 2.0;
                    partial_x += dx_div_dy * sigma_xy_tl / 2.0;

                    partial_y1 -= dy_div_dx * sigma_xy_tl / 2.0;
                    partial_y1 -= sigma_y_tl / 2.0;

                } else {
                    if(has_im1_j_flux1_sx) {
                        real_cpu local_sigma_x_tl = 2.0 / (1.0 / neighbours[LEFT]->sigma.x + 1.0 / grid_cell->sigma.x);
                        partial_x += local_sigma_x_tl;

                        if(neighbours[LEFT]->sigma.xy > 0 && grid_cell->sigma.xy > 0) {
                            real_cpu local_sigma_xy_tl =
                                2.0 / (1.0 / neighbours[LEFT]->sigma.xy + 1.0 / grid_cell->sigma.xy);
                            partial_y1 -= dy_div_dx * local_sigma_xy_tl;
                        }
                    }

                    if(has_im1_j_flux1_sxy) {
                        if(!neighbours[TOP]) {
                            if(neighbours[LEFT]->sigma.xy > 0 &&  neighbours[TOP_LEFT]->sigma.xy > 0) {
                                real_cpu local_sigma_xy_tl =
                                    2.0 / (1.0 / neighbours[LEFT]->sigma.xy + 1.0 / neighbours[TOP_LEFT]->sigma.xy);
                                partial_x += dx_div_dy * local_sigma_xy_tl;
                            }

                            real_cpu local_sigma_y_tl =
                                2.0 / (1.0 / neighbours[LEFT]->sigma.y + 1.0 / neighbours[TOP_LEFT]->sigma.y);
                            partial_y1 -= local_sigma_y_tl;
                        }
                    }
                }

                if(has_all_n_for_im1_j_flux2) {
                    partial_x += sigma_x_dl / 2.0;
                    partial_x -= dx_div_dy * sigma_xy_dl / 2.0;

                    partial_y2 += dy_div_dx * sigma_xy_dl / 2.0;
                    partial_y2 -= sigma_y_dl / 2.0;

                } else {
                    if(has_im1_j_flux2_sx) {
                        real_cpu local_sigma_x_dl = 2.0 / (1.0 / neighbours[LEFT]->sigma.x + 1.0 / grid_cell->sigma.x);
                        partial_x += local_sigma_x_dl;

                        if(neighbours[LEFT]->sigma.xy > 0 &&  grid_cell->sigma.xy > 0) {
                            real_cpu local_sigma_xy_dl =
                                2.0 / (1.0 / neighbours[LEFT]->sigma.xy + 1.0 / grid_cell->sigma.xy);
                            partial_y2 += dy_div_dx * local_sigma_xy_dl;
                        }
                    }
                    if(has_im1_j_flux2_sxy) {
                        if(!neighbours[DOWN]) {
                            if(neighbours[LEFT]->sigma.xy > 0 &&  neighbours[DOWN_LEFT]->sigma.xy > 0) {
                                real_cpu local_sigma_xy_dl =
                                    2.0 / (1.0 / neighbours[LEFT]->sigma.xy + 1.0 / neighbours[DOWN_LEFT]->sigma.xy);
                                partial_x -= dx_div_dy * local_sigma_xy_dl;
                            }

                            real_cpu local_sigma_y_dl =
                                2.0 / (1.0 / neighbours[LEFT]->sigma.y + 1.0 / neighbours[DOWN_LEFT]->sigma.y);
                            partial_y2 -= local_sigma_y_dl;
                        }
                    }
                }

                if(has_im1_j_flux1 && has_im1_j_flux2) {
                    partial_x /= 2.0;
                }

                if(has_i_jp1_flux1 && has_i_jp1_flux2) {
                    partial_y1 /= 2.0;
                }

                if(has_i_jm1_flux1 && has_i_jm1_flux2) {
                    partial_y2 /= 2.0;
                }

                new_element.value += -(partial_x + partial_y1 + partial_y2);

                break;

            case FRONT_TOP:
                break;
            case FRONT_DOWN:
                break;
            case FRONT_RIGHT:
                break;
            case FRONT_LEFT:
                break;
            case FRONT_TOP_RIGHT:
                break;
            case FRONT_TOP_LEFT:
                break;
            case FRONT_DOWN_RIGHT:
                break;
            case FRONT_DOWN_LEFT:
                break;
            case BACK_TOP:
                break;
            case BACK_DOWN:
                break;
            case BACK_RIGHT:
                break;
            case BACK_LEFT:
                break;
            case BACK_TOP_RIGHT:
                break;
            case BACK_TOP_LEFT:
                break;
            case BACK_DOWN_RIGHT:
                break;
            case BACK_DOWN_LEFT:
                break;

            case TOP_RIGHT: // i+1, j+1

                partial_x = 0.0;
                partial_y1 = 0.0;
                if(has_all_n_for_ip1_j_flux1) {
                    //[TOP_RIGHT];
                    partial_x += sigma_x_tr / 2.0;
                    partial_x += dx_div_dy * sigma_xy_tr / 2.0;

                    partial_y1 += dy_div_dx * sigma_xy_tr / 2.0;
                    partial_y1 += sigma_y_tr / 2.0;

                } else {
                    if(has_ip1_j_flux1_sx && neighbours[RIGHT] == NULL) {
                        real_cpu local_sigma_x_tr =
                            2.0 / (1.0 / neighbours[TOP_RIGHT]->sigma.x + 1.0 / neighbours[TOP]->sigma.x);
                        partial_x += local_sigma_x_tr;

                        if(neighbours[TOP_RIGHT]->sigma.xy > 0 &&  neighbours[TOP]->sigma.xy > 0) {
                            real_cpu local_sigma_xy_tr =
                                2.0 / (1.0 / neighbours[TOP_RIGHT]->sigma.xy + 1.0 / neighbours[TOP]->sigma.xy);
                            partial_y1 += dy_div_dx * local_sigma_xy_tr;
                        }
                    }

                    if(has_ip1_j_flux1_sxy && neighbours[TOP] == NULL) {

                        if(neighbours[TOP_RIGHT]->sigma.xy > 0 &&  neighbours[TOP]->sigma.xy > 0) {
                            real_cpu local_sigma_xy_tr =
                                2.0 / (1.0 / neighbours[TOP_RIGHT]->sigma.xy + 1.0 / neighbours[RIGHT]->sigma.xy);
                            partial_x += dx_div_dy * local_sigma_xy_tr;
                        }

                        real_cpu local_sigma_y_tr =
                            2.0 / (1.0 / neighbours[TOP_RIGHT]->sigma.y + 1.0 / neighbours[RIGHT]->sigma.y);
                        partial_y1 += local_sigma_y_tr;
                    }
                }

                if(has_ip1_j_flux1 && has_ip1_j_flux2) {
                    partial_x /= 2.0;
                }

                if(has_i_jp1_flux1 && has_i_jp1_flux2) {
                    partial_y1 /= 2.0;
                }

                new_element.value += -(partial_x + partial_y1);
                break;

            case TOP_LEFT: // i-1, j+1

                partial_x = 0.0;
                partial_y1 = 0.0;
                if(has_all_n_for_im1_j_flux1) {
                    partial_x += sigma_x_tl / 2.0;
                    partial_x -= dx_div_dy * sigma_xy_tl / 2.0;

                    partial_y1 -= dy_div_dx * sigma_xy_tl / 2.0;
                    partial_y1 += sigma_y_tl / 2.0;

                } else {
                    if(has_im1_j_flux1_sx && neighbours[LEFT] == NULL) {
                        real_cpu local_sigma_x_tl =
                            2.0 / (1.0 / neighbours[TOP_LEFT]->sigma.x + 1.0 / neighbours[TOP]->sigma.x);
                        partial_x += local_sigma_x_tl;

                        if(neighbours[TOP_LEFT]->sigma.xy > 0 &&  neighbours[TOP]->sigma.xy > 0) {
                            real_cpu local_sigma_xy_tl =
                                2.0 / (1.0 / neighbours[TOP_LEFT]->sigma.xy + 1.0 / neighbours[TOP]->sigma.xy);
                            partial_y1 -= dy_div_dx * local_sigma_xy_tl;
                        }
                    }

                    if(has_im1_j_flux1_sxy && neighbours[TOP] == NULL) {

                        if(neighbours[TOP_LEFT]->sigma.xy > 0 &&  neighbours[TOP]->sigma.xy > 0) {
                            real_cpu local_sigma_xy_tl =
                                2.0 / (1.0 / neighbours[TOP_LEFT]->sigma.xy + 1.0 / neighbours[LEFT]->sigma.xy);
                            partial_x -= dx_div_dy * local_sigma_xy_tl;
                        }

                        real_cpu local_sigma_y_tl =
                            2.0 / (1.0 / neighbours[TOP_LEFT]->sigma.y + 1.0 / neighbours[LEFT]->sigma.y);
                        partial_y1 += local_sigma_y_tl;
                    }
                }

                if(has_im1_j_flux1 && has_im1_j_flux2) {
                    partial_x /= 2.0;
                }

                if(has_i_jp1_flux1 && has_i_jp1_flux2) {
                    partial_y1 /= 2.0;
                }

                new_element.value += -(partial_x + partial_y1);
                break;

            case DOWN_RIGHT: // i+1, j-1

                partial_x = 0.0;
                partial_y1 = 0.0;

                if(has_all_n_for_ip1_j_flux2) {
                    partial_x += sigma_x_dr / 2.0;
                    partial_x -= dx_div_dy * sigma_xy_dr / 2.0;

                    partial_y1 -= dy_div_dx * sigma_xy_dr / 2.0;
                    partial_y1 += sigma_y_dr / 2.0;

                } else {
                    if(has_ip1_j_flux2_sx && neighbours[RIGHT] == NULL) {
                        real_cpu local_sigma_x_dr =
                            2.0 / (1.0 / neighbours[DOWN_RIGHT]->sigma.x + 1.0 / neighbours[DOWN]->sigma.x);
                        partial_x += local_sigma_x_dr;

                        if(neighbours[DOWN_RIGHT]->sigma.xy > 0 &&  neighbours[DOWN]->sigma.xy > 0) {
                            real_cpu local_sigma_xy_dr =
                                2.0 / (1.0 / neighbours[DOWN_RIGHT]->sigma.xy + 1.0 / neighbours[DOWN]->sigma.xy);
                            partial_y1 -= dy_div_dx * local_sigma_xy_dr;
                        }
                    }

                    if(has_ip1_j_flux2_sxy && neighbours[DOWN] == NULL) {
                        if(neighbours[DOWN_RIGHT]->sigma.xy > 0 &&  neighbours[RIGHT]->sigma.xy > 0) {
                            real_cpu local_sigma_xy_dr =
                                2.0 / (1.0 / neighbours[DOWN_RIGHT]->sigma.xy + 1.0 / neighbours[RIGHT]->sigma.xy);
                            partial_x -= dx_div_dy * local_sigma_xy_dr;
                        }

                        real_cpu local_sigma_y_dr =
                            2.0 / (1.0 / neighbours[DOWN_RIGHT]->sigma.y + 1.0 / neighbours[RIGHT]->sigma.y);
                        partial_y1 += local_sigma_y_dr;
                    }
                }

                if(has_ip1_j_flux1 && has_ip1_j_flux2) {
                    partial_x /= 2.0;
                }

                if(has_i_jm1_flux1 && has_i_jm1_flux2) {
                    partial_y1 /= 2.0;
                }

                new_element.value += -(partial_x + partial_y1);

                break;
            case DOWN_LEFT: // i-1, j-1

                partial_x = 0.0;
                partial_y1 = 0.0;

                if(has_all_n_for_im1_j_flux2) {
                    partial_x += sigma_x_dl / 2.0;
                    partial_x += dx_div_dy * sigma_xy_dl / 2.0;

                    partial_y1 += dy_div_dx * sigma_xy_dl / 2.0;
                    partial_y1 += sigma_y_dl / 2.0;

                } else {
                    if(has_im1_j_flux2_sx && neighbours[LEFT] == NULL) {
                        real_cpu local_sigma_x_dl =
                            2.0 / (1.0 / neighbours[DOWN_LEFT]->sigma.x + 1.0 / neighbours[DOWN]->sigma.x);
                        partial_x += local_sigma_x_dl;

                        real_cpu local_sigma_xy_dl =
                            2.0 / (1.0 / neighbours[DOWN_LEFT]->sigma.xy + 1.0 / neighbours[DOWN]->sigma.xy);
                        partial_y1 += dy_div_dx * local_sigma_xy_dl;
                    }

                    if(has_im1_j_flux2_sxy && neighbours[DOWN] == NULL) {

                        if(neighbours[DOWN_LEFT]->sigma.xy > 0 &&  neighbours[LEFT]->sigma.xy > 0) {
                            real_cpu local_sigma_xy_dl =
                                2.0 / (1.0 / neighbours[DOWN_LEFT]->sigma.xy + 1.0 / neighbours[LEFT]->sigma.xy);
                            partial_x += dx_div_dy * local_sigma_xy_dl;
                        }

                        real_cpu local_sigma_y_dl =
                            2.0 / (1.0 / neighbours[DOWN_LEFT]->sigma.y + 1.0 / neighbours[LEFT]->sigma.y);
                        partial_y1 += local_sigma_y_dl;
                    }
                }

                if(has_im1_j_flux1 && has_im1_j_flux2) {
                    partial_x /= 2.0;
                }

                if(has_i_jm1_flux1 && has_i_jm1_flux2) {
                    partial_y1 /= 2.0;
                }

                new_element.value += -(partial_x + partial_y1);
                break;
            default:
                break;
            }

            new_element.column = neighbours[direction]->grid_position;
            new_element.cell = neighbours[direction];
            arrput(grid_cell->elements, new_element);
        }
    }
}