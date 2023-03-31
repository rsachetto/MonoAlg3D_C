#include "gui_mesh_helpers.h"
#include "gui_colors.h"
#include "raylib_ext.h"
#include <float.h>

struct mesh_info *new_mesh_info() {
    struct mesh_info *mesh_info = (struct mesh_info *)malloc(sizeof(struct mesh_info));
    mesh_info->center_calculated = false;
    return mesh_info;
}

#define SET_MAX(coord)                                                                                                                                         \
    do {                                                                                                                                                       \
        result.coord = mesh_max.coord;                                                                                                                         \
        if(mesh_max.coord != mesh_min.coord)                                                                                                                   \
            result.coord = (mesh_max.coord - mesh_min.coord) / 2.0f;                                                                                           \
    } while(0)

#define SET_MESH_MIN_MAX(vec1, vec2, vec3)                                                                                                                     \
    do {                                                                                                                                                       \
        mesh_info->vec1.x = vec2.x + vec3.x / 2.0f;                                                                                                            \
        mesh_info->vec1.y = vec2.y + vec3.y / 2.0f;                                                                                                            \
        mesh_info->vec1.z = vec2.z + vec3.z / 2.0f;                                                                                                            \
    } while(0)

static inline Vector3 get_max_min(Vector3 mesh_max, Vector3 mesh_min, Vector3 mesh_max_d, Vector3 mesh_min_d, struct mesh_info *mesh_info) {

    Vector3 result = V3_SAME(0.0f);

    SET_MAX(x);
    SET_MAX(y);
    SET_MAX(z);

    mesh_info->center_calculated = true;

    SET_MESH_MIN_MAX(max_size, mesh_max, mesh_max_d);
    SET_MESH_MIN_MAX(min_size, mesh_min, mesh_min_d);

    return result;
}

#define SET_MAX_MIN(coord)                                                                                                                                     \
    if(mesh_center.coord > mesh_max.coord) {                                                                                                                   \
        mesh_max.coord = mesh_center.coord;                                                                                                                    \
        mesh_max_d.coord = d.coord;                                                                                                                            \
    } else if(mesh_center.coord < mesh_min.coord) {                                                                                                            \
        mesh_min.coord = mesh_center.coord;                                                                                                                    \
        mesh_min_d.coord = d.coord;                                                                                                                            \
    }

#define ITERATE_OVER_CELLS()                                                                                                                                   \
    for(uint32_t i = 0; i < n_active; i++) {                                                                                                                   \
        grid_cell = ac[i];                                                                                                                                     \
        d.x = grid_cell->discretization.x;                                                                                                                     \
        d.y = grid_cell->discretization.y;                                                                                                                     \
        d.z = grid_cell->discretization.z;                                                                                                                     \
        mesh_center.x = grid_cell->center.x;                                                                                                                   \
        mesh_center.y = grid_cell->center.y;                                                                                                                   \
        mesh_center.z = grid_cell->center.z;                                                                                                                   \
        SET_MAX_MIN(x);                                                                                                                                        \
        SET_MAX_MIN(y);                                                                                                                                        \
        SET_MAX_MIN(z);                                                                                                                                        \
    }

Vector3 find_mesh_center(struct grid *grid_to_draw, struct mesh_info *mesh_info) {

    uint32_t n_active = grid_to_draw->num_active_cells;
    struct cell_node **ac = grid_to_draw->active_cells;
    struct cell_node *grid_cell;

    Vector3 mesh_max = V3_SAME(FLT_MIN);
    Vector3 mesh_min = V3_SAME(FLT_MAX);

    Vector3 mesh_max_d = V3_SAME(FLT_MIN);
    Vector3 mesh_min_d = V3_SAME(FLT_MAX);

    Vector3 d;
    Vector3 mesh_center;

    if(ac) {
        ITERATE_OVER_CELLS()
    } else {
        if(grid_to_draw->purkinje) {
            n_active = grid_to_draw->purkinje->number_of_purkinje_cells;
            ac = grid_to_draw->purkinje->purkinje_cells;
            ITERATE_OVER_CELLS()
        }
    }

    return get_max_min(mesh_max, mesh_min, mesh_max_d, mesh_min_d, mesh_info);
}

#define SET_CENTER()                                                                                                                                           \
    do {                                                                                                                                                       \
        mesh_center.x = (float)points[cells[i]].x + d.x / 2.0f;                                                                                                \
        mesh_center.y = (float)points[cells[i]].y + d.y / 2.0f;                                                                                                \
        mesh_center.z = (float)points[cells[i]].z + d.z / 2.0f;                                                                                                \
        SET_MAX_MIN(x);                                                                                                                                        \
        SET_MAX_MIN(y);                                                                                                                                        \
        SET_MAX_MIN(z);                                                                                                                                        \
    } while(0)

Vector3 find_mesh_center_vtk(struct vtk_unstructured_grid *grid_to_draw, struct mesh_info *mesh_info) {


    uint32_t n_active = grid_to_draw->num_cells;

    Vector3 mesh_max = V3_SAME(FLT_MIN);
    Vector3 mesh_min = V3_SAME(FLT_MAX);

    Vector3 mesh_max_d = V3_SAME(FLT_MIN);
    Vector3 mesh_min_d = V3_SAME(FLT_MAX);

    int64_t *cells = grid_to_draw->cells;
    point3d_array points = grid_to_draw->points;

    assert(points);
    uint32_t num_points = grid_to_draw->points_per_cell;

    Vector3 d;
    Vector3 mesh_center;

    if(cells) {
        for(uint32_t i = 0; i < n_active * num_points; i += num_points) {
            d.x = (float)fabs((points[cells[i]].x - points[cells[i + 1]].x));
            d.y = (float)fabs((points[cells[i]].y - points[cells[i + 3]].y));
            d.z = (float)fabs((points[cells[i]].z - points[cells[i + 4]].z));
            SET_CENTER();
        }
    } else if(grid_to_draw->purkinje) {

        n_active = grid_to_draw->purkinje->num_cells;
        cells = grid_to_draw->purkinje->cells;
        points = grid_to_draw->purkinje->points;
        num_points = grid_to_draw->purkinje->points_per_cell;

        for(uint32_t i = 0; i < n_active * num_points; i += num_points) {
            d.x = (float)fabs((points[cells[i]].x - points[cells[i + 1]].x));
            d.y = (float)fabs((points[cells[i]].y - points[cells[i + 1]].y));
            d.z = (float)fabs((points[cells[i]].z - points[cells[i + 1]].z));
            SET_CENTER();
        }
    }

    return get_max_min(mesh_max, mesh_min, mesh_max_d, mesh_min_d, mesh_info);
}

static void check_mouse_over_volume(struct voxel *voxel, struct gui_state *gui_state) {

    Vector3 p_draw = voxel->position_draw;
    Vector3 p_mesh = voxel->position_mesh;

    Vector3 cube_size = voxel->size;

    float csx = cube_size.x;
    float csy = cube_size.y;
    float csz = cube_size.z;

    BoundingBox bb = (BoundingBox){(Vector3){p_draw.x - csx / 2, p_draw.y - csy / 2, p_draw.z - csz / 2},
                                   (Vector3){p_draw.x + csx / 2, p_draw.y + csy / 2, p_draw.z + csz / 2}};


    RayCollision collision_mouse_over = GetRayCollisionBox(gui_state->ray_mouse_over, bb);

    if(collision_mouse_over.hit) {
        if(gui_state->ray_mouse_over_hit_distance > collision_mouse_over.distance) {
            gui_state->current_mouse_over_volume.position_draw = (Vector3){p_mesh.x, p_mesh.y, p_mesh.z};
            gui_state->current_mouse_over_volume.matrix_position = voxel->matrix_position;
            gui_state->ray_mouse_over_hit_distance = collision_mouse_over.distance;
        }
    }

}

static bool check_volume_selection(struct voxel *voxel, struct gui_state *gui_state) {

    Vector3 p_draw = voxel->position_draw;
    Vector3 p_mesh = voxel->position_mesh;

    Vector3 cube_size = voxel->size;

    float csx = cube_size.x;
    float csy = cube_size.y;
    float csz = cube_size.z;

    BoundingBox bb = (BoundingBox){(Vector3){p_draw.x - csx / 2, p_draw.y - csy / 2, p_draw.z - csz / 2},
                                   (Vector3){p_draw.x + csx / 2, p_draw.y + csy / 2, p_draw.z + csz / 2}};

    RayCollision collision = {0};

    if(gui_state->double_clicked) {
        collision = GetRayCollisionBox(gui_state->ray, bb);

        if(collision.hit) {
            if(gui_state->ray_hit_distance > collision.distance) {
                gui_state->current_selected_volume = *voxel;
                gui_state->ray_hit_distance = collision.distance;
            }

        } else {
            collision.hit = false;
        }
    } else {
        action_potential_array aps = (struct action_potential *)hmget(gui_state->ap_graph_config->selected_aps, p_mesh);
        if(aps != NULL) {
            // This is needed when we are using adaptive meshes
            hmput(gui_state->current_selected_volumes, p_mesh, *voxel);
        }
    }

    if(gui_state->found_volume.position_mesh.x >= 0 && gui_state->found_volume.position_mesh.y >= 0 && gui_state->found_volume.position_mesh.z >= 0) {
        if(gui_state->found_volume.position_mesh.x == p_mesh.x && gui_state->found_volume.position_mesh.y == p_mesh.y &&
           gui_state->found_volume.position_mesh.z == p_mesh.z) {
            gui_state->current_selected_volume = *voxel;
            collision.hit = true;
            gui_state->found_volume.position_mesh = (Vector3){-1, -1, -1};
        }
    }

    return collision.hit;
}

static void add_or_remove_selected_to_ap_graph(struct gui_shared_info *gui_config, struct gui_state *gui_state) {

    Vector3 p_mesh = gui_state->current_selected_volume.position_mesh;

    action_potential_array aps = (struct action_potential *)hmget(gui_state->ap_graph_config->selected_aps, p_mesh);

    if(aps == NULL) {
        struct action_potential ap;
        ap.t = gui_config->time;
        ap.v = gui_state->current_selected_volume.v;

        arrput(aps, ap);
        hmput(gui_state->ap_graph_config->selected_aps, p_mesh, aps);

        hmput(gui_state->current_selected_volumes, gui_state->current_selected_volume.position_mesh, gui_state->current_selected_volume);
    } else {
        (void)hmdel(gui_state->ap_graph_config->selected_aps, p_mesh);
        (void)hmdel(gui_state->current_selected_volumes, p_mesh);
    }
}

// we only trace the arrays that were previously added to the ap graph
static void trace_ap(struct gui_state *gui_state, struct voxel *voxel, float t) {

    if(hmlen(gui_state->ap_graph_config->selected_aps) == 0) {
        return;
    } else {
        Vector3 p_mesh = voxel->position_mesh;
        action_potential_array aps = NULL;

        aps = (struct action_potential *)hmget(gui_state->ap_graph_config->selected_aps, p_mesh);
        if(aps != NULL) {
            struct action_potential ap1;
            ap1.t = t;
            ap1.v = voxel->v;

            size_t aps_len = arrlen(aps);

            if(aps_len == 0 || ap1.t > aps[aps_len - 1].t) {
                arrput(aps, ap1);
                hmput(gui_state->ap_graph_config->selected_aps, p_mesh, aps);
            }
        }
    }
}

static void update_selected(bool collision, struct gui_state *gui_state, struct gui_shared_info *gui_config, Color *colors) {

    if(collision) {
        add_or_remove_selected_to_ap_graph(gui_config, gui_state);
    }

    int n = hmlen(gui_state->current_selected_volumes);

    for(int i = 0; i < n; i++) {
        uint32_t idx = gui_state->current_selected_volumes[i].value.draw_index;
        Color c = colors[idx];
        c.a = 0;
        colors[idx] = c;
    }
}

void draw_vtk_unstructured_grid(struct gui_shared_info *gui_config, struct gui_state *gui_state, int grid_mask, struct draw_context *draw_context) {

    static bool discretization_calc = false;

    struct vtk_unstructured_grid *grid_to_draw = gui_config->grid_info.vtk_grid;
    if (!grid_to_draw) return;

    int64_t *cells = grid_to_draw->cells;
    if(!cells) return;

    point3d_array points = grid_to_draw->points;
    if(!points) return;

    f32_array values = grid_to_draw->values;

    bool single_file = (gui_config->final_file_index == 0);

    uint32_t n_active = grid_to_draw->num_cells;
    uint32_t num_points = grid_to_draw->points_per_cell;
    int j = 0;

    struct voxel voxel;

    gui_state->ray_hit_distance = FLT_MAX;
    gui_state->ray_mouse_over_hit_distance = FLT_MAX;

    float min_v = 0;
    float max_v = 1;

    int extra_data_len = arrlen(grid_to_draw->extra_values);

    int current_data_index = gui_state->current_data_index;

    if(gui_config->final_file_index == 0) {
        if(current_data_index == -1) {
            min_v = grid_to_draw->min_v;
            max_v = grid_to_draw->max_v;
        } else if(extra_data_len > 0 && current_data_index < extra_data_len) {
            min_v = grid_to_draw->min_extra_value[current_data_index];
            max_v = grid_to_draw->max_extra_value[current_data_index];
        }
    } else {
        max_v = gui_config->max_v;
        min_v = gui_config->min_v;
    }

    if(min_v == max_v) {
        max_v = max_v + 0.00001f;
    }

    gui_config->min_v = min_v;
    gui_config->max_v = max_v;

    float time = gui_config->time;

    bool collision = false;

    Vector3 n = gui_state->plane_normal;
    Vector3 p = gui_state->plane_point;

    float scale = gui_state->mesh_scale_factor;
    Vector3 mesh_offset = gui_state->mesh_offset;

    bool slicing_mode = gui_state->slicing_mode;

    int count = 0;
    for(uint32_t i = 0; i < n_active * num_points; i += num_points) {

        if(grid_mask != 2) {
            if(grid_to_draw->cell_visibility && !grid_to_draw->cell_visibility[j]) {
                j += 1;
                continue;
            }
        }

        struct point_3d cell_point = points[cells[i]];
        float cell_point_x = (float) cell_point.x;
        float cell_point_y = (float) cell_point.y;
        float cell_point_z = (float) cell_point.z;

        float mesh_center_x, mesh_center_y, mesh_center_z;
        float dx, dy, dz;

        dx = fabs((cell_point_x - points[cells[i + 1]].x));
        dy = fabs((cell_point_y - points[cells[i + 3]].y));
        dz = fabs((cell_point_z - points[cells[i + 4]].z));

        if(!discretization_calc) {
            grid_to_draw->average_discretization.x += dx;
            grid_to_draw->average_discretization.y += dy;
            grid_to_draw->average_discretization.z += dz;
        }

        mesh_center_x = cell_point_x + dx / 2.0f;
        mesh_center_y = cell_point_y + dy / 2.0f;
        mesh_center_z = cell_point_z + dz / 2.0f;

        voxel.position_draw.x = (mesh_center_x - mesh_offset.x) / scale;
        voxel.position_draw.y = (mesh_center_y - mesh_offset.y) / scale;
        voxel.position_draw.z = (mesh_center_z - mesh_offset.z) / scale;

        if(slicing_mode) {
            Vector3 test = Vector3Subtract(voxel.position_draw, p);
            float side = Vector3DotProduct(test, n);

            if(side < 0) {
                j += 1;
                continue;
            }
        }

        voxel.v = 0.0f;

        if(current_data_index == -1) {
            if(values) {
                voxel.v = values[j];
            }
        } else if(extra_data_len > 0 && current_data_index < extra_data_len) {
            voxel.v = grid_to_draw->extra_values[current_data_index][j];
        }

        voxel.size.x = dx / scale;
        voxel.size.y = dy / scale;
        voxel.size.z = dz / scale;
        voxel.position_mesh = (Vector3){mesh_center_x, mesh_center_y, mesh_center_z};

        draw_context->translations[count] = MatrixTranslate(voxel.position_draw.x, voxel.position_draw.y, voxel.position_draw.z);
        draw_context->translations[count] = MatrixMultiply(MatrixScale(voxel.size.x, voxel.size.y, voxel.size.z), draw_context->translations[count]);

        draw_context->colors[count] = get_color((voxel.v - min_v) / (max_v - min_v), gui_state->voxel_alpha, gui_state->current_scale);
        voxel.draw_index = count;

        bool searching = gui_state->found_volume.position_mesh.x >= 0
                         && gui_state->found_volume.position_mesh.y >= 0
                         && gui_state->found_volume.position_mesh.z >= 0;

        if((gui_state->double_clicked || gui_config->adaptive || searching) && !single_file) {
            collision |= check_volume_selection(&voxel, gui_state);
        }

        if(gui_state->ctrl_pressed) {
            check_mouse_over_volume(&voxel, gui_state);
        }

        if(!single_file) {
            trace_ap(gui_state, &voxel, time);
        }

        count++;

        j += 1;
    }

    if(!single_file) {
        update_selected(collision, gui_state, gui_config, draw_context->colors);
    }

    if(!discretization_calc) {
        grid_to_draw->average_discretization.x /= count;
        grid_to_draw->average_discretization.y /= count;
        grid_to_draw->average_discretization.z /= count;
        discretization_calc = true;
    }

    gui_state->found_volume.position_mesh = (Vector3){-1, -1, -1};
    DrawMeshInstancedWithColors(draw_context, grid_mask, count);
}

void draw_vtk_purkinje_network(struct gui_shared_info *gui_config, struct gui_state *gui_state) {

    if(gui_config->grid_info.vtk_grid == NULL || gui_config->grid_info.vtk_grid->purkinje == NULL)
        return;

    struct vtk_unstructured_grid *grid_to_draw = gui_config->grid_info.vtk_grid->purkinje;

    int64_t *cells = grid_to_draw->cells;
    if(!cells)
        return;

    point3d_array points = grid_to_draw->points;
    if(!points)
        return;

    uint32_t n_active = grid_to_draw->num_cells;

    uint32_t num_points = grid_to_draw->points_per_cell;

    float min_v = gui_config->min_v;
    float max_v = gui_config->max_v;

    int j = 0;

    float scale = gui_state->mesh_scale_factor;
    Vector3 mesh_offset = gui_state->mesh_offset;

    for(uint32_t i = 0; i < n_active * num_points; i += num_points) {

        Vector3 start_pos;
        start_pos.x = ((float)points[cells[i]].x - mesh_offset.x) / scale;
        start_pos.y = ((float)points[cells[i]].y - mesh_offset.y) / scale;
        start_pos.z = ((float)points[cells[i]].z - mesh_offset.z) / scale;

        Vector3 end_pos;
        end_pos.x = ((float)points[cells[i + 1]].x - mesh_offset.x) / scale;
        end_pos.y = ((float)points[cells[i + 1]].y - mesh_offset.y) / scale;
        end_pos.z = ((float)points[cells[i + 1]].z - mesh_offset.z) / scale;

        Color c = BLUE;

        if(grid_to_draw->values) {
            float v = grid_to_draw->values[j];
            c = get_color((v - min_v) / (max_v - min_v), gui_state->voxel_alpha, gui_state->current_scale);
        }

        j++;

        DrawLine3D(start_pos, end_pos, c);
    }
}

void draw_alg_mesh(struct gui_shared_info *gui_config, struct gui_state *gui_state, int grid_mask, struct draw_context *draw_context) {

    struct grid *grid_to_draw = gui_config->grid_info.alg_grid;

    if(!grid_to_draw)
        return;

    struct voxel voxel;

    uint32_t n_active = grid_to_draw->num_active_cells;
    struct cell_node **ac = grid_to_draw->active_cells;

    if(ac) {

        float scale = gui_state->mesh_scale_factor;
        Vector3 mesh_offset = gui_state->mesh_offset;

        float offsetx_over_scale = mesh_offset.x / scale;
        float offsety_over_scale = mesh_offset.y / scale;
        float offsetz_over_scale = mesh_offset.z / scale;

        float min_v = gui_config->min_v;
        float max_v = gui_config->max_v;
        float time = gui_config->time;

        gui_state->ray_hit_distance = FLT_MAX;
        gui_state->ray_mouse_over_hit_distance = FLT_MAX;

        int count = 0;

        bool collision = false;

        /*
        Vector3 n = gui_state->plane_normal;
        Vector3 p = gui_state->plane_point;
        */
        for(uint32_t i = 0; i < n_active; i++) {

            if(grid_mask != 2) {
                if(!ac[i]->visible) {
                    continue;
                }
            }

            struct cell_node *grid_cell;

            grid_cell = ac[i];

            voxel.position_draw.x = (float)grid_cell->center.x / scale - offsetx_over_scale;
            voxel.position_draw.y = (float)grid_cell->center.y / scale - offsety_over_scale;
            voxel.position_draw.z = (float)grid_cell->center.z / scale - offsetz_over_scale;

            /*
             * No slicing for simulation visualization for now
            if(gui_state->slicing) {
                Vector3 test = Vector3Subtract(voxel.position_draw, p);
                float side = Vector3DotProduct(test, n);

                if(side < 0) {
                    continue;
                }
            }
            */
            voxel.size.x = (float)grid_cell->discretization.x / scale;
            voxel.size.y = (float)grid_cell->discretization.y / scale;
            voxel.size.z = (float)grid_cell->discretization.z / scale;

            voxel.position_mesh = (Vector3){(float)grid_cell->center.x, (float)grid_cell->center.y, (float)grid_cell->center.z};
            voxel.v = (float)grid_cell->v;
            voxel.matrix_position = grid_cell->grid_position;

            draw_context->translations[count] = MatrixTranslate(voxel.position_draw.x, voxel.position_draw.y, voxel.position_draw.z);
            draw_context->translations[count] = MatrixMultiply(MatrixScale(voxel.size.x, voxel.size.y, voxel.size.z), draw_context->translations[count]);

            draw_context->colors[count] = get_color((voxel.v - min_v) / (max_v - min_v), gui_state->voxel_alpha, gui_state->current_scale);

            voxel.draw_index = count;

            bool searching = gui_state->found_volume.position_mesh.x >= 0
                         && gui_state->found_volume.position_mesh.y >= 0
                         && gui_state->found_volume.position_mesh.z >= 0;

            if(gui_state->double_clicked || gui_config->adaptive || searching) {
                collision |= check_volume_selection(&voxel, gui_state);
            }

            if(gui_state->ctrl_pressed) {
                check_mouse_over_volume(&voxel, gui_state);
            }

            trace_ap(gui_state, &voxel, time);

            count++;
        }

        update_selected(collision, gui_state, gui_config, draw_context->colors);
        DrawMeshInstancedWithColors(draw_context, grid_mask, count);
    }
}

void draw_alg_purkinje_network(struct gui_shared_info *gui_config, struct gui_state *gui_state) {

    if(gui_config->grid_info.alg_grid == NULL || gui_config->grid_info.alg_grid->purkinje == NULL)
        return;

    struct grid_purkinje *grid_to_draw = gui_config->grid_info.alg_grid->purkinje;

    uint32_t n_active = grid_to_draw->num_active_purkinje_cells;
    struct cell_node **ac = grid_to_draw->purkinje_cells;

    if(ac == NULL)
        return;

    float min_v = gui_config->min_v;
    float max_v = gui_config->max_v;

    float scale = gui_state->mesh_scale_factor;
    Vector3 mesh_offset = gui_state->mesh_offset;

    for(uint32_t i = 0; i < n_active - 1; i++) {

        Vector3 start_pos;
        start_pos.x = ((float)ac[i]->center.x - mesh_offset.x) / scale;
        start_pos.y = ((float)ac[i]->center.y - mesh_offset.y) / scale;
        start_pos.z = ((float)ac[i]->center.z - mesh_offset.z) / scale;

        Vector3 end_pos;
        end_pos.x = ((float)ac[i + 1]->center.x - mesh_offset.x) / scale;
        end_pos.y = ((float)ac[i + 1]->center.y - mesh_offset.y) / scale;
        end_pos.z = ((float)ac[i + 1]->center.z - mesh_offset.z) / scale;

        Color c = get_color(((float)ac[i]->v - min_v) / (max_v - min_v), gui_state->voxel_alpha, gui_state->current_scale);
        DrawLine3D(start_pos, end_pos, c);
    }
}

void set_visibility_after_split(struct gui_shared_info *gui_config, struct gui_state *gui_state) {

    struct vtk_unstructured_grid *grid = gui_config->grid_info.vtk_grid;

    if(!grid)
        return;

    int64_t *cells = grid->cells;
    if(!cells)
        return;

    point3d_array points = grid->points;
    if(!points)
        return;

    uint32_t n_active = grid->num_cells;
    uint32_t num_points = grid->points_per_cell;

    Vector3 cube_position;
    Vector3 n = gui_state->plane_normal;
    Vector3 p = gui_state->plane_point;

    if(gui_state->old_cell_visibility == NULL) {
        int n_vis = arrlen(grid->cell_visibility);
        arrsetlen(gui_state->old_cell_visibility, n_vis);
        memcpy(gui_state->old_cell_visibility, grid->cell_visibility, n_vis * sizeof(uint8_t));
    }

    struct cell_hash_entry *cells_hash = NULL;

    for(uint32_t i = 0; i < n_active * num_points; i += num_points) {
        struct point_3d mesh_center;
        struct point_3d discretization;
        const struct point_3d cell_point = points[cells[i]];

        discretization.x = (float)fabs((cell_point.x - points[cells[i + 1]].x));
        discretization.y = (float)fabs((cell_point.y - points[cells[i + 3]].y));
        discretization.z = (float)fabs((cell_point.z - points[cells[i + 4]].z));

        mesh_center.x = (float)cell_point.x + discretization.x / 2.0f;
        mesh_center.y = (float)cell_point.y + discretization.y / 2.0f;
        mesh_center.z = (float)cell_point.z + discretization.z / 2.0f;

        cube_position.x = ((float)mesh_center.x - gui_state->mesh_offset.x) / gui_state->mesh_scale_factor;
        cube_position.y = ((float)mesh_center.y - gui_state->mesh_offset.y) / gui_state->mesh_scale_factor;
        cube_position.z = ((float)mesh_center.z - gui_state->mesh_offset.z) / gui_state->mesh_scale_factor;

        // Testing against the plane
        Vector3 test = Vector3Subtract(cube_position, p);
        float side = Vector3DotProduct(test, n);

        if(side < 0) {
            discretization = (struct point_3d){-1, -1, -1};
        }

        hmput(cells_hash, mesh_center, discretization);
    }

    calc_visibility(&gui_config->grid_info.vtk_grid, cells_hash, n_active);
    hmfree(cells_hash);
}

void reset_grid_visibility(struct gui_shared_info *gui_config, struct gui_state *gui_state) {
    struct vtk_unstructured_grid *grid = gui_config->grid_info.vtk_grid;
    int n_vis = arrlen(gui_state->old_cell_visibility);
    arrsetlen(grid->cell_visibility, n_vis);
    memcpy(grid->cell_visibility, gui_state->old_cell_visibility, n_vis * sizeof(uint8_t));
}

//TODO
void calc_min_max_bounds() {

}
