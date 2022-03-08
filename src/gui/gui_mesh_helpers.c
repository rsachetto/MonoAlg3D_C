#include "gui_mesh_helpers.h"
#include "gui_colors.h"
#include "raylib_ext.h"
#include <float.h>

#define SET_MAX(coord)                                                                                                                                         \
    result.coord = mesh_max.coord;                                                                                                                             \
    if(mesh_max.coord != mesh_min.coord)                                                                                                                       \
        result.coord = (mesh_max.coord - mesh_min.coord) / 2.0f;

#define SET_MESH_MIN_MAX(vec1, vec2, vec3)                                                                                                                     \
    mesh_info->vec1.x = vec2.x + vec3.x / 2.0f;                                                                                                                \
    mesh_info->vec1.y = vec2.y + vec3.y / 2.0f;                                                                                                                \
    mesh_info->vec1.z = vec2.z + vec3.z / 2.0f;

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
    }
    else {
        if(grid_to_draw->purkinje) {
            n_active = grid_to_draw->purkinje->number_of_purkinje_cells;
            ac = grid_to_draw->purkinje->purkinje_cells;
            ITERATE_OVER_CELLS ()
        }

    }

    return get_max_min(mesh_max, mesh_min, mesh_max_d, mesh_min_d, mesh_info);
}

#define SET_CENTER()                                                                                                                                           \
    mesh_center.x = (float)points[cells[i]].x + d.x / 2.0f;                                                                                                    \
    mesh_center.y = (float)points[cells[i]].y + d.y / 2.0f;                                                                                                    \
    mesh_center.z = (float)points[cells[i]].z + d.z / 2.0f;                                                                                                    \
    SET_MAX_MIN(x);                                                                                                                                            \
    SET_MAX_MIN(y);                                                                                                                                            \
    SET_MAX_MIN(z);

Vector3 find_mesh_center_vtk(struct vtk_unstructured_grid *grid_to_draw, struct mesh_info *mesh_info) {
    uint32_t n_active = grid_to_draw->num_cells;

    Vector3 mesh_max = V3_SAME(FLT_MIN);
    Vector3 mesh_min = V3_SAME(FLT_MAX);

    Vector3 mesh_max_d = V3_SAME(FLT_MIN);
    Vector3 mesh_min_d = V3_SAME(FLT_MAX);

    int64_t *cells = grid_to_draw->cells;
    point3d_array points = grid_to_draw->points;

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
    }
    else if(grid_to_draw->purkinje) {

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

static bool check_volume_selection(struct voxel *voxel, struct gui_state *gui_state, real_cpu min_v, real_cpu max_v, real_cpu t) {

    Vector3 p_draw = voxel->position_draw;
    Vector3 p_mesh = voxel->position_mesh;

    Vector3 cube_size = voxel->size;

    RayCollision collision = {0};
    RayCollision collision_mouse_over = {0};

    float csx = cube_size.x;
    float csy = cube_size.y;
    float csz = cube_size.z;

    BoundingBox bb = (BoundingBox){(Vector3){p_draw.x - csx / 2, p_draw.y - csy / 2, p_draw.z - csz / 2},
                                   (Vector3){p_draw.x + csx / 2, p_draw.y + csy / 2, p_draw.z + csz / 2}};

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

    if(gui_state->found_volume.position_mesh.x == p_mesh.x && gui_state->found_volume.position_mesh.y == p_mesh.y &&
       gui_state->found_volume.position_mesh.z == p_mesh.z) {
        gui_state->current_selected_volume = *voxel;
        gui_state->found_volume.position_mesh = (Vector3){-1, -1, -1};
        collision.hit = true;
    }

    collision_mouse_over = GetRayCollisionBox(gui_state->ray_mouse_over, bb);

    if(collision_mouse_over.hit) {
        if(gui_state->ray_mouse_over_hit_distance > collision_mouse_over.distance) {
            gui_state->current_mouse_over_volume.position_draw = (Vector3){p_mesh.x, p_mesh.y, p_mesh.z};
            gui_state->current_mouse_over_volume.matrix_position = voxel->matrix_position;
            gui_state->ray_mouse_over_hit_distance = collision_mouse_over.distance;
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

    Vector3 p_mesh = voxel->position_mesh;

    action_potential_array aps = (struct action_potential *)hmget(gui_state->ap_graph_config->selected_aps, p_mesh);

    struct action_potential ap1;
    ap1.t = t;
    ap1.v = voxel->v;
    size_t aps_len = arrlen(aps);

    if(aps != NULL) {
        if(aps_len == 0 || ap1.t > aps[aps_len - 1].t) {
            arrput(aps, ap1);
            hmput(gui_state->ap_graph_config->selected_aps, p_mesh, aps);
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

void draw_vtk_unstructured_grid(struct gui_shared_info *gui_config, Vector3 mesh_offset, float scale, struct gui_state *gui_state, int grid_mask, Shader shader, Mesh cube) {

    struct vtk_unstructured_grid *grid_to_draw = gui_config->grid_info.vtk_grid;

    if(!grid_to_draw)
        return;

    int64_t *cells = grid_to_draw->cells;
    if(!cells)
        return;

    point3d_array points = grid_to_draw->points;
    if(!points)
        return;

    uint32_t n_active = grid_to_draw->num_cells;

    uint32_t num_points = grid_to_draw->points_per_cell;
    int j = 0;

    struct voxel voxel;

    gui_state->ray_hit_distance = FLT_MAX;
    gui_state->ray_mouse_over_hit_distance = FLT_MAX;

    Matrix *translations = RL_MALLOC(n_active * sizeof(Matrix)); // Locations of instances
    Color *colors = RL_MALLOC(n_active * sizeof(Color));

    int count = 0;

    float min_v = 0;
    float max_v = 1;

    int extra_data_len = arrlen(grid_to_draw->extra_values);

    if(gui_config->final_file_index == 0) {
        if(gui_state->current_data_index == -1) {
            min_v = grid_to_draw->min_v;
            max_v = grid_to_draw->max_v;
        }
        else if(extra_data_len > 0 && gui_state->current_data_index < extra_data_len) {
            min_v = grid_to_draw->min_extra_value[gui_state->current_data_index];
            max_v = grid_to_draw->max_extra_value[gui_state->current_data_index];
        }
    }
    else {
        max_v = gui_config->max_v;
        min_v = gui_config->min_v;
    }

    if(min_v == max_v) {
        max_v = max_v + 0.00001;
    }

    gui_config->min_v = min_v;
    gui_config->max_v = max_v;

    float time = gui_config->time;

    bool collision = false;

    Vector3 n = gui_state->plane_normal;
    Vector3 p = gui_state->plane_point;

    for(uint32_t i = 0; i < n_active * num_points; i += num_points) {

        if(grid_mask != 2 && !gui_state->sliced) {
            if(grid_to_draw->cell_visibility && !grid_to_draw->cell_visibility[j]) {
                j += 1;
                continue;
            }
        }

        float mesh_center_x, mesh_center_y, mesh_center_z;
        float dx, dy, dz;

        dx = (float)fabs((points[cells[i]].x - points[cells[i + 1]].x));
        dy = (float)fabs((points[cells[i]].y - points[cells[i + 3]].y));
        dz = (float)fabs((points[cells[i]].z - points[cells[i + 4]].z));

        mesh_center_x = (float)points[cells[i]].x + dx / 2.0f;
        mesh_center_y = (float)points[cells[i]].y + dy / 2.0f;
        mesh_center_z = (float)points[cells[i]].z + dz / 2.0f;

        voxel.position_draw.x = (mesh_center_x - mesh_offset.x) / scale;
        voxel.position_draw.y = (mesh_center_y - mesh_offset.y) / scale;
        voxel.position_draw.z = (mesh_center_z - mesh_offset.z) / scale;

        //TODO: maybe we do not need to check the plain
        if(gui_state->slicing || gui_state->sliced) {
            Vector3 test = Vector3Subtract(voxel.position_draw, p);
            float side = Vector3DotProduct(test, n);

            if(side < 0) {
                j += 1;
                continue;
            }
        }

        Color c = BLUE;
        voxel.v = 0.0;

        if(gui_state->current_data_index == -1) {
            if(grid_to_draw->values) {
                voxel.v = grid_to_draw->values[j];
            }
        }
        else if(extra_data_len > 0 && gui_state->current_data_index < extra_data_len) {
            voxel.v = grid_to_draw->extra_values[gui_state->current_data_index][j];
        }

        c = get_color((voxel.v - min_v) / (max_v - min_v), gui_state->voxel_alpha, gui_state->current_scale);

        voxel.size.x = dx / scale;
        voxel.size.y = dy / scale;
        voxel.size.z = dz / scale;
        voxel.position_mesh = (Vector3){mesh_center_x, mesh_center_y, mesh_center_z};

        translations[count] = MatrixTranslate(voxel.position_draw.x, voxel.position_draw.y, voxel.position_draw.z);
        translations[count] = MatrixMultiply(MatrixScale(voxel.size.x, voxel.size.y, voxel.size.z), translations[count]);

        colors[count] = c;
        voxel.draw_index = count;

        collision |= check_volume_selection(&voxel, gui_state, min_v, max_v, time);
        trace_ap(gui_state, &voxel, time);

        count++;

        j += 1;
    }

    update_selected(collision, gui_state, gui_config, colors);

    DrawMeshInstancedWithColors(cube, shader, colors, translations, grid_mask, count);
    free(translations);
    free(colors);
}

void draw_vtk_purkinje_network(struct gui_shared_info *gui_config, Vector3 mesh_offset, float scale, struct gui_state *gui_state, int grid_mask) {

    if(gui_config->grid_info.vtk_grid == NULL || gui_config->grid_info.vtk_grid->purkinje == NULL) return;

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

    for(uint32_t i = 0; i < n_active * num_points; i += num_points) {

        Vector3 start_pos;
        start_pos.x = (points[cells[i]].x - mesh_offset.x)/scale;
        start_pos.y = (points[cells[i]].y - mesh_offset.y)/scale;
        start_pos.z = (points[cells[i]].z - mesh_offset.z)/scale;

        Vector3 end_pos;
        end_pos.x = (points[cells[i+1]].x - mesh_offset.x)/scale;
        end_pos.y = (points[cells[i+1]].y - mesh_offset.y)/scale;
        end_pos.z = (points[cells[i+1]].z - mesh_offset.z)/scale;

        float v = 0.0;
        Color c = BLUE;

        if(grid_to_draw->values) {
            v = grid_to_draw->values[j];
            c = get_color((v - min_v) / (max_v - min_v), gui_state->voxel_alpha, gui_state->current_scale);
        }

        j++;

        DrawLine3D(start_pos, end_pos, c);
    }
}

void draw_alg_mesh(struct gui_shared_info *gui_config, Vector3 mesh_offset, float scale, struct gui_state *gui_state, int grid_mask, Shader shader, Mesh cube) {

    struct grid *grid_to_draw = gui_config->grid_info.alg_grid;

    if(!grid_to_draw)
        return;

    struct voxel voxel;

    uint32_t n_active = grid_to_draw->num_active_cells;
    struct cell_node **ac = grid_to_draw->active_cells;

    float offsetx_over_scale = mesh_offset.x / scale;
    float offsety_over_scale = mesh_offset.y / scale;
    float offsetz_over_scale = mesh_offset.z / scale;

    float min_v = gui_config->min_v;
    float max_v = gui_config->max_v;
    float time = gui_config->time;

    gui_state->ray_hit_distance = FLT_MAX;
    gui_state->ray_mouse_over_hit_distance = FLT_MAX;

    Matrix *translations = RL_MALLOC(n_active * sizeof(Matrix)); // Locations of instances
    Color *colors = RL_MALLOC(n_active * sizeof(Color));

    if(ac) {

        int count = 0;

        bool collision = false;

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

            voxel.size.x = (float)grid_cell->discretization.x / scale;
            voxel.size.y = (float)grid_cell->discretization.y / scale;
            voxel.size.z = (float)grid_cell->discretization.z / scale;

            voxel.position_mesh = (Vector3){(float)grid_cell->center.x, (float)grid_cell->center.y, (float)grid_cell->center.z};
            voxel.v = (float)grid_cell->v;
            voxel.matrix_position = grid_cell->grid_position;

            translations[count] = MatrixTranslate(voxel.position_draw.x, voxel.position_draw.y, voxel.position_draw.z);
            translations[count] = MatrixMultiply(MatrixScale(voxel.size.x, voxel.size.y, voxel.size.z), translations[count]);

            colors[count] = get_color((voxel.v - min_v) / (max_v - min_v), gui_state->voxel_alpha, gui_state->current_scale);

            voxel.draw_index = count;
            collision |= check_volume_selection(&voxel, gui_state, min_v, max_v, time);
            trace_ap(gui_state, &voxel, time);

            count++;
        }

        update_selected(collision, gui_state, gui_config, colors);
        DrawMeshInstancedWithColors(cube, shader, colors, translations, grid_mask, count);

        free(translations);
        free(colors);
    }
}

void draw_alg_purkinje_network(struct gui_shared_info *gui_config, Vector3 mesh_offset, float scale, struct gui_state *gui_state, int grid_mask) {

    if(gui_config->grid_info.alg_grid == NULL || gui_config->grid_info.alg_grid->purkinje == NULL) return;

    struct grid_purkinje* grid_to_draw = gui_config->grid_info.alg_grid->purkinje;

    uint32_t n_active = grid_to_draw->num_active_purkinje_cells;
    struct cell_node **ac = grid_to_draw->purkinje_cells;

    if(ac == NULL) return;

    float min_v = gui_config->min_v;
    float max_v = gui_config->max_v;


    for(uint32_t i = 0; i < n_active-1; i++) {

        Vector3 start_pos;
        start_pos.x = (ac[i]->center.x - mesh_offset.x)/scale;
        start_pos.y = (ac[i]->center.y - mesh_offset.y)/scale;
        start_pos.z = (ac[i]->center.z - mesh_offset.z)/scale;

        Vector3 end_pos;
        end_pos.x = (ac[i+1]->center.x - mesh_offset.x)/scale;
        end_pos.y = (ac[i+1]->center.y - mesh_offset.y)/scale;
        end_pos.z = (ac[i+1]->center.z - mesh_offset.z)/scale;

        Color c = get_color((ac[i]->v - min_v) / (max_v - min_v), gui_state->voxel_alpha, gui_state->current_scale);
        DrawLine3D(start_pos, end_pos, c);
    }

}
