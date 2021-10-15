//
// Created by sachetto on 11/11/17.
//

#include "gui.h"
#include <float.h>
#include <string.h>

#include "../3dparty/stb_ds.h"
#include "../3dparty/tinyfiledialogs/tinyfiledialogs.h"
#include "../utils/file_utils.h"
#include "color_maps.h"

#include "../3dparty/raylib/src/camera.h"

#define RAYGUI_IMPLEMENTATION
#define RAYGUI_SUPPORT_RICONS
#include "../3dparty/raylib/src/extras/raygui.h"

#include "../3dparty/raylib/src/external/glad.h"
#include "../3dparty/raylib/src/rlgl.h"

#define GUI_TEXTBOX_EXTENDED_IMPLEMENTATION
#include "../3dparty/raylib/src/extras/gui_textbox_extended.h"

#include "raylib_ext.h"
#undef RAYGUI_IMPLEMENTATION

static void set_camera_params(Camera3D *camera, bool set_mode) {
    camera->position = (Vector3){0.1f, 0.1f, 12.0f};
    camera->target = (Vector3){0.f, 0.f, 0.f};
    camera->up = (Vector3){0.0f, 1.0f, 0.0f};
    camera->fovy = 45.0f;

    if(set_mode) {
        camera->projection = CAMERA_PERSPECTIVE;
        SetCameraMode(*camera, CAMERA_FREE);
    }
}

static struct gui_state *new_gui_state_with_font_sizes(float font_size_small, float font_size_big, float ui_scale) {

    MaximizeWindow();

    struct gui_state *gui_state = CALLOC_ONE_TYPE(struct gui_state);
    gui_state->current_selected_volumes = NULL;

    gui_state->ui_scale = ui_scale;

    set_camera_params(&(gui_state->camera), true);

    gui_state->font = LoadFont("res/Roboto.ttf");

    gui_state->font_size_small = font_size_small * ui_scale;
    gui_state->font_size_big = font_size_big * ui_scale;

    gui_state->font_spacing_big = gui_state->font_size_big / (float)gui_state->font.baseSize;
    gui_state->font_spacing_small = gui_state->font_size_small / (float)gui_state->font.baseSize;

    gui_state->handle_keyboard_input = true;

    gui_state->show_ap = true;
    gui_state->scale.window.show = true;

    gui_state->help_box.window.show = false;

    gui_state->end_info_box.window.show = true;
    gui_state->mesh_info_box.window.show = true;
    gui_state->ray.position = (Vector3){FLT_MAX, FLT_MAX, FLT_MAX};
    gui_state->ray.direction = (Vector3){FLT_MAX, FLT_MAX, FLT_MAX};
    gui_state->current_mouse_over_volume.position_draw = (Vector3){-1, -1, -1};
    gui_state->found_volume.position_mesh = (Vector3){-1, -1, -1};
    gui_state->mouse_pos = (Vector2){0, 0};
    gui_state->current_scale = 0;
    gui_state->voxel_alpha = 255;
    gui_state->scale_alpha = 255;
    gui_state->mouse_timer = -1;
    gui_state->search_window.move = false;

    gui_state->ap_graph_config = MALLOC_ONE_TYPE(struct ap_graph_config);
    gui_state->ap_graph_config->selected_ap_point = (Vector2){FLT_MAX, FLT_MAX};
    gui_state->ap_graph_config->selected_point_for_apd1 = (Vector2){FLT_MAX, FLT_MAX};
    gui_state->ap_graph_config->selected_point_for_apd2 = (Vector2){FLT_MAX, FLT_MAX};
    gui_state->ap_graph_config->selected_aps = NULL;
    gui_state->ap_graph_config->graph.drag = false;
    gui_state->ap_graph_config->graph.move = false;

    gui_state->current_window_width = GetScreenWidth();
    gui_state->current_window_height = GetScreenHeight();

    gui_state->ap_graph_config->graph.bounds.height = 300.0f * ui_scale;
    gui_state->ap_graph_config->graph.bounds.width = 690.0f * ui_scale;

    gui_state->ap_graph_config->graph.bounds.x = 10;
    gui_state->ap_graph_config->graph.bounds.y = (float)gui_state->current_window_height - gui_state->ap_graph_config->graph.bounds.height - 90;

    hmdefault(gui_state->ap_graph_config->selected_aps, NULL);

    gui_state->search_window.bounds.width = 220;
    gui_state->search_window.bounds.height = 100;

    gui_state->scale.window.bounds.x = (float)gui_state->current_window_width - 30.0f * ui_scale;
    gui_state->scale.window.bounds.y = (float)gui_state->current_window_height / 1.5f;

    gui_state->scale.window.bounds.width = 20;
    gui_state->scale.window.bounds.height = 0;
    gui_state->scale.calc_bounds = true;

    gui_state->controls_window.show = false;

    gui_state->controls_window.bounds.width = 5.0*32.0 + 5*4;
    gui_state->controls_window.bounds.height = 38.0f + WINDOW_STATUSBAR_HEIGHT;

    gui_state->show_coordinates = true;
    gui_state->double_clicked = false;

    return gui_state;
}

static struct mesh_info *new_mesh_info() {
    struct mesh_info *mesh_info = (struct mesh_info *)malloc(sizeof(struct mesh_info));
    mesh_info->center_calculated = false;
    return mesh_info;
}

static Color get_color(real_cpu value, int alpha, int current_scale) {

    int idx1;
    int idx2;
    real_cpu fract_between = 0;

    if(value <= 0) {
        idx1 = idx2 = 0;
    } else if(value >= 1) {
        idx1 = idx2 = NUM_COLORS - 1;
    } else {
        value = value * (NUM_COLORS - 1);
        idx1 = (int)floor(value); // Our desired color will be after this index.
        idx2 = idx1 + 1;          // ... and before this index (inclusive).
        fract_between = value - (real_cpu)idx1;
    }

    real_cpu color_idx1_0 = color_scales[current_scale][idx1][0];
    real_cpu color_idx1_1 = color_scales[current_scale][idx1][1];
    real_cpu color_idx1_2 = color_scales[current_scale][idx1][2];

    real_cpu color_idx2_0 = color_scales[current_scale][idx2][0];
    real_cpu color_idx2_1 = color_scales[current_scale][idx2][1];
    real_cpu color_idx2_2 = color_scales[current_scale][idx2][2];

    Color result;

    result.r = (unsigned char)(((color_idx2_0 - color_idx1_0) * fract_between + color_idx1_0) * 255);
    result.g = (unsigned char)(((color_idx2_1 - color_idx1_1) * fract_between + color_idx1_1) * 255);
    result.b = (unsigned char)(((color_idx2_2 - color_idx1_2) * fract_between + color_idx1_2) * 255);
    result.a = alpha;

    return result;
}

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

static Vector3 find_mesh_center(struct grid *grid_to_draw, struct mesh_info *mesh_info) {

    uint32_t n_active = grid_to_draw->num_active_cells;
    struct cell_node **ac = grid_to_draw->active_cells;
    struct cell_node *grid_cell;

    Vector3 mesh_max = V3_SAME(FLT_MIN);
    Vector3 mesh_min = V3_SAME(FLT_MAX);

    Vector3 mesh_max_d = V3_SAME(FLT_MIN);
    Vector3 mesh_min_d = V3_SAME(FLT_MAX);

    if(ac) {

        Vector3 d;
        Vector3 mesh_center;

        for(uint32_t i = 0; i < n_active; i++) {
            grid_cell = ac[i];

            d.x = grid_cell->discretization.x;
            d.y = grid_cell->discretization.y;
            d.z = grid_cell->discretization.z;

            mesh_center.x = grid_cell->center.x;
            mesh_center.y = grid_cell->center.y;
            mesh_center.z = grid_cell->center.z;

            SET_MAX_MIN(x);
            SET_MAX_MIN(y);
            SET_MAX_MIN(z);
        }
    }

    return get_max_min(mesh_max, mesh_min, mesh_max_d, mesh_min_d, mesh_info);
}

static Vector3 find_mesh_center_vtk(struct vtk_unstructured_grid *grid_to_draw, struct mesh_info *mesh_info) {
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

    for(uint32_t i = 0; i < n_active * num_points; i += num_points) {

        d.x = (float)fabs((points[cells[i]].x - points[cells[i + 1]].x));
        d.y = (float)fabs((points[cells[i]].y - points[cells[i + 3]].y));
        d.z = (float)fabs((points[cells[i]].z - points[cells[i + 4]].z));

        mesh_center.x = (float)points[cells[i]].x + d.x / 2.0f;
        mesh_center.y = (float)points[cells[i]].y + d.y / 2.0f;
        mesh_center.z = (float)points[cells[i]].z + d.z / 2.0f;

        SET_MAX_MIN(x);
        SET_MAX_MIN(y);
        SET_MAX_MIN(z);
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

static void draw_vtk_unstructured_grid(struct gui_shared_info *gui_config, Vector3 mesh_offset, float scale, struct gui_state *gui_state, Shader shader,
                                       Mesh cube, int grid_mask) {

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

    float min_v = gui_config->min_v;
    float max_v = gui_config->max_v;
    float time = gui_config->time;

    bool collision = false;

    for(uint32_t i = 0; i < n_active * num_points; i += num_points) {

        if(grid_to_draw->cell_visibility && !grid_to_draw->cell_visibility[j]) {
            j += 1;
            continue;
        }

        float mesh_center_x, mesh_center_y, mesh_center_z;
        float dx, dy, dz;

        dx = (float)fabs((points[cells[i]].x - points[cells[i + 1]].x));
        dy = (float)fabs((points[cells[i]].y - points[cells[i + 3]].y));
        dz = (float)fabs((points[cells[i]].z - points[cells[i + 4]].z));

        mesh_center_x = (float)points[cells[i]].x + dx / 2.0f;
        mesh_center_y = (float)points[cells[i]].y + dy / 2.0f;
        mesh_center_z = (float)points[cells[i]].z + dz / 2.0f;

        voxel.v = grid_to_draw->values[j];

        voxel.position_draw.x = (mesh_center_x - mesh_offset.x) / scale;
        voxel.position_draw.y = (mesh_center_y - mesh_offset.y) / scale;
        voxel.position_draw.z = (mesh_center_z - mesh_offset.z) / scale;

        voxel.size.x = dx / scale;
        voxel.size.y = dy / scale;
        voxel.size.z = dz / scale;
        voxel.position_mesh = (Vector3){mesh_center_x, mesh_center_y, mesh_center_z};

        translations[count] = MatrixTranslate(voxel.position_draw.x, voxel.position_draw.y, voxel.position_draw.z);
        translations[count] = MatrixMultiply(MatrixScale(voxel.size.x, voxel.size.y, voxel.size.z), translations[count]);

        colors[count] = get_color((voxel.v - min_v) / (max_v - min_v), gui_state->voxel_alpha, gui_state->current_scale);

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

static void draw_alg_mesh(struct gui_shared_info *gui_config, Vector3 mesh_offset, float scale, struct gui_state *gui_state, Shader shader, Mesh cube,
                          int grid_mask) {

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

            if(!ac[i]->visible) {
                continue;
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

static void draw_ap_graph(struct gui_state *gui_state, struct gui_shared_info *gui_config) {

    int n = hmlen(gui_state->ap_graph_config->selected_aps);

    const float graph_w = gui_state->ap_graph_config->graph.bounds.width;
    const float graph_h = gui_state->ap_graph_config->graph.bounds.height;

    if(gui_state->ap_graph_config->graph.bounds.x + graph_w > (float)gui_state->current_window_width || gui_state->ap_graph_config->graph.bounds.x < 0) {
        gui_state->ap_graph_config->graph.bounds.x = 10;
    }

    if(gui_state->ap_graph_config->graph.bounds.y + graph_h > (float)gui_state->current_window_height) {
        gui_state->ap_graph_config->graph.bounds.y = (float)gui_state->current_window_height - graph_h - 90.0f;
    }

    if(gui_state->ap_graph_config->graph.bounds.y < 0) {
        gui_state->ap_graph_config->graph.bounds.y = 0;
    }

    const float graph_x = gui_state->ap_graph_config->graph.bounds.x;
    const float graph_y = gui_state->ap_graph_config->graph.bounds.y;

    static const Color colors[] = {DARKGRAY, GOLD,     ORANGE, PINK,   RED,        MAROON, GREEN,     LIME,  DARKGREEN,
                                   BLUE,     DARKBLUE, PURPLE, VIOLET, DARKPURPLE, BROWN,  DARKBROWN, BLACK, MAGENTA};

    int num_colors = SIZEOF(colors);

    Font font = gui_state->font;

    float font_size_big = gui_state->font_size_big - 6;
    float font_size_small = gui_state->font_size_small - 6;

    float spacing_big = font_size_big / (float)font.baseSize;
    float spacing_small = font_size_small / (float)font.baseSize;

    Rectangle graph_window = gui_state->ap_graph_config->graph.bounds;
    graph_window.y -= WINDOW_STATUSBAR_HEIGHT;
    graph_window.height += WINDOW_STATUSBAR_HEIGHT;
    gui_state->show_ap = !GuiWindowBox(graph_window, TextFormat("Selected APs (%i)", n));

    gui_state->ap_graph_config->drag_graph_button = (Rectangle){(graph_x + graph_w) - 16.0f, graph_y + graph_h - 16.0f, 16.0f, 16.0f};

    GuiButton(gui_state->ap_graph_config->drag_graph_button, "#71#");

    char tmp[TMP_SIZE];

    snprintf(tmp, TMP_SIZE, "%.2lf", gui_config->final_time);
    Vector2 text_width = MeasureTextEx(font, tmp, font_size_small, spacing_small);

    snprintf(tmp, TMP_SIZE, "%.2lf", gui_config->min_v);
    Vector2 text_width_y = MeasureTextEx(font, tmp, font_size_small, spacing_small);

    gui_state->ap_graph_config->min_x = graph_x + 2.6f * text_width_y.x;
    gui_state->ap_graph_config->max_x = graph_x + graph_w - text_width.x;

    gui_state->ap_graph_config->min_y = graph_y + graph_h - text_width.y;
    gui_state->ap_graph_config->max_y = graph_y + 20.0f; // This is actually the smallest allowed y

    Vector2 p1, p2;

    uint32_t num_ticks;
    float tick_ofsset = 10;
    num_ticks = (int)(gui_config->final_time / tick_ofsset);

    if(num_ticks < MIN_HORIZONTAL_TICKS) {
        num_ticks = MIN_HORIZONTAL_TICKS;
        tick_ofsset = gui_config->final_time / (float)num_ticks;
    } else if(num_ticks > MAX_HORIZONTAL_TICKS) {
        num_ticks = MAX_HORIZONTAL_TICKS;
        tick_ofsset = gui_config->final_time / (float)num_ticks;
    }

    char *time_text;
    char *time_template;
    bool steps = false;

    if(gui_config->dt == 0) {
        time_text = "Time (steps)";
        time_template = "%d";
        steps = true;
    } else {
        time_text = "Time (ms)";
        time_template = "%.2lf";
    }

    // Draw x label
    text_width = MeasureTextEx(font, time_text, font_size_big, spacing_big);

    gui_state->ap_graph_config->min_y -= text_width.y * 1.5f;

    float min_x = gui_state->ap_graph_config->min_x;
    float min_y = gui_state->ap_graph_config->min_y;

    float max_x = gui_state->ap_graph_config->max_x;
    float max_y = gui_state->ap_graph_config->max_y;

    Vector2 text_position = (Vector2){graph_x + gui_state->ap_graph_config->graph.bounds.width / 2.0f - text_width.x / 2.0f, min_y + text_width.y};

    DrawTextEx(font, time_text, text_position, font_size_big, spacing_big, BLACK);

    float time = 0.0f;

    // Draw horizontal ticks (t)
    for(uint32_t t = 0; t <= num_ticks; t++) {

        p1.x = Remap(time, 0.0f, gui_config->final_time, min_x, max_x);
        p1.y = min_y - 5;

        p2.x = p1.x;
        p2.y = min_y + 5;

        if(!(t % 2)) {

            float remaped_data = Remap(p1.x, min_x, max_x, 0.0f, gui_config->final_time);

            if(steps) {
                snprintf(tmp, TMP_SIZE, time_template, (int)remaped_data);
            } else {
                snprintf(tmp, TMP_SIZE, time_template, remaped_data);
            }

            text_width = MeasureTextEx(font, tmp, font_size_small, spacing_small);
            DrawTextEx(font, tmp, (Vector2){p1.x - text_width.x / 2.0f, p1.y + 10}, font_size_small, spacing_small, RED);
        }

        DrawLineV(p1, p2, RED);

        DrawLineV(p1, (Vector2){p1.x, gui_state->ap_graph_config->max_y}, LIGHTGRAY);
        time += tick_ofsset;
    }

    tick_ofsset = 10;
    num_ticks = (uint32_t)((gui_config->max_v - gui_config->min_v) / tick_ofsset);

    if(num_ticks < MIN_VERTICAL_TICKS) {
        num_ticks = MIN_VERTICAL_TICKS;
        tick_ofsset = (gui_config->max_v - gui_config->min_v) / (float)num_ticks;
    } else if(num_ticks > MAX_VERTICAL_TICKS) {
        num_ticks = MAX_VERTICAL_TICKS;
        tick_ofsset = (gui_config->max_v - gui_config->min_v) / (float)num_ticks;
    }

    // Draw y label
    {
        char *label;
        label = "Vm (mV)";
        text_width = MeasureTextEx(font, label, font_size_big, spacing_big);

        text_position.x = gui_state->ap_graph_config->graph.bounds.x + 15;
        text_position.y = min_y - ((min_y - max_y) / 2.0f) + (text_width.x / 2.0f);

        DrawTextPro(font, label, text_position, (Vector2){0, 0}, -90, font_size_big, spacing_big, BLACK);
    }

    float v = (float)gui_config->min_v;
    snprintf(tmp, TMP_SIZE, "%.2lf", v);
    Vector2 max_w = MeasureTextEx(font, tmp, font_size_small, spacing_small);

    // Draw vertical ticks (Vm)
    for(uint32_t t = 0; t <= num_ticks; t++) {

        p1.x = (float)graph_x + 5.0f;
        p1.y = Remap(v, gui_config->min_v, gui_config->max_v, min_y, gui_state->ap_graph_config->max_y);

        snprintf(tmp, TMP_SIZE, "%.2lf", Remap(p1.y, min_y, gui_state->ap_graph_config->max_y, gui_config->min_v, gui_config->max_v));
        text_width = MeasureTextEx(font, tmp, font_size_small, spacing_small);

        DrawTextEx(font, tmp, (Vector2){p1.x + (max_w.x - text_width.x / 2) + 20, p1.y - text_width.y / 2.0f}, font_size_small, spacing_small, RED);

        p1.x = min_x - 5.0f;
        p2.x = p1.x + 10.0f;
        p2.y = p1.y;

        // Grid line
        DrawLineV(p1, (Vector2){gui_state->ap_graph_config->max_x, p1.y}, LIGHTGRAY);

        // Tick line
        DrawLineV(p1, p2, RED);

        v += tick_ofsset;
    }

    // Draw vertical line
    {
        p1.x = min_x;
        p1.y = min_y + 5;

        p2.x = p1.x;
        p2.y = max_y;
        DrawLineV(p1, p2, RED);
    }

    // Draw horizontal line
    {
        p1.x = min_x;
        p1.y = min_y;

        p2.x = max_x;
        p2.y = p1.y;
        DrawLineV(p1, p2, RED);
    }

    struct action_potential *aps;

    // Draw the function
    for(int j = 0; j < n; j++) {

        aps = (struct action_potential *)gui_state->ap_graph_config->selected_aps[j].value;
        int c = arrlen(aps);
        int step = 1;

        if(c > 0) {
            Color line_color = colors[j % num_colors];
            for(int i = 0; i < c; i += step) {

                if(aps[i].t <= gui_config->time) {
                    if(i + step < c) {

                        p1.x = Remap(aps[i].t, 0.0f, gui_config->final_time, min_x, max_x);
                        p1.y = Remap(aps[i].v, gui_config->min_v, gui_config->max_v, min_y, max_y);

                        p2.x = Remap(aps[i + step].t, 0.0f, gui_config->final_time, min_x, max_x);
                        p2.y = Remap(aps[i + step].v, gui_config->min_v, gui_config->max_v, min_y, max_y);

                        if(aps[i + step].v > gui_config->max_v)
                            gui_config->max_v = aps[i + step].v;
                        if(aps[i + step].v < gui_config->min_v)
                            gui_config->min_v = aps[i + step].v;

                        DrawLineEx(p1, p2, 2.0f, line_color);
                    }
                }
            }
        }
    }

    // Draw AP coordinates over mouse cursor
    if(!gui_state->ap_graph_config->graph.drag && gui_state->ap_graph_config->selected_ap_point.x != FLT_MAX &&
       gui_state->ap_graph_config->selected_ap_point.y != FLT_MAX && gui_state->mouse_pos.x < max_x && gui_state->mouse_pos.x > min_x &&
       gui_state->mouse_pos.y < min_y && gui_state->mouse_pos.y > max_y) {

        char *tmp_point = "%.2lf, %.2lf";
        snprintf(tmp, TMP_SIZE, tmp_point, gui_state->ap_graph_config->selected_ap_point.x, gui_state->ap_graph_config->selected_ap_point.y);
        text_width = MeasureTextEx(font, tmp, font_size_small, spacing_small);
        DrawTextEx(font, tmp, (Vector2){gui_state->mouse_pos.x - text_width.x / 2, gui_state->mouse_pos.y - text_width.y}, font_size_small, spacing_small,
                   BLACK);
    }

    if(gui_state->ap_graph_config->selected_point_for_apd1.x != FLT_MAX && gui_state->ap_graph_config->selected_point_for_apd1.y != FLT_MAX) {
        DrawCircleV(gui_state->ap_graph_config->selected_point_for_apd1, 4, RED);
    }

    if(gui_state->ap_graph_config->selected_point_for_apd2.x != FLT_MAX && gui_state->ap_graph_config->selected_point_for_apd2.y != FLT_MAX) {
        DrawCircleV(gui_state->ap_graph_config->selected_point_for_apd2, 4, RED);
        DrawLineV(gui_state->ap_graph_config->selected_point_for_apd1, gui_state->ap_graph_config->selected_point_for_apd2, RED);

        float t1 = Remap(gui_state->ap_graph_config->selected_point_for_apd1.x, min_x, max_x, 0.0f, gui_config->final_time);
        float t2 = Remap(gui_state->ap_graph_config->selected_point_for_apd2.x, min_x, max_x, 0.0f, gui_config->final_time);

        char *tmp_point = "dt = %.4lf";
        snprintf(tmp, TMP_SIZE, tmp_point, fabsf(t2 - t1));
        text_width = MeasureTextEx(font, tmp, font_size_small, spacing_small);

        float x = fminf(gui_state->ap_graph_config->selected_point_for_apd1.x, gui_state->ap_graph_config->selected_point_for_apd2.x);

        DrawTextEx(font, tmp, (Vector2){x + text_width.x / 2.0f, gui_state->ap_graph_config->selected_point_for_apd1.y - text_width.y}, font_size_small,
                   spacing_small, BLACK);
    }
}

static inline void move_rect(Vector2 new_pos, Rectangle *rect) {

    float new_x = new_pos.x;
    float new_y = new_pos.y;

    if(new_x > 1 && new_x + rect->width < (float)GetScreenWidth()) {
        rect->x = new_x;
    }
    if(new_y > 1 && new_y + rect->height < (float)GetScreenHeight()) {
        rect->y = new_y;
    }
}

static inline void drag_box(Vector2 mouse_pos, Rectangle *box) {
    float new_x = mouse_pos.x + box->width / 2;
    float new_y = mouse_pos.y - box->height / 2;
    move_rect((Vector2){new_x, new_y}, box);
}

static void check_window_bounds(Rectangle *box, float current_window_width, float current_window_height) {

    if(box->x + box->width > current_window_width)
        move_rect((Vector2){current_window_width - box->width - 10, box->y}, box);

    if(box->y + box->height > current_window_height)
        move_rect((Vector2){box->x, current_window_height - box->height - 10}, box);
}

static void draw_scale(float min_v, float max_v, struct gui_state *gui_state, bool int_scale) {

    float scale_width = 20 * gui_state->ui_scale;
    check_window_bounds(&(gui_state->scale.window.bounds), (float)gui_state->current_window_width, (float)gui_state->current_window_height);

    float spacing_small = gui_state->font_spacing_small;
    float spacing_big = gui_state->font_spacing_big;

    uint32_t num_ticks;
    float tick_ofsset = 12.0f;

    if(!int_scale) {
        num_ticks = (uint32_t)((max_v - min_v) / tick_ofsset);

        if(num_ticks < 5) {
            num_ticks = 5;
            tick_ofsset = (max_v - min_v) / (float)num_ticks;
        }
    } else {
        num_ticks = 0;
        for(uint32_t i = (uint32_t)min_v; i < (uint32_t)max_v; i++) {
            num_ticks++;
        }
        tick_ofsset = (max_v - min_v) / (float)num_ticks;
    }

    char tmp[TMP_SIZE];

    float v = max_v;
    snprintf(tmp, TMP_SIZE, "%.2lf", v);
    Vector2 max_w = MeasureTextEx(gui_state->font, tmp, gui_state->font_size_small, spacing_small);

    Vector2 p1, p2, width;

    float scale_rec_height = 30.0f * gui_state->ui_scale;
    Color color;

    if(gui_state->scale.calc_bounds) {
        gui_state->scale.window.bounds.height += max_w.y;
        for(int t = 0; t <= num_ticks; t++) {
            gui_state->scale.window.bounds.height += scale_rec_height;
        }

        gui_state->scale.window.bounds.y -= gui_state->scale.window.bounds.height;

        gui_state->scale.calc_bounds = false;
    }

    Vector2 scale_bounds = (Vector2) {(float) gui_state->scale.window.bounds.x, (float)gui_state->scale.window.bounds.y};

    if(!int_scale) {
        width = MeasureTextEx(gui_state->font, "Vm", gui_state->font_size_big, spacing_big);
        float diff = (float)scale_width - width.x;

        p1.x = scale_bounds.x + (diff / 2.0f);
        p1.y = scale_bounds.y - width.y;
        DrawTextEx(gui_state->font, "Vm", p1, gui_state->font_size_big, spacing_big, BLACK);
    }

    float initial_y = scale_bounds.y;

    for(uint32_t t = 0; t <= num_ticks; t++) {
        p1.x = scale_bounds.x - 50.0f * gui_state->ui_scale;
        p1.y = initial_y + (float)scale_rec_height / 2.0f;

        snprintf(tmp, TMP_SIZE, "%.2lf", v);
        width = MeasureTextEx(gui_state->font, tmp, gui_state->font_size_small, spacing_small);

        DrawTextEx(gui_state->font, tmp, (Vector2){p1.x + (max_w.x - width.x), p1.y - width.y / 2.0f}, gui_state->font_size_small, spacing_small, BLACK);

        p1.x = p1.x + max_w.x + 2.5f;
        p2.x = p1.x + 10.0f;
        p2.y = p1.y;

        DrawLineV(p1, p2, BLACK);
        color = get_color((v - min_v) / (max_v - min_v), gui_state->scale_alpha, gui_state->current_scale);

        DrawRectangle((int)scale_bounds.x, (int)initial_y, (int)scale_width, (int)scale_rec_height, color);
        initial_y += scale_rec_height;
        v -= tick_ofsset;
    }
}

static void draw_control_window(struct gui_state *gui_state, struct gui_shared_info *gui_config) {

    check_window_bounds(&gui_state->controls_window.bounds, gui_state->current_window_width, gui_state->current_window_height);
    gui_state->controls_window.show = !GuiWindowBox(gui_state->controls_window.bounds, "Controls");

    Rectangle button_pos = (Rectangle){0, 0, 32, 32};

    button_pos.x = gui_state->controls_window.bounds.x + 2.0;
    button_pos.y = gui_state->controls_window.bounds.y + WINDOW_STATUSBAR_HEIGHT + 3.0;

    GuiButton(button_pos, "#71#");

    button_pos.x += button_pos.width + 4.0;
    GuiButton(button_pos, "#71#");

    button_pos.x += button_pos.width + 4.0;
    GuiButton(button_pos, "#71#");

    button_pos.x += button_pos.width + 4.0;
    GuiButton(button_pos, "#71#");

    button_pos.x += button_pos.width + 4.0;
    GuiButton(button_pos, "#71#");


}

static void draw_box(struct gui_text_window *box, struct gui_state *gui_state, float text_offset) {


    check_window_bounds(&box->window.bounds, (float)gui_state->current_window_width, (float)gui_state->current_window_height);

    float text_x = box->window.bounds.x + 20;
    float text_y = box->window.bounds.y + 10 + WINDOW_STATUSBAR_HEIGHT;

    box->window.show = !GuiWindowBox(box->window.bounds, box->title);

    int num_lines = box->num_lines;

    for(int i = 0; i < num_lines; i++) {
        DrawTextEx(gui_state->font, box->lines[i], (Vector2){text_x, text_y}, gui_state->font_size_small, gui_state->font_spacing_small, BLACK);
        text_y += text_offset;
    }
}

static inline void configure_end_info_box_strings(struct gui_shared_info *gui_config, char ***info_string) {

    char tmp[TMP_SIZE];

    int index = 0;

    snprintf(tmp, TMP_SIZE, "Resolution Time: %ld s", gui_config->solver_time / 1000 / 1000);
    (*(info_string))[index++] = strdup(tmp);

    snprintf(tmp, TMP_SIZE, "ODE Total Time: %ld s", gui_config->ode_total_time / 1000 / 1000);
    (*(info_string))[index++] = strdup(tmp);

    snprintf(tmp, TMP_SIZE, "CG Total Time: %ld s", gui_config->cg_total_time / 1000 / 1000);
    (*(info_string))[index++] = strdup(tmp);

    snprintf(tmp, TMP_SIZE, "Mat time: %ld s", gui_config->total_mat_time / 1000 / 1000);
    (*(info_string))[index++] = strdup(tmp);

    snprintf(tmp, TMP_SIZE, "Refine time: %ld s", gui_config->total_ref_time / 1000 / 1000);
    (*(info_string))[index++] = strdup(tmp);

    snprintf(tmp, TMP_SIZE, "Unrefine time: %ld s", gui_config->total_deref_time / 1000 / 1000);
    (*(info_string))[index++] = strdup(tmp);

    snprintf(tmp, TMP_SIZE, "Write time: %ld s", gui_config->total_write_time / 1000 / 1000);
    (*(info_string))[index++] = strdup(tmp);

    snprintf(tmp, TMP_SIZE, "CG Total Iterations: %ld", gui_config->total_cg_it);
    (*(info_string))[index++] = strdup(tmp);

    snprintf(tmp, TMP_SIZE, "Final Time: %.3lf ms", gui_config->time);
    (*(info_string))[index] = strdup(tmp);
}

static inline bool configure_mesh_info_box_strings(struct gui_state *gui_state, struct gui_shared_info *gui_config, char ***info_string, int draw_type,
                                                   struct mesh_info *mesh_info) {

    if(!gui_config->grid_info.alg_grid && !gui_config->grid_info.vtk_grid)
        return false;

    char tmp[TMP_SIZE];

    int index = 0;

    uint32_t n_active = 0;

    if(draw_type == DRAW_SIMULATION) {
        n_active = gui_config->grid_info.alg_grid->num_active_cells;
    } else {
        n_active = gui_config->grid_info.vtk_grid->num_cells;
    }

    snprintf(tmp, TMP_SIZE, " - Num. of Volumes: %u", n_active);
    (*(info_string))[index++] = strdup(tmp);

    snprintf(tmp, TMP_SIZE, " - Max X: %f", mesh_info->max_size.x);
    (*(info_string))[index++] = strdup(tmp);

    snprintf(tmp, TMP_SIZE, " - Max Y: %f", mesh_info->max_size.y);
    (*(info_string))[index++] = strdup(tmp);

    snprintf(tmp, TMP_SIZE, " - Max Z: %f", mesh_info->max_size.z);
    (*(info_string))[index++] = strdup(tmp);

    snprintf(tmp, TMP_SIZE, " - Min X: %f", mesh_info->min_size.x);
    (*(info_string))[index++] = strdup(tmp);

    snprintf(tmp, TMP_SIZE, " - Min Y: %f", mesh_info->min_size.y);
    (*(info_string))[index++] = strdup(tmp);

    snprintf(tmp, TMP_SIZE, " - Min Z: %f", mesh_info->min_size.z);
    (*(info_string))[index++] = strdup(tmp);

    if(draw_type == DRAW_SIMULATION) {
        if(gui_config->paused) {
            snprintf(tmp, TMP_SIZE, "Simulation paused: %.3lf of %.3lf ms", gui_config->time, gui_config->final_time);
        } else if(gui_config->simulating) {
            snprintf(tmp, TMP_SIZE, "Simulation running: %.3lf of %.3lf ms", gui_config->time, gui_config->final_time);

        } else {
            snprintf(tmp, TMP_SIZE, "Simulation finished: %.3lf of %.3lf ms", gui_config->time, gui_config->final_time);
        }
    } else {
        if(gui_config->paused) {
            if(gui_config->dt == 0) {
                snprintf(tmp, TMP_SIZE, "Visualization paused: %d of %d steps", (int)gui_config->time, (int)gui_config->final_time);
            } else {
                snprintf(tmp, TMP_SIZE, "Visualization paused: %.3lf of %.3lf ms", gui_config->time, gui_config->final_time);
            }
        } else if(gui_config->simulating) {
            if(gui_config->dt == 0) {
                snprintf(tmp, TMP_SIZE, "Visualization running: %d of %d steps", (int)gui_config->time, (int)gui_config->final_time);
            } else {
                snprintf(tmp, TMP_SIZE, "Visualization running: %.3lf of %.3lf ms", gui_config->time, gui_config->final_time);
            }

        } else {
            if(gui_config->dt == 0) {
                snprintf(tmp, TMP_SIZE, "Visualization finished: %d of %d steps", (int)gui_config->time, (int)gui_config->final_time);
            } else {
                snprintf(tmp, TMP_SIZE, "Visualization finished: %.3lf of %.3lf ms", gui_config->time, gui_config->final_time);
            }
        }
    }

    Vector2 wider_text = MeasureTextEx(gui_state->font, tmp, gui_state->font_size_small, gui_state->font_spacing_small);

    float box_w = wider_text.x * 1.08f;

    gui_state->mesh_info_box.window.bounds.width = box_w;

    (*(info_string))[index] = strdup(tmp);

    return true;
}

static bool draw_selection_box(struct gui_state *gui_state) {

#define CENTER_X "Center X"
#define CENTER_Y "Center Y"
#define CENTER_Z "Center Z"

    Vector2 text_box_size = MeasureTextEx(gui_state->font, CENTER_X, gui_state->font_size_small, gui_state->font_spacing_small);

    float text_box_y_dist = text_box_size.y * 2.5f;
    float label_box_y_dist = 30;
    float x_off = 10;

    static char center_x_text[128] = {0};
    static char center_y_text[128] = {0};
    static char center_z_text[128] = {0};

    float pos_x = gui_state->search_window.bounds.x;
    float pos_y = gui_state->search_window.bounds.y;

    float box_pos = pos_x + x_off;
    gui_state->search_window.bounds.width = text_box_size.x * 3.5f;
    gui_state->search_window.bounds.height = (text_box_size.y + text_box_y_dist) * 1.6f;

    bool window_closed = GuiWindowBox(gui_state->search_window.bounds, "Enter the center of the cell");

    DrawTextEx(gui_state->font, CENTER_X, (Vector2){box_pos, pos_y + label_box_y_dist}, gui_state->font_size_small, gui_state->font_spacing_small, BLACK);

    GuiTextBoxEx((Rectangle){box_pos, pos_y + text_box_y_dist, text_box_size.x, text_box_size.y}, center_x_text, SIZEOF(center_x_text) - 1, true);

    box_pos = pos_x + text_box_size.x + 2 * x_off;
    DrawTextEx(gui_state->font, CENTER_Y, (Vector2){box_pos, pos_y + label_box_y_dist}, gui_state->font_size_small, gui_state->font_spacing_small, BLACK);
    GuiTextBoxEx((Rectangle){box_pos, pos_y + text_box_y_dist, text_box_size.x, text_box_size.y}, center_y_text, SIZEOF(center_y_text) - 1, true);

    box_pos = pos_x + 2 * text_box_size.x + 3 * x_off;
    DrawTextEx(gui_state->font, CENTER_Z, (Vector2){box_pos, pos_y + label_box_y_dist}, gui_state->font_size_small, gui_state->font_spacing_small, BLACK);
    GuiTextBoxEx((Rectangle){box_pos, pos_y + text_box_y_dist, text_box_size.x, text_box_size.y}, center_z_text, SIZEOF(center_z_text) - 1, true);

    bool btn_ok_clicked =
        GuiButton((Rectangle){pos_x + text_box_size.x + 2 * x_off, pos_y + (text_box_size.y + text_box_y_dist) * 1.2f, text_box_size.x, text_box_size.y}, "OK");

    if(btn_ok_clicked) {
        gui_state->found_volume.position_mesh.x = strtof(center_x_text, NULL);
        gui_state->found_volume.position_mesh.y = strtof(center_y_text, NULL);
        gui_state->found_volume.position_mesh.z = strtof(center_z_text, NULL);
    }

    return window_closed || btn_ok_clicked;
}

static inline void reset_ui(struct gui_state *gui_state) {

    gui_state->help_box.window.bounds.x = 10;
    gui_state->help_box.window.bounds.y = 10;

    gui_state->ap_graph_config->graph.bounds.height = 300.0f * gui_state->ui_scale;
    gui_state->ap_graph_config->graph.bounds.width = 690.0f * gui_state->ui_scale;

    gui_state->ap_graph_config->graph.bounds.x = 10;
    gui_state->ap_graph_config->graph.bounds.y = (float)gui_state->current_window_height - gui_state->ap_graph_config->graph.bounds.height - 90;

    gui_state->search_window.bounds.width = 220;
    gui_state->search_window.bounds.height = 100;

    gui_state->mesh_info_box.window.bounds.x = (float)gui_state->current_window_width - gui_state->mesh_info_box.window.bounds.width - 10;
    gui_state->mesh_info_box.window.bounds.y = 10.0f;

    gui_state->end_info_box.window.bounds.x = gui_state->mesh_info_box.window.bounds.x - gui_state->mesh_info_box.window.bounds.width - 10;
    gui_state->end_info_box.window.bounds.y = gui_state->mesh_info_box.window.bounds.y;

    gui_state->scale.window.bounds.x = (float)gui_state->current_window_width - 30.0f * gui_state->ui_scale;
    gui_state->scale.window.bounds.y = (float)gui_state->current_window_height / 1.5f;

    gui_state->scale.window.bounds.width = 20;
    gui_state->scale.window.bounds.height = 0;
    gui_state->scale.calc_bounds = true;

    gui_state->controls_window.bounds.x = (float)gui_state->current_window_width/2.0;
    gui_state->controls_window.bounds.y = 10;

}

static void reset(struct gui_shared_info *gui_config, struct gui_state *gui_state, bool full_reset) {

    if(!gui_config->paused) {
        return;
    }

    gui_state->voxel_alpha = 255;

    for(size_t i = 0; i < hmlen(gui_state->ap_graph_config->selected_aps); i++) {
        arrsetlen(gui_state->ap_graph_config->selected_aps[i].value, 0);
    }

    if(gui_config->paused) {
        omp_unset_lock(&gui_config->sleep_lock);
        gui_config->paused = false;
    }

    gui_config->restart = true;
    gui_config->grid_info.alg_grid = NULL;
    gui_config->grid_info.vtk_grid = NULL;

    gui_state->ray.position = (Vector3){FLT_MAX, FLT_MAX, FLT_MAX};
    gui_state->ray.direction = (Vector3){FLT_MAX, FLT_MAX, FLT_MAX};
    omp_unset_lock(&gui_config->sleep_lock);
    gui_state->current_selected_volume.position_draw = (Vector3){FLT_MAX, FLT_MAX, FLT_MAX};
    gui_state->current_mouse_over_volume.position_draw = (Vector3){-1, -1, -1};
    gui_state->ap_graph_config->selected_point_for_apd1 = (Vector2){FLT_MAX, FLT_MAX};
    gui_state->ap_graph_config->selected_point_for_apd2 = (Vector2){FLT_MAX, FLT_MAX};

    if(full_reset) {

        for(size_t i = 0; i < hmlen(gui_state->ap_graph_config->selected_aps); i++) {
            arrfree(gui_state->ap_graph_config->selected_aps[i].value);
        }

        hmfree(gui_state->ap_graph_config->selected_aps);
        gui_state->ap_graph_config->selected_aps = NULL;
        hmdefault(gui_state->ap_graph_config->selected_aps, NULL);

        hmfree(gui_state->current_selected_volumes);
        gui_state->current_selected_volumes = NULL;

        set_camera_params(&(gui_state->camera), false);

        gui_state->scale_alpha = 255;

        reset_ui(gui_state);

        gui_state->show_coordinates = true;
    }
}

static void handle_keyboard_input(struct gui_shared_info *gui_config, struct gui_state *gui_state) {

    if(gui_config->paused) {

        if(IsKeyDown(KEY_RIGHT_CONTROL) || IsKeyDown((KEY_LEFT_CONTROL))) {
            if(IsKeyPressed(KEY_S)) {
                char const *filter[1] = {"*.vtu"};

                const char *save_path = tinyfd_saveFileDialog("Save VTK file", gui_config->input, 1, filter, "vtu files");

                if(save_path) {
                    save_vtk_unstructured_grid_as_vtu_compressed(gui_config->grid_info.vtk_grid, save_path, 6);
                    log_info("Saved vtk file as %s\n", save_path);
                }

                return;
            }
        }

        if(IsKeyPressed(KEY_RIGHT) || IsKeyDown(KEY_UP)) {
            gui_config->advance_or_return = 1;
            omp_unset_lock(&gui_config->sleep_lock);
            return;
        }

        if(gui_config->draw_type == DRAW_FILE) {
            // Return one step only works on file visualization...
            if(IsKeyPressed(KEY_LEFT) || IsKeyDown(KEY_DOWN)) {
                gui_config->advance_or_return = -1;
                omp_unset_lock(&gui_config->sleep_lock);
                return;
            }
        }
    }

    if(IsKeyDown(KEY_RIGHT_CONTROL) || IsKeyDown((KEY_LEFT_CONTROL))) {
        if(IsKeyPressed(KEY_F)) {
            gui_state->show_selection_box = true;
            gui_state->search_window.bounds.x = (float)GetScreenWidth() / 2.0f - gui_state->search_window.bounds.width;
            gui_state->search_window.bounds.y = (float)GetScreenHeight() / 2.0f - gui_state->search_window.bounds.height;
            return;
        }
    }

    if(IsKeyPressed(KEY_Q)) {
        gui_state->scale.window.show = !gui_state->scale.window.show;
        return;
    }

    if(IsKeyDown(KEY_A)) {
        if(gui_state->voxel_alpha - 1 >= 1) {
            gui_state->voxel_alpha = gui_state->voxel_alpha - 1;
        }
        return;
    }

    if(IsKeyDown(KEY_Z)) {
        if(gui_state->voxel_alpha + 1 <= 255) {
            gui_state->voxel_alpha = gui_state->voxel_alpha + 1;
        }
        return;
    }

    if(IsKeyPressed(KEY_G)) {
        gui_state->draw_grid_only = !gui_state->draw_grid_only;
        return;
    }

    if(IsKeyPressed(KEY_PERIOD)) {
        gui_state->current_scale = (gui_state->current_scale + 1) % NUM_SCALES;
        return;
    }

    if(IsKeyPressed(KEY_COMMA)) {
        if(gui_state->current_scale - 1 >= 0) {
            gui_state->current_scale = (gui_state->current_scale - 1);
        } else {
            gui_state->current_scale = NUM_SCALES - 1;
        }
        return;
    }

    if(IsKeyPressed(KEY_L)) {
        gui_state->draw_grid_lines = !gui_state->draw_grid_lines;
        return;
    }

    if(IsKeyPressed(KEY_SPACE)) {
        gui_config->paused = !gui_config->paused;
        return;
    }

    if(IsKeyPressed(KEY_R)) {
        bool full_reset = false;

        if(IsKeyDown(KEY_LEFT_ALT)) {
            full_reset = true;
        }

        reset(gui_config, gui_state, full_reset);
        return;
    }

    if(IsKeyPressed(KEY_X)) {
        gui_state->show_ap = !gui_state->show_ap;
        return;
    }

    if(IsKeyPressed(KEY_C)) {
        gui_state->scale.window.show = gui_state->c_pressed;
        gui_state->show_ap = gui_state->c_pressed;
        gui_state->end_info_box.window.show = gui_state->c_pressed;
        gui_state->mesh_info_box.window.show = gui_state->c_pressed;
        gui_state->c_pressed = !gui_state->c_pressed;
        gui_state->show_coordinates = !gui_state->show_coordinates;
        return;
    }

    if(IsKeyPressed(KEY_H)) {
        gui_state->help_box.window.show = !gui_state->help_box.window.show;
        return;
    }

    if(gui_config->draw_type == DRAW_FILE) {

        if(IsKeyPressed(KEY_O)) {

            gui_config->paused = true;

            char *buf = get_current_directory();

            char const *tmp = tinyfd_selectFolderDialog("Select a directory", buf);
            if(tmp) {
                gui_config->input = strdup(tmp);
            } else {
                gui_config->input = NULL;
            }

            free(buf);

            if(gui_config->input) {
                reset(gui_config, gui_state, true);
                free(gui_config->error_message);
                gui_config->error_message = strdup("Loading Mesh...");
            }

            return;
        }

        if(IsKeyPressed(KEY_F)) {

            gui_config->paused = true;

            char *buf = get_current_directory();

            char const *filter[4] = {"*.pvd", "*.acm", "*.vtk", "*.vtu"};

            char const *tmp = tinyfd_openFileDialog("Select a simulation file", buf, 4, filter, "simulation result (pvd, vtk, vtu or acm)", 0);

            if(tmp) {
                gui_config->input = strdup(tmp);
            } else {
                gui_config->input = NULL;
            }

            free(buf);

            if(tmp) {
                reset(gui_config, gui_state, true);
                free(gui_config->error_message);
                gui_config->error_message = strdup("Loading Mesh...");
            }
            return;
        }
    }
}

static void handle_input(struct gui_shared_info *gui_config, struct mesh_info *mesh_info, struct gui_state *gui_state) {

    if(gui_state->handle_keyboard_input) {
        handle_keyboard_input(gui_config, gui_state);
    }

    gui_state->mouse_pos = GetMousePosition();
    gui_state->ray_mouse_over = GetMouseRay(GetMousePosition(), gui_state->camera);

    if(IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {

        gui_state->ray = GetMouseRay(GetMousePosition(), gui_state->camera);

        if(!gui_state->show_selection_box) {
            if(gui_state->mouse_timer == -1) {
                gui_state->double_clicked = false;
                gui_state->mouse_timer = GetTime();
            } else {
                double delay = GetTime() - gui_state->mouse_timer;
                if(delay < DOUBLE_CLICK_DELAY) {
                    gui_state->double_clicked = true;
                    gui_state->mouse_timer = -1;
                } else {
                    gui_state->mouse_timer = -1;
                    gui_state->double_clicked = false;
                }
            }
        }

        if(CheckCollisionPointRec(gui_state->mouse_pos,
                                  (Rectangle){gui_state->search_window.bounds.x, gui_state->search_window.bounds.y, gui_state->search_window.bounds.width - 18, WINDOW_STATUSBAR_HEIGHT})) {
            gui_state->search_window.move = true;
        } else if(CheckCollisionPointRec(gui_state->mouse_pos, gui_state->ap_graph_config->drag_graph_button)) {
            gui_state->ap_graph_config->graph.drag = true;
        } else if(CheckCollisionPointRec(gui_state->mouse_pos,
                                         (Rectangle){gui_state->ap_graph_config->graph.bounds.x, gui_state->ap_graph_config->graph.bounds.y - WINDOW_STATUSBAR_HEIGHT,
                                                     gui_state->ap_graph_config->graph.bounds.width - 18, WINDOW_STATUSBAR_HEIGHT})) {
            gui_state->ap_graph_config->graph.move = true;
        } else if(CheckCollisionPointRec(gui_state->mouse_pos, (Rectangle){gui_state->help_box.window.bounds.x, gui_state->help_box.window.bounds.y,
                                                                           gui_state->help_box.window.bounds.width - 18, WINDOW_STATUSBAR_HEIGHT})) {
            gui_state->help_box.window.move = true;
        } else if(CheckCollisionPointRec(gui_state->mouse_pos, (Rectangle){gui_state->mesh_info_box.window.bounds.x, gui_state->mesh_info_box.window.bounds.y,
                                                                           gui_state->mesh_info_box.window.bounds.width - 18, WINDOW_STATUSBAR_HEIGHT})) {
            gui_state->mesh_info_box.window.move = true;
        } else if(CheckCollisionPointRec(gui_state->mouse_pos, (Rectangle){gui_state->end_info_box.window.bounds.x, gui_state->end_info_box.window.bounds.y,
                                                                           gui_state->end_info_box.window.bounds.width - 18, WINDOW_STATUSBAR_HEIGHT})) {
            gui_state->end_info_box.window.move = true;
        } else if(CheckCollisionPointRec(gui_state->mouse_pos, (Rectangle){gui_state->scale.window.bounds.x, gui_state->scale.window.bounds.y, gui_state->scale.window.bounds.width,
                                                                           gui_state->scale.window.bounds.height})) {
            gui_state->scale.window.move = true;
        }

    } else if(IsMouseButtonPressed(MOUSE_RIGHT_BUTTON)) {
        if(hmlen(gui_state->ap_graph_config->selected_aps) && gui_state->show_ap) {
            if(CheckCollisionPointRec(gui_state->mouse_pos, gui_state->ap_graph_config->graph.bounds)) {
                if(gui_state->ap_graph_config->selected_point_for_apd1.x == FLT_MAX && gui_state->ap_graph_config->selected_point_for_apd1.y == FLT_MAX) {
                    gui_state->ap_graph_config->selected_point_for_apd1.x = gui_state->mouse_pos.x;
                    gui_state->ap_graph_config->selected_point_for_apd1.y = gui_state->mouse_pos.y;
                } else {
                    if(gui_state->ap_graph_config->selected_point_for_apd2.x == FLT_MAX && gui_state->ap_graph_config->selected_point_for_apd2.y == FLT_MAX) {
                        gui_state->ap_graph_config->selected_point_for_apd2.x = gui_state->mouse_pos.x;
                        gui_state->ap_graph_config->selected_point_for_apd2.y = gui_state->ap_graph_config->selected_point_for_apd1.y;
                    } else {
                        gui_state->ap_graph_config->selected_point_for_apd1.x = gui_state->mouse_pos.x;
                        gui_state->ap_graph_config->selected_point_for_apd1.y = gui_state->mouse_pos.y;

                        gui_state->ap_graph_config->selected_point_for_apd2.x = FLT_MAX;
                        gui_state->ap_graph_config->selected_point_for_apd2.y = FLT_MAX;
                    }
                }
            }
        }
    }

    if(gui_state->search_window.move) {
        gui_state->search_window.bounds.x = (gui_state->mouse_pos.x) - (gui_state->search_window.bounds.width - 18.0f) / 2.0f;
        gui_state->search_window.bounds.y = (gui_state->mouse_pos.y) - WINDOW_STATUSBAR_HEIGHT / 2.0f;

        if(IsMouseButtonReleased(MOUSE_LEFT_BUTTON)) {
            gui_state->search_window.move = false;
        }
    } else if(gui_state->ap_graph_config->graph.drag) {

        float new_heigth = gui_state->mouse_pos.y - gui_state->ap_graph_config->graph.bounds.y;

        if(new_heigth > 100) {
            gui_state->ap_graph_config->graph.bounds.height = new_heigth;
        }

        float new_width = gui_state->mouse_pos.x - gui_state->ap_graph_config->graph.bounds.x;

        if(new_width > 200) {
            gui_state->ap_graph_config->graph.bounds.width = new_width;
        }

        gui_state->ap_graph_config->selected_point_for_apd1.x = FLT_MAX;
        gui_state->ap_graph_config->selected_point_for_apd1.y = FLT_MAX;
        gui_state->ap_graph_config->selected_point_for_apd2.x = FLT_MAX;
        gui_state->ap_graph_config->selected_point_for_apd2.y = FLT_MAX;

        if(IsMouseButtonReleased(MOUSE_LEFT_BUTTON)) {
            gui_state->ap_graph_config->graph.drag = false;
        }
    } else if(gui_state->ap_graph_config->graph.move) {
        move_rect((Vector2){(gui_state->mouse_pos.x) - (gui_state->ap_graph_config->graph.bounds.width - 18.0f) / 2.0f,
                            (gui_state->mouse_pos.y) + WINDOW_STATUSBAR_HEIGHT / 2.0f},
                  &gui_state->ap_graph_config->graph.bounds);
        if(IsMouseButtonReleased(MOUSE_LEFT_BUTTON))
            gui_state->ap_graph_config->graph.move = false;

        gui_state->ap_graph_config->selected_point_for_apd1.x = FLT_MAX;
        gui_state->ap_graph_config->selected_point_for_apd1.y = FLT_MAX;
        gui_state->ap_graph_config->selected_point_for_apd2.x = FLT_MAX;
        gui_state->ap_graph_config->selected_point_for_apd2.y = FLT_MAX;
    } else if(gui_state->help_box.window.move) {
        move_rect(
            (Vector2){(gui_state->mouse_pos.x) - (gui_state->help_box.window.bounds.width - 18.0f) / 2.0f, (gui_state->mouse_pos.y) - WINDOW_STATUSBAR_HEIGHT / 2.0f},
            &gui_state->help_box.window.bounds);
        if(IsMouseButtonReleased(MOUSE_LEFT_BUTTON))
            gui_state->help_box.window.move = false;
    } else if(gui_state->mesh_info_box.window.move) {
        move_rect((Vector2){(gui_state->mouse_pos.x) - (gui_state->mesh_info_box.window.bounds.width - 18.0f) / 2.0f,
                            (gui_state->mouse_pos.y) - WINDOW_STATUSBAR_HEIGHT / 2.0f},
                  &gui_state->mesh_info_box.window.bounds);
        if(IsMouseButtonReleased(MOUSE_LEFT_BUTTON))
            gui_state->mesh_info_box.window.move = false;

    } else if(gui_state->scale.window.move) {
        drag_box(gui_state->mouse_pos, &gui_state->scale.window.bounds);
        if(IsMouseButtonReleased(MOUSE_LEFT_BUTTON))
            gui_state->scale.window.move = false;
    } else if(gui_state->end_info_box.window.move) {
        move_rect((Vector2){(gui_state->mouse_pos.x) - (gui_state->end_info_box.window.bounds.width - 18.0f) / 2.0f,
                            (gui_state->mouse_pos.y) - WINDOW_STATUSBAR_HEIGHT / 2.0f},
                  &gui_state->end_info_box.window.bounds);
        if(IsMouseButtonReleased(MOUSE_LEFT_BUTTON))
            gui_state->end_info_box.window.move = false;
    }

    if(hmlen(gui_state->ap_graph_config->selected_aps) && gui_state->show_ap) {
        float t = Remap(gui_state->mouse_pos.x, gui_state->ap_graph_config->min_x, gui_state->ap_graph_config->max_x, 0.0f, gui_config->final_time);
        float v = Remap(gui_state->mouse_pos.y, gui_state->ap_graph_config->min_y, gui_state->ap_graph_config->max_y, gui_config->min_v, gui_config->max_v);
        if(CheckCollisionPointRec(gui_state->mouse_pos, gui_state->ap_graph_config->graph.bounds)) {
            gui_state->ap_graph_config->selected_ap_point.x = t;
            gui_state->ap_graph_config->selected_ap_point.y = v;
        } else {
            gui_state->ap_graph_config->selected_ap_point.x = FLT_MAX;
            gui_state->ap_graph_config->selected_ap_point.y = FLT_MAX;
        }
    }
}

static void configure_info_boxes_sizes(struct gui_state *gui_state, int help_box_lines, int mesh_info_box_lines, int end_info_box_lines, float box_w,
                                       float text_offset) {

    float margin = 25.0f;

    gui_state->help_box.window.bounds.width = box_w;
    gui_state->help_box.window.bounds.height = (text_offset * (float)help_box_lines) + margin;
    gui_state->help_box.num_lines = help_box_lines;

    box_w = box_w - 100;
    gui_state->mesh_info_box.window.bounds.width = box_w;
    gui_state->mesh_info_box.window.bounds.height = (text_offset * (float)mesh_info_box_lines) + margin;
    gui_state->mesh_info_box.num_lines = mesh_info_box_lines;

    gui_state->mesh_info_box.window.bounds.x = (float)gui_state->current_window_width - box_w - margin;
    gui_state->mesh_info_box.window.bounds.y = 10.0f;

    gui_state->end_info_box.window.bounds.width = box_w;
    gui_state->end_info_box.window.bounds.height = (text_offset * (float)end_info_box_lines) + margin;
    gui_state->end_info_box.num_lines = end_info_box_lines;

    gui_state->end_info_box.window.bounds.x = gui_state->mesh_info_box.window.bounds.x - gui_state->mesh_info_box.window.bounds.width - margin;
    gui_state->end_info_box.window.bounds.y = gui_state->mesh_info_box.window.bounds.y;
}

static void draw_coordinates(struct gui_state *gui_state) {

    const float line_size = 1.0f;
    const float arrow_offset = 0.1f;
    static bool first_draw = true;

    if(first_draw) {
        gui_state->coordinates_cube = (Vector3){-(line_size / 2.0f) + 0.5f, -5.0f + 0.5f, -5.5f};
        first_draw = false;
    }

    Vector3 start_pos = (Vector3){gui_state->coordinates_cube.x - 0.5f, gui_state->coordinates_cube.y - 0.5f, gui_state->coordinates_cube.z - 0.5f};
    Vector3 end_pos = (Vector3){start_pos.x + line_size, start_pos.y, gui_state->coordinates_cube.z - 0.5f};

    DrawLine3D(start_pos, end_pos, RED);
    DrawLine3D((Vector3){end_pos.x - arrow_offset, end_pos.y + arrow_offset, end_pos.z}, end_pos, RED);
    DrawLine3D((Vector3){end_pos.x - arrow_offset, end_pos.y - arrow_offset, end_pos.z}, end_pos, RED);

    gui_state->coordinates_label_x_position = GetWorldToScreen(end_pos, gui_state->camera);

    end_pos = (Vector3){start_pos.x, start_pos.y + line_size, end_pos.z};

    DrawLine3D(start_pos, end_pos, GREEN);
    DrawLine3D((Vector3){end_pos.x - arrow_offset, end_pos.y - arrow_offset, end_pos.z}, end_pos, GREEN);
    DrawLine3D((Vector3){end_pos.x + arrow_offset, end_pos.y - arrow_offset, end_pos.z}, end_pos, GREEN);

    gui_state->coordinates_label_y_position = GetWorldToScreen((Vector3){end_pos.x, end_pos.y + 0.4f, end_pos.z}, gui_state->camera);

    end_pos = (Vector3){start_pos.x, start_pos.y, start_pos.z + line_size};

    DrawLine3D(start_pos, end_pos, DARKBLUE);
    DrawLine3D((Vector3){end_pos.x - arrow_offset, end_pos.y, end_pos.z - arrow_offset}, end_pos, DARKBLUE);
    DrawLine3D((Vector3){end_pos.x + arrow_offset, end_pos.y, end_pos.z - arrow_offset}, end_pos, DARKBLUE);

    gui_state->coordinates_label_z_position = GetWorldToScreen(end_pos, gui_state->camera);
}

void init_and_open_gui_window(struct gui_shared_info *gui_config) {

    const int end_info_box_lines = 9;
    const int mesh_info_box_lines = 8;

    omp_set_lock(&gui_config->sleep_lock);

    SetConfigFlags(FLAG_WINDOW_RESIZABLE | FLAG_MSAA_4X_HINT);
    char *window_title = NULL;

    int draw_type = gui_config->draw_type;

    if(draw_type == DRAW_SIMULATION) {
        size_t max = strlen(gui_config->config_name) + strlen("Simulation visualization - ") + 2;
        window_title = (char *)malloc(max);
        snprintf(window_title, max + 1, "Simulation visualization - %s", gui_config->config_name);
    } else {
        window_title = strdup("Monoalg3D Simulator");
    }

    InitWindow(0, 0, window_title);

    if(gui_config->ui_scale == 0.0) {
        Vector2 ui_scale = GetWindowScaleDPI();
        gui_config->ui_scale = ui_scale.x;
    }
    const int font_size_small = 16;
    const int font_size_big = 20;

    struct gui_state *gui_state = new_gui_state_with_font_sizes((float)font_size_small, (float)font_size_big, gui_config->ui_scale);

    free(window_title);

    SetTargetFPS(60);
    Image icon = LoadImage("res/icon.png");

    if(icon.data) {
        SetWindowIcon(icon);
    }

    UnloadImage(icon);

    float scale = 1.0f;

    Vector3 mesh_offset = (Vector3){0, 0, 0};

    const char *help_box_strings[] = {" - Mouse Wheel to Zoom in-out",
                                      " - Mouse Wheel Pressed to Pan",
                                      " - Alt + Mouse Wheel Pressed to Rotate",
                                      " - Alt + Ctrl + Mouse Wheel Pressed for Smooth Zoom",
                                      " - Ctrl + F to search a cell based on it's center",
                                      " - G to only draw the grid lines",
                                      " - L to enable or disable the grid lines",
                                      " - R to restart simulation (only works when paused)",
                                      " - Alt + R to restart simulation and the box positions",
                                      " - X to show/hide AP visualization",
                                      " - Q to show/hide scale",
                                      " - C to show/hide everything except grid",
                                      " - F to open a simulation file",
                                      " - O to open a simulation directory",
                                      " - F12 to take a screenshot",
                                      " - CTRL + F12 to start/stop recording the screen",
                                      " - . or , to change color scales",
                                      " - Right arrow to advance one dt when paused",
                                      " - Hold up arrow to advance time when paused",
                                      " - Double click on a volume to show the AP",
                                      " - Space to start or pause simulation"};

    int help_box_lines = SIZEOF(help_box_strings);

    float wider_text_w = 0;
    float wider_text_h = 0;

    for(int i = 0; i < help_box_lines; i++) {
        Vector2 tmp = MeasureTextEx(gui_state->font, help_box_strings[i], gui_state->font_size_small, gui_state->font_spacing_small);
        if(tmp.x > wider_text_w) {
            wider_text_w = tmp.x;
            wider_text_h = tmp.y;
        }
    }

    gui_state->help_box.lines = (char **)help_box_strings;
    gui_state->help_box.title = "Default controls";

    float text_offset = 1.5f * wider_text_h;
    float box_w = wider_text_w * 1.08f;

    gui_state->end_info_box.lines = (char **)malloc(sizeof(char *) * end_info_box_lines);
    gui_state->end_info_box.title = "Simulation time";

    gui_state->mesh_info_box.title = "Mesh information";
    gui_state->mesh_info_box.lines = (char **)malloc(sizeof(char *) * mesh_info_box_lines);
    configure_info_boxes_sizes(gui_state, help_box_lines, mesh_info_box_lines, end_info_box_lines, box_w, text_offset);

    Vector2 error_message_width;
    struct mesh_info *mesh_info = new_mesh_info();
    bool end_info_box_strings_configured = false;

    Shader shader = LoadShader("res/instanced_vertex_shader.vs", "res/fragment_shader.fs");

    shader.locs[SHADER_LOC_MATRIX_MVP] = GetShaderLocation(shader, "mvp");
    shader.locs[SHADER_LOC_MATRIX_MODEL] = GetShaderLocationAttrib(shader, "instanceTransform");
    shader.locs[SHADER_LOC_VERTEX_COLOR] = GetShaderLocationAttrib(shader, "color");
    Mesh cube = GenMeshCube(1.0f, 1.0f, 1.0f);

    int grid_mask = 0;

    while(!WindowShouldClose()) {

        if(gui_state->draw_grid_only) {
            grid_mask = 2;
        } else if(gui_state->draw_grid_lines) {
            grid_mask = 1;
        } else {
            grid_mask = 0;
        }

        UpdateCamera(&(gui_state->camera));

        BeginDrawing();

        if(IsWindowResized()) {
            gui_state->current_window_width = GetScreenWidth();
            gui_state->current_window_height = GetScreenHeight();
            reset_ui(gui_state);
        }

        gui_state->handle_keyboard_input = !gui_state->show_selection_box;

        handle_input(gui_config, mesh_info, gui_state);
        if(gui_config->grid_info.loaded) {

            omp_set_lock(&gui_config->draw_lock);

            if(draw_type == DRAW_FILE) {
                if(gui_config->grid_info.file_name) {
                    size_t max = strlen(gui_config->grid_info.file_name) + strlen("Visualizing file - ") + 2;
                    window_title = (char *)malloc(max);
                    snprintf(window_title, max + 1, "Visualizing file - %s", gui_config->grid_info.file_name);
                    SetWindowTitle(window_title);
                    free(window_title);
                }
            }

            ClearBackground(GRAY);

            BeginMode3D(gui_state->camera);

            if(!mesh_info->center_calculated) {
                if(draw_type == DRAW_SIMULATION) {
                    mesh_offset = find_mesh_center(gui_config->grid_info.alg_grid, mesh_info);
                } else {
                    mesh_offset = find_mesh_center_vtk(gui_config->grid_info.vtk_grid, mesh_info);
                }
                scale = fmaxf(mesh_offset.x, fmaxf(mesh_offset.y, mesh_offset.z)) / 1.8f;
            }

            if(draw_type == DRAW_SIMULATION) {
                draw_alg_mesh(gui_config, mesh_offset, scale, gui_state, shader, cube, grid_mask);
            } else if(draw_type == DRAW_FILE) {
                draw_vtk_unstructured_grid(gui_config, mesh_offset, scale, gui_state, shader, cube, grid_mask);
            }

            gui_state->double_clicked = false;
            if(gui_state->show_coordinates) {
                draw_coordinates(gui_state);
            }

            EndMode3D();

            // We finished drawing everything that depends on the mesh being loaded
            omp_unset_lock(&gui_config->draw_lock);

            if(gui_state->controls_window.show) {
                draw_control_window(gui_state, gui_config);
            }

            if(gui_state->show_coordinates) {
                DrawText("x", (int)gui_state->coordinates_label_x_position.x, (int)gui_state->coordinates_label_x_position.y, (int)gui_state->font_size_big,
                         RED);
                DrawText("y", (int)gui_state->coordinates_label_y_position.x, (int)gui_state->coordinates_label_y_position.y, (int)gui_state->font_size_big,
                         GREEN);
                DrawText("z", (int)gui_state->coordinates_label_z_position.x, (int)gui_state->coordinates_label_z_position.y, (int)gui_state->font_size_big,
                         DARKBLUE);
            }

            if(gui_state->mesh_info_box.window.show) {
                bool configured = configure_mesh_info_box_strings(gui_state, gui_config, &gui_state->mesh_info_box.lines, draw_type, mesh_info);

                if(configured) {
                    draw_box(&gui_state->mesh_info_box, gui_state, text_offset);

                    for(int i = 0; i < mesh_info_box_lines; i++) {
                        free(gui_state->mesh_info_box.lines[i]);
                    }
                }
            }

            if(gui_state->scale.window.show) {
                draw_scale(gui_config->min_v, gui_config->max_v, gui_state, gui_config->int_scale);
            }

            if(hmlen(gui_state->ap_graph_config->selected_aps) && gui_state->show_ap) {
                draw_ap_graph(gui_state, gui_config);
            }

            if(gui_state->help_box.window.show) {
                draw_box(&gui_state->help_box, gui_state, text_offset);
            }

            if(!gui_config->simulating) {
                if(draw_type == DRAW_SIMULATION) {
                    if(gui_state->end_info_box.window.show) {
                        if(!end_info_box_strings_configured) {
                            configure_end_info_box_strings(gui_config, &gui_state->end_info_box.lines);
                            end_info_box_strings_configured = true;
                        }
                        draw_box(&gui_state->end_info_box, gui_state, text_offset);
                    }
                }
            }

            if(!gui_config->paused) {
                omp_unset_lock(&gui_config->sleep_lock);
            }

            if(gui_state->show_selection_box) {
                gui_state->show_selection_box = !draw_selection_box(gui_state);
            }

        } else {

            ClearBackground(GRAY);
            float spacing = gui_state->font_spacing_big;
            static Color c = RED;

            if(!gui_config->error_message) {
                gui_config->error_message = strdup("Loading Mesh...");
                c = WHITE;
            }

            error_message_width = MeasureTextEx(gui_state->font, gui_config->error_message, gui_state->font_size_big, spacing);

            int posx = GetScreenWidth() / 2 - (int)error_message_width.x / 2;
            int posy = GetScreenHeight() / 2 - 50;

            int rec_width = (int)(error_message_width.x) + 50;
            int rec_height = (int)(error_message_width.y) + 2;

            DrawRectangle(posx, posy, rec_width, rec_height, c);
            DrawRectangleLines(posx, posy, rec_width, rec_height, BLACK);

            // This should not happen... but it does....
            if(gui_config->error_message) {
                DrawTextEx(gui_state->font, gui_config->error_message, (Vector2){(float)posx + ((float)rec_width - error_message_width.x) / 2, (float)posy},
                           gui_state->font_size_big, gui_state->font_spacing_big, BLACK);
            }
        }

        // Draw FPS
        int fps = GetFPS();
        const char *text = TextFormat("%2i FPS - Frame Time %lf", fps, GetFrameTime());
        Vector2 text_size = MeasureTextEx(gui_state->font, text, gui_state->font_size_big, gui_state->font_spacing_big);

        DrawTextEx(gui_state->font, text,
                   (Vector2){((float)gui_state->current_window_width - text_size.x - 10.0f), ((float)gui_state->current_window_height - text_size.y - 30)},
                   gui_state->font_size_big, gui_state->font_spacing_big, BLACK);

        text_size = MeasureTextEx(gui_state->font, "Press H to show/hide the help box", gui_state->font_size_big, gui_state->font_spacing_big);

        DrawTextEx(gui_state->font, "Press H to show/hide the help box", (Vector2){10.0f, ((float)gui_state->current_window_height - text_size.y - 30.0f)},
                   gui_state->font_size_big, gui_state->font_spacing_big, BLACK);

        float upper_y = text_size.y + 30;

        if(gui_state->current_mouse_over_volume.position_draw.x != -1) {

            Vector2 info_pos;

            if(gui_config->draw_type == DRAW_SIMULATION) {
                text = TextFormat("Mouse is on Volume: %.2lf, %.2lf, %.2lf with grid position %i", gui_state->current_mouse_over_volume.position_draw.x,
                                  gui_state->current_mouse_over_volume.position_draw.y, gui_state->current_mouse_over_volume.position_draw.z,
                                  gui_state->current_mouse_over_volume.matrix_position + 1);

            } else {

                text = TextFormat("Mouse is on Volume: %.2lf, %.2lf, %.2lf", gui_state->current_mouse_over_volume.position_draw.x,
                                  gui_state->current_mouse_over_volume.position_draw.y, gui_state->current_mouse_over_volume.position_draw.z);
            }

            text_size = MeasureTextEx(gui_state->font, text, gui_state->font_size_big, gui_state->font_spacing_big);
            info_pos =
                (Vector2){((float)gui_state->current_window_width - text_size.x - 10), ((float)gui_state->current_window_height - text_size.y - upper_y)};
            DrawTextEx(gui_state->font, text, info_pos, gui_state->font_size_big, gui_state->font_spacing_big, BLACK);
        }

        EndDrawing();
    }

    gui_config->exit = true;
    free(mesh_info);

    if(end_info_box_strings_configured) {
        for(int i = 0; i < end_info_box_lines; i++) {
            free(gui_state->end_info_box.lines[i]);
        }
    }

    free(gui_state->end_info_box.lines);

    hmfree(gui_state->ap_graph_config->selected_aps);

    free(gui_state->ap_graph_config);
    free(gui_state);

    omp_unset_lock(&gui_config->draw_lock);
    omp_unset_lock(&gui_config->sleep_lock);

    CloseWindow();
}
