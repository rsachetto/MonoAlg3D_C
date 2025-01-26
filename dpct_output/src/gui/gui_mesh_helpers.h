#ifndef __GUI_MESH_HELPERS_H
#define __GUI_MESH_HELPERS_H 

#include "gui.h"
#include "raylib_ext.h"

Vector3 find_mesh_center(struct grid *grid_to_draw, struct mesh_info *mesh_info);
Vector3 find_mesh_center_vtk(struct vtk_unstructured_grid *grid_to_draw, struct mesh_info *mesh_info);
struct mesh_info *new_mesh_info();
void draw_alg_mesh(struct gui_shared_info *gui_config, struct gui_state *gui_state, int grid_mask, struct draw_context *draw_context);
void draw_vtk_unstructured_grid(struct gui_shared_info *gui_config, struct gui_state *gui_state, int grid_mask, struct draw_context *draw_context);
void draw_alg_purkinje_network(struct gui_shared_info *gui_config, struct gui_state *gui_state);
void draw_vtk_purkinje_network(struct gui_shared_info *gui_config, struct gui_state *gui_state);
void set_visibility_after_split(struct gui_shared_info *gui_config, struct gui_state *gui_state);
void reset_grid_visibility(struct gui_shared_info *gui_config, struct gui_state *gui_state);

#endif /* __GUI_MESH_HELPERS_H */
