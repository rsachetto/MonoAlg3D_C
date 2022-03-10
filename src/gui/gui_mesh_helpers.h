#ifndef __GUI_MESH_HELPERS_H
#define __GUI_MESH_HELPERS_H 

#include "gui.h"

Vector3 find_mesh_center(struct grid *grid_to_draw, struct mesh_info *mesh_info);
Vector3 find_mesh_center_vtk(struct vtk_unstructured_grid *grid_to_draw, struct mesh_info *mesh_info);
void draw_alg_mesh(struct gui_shared_info *gui_config, struct gui_state *gui_state, int grid_mask, Shader shader, Mesh cube);
void draw_alg_purkinje_network(struct gui_shared_info *gui_config, struct gui_state *gui_state, int grid_mask);
void draw_vtk_unstructured_grid(struct gui_shared_info *gui_config, struct gui_state *gui_state, int grid_mask, Shader shader, Mesh cube);
void draw_vtk_purkinje_network(struct gui_shared_info *gui_config, struct gui_state *gui_state, int grid_mask);
void set_visibility_after_split(struct gui_shared_info *gui_config, struct gui_state *gui_state);
void reset_grid_visibility(struct gui_shared_info *gui_config, struct gui_state *gui_state);

#endif /* __GUI_MESH_HELPERS_H */
