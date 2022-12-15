#ifndef __MESH_INFO_DATA_H
#define __MESH_INFO_DATA_H

#include <stdbool.h>
#include "../common_types/common_types.h"

struct default_fibrotic_mesh_info {
    bool fibrotic;
    bool border_zone;
    int tissue_type;
};

struct dti_mesh_info {

	enum dti_transmurality_labels {
		DTI_FAST_ENDO,
		DTI_ENDO,
		DTI_MID,
		DTI_EPI
	} dti_transmurality_labels;
	float transmurality;
	float base_apex_heterogeneity;
};

#define FIBROTIC_INFO(grid_cell) (struct default_fibrotic_mesh_info *)(grid_cell)->mesh_extra_info
#define FIBROTIC(grid_cell) (FIBROTIC_INFO(grid_cell))->fibrotic
#define BORDER_ZONE(grid_cell) (FIBROTIC_INFO(grid_cell))->border_zone
#define TISSUE_TYPE(grid_cell) (FIBROTIC_INFO(grid_cell))->tissue_type

#define INITIALIZE_FIBROTIC_INFO(grid_cell)                                                                            \
    do {                                                                                                               \
        size_t __size__ = sizeof (struct default_fibrotic_mesh_info);                                                  \
        (grid_cell)->mesh_extra_info = malloc (__size__);                                                              \
        (grid_cell)->mesh_extra_info_size = __size__;                                                                  \
        FIBROTIC ((grid_cell)) = false;                                                                                \
        BORDER_ZONE (grid_cell) = false;                                                                               \
        TISSUE_TYPE ((grid_cell)) = 0;                                                                                 \
} while (0)

#define INITIALIZE_DTI_MESH_INFO(grid_cell) ALLOCATE_MESH_INFO(grid_cell, dti_mesh_info);   
#define DTI_MESH_INFO(grid_cell) (struct dti_mesh_info *)grid_cell->mesh_extra_info
#define DTI_MESH_TRANSMURALITY_LABELS(grid_cell) (DTI_MESH_INFO(grid_cell))->dti_transmurality_labels
#define DTI_MESH_TRANSMURALITY(grid_cell) (DTI_MESH_INFO(grid_cell))->transmurality
#define DTI_MESH_BASE_APEX_HETEROGENEITY(grid_cell) (DTI_MESH_INFO(grid_cell))->base_apex_heterogeneity

#endif /* __MESH_INFO_DATA_H */
