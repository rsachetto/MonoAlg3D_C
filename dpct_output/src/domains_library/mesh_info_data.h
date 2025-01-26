#ifndef __MESH_INFO_DATA_H
#define __MESH_INFO_DATA_H

#include <stdbool.h>
#include "../common_types/common_types.h"

struct default_fibrotic_mesh_info {
    bool fibrotic;
    bool border_zone;
    int tissue_type;
};

// TODO: Move this struct and its macros to "custom_mesh_info_data.h"!
struct dti_mesh_info {
		float dti_transmurality_labels;
		float transmurality;
    	float apicobasal;
		float base_apex_heterogeneity;
    	float sf_IKs;
    	int fast_endo;
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
#define DTI_MESH_APICOBASAL(grid_cell) (DTI_MESH_INFO(grid_cell))->apicobasal
#define DTI_MESH_FAST_ENDO(grid_cell) (DTI_MESH_INFO(grid_cell))->fast_endo
#define DTI_MESH_SCALE_FACTOR_IKS(grid_cell) (DTI_MESH_INFO(grid_cell))->sf_IKs

#endif /* __MESH_INFO_DATA_H */
