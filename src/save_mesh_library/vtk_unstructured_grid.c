//
// Created by sachetto on 30/10/18.
//

#include "vtk_unstructured_grid.h"
#include "../alg/cell/cell.h"
#include "../string/sds.h"
#include "data_utils.h"
#include <inttypes.h>
#include <math.h>
#include <stdint.h>

struct vtk_unstructured_grid *new_vtk_unstructured_grid() {
    struct vtk_unstructured_grid *grid = (struct vtk_unstructured_grid *)malloc(sizeof(struct vtk_unstructured_grid));

    grid->cells = NULL;
    grid->values = NULL;
    grid->points = NULL;

    grid->num_cells = 0;
    grid->num_points = 0;
    grid->points_per_cell = 8;
    grid->cell_type = 12;

    return grid;
}

void free_vtk_unstructured_grid(struct vtk_unstructured_grid *vtk_grid) {
    if(vtk_grid) {
        sb_free(vtk_grid->cells);
        sb_free(vtk_grid->values);
        sb_free(vtk_grid->points);
    }
}

void new_vtk_unstructured_grid_from_alg_grid(struct vtk_unstructured_grid **vtk_grid, struct grid *grid, bool clip_with_plain,
                                                                      float *plain_coordinates, bool clip_with_bounds,
                                                                      float *bounds, bool read_only_values) {

    static bool mesh_already_loaded =  false;

    if(grid == NULL) {
        return;
    }

    if(!read_only_values) {
        *vtk_grid = new_vtk_unstructured_grid();
    }
    else {
        if(!(*vtk_grid) && mesh_already_loaded) {
            fprintf(stderr,
                    "Function new_vtk_unstructured_grid_from_alg_grid can only be called with read_only_values if the grid is already loaded");
            exit(EXIT_FAILURE);
        }

        if(mesh_already_loaded) {
            assert(*vtk_grid);
            sb_free((*vtk_grid)->values);
            (*vtk_grid)->values = NULL;
        }
        else {
            *vtk_grid = new_vtk_unstructured_grid();
        }
    }

    float min_x = 0.0;
    float min_y = 0.0;
    float min_z = 0.0;
    float max_x = 0.0;
    float max_y = 0.0;
    float max_z = 0.0;

    float p0[3] = {0, 0, 0};
    float n[3] = {0, 0, 0};

    if(!plain_coordinates) {
        clip_with_plain = false;
    } else {
        p0[0] = plain_coordinates[0];
        p0[1] = plain_coordinates[1];
        p0[2] = plain_coordinates[2];

        n[0] = plain_coordinates[3];
        n[0] = plain_coordinates[4];
        n[0] = plain_coordinates[5];
    }

    if(!bounds) {
        clip_with_bounds = false;
    } else {
        min_x = bounds[0];
        min_y = bounds[1];
        min_z = bounds[2];
        max_x = bounds[3];
        max_y = bounds[4];
        max_z = bounds[5];
    }

    struct cell_node *grid_cell = grid->first_cell;

    float center_x, center_y, center_z;
    //double v;

    struct point_3d aux1;
    struct point_3d aux2;
    struct point_3d aux3;
    struct point_3d aux4;
    struct point_3d aux5;
    struct point_3d aux6;
    struct point_3d aux7;
    struct point_3d aux8;

    uint32_t id = 0;
    uint32_t num_cells = 0;

    float l = sqrtf(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    float A = n[0] / l;
    float B = n[1] / l;
    float C = n[2] / l;
    float D = -(n[0] * p0[0] + n[1] * p0[1] + n[2] * p0[2]);

    double side;
    struct point_hash *hash = point_hash_create();

    while(grid_cell != 0) {

        if(grid_cell->active) {

            center_x = grid_cell->center_x;
            center_y = grid_cell->center_y;
            center_z = grid_cell->center_z;

            if(clip_with_plain) {
                side = A * center_x + B * center_y + C * center_z + D;
                if(side < 0) {
                    grid_cell = grid_cell->next;
                    continue;
                }
            }

            if(clip_with_bounds) {
                bool ignore_cell = center_x < min_x || center_x > max_x || center_y < min_y || center_y > max_y ||
                                   center_z < min_z || center_z > max_z;

                if(ignore_cell) {
                    grid_cell = grid_cell->next;
                    continue;
                }
            }

            sb_push((*vtk_grid)->values, grid_cell->v);

            if(mesh_already_loaded && read_only_values) {
                grid_cell = grid_cell->next;
                continue;
            }

            float half_face_x = grid_cell->dx / 2.0f;
            float half_face_y = grid_cell->dy / 2.0f;
            float half_face_z = grid_cell->dz / 2.0f;

            aux1.x = center_x - half_face_x;
            aux1.y = center_y - half_face_y;
            aux1.z = center_z - half_face_z;

            aux2.x = center_x + half_face_x;
            aux2.y = center_y - half_face_y;
            aux2.z = center_z - half_face_z;

            aux3.x = center_x + half_face_x;
            aux3.y = center_y + half_face_y;
            aux3.z = center_z - half_face_z;

            aux4.x = center_x - half_face_x;
            aux4.y = center_y + half_face_y;
            aux4.z = center_z - half_face_z;

            aux5.x = center_x - half_face_x;
            aux5.y = center_y - half_face_y;
            aux5.z = center_z + half_face_z;

            aux6.x = center_x + half_face_x;
            aux6.y = center_y - half_face_y;
            aux6.z = center_z + half_face_z;

            aux7.x = center_x + half_face_x;
            aux7.y = center_y + half_face_y;
            aux7.z = center_z + half_face_z;

            aux8.x = center_x - half_face_x;
            aux8.y = center_y + half_face_y;
            aux8.z = center_z + half_face_z;

            if(point_hash_search(hash, aux1) == -1) {
                sb_push((*vtk_grid)->points, aux1);
                point_hash_insert(hash, aux1, id);
                id++;
            }

            if(point_hash_search(hash, aux2) == -1) {
                sb_push((*vtk_grid)->points, aux2);
                point_hash_insert(hash, aux2, id);
                id++;
            }

            if(point_hash_search(hash, aux3) == -1) {
                point_hash_insert(hash, aux3, id);
                sb_push((*vtk_grid)->points, aux3);
                id++;
            }

            if(point_hash_search(hash, aux4) == -1) {
                point_hash_insert(hash, aux4, id);
                sb_push((*vtk_grid)->points, aux4);
                id++;
            }

            if(point_hash_search(hash, aux5) == -1) {
                sb_push((*vtk_grid)->points, aux5);
                point_hash_insert(hash, aux5, id);
                id++;
            }

            if(point_hash_search(hash, aux6) == -1) {
                sb_push((*vtk_grid)->points, aux6);
                point_hash_insert(hash, aux6, id);
                id++;
            }

            if(point_hash_search(hash, aux7) == -1) {
                sb_push((*vtk_grid)->points, aux7);
                point_hash_insert(hash, aux7, id);
                id++;
            }

            if(point_hash_search(hash, aux8) == -1) {
                sb_push((*vtk_grid)->points, aux8);
                point_hash_insert(hash, aux8, id);
                id++;
            }

            sb_push((*vtk_grid)->cells, point_hash_search(hash, aux1));
            sb_push((*vtk_grid)->cells, point_hash_search(hash, aux2));
            sb_push((*vtk_grid)->cells, point_hash_search(hash, aux3));
            sb_push((*vtk_grid)->cells, point_hash_search(hash, aux4));
            sb_push((*vtk_grid)->cells, point_hash_search(hash, aux5));
            sb_push((*vtk_grid)->cells, point_hash_search(hash, aux6));
            sb_push((*vtk_grid)->cells, point_hash_search(hash, aux7));
            sb_push((*vtk_grid)->cells, point_hash_search(hash, aux8));
            num_cells++;
        }

        grid_cell = grid_cell->next;
    }

    if(!mesh_already_loaded) {
        (*vtk_grid)->num_cells = num_cells;
        (*vtk_grid)->num_points = id;

        if(read_only_values)
            mesh_already_loaded = true;
    }
}

sds create_common_vtu_header(bool compressed, int num_points, int num_cells) {

    sds header = sdsempty();

    if(compressed) {
        header = sdscat(header, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" "
                                "header_type=\"UInt64\" compressor=\"vtkZLibDataCompressor\">\n");
    } else {
        header = sdscat(
            header,
            "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
    }

    header = sdscat(header, "  <UnstructuredGrid>\n");

    header = sdscatprintf(header, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", num_points, num_cells);

    header = sdscat(header, "      <CellData Scalars=\"Scalars_\">\n");

    return header;
}

void save_vtk_unstructured_grid_as_vtu(struct vtk_unstructured_grid *vtk_grid, char *filename, bool binary) {

    size_t offset = 0;

    sds file_content = create_common_vtu_header(false, vtk_grid->num_points, vtk_grid->num_cells);

    if(binary) {
        file_content = sdscat(
            file_content,
            "        <DataArray type=\"Float32\" Name=\"Scalars_\" format=\"appended\" offset=\"0\">\n"); // First
                                                                                                          // offset is
                                                                                                          // always 0

    } else {
        file_content =
            sdscat(file_content, "        <DataArray type=\"Float32\" Name=\"Scalars_\" format=\"ascii\">\n");
    }

    if(!binary) {
        size_t num_values = sb_count(vtk_grid->values);

        for(int i = 0; i < num_values; i++) {
            file_content = sdscatprintf(file_content, "     %lf ", vtk_grid->values[i]);
        }
    }

    file_content = sdscat(file_content, "        </DataArray>\n");
    file_content = sdscat(file_content, "      </CellData>\n");

    offset = (vtk_grid->num_cells * 4) + 8;

    file_content = sdscat(file_content, "      <Points>\n");

    if(binary) {
        file_content = sdscatprintf(file_content,
                                    "        <DataArray type=\"Float32\" Name=\"Points\" "
                                    "NumberOfComponents=\"3\" format=\"appended\" offset=\"%zu\">\n",
                                    offset);

    } else {
        file_content =
            sdscat(file_content,
                   "        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    }

    if(!binary) {
        int num_points = sb_count(vtk_grid->points);
        for(int i = 0; i < num_points; i++) {
            struct point_3d p = vtk_grid->points[i];
            file_content = sdscatprintf(file_content, "%lf %lf %lf\n", p.x, p.y, p.z);
        }
    }

    file_content = sdscat(file_content, "        </DataArray>\n");
    file_content = sdscat(file_content, "      </Points>\n");

    file_content = sdscat(file_content, "      <Cells>\n");

    offset += (vtk_grid->num_points * 4 * 3) + 8; // 3*32 bits float for each point

    if(binary) {
        file_content = sdscatprintf(
            file_content,
            "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"appended\" offset=\"%zu\">\n", offset);

    } else {
        file_content =
            sdscat(file_content, "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    }

    int points_per_cell = vtk_grid->points_per_cell;
    int cell_type = vtk_grid->cell_type;

    int num_cells = vtk_grid->num_cells;

    if(!binary) {
        for(int i = 0; i < num_cells; i++) {
            file_content = sdscat(file_content, "     ");
            for(int j = 0; j < points_per_cell; j++) {
                file_content = sdscatprintf(file_content, "%d ", vtk_grid->cells[points_per_cell * i + j]);
            }

            file_content = sdscat(file_content, "\n");
        }
    }

    file_content = sdscat(file_content, "        </DataArray>\n");

    offset += (vtk_grid->num_cells * 8 * points_per_cell) + 8; // 64 bits

    if(binary) {
        file_content = sdscatprintf(
            file_content, "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"appended\" offset=\"%zu\">\n",
            offset);
    } else {
        file_content = sdscat(file_content, "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    }

    {
        int offset_local = 8;

        if(!binary) {
            for(int i = 0; i < num_cells; i++) {
                file_content = sdscat(file_content, "     ");
                file_content = sdscatprintf(file_content, "%d ", offset_local);
                offset_local += 8;
                file_content = sdscat(file_content, "\n");
            }
        }
    }

    file_content = sdscat(file_content, "        </DataArray>\n");

    offset += (vtk_grid->num_cells * 8) + 8; // 64 bits

    if(binary) {
        file_content = sdscatprintf(
            file_content, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\"%zu\">\n",
            offset);
    } else {
        file_content = sdscat(file_content, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    }

    if(!binary) {

        for(int i = 0; i < num_cells; i++) {
            file_content = sdscatprintf(file_content, "    %d\n", cell_type);
        }
    }

    file_content = sdscat(file_content, "        </DataArray>\n");

    file_content = sdscat(file_content, "      </Cells>\n");

    file_content = sdscat(file_content, "    </Piece>\n");
    file_content = sdscat(file_content, "  </UnstructuredGrid>\n");

    size_t size_until_now = 0;

    if(binary) {
        file_content = sdscat(file_content, "  <AppendedData encoding=\"raw\">\n   _");

        size_until_now = sdslen(file_content);

        // scalars
        uint64_t block_size = sizeof(float) * vtk_grid->num_cells;
        file_content = sdscatlen(file_content, &block_size, sizeof(uint64_t));
        file_content = sdscatlen(file_content, vtk_grid->values, (size_t)block_size);
        size_until_now += (sizeof(float) * vtk_grid->num_cells + sizeof(uint64_t));

        // Points
        block_size = sizeof(struct point_3d) * vtk_grid->num_points;
        file_content = sdscatlen(file_content, &block_size, sizeof(uint64_t));
        file_content = sdscatlen(file_content, vtk_grid->points, (size_t)block_size);
        size_until_now += (sizeof(struct point_3d) * vtk_grid->num_points + sizeof(uint64_t));

        // connectivity
        block_size = num_cells * points_per_cell * sizeof(int64_t);
        file_content = sdscatlen(file_content, &block_size, sizeof(uint64_t));
        size_until_now += sizeof(uint64_t);

        for(int i = 0; i < num_cells; i++) {
            for(int j = 0; j < points_per_cell; j++) {
                int64_t aux = (int64_t)vtk_grid->cells[points_per_cell * i + j];
                file_content = sdscatlen(file_content, &aux, sizeof(int64_t));
                size_until_now += sizeof(int64_t);
            }
        }

        // offsets
        block_size = num_cells * sizeof(int64_t);
        file_content = sdscatlen(file_content, &block_size, sizeof(uint64_t));
        size_until_now += sizeof(uint64_t);

        int64_t offset_local = 8;
        for(int i = 0; i < num_cells; i++) {
            file_content = sdscatlen(file_content, &offset_local, sizeof(int64_t));
            offset_local += 8;
            size_until_now += sizeof(int64_t);
        }

        // types
        block_size = num_cells * sizeof(uint8_t);
        file_content = sdscatlen(file_content, &block_size, sizeof(uint64_t));
        size_until_now += sizeof(uint64_t);

        for(int i = 0; i < num_cells; i++) {
            uint8_t aux = (uint8_t)cell_type;
            file_content = sdscatlen(file_content, &aux, sizeof(uint8_t));
            size_until_now += sizeof(uint8_t);
        }

        file_content = sdscat(file_content, "\n  </AppendedData>\n");
    }

    file_content = sdscat(file_content, "</VTKFile>\n");

    size_until_now += 29;

    FILE *output_file = NULL;

    if(binary) {
        output_file = fopen(filename, "wb");
        fwrite(file_content, size_until_now, 1, output_file);
    } else {
        output_file = fopen(filename, "w");
        fprintf(output_file, "%s", file_content);
    }

    sdsfree(file_content);
    fclose(output_file);
}

//TODO: for non adaptive meshes we need to compress only the values data array. The other arrays remain the same. So we need to save them somewhere
void save_vtk_unstructured_grid_as_vtu_compressed(struct vtk_unstructured_grid *vtk_grid, char *filename, int compression_level) {

    sds first_file_part = create_common_vtu_header(true, vtk_grid->num_points, vtk_grid->num_cells);

    size_t offset = 0;

    first_file_part = sdscat(
        first_file_part, "        <DataArray type=\"Float32\" Name=\"Scalars_\" format=\"appended\" offset=\"0\" />\n");

    first_file_part = sdscat(first_file_part, "      </CellData>\n");
    first_file_part = sdscat(first_file_part, "      <Points>\n");

    sds points_array_header = sdsnew("        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" "
                                     "format=\"appended\" offset=\"%zu\" />\n");

    sds points_array_header_end = sdsnew("      </Points>\n");
    points_array_header_end = sdscat(points_array_header_end, "      <Cells>\n");

    sds connectivity_array_header =
        sdsnew("        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"appended\" offset=\"%zu\"/>\n");

    sds offsets_array_header =
        sdsnew("        <DataArray type=\"Int64\" Name=\"offsets\" format=\"appended\" offset=\"%zu\"/>\n");

    sds types_array_header =
        sdsnew("        <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\"%zu\"/>\n");


    sds data_end = sdsnew("      </Cells>\n");

    data_end = sdscat(data_end, "    </Piece>\n");
    data_end = sdscat(data_end, "  </UnstructuredGrid>\n");

    sds appended_begin = sdsnew("  <AppendedData encoding=\"raw\">\n   _");

    FILE *output_file = NULL;

    output_file = fopen(filename, "wb");

    fwrite(first_file_part, sdslen(first_file_part), 1, output_file);
    sdsfree(first_file_part);
    {

        size_t data_size;

        size_t data_size_after_compression_for_values;
        size_t num_block_for_values;
        size_t block_size_uncompressed_for_values;
        size_t *block_sizes_compressed_for_values;
        size_t last_block_size_for_values;

        unsigned char *data_to_compress;

        data_size = vtk_grid->num_cells * 4; // 32 bit float for each cell

        data_to_compress = (unsigned char *)vtk_grid->values;
        unsigned char *compressed_data_for_values;

        calculate_blocks_and_compress_data(data_size, &data_size_after_compression_for_values, data_to_compress,
                                           &(compressed_data_for_values), &num_block_for_values,
                                           &block_size_uncompressed_for_values, &block_sizes_compressed_for_values,
                                           &last_block_size_for_values, compression_level);

        offset = 8 + 8 + 8 + (num_block_for_values * 8) + data_size_after_compression_for_values;

        sds tmp = sdscatprintf(sdsempty(), points_array_header, offset);
        fwrite(tmp, sdslen(tmp), 1, output_file);
        sdsfree(tmp);
        sdsfree(points_array_header);

        fwrite(points_array_header_end, sdslen(points_array_header_end), 1, output_file);
        sdsfree(points_array_header_end);

        size_t data_size_after_compression_for_points;
        size_t num_block_for_points;
        size_t block_size_uncompressed_for_points;
        size_t *block_sizes_compressed_for_points;
        size_t last_block_size_for_points;

        data_size = vtk_grid->num_points * 4 * 3; // 3 points, with 32 bit float each point

        data_to_compress = (unsigned char *)vtk_grid->points;

        unsigned char *compressed_data_for_points;

        calculate_blocks_and_compress_data(data_size, &data_size_after_compression_for_points, data_to_compress,
                                           &(compressed_data_for_points), &num_block_for_points,
                                           &block_size_uncompressed_for_points, &block_sizes_compressed_for_points,
                                           &last_block_size_for_points, compression_level);

        offset += 8 + 8 + 8 + (num_block_for_points * 8) + data_size_after_compression_for_points;

        tmp = sdscatprintf(sdsempty(), connectivity_array_header, offset);
        fwrite(tmp, sdslen(tmp), 1, output_file);
        sdsfree(tmp);
        sdsfree(connectivity_array_header);

        size_t data_size_after_compression_for_connectivity;
        size_t num_block_for_connectivity;
        size_t block_size_uncompressed_for_connectivity;
        size_t *block_sizes_compressed_for_connectivity;
        size_t last_block_size_for_connectivity;

        sds aux_data = sdsempty();
        for(int i = 0; i < vtk_grid->num_cells; i++) {
            for(int j = 0; j < vtk_grid->points_per_cell; j++) {
                int64_t aux = (int64_t)vtk_grid->cells[vtk_grid->points_per_cell * i + j];
                aux_data = sdscatlen(aux_data, &aux, sizeof(int64_t));
            }
        }

        data_size = vtk_grid->num_cells * vtk_grid->points_per_cell * sizeof(int64_t);

        data_to_compress = (unsigned char *)aux_data;
        unsigned char *compressed_data_for_connectivity;

        calculate_blocks_and_compress_data(data_size, &data_size_after_compression_for_connectivity, data_to_compress,
                                           &compressed_data_for_connectivity, &num_block_for_connectivity,
                                           &block_size_uncompressed_for_connectivity,
                                           &block_sizes_compressed_for_connectivity, &last_block_size_for_connectivity, compression_level);

        offset += 8 + 8 + 8 + (num_block_for_connectivity * 8) + data_size_after_compression_for_connectivity;

        tmp = sdscatprintf(sdsempty(), offsets_array_header, offset);
        fwrite(tmp, sdslen(tmp), 1, output_file);
        sdsfree(tmp);
        sdsfree(offsets_array_header);

        sdsfree(aux_data);

        // offsets
        size_t data_size_after_compression_for_offsets;
        size_t num_block_for_offsets;
        size_t block_size_uncompressed_for_offsets;
        size_t *block_sizes_compressed_for_offsets;
        size_t last_block_size_for_offsets;

        aux_data = sdsempty();
        int64_t offset_local = vtk_grid->points_per_cell;
        for(int i = 0; i < vtk_grid->num_cells; i++) {
            aux_data = sdscatlen(aux_data, &offset_local, sizeof(int64_t));
            offset_local += vtk_grid->points_per_cell;
        }

        data_size = vtk_grid->num_cells * sizeof(int64_t);

        data_to_compress = (unsigned char *)aux_data;
        unsigned char *compressed_data_for_offsets;

        calculate_blocks_and_compress_data(data_size, &data_size_after_compression_for_offsets, data_to_compress,
                                           &compressed_data_for_offsets, &num_block_for_offsets,
                                           &block_size_uncompressed_for_offsets, &block_sizes_compressed_for_offsets,
                                           &last_block_size_for_offsets, compression_level);

        offset += 8 + 8 + 8 + (num_block_for_offsets * 8) + data_size_after_compression_for_offsets;

        tmp = sdscatprintf(sdsempty(), types_array_header, offset);
        fwrite(tmp, sdslen(tmp), 1, output_file);
        sdsfree(tmp);
        sdsfree(types_array_header);

        sdsfree(aux_data);

        // types
        size_t data_size_after_compression_for_types;
        size_t num_block_for_types;
        size_t block_size_uncompressed_for_types;
        size_t *block_sizes_compressed_for_types;
        size_t last_block_size_for_types;

        aux_data = sdsempty();
        for(int i = 0; i < vtk_grid->num_cells; i++) {
            uint8_t aux = (uint8_t)vtk_grid->cell_type;
            aux_data = sdscatlen(aux_data, &aux, sizeof(uint8_t));
        }

        data_size = vtk_grid->num_cells * sizeof(uint8_t);

        data_to_compress = (unsigned char *)aux_data;
        unsigned char *compressed_data_for_types;

        calculate_blocks_and_compress_data(data_size, &data_size_after_compression_for_types, data_to_compress,
                                           &compressed_data_for_types, &num_block_for_types,
                                           &block_size_uncompressed_for_types, &block_sizes_compressed_for_types,
                                           &last_block_size_for_types, compression_level);

        sdsfree(aux_data);

        fwrite(data_end, sdslen(data_end), 1, output_file);
        sdsfree(data_end);

        fwrite(appended_begin, sdslen(appended_begin), 1, output_file);
        sdsfree(appended_begin);

        // Now we can save the compressed data
        //Values
        fwrite(&num_block_for_values, sizeof(uint64_t), 1, output_file);
        fwrite(&block_size_uncompressed_for_values, sizeof(uint64_t), 1, output_file);
        fwrite(&last_block_size_for_values, sizeof(uint64_t), 1, output_file);
        fwrite(block_sizes_compressed_for_values, sizeof(uint64_t), num_block_for_values, output_file);
        fwrite(compressed_data_for_values, data_size_after_compression_for_values, 1, output_file);
        free(compressed_data_for_values);
        free(block_sizes_compressed_for_values);
        
        //Points
        fwrite(&num_block_for_points, sizeof(uint64_t), 1, output_file);
        fwrite(&block_size_uncompressed_for_points, sizeof(uint64_t), 1, output_file);
        fwrite(&last_block_size_for_points, sizeof(uint64_t), 1, output_file);
        fwrite(block_sizes_compressed_for_points, sizeof(uint64_t), num_block_for_points, output_file);
        fwrite(compressed_data_for_points, data_size_after_compression_for_points, 1, output_file);

        free(compressed_data_for_points);
        free(block_sizes_compressed_for_points);

        //connectivity
        fwrite(&num_block_for_connectivity, sizeof(uint64_t), 1, output_file);
        fwrite(&block_size_uncompressed_for_connectivity, sizeof(uint64_t), 1, output_file);
        fwrite(&last_block_size_for_connectivity, sizeof(uint64_t), 1, output_file);
        fwrite(block_sizes_compressed_for_connectivity, sizeof(uint64_t), num_block_for_connectivity, output_file);
        fwrite(compressed_data_for_connectivity, data_size_after_compression_for_connectivity, 1, output_file);

        free(compressed_data_for_connectivity);
        free(block_sizes_compressed_for_connectivity);

        //offsets
        fwrite(&num_block_for_offsets, sizeof(uint64_t), 1, output_file);
        fwrite(&block_size_uncompressed_for_offsets, sizeof(uint64_t), 1, output_file);
        fwrite(&last_block_size_for_offsets, sizeof(uint64_t), 1, output_file);
        fwrite(block_sizes_compressed_for_offsets, sizeof(uint64_t), num_block_for_offsets, output_file);
        fwrite(compressed_data_for_offsets, data_size_after_compression_for_offsets, 1, output_file);

        free(compressed_data_for_offsets);
        free(block_sizes_compressed_for_offsets);

        //types
        fwrite(&num_block_for_types, sizeof(uint64_t), 1, output_file);
        fwrite(&block_size_uncompressed_for_types, sizeof(uint64_t), 1, output_file);
        fwrite(&last_block_size_for_types, sizeof(uint64_t), 1, output_file);
        fwrite(block_sizes_compressed_for_types, sizeof(uint64_t), num_block_for_types, output_file);
        fwrite(compressed_data_for_types, data_size_after_compression_for_types, 1, output_file);

        free(compressed_data_for_types);
        free(block_sizes_compressed_for_types);
    }

    sds appended_end = sdsnew("\n  </AppendedData>\n");
    appended_end = sdscat(appended_end, "</VTKFile>\n");

    fwrite(appended_end, sdslen(appended_end), 1, output_file);
    sdsfree(appended_end);

    fclose(output_file);
}

void save_vtk_unstructured_grid_as_legacy_vtk(struct vtk_unstructured_grid *vtk_grid, char *filename, bool binary) {

    sds file_content = sdsempty();

    file_content = sdscat(file_content, "# vtk DataFile Version 4.2\n");
    file_content = sdscat(file_content, "vtk output\n");

    if(binary) {
        file_content = sdscat(file_content, "BINARY\n");
    } else {
        file_content = sdscat(file_content, "ASCII\n");
    }

    file_content = sdscat(file_content, "DATASET UNSTRUCTURED_GRID\n");
    file_content = sdscatprintf(file_content, "POINTS %d float\n", vtk_grid->num_points);

    size_t size_until_now = sdslen(file_content);

    int num_points = sb_count(vtk_grid->points);
    for(int i = 0; i < num_points; i++) {
        struct point_3d p = vtk_grid->points[i];
        if(binary) {
            file_content = write_binary_point(file_content, &p);
            size_until_now += 3 * sizeof(int);
        } else {
            file_content = sdscatprintf(file_content, "%lf %lf %lf\n", p.x, p.y, p.z);
        }
    }

    int num_cells = vtk_grid->num_cells;

    {
        sds tmp = sdscatprintf(sdsempty(), "\nCELLS %d %d\n", num_cells, 9 * num_cells);

        size_until_now += sdslen(tmp);

        file_content = sdscatsds(file_content, tmp);

        sdsfree(tmp);
    }

    int points_per_cell = 8;
    int cell_type = 12;

    int points_per_cell_swapped = invert_bytes(points_per_cell);
    int cell_type_swapped = invert_bytes(cell_type);

    for(int i = 0; i < num_cells; i++) {
        if(binary) {
            file_content = sdscatlen(file_content, &points_per_cell_swapped, sizeof(int));
            size_until_now += sizeof(int);

        } else {
            file_content = sdscatprintf(file_content, "%d ", points_per_cell);
        }

        for(int j = 0; j < points_per_cell; j++) {
            if(binary) {
                int aux = invert_bytes(vtk_grid->cells[points_per_cell * i + j]);
                file_content = sdscatlen(file_content, &aux, sizeof(int));
                size_until_now += sizeof(int);
            } else {
                file_content = sdscatprintf(file_content, "%d ", vtk_grid->cells[points_per_cell * i + j]);
            }
        }

        if(!binary)
            file_content = sdscat(file_content, "\n");
    }

    {
        sds tmp = sdscatprintf(sdsempty(), "\nCELL_TYPES %d\n", num_cells);

        size_until_now += sdslen(tmp);

        file_content = sdscatsds(file_content, tmp);

        sdsfree(tmp);
    }

    for(int i = 0; i < num_cells; i++) {
        if(binary) {
            file_content = sdscatlen(file_content, &cell_type_swapped, sizeof(int));
            size_until_now += sizeof(int);
        } else {
            file_content = sdscatprintf(file_content, "%d\n", cell_type);
        }
    }

    {
        sds tmp = sdscatprintf(sdsempty(), "\nCELL_DATA %d\n", num_cells);
        tmp = sdscat(tmp, "SCALARS Scalars_ float\n");
        tmp = sdscat(tmp, "LOOKUP_TABLE default\n");

        size_until_now += sdslen(tmp);

        file_content = sdscatsds(file_content, tmp);
        sdsfree(tmp);
    }

    size_t num_values = sb_count(vtk_grid->values);

    for(int i = 0; i < num_values; i++) {
        if(binary) {
            int aux = invert_bytes(*((int *)&(vtk_grid->values[i])));
            file_content = sdscatlen(file_content, &aux, sizeof(int));
            size_until_now += sizeof(int);
        } else {
            file_content = sdscatprintf(file_content, "%lf ", vtk_grid->values[i]);
        }
    }

    {
        sds tmp = sdscat(sdsempty(), "\nMETADATA\n");
        tmp = sdscat(tmp, "INFORMATION 0\n\n");
        size_until_now += sdslen(tmp);
        file_content = sdscatsds(file_content, tmp);
        sdsfree(tmp);
    }

    FILE *output_file = fopen(filename, "w");
    if(binary) {
        fwrite(file_content, size_until_now, 1, output_file);
    } else {
        fprintf(output_file, "%s", file_content);
    }

    sdsfree(file_content);
    fclose(output_file);
}