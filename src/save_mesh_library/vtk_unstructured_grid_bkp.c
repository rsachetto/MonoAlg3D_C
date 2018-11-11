//
// Created by sachetto on 30/10/18.
//

#include "vtk_unstructured_grid.h"
#include "../alg/cell/cell.h"
#include "../string/sds.h"
#include <math.h>
#include <inttypes.h>
#include "data_utils.h"


struct vtk_unstructured_grid *new_vtk_unstructured_grid() {
    struct vtk_unstructured_grid *grid = (struct vtk_unstructured_grid *)malloc(sizeof(struct vtk_unstructured_grid));

    grid->num_cells = 0;
    grid->num_points = 0;
    grid->cells = NULL;
    grid->values = NULL;
    grid->points = NULL;

    return grid;
}

void free_vtk_unstructured_grid(struct vtk_unstructured_grid *vtk_grid) {
    if(vtk_grid) {
        sb_free(vtk_grid->cells);
        sb_free(vtk_grid->values);
        sb_free(vtk_grid->points);
    }

}

struct vtk_unstructured_grid *new_vtk_unstructured_grid_from_alg_grid(struct grid *grid, bool clip_with_plain,
                                                                      float *plain_coordinates, bool clip_with_bounds,
                                                                      float *bounds) {

    if(grid == NULL)
        return NULL;

    struct vtk_unstructured_grid *vtk_grid = new_vtk_unstructured_grid();

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

    float center_x, center_y, center_z, half_face;
    double v;

    struct point_3d aux1;
    struct point_3d aux2;
    struct point_3d aux3;
    struct point_3d aux4;
    struct point_3d aux5;
    struct point_3d aux6;
    struct point_3d aux7;
    struct point_3d aux8;

    int id = 0;
    int num_cells = 0;

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

            v = grid_cell->v;
            half_face = grid_cell->half_face_length;

            sb_push(vtk_grid->values, v);

            aux1.x = center_x - half_face;
            aux1.y = center_y - half_face;
            aux1.z = center_z - half_face;

            aux2.x = center_x + half_face;
            aux2.y = center_y - half_face;
            aux2.z = center_z - half_face;

            aux3.x = center_x + half_face;
            aux3.y = center_y + half_face;
            aux3.z = center_z - half_face;

            aux4.x = center_x - half_face;
            aux4.y = center_y + half_face;
            aux4.z = center_z - half_face;

            aux5.x = center_x - half_face;
            aux5.y = center_y - half_face;
            aux5.z = center_z + half_face;

            aux6.x = center_x + half_face;
            aux6.y = center_y - half_face;
            aux6.z = center_z + half_face;

            aux7.x = center_x + half_face;
            aux7.y = center_y + half_face;
            aux7.z = center_z + half_face;

            aux8.x = center_x - half_face;
            aux8.y = center_y + half_face;
            aux8.z = center_z + half_face;

            if(point_hash_search(hash, aux1) == -1) {
                sb_push(vtk_grid->points, aux1);
                point_hash_insert(hash, aux1, id);
                id++;
            }

            if(point_hash_search(hash, aux2) == -1) {
                sb_push(vtk_grid->points, aux2);
                point_hash_insert(hash, aux2, id);
                id++;
            }

            if(point_hash_search(hash, aux3) == -1) {
                point_hash_insert(hash, aux3, id);
                sb_push(vtk_grid->points, aux3);
                id++;
            }

            if(point_hash_search(hash, aux4) == -1) {
                point_hash_insert(hash, aux4, id);
                sb_push(vtk_grid->points, aux4);
                id++;
            }

            if(point_hash_search(hash, aux5) == -1) {
                sb_push(vtk_grid->points, aux5);
                point_hash_insert(hash, aux5, id);
                id++;
            }

            if(point_hash_search(hash, aux6) == -1) {
                sb_push(vtk_grid->points, aux6);
                point_hash_insert(hash, aux6, id);
                id++;
            }

            if(point_hash_search(hash, aux7) == -1) {
                sb_push(vtk_grid->points, aux7);
                point_hash_insert(hash, aux7, id);
                id++;
            }

            if(point_hash_search(hash, aux8) == -1) {
                sb_push(vtk_grid->points, aux8);
                point_hash_insert(hash, aux8, id);
                id++;
            }

            sb_push(vtk_grid->cells, point_hash_search(hash, aux1));
            sb_push(vtk_grid->cells, point_hash_search(hash, aux2));
            sb_push(vtk_grid->cells, point_hash_search(hash, aux3));
            sb_push(vtk_grid->cells, point_hash_search(hash, aux4));
            sb_push(vtk_grid->cells, point_hash_search(hash, aux5));
            sb_push(vtk_grid->cells, point_hash_search(hash, aux6));
            sb_push(vtk_grid->cells, point_hash_search(hash, aux7));
            sb_push(vtk_grid->cells, point_hash_search(hash, aux8));
            num_cells++;
        }

        grid_cell = grid_cell->next;
    }

    vtk_grid->num_cells = num_cells;
    vtk_grid->num_points = id;

    return vtk_grid;
}

void save_vtk_unstructured_grid_as_vtu(struct vtk_unstructured_grid *vtk_grid, char *filename, bool binary, bool compress, int level) {

    sds *file_headers = NULL;
    sds file_content = sdsempty();

    sb_push(file_headers, sdsempty());

    file_content = sdscat(
        file_content,
        "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");

    file_content = sdscat(file_content, "  <UnstructuredGrid>\n");

    file_content = sdscatprintf(file_content, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",
                                vtk_grid->num_points, vtk_grid->num_cells);

    file_content = sdscat(file_content, "      <CellData Scalars=\"Scalars_\">\n");

    int offset = 0;

    if(binary) {
        file_content = sdscatprintf(file_content, "        <DataArray type=\"Float32\" Name=\"Scalars_\" format=\"appended\" offset=\"%d\">\n", offset);
    }
    else {
        file_content = sdscat(file_content, "        <DataArray type=\"Float32\" Name=\"Scalars_\" format=\"ascii\">\n");
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
        file_content = sdscatprintf(
            file_content,
            "        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\">\n", offset);
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

    offset += (vtk_grid->num_points * 4 * 3) + 8; //3*32 bits float for each point

    if(binary) {
        file_content =
            sdscatprintf(file_content, "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"appended\" offset=\"%d\"/>\n", offset);
    } else {
        file_content =
            sdscat(file_content, "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    }

    int points_per_cell = 8;
    int cell_type = 12;

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

    if(!binary)
        file_content = sdscat(file_content, "        </DataArray>\n");

    offset += (vtk_grid->num_cells * 8 * points_per_cell) + 8; //64 bits

    if(binary) {
        file_content =
            sdscatprintf(file_content, "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"appended\" offset=\"%d\"/>\n", offset);
    } else {
        file_content = sdscat(file_content, "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    }

    {
        int offset_local = 8;

        if (!binary) {
            for (int i = 0; i < num_cells; i++) {
                file_content = sdscat(file_content, "     ");
                file_content = sdscatprintf(file_content, "%d ", offset_local);
                offset_local += 8;
                file_content = sdscat(file_content, "\n");
            }
        }
    }

    if(!binary)
        file_content = sdscat(file_content, "        </DataArray>\n");

    offset += (vtk_grid->num_cells * 8) + 8; //64 bits

    if(binary) {
        file_content = sdscatprintf(file_content, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\"%d\"/>\n", offset);
    } else {
        file_content = sdscat(file_content, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    }

    if(!binary) {

        for(int i = 0; i < num_cells; i++) {
            file_content = sdscatprintf(file_content, "    %d\n", cell_type);
        }
    }

    offset += (vtk_grid->num_cells) + 8; //1 byte

    if(!binary)
        file_content = sdscat(file_content, "        </DataArray>\n");

    file_content = sdscat(file_content, "      </Cells>\n");

    file_content = sdscat(file_content, "    </Piece>\n");
    file_content = sdscat(file_content, "  </UnstructuredGrid>\n");

    size_t size_until_now = 0;

    if(binary) {
        file_content = sdscat(file_content, "  <AppendedData encoding=\"raw\">\n   _");

        size_until_now = sdslen(file_content);

        //scalars
        uint64_t block_size = sizeof(float)*vtk_grid->num_cells;
        file_content = sdscatlen(file_content, &block_size, sizeof(uint64_t));
        file_content = sdscatlen(file_content, vtk_grid->values, (size_t) block_size);
        size_until_now += (sizeof(float)*vtk_grid->num_cells + sizeof(uint64_t));


        //Points
        block_size = sizeof(struct point_3d)*vtk_grid->num_points;
        file_content = sdscatlen(file_content, &block_size, sizeof(uint64_t));
        file_content = sdscatlen(file_content, vtk_grid->points, (size_t) block_size);
        size_until_now += (sizeof(struct point_3d)*vtk_grid->num_points + sizeof(uint64_t));

        //connectivity
        block_size = num_cells*points_per_cell*sizeof(int64_t);
        file_content = sdscatlen(file_content, &block_size, sizeof(uint64_t));
        size_until_now += sizeof(uint64_t);

        for(int i = 0; i < num_cells; i++) {
            for(int j = 0; j < points_per_cell; j++) {
                int64_t aux = (int64_t ) vtk_grid->cells[points_per_cell * i + j];
                file_content = sdscatlen(file_content, &aux, sizeof(int64_t));
                size_until_now += sizeof(int64_t);
            }
        }

        //offsets
        block_size = num_cells*sizeof(int64_t);
        file_content = sdscatlen(file_content, &block_size, sizeof(uint64_t));
        size_until_now += sizeof(uint64_t);

        int64_t offset_local = 8;
        for(int i = 0; i < num_cells; i++) {
            file_content = sdscatlen(file_content, &offset_local, sizeof(int64_t));
            offset_local += 8;
            size_until_now += sizeof(int64_t);
        }

        //types
        block_size = num_cells*sizeof(uint8_t);
        file_content = sdscatlen(file_content, &block_size, sizeof(uint64_t));
        size_until_now += sizeof(uint64_t);

        for(int i = 0; i < num_cells; i++) {
            uint8_t aux = (uint8_t )cell_type;
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
    }
    else {
        output_file = fopen(filename, "w");
        fprintf(output_file, "%s", file_content);
    }

    sdsfree(file_content);
    fclose(output_file);
}

void save_vtk_unstructured_grid_as_legacy_vtk(struct vtk_unstructured_grid *vtk_grid, char *filename, bool binary) {

    sds file_content = sdsempty();

    file_content = sdscat(file_content,"# vtk DataFile Version 4.2\n");
    file_content = sdscat(file_content,"vtk output\n");

    if(binary) {
        file_content = sdscat(file_content,"BINARY\n");
    } else {
        file_content = sdscat(file_content,"ASCII\n");
    }

    file_content = sdscat(file_content,"DATASET UNSTRUCTURED_GRID\n");
    file_content = sdscatprintf(file_content,"POINTS %d float\n", vtk_grid->num_points);

    size_t size_until_now = sdslen(file_content);

    int num_points = sb_count(vtk_grid->points);
    for(int i = 0; i < num_points; i++) {
        struct point_3d p = vtk_grid->points[i];
        if(binary) {
            file_content = write_binary_point(file_content, &p);
            size_until_now += 3*sizeof(int);
        }
        else {
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
    }
    else {
        fprintf(output_file, "%s", file_content);
    }

    sdsfree(file_content);
    fclose(output_file);
}