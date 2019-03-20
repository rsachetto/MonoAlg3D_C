//
// Created by bergolho on 25/01/19.
//

#include "vtk_polydata_grid.h"
#include "../alg/cell/cell.h"
#include "../string/sds.h"
#include "data_utils.h"
#include <inttypes.h>
#include <math.h>
#include <stdint.h>

#include "../common_types/common_types.h"
#include "../single_file_libraries/stb_ds.h"

struct vtk_polydata_grid *new_vtk_polydata_grid ()
{
    struct vtk_polydata_grid *grid = (struct vtk_polydata_grid *)malloc(sizeof(struct vtk_polydata_grid));

    grid->values = NULL;
    grid->points = NULL;
    grid->lines = NULL;

    grid->num_points = 0;
    grid->num_lines = 0;

    return grid;
}

sds create_common_vtp_header(bool compressed, int num_points, int num_lines) 
{

    sds header = sdsempty();

    if(compressed) 
    {
        header = sdscat(header, "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" "
                                "header_type=\"UInt64\" compressor=\"vtkZLibDataCompressor\">\n");
    } 
    else 
    {
        header = sdscat(
            header,
            "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
    }

    header = sdscat(header, "  <PolyData>\n");

    header = sdscatprintf(header, "    <Piece NumberOfPoints=\"%d\" NumberOfLines=\"%d\">\n", num_points, num_lines);

    header = sdscat(header, "      <PointData Scalars=\"Scalars_\">\n");

    return header;
}

void new_vtk_polydata_grid_from_purkinje_grid(struct vtk_polydata_grid **vtk_grid, struct grid *grid, bool clip_with_plain,
                                                                     real_cpu *plain_coordinates, bool clip_with_bounds,
                                                                     real_cpu *bounds, bool read_only_values)
{
    static bool mesh_already_loaded =  false;

    if(grid == NULL) 
    {
        return;
    }

    if(!read_only_values)
    {
        *vtk_grid = new_vtk_polydata_grid();
    }
    else
    {
        if(!(*vtk_grid) && mesh_already_loaded)
        {
            fprintf(stderr,
                    "Function new_vtk_polydata_grid_from_purkinje_grid can only be called with read_only_values if the grid is already loaded");
            exit(EXIT_FAILURE);
        }

        if(mesh_already_loaded)
        {
            assert(*vtk_grid);
            arrfree((*vtk_grid)->values);
            (*vtk_grid)->values = NULL;
        }
        else
        {
            *vtk_grid = new_vtk_polydata_grid();
        }
    }

    struct cell_node *grid_cell = grid->first_cell;
    struct node *u = grid->the_purkinje_network->list_nodes;

    struct point_3d aux;
    struct line auxl;
    real_cpu center_x, center_y, center_z;

    uint32_t id = 0;
//    uint32_t num_cells = 0;

    struct point_hash_entry *hash = NULL;
    hmdefault(hash, -1);

    while (grid_cell != NULL)
    {
        if (grid_cell->active)
        {
            center_x = grid_cell->center_x;
            center_y = grid_cell->center_y;
            center_z = grid_cell->center_z;

            arrput((*vtk_grid)->values, grid_cell->v);

            // This 'if' statement do not let us re-insert points and lines to the arrays ... =)
            if(mesh_already_loaded && read_only_values)
            {
                grid_cell = grid_cell->next;
                u = u->next;
                continue;
            }

            // Insert the point to the array of points
            aux.x = center_x;
            aux.y = center_y;
            aux.z = center_z;

            // Search for duplicates
            if(hmget(hash, aux) == -1)
            {
                arrput((*vtk_grid)->points, aux);
                hmput(hash, aux, id);
                id++;
            }

            // Insert the edge to the array of lines
            struct edge *v = u->list_edges;
            while (v != NULL)
            {
                auxl.source = u->id;
                auxl.destination = v->id;

                arrput((*vtk_grid)->lines, auxl);

                v = v->next;
            }

        }
        grid_cell = grid_cell->next;
        u = u->next;
    }

    if(!mesh_already_loaded)
    {
        (*vtk_grid)->num_points = id;
        (*vtk_grid)->num_lines = grid->the_purkinje_network->total_edges;

        if(read_only_values)
            mesh_already_loaded = true;
    }

}

void save_vtk_polydata_grid_as_legacy_vtk(struct vtk_polydata_grid *vtk_grid,\
                                        char *filename, bool binary)
{
    sds file_content = sdsempty();

    file_content = sdscat(file_content, "# vtk DataFile Version 4.2\n");
    file_content = sdscat(file_content, "vtk output\n");
    if(binary)
    {
        file_content = sdscat(file_content, "BINARY\n");
    }
    else
    {
        file_content = sdscat(file_content, "ASCII\n");
    }

    file_content = sdscat(file_content, "DATASET POLYDATA\n");
    file_content = sdscatprintf(file_content, "POINTS %d float\n", vtk_grid->num_points);

    size_t size_until_now = sdslen(file_content);

    int num_points = arrlen(vtk_grid->points);
    for(int i = 0; i < num_points; i++)
    {
        struct point_3d p = vtk_grid->points[i];

        if(binary)
        {
            file_content = write_binary_point(file_content, &p);
            size_until_now += 3 * sizeof(int);
        }
        else
        {
            file_content = sdscatprintf(file_content, "%lf %lf %lf\n", p.x, p.y, p.z);
        }
    }

    int num_lines = vtk_grid->num_lines;

    {
        sds tmp = sdscatprintf(sdsempty(), "\nLINES %d %d\n", num_lines, 3 * num_lines);

        size_until_now += sdslen(tmp);

        file_content = sdscatsds(file_content, tmp);
        sdsfree(tmp);
    }

    for (int i = 0; i < num_lines; i++)
    {
        struct line l = vtk_grid->lines[i];

        if (binary)
        {
            file_content = write_binary_line(file_content, &l);
            size_until_now += 3 * sizeof(int);
        }
        else
        {
            file_content = sdscatprintf(file_content, "2 %lu %lu\n",l.source,l.destination);
        }

    }

    int num_values = arrlen(vtk_grid->values);

    {
        sds tmp = sdscatprintf(sdsempty(), "POINT_DATA %d\n", num_values);
        tmp = sdscat(tmp, "SCALARS Vm float\n");
        tmp = sdscat(tmp, "LOOKUP_TABLE default\n");

        size_until_now += sdslen(tmp);

        file_content = sdscatsds(file_content, tmp);
        sdsfree(tmp);
    }

    for (int i = 0; i < num_values; i++)
    {
        if(binary)
        {
            int aux = invert_bytes(*((int *)&(vtk_grid->values[i])));
            file_content = sdscatlen(file_content, &aux, sizeof(int));
            size_until_now += sizeof(int);
        }
        else
        {
            file_content = sdscatprintf(file_content, "%lf ", vtk_grid->values[i]);
        }
    }

    FILE *output_file = NULL;

    if(binary)
    {
        output_file = fopen(filename, "wb");
        fwrite(file_content, size_until_now, 1, output_file);
    }
    else
    {
        output_file = fopen(filename, "w");
        fprintf(output_file, "%s", file_content);
    }

    sdsfree(file_content);
    fclose(output_file);

}

void save_vtk_polydata_grid_as_vtp (struct vtk_polydata_grid *vtk_grid, char *filename, bool binary)
{

    size_t offset = 0;

    sds file_content = create_common_vtp_header(false, vtk_grid->num_points, vtk_grid->num_lines);

    if(binary)
    {
        file_content = sdscat(
            file_content,
            "        <DataArray type=\"Float32\" Name=\"Scalars_\" format=\"appended\" offset=\"0\">\n"); // First
                                                                                                          // offset is
                                                                                                          // always 0

    }
    else
    {
        file_content =
            sdscat(file_content, "        <DataArray type=\"Float32\" Name=\"Scalars_\" format=\"ascii\">\n");
    }

    size_t num_values = arrlen(vtk_grid->values);

    if(!binary)
    {

        for(int i = 0; i < num_values; i++)
        {
            file_content = sdscatprintf(file_content, "     %lf ", vtk_grid->values[i]);
        }
    }

    file_content = sdscat(file_content, "        </DataArray>\n");
    file_content = sdscat(file_content, "      </PointData>\n");

    offset = (num_values * 4) + 8;

    file_content = sdscat(file_content, "      <Points>\n");

    if(binary)
    {
        file_content = sdscatprintf(file_content,
                                    "        <DataArray type=\"Float32\" Name=\"Points\" "
                                    "NumberOfComponents=\"3\" format=\"appended\" offset=\"%zu\">\n",
                                    offset);

    }
    else
    {
        file_content =
            sdscat(file_content,
                   "        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    }

    if(!binary)
    {
        int num_points = arrlen(vtk_grid->points);
        for(int i = 0; i < num_points; i++)
        {
            struct point_3d p = vtk_grid->points[i];
            file_content = sdscatprintf(file_content, "%lf %lf %lf\n", p.x, p.y, p.z);
        }
    }

    file_content = sdscat(file_content, "        </DataArray>\n");
    file_content = sdscat(file_content, "      </Points>\n");

    file_content = sdscat(file_content, "      <Lines>\n");

    offset += (vtk_grid->num_points * 4 * 3) + 8; // 3*32 bits float for each point

    if(binary)
    {
        file_content = sdscatprintf(
            file_content,
            "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"appended\" offset=\"%zu\">\n", offset);

    }
    else
    {
        file_content =
            sdscat(file_content, "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    }

    int num_lines = vtk_grid->num_lines;

    if(!binary)
    {
        for(int i = 0; i < num_lines; i++)
        {
            file_content = sdscatprintf(file_content,"%lu %lu\n",vtk_grid->lines[i].source,vtk_grid->lines[i].destination);
        }
    }

    file_content = sdscat(file_content, "        </DataArray>\n");

    offset += (vtk_grid->num_lines * 2 * 8) + 8; // 64 bits for the line index

    if(binary)
    {
        file_content = sdscatprintf(
            file_content, "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"appended\" offset=\"%zu\">\n",
            offset);
    }
    else
    {
        file_content = sdscat(file_content, "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    }

    {
        int offset_local = 2;

        if(!binary)
        {
            for(int i = 0; i < num_lines; i++)
            {
                file_content = sdscat(file_content, "     ");
                file_content = sdscatprintf(file_content, "%d ", offset_local);
                offset_local += 2;
                file_content = sdscat(file_content, "\n");
            }
        }
    }

    file_content = sdscat(file_content, "        </DataArray>\n");
    file_content = sdscat(file_content, "      </Lines>\n");


    offset += (vtk_grid->num_lines * 8) + 8; // 64 bits

    file_content = sdscat(file_content, "    </Piece>\n");
    file_content = sdscat(file_content, "  </PolyData>\n");

    size_t size_until_now = 0;

    if(binary)
    {
        file_content = sdscat(file_content, "  <AppendedData encoding=\"raw\">\n   _");

        size_until_now = sdslen(file_content);

        // scalars
        uint64_t block_size = sizeof(real_cpu) * vtk_grid->num_points;
        file_content = sdscatlen(file_content, &block_size, sizeof(uint64_t));
        file_content = sdscatlen(file_content, vtk_grid->values, (size_t)block_size);
        size_until_now += (sizeof(real_cpu) * vtk_grid->num_points + sizeof(uint64_t));

        // Points
        block_size = sizeof(struct point_3d) * vtk_grid->num_points;
        file_content = sdscatlen(file_content, &block_size, sizeof(uint64_t));
        file_content = sdscatlen(file_content, vtk_grid->points, (size_t)block_size);
        size_until_now += (sizeof(struct point_3d) * vtk_grid->num_points + sizeof(uint64_t));

        // connectivity
        block_size = vtk_grid->num_lines * 2 * sizeof(uint64_t);
        file_content = sdscatlen(file_content, &block_size, sizeof(uint64_t));
        for(int i = 0; i < num_lines; i++) 
        {
            uint64_t source = (uint64_t)vtk_grid->lines[i].source;
            uint64_t destination = (uint64_t)vtk_grid->lines[i].destination;
            file_content = sdscatlen(file_content, &source, sizeof(uint64_t));
            file_content = sdscatlen(file_content, &destination, sizeof(uint64_t));
            size_until_now += sizeof(uint64_t) * 2;
        }

        // offsets
        block_size = vtk_grid->num_lines * sizeof(int64_t);
        file_content = sdscatlen(file_content, &block_size, sizeof(uint64_t));
        size_until_now += sizeof(int64_t);

        int64_t offset_local = 2;
        for(int i = 0; i < vtk_grid->num_lines; i++) 
        {
            file_content = sdscatlen(file_content, &offset_local, sizeof(int64_t));
            offset_local += 2;
            size_until_now += sizeof(int64_t);
        }

        file_content = sdscat(file_content, "\n  </AppendedData>\n");
    }

    file_content = sdscat(file_content, "</VTKFile>\n");

    size_until_now += 29;

    FILE *output_file = NULL;

    if(binary) 
    {
        output_file = fopen(filename, "wb");
        fwrite(file_content, size_until_now, 1, output_file);
    } 
    else 
    {
        output_file = fopen(filename, "w");
        fprintf(output_file, "%s", file_content);
    }

    sdsfree(file_content);
    fclose(output_file);

}

void save_vtk_polydata_grid_as_vtp_compressed (struct vtk_polydata_grid *vtk_grid, char *filename, int compression_level)
{
    printf("\tIn 'save_vtk_polydata_grid_as_vtp_compressed'\n");

    printf("\tLeaving 'save_vtk_polydata_grid_as_vtp_compressed'\n");
    exit(EXIT_FAILURE);
}



void free_vtk_polydata_grid(struct vtk_polydata_grid *vtk_grid)
{
    if(vtk_grid) 
    {
        arrfree(vtk_grid->lines);
        arrfree(vtk_grid->values);
        arrfree(vtk_grid->points);
        free(vtk_grid);
    }
}