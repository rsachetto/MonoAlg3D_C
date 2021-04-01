//
// Created by sachetto on 01/11/18.
//

#ifndef MONOALG3D_DATA_UTILS_H
#define MONOALG3D_DATA_UTILS_H

#include "../3dparty/miniz/miniz.h"
#include "../3dparty/sds/sds.h"
#include "../common_types/common_types.h"
#include <stdbool.h>


#define NUMBER_OF_POINTS "NumberOfPoints"
#define NUMBER_OF_CELLS  "NumberOfCells"
#define POINTS           "Points"
#define CELLS            "Cells"
#define CELL_TYPES       "Cell_Types"
#define CELL_DATA        "Cell_Data"
#define DATAARRAY        "DataArray"
#define SCALARS          "SCALARS"
#define OFFSET           "offset"
#define NAME             "Name"
#define SCALARS_NAME     "Scalars_"
#define OFFSETS          "offsets"
#define CONNECTIVITY     "connectivity"
#define TYPES            "types"
#define TYPE             "type"
#define APPENDEDDATA     "AppendedData"
#define ENCODING         "encoding"
#define HEADER_TYPE      "header_type"
#define COMPRESSOR       "compressor"
#define FORMAT           "format"
#define ASCII            "ascii"
#define LOOKUP_TABLE     "LOOKUP_TABLE"

struct f32points {
    float x;
    float y;
    float z;
};


struct parser_state {
    char *number_of_points;
    char *number_of_cells;

    char *celldata_ofsset;
    char *points_ofsset;
    char *cells_connectivity_ofsset;
    char *cells_offsets_ofsset;
    char *cells_types_ofsset;
    char *name_value;
    char *format;
    char *point_data_type;
    char *celldata_ascii;
    char *points_ascii;
    char *cells_connectivity_ascii;
    char *encoding_type;
    char *header_type;
    char *base64_content;

    bool in_dataarray;
    bool in_points_array;
    bool compressed;
    bool binary;
    bool ascii;
};


int invert_bytes(int data);
sds write_binary_point(sds output_string, struct point_3d *p);
sds write_binary_line (sds output_string, struct line *l);

size_t uncompress_buffer(unsigned char const* compressed_data,
                         size_t compressed_size,
                         unsigned char* uncompressed_data,
                         size_t uncompressed_size);

static size_t compress_buffer(unsigned char const *uncompressed_data, size_t uncompressed_data_size,
                              unsigned char *compressed_data, size_t compression_space, int level);

void calculate_blocks_and_compress_data(size_t total_data_size_before_compression,
                                        size_t *total_data_size_after_compression, unsigned char *data_to_compress,
                                        unsigned char **compressed_buffer, size_t *num_blocks,
                                        size_t *block_size_uncompressed, size_t **block_sizes_compressed,
                                        size_t *last_block_size, int compression_level);

size_t get_block_sizes_from_compressed_vtu_file(char *raw_data, size_t header_size, uint64_t *num_blocks, uint64_t *block_size_uncompressed, uint64_t *last_block_size, uint64_t  **block_sizes_compressed);
void get_data_block_from_compressed_vtu_file(const char *raw_data, void* values, size_t header_size, uint64_t num_blocks, uint64_t block_size_uncompressed, uint64_t last_block_size, uint64_t  *block_sizes_compressed);

void get_data_block_from_uncompressed_binary_vtu_file(char *raw_data, void* values, size_t header_size);
void read_binary_point(void *source, struct point_3d *p);

#endif // MONOALG3D_DATA_UTILS_H
