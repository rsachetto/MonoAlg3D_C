//
// Created by sachetto on 01/11/18.
//

#ifndef MONOALG3D_DATA_UTILS_H
#define MONOALG3D_DATA_UTILS_H

#include "../3dparty/miniz/miniz.h"
#include "../3dparty/sds/sds.h"
#include "../common_types/common_types.h"
#include <stdbool.h>


#define NUMBER_OF_POINTS     "NumberOfPoints"
#define NUMBER_OF_POINTS_LEN 14

#define NUMBER_OF_CELLS      "NumberOfCells"
#define NUMBER_OF_CELLS_LEN  13

#define POINTS           "Points"
#define POINTS_LEN        6

#define CELLS            "Cells"
#define CELLS_LEN         5

#define CELL_TYPES       "Cell_Types"
#define CELL_TYPES_LEN   10

#define CELL_DATA        "Cell_Data"
#define CELL_DATA_LEN    9

#define DATAARRAY        "DataArray"
#define DATAARRAY_LEN    9

#define SCALARS          "SCALARS"
#define SCALARS_LEN       7

#define OFFSET           "offset"
#define OFFSET_LEN       6

#define NAME             "Name"
#define NAME_LEN         4

#define SCALARS_NAME     "Scalars_"
#define SCALARS_NAME_LEN 8

#define OFFSETS          "offsets"
#define OFFSETS_LEN      7

#define CONNECTIVITY      "connectivity"
#define CONNECTIVITY_LEN  12

#define TYPES            "types"
#define TYPES_LEN        5

#define TYPE             "type"
#define TYPE_LEN         4

#define APPENDEDDATA     "AppendedData"
#define APPENDEDDATA_LEN 12

#define ENCODING         "encoding"
#define ENCODING_LEN     8


#define HEADER_TYPE      "header_type"
#define HEADER_TYPE_LEN  11

#define COMPRESSOR       "compressor"
#define COMPRESSOR_LEN   10

#define FORMAT           "format"
#define FORMAT_LEN       6

#define ASCII            "ascii"
#define ASCII_LEN        5

#define LOOKUP_TABLE     "LOOKUP_TABLE"
#define LOOKUP_TABLE_LEN 12

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
    bool can_handle;
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
void get_data_block_from_compressed_vtu_file(const char *raw_data, void* values, uint64_t num_blocks, uint64_t block_size_uncompressed, uint64_t last_block_size, uint64_t  *block_sizes_compressed);

void get_data_block_from_uncompressed_binary_vtu_file(char *raw_data, void* values, size_t header_size);
void read_binary_point(void *source, struct point_3d *p);

#endif // MONOALG3D_DATA_UTILS_H
