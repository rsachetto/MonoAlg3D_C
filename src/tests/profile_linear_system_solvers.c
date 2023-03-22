//
// Created by sachetto on 07/04/2020.
//

#include "../alg/grid/grid.h"
#include "../config/domain_config.h"
#include "../config/save_mesh_config.h"
#include "../config/config_parser.h"
#include "../utils/file_utils.h"
#include "../3dparty/ini_parser/ini.h"
#include "../3dparty/sds/sds.h"

#include "../3dparty/stb_ds.h"
#include "../config_helpers/config_helpers.h"
#include "../utils/stop_watch.h"
#include "../logger/logger.h"

#include <gdbm.h>

struct elapsed_times {
    double init_time;
    double run_time;
    double end_time;
    double total_time;

} __attribute__((packed));

#include "profile_solvers_common.c"

int main(int argc, char **argv) {

    const long long nruns = 100;

    struct elapsed_times times;
    struct elapsed_times average_times = { 0 };

    if(argc != 2) {
        printf("Usage: ./%s hardware_key", argv[0]);
        exit(EXIT_FAILURE);
    }

    sds hash_key_with_size = sdsnew(argv[1]);
    sds nruns_string = sdsfromlonglong(nruns);
    hash_key_with_size = sdscat(hash_key_with_size, nruns_string);
    sdsfree(nruns_string);

    datum hash_key;
    hash_key.dptr = (char*)hash_key_with_size;
    hash_key.dsize = (int)strlen(hash_key.dptr);

    FILE *A = NULL;
    FILE *B = NULL;
    FILE *X = NULL;

    A = fopen("tests_bin/A4.txt", "r");
    B = fopen("tests_bin/B4.txt", "r");
    X = fopen("tests_bin/X4.txt", "r");

    assert(A);
    assert(B);
    assert(X);

    struct grid *grid = new_grid();
    assert (grid);

    construct_grid_from_file(grid, A, B);

    for(int i = 0; i < nruns; i++) {

        profile_solver(true, "cpu_conjugate_gradient", "init_cpu_conjugate_gradient", NULL, grid, 1, &times);

        average_times.init_time += times.init_time;
        average_times.run_time  += times.run_time;
        average_times.end_time  += times.end_time;

    }

    clean_and_free_grid(grid);
    fclose(A);
    fclose(B);
    fclose(X);

    average_times.init_time /= nruns;
    average_times.run_time  /= nruns;
    average_times.end_time  /= nruns;

    average_times.total_time = average_times.init_time + average_times.run_time + average_times.end_time;

    GDBM_FILE f;

    f = gdbm_open( "./tests_bin/profile_solver_times.gdbm", 4096, GDBM_WRCREAT, 0644, NULL );

    datum content = gdbm_fetch (f, hash_key);

    datum data;
    data.dptr = (char*)(&average_times);
    data.dsize = sizeof(average_times);

    printf("CURRENT RUN\n");
    printf("Avg Init function time: %lf us\n", average_times.init_time);
    printf("Avg Run function time: %lf us\n", average_times.run_time);
    printf("Avg End function time: %lf us\n", average_times.end_time);
    printf("Avg Total time: %lf us\n", average_times.total_time);

    printf("---------------------------------------------------\n");

    if (content.dptr == NULL) {
        printf("\nFirst run in this hardware with nruns = %lld\n", nruns);
        printf("---------------------------------------------------\n");
        gdbm_store(f, hash_key, data, GDBM_INSERT);
    }
    else {
        printf("BEST RUN\n");
        struct elapsed_times *best_run = (struct elapsed_times *)content.dptr;
        printf("Avg Init function time: %lf us\n", best_run->init_time);
        printf("Avg Run function time: %lf us\n", best_run->run_time);
        printf("Avg End function time: %lf us\n", best_run->end_time);
        printf("Avg Total time: %lf us\n", best_run->total_time);

        printf("---------------------------------------------------\n");

        double speedup = (double)best_run->total_time/(double)average_times.total_time;

        if(speedup > 1.0) {
            printf("Current run is %lf x faster than best run. Replacing record.\n", speedup);
            gdbm_store(f, hash_key, data, GDBM_REPLACE);
        }
        else if(speedup < 1.0) {
            printf("Current run is %lf x slower than best run.\n", 1.0/speedup);
        }

        free(best_run);
    }

    sdsfree(hash_key_with_size);
    gdbm_close(f);
}
