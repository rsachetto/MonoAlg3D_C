//
// Created by sachetto on 07/04/2020.
//

#include "../3dparty/ini_parser/ini.h"
#include "../3dparty/sds/sds.h"
#include "../3dparty/stb_ds.h"
#include "../alg/grid/grid.h"
#include "../config/domain_config.h"
#include "../config/save_mesh_config.h"
#include "../utils/stop_watch.h"
#include "../vtk_utils/vtk_unstructured_grid.h"

#include <gdbm.h>

struct elapsed_times {
    double total_time;
} __attribute__((packed));

int main(int argc, char **argv) {

    struct elapsed_times average_times = { 0 };

    if(argc < 2 || argc > 3) {
        printf("Usage: %s hardware_key or %s hardware_key n_runs", argv[0], argv[0]);
        exit(EXIT_FAILURE);
    }

    sds hash_key_with_size = sdsnew(argv[1]);

    long long nruns = 30;

    if(argc == 3) {
        nruns = strtol(argv[2], NULL, 10);
    }

    sds nruns_string = sdsfromlonglong(nruns);
    hash_key_with_size = sdscat(hash_key_with_size, nruns_string);
    sdsfree(nruns_string);

    datum hash_key;
    hash_key.dptr = (char*)hash_key_with_size;
    hash_key.dsize = strlen(hash_key.dptr);

    for(int i = 0; i < nruns; i++) {

        struct stop_watch total_time;
        start_stop_watch(&total_time);
        struct vtk_unstructured_grid *vtk_grid = new_vtk_unstructured_grid_from_file("tests_bin/profile_mesh_load.geo", true);
        set_vtk_grid_visibility(&vtk_grid);
        average_times.total_time += stop_stop_watch(&total_time);

        free_vtk_unstructured_grid(vtk_grid);

    }

    average_times.total_time /= nruns;

    GDBM_FILE f;

    f = gdbm_open( "./tests_bin/profile_en_load_times.gdbm", 4096, GDBM_WRCREAT, 0644, NULL );

    datum content = gdbm_fetch (f, hash_key);

    datum data;
    data.dptr = (char*)(&average_times);
    data.dsize = sizeof(average_times);

    printf("CURRENT RUN\n");
    printf("Avg Total time: %lf μs (%lf ms)\n", average_times.total_time, average_times.total_time/1000.0);

    printf("---------------------------------------------------\n");

    if (content.dptr == NULL) {
        printf("\nFirst run in this hardware with nruns = %lld\n", nruns);
        printf("---------------------------------------------------\n");
        gdbm_store(f, hash_key, data, GDBM_INSERT);
    }
    else {
        printf("BEST RUN\n");
        struct elapsed_times *best_run = (struct elapsed_times *)content.dptr;
        printf("Avg Total time: %lf μs (%lf ms) \n", best_run->total_time, best_run->total_time/1000.0);

        printf("---------------------------------------------------\n");

        double speedup = best_run->total_time/average_times.total_time;

        //10% speedup
        if(speedup > 1.0) {
            printf("Current run is %lfx faster than best run. Replacing record.\n", speedup);
            gdbm_store(f, hash_key, data, GDBM_REPLACE);
        }
        else if(speedup < 1.0) {
            printf("Current run is %lfx slower than best run.\n", 1.0/speedup);
        }

        free(best_run);
    }

    sdsfree(hash_key_with_size);
    gdbm_close(f);
}
