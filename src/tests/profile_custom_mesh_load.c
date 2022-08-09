//
// Created by sachetto on 07/04/2020.
//

#include "../alg/grid/grid.h"
#include "../config/domain_config.h"
#include "../config/save_mesh_config.h"
#include "../3dparty/ini_parser/ini.h"
#include "../3dparty/sds/sds.h"
#include "../3dparty/stb_ds.h"
#include "../utils/stop_watch.h"
#include "../vtk_utils/vtk_unstructured_grid.h"

#include <gdbm.h>

struct elapsed_times {
    double config_time;
    double create_grid_time;
    double order_grid_time;
    double init_time;
    double save_time;
    double end_time;
    double clean_time;
    double total_time;
    double load_mesh_time;
} __attribute__((packed));

//TODO: check for memory leaks with valgrind
int profile_custom_mesh_load(char *discretization, struct elapsed_times *times) {

    set_no_stdout(true);

    struct grid *grid = new_grid();
    struct config *domain_config;

    domain_config = alloc_and_init_config_data();

    shput_dup_value(domain_config->config_data, "maximum_discretization", discretization);
    shput_dup_value(domain_config->config_data, "mesh_file", "meshes/rabheart.alg");

    domain_config->main_function_name = strdup("initialize_grid_with_rabbit_mesh");
    shput_dup_value(domain_config->config_data, "name", "Test custom mesh");

    shput(domain_config->config_data, "side_length_x", strdup(discretization));
    shput(domain_config->config_data, "side_length_y", strdup(discretization));
    shput(domain_config->config_data, "side_length_z", strdup(discretization));

    struct stop_watch config_time;
    start_stop_watch(&config_time);
    init_config_functions(domain_config, "./shared_libs/libdefault_domains.so", "domain");
    times->config_time = stop_stop_watch(&config_time);

    struct stop_watch create_grid_time;
    start_stop_watch(&create_grid_time);
    int success = ((set_spatial_domain_fn*)domain_config->main_function)(domain_config, grid);
    times->create_grid_time = stop_stop_watch(&create_grid_time);

    if(!success) {
        clean_and_free_grid(grid);
        free_config_data(domain_config);
        return 0;
    }

    struct stop_watch order_grid_time;
    start_stop_watch(&order_grid_time);
    order_grid_cells(grid);
    times->order_grid_time = stop_stop_watch(&order_grid_time);

    struct config *save_mesh_config = alloc_and_init_config_data();

    save_mesh_config->init_function_name = strdup("init_save_as_vtk_or_vtu");
    save_mesh_config->main_function_name = strdup("save_as_vtu");
    save_mesh_config->end_function_name = strdup("end_save_as_vtk_or_vtu");
    shput_dup_value(save_mesh_config->config_data, "output_dir", "/tmp");
    shput_dup_value(save_mesh_config->config_data, "print_rate", "1");

    sds file_prefix = sdscatprintf(sdsempty(), "test_custom_mesh_%s", discretization);

    init_config_functions(save_mesh_config, "shared_libs/libdefault_save_mesh.so", "save_result");

    shput(save_mesh_config->config_data, "file_prefix", strdup(file_prefix));
    shput(save_mesh_config->config_data, "compress", strdup("yes"));
    shput(save_mesh_config->config_data, "save_pvd", strdup("no"));

    struct time_info ti = ZERO_TIME_INFO;

    struct stop_watch init_time;
    start_stop_watch(&init_time);
    ((init_save_mesh_fn *)save_mesh_config->init_function)(save_mesh_config);
    times->init_time = stop_stop_watch(&init_time);

    struct stop_watch save_time;
    start_stop_watch(&save_time);
    ((save_mesh_fn *)save_mesh_config->main_function)(&ti, save_mesh_config, grid, NULL, NULL);
    times->save_time = stop_stop_watch(&save_time);

    struct stop_watch end_time;
    start_stop_watch(&end_time);
    ((end_save_mesh_fn *)save_mesh_config->end_function)(save_mesh_config, grid);
    times->end_time = stop_stop_watch(&end_time);

    file_prefix = sdscat(file_prefix, "_it_0.vtu");
    
    sds tmp_mesh = sdsnew("/tmp/");
    tmp_mesh = sdscat(tmp_mesh, file_prefix);    

    struct stop_watch load_time;
    start_stop_watch(&load_time);    
    struct vtk_unstructured_grid *loaded_mesh = new_vtk_unstructured_grid_from_file(tmp_mesh, true);
    times->load_mesh_time= stop_stop_watch(&load_time);
    
    sdsfree(file_prefix);
    sdsfree(tmp_mesh);

    struct stop_watch clean_time;
    start_stop_watch(&clean_time);
    clean_and_free_grid(grid);
    times->clean_time = stop_stop_watch(&clean_time);

    free_vtk_unstructured_grid(loaded_mesh);
    free_config_data(domain_config);

    return 1;

}

int main(int argc, char **argv) {

    const long long nruns = 5;

    struct elapsed_times times;
    struct elapsed_times average_times = { 0 };

    if(argc != 2) {
        printf("Usage: %s hardware_key", argv[0]);
        exit(EXIT_FAILURE);
    }

    sds hash_key_with_size = sdsnew(argv[1]);
    sds nruns_string = sdsfromlonglong(nruns);
    hash_key_with_size = sdscat(hash_key_with_size, nruns_string);
    sdsfree(nruns_string);

    datum hash_key;
    hash_key.dptr = (char*)hash_key_with_size;
    hash_key.dsize = strlen(hash_key.dptr);

    for(int i = 0; i < nruns; i++) {
        profile_custom_mesh_load("500", &times);

        average_times.config_time      += times.config_time;
        average_times.create_grid_time += times.create_grid_time;
        average_times.order_grid_time  += times.order_grid_time;
        average_times.init_time        += times.init_time;
        average_times.save_time        += times.save_time;
        average_times.end_time         += times.end_time;
        average_times.clean_time       += times.clean_time;
        average_times.load_mesh_time   += times.load_mesh_time;

    }

    average_times.config_time      /= nruns;
    average_times.create_grid_time /= nruns;
    average_times.order_grid_time  /= nruns;
    average_times.init_time        /= nruns;
    average_times.save_time        /= nruns;
    average_times.end_time         /= nruns;
    average_times.clean_time       /= nruns;
    average_times.load_mesh_time   /= nruns;

    average_times.total_time = average_times.config_time + average_times.create_grid_time +
                               average_times.order_grid_time + average_times.init_time +
                               average_times.save_time + average_times.end_time +
                               average_times.clean_time + average_times.load_mesh_time;

    GDBM_FILE f;

    f = gdbm_open( "./tests_bin/profile_custom_mesh_times.gdbm", 4096, GDBM_WRCREAT, 0644, NULL );

    datum content = gdbm_fetch (f, hash_key);

    datum data;
    data.dptr = (char*)(&average_times);
    data.dsize = sizeof(average_times);

    printf("CURRENT RUN\n");
    printf("Avg Config time: %lf μs\n", average_times.config_time);
    printf("Avg Create grid time: %lf μs\n", average_times.create_grid_time);
    printf("Avg Order grid time: %lf μs\n", average_times.order_grid_time);
    printf("Avg Init save function time: μs %lf\n", average_times.init_time);
    printf("Avg Save function time: %lf μs\n", average_times.save_time);
    printf("Avg End save function time: %lf μs\n", average_times.end_time);
    printf("Avg Load mesh from file time: %lf μs\n", average_times.load_mesh_time);
    printf("Avg Clean grid time: %lf μs\n", average_times.clean_time);
    printf("Avg Total time: %lf μs\n", average_times.total_time);

    printf("---------------------------------------------------\n");

    if (content.dptr == NULL) {
        printf("\nFirst run in this hardware with nruns = %lld\n", nruns);
        printf("---------------------------------------------------\n");
        gdbm_store(f, hash_key, data, GDBM_INSERT);
    }
    else {
        printf("BEST RUN\n");
        struct elapsed_times *best_run = (struct elapsed_times *)content.dptr;
        printf("Avg Config time: %lf μs\n", best_run->config_time);
        printf("Avg Create grid time: %lf μs\n", best_run->create_grid_time);
        printf("Avg Order grid time: %lf μs\n", best_run->order_grid_time);
        printf("Avg Init save function time: %lf μs\n", best_run->init_time);
        printf("Avg Save function time: %lf μs\n", best_run->save_time);
        printf("Avg End save function time: %lf μs\n", best_run->end_time);
        printf("Avg Load mesh from file time: %lf μs\n", average_times.load_mesh_time);
        printf("Avg Clean grid time: %lf μs\n", best_run->clean_time);
        printf("Avg Total time: %lf μs\n", best_run->total_time);

        printf("---------------------------------------------------\n");

        double speedup = best_run->total_time/average_times.total_time;

        //10% speedup
        if(speedup > 1.0) {
            printf("Current run is %lfx faster than best run. Replacing record.\n", speedup);
            gdbm_store(f, hash_key, data, GDBM_REPLACE);
        }
        else if(speedup < 1.0) {
            printf("Current run is %lfx slower than best run.\n", speedup);
        }

        free(best_run);
    }

    sdsfree(hash_key_with_size);
    gdbm_close(f);
}
