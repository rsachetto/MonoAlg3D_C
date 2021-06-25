//
// Created by sachetto on 07/04/2020.
//

#include "../3dparty/sds/sds.h"
#include "../3dparty/stb_ds.h"
#include "../logger/logger.h"
#include "../utils/stop_watch.h"
#include "common.h"
#include <gdbm.h>

struct elapsed_times {
    double config_time;
    double simulation_time;
    double total_time;
} __attribute__((packed));

int main(int argc, char **argv) {

    const long long nruns = 1;

    struct elapsed_times times;

    if(argc != 2) {
        printf("Usage: %s hardware_key", argv[0]);
        exit(EXIT_FAILURE);
    }

    set_no_stdout(true);

    struct stop_watch config_time;
    start_stop_watch(&config_time);

    char *out_dir_gpu_no_precond = "tests_bin/circle_cg_gpu_no_precond";

    struct user_options *options = load_options_from_file("example_configs/plain_mesh_with_fibrosis_and_border_zone_inside_circle_example_2cm.ini");
    options->final_time = 100.0;

    free(options->save_mesh_config->main_function_name);

    options->save_mesh_config->main_function_name = strdup("save_as_text_or_binary");
    shput_dup_value(options->save_mesh_config->config_data, "print_rate", "50");
    shput_dup_value(options->save_mesh_config->config_data, "file_prefix", "V");

    shput_dup_value(options->domain_config->config_data, strdup("seed"), "150");

    free(options->linear_system_solver_config->main_function_name);
    free(options->linear_system_solver_config->init_function_name);
    free(options->linear_system_solver_config->end_function_name);

    options->linear_system_solver_config->main_function_name = strdup("conjugate_gradient");
    options->linear_system_solver_config->init_function_name = strdup("init_conjugate_gradient");
    options->linear_system_solver_config->end_function_name = strdup("end_conjugate_gradient");

    shput_dup_value(options->linear_system_solver_config->config_data, "use_preconditioner", "false");

#ifdef COMPILE_CUDA
    //solve odes also in the GPU
    shput_dup_value(options->linear_system_solver_config->config_data, "use_gpu", "true");
    options->gpu = true;
#else
    shput_dup_value(options->linear_system_solver_config->config_data, "use_gpu", "false");
    options->gpu = false;
#endif

    times.config_time = stop_stop_watch(&config_time);

    struct stop_watch simulation_time;
    start_stop_watch(&simulation_time);

    int success = 1;
    success = run_simulation_with_config(options, out_dir_gpu_no_precond);
    times.simulation_time = stop_stop_watch(&simulation_time);

    free_user_options(options);

    if(!success) {
        fprintf(stderr, "Error running simulation!\n");
        return EXIT_FAILURE;
    }

    sds hash_key_with_size = sdsnew(argv[1]);
    sds nruns_string = sdsfromlonglong(nruns);
    hash_key_with_size = sdscat(hash_key_with_size, nruns_string);
    sdsfree(nruns_string);

    datum hash_key;
    hash_key.dptr = (char*)hash_key_with_size;
    hash_key.dsize = strlen(hash_key.dptr);

    times.total_time = times.config_time + times.simulation_time;

    GDBM_FILE f;

    f = gdbm_open( "./tests_bin/profile_all_gpu_simulation_times.gdbm", 4096, GDBM_WRCREAT, 0644, NULL );

    datum content = gdbm_fetch (f, hash_key);

    datum data;
    data.dptr = (char*)(&times);
    data.dsize = sizeof(times);

    printf("CURRENT RUN\n");
    printf("Config time: %lf μs\n",     times.config_time);
    printf("Simulation time: %lf μs\n", times.simulation_time);
    printf("Total time: %lf μs\n",  times.total_time);

    printf("---------------------------------------------------\n");

    if (content.dptr == NULL) {
        printf("\nFirst run in this hardware with nruns = %lld\n", nruns);
        printf("---------------------------------------------------\n");
        gdbm_store(f, hash_key, data, GDBM_INSERT);
    }
    else {
        printf("BEST RUN\n");
        struct elapsed_times *best_run = (struct elapsed_times *)content.dptr;
        printf("Config time: %lf μs\n",     best_run->config_time);
        printf("Simulation time: %lf μs\n", best_run->config_time);
        printf("Total time: %lf μs\n",      best_run->total_time);
        printf("---------------------------------------------------\n");

        double speedup = best_run->total_time/times.total_time;

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
