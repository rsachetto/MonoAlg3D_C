# Author: Lucas Berg
# This script is responsible for generating the configuration files for the Purkinje coupling experiment
# of the 2020 Scientific Reports paper

import sys

def write_config_file (scenario,individual):
    if (scenario == 1):
        config_filename = "outputs/elnaz_purkinje_coupled_sc1.ini"
        output_dir = "./outputs/elnaz_purkinje_coupled_sc1"
        ode_library_file = "shared_libs/libtentusscher_epi_2004_S1.so"
    else:
        config_filename = "outputs/elnaz_purkinje_coupled_sc%u_%u.ini" % (scenario,individual)
        output_dir = "./outputs/elnaz_purkinje_coupled_sc%u_%u" % (scenario,individual)
        ode_library_file = "shared_libs/libtentusscher_epi_2004_S%u_%u.so" % (scenario,individual)

    file = open(config_filename,"w")

    # Write [main] section
    file.write("[main]\n")
    file.write("num_threads = 6\n")
    file.write("dt_pde = 0.02\n")
    file.write("simulation_time = 2000.0\n")
    file.write("abort_on_no_activity = false\n")
    file.write("use_adaptivity = false\n")
    file.write("calc_activation_time = true\n")
    file.write("print_apd_map = true\n")
    file.write("\n")

    # Write [update_monodomain] section
    file.write("[update_monodomain]\n")
    file.write("main_function = update_monodomain_default\n")
    #file.write("library_file=/home/berg/MonoAlg3D_C/shared_libs/libdefault_update_monodomain.so\n")
    file.write("\n")

    # Write [save_result] section
    file.write("[save_result]\n")
    file.write("print_rate = 100001\n")
    file.write("output_dir = %s\n" % (output_dir))
    file.write("main_function = save_as_vtu_tissue_coupled_vtp_purkinje\n")
    file.write("save_pvd = true\n")
    file.write("file_prefix = V_Tissue\n")
    file.write("file_prefix_purkinje = V_Purkinje\n")
    file.write("binary = false\n")
    file.write("compress = false\n")
    #file.write("library_file=/home/berg/MonoAlg3D_C/shared_libs/libdefault_save_mesh.so\n")
    file.write("\n")

    # Write [assembly_matrix] section
    file.write("[assembly_matrix]\n")
    file.write("init_function = set_initial_conditions_coupled_fvm\n")
    file.write("sigma_x = 0.0000176\n")
    file.write("sigma_y = 0.0000134\n")
    file.write("sigma_z = 0.0000176\n")
    file.write("sigma_purkinje = 0.0004\n")
    file.write("main_function = purkinje_coupled_endocardium_assembly_matrix\n")
    file.write("library_file = shared_libs/libpurkinje_coupled_matrix_assembly.so\n")
    #file.write("library_file=/home/berg/MonoAlg3D_C/shared_libs/libpurkinje_coupled_matrix_assembly.so\n")
    file.write("\n")

    # Write [linear_system_solver] section
    file.write("[linear_system_solver]\n")
    file.write("tolerance = 1e-16\n")
    file.write("use_preconditioner = no\n")
    file.write("max_iterations = 200\n")
    file.write("main_function = conjugate_gradient\n")
    file.write("library_file = shared_libs/libdefault_linear_system_solver.so\n")
    #file.write("library_file = /home/berg/MonoAlg3D_C/shared_libs/libdefault_linear_system_solver.so\n")
    file.write("\n")

    # Write [alg] section
    file.write("[alg]\n")
    file.write("refinement_bound = 0.11\n")
    file.write("derefinement_bound = 0.10\n")
    file.write("refine_each = 1\n")
    file.write("derefine_each = 1\n")
    file.write("\n")

    # Write [domain] section
    file.write("[domain]\n")
    file.write("name = Plain Mesh\n")
    file.write("num_layers = 1\n")
    file.write("start_dx = 200.0\n")
    file.write("start_dy = 200.0\n")
    file.write("start_dz = 200.0\n")
    file.write("side_length=20000\n")
    file.write("main_function = initialize_grid_with_square_mesh\n")
    #file.write("library_file = /home/berg/MonoAlg3D_C/shared_libs/libdefault_domains.so.so\n")    
    file.write("\n")

    # Write [ode_solver] section
    file.write("[ode_solver]\n")
    file.write("dt_ode = 0.02\n")
    file.write("use_gpu = no\n")
    file.write("gpu_id = 0\n")
    file.write("library_file = %s\n" % (ode_library_file))
    #file.write("library_file = /home/berg/MonoAlg3D_C/shared_libs/libten_tusscher_2006.so\n")
    file.write("\n")

    # Write [purkinje] section
    file.write("[purkinje]\n")
    file.write("name = Simple Purkinje\n")
    file.write("start_discretization = 100.0\n")
    file.write("start_dx = 100.0\n")
    file.write("rpmj = 1.0e+03\n")
    file.write("pmj_scale = 0.15\n")
    file.write("retro_propagation = false\n")
    file.write("network_file = networks/lucas-fractal-2.vtk\n")
    file.write("main_function = initialize_purkinje_with_custom_mesh\n")
    file.write("library_file = shared_libs/libdefault_purkinje.so\n")
    #file.write("library_file = /home/berg/MonoAlg3D_C/shared_libs/libdefault_purkinje.so\n")
    file.write("\n")

    # Write [purkinje_ode_solver] section
    file.write("[purkinje_ode_solver]\n")
    file.write("dt_ode = 0.01\n")
    file.write("use_gpu = no\n")
    file.write("gpu_id = 0\n")
    file.write("library_file = shared_libs/libli_rudy_2011.so\n")
    #file.write("library_file = /home/berg/MonoAlg3D_C/shared_libs/libstewart_aslanidi_noble_2009.so\n")
    file.write("\n")

    # Write [stimulus] section
    file.write("[stim_purkinje_s1]\n")
    file.write("start = 0.0\n")
    file.write("duration = 2.0\n")
    file.write("period = 1000.0\n")
    file.write("current = -80.0\n")
    file.write("id_limit = 10\n")
    file.write("main_function = stim_if_id_less_than\n")
    #file.write("library_file = /home/berg/MonoAlg3D_C/shared_libs/libdefault_stimuli.so")
    file.write("\n")

    file.close()


def main():
	
    num_scenarios = 2
    num_individuals = 20

    write_config_file(1,1)
    for scenario in range(2,num_scenarios+2):
        for individual in range(1,num_individuals+1):
            print("%u %u" % (scenario,individual))
            write_config_file(scenario,individual)

if __name__ == "__main__":
	main()
