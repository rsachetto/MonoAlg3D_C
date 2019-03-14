//
// Created by sachetto on 12/10/17.
//

#ifndef MONOALG3D_INI_FILE_HEADERS_H
#define MONOALG3D_INI_FILE_HEADERS_H

#define ODE_SECTION "ode_solver"
#define MAIN_SECTION "main"
#define ALG_SECTION "alg"
#define STIM_SECTION "stim"
#define DOMAIN_SECTION "domain"
#define PURKINJE_SECTION "purkinje"
#define SAVE_RESULT_SECTION "save_result"
#define SAVE_STATE_SECTION "save_state" 
#define RESTORE_STATE_SECTION "restore_state" 
#define MATRIX_ASSEMBLY_SECTION "assembly_matrix"
#define UPDATE_MONODOMAIN_SECTION "update_monodomain"
#define LINEAR_SYSTEM_SOLVER_SECTION "linear_system_solver"
#define EXTRA_DATA_SECTION "extra_data"
#define BATCH_SECTION "batch"
#define MODIFICATION_SECTION "modify"

#define MATCH_SECTION_AND_NAME(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
#define MATCH_SECTION(s) strcmp(section, s) == 0
#define MATCH_NAME(v) strcmp(name, v) == 0
#define SECTION_STARTS_WITH(s) strncmp(section, s, strlen(s)) == 0
#define A_STARTS_WITH_B(a ,b) strncmp(a, b, strlen(b)) == 0


#endif //MONOALG3D_INI_FILE_HEADERS_H
