//
// Created by sachetto on 12/10/17.
//

#ifndef MONOALG3D_INI_FILE_HEADERS_H
#define MONOALG3D_INI_FILE_HEADERS_H

#define ODE_SECTION "ode_solver"
#define MAIN_SECTION "main"
#define ALG_SECTION "alg"
#define STIM_SECTION "stim"

#define MATCH_SECTION_VALUE_AND_NAME(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
#define MATCH_NAME(v) strcmp(name, v) == 0
#define SECTION_STARTS_WITH(s) strncmp(section, s, strlen(s)) == 0


#endif //MONOALG3D_INI_FILE_HEADERS_H
