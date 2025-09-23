#include "../common_types/common_types.h"
#include "../config/config_parser.h"

struct changed_parameters {
    char *section;
    char *name;
    union {
        string_array values;
        char *value;
    };
};

struct simulation {
    uint32_t run_number;
    struct changed_parameters *parameters;
};

#ifdef __cplusplus
extern "C" {
#endif

void configure_new_parameters(struct changed_parameters *changed, struct user_options *options);
void free_current_simulation_resources(struct monodomain_solver *monodomain_solver, struct ode_solver *ode_solver, struct grid *the_grid);
struct changed_parameters parse_range_or_list_values(char *directive_rhs, char *directive_lhs);
struct simulation *generate_all_simulations(struct string_hash_entry *modify_directives, int num_sims);
void print_simulations(struct simulation *all_simulations);

#ifdef __cplusplus
}
#endif
