////
//// Created by sachetto on 06/10/17.
////
#include <criterion/criterion.h>
#include <criterion/internal/assert.h>
#include <signal.h>

#include "../alg/grid/grid.h"
#include "../config/linear_system_solver_config.h"
#include "../config/domain_config.h"
#include "../config/save_mesh_config.h"
#include "../config/config_parser.h"
#include "../utils/file_utils.h"
#include "../3dparty/ini_parser/ini.h"
#include "../3dparty/sds/sds.h"
#include "../logger/logger.h"

#include "../3dparty/stb_ds.h"
#include "../config_helpers/config_helpers.h"

#include "../3dparty/fast_double_parser.h"

Test (utils, arr_int) {

    int *v = NULL;

    arrsetcap(v, 1);

    cr_assert_eq(arrlen(v), 0);
    cr_assert_geq(arrcap(v), 1);

    arrput(v, 0);
    arrput(v, 1);
    arrput(v, 2);

    cr_assert_eq(arrlen(v), 3);

    cr_assert_eq(v[0], 0);
    cr_assert_eq(v[1], 1);
    cr_assert_eq(v[2], 2);
}

Test (utils, arr_float) {

    float *v = NULL;

    arrsetcap(v, 1);

    cr_assert_eq(arrlen(v), 0);
    cr_assert_geq(arrcap(v), 1);

    arrput(v, 0);
    arrput(v, 0.5);
    arrput(v, 2.5);

    cr_assert_eq(arrlen(v), 3);

    cr_assert_float_eq(v[0], 0.0, 1e-10);
    cr_assert_float_eq(v[1], 0.5, 1e-10);
    cr_assert_float_eq(v[2], 2.5, 1e-10);
}

Test (utils, arr_double) {

    real_cpu *v = NULL;

    arrsetcap(v, 1);

    cr_assert_eq(arrlen(v), 0);
    cr_assert_geq(arrcap(v), 1);

    arrput(v, 0);
    arrput(v, 0.5);
    arrput(v, 2.5);

    cr_assert_eq(arrlen(v), 3);

    cr_assert_float_eq(v[0], 0.0, 1e-10);
    cr_assert_float_eq(v[1], 0.5, 1e-10);
    cr_assert_float_eq(v[2], 2.5, 1e-10);

}

Test (utils, arr_element) {

    struct element *v = NULL;

    arrsetcap(v, 1);

    cr_assert_eq(arrlen(v), 0);
    cr_assert_geq(arrcap(v), 1);

    struct cell_node *c = new_cell_node();
#ifdef ENABLE_DDM
    struct element a = {'a', 0, 0, 1, c};
#else
    struct element a = {0, 0, 1, c};
#endif
    arrput(v, a);

    a.column = 2;
    a.value = -2.2;
    arrput(v, a);

    a.column = 3;
    a.value = 3.5;
    arrput(v, a);

    cr_assert_eq(arrlen(v), 3);

    cr_assert_eq(v[0].column, 1);
    cr_assert_float_eq(v[0].value, 0.0, 1e-10);
    cr_assert_eq(v[0].cell, c);

    cr_assert_eq(v[1].column, 2);
    cr_assert_float_eq(v[1].value, -2.2, 1e-10);
    cr_assert_eq(v[1].cell, c);

    cr_assert_eq(v[2].column, 3);
    cr_assert_float_eq(v[2].value, 3.5, 1e-10);
    cr_assert_eq(v[2].cell, c);

    struct element b = arrpop(v);
    cr_assert_eq(arrlen(v), 2);

    cr_assert_eq(b.column, 3);

    free(c);
}

static void test_file(char *file) {
{
        size_t size;
        char *values = read_entire_file_with_mmap(file, &size);
        char *fast_floats = values;
        char *strtod_floats = values;
        while(*fast_floats) {
            double n_fast;
            const char *end = parse_number(fast_floats, &n_fast);

            if(!end) break;

            fast_floats += (end - fast_floats + 1);

            double n_strtod = 0;
            char *end2;
            n_strtod = strtod(strtod_floats, &end2);

            strtod_floats += (end2 - strtod_floats + 1);

            cr_assert_float_eq(n_fast, n_strtod, 0.00001);
        }

        munmap(values, size);
    }
}

Test (float_parser, fast_parser) {
    test_file("./tests_bin/float_pos.txt");
    test_file("./tests_bin/float_neg.txt");
    test_file("./tests_bin/float_mixed.txt");
}
