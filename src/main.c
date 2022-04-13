#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include "gram_schmidt.h"
#define MAX_TEST_DATA (3)
#define MEPS_MAIN (1e-10)

struct data_t {
    int32_t dim;
    int32_t vec_num;
    double *vecs;
    const double *ans;
};

double vec1[9] = {
    1.0, 1.0, 0.0,
    1.0, 0.0, 2.0,
    2.0, 1.0, 3.0,
};
double vec2[12] = {
    1.0, 1.0, -2.0, 2.0,
    0.0, 1.0, -1.0, 0.0,
    3.0, 5.0, -2.0, 1.0,
};
double vec3[20] = {
    1.0,  2.0,  1.0,  3.0, -1.0,
    0.0, -2.0, -3.0,  3.0, -2.0,
    2.0, -3.0,  0.0, -4.0,  3.0,
    3.0,  1.0,  2.0,  3.0, -4.0,
};
const double ans1[9] = {
    sqrt(2)/2,  sqrt(2)/2,         0.0,
    sqrt(2)/6, -sqrt(2)/6, 2*sqrt(2)/3,
    -2.0/3.0,    2.0/3.0,     1.0/3.0,
};
const double ans2[12] = {
    sqrt(10)/10,     sqrt(10)/10,     -sqrt(10)/5,      sqrt(10)/5,
    -3*sqrt(110)/110, 7*sqrt(110)/110, -2*sqrt(110)/55, -3*sqrt(110)/55,
    26*sqrt(165)/495,  4*sqrt(165)/99,  4*sqrt(165)/99,  -sqrt(165)/165,
};
const double ans3[20] = {
    0.25, 0.5, 0.25, 0.75, -0.25,
    -0.05, -0.5, -0.65, 0.45, -0.35,
    89*sqrt(5334)/7620, -9*sqrt(5334)/1778, 7*sqrt(5334)/7620, 45*sqrt(5334)/53340, 87*sqrt(5334)/17780,
    334*sqrt(347494098)/24821007, -626*sqrt(347494098)/57915683, 913*sqrt(347494098)/49642014,-6257*sqrt(347494098)/347494098,-2536*sqrt(347494098)/57915683,
};

struct data_t TEST_DATA[MAX_TEST_DATA] = {
    {3, 3, vec1, ans1},
    {4, 3, vec2, ans2},
    {5, 4, vec3, ans3},
};

static double calculate_error(const struct data_t *target) {
    int32_t i, j;
    int32_t dim, vec_num;
    double *estimated_vec;
    const double *true_vec;
    double sum = 0.0;
    dim = target->dim;
    vec_num = target->vec_num;

    for (i = 0; i < vec_num; i++) {
        estimated_vec = &(target->vecs[i * dim]);
        true_vec = &(target->ans[i * dim]);

        for (j = 0; j < dim; j++) {
            sum += fabs(estimated_vec[j] - true_vec[j]);
        }
    }

    return sum;
}

static void print_vec(const struct data_t *target) {
    int32_t i, j;
    int32_t dim, vec_num;
    double *estimated_vec;
    const double *true_vec;
    dim = target->dim;
    vec_num = target->vec_num;

    for (i = 0; i < vec_num; i++) {
        estimated_vec = &(target->vecs[i * dim]);
        true_vec = &(target->ans[i * dim]);
        printf("[%03d]\n", i);

        for (j = 0; j < dim; j++) {
            printf("    %+.5f (%.13e)\n", estimated_vec[j], fabs(estimated_vec[j] - true_vec[j]));
        }
        printf("\n");
    }
}

int main(int argc, char **argv) {
    struct data_t target;
    double err;
    int32_t ret;
    int i;

    for (i = 0; i < (int)MAX_TEST_DATA; i++) {
        memcpy(&target, &TEST_DATA[i], sizeof(struct data_t));
        printf("Dim: %d, Num: %d\n", target.dim, target.vec_num);
        ret = GS_orthonormalization(target.dim, target.vec_num, target.vecs);
        err = calculate_error((const struct data_t *)&target);
        print_vec((const struct data_t *)&target);
        printf("error: %.13e\n", err);
        printf("\n");

        if ((int32_t)GS_OK != ret) {
            goto EXIT_MAIN;
        }
    }
EXIT_MAIN:

    return 0;
}
