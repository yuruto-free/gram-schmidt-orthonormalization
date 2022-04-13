#include <stddef.h>
#include <stdint.h>
#include <malloc.h>
#include <math.h>
#include "gram_schmidt.h"

#define RETURN_OK_GS (0)
#define RETURN_NG_GS (1)
#define MEPS_GS (1e-10)

/*
 * \brief ノルムの計算
 * \param[in]  dim  次元数
 * \param[in]  vec  ベクトル
 * \param[out] norm ノルムの値
 * \retval RETURN_OK_GS 正常終了
 * \retval RETURN_NG_GS 異常終了
*/
static int32_t euclid_norm(int32_t dim, double *vec, double *norm);
/*
 * \brief 内積の計算
 * \param[in]  dim        次元数
 * \param[in]  *left_vec  ベクトル
 * \param[in]  *right_vec ベクトル
 * \param[out] *dot       内積値
 * \retval RETURN_OK_GS 正常終了
 * \retval RETURN_NG_GS 異常終了
*/
static int32_t inner_product(int32_t dim, double *left_vec, double *right_vec, double *dot);
/*
 * \brief ベクトルの差の計算
 * \param[in]    dim     次元数
 * \param[in]    coef    係数
 * \param[inout] *target 引かれる数となるベクトル
 * \param[in]    *base   引く数となるベクトル
 * \retval RETURN_OK_GS 正常終了
 * \retval RETURN_NG_GS 異常終了
*/
static int32_t sub_vector(int32_t dim, double coef, double *target, double *base);
/*
 * \brief 直交するベクトルの計算
 * \param[in]    dim         次元数
 * \param[in]    current_idx 計算対象のインデックス
 * \param[inout] *dots       内積値（作業用）
 * \param[inout] vecs        ベクトル（dim * vec_num、正規直交化するベクトル：行ベクトル形式）
 * \retval RETURN_OK_GS 正常終了
 * \retval RETURN_NG_GS 異常終了
*/
static int32_t orthogonal_vector(int32_t dim, int32_t current_idx, double *dots, double *vecs);
/*
 * \brief スケーリング
 * \param[in]    dim  次元数
 * \param[in]    coef スケーリング値
 * \param[inout] vec  ベクトル
 * \retval RETURN_OK_GS 正常終了
 * \retval RETURN_NG_GS 異常終了
*/
static int32_t scaling(int32_t dim, double coef, double *vec);

int32_t GS_orthonormalization(int32_t dim, int32_t vec_num, double *vecs) {
    int32_t ret = (int32_t)GS_NG;
    int32_t func_val;
    int32_t k;
    double norm;
    double *target;
    double *work;

    if (NULL != vecs) {
        work = (double *)malloc(sizeof(double) * (vec_num - 1));
        if (NULL == work) {
            goto EXIT_GS_ORTH;
        }

        target = &vecs[0];
        /* e_{0}のノルムを計算 */
        func_val = euclid_norm(dim, target, &norm);
        if (((int32_t)RETURN_OK_GS != func_val) || (fabs(norm) < (double)MEPS_GS)) {
            goto EXIT_GS_ORTH;
        }
        /* e_{0}を計算 */
        func_val = scaling(dim, 1.0/norm, target);
        if ((int32_t)RETURN_OK_GS != func_val) {
            goto EXIT_GS_ORTH;
        }
        /* e_{1}からe_{n-1}を計算 */
        for (k = 1; k < vec_num; k++) {
            target = &vecs[k * dim];
            /* 直交ベクトルf_{k}を計算 */
            func_val = orthogonal_vector(dim, k, work, vecs);
            if ((int32_t)RETURN_OK_GS != func_val) {
                goto EXIT_GS_ORTH;
            }
            /* f_{k}のノルムを計算 */
            func_val = euclid_norm(dim, target, &norm);
            if (((int32_t)RETURN_OK_GS != func_val) || (fabs(norm) < (double)MEPS_GS)) {
                goto EXIT_GS_ORTH;
            }
            /* e_{k}を計算 */
            func_val = scaling(dim, 1.0/norm, target);
            if ((int32_t)RETURN_OK_GS != func_val) {
                goto EXIT_GS_ORTH;
            }
        }
        free(work);
        ret = (int32_t)GS_OK;
    }
EXIT_GS_ORTH:

    return ret;
}

static int32_t euclid_norm(int32_t dim, double *vec, double *norm) {
    int32_t ret = (int32_t)RETURN_NG_GS;
    int32_t i;
    double sum;

    if ((NULL != vec) && (NULL != norm)) {
        sum = 0.0;

        for (i = 0; i < dim; i++) {
            sum += vec[i] * vec[i];
        }
        (*norm) = sqrt(sum);
        ret = (int32_t)RETURN_OK_GS;
    }

    return ret;
}

static int32_t inner_product(int32_t dim, double *left_vec, double *right_vec, double *dot) {
    int32_t ret = (int32_t)RETURN_NG_GS;
    int32_t i;
    double sum;

    if ((NULL != left_vec) && (NULL != right_vec) && (NULL != dot)) {
        sum = 0.0;

        for (i = 0; i < dim; i++) {
            sum += left_vec[i] * right_vec[i];
        }
        (*dot) = sum;
        ret = (int32_t)RETURN_OK_GS;
    }

    return ret;
}

static int32_t sub_vector(int32_t dim, double coef, double *target, double *base) {
    int32_t ret = (int32_t)RETURN_NG_GS;
    int32_t i;

    if ((NULL != target) && (NULL != base)) {
        for (i = 0; i < dim; i++) {
            target[i] -= coef * base[i];
        }
        ret = (int32_t)RETURN_OK_GS;
    }

    return ret;
}

static int32_t orthogonal_vector(int32_t dim, int32_t current_idx, double *dots, double *vecs) {
    int32_t ret = (int32_t)RETURN_NG_GS;
    int32_t i;
    int32_t func_val;
    double *target;

    if ((NULL != dots) && (NULL != vecs)) {
        target = &vecs[current_idx * dim];

        /* 内積を計算 */
        for (i = 0; i < current_idx; i++) {
            func_val = inner_product(dim, target, &vecs[i * dim], &dots[i]);

            if ((int32_t)RETURN_OK_GS != func_val) {
                goto EXIT_ORTH_VEC;
            }
        }
        /* ベクトルの差を計算 */
        for (i = 0; i < current_idx; i++) {
            func_val = sub_vector(dim, dots[i], target, &vecs[i * dim]);

            if ((int32_t)RETURN_OK_GS != func_val) {
                goto EXIT_ORTH_VEC;
            }
        }
        ret = (int32_t)RETURN_OK_GS;
    }
EXIT_ORTH_VEC:

    return ret;
}

static int32_t scaling(int32_t dim, double coef, double *vec) {
    int32_t ret = (int32_t)RETURN_NG_GS;
    int32_t i;

    if (NULL != vec) {
        for (i = 0; i < dim; i++) {
            vec[i] *= coef;
        }
        ret = (int32_t)RETURN_OK_GS;
    }

    return ret;
}
