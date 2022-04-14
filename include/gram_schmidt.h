#ifndef GRAM_SCHMIDT_H__
#define GRAM_SCHMIDT_H__

#include <stdint.h>

#define GS_OK (0)
#define GS_NG (1)

/*
 * \brief グラム・シュミットの正規直交化法による正規直交基底への変換
 * \param[in]    dim  次元数
 * \param[in]    num  正規直交化するベクトルの数
 * \param[inout] vecs ベクトル
 *                    サイズ：dim * vec_num
 *                    正規直交化するベクトル：行ベクトル形式
 * \retval       GS_OK   正常終了
 * \retval       GS_NG   異常終了
 */
int32_t GS_orthonormalization(int32_t dim, int32_t num, double *vecs);

#endif
