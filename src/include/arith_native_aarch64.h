/* Copyright (c) The mldsa-native project authors
 * SPDX-License-Identifier: Apache-2.0 OR ISC OR MIT
 */

#ifndef MONTGOMERY_OPS_ARITH_NATIVE_AARCH64_H
#define MONTGOMERY_OPS_ARITH_NATIVE_AARCH64_H

#include "../../common.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  int32_t coeffs[MLDSA_N];
} MLD_ALIGN mld_poly;

typedef struct
{
  mld_poly vec[4];
} MLD_ALIGN montgomery_polyvecl_l4;

typedef struct
{
  mld_poly vec[5];
} MLD_ALIGN montgomery_polyvecl_l5;

typedef struct
{
  mld_poly vec[7];
} MLD_ALIGN montgomery_polyvecl_l7;

#define mld_poly_pointwise_montgomery_c \
  MLD_NAMESPACE(poly_pointwise_montgomery_c)
void mld_poly_pointwise_montgomery_c(mld_poly *c, const mld_poly *a,
                                     const mld_poly *b);

#define mld_polyvecl_pointwise_acc_montgomery_l4_c \
  MLD_NAMESPACE(polyvecl_pointwise_acc_montgomery_l4_c)
void mld_polyvecl_pointwise_acc_montgomery_l4_c(mld_poly *w,
                                                const montgomery_polyvecl_l4 *u,
                                                const montgomery_polyvecl_l4 *v);

#define mld_polyvecl_pointwise_acc_montgomery_l5_c \
  MLD_NAMESPACE(polyvecl_pointwise_acc_montgomery_l5_c)
void mld_polyvecl_pointwise_acc_montgomery_l5_c(mld_poly *w,
                                                const montgomery_polyvecl_l5 *u,
                                                const montgomery_polyvecl_l5 *v);

#define mld_polyvecl_pointwise_acc_montgomery_l7_c \
  MLD_NAMESPACE(polyvecl_pointwise_acc_montgomery_l7_c)
void mld_polyvecl_pointwise_acc_montgomery_l7_c(mld_poly *w,
                                                const montgomery_polyvecl_l7 *u,
                                                const montgomery_polyvecl_l7 *v);

#define mld_poly_pointwise_montgomery_asm \
  MLD_NAMESPACE(poly_pointwise_montgomery_asm)
void mld_poly_pointwise_montgomery_asm(int32_t *out, const int32_t *in0,
                                       const int32_t *in1);

#define mld_poly_pointwise_montgomery_asm_sve2 \
  MLD_NAMESPACE(poly_pointwise_montgomery_asm_sve2)
void mld_poly_pointwise_montgomery_asm_sve2(int32_t *out, const int32_t *in0,
                                            const int32_t *in1);

#define mld_polyvecl_pointwise_acc_montgomery_l4_asm \
  MLD_NAMESPACE(polyvecl_pointwise_acc_montgomery_l4_asm)
void mld_polyvecl_pointwise_acc_montgomery_l4_asm(int32_t *w,
                                                  const int32_t *u,
                                                  const int32_t *v);

#define mld_polyvecl_pointwise_acc_montgomery_l4_asm_sve2 \
  MLD_NAMESPACE(polyvecl_pointwise_acc_montgomery_l4_asm_sve2)
void mld_polyvecl_pointwise_acc_montgomery_l4_asm_sve2(int32_t *w,
                                                       const int32_t *u,
                                                       const int32_t *v);

#define mld_polyvecl_pointwise_acc_montgomery_l5_asm \
  MLD_NAMESPACE(polyvecl_pointwise_acc_montgomery_l5_asm)
void mld_polyvecl_pointwise_acc_montgomery_l5_asm(int32_t *w,
                                                  const int32_t *u,
                                                  const int32_t *v);

#define mld_polyvecl_pointwise_acc_montgomery_l5_asm_sve2 \
  MLD_NAMESPACE(polyvecl_pointwise_acc_montgomery_l5_asm_sve2)
void mld_polyvecl_pointwise_acc_montgomery_l5_asm_sve2(int32_t *w,
                                                       const int32_t *u,
                                                       const int32_t *v);

#define mld_polyvecl_pointwise_acc_montgomery_l7_asm \
  MLD_NAMESPACE(polyvecl_pointwise_acc_montgomery_l7_asm)
void mld_polyvecl_pointwise_acc_montgomery_l7_asm(int32_t *w,
                                                  const int32_t *u,
                                                  const int32_t *v);

#define mld_polyvecl_pointwise_acc_montgomery_l7_asm_sve2 \
  MLD_NAMESPACE(polyvecl_pointwise_acc_montgomery_l7_asm_sve2)
void mld_polyvecl_pointwise_acc_montgomery_l7_asm_sve2(int32_t *w,
                                                       const int32_t *u,
                                                       const int32_t *v);

#ifdef __cplusplus
}
#endif

#endif /* MONTGOMERY_OPS_ARITH_NATIVE_AARCH64_H */
