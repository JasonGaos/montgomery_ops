/* Copyright (c) The mldsa-native project authors
 * SPDX-License-Identifier: Apache-2.0 OR ISC OR MIT
 */

#include "include/arith_native_aarch64.h"

#include "montgomery_reduce.h"

#define MLDSA_L 7
#define mld_polyvecl montgomery_polyvecl_l7
#define mld_polyvecl_pointwise_acc_montgomery_c \
  mld_polyvecl_pointwise_acc_montgomery_l7_c

#include "polyvecl_pointwise_acc_montgomery_ref.inc"

#undef mld_polyvecl_pointwise_acc_montgomery_c
#undef mld_polyvecl
#undef MLDSA_L
