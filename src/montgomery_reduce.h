/* Copyright (c) The mldsa-native project authors
 * SPDX-License-Identifier: Apache-2.0 OR ISC OR MIT
 */

#ifndef MONTGOMERY_OPS_REDUCE_H
#define MONTGOMERY_OPS_REDUCE_H

#include "montgomery_ref_compat.h"

/* check-magic: -4186625 == pow(2,32,MLDSA_Q) */
#define MLD_MONT -4186625

MLD_MUST_CHECK_RETURN_VALUE
static MLD_INLINE int32_t mld_montgomery_reduce(int64_t a)
__contract__(
  requires(a > -(((int64_t)1 << 31) * MLDSA_Q) &&
           a < (((int64_t)1 << 31) * MLDSA_Q))
  ensures(return_value > -MLDSA_Q && return_value < MLDSA_Q))
{
  const uint64_t QINV = 58728449;
  const uint32_t a_reduced = mld_cast_int64_to_uint32(a);
  const uint32_t a_inverted = (a_reduced * QINV) & UINT32_MAX;
  const int32_t t = mld_cast_uint32_to_int32(a_inverted);

  int64_t r;

  mld_assert(a < +(INT64_MAX - (((int64_t)1 << 31) * MLDSA_Q)) &&
             a > -(INT64_MAX - (((int64_t)1 << 31) * MLDSA_Q)));

  r = a - (int64_t)t * MLDSA_Q;
  r = r >> 32;

  return (int32_t)r;
}

#endif /* MONTGOMERY_OPS_REDUCE_H */
