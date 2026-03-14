/* Copyright (c) The mldsa-native project authors
 * SPDX-License-Identifier: Apache-2.0 OR ISC OR MIT
 */

#include "include/arith_native_aarch64.h"

#include "montgomery_reduce.h"

MLD_STATIC_TESTABLE void mld_poly_pointwise_montgomery_c(mld_poly *c,
                                                         const mld_poly *a,
                                                         const mld_poly *b)
__contract__(
  requires(memory_no_alias(a, sizeof(mld_poly)))
  requires(memory_no_alias(b, sizeof(mld_poly)))
  requires(memory_no_alias(c, sizeof(mld_poly)))
  requires(array_abs_bound(a->coeffs, 0, MLDSA_N, MLD_NTT_BOUND))
  requires(array_abs_bound(b->coeffs, 0, MLDSA_N, MLD_NTT_BOUND))
  assigns(memory_slice(c, sizeof(mld_poly)))
  ensures(array_abs_bound(c->coeffs, 0, MLDSA_N, MLDSA_Q))
)
{
  unsigned int i;
  mld_assert_abs_bound(a->coeffs, MLDSA_N, MLD_NTT_BOUND);
  mld_assert_abs_bound(b->coeffs, MLDSA_N, MLD_NTT_BOUND);

  for (i = 0; i < MLDSA_N; ++i)
  __loop__(
    invariant(i <= MLDSA_N)
    invariant(array_abs_bound(c->coeffs, 0, i, MLDSA_Q))
  )
  {
    c->coeffs[i] = mld_montgomery_reduce((int64_t)a->coeffs[i] * b->coeffs[i]);
  }
  mld_assert_abs_bound(c->coeffs, MLDSA_N, MLDSA_Q);
}
