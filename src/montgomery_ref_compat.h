/* Copyright (c) The mldsa-native project authors
 * Copyright (c) The mlkem-native project authors
 * SPDX-License-Identifier: Apache-2.0 OR ISC OR MIT
 */

#ifndef MONTGOMERY_OPS_REF_COMPAT_H
#define MONTGOMERY_OPS_REF_COMPAT_H

#include "../common.h"

#include <limits.h>
#include <stdint.h>

#define __contract__(x)
#define __loop__(x)
#define cassert(x)

#define mld_assert(val) \
  do                    \
  {                     \
  } while (0)

#define mld_assert_bound(ptr, len, value_lb, value_ub) \
  do                                                   \
  {                                                    \
  } while (0)

#define mld_assert_abs_bound(ptr, len, value_abs_bd) \
  do                                                 \
  {                                                  \
  } while (0)

#define mld_assert_bound_2d(ptr, len0, len1, value_lb, value_ub) \
  do                                                             \
  {                                                              \
  } while (0)

#define mld_assert_abs_bound_2d(ptr, len0, len1, value_abs_bd) \
  do                                                           \
  {                                                            \
  } while (0)

MLD_MUST_CHECK_RETURN_VALUE
static MLD_ALWAYS_INLINE int32_t mld_cast_uint32_to_int32(uint32_t x)
{
  return (int32_t)x;
}

MLD_MUST_CHECK_RETURN_VALUE
static MLD_ALWAYS_INLINE uint32_t mld_cast_int64_to_uint32(int64_t x)
{
  return (uint32_t)(x & (int64_t)UINT32_MAX);
}

#endif /* MONTGOMERY_OPS_REF_COMPAT_H */
