/* Copyright (c) The mldsa-native project authors
 * SPDX-License-Identifier: Apache-2.0 OR ISC OR MIT
 */

#ifndef MONTGOMERY_OPS_COMMON_H
#define MONTGOMERY_OPS_COMMON_H

#ifndef __ASSEMBLER__
#include <stddef.h>
#include <stdint.h>
#endif

#define MLD_ARITH_BACKEND_AARCH64
#define MLD_CONFIG_MULTILEVEL_WITH_SHARED

#define MLDSA_N 256
#define MLDSA_Q 8380417
#define MLDSA_Q_HALF ((MLDSA_Q + 1) / 2)
#define MLD_NTT_BOUND (9 * MLDSA_Q)

#if !defined(MLD_INLINE)
#if defined(_MSC_VER)
#define MLD_INLINE __inline
#elif defined(inline) || \
    (defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L)
#define MLD_INLINE inline
#elif defined(__GNUC__) || defined(__clang__)
#define MLD_INLINE __attribute__((unused))
#else
#define MLD_INLINE
#endif
#endif

#if !defined(MLD_ALWAYS_INLINE)
#if defined(_MSC_VER)
#define MLD_ALWAYS_INLINE __forceinline
#elif (defined(__GNUC__) || defined(__clang__)) && \
    (defined(inline) ||                            \
     (defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L))
#define MLD_ALWAYS_INLINE MLD_INLINE __attribute__((always_inline))
#else
#define MLD_ALWAYS_INLINE MLD_INLINE
#endif
#endif

#if defined(__GNUC__) || defined(__clang__)
#define MLD_MUST_CHECK_RETURN_VALUE __attribute__((warn_unused_result))
#else
#define MLD_MUST_CHECK_RETURN_VALUE
#endif

#ifndef MLD_STATIC_TESTABLE
#define MLD_STATIC_TESTABLE
#endif

#define MLD_DEFAULT_ALIGN 32
#if defined(__GNUC__)
#define MLD_ALIGN __attribute__((aligned(MLD_DEFAULT_ALIGN)))
#elif defined(_MSC_VER)
#define MLD_ALIGN __declspec(align(MLD_DEFAULT_ALIGN))
#else
#define MLD_ALIGN
#endif

#define MLD_CONCAT_(x1, x2) x1##x2
#define MLD_CONCAT(x1, x2) MLD_CONCAT_(x1, x2)
#define MLD_CONFIG_NAMESPACE_PREFIX mldsa
#define MLD_NAMESPACE_PREFIX MLD_CONCAT(MLD_CONFIG_NAMESPACE_PREFIX, _)
#define MLD_NAMESPACE(sym) MLD_CONCAT(MLD_NAMESPACE_PREFIX, sym)

#if defined(__APPLE__)
#define MLD_ASM_NAMESPACE(sym) MLD_CONCAT(_, MLD_NAMESPACE(sym))
#else
#define MLD_ASM_NAMESPACE(sym) MLD_NAMESPACE(sym)
#endif

#define MLD_ASM_FN_SYMBOL(sym) MLD_ASM_NAMESPACE(sym) :

#if defined(__ELF__)
#define MLD_ASM_FN_SIZE(sym) \
  .size MLD_ASM_NAMESPACE(sym), .- MLD_ASM_NAMESPACE(sym)
#else
#define MLD_ASM_FN_SIZE(sym)
#endif

#endif /* MONTGOMERY_OPS_COMMON_H */
