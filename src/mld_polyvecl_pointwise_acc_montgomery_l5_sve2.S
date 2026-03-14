/* Copyright (c) The mldsa-native project authors
 * SPDX-License-Identifier: Apache-2.0 OR ISC OR MIT
 */

#include "../common.h"
#if defined(MLD_ARITH_BACKEND_AARCH64) && !defined(MLD_CONFIG_MULTILEVEL_NO_SHARED) && \
           (defined(MLD_CONFIG_MULTILEVEL_WITH_SHARED) || MLDSA_L == 5)

#if defined(__ELF__)
.section .note.GNU-stack,"",@progbits
#endif

.text
.global MLD_ASM_NAMESPACE(polyvecl_pointwise_acc_montgomery_l5_asm_sve2)
.balign 4
MLD_ASM_FN_SYMBOL(polyvecl_pointwise_acc_montgomery_l5_asm_sve2)
    // TODO: replace this placeholder with a real SVE2 implementation.
    movz w3, #0xe001
    movk w3, #0x7f, lsl #16
    str  w3, [x0]
    ret

MLD_ASM_FN_SIZE(polyvecl_pointwise_acc_montgomery_l5_asm_sve2)

#endif /* MLD_ARITH_BACKEND_AARCH64 && !MLD_CONFIG_MULTILEVEL_NO_SHARED && \
          (MLD_CONFIG_MULTILEVEL_WITH_SHARED || MLDSA_L == 5) */
