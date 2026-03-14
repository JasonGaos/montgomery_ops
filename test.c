/* Copyright (c) The mldsa-native project authors
 * SPDX-License-Identifier: Apache-2.0 OR ISC OR MIT
 */

#include "src/include/arith_native_aarch64.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NUM_RANDOM_TESTS 32

static uint32_t g_rng_state = 1u;

static uint32_t next_u32(void)
{
  uint32_t x = g_rng_state;
  x ^= x << 13;
  x ^= x >> 17;
  x ^= x << 5;
  g_rng_state = x;
  return x;
}

static void generate_i32_array_zeros(int32_t *data, size_t len)
{
  memset(data, 0, len * sizeof(int32_t));
}

static void generate_i32_array_single(int32_t *data, size_t len, size_t pos,
                                      int32_t value)
{
  memset(data, 0, len * sizeof(int32_t));
  data[pos] = value;
}

static void generate_i32_array_ranged(int32_t *data, size_t len, int min_incl,
                                      int max_excl)
{
  size_t i;
  uint32_t range = (uint32_t)(max_excl - min_incl);

  for (i = 0; i < len; i++)
  {
    data[i] = (int32_t)((unsigned)min_incl + (next_u32() % range));
  }
}

static int compare_i32_arrays(const int32_t *got, const int32_t *expected,
                              size_t len, const char *label)
{
  size_t i;

  for (i = 0; i < len; i++)
  {
    if (got[i] != expected[i])
    {
      fprintf(stderr,
              "FAIL: %s mismatch at index %zu: got=%d expected=%d\n",
              label, i, got[i], expected[i]);
      return 1;
    }
  }

  printf("PASS: %s\n", label);
  return 0;
}

static int run_poly_pointwise_case(const mld_poly *a, const mld_poly *b,
                                   const char *name)
{
  mld_poly ref;
  mld_poly neon;
  mld_poly sve2;
  char label[128];
  int errors = 0;

  generate_i32_array_zeros(ref.coeffs, MLDSA_N);
  generate_i32_array_zeros(neon.coeffs, MLDSA_N);
  generate_i32_array_zeros(sve2.coeffs, MLDSA_N);

  mld_poly_pointwise_montgomery_c(&ref, a, b);
  mld_poly_pointwise_montgomery_asm(neon.coeffs, a->coeffs, b->coeffs);
  mld_poly_pointwise_montgomery_asm_sve2(sve2.coeffs, a->coeffs, b->coeffs);

  snprintf(label, sizeof(label), "%s: NEON vs C", name);
  errors += compare_i32_arrays(neon.coeffs, ref.coeffs, MLDSA_N, label);

  snprintf(label, sizeof(label), "%s: SVE2 vs C", name);
  errors += compare_i32_arrays(sve2.coeffs, ref.coeffs, MLDSA_N, label);

  return errors;
}

static int test_poly_pointwise_montgomery(void)
{
  mld_poly a;
  mld_poly b;
  int errors = 0;
  int i;
  int pos;

  generate_i32_array_zeros(a.coeffs, MLDSA_N);
  generate_i32_array_zeros(b.coeffs, MLDSA_N);
  errors += run_poly_pointwise_case(&a, &b, "pointwise_montgomery_zeros");

  for (pos = 0; pos < MLDSA_N; pos += MLDSA_N / 8)
  {
    generate_i32_array_single(a.coeffs, MLDSA_N, (size_t)pos, 1);
    generate_i32_array_single(b.coeffs, MLDSA_N, (size_t)pos, 1);
    errors += run_poly_pointwise_case(&a, &b, "pointwise_montgomery_single");
  }

  for (i = 0; i < NUM_RANDOM_TESTS; i++)
  {
    generate_i32_array_ranged(a.coeffs, MLDSA_N, -(MLD_NTT_BOUND - 1),
                              MLD_NTT_BOUND);
    generate_i32_array_ranged(b.coeffs, MLDSA_N, -(MLD_NTT_BOUND - 1),
                              MLD_NTT_BOUND);
    errors += run_poly_pointwise_case(&a, &b, "pointwise_montgomery_random");
  }

  return errors;
}

#define DEFINE_ACC_CASE(L, vec_type, c_fn, neon_fn, sve2_fn)                     \
static int run_polyvecl_pointwise_acc_montgomery_l##L##_case(                    \
    const vec_type *u, const vec_type *v, const char *name)                      \
{                                                                                 \
  mld_poly ref;                                                                   \
  mld_poly neon;                                                                  \
  mld_poly sve2;                                                                  \
  char label[128];                                                                \
  int errors = 0;                                                                 \
                                                                                  \
  generate_i32_array_zeros(ref.coeffs, MLDSA_N);                                  \
  generate_i32_array_zeros(neon.coeffs, MLDSA_N);                                 \
  generate_i32_array_zeros(sve2.coeffs, MLDSA_N);                                 \
                                                                                  \
  c_fn(&ref, u, v);                                                               \
  neon_fn(neon.coeffs, (const int32_t *)u->vec, (const int32_t *)v->vec);        \
  sve2_fn(sve2.coeffs, (const int32_t *)u->vec, (const int32_t *)v->vec);        \
                                                                                  \
  snprintf(label, sizeof(label), "%s: NEON vs C", name);                         \
  errors += compare_i32_arrays(neon.coeffs, ref.coeffs, MLDSA_N, label);         \
                                                                                  \
  snprintf(label, sizeof(label), "%s: SVE2 vs C", name);                         \
  errors += compare_i32_arrays(sve2.coeffs, ref.coeffs, MLDSA_N, label);         \
                                                                                  \
  return errors;                                                                  \
}                                                                                 \
                                                                                  \
static int test_polyvecl_pointwise_acc_montgomery_l##L(void)                     \
{                                                                                 \
  vec_type u;                                                                     \
  vec_type v;                                                                     \
  int errors = 0;                                                                 \
  int i;                                                                          \
                                                                                  \
  generate_i32_array_zeros((int32_t *)&u, (size_t)(L) * MLDSA_N);                \
  generate_i32_array_zeros((int32_t *)&v, (size_t)(L) * MLDSA_N);                \
  errors += run_polyvecl_pointwise_acc_montgomery_l##L##_case(                   \
      &u, &v, "polyvecl_acc_l" #L "_zeros");                                     \
                                                                                  \
  for (i = 0; i < NUM_RANDOM_TESTS; i++)                                          \
  {                                                                               \
    generate_i32_array_ranged((int32_t *)&u, (size_t)(L) * MLDSA_N, 0, MLDSA_Q); \
    generate_i32_array_ranged((int32_t *)&v, (size_t)(L) * MLDSA_N,              \
                              -MLD_NTT_BOUND + 1, MLD_NTT_BOUND);                \
    errors += run_polyvecl_pointwise_acc_montgomery_l##L##_case(                 \
        &u, &v, "polyvecl_acc_l" #L "_random");                                  \
  }                                                                               \
                                                                                  \
  return errors;                                                                  \
}

DEFINE_ACC_CASE(4, montgomery_polyvecl_l4, mld_polyvecl_pointwise_acc_montgomery_l4_c,
                mld_polyvecl_pointwise_acc_montgomery_l4_asm,
                mld_polyvecl_pointwise_acc_montgomery_l4_asm_sve2)
DEFINE_ACC_CASE(5, montgomery_polyvecl_l5, mld_polyvecl_pointwise_acc_montgomery_l5_c,
                mld_polyvecl_pointwise_acc_montgomery_l5_asm,
                mld_polyvecl_pointwise_acc_montgomery_l5_asm_sve2)
DEFINE_ACC_CASE(7, montgomery_polyvecl_l7, mld_polyvecl_pointwise_acc_montgomery_l7_c,
                mld_polyvecl_pointwise_acc_montgomery_l7_asm,
                mld_polyvecl_pointwise_acc_montgomery_l7_asm_sve2)

int main(void)
{
  int errors = 0;

  printf("Running montgomery_ops tests.\n");
  printf("Default SVE2 placeholders are expected to fail until implemented.\n");

  errors += test_poly_pointwise_montgomery();
  errors += test_polyvecl_pointwise_acc_montgomery_l4();
  errors += test_polyvecl_pointwise_acc_montgomery_l5();
  errors += test_polyvecl_pointwise_acc_montgomery_l7();

  if (errors != 0)
  {
    fprintf(stderr, "montgomery_ops completed with %d mismatches.\n", errors);
    return 1;
  }

  printf("montgomery_ops completed without mismatches.\n");
  return 0;
}
