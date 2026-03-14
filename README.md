# montgomery_ops

`montgomery_ops` is a small standalone sandbox for the Montgomery-related
polynomial operations used by `mldsa-native`.

It exists to make SVE2 development easier without building the full ML-DSA
tree. The codebase keeps three implementations side by side:

- Exact copied C reference logic from `mldsa-native/src`
- Copied AArch64 NEON assembly from `dev/aarch64_opt/src`
- Placeholder AArch64 SVE2 assembly with `_sve2` suffix

The default build links all three paths into one test binary.

## Scope

This sandbox currently covers four functions:

- `poly_pointwise_montgomery`
- `polyvecl_pointwise_acc_montgomery_l4`
- `polyvecl_pointwise_acc_montgomery_l5`
- `polyvecl_pointwise_acc_montgomery_l7`

The public declarations live in `src/include/arith_native_aarch64.h`.

## What Each Function Does

- `mld_poly_pointwise_montgomery_{c,asm,asm_sve2}`
  Performs coefficient-wise multiplication of two NTT-domain polynomials and
  applies Montgomery reduction to each product.

- `mld_polyvecl_pointwise_acc_montgomery_l4_{c,asm,asm_sve2}`
  Computes the Montgomery-domain inner product of two length-4 polynomial
  vectors and accumulates the result into one output polynomial.

- `mld_polyvecl_pointwise_acc_montgomery_l5_{c,asm,asm_sve2}`
  Same operation for length-5 polynomial vectors.

- `mld_polyvecl_pointwise_acc_montgomery_l7_{c,asm,asm_sve2}`
  Same operation for length-7 polynomial vectors.

The explicit `l4`, `l5`, and `l7` names let all parameter variants coexist in
one standalone binary.

## Source Provenance

- C reference logic comes from upstream `mldsa/src/poly.c`,
  `mldsa/src/polyvec.c`, `mldsa/src/reduce.h`, and `mldsa/src/ct.h`
- NEON assembly comes from upstream `dev/aarch64_opt/src`
- Minimal compatibility shims are local to this sandbox

The goal is to preserve upstream logic, not to re-derive the algorithms.

## Repository Layout

- `common.h`
  Minimal standalone configuration and namespace macros

- `src/include/arith_native_aarch64.h`
  Public types and entrypoint declarations

- `src/*_ref.c`, `src/*_ref.inc`, `src/montgomery_reduce.h`,
  `src/montgomery_ref_compat.h`
  Copied C reference path plus the smallest compatibility layer needed to
  compile it standalone

- `src/*.S`
  Copied NEON assembly and placeholder SVE2 assembly

- `test.c`
  Functional comparison harness for C vs NEON vs SVE2

## Build Requirements

- Linux
- AArch64
- A compiler that accepts `-march=armv8-a+sve2`
- CMake 3.16 or newer

The CMake project is intentionally strict and fails during configure on
unsupported hosts. This is by design because the default build must include the
SVE2 objects.

## Build

```sh
cmake -S . -B build
cmake --build build -j
```

## Test

```sh
cmake --build build --target run
```

Or run the binary directly:

```sh
./build/test_montgomery_ops
```

## Current Expected Test Result

Today the default test run is expected to fail.

- C vs NEON should pass
- C vs SVE2 should fail
- The test binary should exit non-zero

This is intentional. The `_sve2` files are placeholder implementations that
write a deterministic impossible coefficient so the harness clearly shows where
real SVE2 work is still missing.

## Next Step

Replace the placeholder `_sve2` assembly files with real implementations while
keeping the existing test harness unchanged. Once the SVE2 functions match the
C reference and NEON behavior, the same binary becomes the correctness check.
