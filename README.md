# montgomery_ops

`montgomery_ops` is a small standalone sandbox for the Montgomery-related
polynomial operations used by `mldsa-native`.

It exists to make SVE2 development easier without building the full ML-DSA
tree. The codebase keeps three implementations side by side:

- Exact copied C reference logic from `mldsa-native/src`
- Copied AArch64 NEON assembly from `mldsa/src/native/aarch64/src`
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

## Montgomery NEON 学习笔记

这一节只讨论仓库里现有的 4 个 NEON 实现:

- `src/pointwise_montgomery.S`
- `src/mld_polyvecl_pointwise_acc_montgomery_l4.S`
- `src/mld_polyvecl_pointwise_acc_montgomery_l5.S`
- `src/mld_polyvecl_pointwise_acc_montgomery_l7.S`

### 4 个函数的共同骨架

这 4 个函数的主体结构其实完全一致, 差别只在于每个系数块要累加多少个多项式乘积:

- `poly_pointwise_montgomery`
  每个输出系数只做一次乘法, 然后立刻做一次 Montgomery reduction。
- `polyvecl_pointwise_acc_montgomery_l4`
  先做 1 次 `pmull`, 再做 3 次 `pmlal`, 最后统一 reduction。
- `polyvecl_pointwise_acc_montgomery_l5`
  先做 1 次 `pmull`, 再做 4 次 `pmlal`, 最后统一 reduction。
- `polyvecl_pointwise_acc_montgomery_l7`
  先做 1 次 `pmull`, 再做 6 次 `pmlal`, 最后统一 reduction。

每次循环都会并行处理 16 个 32-bit 系数:

- 4 个 `q` 寄存器
- 每个 `q` 寄存器装 4 个 `int32_t`
- 所以一次 loop 是 `4 x 4 = 16` 个系数

`count` 初值是 `MLDSA_N / 4`, 但每轮循环最后做的是 `subs count, count, #4`。原因不是每轮只处理 4 个系数, 而是每轮处理了 4 组 `q` 寄存器, 也就是 16 个系数, 因而一次减掉 4 个 "`q` 组"。

### `load_polys` 宏到底在做什么

`l4/l5/l7` 三个内积版本共享同一个 `load_polys` 宏:

```asm
.macro load_polys p0, p1, p2, p3, ptr, idx
.if \idx == 0
    ldr \p1, [\ptr, #1*16]
    ldr \p2, [\ptr, #2*16]
    ldr \p3, [\ptr, #3*16]
    ldr \p0, [\ptr], #4*16
.else
    ldr \p0, [\ptr, #(1024*\idx-4*16)]
    ldr \p1, [\ptr, #(1024*\idx-3*16)]
    ldr \p2, [\ptr, #(1024*\idx-2*16)]
    ldr \p3, [\ptr, #(1024*\idx-1*16)]
.endif
.endm
```

这里的关键点是:

- 一个多项式有 `256` 个系数。
- 每个系数是 `int32_t`, 所以一个多项式占 `256 x 4 = 1024` 字节。
- `idx == 0` 时, `ldr q0, [ptr], #64` 在取完当前多项式的当前 16 个系数后, 顺手把指针前移 64 字节。
- `idx > 0` 时, `1024*idx - 64` 这类偏移, 表示 "从已经前移过 64 字节的 `ptr` 回退 64 字节, 再跳到第 `idx` 个多项式的同一系数块"。

换句话说, `load_polys` 不是顺序扫完整个多项式向量, 而是在同一个系数块位置上, 在不同多项式之间横向取数据, 这样才能做点积累加。

### `montgomery_reduce_long` 宏的真实含义

4 个 NEON 文件都共用同一个 reduction 宏:

```asm
.macro montgomery_reduce_long res, inl, inh
    uzp1   \res\().4s, \inl\().4s, \inh\().4s
    mul    \res\().4s, \res\().4s, modulus_twisted.4s
    smlsl  \inl\().2d, \res\().2s, modulus.2s
    smlsl2 \inh\().2d, \res\().4s, modulus.4s
    uzp2   \res\().4s, \inl\().4s, \inh\().4s
.endm
```

它对应的标量公式就是:

```c
uint32_t t = (uint32_t)a * QINV;
int64_t r = a - (int64_t)(int32_t)t * Q;
return (int32_t)(r >> 32);
```

按指令拆开看:

- `uzp1 res.4s, inl.4s, inh.4s`
  从 4 个 64-bit 积里抽出每个积的低 32 位, 形成 `[a0_lo, a1_lo, a2_lo, a3_lo]`。
- `mul res.4s, res.4s, modulus_twisted.4s`
  做 32-bit 逐 lane 乘法, 只保留低 32 位, 正好等价于模 `2^32` 的乘法, 也就是 `t = (uint32_t)a * QINV`。
- `smlsl` / `smlsl2`
  以 64-bit lane 为目标累加器, 做 `in -= (int64_t)t * Q`。前者处理前两项, 后者处理后两项。
- `uzp2 res.4s, inl.4s, inh.4s`
  此时每个 64-bit lane 的低 32 位已经被消成 0, 取高 32 位就等价于 `r >> 32` 的结果。

也就是说, 这个宏的核心不是 "重新打包数据", 而是:

1. 抽出 64-bit 积的低 32 位做 `t`
2. 用 `t * Q` 消掉 64-bit 积的低半部分
3. 直接把剩下的高 32 位当作 Montgomery reduction 结果

### 为什么可以先累加再统一 reduction

每个文件里的注释都在跟踪边界:

- `poly_pointwise_montgomery` 的单次乘积满足 `81q^2 < qR/2`
- `l4` 最终累加满足 `36q^2 < qR/2`
- `l5` 最终累加满足 `45q^2 < qR/2`
- `l7` 最终累加满足 `63q^2 < qR/2`

这里 `R = 2^32`。这些上界说明 64-bit 累加仍然落在 Montgomery reduction 允许的范围内, 所以内积版本可以把多个乘积先加完, 最后再统一做一次 reduction, 不需要每加一次就 reduce 一次。

### 指令速查: NEON 到 SVE / SVE2

| NEON 指令 | 在这 4 个函数里的作用 | SVE 写法 | SVE2 写法 | 迁移时要注意什么 |
| --- | --- | --- | --- | --- |
| `movz` / `movk` | 先在 `wtmp` 里拼出 `Q` 和 `QINV` 常量 | 同一条 A64 标量指令 | 同一条 A64 标量指令 | 这两条不是 SIMD 指令, SVE/SVE2 下照用 |
| `dup vX.4s, wtmp` | 把 `Q` 或 `QINV` 广播到所有 32-bit lane | `dup zX.s, wtmp` | `dup zX.s, wtmp` | SVE 广播长度由硬件 `VL` 决定 |
| `ldr qX, [ptr, ...]` / `str qX, [ptr, ...]` | 读写 4 个 32-bit 系数 | `ld1w { zX.s }, p0/z, [ptr, ...]` / `st1w { zX.s }, p0, [ptr, ...]` | 同 SVE | SVE 需要 predicate; VLA 代码通常配合 `whilelt` 使用 |
| `smull` / `smull2` | 32-bit 乘 32-bit, 扩成 64-bit; `smull` 取低半, `smull2` 取高半 | 没有单条等价指令, 一般写成 `sunpk*` 加 `mul` | `smullb` / `smullt` | SVE2 的 `b/t` 处理的是交错的 even/odd 元素, 不是 NEON 的 low/high half |
| `smlal` / `smlal2` | 在已有 64-bit 累加器上继续加乘积 | `sunpk*` 加 `mla` | `smlalb` / `smlalt` | 和上一条一样, SVE2 的 lane 次序有变化 |
| `uzp1` | 从 64-bit lane 里抽低 32 位 | `uzp1 zT.s, zA.s, zB.s` | 同 SVE | 如果前面用了 `smullb/smullt`, 常常要先用 `trn1/trn2` 或 `zip1/zip2` 恢复 lane 顺序 |
| `mul vT.4s, vT.4s, vQinv.4s` | 计算 `t = low32(a) * QINV (mod 2^32)` | `mul zT.s, p0/m, zT.s, zQinv.s` | 同 SVE | predicated `mul` 的目标寄存器需要同时作为第一个源操作数 |
| `smlsl` / `smlsl2` | 做 `acc -= (int64_t)t * Q` | 把 `t` 和 `Q` 先扩成 64-bit, 再用 `mls` | `smlslb` / `smlslt` | SVE2 仍然有 even/odd 次序问题 |
| `uzp2` | 从 64-bit lane 里抽高 32 位, 也就是 reduction 输出 | `uzp2 zT.s, zA.s, zB.s` | 同 SVE | 前提是低 32 位已经被 `t * Q` 消成 0 |
| `subs` + `cbnz` | 固定 128-bit 宽度的计数循环 | 常见写法是 `whilelt` + `incw` + 条件分支 | 同 SVE | 真正的 SVE 版本最好写成 VLA, 不要把 `q` 寄存器个数硬编码成 4 |

### 关键指令的具体用法

#### `movz` / `movk` / `dup`

这组指令负责准备常量:

```asm
movz wtmp, #0xe001
movk wtmp, #0x7f, lsl #16
dup  modulus.4s, wtmp
```

含义是:

- 先用 `movz`/`movk` 在通用寄存器里拼出 32-bit 常量
- 再用 `dup` 广播到向量寄存器

移到 SVE/SVE2 时, 典型写法就是:

```asm
movz w3, #0xe001
movk w3, #0x7f, lsl #16
dup  z0.s, w3
```

#### `ldr q` / `str q`

NEON 版本默认 128-bit 固定宽度, 所以一次装 4 个 `int32_t`:

```asm
ldr q_a_1, [a_ptr, #1*16]
ldr q_a_2, [a_ptr, #2*16]
ldr q_a_3, [a_ptr, #3*16]
ldr q_a_0, [a_ptr], #4*16
```

SVE/SVE2 不再假设固定 4 lanes, 一般写成:

```asm
whilelt p0.s, x_idx, x_n
ld1w { z0.s }, p0/z, [x1, x_idx, lsl #2]
st1w { z0.s }, p0,   [x0, x_idx, lsl #2]
incw x_idx
```

用法要点:

- `whilelt p0.s, x_idx, x_n` 生成当前有效 lanes 的谓词
- `ld1w ... p0/z` 只装活跃 lanes, 非活跃 lanes 置 0
- `st1w ... p0` 只写活跃 lanes
- `incw x_idx` 让 `x_idx` 按当前硬件的 32-bit lane 数前进

#### `smull` / `smull2`

NEON 宏:

```asm
smull  c_0_lo.2d, a_0.2s, b_0.2s
smull2 c_0_hi.2d, a_0.4s, b_0.4s
```

用法含义:

- `smull` 处理 `a_0` / `b_0` 的低半部分 2 个 32-bit lanes
- `smull2` 处理高半部分 2 个 32-bit lanes
- 结果扩成 64-bit lane, 供后面的 Montgomery reduction 使用

SVE 没有单条 32->64 的 signed multiply-long 对应物, 典型写法是先扩再乘:

```asm
sunpklo zA0.d, zA.s
sunpkhi zA1.d, zA.s
sunpklo zB0.d, zB.s
sunpkhi zB1.d, zB.s
mul     zP0.d, p1/m, zA0.d, zB0.d
mul     zP1.d, p1/m, zA1.d, zB1.d
```

这里的 `p1` 应该理解成与 `.d` 结果匹配的 predicate, 也就是 widened 之后那一组 64-bit lanes 的有效掩码。

SVE2 则有更接近的写法:

```asm
smullb zP_even.d, zA.s, zB.s
smullt zP_odd.d,  zA.s, zB.s
```

但这里必须记住一个最重要的差异:

- NEON 的 widen-long 看的是 low half / high half
- SVE2 的 `*b` / `*t` 看的是交错的 even / odd 元素

所以 SVE2 虽然有近似一对一的指令名, 但 lane 排列不是直接兼容的。

#### `smlal` / `smlal2`

NEON 用它们把内积的后续乘积加到 64-bit 累加器里:

```asm
smlal  c_0_lo.2d, a_0.2s, b_0.2s
smlal2 c_0_hi.2d, a_0.4s, b_0.4s
```

SVE 里没有单条 signed multiply-accumulate-long, 通常还是先扩再 `mla`:

```asm
mla zP0.d, p1/m, zA0.d, zB0.d
mla zP1.d, p1/m, zA1.d, zB1.d
```

SVE2 的最近对应物是:

```asm
smlalb zP_even.d, zA.s, zB.s
smlalt zP_odd.d,  zA.s, zB.s
```

用途和 `smullb/smullt` 一样, 只是把结果加到已有 64-bit 累加器而不是覆盖。

#### `uzp1`

在这几个 Montgomery 函数里, `uzp1` 不是普通的重排技巧, 而是 reduction 的第一步:

```asm
uzp1 res.4s, inl.4s, inh.4s
```

它把 4 个 64-bit lane 的低 32 位抽出来, 正好得到做 `t = low32(a) * QINV` 所需的 `t` 输入。

SVE/SVE2 仍然可以直接用同名指令:

```asm
uzp1 zT.s, zAcc0.s, zAcc1.s
```

如果前面乘法用了 SVE2 的 `smullb/smullt`, 常见写法会先用 `trn1` / `trn2` 把 even/odd 结果重新排成 Montgomery reduction 需要的顺序, 再做 `uzp1`。

#### `mul`

这条在 reduction 宏里很关键:

```asm
mul res.4s, res.4s, modulus_twisted.4s
```

这里并不是要得到完整的 64-bit 乘积, 恰恰相反, 它只需要保留 32-bit 低半部分, 因为 Montgomery reduction 只关心模 `2^32` 的结果。

SVE/SVE2 写法:

```asm
mul zT.s, p0/m, zT.s, zQinv.s
```

要点是:

- 目标寄存器 `zT` 同时是第一个源操作数
- `p0/m` 表示只更新活跃 lanes, 非活跃 lanes 保持原值

#### `smlsl` / `smlsl2`

这两条负责做 `acc -= t * Q`:

```asm
smlsl  inl.2d, res.2s, modulus.2s
smlsl2 inh.2d, res.4s, modulus.4s
```

SVE 版通常要先把 32-bit 的 `t` 和 `Q` 扩成 64-bit, 再做 `mls`:

```asm
sunpklo zT0.d, zT.s
sunpkhi zT1.d, zT.s
mls     zAcc0.d, p1/m, zT0.d, zQ0.d
mls     zAcc1.d, p1/m, zT1.d, zQ1.d
```

SVE2 有直接的 long subtract:

```asm
smlslb zAcc_even.d, zT.s, zQ.s
smlslt zAcc_odd.d,  zT.s, zQ.s
```

它们分别处理 even 和 odd 的 32-bit 元素并扩成 64-bit 后做减法。

#### `uzp2`

`uzp2` 是 reduction 的收尾:

```asm
uzp2 res.4s, inl.4s, inh.4s
```

此时每个 64-bit lane 的低 32 位已经是 0, 所以 `uzp2` 取出的高 32 位就是最终的 Montgomery output。

SVE/SVE2 仍然可以直接写:

```asm
uzp2 zR.s, zAcc0.s, zAcc1.s
```

如果 64-bit 累加器仍处在 even/odd 排列, 同样需要先做一次 `trn1` / `trn2` 或 `zip1` / `zip2` 来恢复期望顺序。

#### `subs` / `cbnz` 与 SVE 的 loop 习惯

NEON 版本最后是固定宽度循环:

```asm
subs count, count, #4
cbnz count, poly_pointwise_montgomery_loop_start
```

SVE/SVE2 更常见的写法是 predicate-driven loop:

```asm
mov     x_idx, #0
whilelt p0.s, x_idx, x_n
loop:
    ...
    incw   x_idx
    whilelt p0.s, x_idx, x_n
    b.mi   loop
```

这种写法的好处是:

- 自动适配不同 `VL`
- 天然处理尾部不足一个完整向量的情况
- 不需要再人为维护 "一次处理几个 `q` 寄存器" 这种固定宽度计数

### 迁移到 SVE / SVE2 时最容易踩的坑

- 最大的坑不是指令名, 而是 lane 顺序。Arm 的 SVE2 入门文档明确指出, transformed widen/narrow/pairwise 指令按 even/odd 交错元素工作, 而 Neon 按 low/high half 工作。
- 因此, `smull -> smullb`、`smull2 -> smullt` 这种机械替换在语义上并不完全等价。
- 如果要严格复刻当前 NEON 宏的中间结果布局, 通常要在 `smullb/smullt` 或 `smlalb/smlalt` 之后增加 `trn1` / `trn2` 或 `zip1` / `zip2`。
- 如果只追求最终数值正确, 也可以接受 even/odd 布局, 但整个 reduction 链条都必须一致地按 even/odd 去设计, 不能中途再按 NEON 的 low/high 假设去取 lane。
- 纯 SVE 没有直接的 32->64 signed multiply-long 和 multiply-accumulate-long 指令, 所以最稳妥的思路是 "先 `sunpk*`, 再 `mul` / `mla` / `mls`"。

### 一句话总结

这 4 个 NEON 实现的本质都是:

1. 按 16 个系数一组装载
2. 在 64-bit lane 中做乘法或点积累加
3. 用 `uzp1 + mul + smlsl/smlsl2 + uzp2` 完成 Montgomery reduction
4. 写回 32-bit 结果

从 NEON 迁到 SVE/SVE2 时, 真正要守住的不是某一条指令的名字, 而是:

- 64-bit 中间值的 lane 排列
- `low32(a)` 和 `high32(r)` 被抽取时的次序
- SVE2 widen 指令的 even/odd 语义和 NEON low/high 语义并不相同

## Source Provenance

- C reference logic comes from upstream `mldsa/src/poly.c`,
  `mldsa/src/polyvec.c`, `mldsa/src/reduce.h`, and `mldsa/src/ct.h`
- NEON assembly comes from upstream `mldsa/src/native/aarch64/src`
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
