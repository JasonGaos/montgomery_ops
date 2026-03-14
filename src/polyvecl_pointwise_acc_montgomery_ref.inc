MLD_STATIC_TESTABLE void mld_polyvecl_pointwise_acc_montgomery_c(
    mld_poly *w, const mld_polyvecl *u, const mld_polyvecl *v)
__contract__(
  requires(memory_no_alias(w, sizeof(mld_poly)))
  requires(memory_no_alias(u, sizeof(mld_polyvecl)))
  requires(memory_no_alias(v, sizeof(mld_polyvecl)))
  requires(forall(l0, 0, MLDSA_L,
                  array_bound(u->vec[l0].coeffs, 0, MLDSA_N, 0, MLDSA_Q)))
  requires(forall(l1, 0, MLDSA_L,
    array_abs_bound(v->vec[l1].coeffs, 0, MLDSA_N, MLD_NTT_BOUND)))
  assigns(memory_slice(w, sizeof(mld_poly)))
  ensures(array_abs_bound(w->coeffs, 0, MLDSA_N, MLDSA_Q))
)
{
  unsigned int i, j;
  mld_assert_bound_2d(u->vec, MLDSA_L, MLDSA_N, 0, MLDSA_Q);
  mld_assert_abs_bound_2d(v->vec, MLDSA_L, MLDSA_N, MLD_NTT_BOUND);
  for (i = 0; i < MLDSA_N; i++)
  __loop__(
    assigns(i, j, memory_slice(w, sizeof(mld_poly)))
    invariant(i <= MLDSA_N)
    invariant(array_abs_bound(w->coeffs, 0, i, MLDSA_Q))
  )
  {
    int64_t t = 0;
    int32_t r;
    for (j = 0; j < MLDSA_L; j++)
    __loop__(
      assigns(j, t)
      invariant(j <= MLDSA_L)
      invariant(t >= -(int64_t)j*(MLDSA_Q - 1)*(MLD_NTT_BOUND - 1))
      invariant(t <= (int64_t)j*(MLDSA_Q - 1)*(MLD_NTT_BOUND - 1))
    )
    {
      t += (int64_t)u->vec[j].coeffs[i] * v->vec[j].coeffs[i];
    }

    r = mld_montgomery_reduce(t);
    w->coeffs[i] = r;
  }

  mld_assert_abs_bound(w->coeffs, MLDSA_N, MLDSA_Q);
}
