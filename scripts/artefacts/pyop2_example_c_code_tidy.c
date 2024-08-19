void pyop2_kernel(int32_t const start, int32_t const end,
                  double const *__restrict__ dat0,
                  double *__restrict__ dat1,
                  int32_t const *__restrict__ map0)
{
  double t0[3 * 2];
  double t1[3 * 2];

  for (int32_t n = start; n <= -1 + end; ++n)
  {
    for (int32_t i10 = 0; i10 <= 2; ++i10)
      for (int32_t i11 = 0; i11 <= 1; ++i11)
        t0[2 * i10 + i11] = dat0[2 * map0[3 * n + i10] + i11];
    for (int32_t i13 = 0; i13 <= 2; ++i13)
      for (int32_t i14 = 0; i14 <= 1; ++i14)
        t1[2 * i13 + i14] = 0.0;
    local_kernel(&(t0[0]), &(t1[0]));
    for (int32_t i7 = 0; i7 <= 1; ++i7)
      for (int32_t i6 = 0; i6 <= 2; ++i6)
        dat1[2 * map0[3 * n + i6] + i7] = t1[2 * i6 + i7];
  }
}
