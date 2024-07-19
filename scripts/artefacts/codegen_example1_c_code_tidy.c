void pyop3_kernel(double *__restrict__ array_0)
{
  for (int32_t i_1 = 0; i_1 <= 1; ++i_1)
    for (int32_t i_0 = 0; i_0 <= 2; ++i_0)
      array_0[1 + 6 * i_0 + i_1] = 666.0;
}
