void pyop3_kernel(double const *__restrict__ array_0,
                  int64_t const *__restrict__ array_1,
                  double *__restrict__ array_2)
{
  double t_0[6];
  double t_1[1];

  for (int32_t i_0 = 0; i_0 <= 4; ++i_0)
  {
    for (int32_t i_2 = 0; i_2 <= 2; ++i_2)
      for (int32_t i_1 = 0; i_1 <= 1; ++i_1)
        t_0[i_1*3+i_2] = array_0[array_1[2*i_0+i_1]*3+i_2];
    t_1[0] = 0.0;
    kernel(&(t_0[0]), &(t_1[0]));
    array_2[i_0] = array_2[i_0] + t_1[0];
  }
}
