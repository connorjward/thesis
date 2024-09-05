void pyop3_loop(...)
{
  int32_t p_0;
  double t_0[1];
  double t_1[1];

  for (int32_t i_0 = 0; i_0 <= end; ++i_0)
  {
    p_0 = array_0[i_0];
    for (int32_t i_1 = 0; i_1 <= -1 + p_0; ++i_1)
    {
      t_0[0] = array_3[array_4[array_1[array_2[i_0] + i_1]]];
      t_1[0] = array_5[array_6[i_0]];
      max_kernel(&(t_0[0]), &(t_1[0]));
      array_5[array_6[i_0]] = t_1[0];
    }
  }
}
