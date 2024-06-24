void loopy_kernel(double const *x,
                  double *y,
                  int64_t const n)
{
  for (int32_t i = 0; i <= -1 + n; ++i)
    y[i] = 2.0 * x[i];
}
