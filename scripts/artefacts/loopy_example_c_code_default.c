#include <stdint.h>

void loopy_kernel(double const *__restrict__ x, double *__restrict__ y, int64_t const n)
{
  for (int32_t i = 0; i <= -1 + n; ++i)
    y[i] = 2.0 * x[i];
}