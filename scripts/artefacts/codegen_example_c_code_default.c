#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>


int32_t cmpfunc(const void * a, const void * b) {
   return ( *(int32_t*)a - *(int32_t*)b );
}
#include <petsc.h>

static void loopy_kernel(double const *__restrict__ x, double *__restrict__ y);
static void loopy_kernel(double const *__restrict__ x, double *__restrict__ y)
{
  for (int32_t i = 0; i <= 5; ++i)
    y[0] = y[0] + x[i];
}

void mykernel(double const *__restrict__ array_0, int64_t const *__restrict__ array_1, double *__restrict__ array_2)
{
  int32_t j_0;
  int32_t j_1;
  int32_t j_2;
  double t_0[6];
  double t_1[1];

  for (int32_t i_0 = 0; i_0 <= 4; ++i_0)
  {
    for (int32_t i_2 = 0; i_2 <= 2; ++i_2)
      for (int32_t i_1 = 0; i_1 <= 1; ++i_1)
      {
        j_0 = i_1 * 3 + i_2;
        t_0[j_0] = array_0[array_1[2 * i_0 + i_1] * 3 + i_2];
      }
    j_1 = 0;
    t_1[j_1] = 0.0;
    loopy_kernel(&(t_0[0]), &(t_1[0]));
    j_2 = 0;
    array_2[i_0] = array_2[i_0] + t_1[j_2];
  }

}