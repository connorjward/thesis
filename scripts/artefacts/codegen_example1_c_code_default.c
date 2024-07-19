#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>


int32_t cmpfunc(const void * a, const void * b) {
   return ( *(int32_t*)a - *(int32_t*)b );
}
#include <petsc.h>

void mykernel(double *__restrict__ array_0)
{
  for (int32_t i_1 = 0; i_1 <= 1; ++i_1)
    for (int32_t i_0 = 0; i_0 <= 2; ++i_0)
      array_0[1 + 6 * i_0 + i_1] = 666.0;

}