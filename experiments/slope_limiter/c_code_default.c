#include <stdbool.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>


int32_t cmpfunc(const void * a, const void * b) {
   return ( *(int32_t*)a - *(int32_t*)b );
}
#include <petsc.h>
#include <petsclog.h>

// Prepare a dummy event so that things compile. This is overwridden using
// the object file.
PetscLogEvent id_pyop3_loop = -1;

static void loopy_kernel(double const *__restrict__ in, double *__restrict__ out);
static void loopy_kernel(double const *__restrict__ in, double *__restrict__ out)
{
  out[0] = ((out[0] > in[0]) ? out[0] : in[0]);
}

void pyop3_loop(int32_t const *__restrict__ array_0, int32_t const *__restrict__ array_1, int32_t const *__restrict__ array_2, double const *__restrict__ array_3, int32_t const *__restrict__ array_4, double *__restrict__ array_5, int32_t const *__restrict__ array_6)
{
  int32_t j_0;
  int32_t j_1;
  int32_t j_2;
  int32_t p_0;
  double t_0[1];
  double t_1[1];

  PetscLogEventBegin(id_pyop3_loop, 0, 0, 0, 0);
  for (int32_t i_0 = 0; i_0 <= 3; ++i_0)
  {
    p_0 = array_0[i_0];
    for (int32_t i_1 = 0; i_1 <= -1 + p_0; ++i_1)
    {
      j_0 = 0;
      t_0[j_0] = array_3[array_4[array_1[array_2[i_0] + i_1]] + 0];
      j_1 = 0;
      t_1[j_1] = array_5[array_6[i_0] + 0];
      loopy_kernel(&(t_0[0]), &(t_1[0]));
      j_2 = 0;
      array_5[array_6[i_0] + 0] = t_1[j_2];
    }
  }

  PetscLogEventEnd(id_pyop3_loop, 0, 0, 0, 0);
}
