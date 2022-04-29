#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

// Make sure these dimensions are consistent with the source tables!!!
#define DIM_POP2_NIONTBL_R (int) 13
#define DIM_POP2_NIONTBL_C (int) 2

// Source tables to read and interpolate
#define NION_POP2_SED_FILENAME (const char *) "../External_tables/stellar_spectra_bpass/nion_continuous.dat"

static double rdata_nion_Z[DIM_POP2_NIONTBL_R], rdata_nion_N[DIM_POP2_NIONTBL_R];
static gsl_spline *spline_nion;
static gsl_interp_accel *xacc_nion;

int init_pop2_nion_table() {
  int i, ii;
  float nion_table_pop2[DIM_POP2_NIONTBL_R][DIM_POP2_NIONTBL_C];
  FILE *F;

  if (!(F = fopen(NION_POP2_SED_FILENAME, "r"))) {
    fprintf(stderr, "Unable to open file: nion_continuous.dat for reading\nAborting\n");
    return -1;
  }
  for (i = 0; i < DIM_POP2_NIONTBL_R; i++) {
    fscanf(F, "%*[^\n]\n");
    fscanf(F, "%e %e", &nion_table_pop2[i][0], &nion_table_pop2[i][1]);
  }
  fclose(F);

  // Populate Z axis
  for (ii = 0; ii < DIM_POP2_NIONTBL_R; ii++) {
    rdata_nion_Z[ii] = nion_table_pop2[ii][0];
    rdata_nion_N[ii] = nion_table_pop2[ii][1];
  }

  const gsl_interp_type *T = gsl_interp_linear;

  size_t nx = sizeof(rdata_nion_Z) / sizeof(rdata_nion_Z[0]);

  spline_nion = gsl_spline_alloc(T, nx);
  xacc_nion = gsl_interp_accel_alloc();

  gsl_spline_init(spline_nion, rdata_nion_Z, rdata_nion_N, nx);
}


float get_nion(float metallicity)
{
  return gsl_spline_eval(spline_nion, metallicity, xacc_nion);
}
