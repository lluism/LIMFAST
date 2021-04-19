#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#define DIM_POP2_METALLICITY_SAMPLES (int) 41
#define DIM_POP2_ION_PARAM_SAMPLES (int) 4
#define DIM_LUM (int) 164

#define LINELUM_POP2_FILENAME (const char *) "../External_tables/Ha_PopII_table_cont.dat"

static double rdata_Z[DIM_POP2_METALLICITY_SAMPLES], rdata_U[DIM_POP2_ION_PARAM_SAMPLES], rdata_LUM[DIM_POP2_METALLICITY_SAMPLES][DIM_POP2_ION_PARAM_SAMPLES], rdata_LUM_flat[DIM_LUM];

static gsl_spline2d *spline_lum;
static gsl_interp_accel *xacc_lum;
static gsl_interp_accel *yacc_lum;


int init_pop2_tbl(float rdata_LUM[DIM_POP2_METALLICITY_SAMPLES][DIM_POP2_ION_PARAM_SAMPLES]) {
  int i, ii, jj, kk;
  FILE *F;
  //initialize Halpha table
  if (!(F = fopen(LINELUM_POP2_FILENAME, "r"))) {
    fprintf(stderr, "Unable to open file: Halpha_table_pop2.dat for reading\nAborting\n");
    return -1;
  }
  for (i = 0; i < DIM_POP2_METALLICITY_SAMPLES; i++) {
    fscanf(F, "%*[^\n]\n");
    fscanf(F, "%e %e %e %e", &rdata_LUM[i][0], &rdata_LUM[i][1], &rdata_LUM[i][2], &rdata_LUM[i][3]);
  }
  fclose(F);
  // Populate Z axis
  for (ii = 0; ii < DIM_POP2_METALLICITY_SAMPLES; ii++) {
    rdata_Z[ii] = -3.+0.1*ii;
  }
  // Populate U axis
  for (jj = 0; jj < DIM_POP2_ION_PARAM_SAMPLES; jj++) {
    rdata_U[jj] = -4.+1.*jj;
  }
  for (jj = 0; jj < DIM_POP2_ION_PARAM_SAMPLES; jj++) {
    for (ii = 0; ii < DIM_POP2_METALLICITY_SAMPLES; ii++) {
      int kk = ii + jj*DIM_POP2_METALLICITY_SAMPLES;
      rdata_LUM_flat[kk] = rdata_LUM[ii][jj];
    }
  }

  const gsl_interp2d_type *T = gsl_interp2d_bilinear;

  size_t nx = sizeof(rdata_Z) / sizeof(rdata_Z[0]);
  size_t ny = sizeof(rdata_U) / sizeof(rdata_U[0]);

  spline_lum = gsl_spline2d_alloc(T, nx, ny);
  xacc_lum = gsl_interp_accel_alloc();
  yacc_lum = gsl_interp_accel_alloc();

  gsl_spline2d_init(spline_lum, rdata_Z, rdata_U, rdata_LUM_flat, nx, ny);

  return 0;
}

float get_pop2_lum(float logZ, float logU)
{
  return gsl_spline2d_eval(spline_lum, logZ, logU, xacc_lum, yacc_lum);
}


//int main(int argc, char const *argv[])
//{
//  float Halpha_table_2[DIM_POP2_METALLICITY_SAMPLES][DIM_POP2_ION_PARAM_SAMPLES];
//  init_pop2_tbl(Halpha_table_2);
//  //printf("CHECK: %e\n", get_mar_base(10.0, 1.520e9));
//  printf("CHECK CHECK: %f\n", Halpha_table_2[0][0]);
//  printf("CHECK CHECK: %f, %f, %f\n", rdata_Z[20], rdata_U[2], Halpha_table_2[20][2]);
//  printf("CHECK CHECK: %f\n", get_pop2_lum(-1., -2.));
//  return 0;
//}
