#include <math.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <sys/wait.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"
#include "../Parameter_files/HEAT_PARAMS.H"
#include "../Parameter_files/SOURCES.H"

/*
  Program compute_line_fields reads previously computed boxes for metallicity and stellar mass
  and uses a look up table of SEDs and emissivities to compute radiation fields of various
  emission lines for each cell in the box.
*/

//make sure these are the same high and low as the selected driver file

//#define POP2_METALLICITY_SAMPLES (int) 24
#define DIM_POP2_METALLICITY (int) 41
#define DIM_POP2_ION_PARAM (int) 4
#define DIM_LUM (int) 164
#define POP3_ION_SAMPLES (int) 3
#define POP3_COLUMNS (int) 8

//#define HALPHA_POP2_FILENAME (const char *) "../External_tables/Halpha_table_pop2.dat"
//#define LYA_POP2_FILENAME (const char *) "../External_tables/Lya_table_pop2.dat"
#define LYA_POP2_TABLE_FILENAME (const char *) "../External_tables/Lya_PopII_table_cont.dat"
#define HA_POP2_TABLE_FILENAME (const char *) "../External_tables/Ha_PopII_table_cont.dat"
#define O3_POP2_TABLE_FILENAME (const char *) "../External_tables/O3_PopII_table_cont.dat"
#define O2_POP2_TABLE_FILENAME (const char *) "../External_tables/O2_PopII_table_cont.dat"
//#define POP3_FILENAME (const char *) "../External_tables/pop3_luminosity.dat"

static double rdata_Z[DIM_POP2_METALLICITY], rdata_U[DIM_POP2_ION_PARAM];
static double rdata_LUM_Lya_flat[DIM_LUM], rdata_LUM_Ha_flat[DIM_LUM], rdata_LUM_O3_flat[DIM_LUM], rdata_LUM_O2_flat[DIM_LUM];

static gsl_spline2d *spline_Lya_lum, *spline_Ha_lum, *spline_O3_lum, *spline_O2_lum;
static gsl_interp_accel *xacc_lum;
static gsl_interp_accel *yacc_lum;

int init_pop2_full_tables(float Lya_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM],
                     float Ha_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM],
                     float O3_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM],
                     float O2_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM]) {
  int i, ii, jj, kk;
  FILE *F;

  //initialize Lya table
  if (!(F = fopen(LYA_POP2_TABLE_FILENAME, "r"))) {
    fprintf(stderr, "Unable to open file: Lya_PopII_table_cont.dat for reading\nAborting\n");
    return -1;
  }
  for (i = 0; i < DIM_POP2_METALLICITY; i++) {
    fscanf(F, "%*[^\n]\n");
    fscanf(F, "%e %e %e %e", &Lya_table_pop2[i][0], &Lya_table_pop2[i][1], &Lya_table_pop2[i][2], &Lya_table_pop2[i][3]);
  }
  fclose(F);

  //initialize Ha table
  if (!(F = fopen(HA_POP2_TABLE_FILENAME, "r"))) {
    fprintf(stderr, "Unable to open file: Ha_PopII_table_cont.dat for reading\nAborting\n");
    return -1;
  }
  for (i = 0; i < DIM_POP2_METALLICITY; i++) {
    fscanf(F, "%*[^\n]\n");
    fscanf(F, "%e %e %e %e", &Ha_table_pop2[i][0], &Ha_table_pop2[i][1], &Ha_table_pop2[i][2], &Ha_table_pop2[i][3]);
  }
  fclose(F);

  //initialize O3 table
  if (!(F = fopen(O3_POP2_TABLE_FILENAME, "r"))) {
    fprintf(stderr, "Unable to open file: O3_PopII_table_cont.dat for reading\nAborting\n");
    return -1;
  }
  for (i = 0; i < DIM_POP2_METALLICITY; i++) {
    fscanf(F, "%*[^\n]\n");
    fscanf(F, "%e %e %e %e", &O3_table_pop2[i][0], &O3_table_pop2[i][1], &O3_table_pop2[i][2], &O3_table_pop2[i][3]);
  }
  fclose(F);

  //initialize O2 table
  if (!(F = fopen(O2_POP2_TABLE_FILENAME, "r"))) {
    fprintf(stderr, "Unable to open file: O2_PopII_table_cont.dat for reading\nAborting\n");
    return -1;
  }
  for (i = 0; i < DIM_POP2_METALLICITY; i++) {
    fscanf(F, "%*[^\n]\n");
    fscanf(F, "%e %e %e %e", &O2_table_pop2[i][0], &O2_table_pop2[i][1], &O2_table_pop2[i][2], &O2_table_pop2[i][3]);
  }
  fclose(F);

  // Populate Z axis
  for (ii = 0; ii < DIM_POP2_METALLICITY; ii++) {
    rdata_Z[ii] = -3.+0.1*ii;
  }
  // Populate U axis
  for (jj = 0; jj < DIM_POP2_ION_PARAM; jj++) {
    rdata_U[jj] = -4.+1.*jj;
  }
  // Populate flat luminosity array
  for (jj = 0; jj < DIM_POP2_ION_PARAM; jj++) {
    for (ii = 0; ii < DIM_POP2_METALLICITY; ii++) {
      int kk = ii + jj*DIM_POP2_METALLICITY;
      rdata_LUM_Lya_flat[kk] = Lya_table_pop2[ii][jj];
      rdata_LUM_Ha_flat[kk] = Ha_table_pop2[ii][jj];
      rdata_LUM_O3_flat[kk] = O3_table_pop2[ii][jj];
      rdata_LUM_O2_flat[kk] = O2_table_pop2[ii][jj];
    }
  }

  const gsl_interp2d_type *T = gsl_interp2d_bilinear;

  size_t nx = sizeof(rdata_Z) / sizeof(rdata_Z[0]);
  size_t ny = sizeof(rdata_U) / sizeof(rdata_U[0]);

  spline_Lya_lum = gsl_spline2d_alloc(T, nx, ny);
  spline_Ha_lum = gsl_spline2d_alloc(T, nx, ny);
  spline_O3_lum = gsl_spline2d_alloc(T, nx, ny);
  spline_O2_lum = gsl_spline2d_alloc(T, nx, ny);
  xacc_lum = gsl_interp_accel_alloc();
  yacc_lum = gsl_interp_accel_alloc();

  gsl_spline2d_init(spline_Lya_lum, rdata_Z, rdata_U, rdata_LUM_Lya_flat, nx, ny);
  gsl_spline2d_init(spline_Ha_lum, rdata_Z, rdata_U, rdata_LUM_Ha_flat, nx, ny);
  gsl_spline2d_init(spline_O3_lum, rdata_Z, rdata_U, rdata_LUM_O3_flat, nx, ny);
  gsl_spline2d_init(spline_O2_lum, rdata_Z, rdata_U, rdata_LUM_O2_flat, nx, ny);

  return 0;
}

float get_Lya_luminosity(float Lya_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM], float logZ, int logU) {
  return gsl_spline2d_eval(spline_Lya_lum, logZ, logU, xacc_lum, yacc_lum);
}

float get_Ha_luminosity(float Ha_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM], float logZ, int logU) {
  return gsl_spline2d_eval(spline_Ha_lum, logZ, logU, xacc_lum, yacc_lum);
}

float get_O3_luminosity(float O3_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM], float logZ, int logU) {
  return gsl_spline2d_eval(spline_O3_lum, logZ, logU, xacc_lum, yacc_lum);
}

float get_O2_luminosity(float O2_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM], float logZ, int logU) {
  return gsl_spline2d_eval(spline_O2_lum, logZ, logU, xacc_lum, yacc_lum);
}

/*
int main(int argc, char const *argv[])
{
  float Lya_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM];
  float Ha_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM];
  float O3_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM];
  float O2_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM];

  init_pop2_tables(Lya_table_pop2, Ha_table_pop2, O3_table_pop2, O2_table_pop2);

  printf("CHECK Lya true value: %f, %f, %f\n", rdata_Z[20], rdata_U[2], Lya_table_pop2[20][2]);
  printf("CHECK Lya interp value: %f\n", get_luminosity(Lya_table_pop2, -1., -2.));

  printf("CHECK Ha true value: %f, %f, %f\n", rdata_Z[20], rdata_U[2], Ha_table_pop2[20][2]);
  printf("CHECK Ha interp value: %f\n", get_luminosity(Ha_table_pop2, -1., -2.));

  printf("CHECK O3 true value: %f, %f, %f\n", rdata_Z[20], rdata_U[2], O3_table_pop2[20][2]);
  printf("CHECK O3 interp value: %f\n", get_luminosity(O3_table_pop2, -1., -2.));

  printf("CHECK O2 true value: %f, %f, %f\n", rdata_Z[20], rdata_U[2], O2_table_pop2[20][2]);
  printf("CHECK O2 interp value: %f\n", get_luminosity(O2_table_pop2, -1., -2.));

  return 0;
}
*/
