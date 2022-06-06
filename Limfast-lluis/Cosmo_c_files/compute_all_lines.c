#include <math.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <sys/wait.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

/*
  Program compute_all_lines reads separately computed look-up tables of line luminosities
  and interpolates them as a function of metallicity and ionization parameter. The resulting
  L(Z, U) relation is supplied to the SFR(M, z) or M_*(M, z) relation to compute the box of
  line intensity.
*/

// Make sure these dimensions are consistent with the source tables!!!
#define DIM_POP2_METALLICITY (int) 12
#define DIM_POP2_ION_PARAM (int) 7
#define DIM_LUM (int) 84
#define POP3_ION_SAMPLES (int) 3
#define POP3_COLUMNS (int) 8

// Source tables to read and interpolate
#define LYA_POP2_TABLE_FILENAME (const char *) "../External_tables/line_tables_bpass/LlineTbl_POP2_HI1216A_Continuous_ND.dat"
#define HA_POP2_TABLE_FILENAME (const char *) "../External_tables/line_tables_bpass/LlineTbl_POP2_HI6563A_Continuous_ND.dat"
#define HB_POP2_TABLE_FILENAME (const char *) "../External_tables/line_tables_bpass/LlineTbl_POP2_HI4861A_Continuous_ND.dat"
#define HE2_POP2_TABLE_FILENAME (const char *) "../External_tables/line_tables_bpass/LlineTbl_POP2_HeII1640A_Continuous_ND.dat"
#define O3_POP2_TABLE_FILENAME (const char *) "../External_tables/line_tables_bpass/LlineTbl_POP2_OIII5007A_Continuous_ND.dat"
#define O2S_POP2_TABLE_FILENAME (const char *) "../External_tables/line_tables_bpass/LlineTbl_POP2_OII3726A_Continuous_ND.dat"
#define O2L_POP2_TABLE_FILENAME (const char *) "../External_tables/line_tables_bpass/LlineTbl_POP2_OII3729A_Continuous_ND.dat"
#define C2_POP2_TABLE_FILENAME (const char *) "../External_tables/line_tables_bpass/LlineTbl_POP2_CII157m.dat"
#define CO10_POP2_TABLE_FILENAME (const char *) "../External_tables/line_tables_bpass/LlineTbl_POP2_CO2600m.dat"


double Z_VALUES[]={9.999e-6, 1.0e-4, 1.0e-3, 2.0e-3, 4.0e-3, 6.0e-3, 8.0e-3, 1.0e-2, 1.4e-2, 2.0e-2, 3.0e-2, 4.001e-2};
static double rdata_Z[DIM_POP2_METALLICITY], rdata_U[DIM_POP2_ION_PARAM];
static double rdata_LUM_LYA_flat[DIM_LUM], rdata_LUM_HA_flat[DIM_LUM], rdata_LUM_HB_flat[DIM_LUM], rdata_LUM_HE2_flat[DIM_LUM],
              rdata_LUM_O3_flat[DIM_LUM], rdata_LUM_O2S_flat[DIM_LUM], rdata_LUM_O2L_flat[DIM_LUM], rdata_LUM_C2_flat[DIM_LUM],
              rdata_LUM_CO10_flat[DIM_LUM];

static gsl_spline2d *spline_lya_lum, *spline_ha_lum, *spline_hb_lum, *spline_he2_lum, *spline_o3_lum,
                    *spline_o2s_lum, *spline_o2l_lum, *spline_c2_lum, *spline_co10_lum;
static gsl_interp_accel *xacc_lum;
static gsl_interp_accel *yacc_lum;

int init_pop2_full_tables(float lya_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM],
                          float ha_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM],
                          float hb_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM],
                          float he2_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM],
                          float o3_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM],
                          float o2s_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM],
                          float o2l_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM],
                          float c2_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM],
                          float co10_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM]) {
  int i, ii, jj, kk;
  FILE *F;

  //initialize Lya table
  if (!(F = fopen(LYA_POP2_TABLE_FILENAME, "r"))) {
    fprintf(stderr, "Unable to open file: LlineTbl_POP2_HI1216A_Continuous_ND.dat for reading\nAborting\n");
    return -1;
  }
  for (i = 0; i < DIM_POP2_METALLICITY; i++) {
    fscanf(F, "%*[^\n]\n");
    fscanf(F, "%e %e %e %e %e %e %e", &lya_table_pop2[i][0], &lya_table_pop2[i][1], &lya_table_pop2[i][2], &lya_table_pop2[i][3],
                                      &lya_table_pop2[i][4], &lya_table_pop2[i][5], &lya_table_pop2[i][6]);
  }
  fclose(F);

  //initialize Ha table
  if (!(F = fopen(HA_POP2_TABLE_FILENAME, "r"))) {
    fprintf(stderr, "Unable to open file: LlineTbl_POP2_HI6563A_Continuous_ND.dat for reading\nAborting\n");
    return -1;
  }
  for (i = 0; i < DIM_POP2_METALLICITY; i++) {
    fscanf(F, "%*[^\n]\n");
    fscanf(F, "%e %e %e %e %e %e %e", &ha_table_pop2[i][0], &ha_table_pop2[i][1], &ha_table_pop2[i][2], &ha_table_pop2[i][3],
                                      &ha_table_pop2[i][4], &ha_table_pop2[i][5], &ha_table_pop2[i][6]);
  }
  fclose(F);

  //initialize Hb table
  if (!(F = fopen(HB_POP2_TABLE_FILENAME, "r"))) {
    fprintf(stderr, "Unable to open file: LlineTbl_POP2_HI4861A_Continuous_ND.dat for reading\nAborting\n");
    return -1;
  }
  for (i = 0; i < DIM_POP2_METALLICITY; i++) {
    fscanf(F, "%*[^\n]\n");
    fscanf(F, "%e %e %e %e %e %e %e", &hb_table_pop2[i][0], &hb_table_pop2[i][1], &hb_table_pop2[i][2], &hb_table_pop2[i][3],
                                      &hb_table_pop2[i][4], &hb_table_pop2[i][5], &hb_table_pop2[i][6]);
  }
  fclose(F);

  //initialize HeII table
  if (!(F = fopen(HE2_POP2_TABLE_FILENAME, "r"))) {
    fprintf(stderr, "Unable to open file: LlineTbl_POP2_HeII1640A_Continuous_ND.dat for reading\nAborting\n");
    return -1;
  }
  for (i = 0; i < DIM_POP2_METALLICITY; i++) {
    fscanf(F, "%*[^\n]\n");
    fscanf(F, "%e %e %e %e %e %e %e", &he2_table_pop2[i][0], &he2_table_pop2[i][1], &he2_table_pop2[i][2], &he2_table_pop2[i][3],
                                      &he2_table_pop2[i][4], &he2_table_pop2[i][5], &he2_table_pop2[i][6]);
  }
  fclose(F);

  //initialize OIII table
  if (!(F = fopen(O3_POP2_TABLE_FILENAME, "r"))) {
    fprintf(stderr, "Unable to open file: LlineTbl_POP2_OIII5007A_Continuous_ND.dat for reading\nAborting\n");
    return -1;
  }
  for (i = 0; i < DIM_POP2_METALLICITY; i++) {
    fscanf(F, "%*[^\n]\n");
    fscanf(F, "%e %e %e %e %e %e %e", &o3_table_pop2[i][0], &o3_table_pop2[i][1], &o3_table_pop2[i][2], &o3_table_pop2[i][3],
                                      &o3_table_pop2[i][4], &o3_table_pop2[i][5], &o3_table_pop2[i][6]);
  }
  fclose(F);

  //initialize OII (3726) table
  if (!(F = fopen(O2S_POP2_TABLE_FILENAME, "r"))) {
    fprintf(stderr, "Unable to open file: LlineTbl_POP2_OII3726A_Continuous_ND.dat for reading\nAborting\n");
    return -1;
  }
  for (i = 0; i < DIM_POP2_METALLICITY; i++) {
    fscanf(F, "%*[^\n]\n");
    fscanf(F, "%e %e %e %e %e %e %e", &o2s_table_pop2[i][0], &o2s_table_pop2[i][1], &o2s_table_pop2[i][2], &o2s_table_pop2[i][3],
                                      &o2s_table_pop2[i][4], &o2s_table_pop2[i][5], &o2s_table_pop2[i][6]);
  }
  fclose(F);

  //initialize OII (3729) table
  if (!(F = fopen(O2L_POP2_TABLE_FILENAME, "r"))) {
    fprintf(stderr, "Unable to open file: LlineTbl_POP2_OII3729A_Continuous_ND.dat for reading\nAborting\n");
    return -1;
  }
  for (i = 0; i < DIM_POP2_METALLICITY; i++) {
    fscanf(F, "%*[^\n]\n");
    fscanf(F, "%e %e %e %e %e %e %e", &o2l_table_pop2[i][0], &o2l_table_pop2[i][1], &o2l_table_pop2[i][2], &o2l_table_pop2[i][3],
                                      &o2l_table_pop2[i][4], &o2l_table_pop2[i][5], &o2l_table_pop2[i][6]);
  }
  fclose(F);

  //initialize CII table
  if (!(F = fopen(C2_POP2_TABLE_FILENAME, "r"))) {
    fprintf(stderr, "Unable to open file: LlineTbl_POP2_CII157m.dat for reading\nAborting\n");
    return -1;
  }
  for (i = 0; i < DIM_POP2_METALLICITY; i++) {
    fscanf(F, "%*[^\n]\n");
    fscanf(F, "%e %e %e %e %e %e %e", &c2_table_pop2[i][0], &c2_table_pop2[i][1], &c2_table_pop2[i][2], &c2_table_pop2[i][3],
                                      &c2_table_pop2[i][4], &c2_table_pop2[i][5], &c2_table_pop2[i][6]);
  }
  fclose(F);

  //initialize CO(1-0) table
  if (!(F = fopen(CO10_POP2_TABLE_FILENAME, "r"))) {
    fprintf(stderr, "Unable to open file: LlineTbl_POP2_CO2600m.dat for reading\nAborting\n");
    return -1;
  }
  for (i = 0; i < DIM_POP2_METALLICITY; i++) {
    fscanf(F, "%*[^\n]\n");
    fscanf(F, "%e %e %e %e %e %e %e", &co10_table_pop2[i][0], &co10_table_pop2[i][1], &co10_table_pop2[i][2], &co10_table_pop2[i][3],
                                      &co10_table_pop2[i][4], &co10_table_pop2[i][5], &co10_table_pop2[i][6]);
  }
  fclose(F);

  // Populate Z axis
  for (ii = 0; ii < DIM_POP2_METALLICITY; ii++) {
    //rdata_Z[ii] = -3.+0.1*ii;
    rdata_Z[ii] = Z_VALUES[ii];
  }
  // Populate U axis
  for (jj = 0; jj < DIM_POP2_ION_PARAM; jj++) {
    rdata_U[jj] = -4.+0.5*jj;
  }
  // Populate flat luminosity array
  for (jj = 0; jj < DIM_POP2_ION_PARAM; jj++) {
    for (ii = 0; ii < DIM_POP2_METALLICITY; ii++) {
      int kk = ii + jj*DIM_POP2_METALLICITY;
      rdata_LUM_LYA_flat[kk] = lya_table_pop2[ii][jj];
      rdata_LUM_HA_flat[kk] = ha_table_pop2[ii][jj];
      rdata_LUM_HB_flat[kk] = hb_table_pop2[ii][jj];
      rdata_LUM_HE2_flat[kk] = he2_table_pop2[ii][jj];
      rdata_LUM_O3_flat[kk] = o3_table_pop2[ii][jj];
      rdata_LUM_O2S_flat[kk] = o2s_table_pop2[ii][jj];
      rdata_LUM_O2L_flat[kk] = o2l_table_pop2[ii][jj];
      rdata_LUM_C2_flat[kk] = c2_table_pop2[ii][jj];
      rdata_LUM_CO10_flat[kk] = co10_table_pop2[ii][jj];
    }
  }

  const gsl_interp2d_type *T = gsl_interp2d_bilinear;

  size_t nx = sizeof(rdata_Z) / sizeof(rdata_Z[0]);
  size_t ny = sizeof(rdata_U) / sizeof(rdata_U[0]);

  spline_lya_lum = gsl_spline2d_alloc(T, nx, ny);
  spline_ha_lum = gsl_spline2d_alloc(T, nx, ny);
  spline_hb_lum = gsl_spline2d_alloc(T, nx, ny);
  spline_he2_lum = gsl_spline2d_alloc(T, nx, ny);
  spline_o3_lum = gsl_spline2d_alloc(T, nx, ny);
  spline_o2s_lum = gsl_spline2d_alloc(T, nx, ny);
  spline_o2l_lum = gsl_spline2d_alloc(T, nx, ny);
  spline_c2_lum = gsl_spline2d_alloc(T, nx, ny);
  spline_co10_lum = gsl_spline2d_alloc(T, nx, ny);
  xacc_lum = gsl_interp_accel_alloc();
  yacc_lum = gsl_interp_accel_alloc();

  gsl_spline2d_init(spline_lya_lum, rdata_Z, rdata_U, rdata_LUM_LYA_flat, nx, ny);
  gsl_spline2d_init(spline_ha_lum, rdata_Z, rdata_U, rdata_LUM_HA_flat, nx, ny);
  gsl_spline2d_init(spline_hb_lum, rdata_Z, rdata_U, rdata_LUM_HB_flat, nx, ny);
  gsl_spline2d_init(spline_he2_lum, rdata_Z, rdata_U, rdata_LUM_HE2_flat, nx, ny);
  gsl_spline2d_init(spline_o3_lum, rdata_Z, rdata_U, rdata_LUM_O3_flat, nx, ny);
  gsl_spline2d_init(spline_o2s_lum, rdata_Z, rdata_U, rdata_LUM_O2S_flat, nx, ny);
  gsl_spline2d_init(spline_o2l_lum, rdata_Z, rdata_U, rdata_LUM_O2L_flat, nx, ny);
  gsl_spline2d_init(spline_c2_lum, rdata_Z, rdata_U, rdata_LUM_C2_flat, nx, ny);
  gsl_spline2d_init(spline_co10_lum, rdata_Z, rdata_U, rdata_LUM_CO10_flat, nx, ny);

  return 0;
}

float get_lya_luminosity(float Z, int logU) {
  return gsl_spline2d_eval(spline_lya_lum, Z, logU, xacc_lum, yacc_lum);
}

float get_ha_luminosity(float Z, int logU) {
  return gsl_spline2d_eval(spline_ha_lum, Z, logU, xacc_lum, yacc_lum);
}

float get_hb_luminosity(float Z, int logU) {
  return gsl_spline2d_eval(spline_hb_lum, Z, logU, xacc_lum, yacc_lum);
}

float get_he2_luminosity(float Z, int logU) {
  return gsl_spline2d_eval(spline_he2_lum, Z, logU, xacc_lum, yacc_lum);
}

float get_o3_luminosity(float Z, int logU) {
  return gsl_spline2d_eval(spline_o3_lum, Z, logU, xacc_lum, yacc_lum);
}

float get_o2s_luminosity(float Z, int logU) {
  return gsl_spline2d_eval(spline_o2s_lum, Z, logU, xacc_lum, yacc_lum);
}

float get_o2l_luminosity(float Z, int logU) {
  return gsl_spline2d_eval(spline_o2l_lum, Z, logU, xacc_lum, yacc_lum);
}

float get_c2_luminosity(float Z, int logU) {
  return gsl_spline2d_eval(spline_c2_lum, Z, logU, xacc_lum, yacc_lum);
}

float get_co10_luminosity(float Z, int logU) {
  return gsl_spline2d_eval(spline_co10_lum, Z, logU, xacc_lum, yacc_lum);
}


/*int main(int argc, char const *argv[])
{
  float lya_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM];
  float ha_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM];
  float hb_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM];
  float he2_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM];
  float o3_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM];
  float o2s_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM];
  float o2l_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM];
  float c2_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM];
  float co10_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM];

  init_pop2_full_tables(lya_table_pop2, ha_table_pop2, hb_table_pop2, he2_table_pop2,
                        o3_table_pop2, o2s_table_pop2, o2l_table_pop2, c2_table_pop2, co10_table_pop2);


  printf("CHECK Lya true value: %f, %f, %f\n", rdata_Z[0], rdata_U[2], lya_table_pop2[0][2]);
  printf("CHECK Lya interp value: %f\n", get_lya_luminosity(0.001, -3.));

  return 0;
}*/

