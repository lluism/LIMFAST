#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#define FILE_TO_READ             "../External_tables/hmf/hmf_ST_planck_TTTEEE_lowl_lowE_best_logM_1400_4-18_z_1201_0-60.hdf5"
#define DATASET_REDSHIFT         "tab_z" // shape: 1203
#define DATASET_MASS             "tab_M" // shape: 1400
#define DATASET_MAR              "tab_MAR" // shape: 1203 x 1400
#define DIM_REDSHIFT             1203 // Redshift really goes from 0 to 60.1, so 60.1/0.05 + 1 = 1203 -GS
#define DIM_MASS                 1400
#define DIM_MAR                  1684200

static double rdata_z[DIM_REDSHIFT], rdata_m[DIM_MASS], rdata_mar[DIM_REDSHIFT][DIM_MASS], rdata_mar_flat[DIM_MAR];

static gsl_spline2d *spline_mar;
static gsl_interp_accel *xacc_mar;
static gsl_interp_accel *yacc_mar;

float get_mar_base(float redshift, float halo_mass)
{
  hid_t         file, space, dset_z, dset_m, dset_mar, dcpl;
  herr_t        status_z, status_m, status_mar;
  H5D_layout_t  layout;
  int           i, j;
  const gsl_interp2d_type *T = gsl_interp2d_bilinear;

  file = H5Fopen (FILE_TO_READ, H5F_ACC_RDONLY, H5P_DEFAULT);
  dset_z = H5Dopen (file, DATASET_REDSHIFT, H5P_DEFAULT);
  dset_m = H5Dopen (file, DATASET_MASS, H5P_DEFAULT);
  dset_mar = H5Dopen (file, DATASET_MAR, H5P_DEFAULT);

  /*
   * Retrieve the dataset creation property list, and print the storage layout
   */
  //dcpl = H5Dget_create_plist (dset);
  //layout = H5Pget_layout (dcpl);
  //printf("Storage layout for %s is: ", DATASET_TO_READ);

  status_z = H5Dread (dset_z, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata_z);
  status_m = H5Dread (dset_m, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata_m);
  status_mar = H5Dread (dset_mar, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata_mar);

  for (j=0; j<DIM_MASS; j++)
    {
    for (i=0; i<DIM_REDSHIFT; i++)
      {
      int k = i + j*DIM_REDSHIFT;
      rdata_mar_flat[k] = rdata_mar[i][j];
      }
    }

  //printf("%s:\n", DATASET_TO_READ);
  //printf("Entries of z array: %.3f, %.3f, %.3f\n", rdata_z[200], rdata_z[200], rdata_z[200]);
  //printf("Entries of m array: %.3f, %.3f, %.3f\n", rdata_m[500], rdata_m[501], rdata_m[502]);
  //printf("Entries of mar array: %.3f, %.3f, %.3f\n", rdata_mar[200][500], rdata_mar[200][501], rdata_mar[200][502]);

  size_t nx = sizeof(rdata_z) / sizeof(rdata_z[0]);
  size_t ny = sizeof(rdata_m) / sizeof(rdata_m[0]);

  gsl_spline2d *spline_mar = gsl_spline2d_alloc(T, nx, ny);
  gsl_interp_accel *xacc_mar = gsl_interp_accel_alloc();
  gsl_interp_accel *yacc_mar = gsl_interp_accel_alloc();

  gsl_spline2d_init(spline_mar, rdata_z, rdata_m, rdata_mar_flat, nx, ny);

  //status = H5Pclose (dcpl);
  status_z = H5Fclose (file);
  status_z = H5Dclose (dset_z);
  status_m = H5Dclose (dset_m);
  status_mar = H5Dclose (dset_mar);

  return gsl_spline2d_eval(spline_mar, redshift, halo_mass, xacc_mar, yacc_mar);
}


void init_mar_tbl()
{
  hid_t         file, space, dset_z, dset_m, dset_mar, dcpl;
  herr_t        status_z, status_m, status_mar;
  H5D_layout_t  layout;
  int           i, j;

  file = H5Fopen (FILE_TO_READ, H5F_ACC_RDONLY, H5P_DEFAULT);
  dset_z = H5Dopen (file, DATASET_REDSHIFT, H5P_DEFAULT);
  dset_m = H5Dopen (file, DATASET_MASS, H5P_DEFAULT);
  dset_mar = H5Dopen (file, DATASET_MAR, H5P_DEFAULT);

  /*
   * Retrieve the dataset creation property list, and print the storage layout
   */
  //dcpl = H5Dget_create_plist (dset);
  //layout = H5Pget_layout (dcpl);
  //printf("Storage layout for %s is: ", DATASET_TO_READ);

  status_z = H5Dread (dset_z, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata_z);
  status_m = H5Dread (dset_m, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata_m);
  status_mar = H5Dread (dset_mar, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata_mar);

  for (j=0; j<DIM_MASS; j++)
    {
    for (i=0; i<DIM_REDSHIFT; i++)
      {
      int k = i + j*DIM_REDSHIFT;
      rdata_mar_flat[k] = rdata_mar[i][j];
      }
    }

  //printf("%s:\n", DATASET_TO_READ);
  //printf("Entries of z array: %.3f, %.3f, %.3f\n", rdata_z[200], rdata_z[200], rdata_z[200]);
  //printf("Entries of m array: %.3f, %.3f, %.3f\n", rdata_m[500], rdata_m[501], rdata_m[502]);
  //printf("Entries of mar array: %.3f, %.3f, %.3f\n", rdata_mar[200][500], rdata_mar[200][501], rdata_mar[200][502]);

  //status = H5Pclose (dcpl);
  status_z = H5Fclose (file);
  status_z = H5Dclose (dset_z);
  status_m = H5Dclose (dset_m);
  status_mar = H5Dclose (dset_mar);

  const gsl_interp2d_type *T = gsl_interp2d_bilinear;

  size_t nx = sizeof(rdata_z) / sizeof(rdata_z[0]);
  size_t ny = sizeof(rdata_m) / sizeof(rdata_m[0]);

  spline_mar = gsl_spline2d_alloc(T, nx, ny);
  xacc_mar = gsl_interp_accel_alloc();
  yacc_mar = gsl_interp_accel_alloc();

  gsl_spline2d_init(spline_mar, rdata_z, rdata_m, rdata_mar_flat, nx, ny);
}



float get_mar(float redshift, float halo_mass)
{
  return gsl_spline2d_eval(spline_mar, redshift, halo_mass, xacc_mar, yacc_mar);
}


//
//int main(int argc, char const *argv[])
//{
//  init_mar_tbl();
//  printf("CHECK: %e\n", get_mar_base(10.0, 1.520e9));
//  printf("CHECK: %e\n", get_mar(10.0, 1.520e9));
//  return 0;
//}
