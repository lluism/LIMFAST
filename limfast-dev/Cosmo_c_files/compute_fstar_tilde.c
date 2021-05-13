#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#define FILE_TO_READf             "../External_tables/IntegratedSFE_ZSTART_30.0.hdf5"
#define DATASET_REDSHIFTf         "tab_z" // shape: 1501
#define DATASET_MASSf             "tab_M" // shape: 701
#define DATASET_MARf              "tab_f" // shape: 1501 x 701
#define DIM_REDSHIFTf             1501 // 
#define DIM_MASSf                 701
#define DIM_MARf                  1052201




static double rdata_zf[DIM_REDSHIFTf], rdata_mf[DIM_MASSf], rdata_marf[DIM_REDSHIFTf][DIM_MASSf], rdata_mar_flatf[DIM_MARf];

static gsl_spline2d *spline_ft;
static gsl_interp_accel *xacc_ft;
static gsl_interp_accel *yacc_ft;

float get_f_base(float redshift, float halo_mass)
{
  hid_t         file, space, dset_zf, dset_mf, dset_marf, dcpl;
  herr_t        status_zf, status_mf, status_marf;
  H5D_layout_t  layout;
  int           i, j;
  const gsl_interp2d_type *T = gsl_interp2d_bilinear;

  file = H5Fopen (FILE_TO_READf, H5F_ACC_RDONLY, H5P_DEFAULT);
  dset_zf = H5Dopen (file, DATASET_REDSHIFTf, H5P_DEFAULT);
  dset_mf = H5Dopen (file, DATASET_MASSf, H5P_DEFAULT);
  dset_marf = H5Dopen (file, DATASET_MARf, H5P_DEFAULT);

  /*
   * Retrieve the dataset creation property list, and print the storage layout
   */
  //dcpl = H5Dget_create_plist (dset);
  //layout = H5Pget_layout (dcpl);
  //printf("Storage layout for %s is: ", DATASET_TO_READ);

  status_zf = H5Dread (dset_zf, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata_zf);
  status_mf = H5Dread (dset_mf, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata_mf);
  status_marf = H5Dread (dset_marf, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata_marf);

  for (j=0; j<DIM_MASSf; j++)
    {
    for (i=0; i<DIM_REDSHIFTf; i++)
      {
      int k = i + j*DIM_REDSHIFTf;
      rdata_mar_flatf[k] = rdata_marf[i][j];
      }
    }

  //printf("%s:\n", DATASET_TO_READ);
  //printf("Entries of z array: %.3f, %.3f, %.3f\n", rdata_z[200], rdata_z[200], rdata_z[200]);
  //printf("Entries of m array: %.3f, %.3f, %.3f\n", rdata_m[500], rdata_m[501], rdata_m[502]);
  //printf("Entries of mar array: %.3f, %.3f, %.3f\n", rdata_mar[200][500], rdata_mar[200][501], rdata_mar[200][502]);

  size_t nx = sizeof(rdata_zf) / sizeof(rdata_zf[0]);
  size_t ny = sizeof(rdata_mf) / sizeof(rdata_mf[0]);

  gsl_spline2d *spline_ft = gsl_spline2d_alloc(T, nx, ny);
  gsl_interp_accel *xacc_ft = gsl_interp_accel_alloc();
  gsl_interp_accel *yacc_ft = gsl_interp_accel_alloc();

  gsl_spline2d_init(spline_ft, rdata_zf, rdata_mf, rdata_mar_flatf, nx, ny);

  //status = H5Pclose (dcpl);
  status_zf = H5Fclose (file);
  status_zf = H5Dclose (dset_zf);
  status_mf = H5Dclose (dset_mf);
  status_marf = H5Dclose (dset_marf);

  return gsl_spline2d_eval(spline_ft, redshift, halo_mass, xacc_ft, yacc_ft);
}


void init_f_tbl()
{
  hid_t         file, space, dset_zf, dset_mf, dset_marf, dcpl;
  herr_t        status_zf, status_mf, status_marf;
  H5D_layout_t  layout;
  int           i, j;

  file = H5Fopen (FILE_TO_READf, H5F_ACC_RDONLY, H5P_DEFAULT);
  dset_zf = H5Dopen (file, DATASET_REDSHIFTf, H5P_DEFAULT);
  dset_mf = H5Dopen (file, DATASET_MASSf, H5P_DEFAULT);
  dset_marf = H5Dopen (file, DATASET_MARf, H5P_DEFAULT);

  /*
   * Retrieve the dataset creation property list, and print the storage layout
   */
  //dcpl = H5Dget_create_plist (dset);
  //layout = H5Pget_layout (dcpl);
  //printf("Storage layout for %s is: ", DATASET_TO_READ);

  status_zf = H5Dread (dset_zf, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata_zf);
  status_mf = H5Dread (dset_mf, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata_mf);
  status_marf = H5Dread (dset_marf, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata_marf);

  for (j=0; j<DIM_MASSf; j++)
    {
    for (i=0; i<DIM_REDSHIFTf; i++)
      {
      int k = i + j*DIM_REDSHIFTf;
      rdata_mar_flatf[k] = rdata_marf[i][j];
      }
    }

  //printf("%s:\n", DATASET_TO_READ);
  //printf("Entries of z array: %.3f, %.3f, %.3f\n", rdata_z[200], rdata_z[200], rdata_z[200]);
  //printf("Entries of m array: %.3f, %.3f, %.3f\n", rdata_m[500], rdata_m[501], rdata_m[502]);
  //printf("Entries of mar array: %.3f, %.3f, %.3f\n", rdata_mar[200][500], rdata_mar[200][501], rdata_mar[200][502]);

  //status = H5Pclose (dcpl);
  status_zf = H5Fclose (file);
  status_zf = H5Dclose (dset_zf);
  status_mf = H5Dclose (dset_mf);
  status_marf = H5Dclose (dset_marf);

  const gsl_interp2d_type *T = gsl_interp2d_bilinear;

  size_t nx = sizeof(rdata_zf) / sizeof(rdata_zf[0]);
  size_t ny = sizeof(rdata_mf) / sizeof(rdata_mf[0]);

  spline_ft = gsl_spline2d_alloc(T, nx, ny);
  xacc_ft = gsl_interp_accel_alloc();
  yacc_ft = gsl_interp_accel_alloc();

  gsl_spline2d_init(spline_ft, rdata_zf, rdata_mf, rdata_mar_flatf, nx, ny);
}



float get_ftilde(float redshift, float halo_mass)
{
  return gsl_spline2d_eval(spline_ft, redshift, halo_mass, xacc_ft, yacc_ft);
}


//
//int main(int argc, char const *argv[])
//{
//  init_mar_tbl();
//  printf("CHECK: %e\n", get_mar_base(10.0, 1.520e9));
//  printf("CHECK: %e\n", get_mar(10.0, 1.520e9));
//  return 0;
//}
