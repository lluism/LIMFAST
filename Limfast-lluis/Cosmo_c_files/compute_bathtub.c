#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#define FILE_TO_READ_BT             "../External_tables/bathtub/tbl_z_251_5-30_logM_901_7_16_momentum.hdf5"
#define DATASET_REDSHIFT_BT         "z"       // redshift, shape: 251
#define DATASET_HALOMASS_BT         "mh"      // halo mass, shape: 801
#define DATASET_SFR_BT              "sfr"     // SFR, shape: 801 x 251
#define DATASET_Z_BT                "Z"       // metallicity, shape: 801 x 251
#define DATASET_MG_BT               "Mg"      // gas mass, shape: 801 x 251
#define DATASET_MS_BT               "Ms"      // stellar mass, shape: 801 x 251
#define DATASET_FSTAR_BT            "fstar"   // stellar mass, shape: 801 x 251
#define DIM_REDSHIFT_BT             251
#define DIM_HALOMASS_BT             901
#define DIM_SFR_BT                  226151
#define DIM_Z_BT                    226151
#define DIM_MG_BT                   226151
#define DIM_MS_BT                   226151
#define DIM_FSTAR_BT                226151

static double rdata_z_bt[DIM_REDSHIFT_BT], rdata_m_bt[DIM_HALOMASS_BT],
              rdata_sfr_bt[DIM_REDSHIFT_BT][DIM_HALOMASS_BT], rdata_sfr_flat_bt[DIM_SFR_BT],
              rdata_Z_bt[DIM_REDSHIFT_BT][DIM_HALOMASS_BT], rdata_Z_flat_bt[DIM_Z_BT],
              rdata_Mg_bt[DIM_REDSHIFT_BT][DIM_HALOMASS_BT], rdata_Mg_flat_bt[DIM_MG_BT],
              rdata_Ms_bt[DIM_REDSHIFT_BT][DIM_HALOMASS_BT], rdata_Ms_flat_bt[DIM_MS_BT],
              rdata_fstar_bt[DIM_REDSHIFT_BT][DIM_HALOMASS_BT], rdata_fstar_flat_bt[DIM_FSTAR_BT];

// For SFR
static gsl_spline2d *spline_sfr_bt;
static gsl_interp_accel *xacc_sfr_bt;
static gsl_interp_accel *yacc_sfr_bt;
// For metallicity
static gsl_spline2d *spline_Z_bt;
static gsl_interp_accel *xacc_Z_bt;
static gsl_interp_accel *yacc_Z_bt;
// For gas mass
static gsl_spline2d *spline_Mg_bt;
static gsl_interp_accel *xacc_Mg_bt;
static gsl_interp_accel *yacc_Mg_bt;
// For stellar mass
static gsl_spline2d *spline_Ms_bt;
static gsl_interp_accel *xacc_Ms_bt;
static gsl_interp_accel *yacc_Ms_bt;
// For fstar
static gsl_spline2d *spline_fstar_bt;
static gsl_interp_accel *xacc_fstar_bt;
static gsl_interp_accel *yacc_fstar_bt;

void init_galaxymodel_tbl()
{
  hid_t         file, space, dset_z, dset_m, dset_sfr, dset_Z, dset_Mg, dset_Ms, dset_fstar, dcpl;
  herr_t        status_z, status_m, status_sfr, status_Z, status_Mg, status_Ms, status_fstar;
  H5D_layout_t  layout;
  int           i, j;

  file = H5Fopen (FILE_TO_READ_BT, H5F_ACC_RDONLY, H5P_DEFAULT);
  dset_z = H5Dopen (file, DATASET_REDSHIFT_BT, H5P_DEFAULT);
  dset_m = H5Dopen (file, DATASET_HALOMASS_BT, H5P_DEFAULT);
  dset_sfr = H5Dopen (file, DATASET_SFR_BT, H5P_DEFAULT);
  dset_Z = H5Dopen (file, DATASET_Z_BT, H5P_DEFAULT);
  dset_Mg = H5Dopen (file, DATASET_MG_BT, H5P_DEFAULT);
  dset_Ms = H5Dopen (file, DATASET_MS_BT, H5P_DEFAULT);
  dset_fstar = H5Dopen (file, DATASET_FSTAR_BT, H5P_DEFAULT);

  /*
   * Retrieve the dataset creation property list, and print the storage layout
   */
  //dcpl = H5Dget_create_plist (dset);
  //layout = H5Pget_layout (dcpl);
  //printf("Storage layout for %s is: ", DATASET_TO_READ);

  status_z = H5Dread (dset_z, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata_z_bt);
  status_m = H5Dread (dset_m, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata_m_bt);
  status_sfr = H5Dread (dset_sfr, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata_sfr_bt);
  status_Z = H5Dread (dset_Z, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata_Z_bt);
  status_Mg = H5Dread (dset_Mg, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata_Mg_bt);
  status_Ms = H5Dread (dset_Ms, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata_Ms_bt);
  status_fstar = H5Dread (dset_fstar, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata_fstar_bt);

  for (j=0; j<DIM_HALOMASS_BT; j++)
    {
    for (i=0; i<DIM_REDSHIFT_BT; i++)
      {
      int k = i + j*DIM_REDSHIFT_BT;
      rdata_sfr_flat_bt[k] = rdata_sfr_bt[i][j];
      rdata_Z_flat_bt[k] = rdata_Z_bt[i][j];
      rdata_Mg_flat_bt[k] = rdata_Mg_bt[i][j];
      rdata_Ms_flat_bt[k] = rdata_Ms_bt[i][j];
      rdata_fstar_flat_bt[k] = rdata_fstar_bt[i][j];
      }
    }

  //status = H5Pclose (dcpl);
  status_z = H5Fclose (file);
  status_z = H5Dclose (dset_z);
  status_m = H5Dclose (dset_m);
  status_sfr = H5Dclose (dset_sfr);
  status_Z = H5Dclose (dset_Z);
  status_Mg = H5Dclose (dset_Mg);
  status_Ms = H5Dclose (dset_Ms);
  status_fstar = H5Dclose (dset_fstar);

  const gsl_interp2d_type *T = gsl_interp2d_bilinear;

  size_t nx = sizeof(rdata_z_bt) / sizeof(rdata_z_bt[0]);
  size_t ny = sizeof(rdata_m_bt) / sizeof(rdata_m_bt[0]);

  spline_sfr_bt = gsl_spline2d_alloc(T, nx, ny);
  xacc_sfr_bt = gsl_interp_accel_alloc();
  yacc_sfr_bt = gsl_interp_accel_alloc();
  spline_Z_bt = gsl_spline2d_alloc(T, nx, ny);
  xacc_Z_bt = gsl_interp_accel_alloc();
  yacc_Z_bt = gsl_interp_accel_alloc();
  spline_Mg_bt = gsl_spline2d_alloc(T, nx, ny);
  xacc_Mg_bt = gsl_interp_accel_alloc();
  yacc_Mg_bt = gsl_interp_accel_alloc();
  spline_Ms_bt = gsl_spline2d_alloc(T, nx, ny);
  xacc_Ms_bt = gsl_interp_accel_alloc();
  yacc_Ms_bt = gsl_interp_accel_alloc();
  spline_fstar_bt = gsl_spline2d_alloc(T, nx, ny);
  xacc_fstar_bt = gsl_interp_accel_alloc();
  yacc_fstar_bt = gsl_interp_accel_alloc();

  gsl_spline2d_init(spline_sfr_bt, rdata_z_bt, rdata_m_bt, rdata_sfr_flat_bt, nx, ny);
  gsl_spline2d_init(spline_Z_bt, rdata_z_bt, rdata_m_bt, rdata_Z_flat_bt, nx, ny);
  gsl_spline2d_init(spline_Mg_bt, rdata_z_bt, rdata_m_bt, rdata_Mg_flat_bt, nx, ny);
  gsl_spline2d_init(spline_Ms_bt, rdata_z_bt, rdata_m_bt, rdata_Ms_flat_bt, nx, ny);
  gsl_spline2d_init(spline_fstar_bt, rdata_z_bt, rdata_m_bt, rdata_fstar_flat_bt, nx, ny);
}



float get_sfr(float redshift, float halo_mass)
{
  return gsl_spline2d_eval(spline_sfr_bt, redshift, halo_mass, xacc_sfr_bt, yacc_sfr_bt);
}

float get_Z(float redshift, float halo_mass)
{
  return gsl_spline2d_eval(spline_Z_bt, redshift, halo_mass, xacc_Z_bt, yacc_Z_bt);
}

float get_Mg(float redshift, float halo_mass)
{
  return gsl_spline2d_eval(spline_Mg_bt, redshift, halo_mass, xacc_Mg_bt, yacc_Mg_bt);
}

float get_Ms(float redshift, float halo_mass)
{
  return gsl_spline2d_eval(spline_Ms_bt, redshift, halo_mass, xacc_Ms_bt, yacc_Ms_bt);
}

float get_fstar(float redshift, float halo_mass)
{
  return gsl_spline2d_eval(spline_fstar_bt, redshift, halo_mass, xacc_fstar_bt, yacc_fstar_bt);
}



//int main(int argc, char const *argv[])
//{
//  init_sfr_tbl();
//  printf("CHECK: %e\n", get_sfr(7.0, 1.0e12));
//  return 0;
//}
