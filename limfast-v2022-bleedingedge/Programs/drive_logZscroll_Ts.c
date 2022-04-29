#include <math.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <sys/wait.h>

#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"
#include "../Parameter_files/HEAT_PARAMS.H"
#include "../Parameter_files/SOURCES.H"

/*
  Program DRIVE_ZSCROLL.C scrolls through the redshifts defined in ANAL_PARAMS.H creating halo, velocity, density, and ionization fields
*/

// IMPORTANT NOTE: because of the interpolation range of the bathtub model, ZLOW must not be smaller than 5.0 whereas ZHIGH must not be
// greater than 25.0. -GS
#define ZLOW (float) (14.0)
//#define ZHIGH  Z_HEAT_MAX
#define ZHIGH (float) (15.0)

int main(int argc, char ** argv){
  //float Z, M, M_MIN, nf;
  float M_MIN;
  float Z, M, nf;
  char cmnd[1000];
  FILE *LOG;
  time_t start_time, curr_time;
  int status;


  time(&start_time);

  // make appropriate directories
  system("mkdir ../Log_files");
  system("mkdir ../Boxes");
  system("mkdir ../Output_files");
  system("mkdir ../Output_files/DNDLNM_files");
  system("mkdir ../Output_files/FgtrM_files");
  system("mkdir ../Output_files/Halo_lists");
  system("mkdir ../Output_files/Size_distributions");
  system("mkdir ../Output_files/Deldel_T_power_spec");
  system("mkdir ../Redshift_interpolate_filelists");
  //  system("mkdir ../Lighttravel_filelists");

  // remove some of the previous (astro) files which might conflict with current run
  system("rm ../Boxes/Ts_evolution/*");
  system("rm ../Boxes/Ts_*");
  system("rm ../Boxes/delta_T_*");
  system("rm ../Boxes/xH_*");
  system("rm ../Boxes/Nrec_*");
  system("rm ../Boxes/z_first*");
  system("rm ../Output_files/Deldel_T_power_spec/*");

  init_ps();

  // open log file
  system("rm ../Log_files/*");
  LOG = log_open("../Log_files/drive_logzscroll_log_file");
  if (!LOG){
    fprintf(stderr, "drive_zscroll_log_file.c: Unable to open log file\n Aborting...\n");
    return -1;
  }


  fprintf(stderr, "Calling init to set up the initial conditions\n");
  fprintf(LOG, "Calling init to set up the initial conditions\n");
  system("./init"); // you only need this call once per realization

  Z = ZLOW*1.0001; // match rounding convention from Ts.c

   // call Ts on the lowest redshift
  if (USE_TS_IN_21CM){
    sprintf(cmnd, "./perturb_field %.2f", Z);
    time(&curr_time);
    fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
    fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
    fflush(NULL);
    system(cmnd);

    sprintf(cmnd, "./Ts %.2f", Z);
    time(&curr_time);
    fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
    fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
    fflush(NULL);
    system(cmnd);
  }


  // now go to highest redshift and step downwards
  while (Z < ZHIGH){
    Z = ((1+Z)*ZPRIME_STEP_FACTOR - 1);
  }
  Z = ((1+Z)/ ZPRIME_STEP_FACTOR - 1);
  while (Z >= ZLOW){

    //set the minimum source mass
    //M_MIN = get_M_min_ion(Z);

    if(USE_GENERAL_SOURCES)
    {
      sources src;
      src = defaultSources();
      M_MIN = src.minMass(Z);
      printf("At z = %f, M_MIN = %e MSUN\n", Z, M_MIN);
    }

    // if USE_HALO_FIELD is turned on in ANAL_PARAMS.H, run the halo finder
    if (USE_HALO_FIELD){
      //  the following only depend on redshift, not ionization field
      // find halos
      sprintf(cmnd, "./find_halos %.2f", Z);
      time(&curr_time);
      fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
      fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
      fflush(NULL);
      system(cmnd);


      // shift halos accordig to their linear velocities
      sprintf(cmnd, "./update_halo_pos %.2f", Z);
      time(&curr_time);
      fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
      fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
      fflush(NULL);
      system(cmnd);
    }

    // shift density field and update velocity field
    sprintf(cmnd, "./perturb_field %.2f", Z);
    time(&curr_time);
    fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
    fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
    fflush(NULL);
    system(cmnd);
    // end of solely redshift dependent things, now do ionization stuff


    // if it is the lowest redshift, let's call Ts.c
    if (USE_TS_IN_21CM && (Z > Z_HEAT_MAX) ) { // NEW CONDITIONAL
      //    if (USE_TS_IN_21CM && (fabs(Z-ZLOW)/Z < 0.0002) ){
      sprintf(cmnd, "./Ts %.2f", Z);
      time(&curr_time);
      fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
      fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
      fflush(NULL);
      system(cmnd);
    } // this will create all of the higher z Ts files in Boxes, provided Ts_verbose is turned on
    // in HEAT_PARAMS.H


    // find bubbles
    if (INHOMO_RECO)
      sprintf(cmnd, "./find_HII_bubbles %f %f", Z, (1+Z)*ZPRIME_STEP_FACTOR - 1 );
    else
      sprintf(cmnd, "./find_HII_bubbles %f", Z );
    time(&curr_time);
    fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
    fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
    fflush(NULL);
    status = system(cmnd);
    nf = WEXITSTATUS(status) / 100.0;
    if (nf < 0){
      fprintf(stderr, "find_HII_bubbles exited...\nAborting run...\n");
      fprintf(LOG,  "find_HII_bubbles exited...\nAborting run...\n");
      return -1;
    }

    /*
    // generate size distributions, first ionized bubbles
    switch(FIND_BUBBLE_ALGORITHM){
    case 2:
      if (USE_HALO_FIELD)
        sprintf(cmnd, "./gen_size_distr %06.2f 0 ../Boxes/xH_z%06.2f_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", Z,Z, HII_EFF_FACTOR, EFF_FACTOR_PL_INDEX, HII_FILTER, M_MIN, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
      else
        sprintf(cmnd, "./gen_size_distr %06.2f 0 ../Boxes/xH_z%06.2f_nohalos_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", Z,Z, HII_EFF_FACTOR, EFF_FACTOR_PL_INDEX, HII_FILTER, M_MIN, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
      break;
    default:
      if (USE_HALO_FIELD)
        sprintf(cmnd, "./gen_size_distr %06.2f 0 ../Boxes/sphere_xH_z%06.2f_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", Z, Z,HII_EFF_FACTOR, EFF_FACTOR_PL_INDEX, HII_FILTER, M_MIN, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
      else
        sprintf(cmnd, "./gen_size_distr %06.2f 0 ../Boxes/sphere_xH_nohalos_z%06.2f_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", Z, Z,HII_EFF_FACTOR, EFF_FACTOR_PL_INDEX, HII_FILTER, M_MIN, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
    }
    time(&curr_time);
    fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
    fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
    fflush(NULL);
    system(cmnd);


    // generate size distributions, then neutral regions
    switch(FIND_BUBBLE_ALGORITHM){
    case 2:
      if (USE_HALO_FIELD)
        sprintf(cmnd, "./gen_size_distr %06.2f 1 ../Boxes/xH_z%06.2f_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", Z, Z,HII_EFF_FACTOR, EFF_FACTOR_PL_INDEX, HII_FILTER, M_MIN, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
      else
        sprintf(cmnd, "./gen_size_distr %06.2f 1 ../Boxes/xH_nohalos_z%06.2f_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", Z, Z,HII_EFF_FACTOR, EFF_FACTOR_PL_INDEX, HII_FILTER, M_MIN, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
      break;
    default:
      if (USE_HALO_FIELD)
        sprintf(cmnd, "./gen_size_distr %06.2f 1 ../Boxes/sphere_xH_z%06.2f_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", Z, Z,HII_EFF_FACTOR, EFF_FACTOR_PL_INDEX, HII_FILTER, M_MIN, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
      else
        sprintf(cmnd, "./gen_size_distr %06.2f 1 ../Boxes/sphere_xH_nohalos_z%06.2f_nf*_eff%.1f_effPLindex%.1f_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", Z, Z,HII_EFF_FACTOR, EFF_FACTOR_PL_INDEX, HII_FILTER, M_MIN, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
    }
    time(&curr_time);
    fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
    fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
    fflush(NULL);
    system(cmnd);
    */


    // do temperature map
    switch(FIND_BUBBLE_ALGORITHM){
    case 2:
      if (USE_HALO_FIELD)
	sprintf(cmnd, "./delta_T %06.2f ../Boxes/xH_z%06.2f_nf*_%i_%.0fMpc ../Boxes/Ts_z%06.2f_*_%.0fMpc", Z, Z, HII_DIM, BOX_LEN, Z, BOX_LEN);
      else
	sprintf(cmnd, "./delta_T %06.2f ../Boxes/xH_nohalos_z%06.2f_nf*_%i_%.0fMpc ../Boxes/Ts_z%06.2f_*_%.0fMpc", Z, Z, HII_DIM, BOX_LEN, Z, BOX_LEN);
      break;
    default:
      if (USE_HALO_FIELD)
	sprintf(cmnd, "./delta_T %06.2f ../Boxes/sphere_xH_z%06.2f_nf*_%i_%.0fMpc ../Boxes/Ts_z%06.2f_*_%.0fMpc", Z, Z, HII_DIM, BOX_LEN, Z, BOX_LEN);
      else
	sprintf(cmnd, "./delta_T %06.2f ../Boxes/sphere_xH_nohalos_z%06.2f_nf*_%i_%.0fMpc ../Boxes/Ts_z%06.2f_*_%.0fMpc", Z, Z, HII_DIM, BOX_LEN, Z, BOX_LEN);
      break;
    }
    time(&curr_time);
    fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
    fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
    fflush(NULL);
    system(cmnd);

    // find bubbles
    if (RSD_FOR_LIM) {
      // run RSD for CII158m
      sprintf(cmnd, "./RSD %06.2f ../Boxes/I_CII158m_z%06.2f_*_%i_%.0fMpc CII158m", Z, Z, HII_DIM, BOX_LEN);
      time(&curr_time);
      fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
      fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
      fflush(NULL);
      system(cmnd);
      // run RSD for CO2600m
      sprintf(cmnd, "./RSD %06.2f ../Boxes/I_CO2600m_z%06.2f_*_%i_%.0fMpc CO2600m", Z, Z, HII_DIM, BOX_LEN);
      time(&curr_time);
      fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
      fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
      fflush(NULL);
      system(cmnd);
      // run RSD for OIII88m
      sprintf(cmnd, "./RSD %06.2f ../Boxes/I_OIII88m_z%06.2f_*_%i_%.0fMpc OIII88m", Z, Z, HII_DIM, BOX_LEN);
      time(&curr_time);
      fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
      fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
      fflush(NULL);
      system(cmnd);
      // run RSD for HI6563A
      sprintf(cmnd, "./RSD %06.2f ../Boxes/I_HI6563A_z%06.2f_*_%i_%.0fMpc HI6563A", Z, Z, HII_DIM, BOX_LEN);
      time(&curr_time);
      fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
      fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
      fflush(NULL);
      system(cmnd);
      // run RSD for HI4861A
      sprintf(cmnd, "./RSD %06.2f ../Boxes/I_HI4861A_z%06.2f_*_%i_%.0fMpc HI4861A", Z, Z, HII_DIM, BOX_LEN);
      time(&curr_time);
      fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
      fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
      fflush(NULL);
      system(cmnd);
      // run RSD for HI1216A
      sprintf(cmnd, "./RSD %06.2f ../Boxes/I_HI1216A_z%06.2f_*_%i_%.0fMpc HI1216A", Z, Z, HII_DIM, BOX_LEN);
      time(&curr_time);
      fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
      fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
      fflush(NULL);
      system(cmnd);
      // run RSD for HI1216AdIGM
      sprintf(cmnd, "./RSD %06.2f ../Boxes/I_HI1216A_z%06.2f_*_%i_%.0fMpc HI1216AdIGM", Z, Z, HII_DIM, BOX_LEN);
      time(&curr_time);
      fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
      fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, -difftime(start_time, curr_time)/60.0);
      fflush(NULL);
      system(cmnd);
    }

    fprintf(stderr, "*************************************\n");
    fflush(NULL);

    // update the redshift value according to the logarithmic stepping in the Ts.c routine
    Z = ((1+Z)/ZPRIME_STEP_FACTOR - 1);
  }


  // Create lightcone boxes from the coeval cubes
  sprintf(cmnd, "ls ../Boxes/xH_*%i_%.0fMpc > ../Redshift_interpolate_filelists/xH_%i_%.0fMpc",
	  HII_DIM, BOX_LEN, HII_DIM, BOX_LEN);
  system(cmnd);
  sprintf(cmnd, "redshift_interpolate_boxes 0 ../Redshift_interpolate_filelists/xH_%i_%.0fMpc", HII_DIM, BOX_LEN);
  system(cmnd);
  fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
  fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
  fflush(NULL);

  sprintf(cmnd, "ls ../Boxes/delta_T_*%i_%.0fMpc > ../Redshift_interpolate_filelists/delta_T_%i_%.0fMpc",
	  HII_DIM, BOX_LEN, HII_DIM, BOX_LEN);
  system(cmnd);
  sprintf(cmnd, "redshift_interpolate_boxes 0 ../Redshift_interpolate_filelists/delta_T_%i_%.0fMpc", HII_DIM, BOX_LEN);
  system(cmnd);
  fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
  fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
  fflush(NULL);

  if (INHOMO_RECO){
    sprintf(cmnd, "ls ../Boxes/Nrec_*%i_%.0fMpc > ../Redshift_interpolate_filelists/Nrec_%i_%.0fMpc",
	    HII_DIM, BOX_LEN, HII_DIM, BOX_LEN);
    system(cmnd);
    sprintf(cmnd, "redshift_interpolate_boxes 0 ../Redshift_interpolate_filelists/Nrec_%i_%.0fMpc", HII_DIM, BOX_LEN);
    system(cmnd);
    fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
    fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
    fflush(NULL);
  }

  sprintf(cmnd, "./extract_delTps.pl 0.1 ../Output_files/Deldel_T_power_spec/ps_z0* > ../Output_files/Deldel_T_power_spec/Power_k0.1vsRedshift.txt");
  system(cmnd);
  fprintf(stderr, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
  fprintf(LOG, "Now calling: %s, %g min have ellapsed\n", cmnd, difftime(curr_time, start_time)/60.0);
  fflush(NULL);

  fclose(LOG);
  return 0;
}
