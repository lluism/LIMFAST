#include "heating_helper_progs.c"

/*
  This is completed version.
*/

/*
  Program Ts calculates the spin temperature field, according to the perscription outlined in
  Mesinger, Furlanetto, Cen (2010).  The fluctuating component is sourced by the mean EPS collapsed fraction in spherical annyli surrounding each cell.

  Usage: Ts <REDSHIFT> [reload zp redshift]  [<stellar fraction for 10^10 Msun halos> <power law index for stellar fraction halo mass scaling>
       <escape fraction for 10^10 Msun halos> <power law index for escape fraction halo mass scaling>
	   <turn-over scale for the duty cycle of galaxies, in units of halo mass> <Soft band X-ray luminosity>]

  The last optional argument is the z' redshift output from which to reload intermediate
  evolution files in ../Boxes/Ts_evolution/

  Memory usage (in floats)~ (<NUMBER OF FILTER STEPS> + 3) x HII_DIM^3

  Author: Andrei Mesinger
  Date: 9.2.2009
*/


// New in v1.4: To calculate follapse fraction for new parametrization
void init_21cmMC_arrays() {

	int i,j;

    for (i=0; i < NUM_FILTER_STEPS_FOR_Ts; i++){
      FcollLowXray_zpp_spline_acc[i] = gsl_interp_accel_alloc ();
      FcollLowLya_zpp_spline_acc[i] = gsl_interp_accel_alloc ();
      FcollLowXray_zpp_spline[i] = gsl_spline_alloc (gsl_interp_cspline, NSFR_low);
      FcollLowLya_zpp_spline[i] = gsl_spline_alloc (gsl_interp_cspline, NSFR_low);

      second_derivs_FcollXray_zpp[i] = calloc(NSFR_high,sizeof(float));
      second_derivs_FcollLya_zpp[i] = calloc(NSFR_high,sizeof(float));
    }

    xi_SFR = calloc((NGL_SFR+1),sizeof(float));
    wi_SFR = calloc((NGL_SFR+1),sizeof(float));

	zpp_interp_table = calloc(zpp_interp_points, sizeof(float));

	redshift_interp_table = calloc(NUM_FILTER_STEPS_FOR_Ts*Nsteps_zp, sizeof(float)); // New

	log10_overdense_low_table = calloc(NSFR_low,sizeof(double));
	log10_FcollzXray_SFR_low_table = (double **)calloc(NUM_FILTER_STEPS_FOR_Ts*Nsteps_zp,sizeof(double *));
  log10_FcollzLya_SFR_low_table = (double **)calloc(NUM_FILTER_STEPS_FOR_Ts*Nsteps_zp,sizeof(double *)); //New
	for(i=0;i<NUM_FILTER_STEPS_FOR_Ts*Nsteps_zp;i++){  // New
			log10_FcollzXray_SFR_low_table[i] = (double *)calloc(NSFR_low,sizeof(double));
      log10_FcollzLya_SFR_low_table[i] = (double *)calloc(NSFR_low,sizeof(double));
	}

	Overdense_high_table = calloc(NSFR_high,sizeof(float));
	FcollzXray_SFR_high_table = (float **)calloc(NUM_FILTER_STEPS_FOR_Ts*Nsteps_zp,sizeof(float *));
  FcollzLya_SFR_high_table = (float **)calloc(NUM_FILTER_STEPS_FOR_Ts*Nsteps_zp,sizeof(float *)); //New
	for(i=0;i<NUM_FILTER_STEPS_FOR_Ts*Nsteps_zp;i++){  // New
			FcollzXray_SFR_high_table[i] = (float *)calloc(NSFR_high,sizeof(float));
      FcollzLya_SFR_high_table[i] = (float *)calloc(NSFR_high,sizeof(float));
	}
}

void destroy_21cmMC_arrays() {

	int i,j,ithread;

    free(xi_SFR);
    free(wi_SFR);

	free(zpp_interp_table);
	free(redshift_interp_table);

	for(i=0;i<NUM_FILTER_STEPS_FOR_Ts*Nsteps_zp;i++) {
		free(log10_FcollzXray_SFR_low_table[i]);
    free(log10_FcollzLya_SFR_low_table[i]);
	}
	free(log10_FcollzXray_SFR_low_table);
  free(log10_FcollzLya_SFR_low_table);

	free(log10_overdense_low_table);

	for(i=0;i<NUM_FILTER_STEPS_FOR_Ts*Nsteps_zp;i++) {
		free(FcollzXray_SFR_high_table[i]);
    free(FcollzLya_SFR_high_table[i]);
	}
	free(FcollzXray_SFR_high_table);
  free(FcollzLya_SFR_high_table);

	free(Overdense_high_table);

    for (i=0; i < NUM_FILTER_STEPS_FOR_Ts; i++){
      gsl_spline_free (FcollLowXray_zpp_spline[i]);
      gsl_spline_free (FcollLowLya_zpp_spline[i]);
      gsl_interp_accel_free (FcollLowXray_zpp_spline_acc[i]);
      gsl_interp_accel_free (FcollLowLya_zpp_spline_acc[i]);
      free(second_derivs_FcollXray_zpp[i]);
      free(second_derivs_FcollLya_zpp[i]);
    }
}

/* Maximum allowed value for the kinetic temperature.
   Useful to set to avoid some spurious behaviour
   when the code is run with redshift poor resolution
   and very high X-ray heating efficiency */
#define MAX_TK (float) 5e4


int main(int argc, char ** argv){
  fftwf_complex *box, *unfiltered_box;
  fftwf_plan plan;
  unsigned long long ct, sample_ct;
  int R_ct,i,ii,j,k, COMPUTE_Ts, x_e_ct;
  float REDSHIFT, growth_factor_z, R, R_factor, zp, mu_for_Ts, filling_factor_of_HI_zp;
  int ithread;
  float *Tk_box, *x_e_box, *Ts, *Lya_box,  J_star_Lya, dzp, prev_zp, zpp, prev_zpp, prev_R;
  FILE *F, *GLOBAL_EVOL, *OUT;
  char filename[500];
  float dz, zeta_ion_eff, Tk_BC, xe_BC, nu, zprev, zcurr, curr_delNL0[NUM_FILTER_STEPS_FOR_Ts];
  double *evolve_ans, ans[2], dansdz[5], Tk_ave, J_alpha_ave, xalpha_ave, J_alpha_tot, Xheat_ave,
    Xion_ave;
double freq_int_heat_tbl[x_int_NXHII][NUM_FILTER_STEPS_FOR_Ts], freq_int_ion_tbl[x_int_NXHII][NUM_FILTER_STEPS_FOR_Ts], freq_int_lya_tbl[x_int_NXHII][NUM_FILTER_STEPS_FOR_Ts];
  int goodSteps,badSteps;
  int m_xHII_low, m_xHII_high, n_ct, zp_ct;
	double freq_int_heat[NUM_FILTER_STEPS_FOR_Ts], freq_int_ion[NUM_FILTER_STEPS_FOR_Ts], freq_int_lya[NUM_FILTER_STEPS_FOR_Ts];
	double nuprime, fcoll_R, fcoll_R_Lya, Ts_ave;
	float *delNL0[NUM_FILTER_STEPS_FOR_Ts], xHII_call, curr_xalpha;
  double mean_metallicity_collapsed; // Used for SED evolution
	float z, Jalpha, TK, TS, xe, deltax, Lya;
	time_t start_time, curr_time;
 	double J_alpha_threads[NUMCORES], xalpha_threads[NUMCORES], Xheat_threads[NUMCORES], Xion_threads[NUMCORES], lower_int_limit;
	float Splined_Fcollzp_mean, Splined_Fcollzpp_X_mean,ION_EFF_FACTOR,fcoll, fcollLya,Splined_Fcollzpp_Lya_mean; // New in v1.4
	float zp_table; //New in v1.4
	int counter,arr_num; // New in v1.4
  int counter2; // For fixing kinky Lya background
  double Luminosity_conversion_factor;
  int RESTART = 0;
  double Tback;
  // Declaring new variables for the kinky fix -GS
  int n_pts_radii;
  double trial_zpp_min,trial_zpp_max,trial_zpp, weight;
  bool first_radii, first_zero;
  n_pts_radii = 1000;

	// Declaring line luminosity tables -GS
	//float lya_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM],
	//			ha_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM],
	//			hb_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM],
	//			he2_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM],
	//			o3_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM],
	//			o2s_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM],
	//			o2l_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM],
	//			c2_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM],
	//			co10_table_pop2[DIM_POP2_METALLICITY][DIM_POP2_ION_PARAM];


 /**********  BEGIN INITIALIZATION   **************************************/
 //New in v1.4
 if (SHARP_CUTOFF) {
   if (argc == 3){
     RESTART = 1;
     zp = atof(argv[2]);
   }
   else if (argc != 2){
     //fprintf(stderr, "Usage: Ts <REDSHIFT>  [reload zp redshift]\nAborting...\n");
     return -1;
   }
   HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY = 0;
   X_LUMINOSITY = pow(10.,L_X);
   F_STAR10 = STELLAR_BARYON_FRAC;
 }
 else {
   if (argc  == 10) {
     RESTART = 1;
     zp = atof(argv[2]);
     F_STAR10 = atof(argv[3]);
     ALPHA_STAR = atof(argv[4]);
     F_ESC10 = atof(argv[5]);
     ALPHA_ESC = atof(argv[6]);
     M_TURN = atof(argv[7]);
     T_AST = atof(argv[8]);
     X_LUMINOSITY = pow(10.,atof(argv[9]));
   }
   else if (argc == 9) {
     F_STAR10 = atof(argv[2]);
     ALPHA_STAR = atof(argv[3]);
     F_ESC10 = atof(argv[4]);
     ALPHA_ESC = atof(argv[5]);
     M_TURN = atof(argv[6]);
     T_AST = atof(argv[7]);
     X_LUMINOSITY = pow(10.,atof(argv[8]));
   }
   else if (argc == 3) {
     RESTART = 1;
     zp = atof(argv[2]);
     F_STAR10 = STELLAR_BARYON_FRAC;
     ALPHA_STAR = STELLAR_BARYON_PL;
     F_ESC10 = ESC_FRAC;
     ALPHA_ESC = ESC_PL;
     M_TURN = M_TURNOVER;
     T_AST = t_STAR;
     X_LUMINOSITY = pow(10.,L_X);
   }
   else if (argc == 2) {
     F_STAR10 = STELLAR_BARYON_FRAC;
     ALPHA_STAR = STELLAR_BARYON_PL;
     F_ESC10 = ESC_FRAC;
     ALPHA_ESC = ESC_PL;
     M_TURN = M_TURNOVER;
     T_AST = t_STAR;
     X_LUMINOSITY = pow(10.,L_X);
   }
   else {
     //fprintf(stderr, "Usage: Ts <REDSHIFT> [reload zp redshift] [<f_star10> <alpha_star> <f_esc10> <alpha_esc> <M_turn> <t_star> <X_luminosity>] \nAborting...\n");
     return -1;
   }
   HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY = 1;
   ION_EFF_FACTOR = N_GAMMA_UV * F_STAR10 * F_ESC10;
 }
 // Initialize source structure - RM
 sources src;
 src = defaultSources();

 M_MIN = M_TURNOVER;
 REDSHIFT = atof(argv[1]);

init_ps();
// Initialization of galaxy model lookup tables -GS
init_galaxymodel_tbl();
// Initialization of Nion lookup table model -GS
init_pop2_nion_table();
// Initialization of line luminosity lookup tables -GS
//init_pop2_full_tables(lya_table_pop2, ha_table_pop2, hb_table_pop2, he2_table_pop2,
//											o3_table_pop2, o2s_table_pop2, o2l_table_pop2, c2_table_pop2, co10_table_pop2);



// - RM
 if(USE_GENERAL_SOURCES) ION_EFF_FACTOR = 1.0;
//printf("11\n");
 // Set Min Mass if necessary - RM
/*
 if(USE_GENERAL_SOURCES)
 {
     if(src.minMassIII(REDSHIFT) > src.minMass(REDSHIFT)
         || src.minMassIII(REDSHIFT) < 0) M_MIN = src.minMass(REDSHIFT);
     else M_MIN = src.minMassIII(REDSHIFT);
 }*/
if(USE_GENERAL_SOURCES) M_MIN = src.minMass(REDSHIFT) * 50.0; //Multipied by 50 just to test mean PS collapse fraction with fiducial 21cmFAST... this variable doesn't actually matter


 system("mkdir ../Log_files");
 system("mkdir ../Output_files");
 system("mkdir ../Boxes/Ts_evolution/");
 system("mkdir ../Output_files/Ts_outs/");
 system("cp ../Parameter_files/* ../Output_files/Ts_outs/");
 system("cp ../Parameter_files/* ../Boxes/Ts_evolution/");
 omp_set_num_threads(NUMCORES);
 growth_factor_z = dicke(REDSHIFT);


 // open log file
 if (!(LOG = fopen("../Log_files/Ts_log", "w") ) ){
   fprintf(stderr, "Unable to open log file for writting\nAborting...\n");
   return -1;
 }

 // Initialize some interpolation tables
 if (init_heat() < 0){
   fclose(LOG); fclose(GLOBAL_EVOL);
   return -1;
 }

 // check if we are in the really high z regime before the first stars; if so, simple
 if (REDSHIFT > Z_HEAT_MAX){
//(FgtrM(REDSHIFT, FMAX(TtoM(REDSHIFT, X_RAY_Tvir_MIN, mu_for_Ts),  M_MIN_WDM)) < 1e-15 ){
   xe = xion_RECFAST(REDSHIFT,0);
   TK = T_RECFAST(REDSHIFT,0);

   // open input
   sprintf(filename, "../Boxes/updated_smoothed_deltax_z%06.2f_%i_%.0fMpc",
	   REDSHIFT, HII_DIM, BOX_LEN);
   F = fopen(filename, "rb");
   if ( !(F = fopen(filename, "rb") ) ){
     fprintf(stderr, "Error opening file %s for reading.\nAborting...\n", filename);
     fprintf(LOG, "Error opening file %s for reading.\nAborting...\n", filename);
     destruct_heat(); return -1;
   }
   fprintf(stderr, "Opened density file %s for reading\n", filename);
   fprintf(LOG, "Opened density file %s for reading\n", filename);

   // open output
   // New in v1.4
   if (HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY) {
   sprintf(filename, "../Boxes/Ts_z%06.2f_L_X%.1e_alphaX%.1f_f_star%06.4f_alpha_star%06.4f_MturnX%.1e_t_star%06.4f_Pop%i_%i_%.0fMpc", REDSHIFT, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, M_TURN, T_AST, Pop, HII_DIM, BOX_LEN);
   }
   else {
   sprintf(filename, "../Boxes/Ts_z%06.2f_L_X%.1e_alphaX%.1f_MminX%.1e_zetaIon%.2f_Pop%i_%i_%.0fMpc", REDSHIFT, X_LUMINOSITY, X_RAY_SPEC_INDEX, M_MIN, HII_EFF_FACTOR, Pop, HII_DIM, BOX_LEN);
   }
   if (!(OUT=fopen(filename, "wb"))){
     fprintf(stderr, "Ts.c: WARNING: Unable to open output file %s\n", filename);
     fprintf(LOG, "Ts.c: WARNING: Unable to open output file %s\n", filename);
     destruct_heat(); return -1;
   }
   fprintf(stderr, "Opened TS file %s for writting\n", filename);
   fprintf(LOG, "Opened TS file %s for writting\n", filename);

   // read file
   for (i=0; i<HII_DIM; i++){
     for (j=0; j<HII_DIM; j++){
       for (k=0; k<HII_DIM; k++){
	 if (fread(&deltax, sizeof(float), 1, F)!=1){
	  fprintf(stderr, "Error reading-in binary density file\nAborting...\n");
	  fprintf(LOG, "Error reading-in binary density file\nAborting...\n");
	  destruct_heat(); return -1;
	 }

	 // compute the spin temperature
	TS = get_Ts(REDSHIFT, deltax, TK, xe, 0, &curr_xalpha);

	// and print it out
	if (fwrite(&TS, sizeof(float), 1, OUT)!=1){
	  fprintf(stderr, "Ts.c: Write error occured while writting Tk box.\n");
	  fprintf(LOG, "Ts.c: Write error occured while writting Tk box.\n");
	  destruct_heat(); return -1;
	 }

       }
     }
   }

    printf("2\n");

   destruct_heat(); fclose(F); fclose(OUT);
   return 0;
 }


 // open global evolution output file
 // New in v1.4
 if (HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY) {
 sprintf(filename, "../Output_files/Ts_outs/global_evolution_Nsteps%i_zprimestepfactor%.3f_L_X%.1e_alphaX%.1f_f_star10%06.4f_alpha_star%06.4f_f_esc10%06.4f_alpha_esc%06.4f_Mturn%.1e_t_star%06.4f_Pop%i_%i_%.0fMpc", NUM_FILTER_STEPS_FOR_Ts, ZPRIME_STEP_FACTOR, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, T_AST, Pop, HII_DIM, BOX_LEN);
   if (argc == 3 || argc == 9) // restarting
     GLOBAL_EVOL = fopen(filename, "a");
   else
     GLOBAL_EVOL = fopen(filename, "w");
 }
 else {
 sprintf(filename, "../Output_files/Ts_outs/global_evolution_zetaIon%.2f_Nsteps%i_zprimestepfactor%.3f_L_X%.1e_alphaX%.1f_TvirminX%.1e_Pop%i_%i_%.0fMpc", HII_EFF_FACTOR, NUM_FILTER_STEPS_FOR_Ts, ZPRIME_STEP_FACTOR, X_LUMINOSITY, X_RAY_SPEC_INDEX, M_TURN, Pop, HII_DIM, BOX_LEN);
   if (argc > 2) // restarting
     GLOBAL_EVOL = fopen(filename, "a");
   else
     GLOBAL_EVOL = fopen(filename, "w");
 }
 if (!GLOBAL_EVOL){
   fprintf(stderr, "Unable to open global evolution file at %s\nAborting...\n",
	   filename);
   fprintf(LOG, "Unable to open global evolution file at %s\nAborting...\n",
	   filename);
   fclose(LOG);
   return -1;
 }

 // set boundary conditions for the evolution equations->  values of Tk and x_e at Z_HEAT_MAX
 if (XION_at_Z_HEAT_MAX > 0) // user has opted to use his/her own value
   xe_BC = XION_at_Z_HEAT_MAX;
 else// will use the results obtained from recfast
   xe_BC = xion_RECFAST(Z_HEAT_MAX,0);
 if (TK_at_Z_HEAT_MAX > 0)
   Tk_BC = TK_at_Z_HEAT_MAX;
 else
   Tk_BC = T_RECFAST(Z_HEAT_MAX,0);


  /******  Now allocate large arrays  ******/

  // allocate memory for the nonlinear density field and open file
  sprintf(filename, "../Boxes/updated_smoothed_deltax_z%06.2f_%i_%.0fMpc",
	  REDSHIFT, HII_DIM, BOX_LEN);
  F = fopen(filename, "rb");
  if ( !(F = fopen(filename, "rb") ) ){
    fprintf(stderr, "Error opening file %s for reading.\nAborting...\n", filename);
    fprintf(LOG, "Error opening file %s for reading.\nAborting...\n", filename);
    fclose(LOG); fclose(GLOBAL_EVOL);
    destruct_heat();
    return -1;
  }
  if (!(box = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS))){
    fprintf(stderr, "Error in memory allocation for %s\nAborting...\n", filename);
    fprintf(LOG, "Error in memory allocation for %s\nAborting...\n", filename);
    fclose(LOG); fclose(GLOBAL_EVOL); fclose(F);   destruct_heat();
    return -1;
  }
  if (!(unfiltered_box = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS))){
    fprintf(stderr, "Error in memory allocation for %s\nAborting...\n", filename);
    fprintf(LOG, "Error in memory allocation for %s\nAborting...\n", filename);
    fclose(LOG); fclose(GLOBAL_EVOL);fclose(F);   destruct_heat(); fftwf_free(box);
    return -1;
  }
  //fprintf(stderr, "Reading in deltax box\n");
  //fprintf(LOG, "Reading in deltax box\n");
  for (i=0; i<HII_DIM; i++){
    for (j=0; j<HII_DIM; j++){
      for (k=0; k<HII_DIM; k++){
	if (fread((float *)unfiltered_box + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, F)!=1){
	  fprintf(stderr, "Error reading-in binary file %s\nAborting...\n", filename);
	  fprintf(LOG, "Error reading-in binary file %s\nAborting...\n", filename);
	  fftwf_free(box); fclose(GLOBAL_EVOL); fclose(F); fclose(LOG); fftwf_free(unfiltered_box);
	  destruct_heat();
	  return -1;
	}
      }
    }
  }
  fclose(F);


  /*** Transform unfiltered box to k-space to prepare for filtering ***/
  //fprintf(stderr, "begin initial ffts, time=%06.2f min\n", (double)clock()/CLOCKS_PER_SEC/60.0);
  //fprintf(LOG, "begin initial ffts, time=%06.2f min\n", (double)clock()/CLOCKS_PER_SEC/60.0);
  plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)unfiltered_box, (fftwf_complex *)unfiltered_box, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
  fftwf_cleanup();
  // remember to add the factor of VOLUME/TOT_NUM_PIXELS when converting from real space to k-space
  // Note: we will leave off factor of VOLUME, in anticipation of the inverse FFT below
  for (ct=0; ct<HII_KSPACE_NUM_PIXELS; ct++){
    unfiltered_box[ct] /= (float)HII_TOT_NUM_PIXELS;
  }
  //fprintf(stderr, "end initial ffts, time=%06.2f min\n", (double)clock()/CLOCKS_PER_SEC/60.0);
  //fprintf(LOG, "end initial ffts, time=%06.2f min\n", (double)clock()/CLOCKS_PER_SEC/60.0);


  /*** Create the z=0 non-linear density fields smoothed on scale R to be used in computing fcoll ***/
  R = L_FACTOR*BOX_LEN/(float)HII_DIM;
  R_factor = pow(R_XLy_MAX/R, 1/(float)NUM_FILTER_STEPS_FOR_Ts);
  //  R_factor = pow(E, log(HII_DIM)/(float)NUM_FILTER_STEPS_FOR_Ts);
  for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
    R_values[R_ct] = R;
    sigma_atR[R_ct] = sigma_z0(RtoM(R));
    //fprintf(stderr, "Processing scale R= %06.2fMpc, time=%06.2f min\n", R,
	   // (double)clock()/CLOCKS_PER_SEC/60.0);
    //fprintf(LOG, "Processing scale R= %06.2fMpc, time=%06.2f min\n", R,
	   // (double)clock()/CLOCKS_PER_SEC/60.0);
    if (! (delNL0[R_ct] = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS))){
      //fprintf(stderr, "Error in memory allocation\nAborting...\n");
      //fprintf(LOG, "Error in memory allocation\nAborting...\n");
      fclose(LOG); fclose(GLOBAL_EVOL);fftwf_free(box);  fftwf_free(unfiltered_box);
      for(ct=0; ct<R_ct; ct++)
	free(delNL0[ct]);
      destruct_heat();
      return -1;
    }

    // copy over unfiltered box
    memcpy(box, unfiltered_box, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    if (R_ct > 0){ // don't filter on cell size
      HII_filter(box, HEAT_FILTER, R);
    }

    // now fft back to real space
    plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)box, (float *)box, FFTW_ESTIMATE);
    fftwf_execute(plan);

    // copy over the values
    for (i=0; i<HII_DIM; i++){
      for (j=0; j<HII_DIM; j++){
	for (k=0; k<HII_DIM; k++){
	  delNL0[R_ct][HII_R_INDEX(i,j,k)] = *((float *) box + HII_R_FFT_INDEX(i,j,k));
	  if (delNL0[R_ct][HII_R_INDEX(i,j,k)] < -1){ // correct for alliasing in the filtering step
	    delNL0[R_ct][HII_R_INDEX(i,j,k)] = -1+FRACT_FLOAT_ERR;
	  }
	  // and linearly extrapolate to z=0
	  delNL0[R_ct][HII_R_INDEX(i,j,k)] /= growth_factor_z;
	}
      }
    }

    R *= R_factor;
  } //end for loop through the filter scales R

  fftwf_destroy_plan(plan);
  fftwf_cleanup();
  fftwf_free(box); fftwf_free(unfiltered_box);// we don't need this anymore

  // now lets allocate memory for our kinetic temperature and residual neutral fraction boxes
  if (!(Tk_box = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS))){
    fprintf(stderr, "Error in memory allocation for Tk box\nAborting...\n");
    fprintf(LOG, "Error in memory allocation for Tk box\nAborting...\n");
    fclose(LOG);fclose(GLOBAL_EVOL);
    for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
      free(delNL0[R_ct]);
    }
    destruct_heat();
    return -1;
  }
  if (!(x_e_box = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS))){
    fprintf(stderr, "Error in memory allocation for xe box\nAborting...\n");
    fprintf(LOG, "Error in memory allocation for xe box\nAborting...\n");
    fclose(LOG);  free(Tk_box);fclose(GLOBAL_EVOL);
    for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
      free(delNL0[R_ct]);
    }
    destruct_heat();
    return -1;
  }

  // and finally allocate memory for the spin temperature box
  if (!(Ts = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS))){
    fprintf(stderr, "Error in memory allocation for Ts box\nAborting...\n");
    fprintf(LOG, "Error in memory allocation for Ts box\nAborting...\n");
    fclose(LOG);  fclose(GLOBAL_EVOL);free(Tk_box); free(x_e_box);
    for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
      free(delNL0[R_ct]);
    }
    destruct_heat();
    return -1;
  }

  // additionally allocate memory for the Lya radiation box
  if (!(Lya_box = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS))){
    fprintf(stderr, "Error in memory allocation for Lya box\nAborting...\n");
    fprintf(LOG, "Error in memory allocation for Lya box\nAborting...\n");
    fclose(LOG);  fclose(GLOBAL_EVOL);free(Tk_box); free(x_e_box); free(Ts);
    for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
      free(delNL0[R_ct]);
    }
    destruct_heat();
    return -1;
  }

  // and initialize to the boundary values at Z_HEAT_END
  if (!RESTART){ // we are not restarting
    for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
      Tk_box[ct] = Tk_BC;
      x_e_box[ct] = xe_BC;
    }
    x_e_ave = xe_BC;
    Tk_ave = Tk_BC;


    printf("Starting at at z_max=%f, Tk=%f, x_e=%e\n", Z_HEAT_MAX, Tk_ave, x_e_ave);
  }
  else{ // we need to load the evolution files from the intermediate output
    // first Tk
	if (HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY) { // New in v1.4
    sprintf(filename, "../Boxes/Ts_evolution/Tk_zprime%06.2f_L_X%.1e_alphaX%.1f_f_star10_%06.4f_alpha_star%06.4f_f_esc10_%06.4f_alpha_esc%06.4f_Mturn%.1e_t_star%06.4f_Pop%i_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, T_AST, Pop, HII_DIM, BOX_LEN);
	}
	else {
    sprintf(filename, "../Boxes/Ts_evolution/Tk_zprime%06.2f_L_X%.1e_alphaX%.1f_Mmin%.1e_zetaIon%.2f_Pop%i_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, M_MIN, HII_EFF_FACTOR, Pop, HII_DIM, BOX_LEN);
	}
    if (!(F=fopen(filename, "rb"))){
      fprintf(stderr, "Ts.c: WARNING: Unable to open input file %s\nAborting\n", filename);
      fprintf(LOG, "Ts.c: WARNING: Unable to open input file %s\nAborting\n", filename);
      fclose(LOG); fclose(GLOBAL_EVOL); free(Tk_box); free(x_e_box); free(Ts); free(Lya_box);
      for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
	free(delNL0[R_ct]);
      }
      destruct_heat();
      return -1;
    }
    else{
      if (mod_fread(Tk_box, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
	fprintf(stderr, "Ts.c: Write error occured while reading Tk box.\nAborting\n");
	fprintf(LOG, "Ts.c: Write error occured while reading Tk box.\nAborting\n");
	fclose(LOG);  free(Tk_box); free(x_e_box); free(Ts); free(Lya_box);
	for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
	  free(delNL0[R_ct]);
	}
	destruct_heat();
      }
      fclose(F);
    }

    printf("4\n");
    // then xe_neutral
	// New in v1.4
    if (HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY) {
    sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_L_X%.1e_alphaX%.1f_f_star10_%06.4f_alpha_star%06.4f_f_esc10_%06.4f_alpha_esc%06.4f_Mturn%.1e_t_star%06.4f_Pop%i_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, T_AST, Pop, HII_DIM, BOX_LEN);
	}
	else {
    sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_L_X%.1e_alphaX%.1f_Mmin%.1e_zetaIon%.2f_Pop%i_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, M_MIN, HII_EFF_FACTOR, Pop, HII_DIM, BOX_LEN);
	}
      if (!(F=fopen(filename, "rb"))){
      fprintf(stderr, "Ts.c: WARNING: Unable to open output file %s\nAborting\n", filename);
      fprintf(LOG, "Ts.c: WARNING: Unable to open output file %s\nAborting\n", filename);
      fclose(LOG);  free(Tk_box); free(x_e_box); free(Ts); free(Lya_box);
      for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
	free(delNL0[R_ct]);
      }
      destruct_heat();
    }
    else{
      if (mod_fread(x_e_box, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
	fprintf(stderr, "Ts.c: Write error occured while reading xe box.\n");
	fprintf(LOG, "Ts.c: Write error occured while reading xe box.\n");
	fclose(LOG);  free(Tk_box); free(x_e_box); free(Ts); free(Lya_box);
	for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
	  free(delNL0[R_ct]);
	}
	destruct_heat();
      }
      fclose(F);
    }
    Tk_ave = x_e_ave = 0;
    for (box_ct=0; box_ct<HII_TOT_NUM_PIXELS; box_ct++){
      Tk_ave += Tk_box[box_ct];
      x_e_ave += x_e_box[box_ct];
    }
    Tk_ave /= (double) HII_TOT_NUM_PIXELS;


    x_e_ave /= (double) HII_TOT_NUM_PIXELS;
    fprintf(stderr, "Rebooting from z'=%f output. <Tk> = %f. <xe> = %e\n", zp, Tk_ave, x_e_ave);

  }

  /***************    END INITIALIZATION   *********************************/

  /*********  FOR DEBUGGING, set IGM to be homogeneous for testing purposes
    for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++)
    for (box_ct=0; box_ct<HII_TOT_NUM_PIXELS;box_ct++)
      delNL0[R_ct][box_ct] = 0;
    /*  *********/

  // main trapezoidal integral over z' (see eq. ? in Mesinger et al. 2009)
  if (!RESTART){
	Nsteps_zp = 0;
    zp = REDSHIFT*1.0001; //higher for rounding
    while (zp < Z_HEAT_MAX) {
	  Nsteps_zp += 1;
      zp = ((1+zp)*ZPRIME_STEP_FACTOR - 1);
	}
    prev_zp = Z_HEAT_MAX;
  }
  else{
    prev_zp = zp;
  }

  zp = ((1+zp)/ ZPRIME_STEP_FACTOR - 1);
  dzp = zp - prev_zp;
  zp_ct=0;
  COMPUTE_Ts = 0;
  /* New in v1.4: set up interpolation table for computing f_coll(z,delta) */
  if (HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY) {
    init_21cmMC_arrays();

	// Find the highest and lowest redshfit to initialise interpolation of the mean collapse fraction for the global reionization.
    determine_zpp_min = REDSHIFT*0.999;
    for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
        if (R_ct==0){
            prev_zpp = zp;
            prev_R = 0;
        }
        else{
            prev_zpp = zpp_edge[R_ct-1];
            prev_R = R_values[R_ct-1];
        }
        zpp_edge[R_ct] = prev_zpp - (R_values[R_ct] - prev_R)*CMperMPC / drdz(prev_zpp); // cell size
        zpp = (zpp_edge[R_ct]+prev_zpp)*0.5; // average redshift value of shell: z'' + 0.5 * dz''
    }
    determine_zpp_max = zpp*1.001;

	zpp_bin_width = (determine_zpp_max - determine_zpp_min)/((float)zpp_interp_points-1.0);
	for (i=0; i<zpp_interp_points;i++) {
	    zpp_interp_table[i] = determine_zpp_min + (determine_zpp_max - determine_zpp_min)*(float)i/((float)zpp_interp_points-1.0);
	}

	/* initialise interpolation of the mean collapse fraction for global reionization.
	   compute 'FgtrM_st_SFR' corresponding to an array of redshift. */
    initialise_FgtrM_st_SFR_spline(zpp_interp_points,determine_zpp_min, determine_zpp_max, M_TURN, ALPHA_STAR, ALPHA_ESC, F_STAR10, F_ESC10);
    //printf("\n Completed initialise Fcoll_spline Time = %06.2f min \n",(double)clock()/CLOCKS_PER_SEC/60.0);

	/* initialise interpolation of the mean collapse fraction with respect to the X-ray heating.
	   compute 'FgtrM_st_SFR' corresponding to an array of redshift, but assume f_{esc10} = 1 and \alpha_{esc} = 0. */
    initialise_Xray_FgtrM_st_SFR_spline(zpp_interp_points,determine_zpp_min, determine_zpp_max, M_TURN, ALPHA_STAR, F_STAR10);
    //printf("\n Completed initialise Fcoll_spline for X-ray Time = %06.2f min \n",(double)clock()/CLOCKS_PER_SEC/60.0);

    initialise_Lya_FgtrM_st_SFR_spline(zpp_interp_points,determine_zpp_min, determine_zpp_max, M_TURN, ALPHA_STAR, F_STAR10);
    //printf("\n Completed initialise Fcoll_spline for Lya Time = %06.2f min \n",(double)clock()/CLOCKS_PER_SEC/60.0);

	// initialise redshift table corresponding to all the redshifts to initialise interpolation for the conditional mass function.
    zp_table = zp;
	counter = 0;
    for (i=0; i<Nsteps_zp; i++) {
      for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
          if (R_ct==0){
              prev_zpp = zp_table;
              prev_R = 0;
          }
          else{
              prev_zpp = zpp_edge[R_ct-1];
              prev_R = R_values[R_ct-1];
          }
          zpp_edge[R_ct] = prev_zpp - (R_values[R_ct] - prev_R)*CMperMPC / drdz(prev_zpp); // cell size
          zpp = (zpp_edge[R_ct]+prev_zpp)*0.5; // average redshift value of shell: z'' + 0.5 * dz''
		  redshift_interp_table[counter] = zpp;
		  counter += 1;
      }
      prev_zp = zp_table;
	  zp_table = ((1+prev_zp) / ZPRIME_STEP_FACTOR - 1);
    }

	/* generate a table for interpolation of the collapse fraction with respect to the X-ray heating, as functions of
	filtering scale, redshift and overdensity.
	   Note that at a given zp, zpp values depends on the filtering scale R, i.e. f_coll(z(R),delta).
	   Compute the conditional mass function, but assume f_{esc10} = 1 and \alpha_{esc} = 0. */
	initialise_Xray_Fcollz_SFR_Conditional_table(Nsteps_zp,NUM_FILTER_STEPS_FOR_Ts,redshift_interp_table,R_values, M_TURN, ALPHA_STAR, F_STAR10);
	printf("\n Generated the table of conditional mass function = %06.2f min \n",(double)clock()/CLOCKS_PER_SEC/60.0);

  initialise_Lya_Fcollz_SFR_Conditional_table(Nsteps_zp,NUM_FILTER_STEPS_FOR_Ts,redshift_interp_table,R_values, M_TURN, ALPHA_STAR, F_STAR10);
  printf("\n Generated the table of conditional mass function (Lya) = %06.2f min \n",(double)clock()/CLOCKS_PER_SEC/60.0);

  }

  counter = 0;
  while (zp > REDSHIFT){
      Tback = T_background(0, 0, zp);

// - RM
 //if(USE_GENERAL_SOURCES) ION_EFF_FACTOR = ionEff(zp, src);

 //if(USE_GENERAL_SOURCES)
 //{
 //    if(src.minMassIII(zp) > src.minMass(zp)
 //        || src.minMassIII(zp) < 0) M_MIN = src.minMass(zp);
 //    else M_MIN = src.minMassIII(zp);
 //}



	// New in v1.4: initialise interpolation of fcoll over zpp and overdensity.
	if (HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY) {
	  arr_num = NUM_FILTER_STEPS_FOR_Ts*counter; // New
	  for (i=0; i<NUM_FILTER_STEPS_FOR_Ts; i++) {
        gsl_spline_init(FcollLowXray_zpp_spline[i], log10_overdense_low_table, log10_FcollzXray_SFR_low_table[arr_num + i], NSFR_low);
        gsl_spline_init(FcollLowLya_zpp_spline[i], log10_overdense_low_table, log10_FcollzLya_SFR_low_table[arr_num + i], NSFR_low);
        spline(Overdense_high_table-1,FcollzXray_SFR_high_table[arr_num + i]-1,NSFR_high,0,0,second_derivs_FcollXray_zpp[i]-1);
        spline(Overdense_high_table-1,FcollzLya_SFR_high_table[arr_num + i]-1,NSFR_high,0,0,second_derivs_FcollLya_zpp[i]-1);
	  }
	}

    // check if we will next compute the spin temperature (i.e. if this is the final zp step)
    if (Ts_verbose || (((1+zp) / ZPRIME_STEP_FACTOR) < (REDSHIFT+1)) )
      COMPUTE_Ts = 1;
    // check if we are in the really high z regime before the first stars..
	if (HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY) { // New in v1.4
	  FgtrM_st_SFR_z(zp,&(Splined_Fcollzp_mean));
      if ( Splined_Fcollzp_mean < 1e-15 )
        NO_LIGHT = 1;
      else
        NO_LIGHT = 0;
	}
	else {
      if (FgtrM(zp, M_MIN) < 1e-15 )
        NO_LIGHT = 1;
      else
        NO_LIGHT = 0;
	}


	//New in v1.4
	if (HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY) {
	  filling_factor_of_HI_zp = 1 - ION_EFF_FACTOR * Splined_Fcollzp_mean / (1.0 - x_e_ave); // fcoll including f_esc
	}
	else {
	  filling_factor_of_HI_zp = 1 - ION_EFF_FACTOR * FgtrM_st(zp, M_MIN) / (1.0 - x_e_ave);
	}

    if (filling_factor_of_HI_zp > 1) filling_factor_of_HI_zp=1;
    if (filling_factor_of_HI_zp < 0) filling_factor_of_HI_zp=0;

    // let's initialize an array of redshifts (z'') corresponding to the
    // far edge of the dz'' filtering shells
    // and the corresponding minimum halo scale, sigma_Tmin,
    // as well as an array of the frequency integrals
    //fprintf(stderr, "Initializing look-up tables. Time=%06.2f min\n", (double)clock()/CLOCKS_PER_SEC/60.0);
    fprintf(LOG, "Initializing look-up tables. Time=%06.2f min\n", (double)clock()/CLOCKS_PER_SEC/60.0);
    time(&start_time);
    for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
      if (R_ct==0){
	prev_zpp = zp;
	prev_R = 0;
      }
      else{
	prev_zpp = zpp_edge[R_ct-1];
	prev_R = R_values[R_ct-1];
      }

      zpp_edge[R_ct] = prev_zpp - (R_values[R_ct] - prev_R)*CMperMPC / drdz(prev_zpp); // cell size

      zpp = (zpp_edge[R_ct]+prev_zpp)*0.5; // average redshift value of shell: z'' + 0.5 * dz''
	  //if (zpp - redshift_interp_table[arr_num+R_ct] > 1e-3) printf("zpp = %.4f, zpp_array = %.4f\n", zpp, redshift_interp_table[arr_num+R_ct]);
      if(SHARP_CUTOFF) sigma_Tmin[R_ct] =  sigma_z0(M_MIN); // In v1.4 sigma_Tmin doesn't nedd to be an array, just a constant.
      // let's now normalize the total collapse fraction so that the mean is the
      // Sheth-Torman collapse fraction
      fcoll_R = 0;
      fcoll_R_Lya;
      sample_ct=0;
      for (box_ct=0; box_ct<HII_TOT_NUM_PIXELS; box_ct+=(HII_TOT_NUM_PIXELS/1e5+1)){
	sample_ct++;
	// New in v1.4

	if (HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY) {
	  growth_zpp = dicke(zpp);
      //---------- interpolation for fcoll starts ----------
      if (delNL0[R_ct][box_ct]*growth_zpp < 1.5){
        if (delNL0[R_ct][box_ct]*growth_zpp < -1.) {
		  fcoll = 0;
      fcollLya = 0;
        }
        else {
          fcoll = gsl_spline_eval(FcollLowXray_zpp_spline[R_ct], log10(delNL0[R_ct][box_ct]*growth_zpp+1.), FcollLowXray_zpp_spline_acc[R_ct]);
          fcollLya = gsl_spline_eval(FcollLowLya_zpp_spline[R_ct], log10(delNL0[R_ct][box_ct]*growth_zpp+1.), FcollLowLya_zpp_spline_acc[R_ct]);
          fcoll = pow(10., fcoll);
        }
      }
      else {
        if (delNL0[R_ct][box_ct]*growth_zpp < 0.99*Deltac) {
          // Usage of 0.99*Deltac arises due to the fact that close to the critical density, the collapsed fraction becomes a little unstable
          // However, such densities should always be collapsed, so just set f_coll to unity.
          // Additionally, the fraction of points in this regime relative to the entire simulation volume is extremely small.
		  //New
          splint(Overdense_high_table-1,FcollzXray_SFR_high_table[arr_num + R_ct]-1,second_derivs_FcollXray_zpp[R_ct]-1,NSFR_high,delNL0[R_ct][box_ct]*growth_zpp,&(fcoll));
          splint(Overdense_high_table-1,FcollzLya_SFR_high_table[arr_num + R_ct]-1,second_derivs_FcollLya_zpp[R_ct]-1,NSFR_high,delNL0[R_ct][box_ct]*growth_zpp,&(fcollLya));
        }
        else {
		  fcoll = 1.;
      fcollLya = 1.;
        }
      }
      Splined_Fcoll = fcoll;
      Splined_Fcoll_Lya = fcollLya;
      //---------- interpolation for fcoll is done ----------
	  fcoll_R += Splined_Fcoll;
    fcoll_R_Lya += Splined_Fcoll_Lya;
	}
	else {
	  fcoll_R += sigmaparam_FgtrM_bias(zpp, sigma_Tmin[R_ct],
					 delNL0[R_ct][box_ct], sigma_atR[R_ct]);
    }
	  }

      fcoll_R /= (double) sample_ct;

	  if (HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY) {// New in v1.4
	    FgtrM_st_SFR_Xray_z(zpp,&(Splined_Fcollzpp_X_mean));
      FgtrM_st_SFR_Lya_z(zpp,&(Splined_Fcollzpp_Lya_mean));
	    ST_over_PS[R_ct] = Splined_Fcollzpp_X_mean / fcoll_R;
      ST_over_PS_Lya[R_ct] = Splined_Fcollzpp_Lya_mean / fcoll_R;
	  }
	  else {
        ST_over_PS[R_ct] = FgtrM_st(zpp, M_MIN) / fcoll_R;
	  }

      if (DEBUG_ON){
	    if (HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY) {
        printf("zpp = %le M_MIN = %le fcoll_R = %le\n", zpp, M_MIN, fcoll_R);
      printf("ST/PS=%g, mean_ST=%g, mean_ps=%g\n, ratios of mean=%g\n", ST_over_PS[R_ct],
		 Splined_Fcollzpp_X_mean,

	     FgtrM(zpp, M_MIN),
	     Splined_Fcollzpp_X_mean/FgtrM(zpp, M_MIN)
	     );
		}
		else {
      printf("ST/PS=%g, mean_ST=%g, mean_ps=%g\n, ratios of mean=%g\n", ST_over_PS[R_ct],
	     FgtrM_st(zpp, M_MIN),
	     FgtrM(zpp, M_MIN),
	     FgtrM_st(zpp, M_MIN)/FgtrM(zpp, M_MIN)
	     );
		}
      }

      lower_int_limit = FMAX(nu_tau_one(zp, zpp, x_e_ave, filling_factor_of_HI_zp), NU_X_THRESH);


/***************  PARALLELIZED LOOP ******************************************************************/
      // set up frequency integral table for later interpolation for the cell's x_e value
#pragma omp parallel shared(freq_int_heat_tbl, freq_int_ion_tbl, COMPUTE_Ts, freq_int_lya_tbl, zp, R_ct, x_e_ave, x_int_XHII, x_int_Energy, x_int_fheat, x_int_n_Lya, x_int_nion_HI, x_int_nion_HeI, x_int_nion_HeII, lower_int_limit) private(x_e_ct)
{
#pragma omp for
  for (x_e_ct = 0; x_e_ct < x_int_NXHII; x_e_ct++){
    freq_int_heat_tbl[x_e_ct][R_ct] = integrate_over_nu(zp, x_int_XHII[x_e_ct], lower_int_limit, 0);
    freq_int_ion_tbl[x_e_ct][R_ct] = integrate_over_nu(zp, x_int_XHII[x_e_ct], lower_int_limit, 1);
    if (COMPUTE_Ts)
      freq_int_lya_tbl[x_e_ct][R_ct] = integrate_over_nu(zp, x_int_XHII[x_e_ct], lower_int_limit, 2);
  }
} // end omp declaration
/***************  END PARALLELIZED LOOP ******************************************************************/

      // and create the sum over Lya transitions from direct Lyn flux
      sum_lyn[R_ct] = 0;
      for (n_ct=NSPEC_MAX; n_ct>=2; n_ct--){
	if (zpp > zmax(zp, n_ct))
	  continue;

	nuprime = nu_n(n_ct)*(1+zpp)/(1.0+zp);
  mean_metallicity_collapsed = FgtrM_st_ZNUMEPOP2(zpp, M_TURN, ALPHA_STAR, F_STAR10)/FgtrM_st_ZDENOPOP2(zpp, M_TURN, ALPHA_STAR, F_STAR10);
  //printf("CHECK z=%.1f, mean_metallicity_collapsed=%.1e\n", zpp, mean_metallicity_collapsed);

  // remove spectral_emissivity calculation, as it is added into one of the fcoll integrals if USE_GENERAL_SOURCES
  //if(USE_GENERAL_SOURCES) sum_lyn[R_ct] += frecycle(n_ct);
  /*else*/ sum_lyn[R_ct] += frecycle(n_ct) * spectral_emissivity_alt(nuprime, 0, Pop, src.fPopIII(zpp), mean_metallicity_collapsed);
      }
            //printf("We are at zp=%f and R_ct=%d\n", zp, R_ct);
      //printf("Check sum_lyn[R_ct]: %.3e, sum_lyn[R_ct-1]: %.3e, first_radii: %d\n", sum_lyn[R_ct], sum_lyn[R_ct-1], first_radii);
      // Here we adopt the fix of kinky Lya background described in the PR (https://github.com/21cmfast/21cmFAST/pull/230)
      //if(R_ct > 1 && sum_lyn[R_ct]==0.0 && sum_lyn[R_ct-1]>0. && first_radii) {
          // The current zpp for which we are getting zero contribution
        //  trial_zpp_max = (prev_zpp - (R_values[R_ct] - prev_R)*CMperMPC / drdz(prev_zpp)+prev_zpp)*0.5;
          // The zpp for the previous radius for which we had a non-zero contribution
          //trial_zpp_min = (zpp_edge[R_ct-2] - (R_values[R_ct-1] - R_values[R_ct-2])*CMperMPC / drdz(zpp_edge[R_ct-2])+zpp_edge[R_ct-2])*0.5;
          // Split the previous radii and current radii into n_pts_radii smaller radii (redshift) to have fine control of where
          // it transitions from zero to non-zero
          // This is a coarse approximation as it assumes that the linear sampling is a good representation of the different
          // volumes of the shells (from different radii).
          //for(ii=0;ii<n_pts_radii;ii++) {
            //  trial_zpp = trial_zpp_min + (trial_zpp_max - trial_zpp_min)*(float)ii/((float)n_pts_radii-1.);
              //counter2 = 0;
              //for (n_ct=NSPEC_MAX; n_ct>=2; n_ct--){
               //   if (trial_zpp > zmax(zp, n_ct))
                 //     continue;
                  //counter2 += 1;
              //}
              //if(counter2==0&&first_zero) {
               //   first_zero = false;
                //  weight = (float)ii/(float)n_pts_radii;
              //}
          //}
          // Now add a non-zero contribution to the previously zero contribution
          // The amount is the weight, multplied by the contribution from the previous radii
          //sum_lyn[R_ct] = weight * sum_lyn[R_ct-1];
          //first_radii = false;
      //}
    } // end loop over R_ct filter steps
    time(&curr_time);
    //fprintf(stderr, "Finishing initializing look-up tables.  It took %06.2f min on the main thread. Time elapsed (total for all threads)=%06.2f\n", difftime(curr_time, start_time)/60.0, (double)clock()/CLOCKS_PER_SEC/60.0);
    //fprintf(LOG, "Finishing initializing look-up tables.  It took %06.2f min on the main thread. Total time elapsed (total for all threads)=%06.2f\n", difftime(curr_time, start_time)/60.0, (double)clock()/CLOCKS_PER_SEC/60.0);
    fflush(NULL);

    // scroll through each cell and update the temperature and residual ionization fraction
    growth_factor_zp = dicke(zp);
    dgrowth_factor_dzp = ddicke_dz(zp);
    dt_dzp = dtdz(zp);
	// New in v1.4
    // Conversion of the input bolometric luminosity to a ZETA_X, as used to be used in Ts.c
    // Conversion here means the code otherwise remains the same as the original Ts.c
    if(fabs(X_RAY_SPEC_INDEX - 1.0) < 0.000001) {
        Luminosity_conversion_factor = NU_X_THRESH * log( NU_X_BAND_MAX/NU_X_THRESH );
        Luminosity_conversion_factor = 1./Luminosity_conversion_factor;
    }
    else {
        Luminosity_conversion_factor = pow( NU_X_BAND_MAX , 1. - X_RAY_SPEC_INDEX ) - pow( NU_X_THRESH , 1. - X_RAY_SPEC_INDEX ) ;
        Luminosity_conversion_factor = 1./Luminosity_conversion_factor;
        Luminosity_conversion_factor *= pow( NU_X_THRESH, - X_RAY_SPEC_INDEX )*(1 - X_RAY_SPEC_INDEX);
    }
    // Finally, convert to the correct units. NU_over_EV*hplank as only want to divide by eV -> erg (owing to the definition of Luminosity)
    Luminosity_conversion_factor *= (3.1556226e7)/(hplank);
    if(USE_GENERAL_SOURCES)
    {
      const_zp_prefactor = ( X_LUMINOSITY * Luminosity_conversion_factor ) / NU_X_THRESH * C
       * OMb * RHOcrit * pow(CMperMPC, -3) * pow(1+zp, X_RAY_SPEC_INDEX+3);
    }
    else
    {
    const_zp_prefactor = ( X_LUMINOSITY * Luminosity_conversion_factor ) / NU_X_THRESH * C
			 * F_STAR10 * OMb * RHOcrit * pow(CMperMPC, -3) * pow(1+zp, X_RAY_SPEC_INDEX+3);
    }

    /* DEFINE NORMALIZATION FACTORS */
    double norm_lya_input = FgtrM_st_SFR_Lya(zp, M_TURN, ALPHA_STAR, ALPHA_ESC, F_STAR10, F_ESC10, Mlim_Fstar, Mlim_Fesc);
    double norm_xray_input = FgtrM_st_SFR_Xray(zp, M_TURN, ALPHA_STAR, ALPHA_ESC, F_STAR10, F_ESC10, Mlim_Fstar, Mlim_Fesc);

    /*Lya addition -MG
    // open previous lya box
    if (HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY != 0) {
    sprintf(filename, "../Boxes/Lya_z%06.2f_L_X%.1e_alphaX%.1f_f_star%06.4f_alpha_star%06.4f_f_esc%06.4f_alpha_esc%06.4f_Mturn%.1e_t_star%06.4f_Pop%i_%i_%.0fMpc", prev_zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, T_AST, Pop, HII_DIM, BOX_LEN);
    }
  else {
    sprintf(filename, "../Boxes/Lya_z%06.2f_L_X%.1e_alphaX%.1f_TvirminX%.1e_zetaIon%.2f_Pop%i_%i_%.0fMpc", prev_zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, M_TURN, HII_EFF_FACTOR, Pop, HII_DIM, BOX_LEN);
  }
    F = fopen(filename, "rb");
    if ( !(F = fopen(filename, "rb") ) ){
      fprintf(stderr, "Error opening file %s for reading.\n", filename);
      fprintf(LOG, "Error opening file %s for reading.\n", filename);
      for (ct=0; ct<HII_TOT_NUM_PIXELS;ct++)
        Lya_box[ct] = 0;
    }
    else {
      fprintf(stderr, "Loading Lya box with previous redshift values at %s\n", filename);
      fprintf(LOG, "Loading Lya box with previous redshift values at %s\n", filename);
      for (ct=0; ct<HII_TOT_NUM_PIXELS;ct++) {
        if (fread(&Lya, sizeof(float), 1, F)!=1){
        fprintf(stderr, "Error reading previous Lya box\nAborting...\n");
        fprintf(LOG, "Error reading previous Lya box\nAborting...\n");
        destruct_heat(); return -1;
       }
       Lya_box[ct] = Lya;
      }
    }*/
    //for now just keep lya box set to 0 under assumption that evolveInt integrates over z //-MG




    /********  LOOP THROUGH BOX *************/
    //fprintf(stderr, "Looping through box at z'=%f, time elapsed  (total for all threads)= %06.2f min\n", zp, (double)clock()/CLOCKS_PER_SEC/60.0);
    fprintf(LOG, "Looping through box at z'=%f, time elapsed  (total for all threads)= %06.2f min\n", zp, (double)clock()/CLOCKS_PER_SEC/60.0);
    fflush(NULL);
    time(&start_time);
    for (ct=0; ct<NUMCORES; ct++)
      J_alpha_threads[ct] = xalpha_threads[ct] = Xheat_threads[ct] = Xion_threads[ct] = 0;
    /***************  PARALLELIZED LOOP ******************************************************************/
#pragma omp parallel shared(COMPUTE_Ts, Tk_box, x_e_box, x_e_ave, delNL0, freq_int_heat_tbl, freq_int_ion_tbl, freq_int_lya_tbl, zp, dzp, Ts, x_int_XHII, x_int_Energy, x_int_fheat, x_int_n_Lya, x_int_nion_HI, x_int_nion_HeI, x_int_nion_HeII, growth_factor_zp, dgrowth_factor_dzp, NO_LIGHT, zpp_edge, sigma_atR, sigma_Tmin, ST_over_PS, ST_over_PS_Lya, sum_lyn, const_zp_prefactor, M_MIN_at_z, M_MIN_at_zp, dt_dzp, J_alpha_threads, xalpha_threads, Xheat_threads, Xion_threads, Lya_box) private(box_ct, ans, xHII_call, R_ct, curr_delNL0, m_xHII_low, m_xHII_high, freq_int_heat, freq_int_ion, freq_int_lya, dansdz, J_alpha_tot, curr_xalpha)
    {
#pragma omp for

    for (box_ct=0; box_ct<HII_TOT_NUM_PIXELS; box_ct++){
      if (!COMPUTE_Ts && (Tk_box[box_ct] > MAX_TK)) //just leave it alone and go to next value
	continue;

      // set to current values before updating
      ans[0] = x_e_box[box_ct];
      ans[1] = Tk_box[box_ct];

      /*
      if (DEBUG_ON){
	if (isnan(ans[0]))
	  fprintf(stderr, "Problem at cell %llu, x_e=%e\n", box_ct, ans[0]);
	if (isnan(ans[1]))
	  fprintf(stderr, "Problem at cell %llu, Tk=%e\n", box_ct, ans[1]);
      }
      */

      xHII_call = x_e_box[box_ct];

      // Check if ionized fraction is within boundaries; if not, adjust to be within
      if (xHII_call > x_int_XHII[x_int_NXHII-1]*0.999) {
	xHII_call = x_int_XHII[x_int_NXHII-1]*0.999;
      } else if (xHII_call < x_int_XHII[0]) {
	xHII_call = 1.001*x_int_XHII[0];
      }
      //interpolate to correct nu integral value based on the cell's ionization state
      for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
	curr_delNL0[R_ct] = delNL0[R_ct][box_ct];
	m_xHII_low = locate_xHII_index(xHII_call);
	m_xHII_high = m_xHII_low + 1;

	// heat
	freq_int_heat[R_ct] = (freq_int_heat_tbl[m_xHII_high][R_ct] -
			       freq_int_heat_tbl[m_xHII_low][R_ct]) /
	                  (x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
	freq_int_heat[R_ct] *= (xHII_call - x_int_XHII[m_xHII_low]);
	freq_int_heat[R_ct] += freq_int_heat_tbl[m_xHII_low][R_ct];

	// ionization
	freq_int_ion[R_ct] = (freq_int_ion_tbl[m_xHII_high][R_ct] -
			       freq_int_ion_tbl[m_xHII_low][R_ct]) /
	                  (x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
	freq_int_ion[R_ct] *= (xHII_call - x_int_XHII[m_xHII_low]);
	freq_int_ion[R_ct] += freq_int_ion_tbl[m_xHII_low][R_ct];

	// lya
	if (COMPUTE_Ts){
	  freq_int_lya[R_ct] = (freq_int_lya_tbl[m_xHII_high][R_ct] -
				 freq_int_lya_tbl[m_xHII_low][R_ct]) /
	    (x_int_XHII[m_xHII_high] - x_int_XHII[m_xHII_low]);
	  freq_int_lya[R_ct] *= (xHII_call - x_int_XHII[m_xHII_low]);
	  freq_int_lya[R_ct] += freq_int_lya_tbl[m_xHII_low][R_ct];
	}
      }

      /********  finally compute the redshift derivatives *************/
      evolveInt(zp, curr_delNL0, freq_int_heat, freq_int_ion, freq_int_lya,
	  	COMPUTE_Ts, ans, dansdz, arr_num, norm_lya_input, norm_xray_input);//, M_TURN,ALPHA_STAR,F_STAR10,T_AST);

      //update quantities
      x_e_box[box_ct] += dansdz[0] * dzp; // remember dzp is negative
      if (x_e_box[box_ct] > 1) // can do this late in evolution if dzp is too large
	x_e_box[box_ct] = 1 - FRACT_FLOAT_ERR;
      else if (x_e_box[box_ct] < 0)
	x_e_box[box_ct] = 0;
      if (Tk_box[box_ct] < MAX_TK)
	Tk_box[box_ct] += dansdz[1] * dzp;

      if (Tk_box[box_ct]<0){ // spurious bahaviour of the trapazoidalintegrator. generally overcooling in underdensities
	Tk_box[box_ct] = T_cmb*(1+zp);
    //Tk_box[box_ct] = T_background(T_cmb, RADIO_EXCESS_FRAC, zp);
    Tk_box[box_ct] = Tback;
      }

      if (COMPUTE_Ts){
	J_alpha_tot = dansdz[2]; //not really d/dz, but the lya flux
	Ts[box_ct] = get_Ts(zp, curr_delNL0[0]*growth_factor_zp,
	Tk_box[box_ct], x_e_box[box_ct], J_alpha_tot, &curr_xalpha);
  Lya_box[box_ct] = J_alpha_tot * 6 * Ly_alpha_E / pow(1+zp, 4);
	J_alpha_threads[omp_get_thread_num()] += J_alpha_tot;
	xalpha_threads[omp_get_thread_num()] += curr_xalpha;
	Xheat_threads[omp_get_thread_num()] += dansdz[3];
	Xion_threads[omp_get_thread_num()] += dansdz[4];
      }
    }

    } // end parallelization pragma


/***************  END PARALLELIZED LOOP ******************************************************************/
    time(&curr_time);
    //fprintf(stderr, "End scrolling through the box, which took %06.2 min\n", difftime(curr_time, start_time)/60.0);
    fprintf(LOG, "End scrolling through the box, which took %06.2 min\n", difftime(curr_time, start_time)/60.0);
    fflush(NULL);

    // compute new average values
    x_e_ave = 0; Tk_ave = 0; Ts_ave = 0; J_alpha_ave = 0; xalpha_ave = 0; Xheat_ave=0; Xion_ave=0;
    for (box_ct=0; box_ct<HII_TOT_NUM_PIXELS; box_ct++){
      x_e_ave += x_e_box[box_ct];
      Tk_ave += Tk_box[box_ct];
      if (COMPUTE_Ts)
	Ts_ave += Ts[box_ct];
    }
    for (ct=0; ct<NUMCORES; ct++){
      J_alpha_ave += J_alpha_threads[ct];
      xalpha_ave += xalpha_threads[ct];
      Xheat_ave += Xheat_threads[ct];
      Xion_ave += Xion_threads[ct];
    }
    Ts_ave /= (double)HII_TOT_NUM_PIXELS;
    x_e_ave /= (double)HII_TOT_NUM_PIXELS;
    Tk_ave /= (double)HII_TOT_NUM_PIXELS;


    J_alpha_ave /= (double)HII_TOT_NUM_PIXELS;
    xalpha_ave /= (double)HII_TOT_NUM_PIXELS;
    Xheat_ave /= (double)HII_TOT_NUM_PIXELS;
    Xion_ave /= (double)HII_TOT_NUM_PIXELS;
    // write to global evolution file
    fprintf(GLOBAL_EVOL, "%f\t%f\t%f\t%e\t%f\t%f\t%e\t%e\t%e\t%e\n", zp, filling_factor_of_HI_zp, Tk_ave, x_e_ave, Ts_ave, T_cmb*(1+zp), J_alpha_ave, xalpha_ave, Xheat_ave, Xion_ave);
    fflush(NULL);

    // output these intermediate boxes
    if ( Ts_verbose || (++zp_ct >= 10)){ // print every 10th z' evolution step, in case we need to restart
      zp_ct=0;
      //fprintf(stderr, "Writting the intermediate output at zp = %.4f, <Tk>=%f, <x_e>=%e\n", zp, Tk_ave, x_e_ave);
      fprintf(LOG, "Writting the intermediate output at zp = %.4f, <Tk>=%f, <x_e>=%e\n", zp, Tk_ave, x_e_ave);
      fflush(NULL);

      // first Tk
	    // New v1.4
	  if (HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY) {
      sprintf(filename, "../Boxes/Ts_evolution/Tk_zprime%06.2f_L_X%.1e_alphaX%.1f_f_star10_%06.4f_alpha_star%06.4f_f_esc10_%06.4f_alpha_esc%06.4f_Mturn%.1e_t_star%06.4f_Pop%i_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, T_AST, Pop, HII_DIM, BOX_LEN);
	  }
	  else {
      sprintf(filename, "../Boxes/Ts_evolution/Tk_zprime%06.2f_L_X%.1e_alphaX%.1f_Mmin%.1e_zetaIon%.2f_Pop%i_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, M_MIN, HII_EFF_FACTOR, Pop, HII_DIM, BOX_LEN);
	  }
      if (!(F=fopen(filename, "wb"))){
	fprintf(stderr, "Ts.c: WARNING: Unable to open output file %s\n", filename);
	fprintf(LOG, "Ts.c: WARNING: Unable to open output file %s\n", filename);
      }
      else{
	if (mod_fwrite(Tk_box, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
	  fprintf(stderr, "Ts.c: Write error occured while writting Tk box.\n");
	  fprintf(LOG, "Ts.c: Write error occured while writting Tk box.\n");
	}
	fclose(F);
      }
      // then xe_neutral
	    // New in v1.4
	  if (HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY) {
      sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_L_X%.1e_alphaX%.1f_f_star10_%06.4f_alpha_star%06.4f_f_esc10_%06.4f_alpha_esc%06.4f_Mturn%.1e_t_star%06.4f_Pop%i_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, T_AST, Pop, HII_DIM, BOX_LEN);
	  }
	  else {
      sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_L_X%.1e_alphaX%.1f_Mmin%.1e_zetaIon%.2f_Pop%i_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, M_MIN, HII_EFF_FACTOR, Pop, HII_DIM, BOX_LEN);
	  }
      if (!(F=fopen(filename, "wb"))){
	fprintf(stderr, "Ts.c: WARNING: Unable to open output file %s\n", filename);
	fprintf(LOG, "Ts.c: WARNING: Unable to open output file %s\n", filename);
      }
      else{
	if (mod_fwrite(x_e_box, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
	  fprintf(stderr, "Ts.c: Write error occured while writting Tk box.\n");
	  fprintf(LOG, "Ts.c: Write error occured while writting Tk box.\n");
	}
	fclose(F);
      }
    }

    // and the spin temperature if desired
  if ( COMPUTE_Ts ){
    // New in v1.4
    if (HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY != 0) {
    sprintf(filename, "../Boxes/Ts_z%06.2f_L_X%.1e_alphaX%.1f_f_star%06.4f_alpha_star%06.4f_f_esc%06.4f_alpha_esc%06.4f_Mturn%.1e_t_star%06.4f_Pop%i_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, T_AST, Pop, HII_DIM, BOX_LEN);
    }
	else {
    sprintf(filename, "../Boxes/Ts_z%06.2f_L_X%.1e_alphaX%.1f_TvirminX%.1e_zetaIon%.2f_Pop%i_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, M_TURN, HII_EFF_FACTOR, Pop, HII_DIM, BOX_LEN);
	}
      if (!(F=fopen(filename, "wb"))){
	fprintf(stderr, "Ts.c: WARNING: Unable to open output file %s\n", filename);
	fprintf(LOG, "Ts.c: WARNING: Unable to open output file %s\n", filename);
      }
      else{
	if (mod_fwrite(Ts, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
	  fprintf(stderr, "Ts.c: Write error occured while writting Tk box.\n");
	  fprintf(LOG, "Ts.c: Write error occured while writting Tk box.\n");
	}
	fclose(F);
      }



    //Also for Lya box -MG
      if (HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY != 0) {
      sprintf(filename, "../Boxes/Lyabkg_z%06.2f_L_X%.1e_alphaX%.1f_f_star%06.4f_alpha_star%06.4f_f_esc%06.4f_alpha_esc%06.4f_Mturn%.1e_t_star%06.4f_Pop%i_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, T_AST, Pop, HII_DIM, BOX_LEN);
      }
    else {
      sprintf(filename, "../Boxes/Lyabkg_z%06.2f_L_X%.1e_alphaX%.1f_TvirminX%.1e_zetaIon%.2f_Pop%i_%i_%.0fMpc", zp, X_LUMINOSITY, X_RAY_SPEC_INDEX, M_TURN, HII_EFF_FACTOR, Pop, HII_DIM, BOX_LEN);
    }
        if (!(F=fopen(filename, "wb"))){
    fprintf(stderr, "Ts.c: WARNING: Unable to open output file %s\n", filename);
    fprintf(LOG, "Ts.c: WARNING: Unable to open output file %s\n", filename);
        }
        else{
    if (mod_fwrite(Lya_box, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
      fprintf(stderr, "Ts.c: Write error occured while writing Lya box.\n");
      fprintf(LOG, "Ts.c: Write error occured while writing Lya box.\n");
    }
    fclose(F);
        }
    }



    prev_zp = zp;
    zp = ((1+prev_zp) / ZPRIME_STEP_FACTOR - 1);
    dzp = zp - prev_zp;
	counter += 1;
  } // end main integral loop over z'



  //deallocate
  fclose(LOG); fclose(GLOBAL_EVOL); free(Tk_box); free(x_e_box); free(Ts); free(Lya_box);
  if (HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY) {
  	destroy_21cmMC_arrays();
	free_interpolation();
  }
  for (R_ct=0; R_ct<NUM_FILTER_STEPS_FOR_Ts; R_ct++){
    free(delNL0[R_ct]);
  }
  destruct_heat();
  return 0;
}
