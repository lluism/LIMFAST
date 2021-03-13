#include "bubble_helper_progs.c"
#include "heating_helper_progs.c"
#include "compute_line_fields.c" //-MG
#include "../Parameter_files/SOURCES.H"
#include "../Parameter_files/COSMOLOGY.H"

/*
  USAGE: find_HII_bubbles [-p <num of processors>] <redshift> [<previous redshift>]
  [<stellar fraction for 10^10 Msun halos> <power law index for stellar fraction halo mass scaling>
   <escape fraction for 10^10 Msun halos> <power law index for escape fraction halo mass scaling>
   <turn-over scale for the duty cycle of galaxies, in units of halo mass>] [<Soft band X-ray luminosity>]

   final optional parameter is MFP (depricated in v1.4)

  Program FIND_HII_BUBBLES generates the ionization field in Boxes/

  The minimum command line parameter is the redshift and, if the INHOMO_RECO flag is set, also the previous redshift.
  (the previous redshift is only used if computing recombinations)

  If optional parameters are not passed, then the
  values in ANAL_PARAMS.H are used.

  NOTE: the optional argument of thread number including the -p flag, MUST
  be the first two arguments.  If these are omitted, num_threads defaults
  to NUMCORES in INIT_PARAMS.H

  See ANAL_PARAMS.H for updated parameters and algorithms!

  Author: Andrei Mesinger
  Date: 01/10/07

  Support for Alpha power law scaling of zeta with halo mass added by Bradley Greig, 2016
*/

float *Fcoll;
float *Fcoll_uni;
float *Fcoll_fs;
float *Fcoll_sfrxs;
float *mass_box;
float *collapsed_mass;
float *star_mass;
float *SFR;
float *stellar_mass;
float *metallicity;
float *Lya_igm;
float *Lya_sf;
float *Halpha_sf;

void init_21cmMC_arrays() { // defined in Cosmo_c_files/ps.c

    Overdense_spline_SFR = calloc(NSFR_high,sizeof(float)); // New in v1.4
    Overdense_spline_SFR_fs = calloc(NSFR_high,sizeof(float)); //-MG
    Overdense_spline_SFR_sfrxs = calloc(NSFR_high,sizeof(float)); //-LMR
    Fcoll_spline_SFR = calloc(NSFR_high,sizeof(float));
    Fcoll_spline_SFR_fs = calloc(NSFR_high,sizeof(float)); //-MG
    Fcoll_spline_SFR_sfrxs = calloc(NSFR_high,sizeof(float)); //-LMR

    if(USE_GENERAL_SOURCES)
    {
      //printf("initialized zeta spline array\n");
      zeta_spline_SFR = calloc(NSFR_high, sizeof(float));
      zeta_second_derivs_SFR = calloc(NSFR_high, sizeof(float));
    }
    second_derivs_SFR = calloc(NSFR_high,sizeof(float));
    second_derivs_SFR_fs = calloc(NSFR_high,sizeof(float)); //-MG
    second_derivs_SFR_sfrxs = calloc(NSFR_high,sizeof(float)); //-LMR
    xi_SFR = calloc((NGL_SFR+1),sizeof(float));
    wi_SFR = calloc((NGL_SFR+1),sizeof(float));
}

void destroy_21cmMC_arrays() {

    free(Mass_Spline);
    free(Sigma_Spline);
    free(dSigmadm_Spline);
    free(second_derivs_sigma);
    free(second_derivs_dsigma);

    free(Overdense_spline_SFR); // New in v1.4
    free(Overdense_spline_SFR_fs); //-MG
    free(Overdense_spline_SFR_sfrxs); //-LMR
    free(Fcoll_spline_SFR);
    free(Fcoll_spline_SFR_fs); //-MG
    free(Fcoll_spline_SFR_sfrxs); //-LMR
    if(USE_GENERAL_SOURCES)
    {
      free(zeta_spline_SFR);
      free(zeta_second_derivs_SFR);
    }
    free(second_derivs_SFR);
    free(second_derivs_SFR_fs); //-MG
    free(second_derivs_SFR_sfrxs); //-LMR
    free(xi_SFR);
    free(wi_SFR);
}


FILE *LOG;
unsigned long long SAMPLING_INTERVAL = (((unsigned long long)(HII_TOT_NUM_PIXELS/1.0e6)) + 1); //used to sample density field to compute mean collapsed fraction

int parse_arguments(int argc, char ** argv, int * num_th, int * arg_offset, float * F_STAR10,
          float * ALPHA_STAR, float * F_ESC10, float * ALPHA_ESC, float * M_TURN, float * T_AST, double * X_LUMINOSITY, float * MFP,
          float * REDSHIFT, float * PREV_REDSHIFT){

  int min_argc = 2;

  if (INHOMO_RECO)
    min_argc++;

  // check if user wants to specify the multi threading
  if ((argc > min_argc) && (argv[1][0]=='-') && ((argv[1][1]=='p') || (argv[1][1]=='P'))){
    // user specified num proc
    *num_th = atoi(argv[2]);
    *arg_offset = 2;
  }
  else{
    *num_th = NUMCORES;
    *arg_offset = 0;
  }

  // parse the remaining arguments
  if (SHARP_CUTOFF) {
  if (argc == (*arg_offset + min_argc)){
      *MFP = R_BUBBLE_MAX;
      *F_STAR10 = STELLAR_BARYON_FRAC;
      *ALPHA_STAR = STELLAR_BARYON_PL;
      *F_ESC10 = ESC_FRAC;
      *ALPHA_ESC = ESC_PL;
      *M_TURN = M_TURNOVER;
      *T_AST = t_STAR;
    *X_LUMINOSITY = 0;
  }
    else{ return 0;} // format is not allowed
  }
  else {
    if (USE_TS_IN_21CM) {
      if (argc == (*arg_offset + min_argc+7)){
        *F_STAR10 = atof(argv[*arg_offset + min_argc]);
        *ALPHA_STAR = atof(argv[*arg_offset + min_argc+1]);
        *F_ESC10 = atof(argv[*arg_offset + min_argc+2]);
        *ALPHA_ESC = atof(argv[*arg_offset + min_argc+3]);
        *M_TURN = atof(argv[*arg_offset + min_argc+4]); // Input value M_TURN
    *T_AST = atof(argv[*arg_offset + min_argc+5]);
    *X_LUMINOSITY = pow(10.,atof(argv[*arg_offset + min_argc+6]));
      }
      else if (argc == (*arg_offset + min_argc)){ //These parameters give the result which is the same with the default model.
        *F_STAR10 = STELLAR_BARYON_FRAC;         // We need to set the parameters with some standard values later on.
        *ALPHA_STAR = STELLAR_BARYON_PL;
        *F_ESC10 = ESC_FRAC;
        *ALPHA_ESC = ESC_PL;
        *M_TURN = M_TURNOVER;
    *T_AST = t_STAR;
    *X_LUMINOSITY = pow(10.,L_X);
      }
      else{ return 0;} // format is not allowed
  }
    else {
      if (argc == (*arg_offset + min_argc+5)){
        *F_STAR10 = atof(argv[*arg_offset + min_argc]);
        *ALPHA_STAR = atof(argv[*arg_offset + min_argc+1]);
        *F_ESC10 = atof(argv[*arg_offset + min_argc+2]);
        *ALPHA_ESC = atof(argv[*arg_offset + min_argc+3]);
        *M_TURN = atof(argv[*arg_offset + min_argc+4]); // Input value M_TURN
    *T_AST = t_STAR;
    *X_LUMINOSITY = 0;
      }
      else if (argc == (*arg_offset + min_argc)){ //These parameters give the result which is the same with the default model.
        *F_STAR10 = STELLAR_BARYON_FRAC;         // We need to set the parameters with some standard values later on.
        *ALPHA_STAR = STELLAR_BARYON_PL;
        *F_ESC10 = ESC_FRAC;
        *ALPHA_ESC = ESC_PL;
        *M_TURN = M_TURNOVER;
    *T_AST = t_STAR;
    *X_LUMINOSITY = 0;
      }
      else{ return 0;} // format is not allowed
  }
  *MFP = R_BUBBLE_MAX;
  }

  *REDSHIFT = atof(argv[*arg_offset+1]);
  if (INHOMO_RECO){
    *PREV_REDSHIFT = atof(argv[*arg_offset+2]);
    if (*PREV_REDSHIFT <= *REDSHIFT ){
      fprintf(stderr, "find_HII_bubbles: previous redshift must be larger than current redshift!!!\nAborting...\n");
      return -1;
    }
  }
  else
    *PREV_REDSHIFT = *REDSHIFT+0.2; // dummy variable which is not used

   fprintf(stderr, "find_HII_bubbles: command line parameters are as follows\nnum threads=%i, f_star10=%g, alpha_satr=%g, f_esc10=%g, alpha_esc=%g, Mturn=%g, t_star=%g, L_X=%g, z=%g, prev z=%g\n",
    *num_th, *F_STAR10, *ALPHA_STAR, *F_ESC10, *ALPHA_ESC, *M_TURN, *T_AST, *X_LUMINOSITY, *REDSHIFT, *PREV_REDSHIFT);

  return 1;
}



/********** MAIN PROGRAM **********/
int main(int argc, char ** argv){
  char filename[1000], error_message[1000];
  FILE *F = NULL, *pPipe = NULL;
  float REDSHIFT, PREV_REDSHIFT, mass, R, xf, yf, zf, growth_factor, pixel_mass, cell_length_factor, massofscaleR;
  float ave_N_min_cell, ION_EFF_FACTOR, M_MIN;
  int x,y,z, N_min_cell, LAST_FILTER_STEP, num_th, arg_offset, i,j,k;
  unsigned long long ct, ion_ct, sample_ct;
  float f_coll_crit, pixel_volume,  density_over_mean, erfc_num, erfc_denom, erfc_denom_cell, res_xH, Splined_Fcoll, Splined_Fcoll_uni, Splined_Fcoll_fs;
  double coll_mass, cell_mass, val, mean_f_coll_st_uni, ST_over_PS_uni, f_coll_uni; //-MG
  double mean_f_coll_st_fs, ST_over_PS_fs, f_coll_fs;
  double mean_f_coll_st_sfrxs, ST_over_PS_sfrxs, f_coll_sfrxs; //-LMR
  float Splined_Fcoll_sfrxs;
  float metal_val, Halpha_lum, Lya_lum, Lya_emissivity, n_Hii; //-MG
  float Lya_table_pop2[POP2_METALLICITY_SAMPLES+1][POP2_ION_PARAM_SAMPLES+1], Halpha_table_pop2[POP2_METALLICITY_SAMPLES+1][POP2_ION_PARAM_SAMPLES+1]; //-MG
  double table_pop3[POP3_ION_SAMPLES][POP3_COLUMNS]; //-MG
  // - RM
  float Splined_zeta;
  float *xH=NULL, TVIR_MIN, MFP, xHI_from_xrays, std_xrays, *z_re=NULL, *Gamma12=NULL, *mfp=NULL;
  fftwf_complex *M_coll_unfiltered=NULL, *M_coll_filtered=NULL, *deltax_unfiltered=NULL, *deltax_filtered=NULL, *xe_unfiltered=NULL, *xe_filtered=NULL;
  fftwf_complex *N_rec_unfiltered=NULL, *N_rec_filtered=NULL;
  fftwf_plan plan;
  double global_xH, ave_xHI_xrays, ave_den, ST_over_PS, mean_f_coll_st, f_coll, ave_fcoll, dNrec;
  const gsl_rng_type * T=NULL;
  gsl_rng * r=NULL;
  i=0;
  double t_ast, dfcolldt, Gamma_R_prefactor, rec;
  float nua, dnua, temparg, Gamma_R, z_eff;
  float F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, Mlim_Fstar, Mlim_Fesc; //New in v1.4
  double X_LUMINOSITY;
  float global_xH_m, fabs_dtdz, ZSTEP;
  const float dz = 0.01;
  *error_message = '\0';
  //TEST
  double mhalo_i,dmhalo,mhalo_min,mhalo_max,Fesc,Fesc2,Fesc3;
  // - RM

  // Initializtion of source structure - RM
  sources src;
  src = defaultSources();

  //fprintf(stderr, "SFR - - - - - - - - : %f\n", src.fstar(8,1e11));

  int HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY = 0;

  double aveR = 0;
  unsigned long long Rct = 0;

  /* TEST computing time */
  double test1,test2,ratio1,ratio2;
  time_t start, start2, end2, end;
  double time_taken;
  float dum; // TEST
  start = clock();

  /*************************************************************************************/
  /******** BEGIN INITIALIZATION ********/
  /*************************************************************************************/

  // PARSE COMMAND LINE ARGUMENTS
  if(SHARP_CUTOFF){
    if( !parse_arguments(argc, argv, &num_th, &arg_offset, &F_STAR10, &ALPHA_STAR, &F_ESC10,
           &ALPHA_ESC, &M_TURN, &T_AST, &X_LUMINOSITY, &MFP, &REDSHIFT, &PREV_REDSHIFT)){
        fprintf(stderr, "find_HII_bubbles <redshift> [<previous redshift>] \n \
            Aborting...\n                               \
        Check that your inclusion (or not) of [<previous redshift>] is consistent with the INHOMO_RECO flag in ../Parameter_files/ANAL_PARAMS.H\nAborting...\n");
      return -1;
  }
  }
  else {
    HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY = 1;
    if( !parse_arguments(argc, argv, &num_th, &arg_offset, &F_STAR10, &ALPHA_STAR, &F_ESC10,
           &ALPHA_ESC, &M_TURN, &T_AST, &X_LUMINOSITY, &MFP, &REDSHIFT, &PREV_REDSHIFT)){
        fprintf(stderr, "find_HII_bubbles <redshift> [<previous redshift>] \n \
        additional optional arguments: <f_star10> <alpha,star> <f_esc10> <alpha,esc> <M_TURNOVER>] [<t_star>]\n \
        Check that your inclusion (or not) of [<previous redshift>] is consistent with the INHOMO_RECO flag in ../Parameter_files/ANAL_PARAMS.H\n \
   Also check that your inclusion (or not) of [<t_star>] is consistent with the USE_TS_IN_21CM flag in ../Parameter_files/HEAT_PARAMS.H\nAborting...\n");
      return -1;
    }
  }



/*  if ((fabs(ALPHA_STAR) > FRACT_FLOAT_ERR) ||
      (fabs(ALPHA_ESC) > FRACT_FLOAT_ERR) ){ // use the new galaxy parametrization in v1.4 (see ANAL_PARAMS.H)
    HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY = 1;
  }  */

  ZSTEP = PREV_REDSHIFT - REDSHIFT;
  fabs_dtdz = fabs(dtdz(REDSHIFT));
  t_ast = T_AST * t_hubble(REDSHIFT);
  growth_factor = dicke(REDSHIFT);
  pixel_volume = pow(BOX_LEN/(float)HII_DIM, 3);
  pixel_mass = RtoM(L_FACTOR*BOX_LEN/(float)HII_DIM);
  // this parameter choice is sensitive to noise on the cell size, at least for the typical
  // cell sizes in RT simulations.  it probably doesn't matter for larger cell sizes.
  cell_length_factor = L_FACTOR;
  if (USE_HALO_FIELD && (FIND_BUBBLE_ALGORITHM==2)
      && ((BOX_LEN/(float)HII_DIM) < 1)){ // fairly arbitrary length based on 2 runs i did
    cell_length_factor = 1;
  }
  init_ps();
  Fcoll = (float *) malloc(sizeof(float)*HII_TOT_FFT_NUM_PIXELS);
  Fcoll_uni = (float *) malloc(sizeof(float)*HII_TOT_FFT_NUM_PIXELS);
  Fcoll_fs = (float *) malloc(sizeof(float)*HII_TOT_FFT_NUM_PIXELS);
  Fcoll_sfrxs = (float *) malloc(sizeof(float)*HII_TOT_FFT_NUM_PIXELS);
  mass_box = (float *) malloc(sizeof(float)*HII_TOT_FFT_NUM_PIXELS);
  collapsed_mass = (float *) malloc(sizeof(float)*HII_TOT_FFT_NUM_PIXELS);
  star_mass = (float *) malloc(sizeof(float)*HII_TOT_FFT_NUM_PIXELS);
  SFR = (float *) malloc(sizeof(float)*HII_TOT_FFT_NUM_PIXELS);
  stellar_mass = (float *) malloc(sizeof(float)*HII_TOT_FFT_NUM_PIXELS);
  metallicity = (float *) malloc(sizeof(float)*HII_TOT_FFT_NUM_PIXELS);
  Lya_igm = (float *) malloc(sizeof(float)*HII_TOT_FFT_NUM_PIXELS);
  Lya_sf = (float *) malloc(sizeof(double)*HII_TOT_FFT_NUM_PIXELS);
  Halpha_sf = (float *) malloc(sizeof(double)*HII_TOT_FFT_NUM_PIXELS);
  if (INHOMO_RECO) {  init_MHR();}
  if (HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY) {
    init_21cmMC_arrays();
    Mlim_Fstar = Mass_limit_bisection(M_TURN/50., 1e16, ALPHA_STAR, F_STAR10);
    Mlim_Fesc = Mass_limit_bisection(M_TURN/50., 1e16, ALPHA_ESC, F_ESC10);
    ION_EFF_FACTOR = N_GAMMA_UV * F_STAR10 * F_ESC10;
  }
  else
    ION_EFF_FACTOR = N_GAMMA_UV * STELLAR_BARYON_FRAC * ESC_FRAC; // Constant ionizing efficiency parameter.

  // compute average ION_EFF_FACTOR if we are using sources defined in SOURCES.H - RM

  if (USE_GENERAL_SOURCES)
  {
    //ION_EFF_FACTOR = ionEff(REDSHIFT, src);
    //printf("%le\n", ION_EFF_FACTOR);
    ION_EFF_FACTOR = 1;
  }

  // Set the minimum halo mass hosting ionizing source mass.
  // For constant ionizing efficiency parameter M_MIN is set to be M_TURN which is a sharp cut-off.
  // For the new parametrization the number of halos hosting active galaxies (i.e. the duty cycle) is assumed to
  // exponentially decrease below M_TURNOVER Msun, : fduty \propto e^(- M_TURNOVER / M)
  // In this case, we define M_MIN = M_TURN/50.
  if (HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY) {
    M_MIN = M_TURN;
  }
  else {
    M_MIN = M_TURNOVER;
  }

  if(USE_GENERAL_SOURCES) M_MIN = src.minMass(REDSHIFT);

  // Set minimum mass if we are using sources defined in SOURCES.H - RM
  //if (USE_GENERAL_SOURCES)
  //{
  //    if(src.minMassIII(REDSHIFT) > src.minMass(REDSHIFT)
  //       || src.minMassIII(REDSHIFT) < 0) M_MIN = src.minMass(REDSHIFT);
  //    else M_MIN = src.minMassIII(REDSHIFT);
 // }


  // check for WDM
  if (P_CUTOFF && ( M_MIN < M_J_WDM())){
    fprintf(stderr, "The default Jeans mass of %e Msun is smaller than the scale supressed by the effective pressure of WDM.\n", M_MIN);
    M_MIN = M_J_WDM();
    fprintf(stderr, "Setting a new effective Jeans mass from WDM pressure supression of %e Msun\n", M_MIN);
  }
  if (USE_GENERAL_SOURCES) {
    initialiseSplinedSigmaM(M_MIN, 5e16); // set up interpolation table for computing sigma_M
  }
  else {
    initialiseSplinedSigmaM(M_MIN/50.0, 5e16);
  }


  // INITIALIZE THREADS
  if (fftwf_init_threads()==0){
    fprintf(stderr, "find_HII_bubbles: ERROR: problemin itializing fftwf threads\nAborting\n.");
    return -1;
  }
  omp_set_num_threads(num_th);
  fftwf_plan_with_nthreads(num_th);
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);


  // OPEN LOG FILE
  system("mkdir ../Log_files");
  sprintf(filename, "../Log_files/HII_bubble_log_file_%d", getpid());
  LOG = fopen(filename, "w");
  if (!LOG){
    fprintf(stderr, "find_HII_bubbles.c: Error opening log file\n");
  }

  // allocate memory for the neutral fraction box
  xH = (float *) fftwf_malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
  if (!xH){
    strcpy(error_message, "find_HII_bubbles.c: Error allocating memory for xH box\nAborting...\n");
    goto CLEANUP;
  }
  for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){    xH[ct] = 1;  }


  // compute the mean collpased fraction at this redshift
  if (HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY){ // New in v1.4
    // - RM
    //if(USE_GENERAL_SOURCES)
    //{
    //  mean_f_coll_st = FgtrM_st(REDSHIFT, src.minMass(REDSHIFT));
    //}
    //else
    //{
      mean_f_coll_st = FgtrM_st_SFR(REDSHIFT, M_TURN, ALPHA_STAR, ALPHA_ESC, F_STAR10, F_ESC10, Mlim_Fstar, Mlim_Fesc);
    //}
  }
  else {
    mean_f_coll_st = FgtrM_st(REDSHIFT, M_MIN);
  }
  //custom f_coll using all calculations from "else" parts -MG
  mean_f_coll_st_uni = FgtrM_st(REDSHIFT, M_MIN);
  mean_f_coll_st_fs = FgtrM_st_fs(REDSHIFT, M_TURN, ALPHA_STAR, F_STAR10, Mlim_Fstar);
  //fprintf(stderr, "MEAN FCOLL ST FS - - - - - - - - : %e\n", mean_f_coll_st_fs);
  mean_f_coll_st_sfrxs = FgtrM_st_sfrxs(REDSHIFT, M_TURN, ALPHA_STAR, F_STAR10, Mlim_Fstar);
  // - RM
  //X_LUMINOSITY = 3.2e40;
    /**********  CHECK IF WE ARE IN THE DARK AGES ******************************/
    // lets check if we are going to bother with computing the inhmogeneous field at all...
    if ((mean_f_coll_st*ION_EFF_FACTOR < HII_ROUND_ERR)){ // way too small to ionize anything...//New in v1.4
      fprintf(stderr, "The ST mean collapse fraction is %e, which is much smaller than the effective critical collapse fraction of %e\n \
                       I will just declare everything to be neutral\n", mean_f_coll_st, 1./ION_EFF_FACTOR);
      fprintf(LOG, "The ST mean collapse fraction is %e, which is much smaller than the effective critical collapse fraction of %e\n \
                    I will just declare everything to be neutral\n", mean_f_coll_st, 1./ION_EFF_FACTOR);

      if (USE_TS_IN_21CM){ // use the x_e box to set residuals
    if(HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY){ // New in v1.4
        sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_L_X%.1e_alphaX%.1f_f_star10_%06.4f_alpha_star%06.4f_f_esc10_%06.4f_alpha_esc%06.4f_Mturn%.1e_t_star%06.4f_Pop%i_%i_%.0fMpc", REDSHIFT, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, T_AST, Pop, HII_DIM, BOX_LEN);
      ("filename: %s\n",filename);
    }
    else {
      sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_L_X%.1e_alphaX%.1f_Mmin%.1e_zetaIon%.2f_Pop%i_%i_%.0fMpc", REDSHIFT, X_LUMINOSITY, X_RAY_SPEC_INDEX, M_MIN, HII_EFF_FACTOR, Pop, HII_DIM, BOX_LEN);
    }
  if (!(F = fopen(filename, "rb"))){
    fprintf(stderr, "find_HII_bubbles: Unable to open x_e file at %s\nAborting...\n", filename);
    fprintf(LOG, "find_HII_bubbles: Unable to open x_e file at %s\nAborting...\n", filename);
    fclose(LOG); fftwf_free(xH); fftwf_cleanup_threads();
      free(Fcoll);
      free(Fcoll_uni);
      free(Fcoll_fs);
      free(mass_box);
      free(collapsed_mass);
      free(star_mass);
      free(SFR);
      free(stellar_mass);
      free(metallicity);
      free(Lya_igm);
      free(Lya_sf);
      free(Halpha_sf);
    free_ps();  if (HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY) {destroy_21cmMC_arrays();} return -1;
  }
  for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
    if (fread(&xH[ct], sizeof(float), 1, F)!=1){
      strcpy(error_message, "find_HII_bubbles.c: Read error occured while reading xe box.\nAborting...\n");
      goto CLEANUP;
    }
    xH[ct] = 1-xH[ct]; // convert from x_e to xH
    if (xH[ct]<0) xH[ct] = 0; //  should not happen....
    global_xH += xH[ct];
  }
  fclose(F);
  F = NULL;
  global_xH /= (double)HII_TOT_NUM_PIXELS;
      }
      else{
  // find the neutral fraction
  init_heat();
  global_xH = 1 - xion_RECFAST(REDSHIFT, 0);;
  destruct_heat();
  for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
    xH[ct] = global_xH;
  }
      }

      // print out the xH box
      global_xH_m = global_xH;
      switch(FIND_BUBBLE_ALGORITHM){
      case 2:
    if(HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY){
        if (USE_HALO_FIELD)
      sprintf(filename, "../Boxes/xH_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_Mturn%.2e_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, HII_FILTER, MFP, HII_DIM, BOX_LEN);
      else
      sprintf(filename, "../Boxes/xH_nohalos_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_Mturn%.2e_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, HII_FILTER, MFP, HII_DIM, BOX_LEN);
    }
    else{
      if (USE_HALO_FIELD)
        sprintf(filename, "../Boxes/xH_z%06.2f_nf%f_eff%.1f_effPLindex0_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, ION_EFF_FACTOR, HII_FILTER, M_MIN, MFP, HII_DIM, BOX_LEN);
          else
            sprintf(filename, "../Boxes/xH_nohalos_z%06.2f_nf%f_eff%.1f_effPLindex0_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, ION_EFF_FACTOR, HII_FILTER, M_MIN, MFP, HII_DIM, BOX_LEN);
        }
      break;
      default:
    if(HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY){
      if (USE_HALO_FIELD)
      sprintf(filename, "../Boxes/sphere_xH_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_Mturn%.2e_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, HII_FILTER, MFP, HII_DIM, BOX_LEN);
    else
      sprintf(filename, "../Boxes/sphere_xH_nohalos_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_Mturn%.2e_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, HII_FILTER, MFP, HII_DIM, BOX_LEN);
    }
    else{
        if (USE_HALO_FIELD)
          sprintf(filename, "../Boxes/sphere_xH_z%06.2f_nf%f_eff%.1f_effPLindex0_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, ION_EFF_FACTOR, HII_FILTER, M_MIN, MFP, HII_DIM, BOX_LEN);
        else
          sprintf(filename, "../Boxes/sphere_xH_nohalos_z%06.2f_nf%f_eff%.1f_effPLindex0_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, ION_EFF_FACTOR, HII_FILTER, M_MIN, MFP, HII_DIM, BOX_LEN);
    }
      }
      F = fopen(filename, "wb");
      fprintf(LOG, "Neutral fraction is %f\nNow writting xH box at %s\n", global_xH, filename);
      fprintf(stderr, "Neutral fraction is %f\nNow writting xH box at %s\n", global_xH, filename);
      if (mod_fwrite(xH, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
  fprintf(stderr, "find_HII_bubbles.c: Write error occured while writting xH box.\n");
  fprintf(LOG, "find_HII_bubbles.c: Write error occured while writting xH box.\n");
      }
      free_ps(); fclose(F); F = NULL; fclose(LOG); fftwf_free(xH); fftwf_cleanup_threads();
      free(Fcoll);
      free(Fcoll_uni);
      free(Fcoll_sfrxs);
      free(Fcoll_fs);
      free(mass_box);
      free(collapsed_mass);
      free(star_mass);
      free(SFR);
      free(stellar_mass);
      free(metallicity);
      free(Lya_igm);
      free(Lya_sf);
      free(Halpha_sf);
      if (HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY) {destroy_21cmMC_arrays();} return (int) (global_xH * 100);
    }

    /*************   END CHECK TO SEE IF WE ARE STILL IN THE DARK AGES *************/





    // ARE WE INCLUDING PRE-IONIZATION FROM X-RAYS??
    if (USE_TS_IN_21CM){
      xe_unfiltered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
      xe_filtered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
      if (!xe_filtered || !xe_unfiltered){
  strcpy(error_message, "find_HII_bubbles.c: Error allocating memory for xe boxes\nAborting...\n");
  goto CLEANUP;
      }

      // and read-in
    if(HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY){ // New in v1.4
      sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_L_X%.1e_alphaX%.1f_f_star10_%06.4f_alpha_star%06.4f_f_esc10_%06.4f_alpha_esc%06.4f_Mturn%.1e_t_star%06.4f_Pop%i_%i_%.0fMpc", REDSHIFT, X_LUMINOSITY, X_RAY_SPEC_INDEX, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, T_AST, Pop, HII_DIM, BOX_LEN);
    }
    else {
    sprintf(filename, "../Boxes/Ts_evolution/xeneutral_zprime%06.2f_L_X%.1e_alphaX%.1f_Mmin%.1e_zetaIon%.2f_Pop%i_%i_%.0fMpc", REDSHIFT, X_LUMINOSITY, X_RAY_SPEC_INDEX, M_MIN, HII_EFF_FACTOR, Pop, HII_DIM, BOX_LEN);
    }
      if (!(F = fopen(filename, "rb"))){
  strcpy(error_message, "find_HII_bubbles.c: Unable to open x_e file at ");
  strcat(error_message, filename);
  strcat(error_message, "\nAborting...\n");
  goto CLEANUP;
      }
      for (i=0; i<HII_DIM; i++){
  for (j=0; j<HII_DIM; j++){
    for (k=0; k<HII_DIM; k++){
      if (fread((float *)xe_unfiltered + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, F)!=1){
        strcpy(error_message, "find_HII_bubbles.c: Read error occured while reading xe box.\nAborting...\n");
        goto CLEANUP;
      }
    }
  }
      }
      fclose(F);
    F = NULL;
    }


    // ARE WE USING A DISCRETE HALO FIELD (identified in the ICs with find_halos.c and evolved  with update_halos.c)
    if (USE_HALO_FIELD){
      // allocate memory for the smoothed halo field
      M_coll_unfiltered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
      M_coll_filtered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
      if (!M_coll_unfiltered || !M_coll_filtered){
  strcpy(error_message, "find_HII_bubbles.c: Error allocating memory for M_coll boxes\nAborting...\n");
  goto CLEANUP;
      }
      for (ct=0; ct<HII_TOT_FFT_NUM_PIXELS; ct++){    *((float *)M_coll_unfiltered + ct) = 0;  }

      // read in the halo list
      sprintf(filename, "../Output_files/Halo_lists/updated_halos_z%06.2f_%i_%.0fMpc", REDSHIFT, DIM, BOX_LEN);
      F = fopen(filename, "r");
      if (!F){
  strcpy(error_message, "find_HII_bubbles.c: Unable to open halo list file: ");
  strcat(error_message, filename);
  strcat(error_message, "\nAborting...\n");
  goto CLEANUP;
      }
      // now read in all halos above our threshold into the smoothed halo field
      fscanf(F, "%e %f %f %f", &mass, &xf, &yf, &zf);
      while (!feof(F) && (mass>=M_MIN)){
  x = xf*HII_DIM;
  y = yf*HII_DIM;
  z = zf*HII_DIM;
  *((float *)M_coll_unfiltered + HII_R_FFT_INDEX(x, y, z)) += mass;
  fscanf(F, "%e %f %f %f", &mass, &xf, &yf, &zf);
      }
      fclose(F);
    F = NULL;
    } // end of the USE_HALO_FIELD option


    // ALLOCATE AND READ-IN THE EVOLVED DENSITY FIELD
    deltax_unfiltered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    deltax_filtered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
    if (!deltax_unfiltered || !deltax_filtered){
      strcpy(error_message, "find_HII_bubbles.c: Error allocating memory for deltax boxes\nAborting...\n");
      goto CLEANUP;
    }
    fprintf(stderr, "Reading in deltax box\n");
    fprintf(LOG, "Reading in deltax box\n");

    sprintf(filename, "../Boxes/updated_smoothed_deltax_z%06.2f_%i_%.0fMpc", REDSHIFT, HII_DIM, BOX_LEN);
    F = fopen(filename, "rb");
    if (!F){
      strcpy(error_message, "find_HII_bubbles.c: Unable to open file: ");
      strcat(error_message, filename);
      strcat(error_message, "\nAborting...\n");
      goto CLEANUP;
    }
    for (i=0; i<HII_DIM; i++){
      for (j=0; j<HII_DIM; j++){
  for (k=0; k<HII_DIM; k++){
    if (fread((float *)deltax_unfiltered + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, F)!=1){
      strcpy(error_message, "find_HII_bubbles.c: Read error occured while reading deltax box.\n");
      goto CLEANUP;
    }
  }
      }
    }
    fclose(F);
    F = NULL;

    // ALLOCATE AND INITIALIZE ADDITIONAL BOXES NEEDED TO KEEP TRACK OF RECOMBINATIONS (Sobacchi & Mesinger 2014; NEW IN v1.3)
    if (INHOMO_RECO){ //  flag in ANAL_PARAMS.H to determine whether to compute recombinations or not
      z_re = (float *) fftwf_malloc(sizeof(float)*HII_TOT_NUM_PIXELS); // the redshift at which the cell is ionized
      Gamma12 = (float *) fftwf_malloc(sizeof(float)*HII_TOT_NUM_PIXELS);  // stores the ionizing backgroud
      N_rec_unfiltered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS); // cumulative number of recombinations
      N_rec_filtered = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
      if (!z_re || !N_rec_filtered || !N_rec_unfiltered || !Gamma12){
  strcpy(error_message, "find_HII_bubbles.c: Error allocating memory for recombination boxes boxes\nAborting...\n");
  goto CLEANUP;
      }
      sprintf(filename, "../Boxes/z_first_ionization_z%06.2f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc",  PREV_REDSHIFT, HII_FILTER, MFP, HII_DIM, BOX_LEN);
      if (F=fopen(filename, "rb")){  // this is the first call for this run, i.e. at the highest redshift
  //check if some read error occurs
  if (mod_fread(z_re, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
    strcpy(error_message, "find_HII_bubbles.c: Read error occured while reading z_re box!\n");
    goto CLEANUP;
  }
      }
      else{ // this is the first call (highest redshift)
  for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++)
    z_re[ct] = -1.0;
      }
      sprintf(filename, "../Boxes/Nrec_z%06.2f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", PREV_REDSHIFT, HII_FILTER, MFP, HII_DIM, BOX_LEN);
      if (F=fopen(filename, "rb")){ // we had prvious boxes
  //check if some read error occurs
  for (i=0; i<HII_DIM; i++){
    for (j=0; j<HII_DIM; j++){
      for (k=0; k<HII_DIM; k++){
        if (fread((float *)N_rec_unfiltered + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, F)!=1){
    strcpy(error_message, "find_HII_bubbles.c: Read error occured while reading N_rec box!\n");
    goto CLEANUP;
        }
      }
    }
  }
      }
      else{
  fprintf(stderr, "Earliest redshift call, initializing Nrec to zeros\n");
  // initialize N_rec
  for (i=0; i<HII_DIM; i++){
    for (j=0; j<HII_DIM; j++){
      for (k=0; k<HII_DIM; k++){
        *((float *)N_rec_unfiltered + HII_R_FFT_INDEX(i,j,k)) = 0.0;
      }
    }
  }
      }

    } //  end if INHOMO_RECO

    // do the fft to get the k-space M_coll field and deltax field
    fprintf(LOG, "begin initial ffts, clock=%06.2f\n", (double)clock()/CLOCKS_PER_SEC);
    fflush(LOG);
    if (USE_HALO_FIELD){
      plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)M_coll_unfiltered, (fftwf_complex *)M_coll_unfiltered, FFTW_ESTIMATE);
      fftwf_execute(plan);
    }
    if (USE_TS_IN_21CM){
      plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)xe_unfiltered, (fftwf_complex *)xe_unfiltered, FFTW_ESTIMATE);
      fftwf_execute(plan);
    }
    if (INHOMO_RECO){
      plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)N_rec_unfiltered, (fftwf_complex *)N_rec_unfiltered, FFTW_ESTIMATE);
      fftwf_execute(plan);
    }
    plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)deltax_unfiltered, (fftwf_complex *)deltax_unfiltered, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    fftwf_cleanup();
    // remember to add the factor of VOLUME/TOT_NUM_PIXELS when converting from
    //  real space to k-space
    // Note: we will leave off factor of VOLUME, in anticipation of the inverse FFT below
    for (ct=0; ct<HII_KSPACE_NUM_PIXELS; ct++){
      if (USE_HALO_FIELD){  M_coll_unfiltered[ct] /= (double)HII_TOT_NUM_PIXELS;}
      if (USE_TS_IN_21CM){ xe_unfiltered[ct] /= (double)HII_TOT_NUM_PIXELS;}
      if (INHOMO_RECO){ N_rec_unfiltered[ct] /= (double)HII_TOT_NUM_PIXELS; }
      deltax_unfiltered[ct] /= (HII_TOT_NUM_PIXELS+0.0);
    }
    fprintf(LOG, "end initial ffts, clock=%06.2f\n", (double)clock()/CLOCKS_PER_SEC);
    fflush(LOG);


     /*********************** DESTRUCTOR ************************************************/
 CLEANUP: if (*error_message != '\0'){
      fprintf(stderr, error_message);
      fprintf(LOG, error_message);
      fftwf_free(xH);
      fftwf_free(z_re);
      fftwf_free(Gamma12);
      fclose(LOG);
      fftwf_free(deltax_unfiltered);
      fftwf_free(deltax_filtered);
      fftwf_free(M_coll_unfiltered);
      fftwf_free(M_coll_filtered);
      fftwf_cleanup_threads();
      free_ps();
      if (INHOMO_RECO) { free_MHR();}
      fftwf_free(xe_filtered);
      fftwf_free(xe_unfiltered);
      fftwf_free(N_rec_unfiltered);
      fftwf_free(N_rec_filtered);
      free(Fcoll);
      free(Fcoll_uni);
      free(Fcoll_fs);
      free(Fcoll_sfrxs);
      free(mass_box);
      free(collapsed_mass);
      free(star_mass);
      free(SFR);
      free(stellar_mass);
      free(metallicity);
      free(Lya_igm);
      free(Lya_sf);
      free(Halpha_sf);
      if (HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY) {destroy_21cmMC_arrays();}

      return -1;
    }



    /*************************************************************************************/
    /*****************  END OF INITIALIZATION ********************************************/
    /*************************************************************************************/


    /*************************************************************************************/
    /***************** LOOP THROUGH THE FILTER RADII (in Mpc)  ***************************/
    /*************************************************************************************/
    // set the max radius we will use, making sure we are always sampling the same values of radius
    // (this avoids aliasing differences w redshift)
    R=fmax(R_BUBBLE_MIN, (cell_length_factor*BOX_LEN/(float)HII_DIM));
    while (R < fmin(MFP, L_FACTOR*BOX_LEN)) { R*= DELTA_R_HII_FACTOR;}
    R /= DELTA_R_HII_FACTOR;
    LAST_FILTER_STEP = 0;
    while (!LAST_FILTER_STEP && (M_MIN < RtoM(R)) ){

      if ( ((R/DELTA_R_HII_FACTOR) <= (cell_length_factor*BOX_LEN/(float)HII_DIM)) || ((R/DELTA_R_HII_FACTOR) <= R_BUBBLE_MIN) ){
  LAST_FILTER_STEP = 1;
  R = fmax(cell_length_factor*BOX_LEN/(double)HII_DIM, R_BUBBLE_MIN);
      }

      fprintf(LOG, "before memcpy, clock=%06.2f\n", (double)clock()/CLOCKS_PER_SEC);
      fflush(LOG);
      if (USE_HALO_FIELD){    memcpy(M_coll_filtered, M_coll_unfiltered, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);}
      if (USE_TS_IN_21CM){    memcpy(xe_filtered, xe_unfiltered, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);}
      if (INHOMO_RECO){  memcpy(N_rec_filtered, N_rec_unfiltered, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);}
      memcpy(deltax_filtered, deltax_unfiltered, sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);


      // if this scale is not the size of the cells, we need to filter the fields
      if (!LAST_FILTER_STEP || (R > cell_length_factor*BOX_LEN/(double)HII_DIM) ){
  fprintf(LOG, "begin filter, clock=%06.2f\n", (double)clock()/CLOCKS_PER_SEC);
        fflush(LOG);
        if (USE_HALO_FIELD) { HII_filter(M_coll_filtered, HII_FILTER, R);}
        if (USE_TS_IN_21CM) { HII_filter(xe_filtered, HII_FILTER, R);}
  if (INHOMO_RECO) {  HII_filter(N_rec_filtered, HII_FILTER, R);}
        HII_filter(deltax_filtered, HII_FILTER, R);
        fprintf(LOG, "end filter, clock=%06.2f\n", (double)clock()/CLOCKS_PER_SEC);
        fflush(LOG);
      }

      // do the FFT to get back to real space
      fprintf(LOG, "begin fft with R=%f, clock=%06.2f\n", R, (double)clock()/CLOCKS_PER_SEC);
      fflush(LOG);
      if (USE_HALO_FIELD) {
  plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)M_coll_filtered, (float *)M_coll_filtered, FFTW_ESTIMATE);
  fftwf_execute(plan);
      }
      if (USE_TS_IN_21CM) {
  plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)xe_filtered, (float *)xe_filtered, FFTW_ESTIMATE);
  fftwf_execute(plan);
      }
      if (INHOMO_RECO){
  plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)N_rec_filtered, (float *)N_rec_filtered, FFTW_ESTIMATE);
  fftwf_execute(plan);
      }
      plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)deltax_filtered, (float *)deltax_filtered, FFTW_ESTIMATE);
      fftwf_execute(plan);
      fftwf_destroy_plan(plan);
      fftwf_cleanup();
      fprintf(LOG, "end fft with R=%f, clock=%06.2f\n", R, (double)clock()/CLOCKS_PER_SEC);
      fflush(LOG);


      // perform sanity checks to account for aliasing effects
      for (x=0; x<HII_DIM; x++){
  for (y=0; y<HII_DIM; y++){
    for (z=0; z<HII_DIM; z++){

      // delta cannot be less than -1
      *((float *)deltax_filtered + HII_R_FFT_INDEX(x,y,z)) =
        FMAX(*((float *)deltax_filtered + HII_R_FFT_INDEX(x,y,z)) , -1+FRACT_FLOAT_ERR);

      // <N_rec> cannot be less than zero
      if (INHOMO_RECO){
        *((float *)N_rec_filtered + HII_R_FFT_INDEX(x,y,z)) =
    FMAX(*((float *)N_rec_filtered + HII_R_FFT_INDEX(x,y,z)) , 0.0);
      }

      // collapsed mass cannot be less than zero
      if (USE_HALO_FIELD){
        *((float *)M_coll_filtered + HII_R_FFT_INDEX(x,y,z)) =
    FMAX(*((float *)M_coll_filtered + HII_R_FFT_INDEX(x,y,z)) , 0.0);
      }

      // x_e has to be between zero and unity
      if (USE_TS_IN_21CM){
        *((float *)xe_filtered + HII_R_FFT_INDEX(x,y,z)) = FMAX(*((float *)xe_filtered + HII_R_FFT_INDEX(x,y,z)) , 0);
        *((float *)xe_filtered + HII_R_FFT_INDEX(x,y,z)) = FMIN(*((float *)xe_filtered + HII_R_FFT_INDEX(x,y,z)) , 0.999);
      }

    }
  }
      } //  end sanity check in box


      // normalize the analytic collapse fractions if we are operating on the density field
      f_coll = 0.0;
      f_coll_uni = 0.0;
      f_coll_fs = 0.0;
      f_coll_sfrxs = 0.0; 
      massofscaleR = RtoM(R);
      if (!USE_HALO_FIELD){
  fprintf(LOG, "begin f_coll normalization, clock=%06.2f\n", (double)clock()/CLOCKS_PER_SEC);
  fflush(LOG);
  temparg =  2*(pow(sigma_z0(M_MIN), 2) - pow(sigma_z0(massofscaleR), 2) );
  erfc_denom = sqrt(temparg);
  if(HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY) { // New in v1.4
    initialiseGL_FcollSFR(NGL_SFR, M_MIN/50.0, massofscaleR, REDSHIFT);
    initialiseFcollSFR_spline(REDSHIFT, massofscaleR,M_TURN,ALPHA_STAR,ALPHA_ESC,F_STAR10,F_ESC10,Mlim_Fstar,Mlim_Fesc);
    initialiseFcollSFR_spline_fs(REDSHIFT, massofscaleR, M_TURN, ALPHA_STAR, F_STAR10, Mlim_Fstar);
    initialiseFcollSFR_spline_sfrxs(REDSHIFT, massofscaleR, M_TURN, ALPHA_STAR, F_STAR10, Mlim_Fstar); // - LMR
    // - RM
    /*
    if(USE_GENERAL_SOURCES)
    {
      //initialiseGL_zetaSFR(NGL_SFR, M_MIN/50.0, massofscaleR);
      initialisezetaSFR_spline(REDSHIFT, massofscaleR,M_TURN,ALPHA_STAR,ALPHA_ESC,F_STAR10,F_ESC10,Mlim_Fstar,Mlim_Fesc);
    }*/
  }
  for (x=0; x<HII_DIM; x++){
    for (y=0; y<HII_DIM; y++){
      for (z=0; z<HII_DIM; z++){
        if (USE_HALO_FIELD){
    Splined_Fcoll = *((float *)M_coll_filtered + HII_R_FFT_INDEX(x,y,z)) / (massofscaleR*density_over_mean);
    Splined_Fcoll *= (4/3.0)*PI*pow(R,3) / pixel_volume;
    Splined_Fcoll_uni = Splined_Fcoll;
    Splined_Fcoll_fs = Splined_Fcoll;
    // - RM
    /*
    if(USE_GENERAL_SOURCES)
    {
      Splined_zeta = *((float *)M_coll_filtered + HII_R_FFT_INDEX(x,y,z)) / (massofscaleR*density_over_mean);
      Splined_zeta *= (4/3.0)*PI*pow(R,3) / pixel_volume;
    }*/
        }

        else{
          // - RM
          /*
      if(USE_GENERAL_SOURCES) {
        density_over_mean = 1.0 + *((float *)deltax_filtered + HII_R_FFT_INDEX(x,y,z));
        zetaSpline_SFR(density_over_mean - 1,&(Splined_zeta));
      }*/
          density_over_mean = 1.0 + *((float *)deltax_filtered + HII_R_FFT_INDEX(x,y,z));
          //if ( (density_over_mean - 1) < Deltac){ // we are not resolving collapsed structures
            //make new Splined_Fcoll_uni that mimics else statement
            //erfc_num = (Deltac - (density_over_mean-1)) /  growth_factor;
            //Splined_Fcoll_uni = splined_erfc(erfc_num/erfc_denom);
            //if (HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY) { // New in v1.4
              //FcollSpline_SFR(density_over_mean - 1,&(Splined_Fcoll));
              //FcollSpline_SFR_fs(density_over_mean - 1,&(Splined_Fcoll_fs));
            //}
            //else{ // we can assume the classic constant ionizing luminosity to halo mass ratio
              //Splined_Fcoll = Splined_Fcoll_uni;
              //Splined_Fcoll_fs = Splined_Fcoll_uni;
            //}
          //}
          //else { // the entrire cell belongs to a collpased halo...  this is rare...
            //Splined_Fcoll =  1.0;
            //Splined_Fcoll_uni = 1.0;
            //Splined_Fcoll_fs = 1.0;
          //}
          // gonna avoid equating the collapsed fractions to 1 and make calculate for sure in all cases here - LMR
            erfc_num = (Deltac - (density_over_mean-1)) /  growth_factor;
            Splined_Fcoll_uni = splined_erfc(erfc_num/erfc_denom);
            FcollSpline_SFR(density_over_mean - 1,&(Splined_Fcoll));
            FcollSpline_SFR_fs(density_over_mean - 1,&(Splined_Fcoll_fs));
            FcollSpline_SFR_sfrxs(density_over_mean - 1,&(Splined_Fcoll_sfrxs));
            //if (HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY) { // New in v1.4
            //  FcollSpline_SFR(density_over_mean - 1,&(Splined_Fcoll));
            //  FcollSpline_SFR_fs(density_over_mean - 1,&(Splined_Fcoll_fs));
              //FcollSpline_SFR_sfrxs(density_over_mean - 1,&(Splined_Fcoll_sfrxs));
            //}
            //else{ // we can assume the classic constant ionizing luminosity to halo mass ratio
              //Splined_Fcoll = Splined_Fcoll_uni;
              //Splined_Fcoll_fs = Splined_Fcoll_uni;
              //}
        }
        // save the value of the collapsed fraction into the Fcoll array
        Fcoll[HII_R_FFT_INDEX(x,y,z)] = Splined_Fcoll;
          //load separate fcoll_uni box assuming 'else' case for everything
          Fcoll_uni[HII_R_FFT_INDEX(x,y,z)] = Splined_Fcoll_uni;
          Fcoll_fs[HII_R_FFT_INDEX(x,y,z)] = Splined_Fcoll_fs;
          Fcoll_sfrxs[HII_R_FFT_INDEX(x,y,z)] = Splined_Fcoll_sfrxs;
        f_coll += Splined_Fcoll;
          //add to separate f_coll_uni using 'else' case
          f_coll_uni += Splined_Fcoll_uni;
          f_coll_fs += Splined_Fcoll_fs;
          f_coll_sfrxs += Splined_Fcoll_sfrxs;


        // - RM
        /*
        if(USE_GENERAL_SOURCES)
        {
          Fcoll[HII_R_FFT_INDEX(x,y,z)] = Splined_Fcoll * Splined_zeta;
          //if(x == 0 && y == 0 && z == 0) printf("Fcoll, zeta = %le %le\n", Splined_Fcoll, Splined_zeta);
          f_coll += Splined_Fcoll;
        }  */
      }
    }
  } //  end loop through Fcoll box

  f_coll /= (double) HII_TOT_NUM_PIXELS; // ave PS fcoll for this filter scale
    f_coll_uni /= (double) HII_TOT_NUM_PIXELS; //same for uni one -MG
    f_coll_fs /= (double) HII_TOT_NUM_PIXELS;
    f_coll_sfrxs /= (double) HII_TOT_NUM_PIXELS;
  ST_over_PS = mean_f_coll_st/f_coll; // normalization ratio used to adjust the PS conditional collapsed fraction
    ST_over_PS_uni = mean_f_coll_st_uni / f_coll_uni; //same for uni one -MG
    ST_over_PS_fs = mean_f_coll_st_fs / f_coll_fs;
    ST_over_PS_sfrxs = mean_f_coll_st_sfrxs / f_coll_sfrxs;
  fprintf(LOG, "end f_coll normalization if, clock=%06.2f\n", (double)clock()/CLOCKS_PER_SEC);
  fflush(LOG);
      } // end ST fcoll normalization

      //     fprintf(stderr, "Last filter %i, R_filter=%f, fcoll=%f, ST_over_PS=%f, mean normalized fcoll=%f\n", LAST_FILTER_STEP, R, f_coll, ST_over_PS, f_coll*ST_over_PS);


      /****************************************************************************/
      /************  MAIN LOOP THROUGH THE BOX FOR THIS FILTER SCALE **************/
      /****************************************************************************/
      fprintf(LOG, "Start of the main loop through the box for this filter scale, clock=%06.2f\n", (double)clock()/CLOCKS_PER_SEC);
      fflush(LOG);
      // now lets scroll through the filtered box
      rec = ave_xHI_xrays = ave_den = ave_fcoll = std_xrays = 0;
      ion_ct=0;
      xHI_from_xrays = 1;
      Gamma_R_prefactor = pow(1+REDSHIFT, 2) * (R*CMperMPC) * SIGMA_HI * ALPHA_UVB / (ALPHA_UVB+2.75) * N_b0 * ION_EFF_FACTOR / 1.0e-12;

      for (x=0; x<HII_DIM; x++){
  for (y=0; y<HII_DIM; y++){
    for (z=0; z<HII_DIM; z++){

      density_over_mean = 1.0 + *((float *)deltax_filtered + HII_R_FFT_INDEX(x,y,z));

      /*
      if (LAST_FILTER_STEP && (density_over_mean > 3)){
        fprintf(stderr, "%g at index: %i, %i, %i\n", density_over_mean, x,y,z);
        strcpy(error_message, "shs\n");
        goto CLEANUP;
      }
      */

      f_coll = ST_over_PS * Fcoll[HII_R_FFT_INDEX(x,y,z)];
      //make custom f_coll_uni using custom st_over_ps
      f_coll_uni = ST_over_PS_uni * Fcoll_uni[HII_R_FFT_INDEX(x,y,z)];
      //f_coll_sfrxs = ST_over_PS_sfrxs * Fcoll_sfrxs[HII_R_FFT_INDEX(x,y,z)];

      if (INHOMO_RECO){
        dfcolldt = f_coll / t_ast;
        Gamma_R = Gamma_R_prefactor * dfcolldt;
        rec = (*((float *)N_rec_filtered + HII_R_FFT_INDEX(x,y,z))); // number of recombinations per mean baryon
        rec /= density_over_mean; // number of recombinations per baryon inside <R>
      }

      // adjust the denominator of the collapse fraction for the residual electron fraction in the neutral medium
      if (USE_TS_IN_21CM){
        xHI_from_xrays =  1 - (*((float *)xe_filtered + HII_R_FFT_INDEX(x,y,z)));
      }

      // check if fully ionized!
      //printf("F_COLL ZETA = %le\n", f_coll);
      //if ( (f_coll*ION_EFF_FACTOR > xHI_from_xrays*(1.0+rec) && !(USE_GENERAL_SOURCES)) || (f_coll > xHI_from_xrays*(1.0 + rec) && (USE_GENERAL_SOURCES)) ){ //IONIZED!! //New in v1.4
      if ( (f_coll*ION_EFF_FACTOR > xHI_from_xrays*(1.0+rec))) { //IONIZED!! //New in v1.4
        // if this is the first crossing of the ionization barrier for this cell (largest R), record the gamma
        // this assumes photon-starved growth of HII regions...  breaks down post EoR
        if (INHOMO_RECO && (xH[HII_R_INDEX(x, y, z)] > FRACT_FLOAT_ERR) ){
    Gamma12[HII_R_INDEX(x, y, z)] = Gamma_R;
    aveR += R;
    Rct++;
        }

        // keep track of the first time this cell is ionized (earliest time)
        if (INHOMO_RECO && (z_re[HII_R_INDEX(x, y, z)] < 0)){
    z_re[HII_R_INDEX(x, y, z)] = REDSHIFT;
        }


        // FLAG CELL(S) AS IONIZED
        if (FIND_BUBBLE_ALGORITHM == 2) // center method
    xH[HII_R_INDEX(x, y, z)] = 0;
        else if (FIND_BUBBLE_ALGORITHM == 1) // sphere method
    update_in_sphere(xH, HII_DIM, R/BOX_LEN, x/(HII_DIM+0.0), y/(HII_DIM+0.0), z/(HII_DIM+0.0));
        else{
    fprintf(stderr, "Incorrect choice of find bubble algorithm set in ANAL_PARAMS.H.\nSetting center method...");
    fprintf(LOG, "Incorrect choice of find bubble algorithm set in ANAL_PARAMS.H.\nSetting center method...");
    xH[HII_R_INDEX(x, y, z)] = 0;
        }

      } // end ionized

      // If not fully ionized, then assign partial ionizations
      else if (LAST_FILTER_STEP && (xH[HII_R_INDEX(x, y, z)] > TINY)){

        if (!USE_HALO_FIELD){
    ave_N_min_cell = f_coll * pixel_mass * density_over_mean / M_MIN; // ave # of M_MIN halos in cell
    if (ave_N_min_cell < N_POISSON){ // add poissonian fluctuations to the nalo number
      N_min_cell = (int) gsl_ran_poisson(r, ave_N_min_cell);
      f_coll = N_min_cell * M_MIN / (pixel_mass*density_over_mean);
    }
        }

        res_xH = xHI_from_xrays - f_coll * ION_EFF_FACTOR;
        //if(USE_GENERAL_SOURCES) res_xH = xHI_from_xrays - f_coll;
        // and make sure fraction doesn't blow up for underdense pixels
        if (res_xH < 0)
    res_xH = 0;
        else if (res_xH > 1)
    res_xH = 1;

        xH[HII_R_INDEX(x, y, z)] = res_xH;
      } // end partial ionizations at last filtering step

      /*
      if ((x==0) && (y==19) && (z==107)){ //DEBUGGING
        fprintf(stderr, "R=%g Mpc, <Delta>_R=%g, <fcoll>_R=%g, <rec>_V=%g, xH=%g, Gamma12=%g, z_re=%g\n",
          R, density_over_mean, f_coll, rec, xH[HII_R_INDEX(x, y, z)], Gamma12[HII_R_INDEX(x, y, z)], z_re[HII_R_INDEX(x, y, z)]);
        fprintf(LOG, "R=%g Mpc, <Delta>_R=%g, <fcoll>_R=%g, <rec>_V=%g, xH=%g, Gamma12=%g, z_re=%g\n",
          R, density_over_mean, f_coll, rec, xH[HII_R_INDEX(x, y, z)], Gamma12[HII_R_INDEX(x, y, z)], z_re[HII_R_INDEX(x, y, z)]);
        fflush(NULL);
      }
      */

    } // z
  } // y
      } // x




      R /= DELTA_R_HII_FACTOR;
    } // END OF LOOP THROUGH FILTER RADII


    int compute_seds = 1;
    //read in last metallicity box to add to (if applicable) -MG
    sprintf(filename, "../Boxes/metallicity_z%06.2f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", PREV_REDSHIFT, HII_FILTER, MFP, HII_DIM, BOX_LEN);
    F = fopen(filename, "rb");
    if (!(F = fopen(filename, "rb"))){
      fprintf(stderr, "Error opening file %s for reading.\n", filename);
      fprintf(LOG, "Error opening file %s for reading.\n", filename);
      for (ct=0; ct<HII_TOT_FFT_NUM_PIXELS;ct++)
        metallicity[ct] = 0;
      //compute_seds = 0;
    }
    else {
      fprintf(stderr, "Loading metallicity box with previous redshift values at %s\n", filename);
      fprintf(LOG, "Loading metallicity box with previous redshift values at %s\n", filename);
      for (x=0; x<HII_DIM; x++) {
        for (y=0; y<HII_DIM; y++) {
          for (z=0; z<HII_DIM; z++) {
            if (fread((float *)metallicity + HII_R_FFT_INDEX(x,y,z), sizeof(float), 1, F)!=1){
              strcpy(error_message, "find_HII_bubbles.c: Read error occured while reading metallicity box!\n");
              compute_seds = 0;
              goto CLEANUP;
            }
          }
        }
      }
      fclose(F);
    }

    //read in last stellar mass box to add to (if applicable) -MG
    sprintf(filename, "../Boxes/stellarmass_z%06.2f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", PREV_REDSHIFT, HII_FILTER, MFP, HII_DIM, BOX_LEN);
    F = fopen(filename, "rb");
    if (!(F = fopen(filename, "rb"))){
      fprintf(stderr, "Error opening file %s for reading.\n", filename);
      fprintf(LOG, "Error opening file %s for reading.\n", filename);
      for (ct=0; ct<HII_TOT_FFT_NUM_PIXELS;ct++)
        stellar_mass[ct] = 0;
    }
    else {
      fprintf(stderr, "Loading stellar mass box with previous redshift values at %s\n", filename);
      fprintf(LOG, "Loading stellar mass box with previous redshift values at %s\n", filename);
      for (x=0; x<HII_DIM; x++) {
        for (y=0; y<HII_DIM; y++) {
          for (z=0; z<HII_DIM; z++) {
            if (fread((float *)stellar_mass + HII_R_FFT_INDEX(x,y,z), sizeof(float), 1, F)!=1){
              strcpy(error_message, "find_HII_bubbles.c: Read error occured while reading stellar mass box!\n");
              goto CLEANUP;
            }
          }
        }
      }
      fclose(F);
    }

    val = timesince(REDSHIFT, PREV_REDSHIFT) / (60*60*24*365.25); //in yrs, calculating here to hopefully speed things up -MG
    //output box now (for sure uses last fcoll calculation) -MG
    for (x=0; x<HII_DIM; x++){
      for (y=0; y<HII_DIM; y++){
        for (z=0; z<HII_DIM; z++){
          //calculate value of mass, collapsed mass, stellar mass, SFR, and metallicity -MG
          ct = HII_R_FFT_INDEX(x,y,z);
          //density_over_mean = 1.0 + *((float *)deltax_filtered + ct);
          //mass_box[ct] = density_over_mean * OMm * RHOcrit * pixel_volume;
          //collapsed_mass[ct] = ST_over_PS_uni * Fcoll_uni[ct] * mass_box[ct];
          //star_mass[ct] = F_b * ST_over_PS_fs * Fcoll_fs[ct] * mass_box[ct]; //F_b times collapsed mass with Fs taken into account
          //SFR[ct] = star_mass[ct] * (60*60*24*365.25) / t_ast; // adjust for t_ast being in seconds -> years
          SFR[ct] = Fcoll_sfrxs[ct] * ST_over_PS_sfrxs *  OMm * RHOcrit;
          stellar_mass[ct] += SFR[ct] * val;
          //metallicity[ct] += Y_z * SFR[ct] * val / (mean_rho_comoving(REDSHIFT) * pixel_volume * Fcoll_uni[ct] * ST_over_PS_uni);
          //n_Hii = (1 - xH[HII_R_INDEX(x,y,z)]) * N_b0 * pow(1+REDSHIFT, 3); //* density_over_mean; //equation 7 in overleaf
          //Lya_emissivity = f_rec_lya * (double) 4.2e-13 * pow(1+REDSHIFT, 3) * (1 + He_abundance / (2 * H_abundance)) * Ly_alpha_E * pow(n_Hii, 2);
          //Lya_igm[ct] = (C * Lya_emissivity) / (4 * PI * Ly_alpha_HZ * hubble(REDSHIFT));
        }
      }
    }
    //average stellar mass box and print
    double ave_sfr = 0.0;
    /*double mass_avg = 0.0;
    double overdens_avg = 0.0;
    double collmass_avg = 0.0;
    double stellar_mass_avg = 0.0;
    double fcoll_avg = 0.0;
    double fcollfs_avg = 0.0;*/
    for (x=0; x<HII_DIM; x++) {
      for (y=0; y<HII_DIM; y++) {
        for (z=0; z<HII_DIM; z++) {
            ct = HII_R_FFT_INDEX(x,y,z);
            ave_sfr += (double) stellar_mass[ct];
            //mass_avg += mass_box[ct];
            //overdens_avg += 1.0 + *((float *)deltax_filtered + ct);
            //collmass_avg += collapsed_mass[ct];
            //stellar_mass_avg += (double) star_mass[ct];
            //fcoll_avg += (double) Fcoll_uni[ct];
            //fcollfs_avg += (double) Fcoll_fs[ct];
        }
      }
    }
    //ave_sfr /= (BOX_LEN * BOX_LEN * BOX_LEN);
    ave_sfr /= (HII_DIM * HII_DIM * HII_DIM);
    //collmass_avg /= (BOX_LEN * BOX_LEN * BOX_LEN);
    //stellar_mass_avg /= (BOX_LEN * BOX_LEN * BOX_LEN);
    //fcoll_avg /= (HII_DIM * HII_DIM * HII_DIM);
    //fcollfs_avg /= (HII_DIM * HII_DIM * HII_DIM);
    fprintf(stderr, "SMD at z = %f - - - - - - - - : %f\n", REDSHIFT, ave_sfr);
    //fprintf(stderr, "Overdensity avg: %f, Mass density: %f, collmass density: %f\n", overdens_avg, mass_avg, collmass_avg);
    //fprintf(stderr, "Stellar mass density: %f, fcoll: %f, fcoll_fs: %f\n", stellar_mass_avg, fcoll_avg, fcollfs_avg);
    //fprintf(stderr, "ST_over_PS_uni: %f, ST_over_PS_fs: %f\n", ST_over_PS_uni, ST_over_PS_fs);*/

    //if custom SEDs are enabled, use metallicity and starmass boxes to calculate other emission line radiation
    if (USE_CUSTOM_SEDS && compute_seds) {
      if (init_pop2_tables(Lya_table_pop2, Halpha_table_pop2) == 0 && init_pop3_table(table_pop3) == 0) {
        //fprintf(stderr, "Calculating line emissivities...\n");
        //fprintf(stderr, "Lya pop3 luminosity test: %f\n", table_pop3[1][5] * table_pop3[1][7] / table_pop3[1][2]);
        //float avg = 0;
        //float avg_lum = 0;
        for (x=0; x<HII_DIM; x++){
          for (y=0; y<HII_DIM; y++){
            for (z=0; z<HII_DIM; z++){
              //use lookup table to calculate intensity of Lya_sf and Halpha_sf
              ct = HII_R_FFT_INDEX(x,y,z);
              metal_val = (float) log10(metallicity[ct] / 0.014); //units of log(z/z_sun);
              if (metal_val > 0.0) { //for now, limit to 0
                metal_val = 0.0;
              }
              //Halpha_lum = lookup(Halpha_table, metal_val, 3);
              Halpha_lum = get_luminosity(Halpha_table_pop2, table_pop3, metal_val, 3, src.fPopIII(REDSHIFT), 3);
              Halpha_sf[ct] = (C * Halpha_lum * Lsun * stellar_mass[ct] / CMperMPC) / (4 * PI * H_alpha_HZ * hubble(REDSHIFT) * pixel_volume * pow(CMperMPC, 2));
              //Lya_lum = lookup(Lya_table, metal_val, 3);
              Lya_lum = get_luminosity(Lya_table_pop2, table_pop3, metal_val, 3, src.fPopIII(REDSHIFT), 5);
              Lya_sf[ct] = (C * Lya_lum * Lsun * stellar_mass[ct] / CMperMPC) / (4 * PI * Ly_alpha_HZ * hubble(REDSHIFT) * pixel_volume * pow(CMperMPC, 2));
              //avg += Lya_sf[ct];
              //avg_lum += Lya_lum;
              if(Halpha_lum <= 0.0 || Lya_lum <= 0.0) {
                fprintf(stderr, "Impossible luminosity at Metallicity = %f, Halpha Luminosity = %f, Lya Luminosity = %f\n", metal_val, Halpha_lum, Lya_lum);
                fprintf(LOG, "Impossible luminosity at Metallicity = %f, Halpha Luminosity = %f, Lya Luminosity = %f\n", metal_val, Halpha_lum, Lya_lum);
              }
            }
          }
        }
        //fprintf(stderr, "avg lya lum: %f, avg lyasf: %f\n", avg_lum/(HII_DIM*HII_DIM*HII_DIM), avg/(HII_DIM*HII_DIM*HII_DIM));
      }
    }

      // find the neutral fraction
    global_xH = 0;
    for (ct=0; ct<HII_TOT_NUM_PIXELS; ct++){
      global_xH += xH[ct];
    }
    global_xH /= (float)HII_TOT_NUM_PIXELS;


    // update the N_rec field
    if (INHOMO_RECO){
      //fft to get the real N_rec  and delta fields
      plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)N_rec_unfiltered, (float *)N_rec_unfiltered, FFTW_ESTIMATE);
      fftwf_execute(plan);
      plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)deltax_unfiltered, (float *)deltax_unfiltered, FFTW_ESTIMATE);
      fftwf_execute(plan);
      fftwf_destroy_plan(plan);
      fftwf_cleanup();
      for (x=0; x<HII_DIM; x++){
  for (y=0; y<HII_DIM; y++){
    for (z=0; z<HII_DIM; z++){

      density_over_mean = 1.0 + (*((float *)deltax_unfiltered + HII_R_FFT_INDEX(x,y,z)));
      z_eff = (1+REDSHIFT) * pow(density_over_mean, 1.0/3.0) - 1;
      dNrec = splined_recombination_rate(z_eff, Gamma12[HII_R_INDEX(x,y,z)]) * fabs_dtdz * ZSTEP * (1 - xH[HII_R_INDEX(x,y,z)]);
      *((float *)N_rec_unfiltered + HII_R_FFT_INDEX(x,y,z)) += dNrec;

      /*
      if ((x==0) && (y==19) && (z==107)){ //DEBUGGING
        fprintf(stderr, "Delta=%g, z_eff=%g, dNrec=%g, Nrec=%g\n\n", density_over_mean, z_eff, dNrec, *((float *)N_rec_unfiltered + HII_R_FFT_INDEX(x,y,z)));
        fprintf(LOG, "Delta=%g, z_eff=%g, dNrec=%g, Nrec=%g\n\n", density_over_mean, z_eff, dNrec, *((float *)N_rec_unfiltered + HII_R_FFT_INDEX(x,y,z)));
        fflush(NULL);
      }
      */

    }
  }
      }
    }

    /*********************************************************************************************/
    /************  END OF CALCULATION FOR THIS REDSHIFT.  NOW WE WILL PRINT OUT THE BOXES ********/
    /*********************************************************************************************/


    fprintf(stderr, "ave R used to compute Gamma in ionized regions = %g Mpc\n", aveR / (double)Rct);
    fprintf(LOG, "ave R used to compute Gamma in ionized regions = %g Mpc\n", aveR / (double)Rct);


    if (INHOMO_RECO){
      // N_rec box
      sprintf(filename, "../Boxes/Nrec_z%06.2f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT,HII_FILTER, MFP, HII_DIM, BOX_LEN);
      if (!(F = fopen(filename, "wb"))){
      sprintf(error_message, "find_HII_bubbles: ERROR: unable to open file for writting Nrec box!\n");
      goto CLEANUP;
      }
      else {
        for (i=0; i<HII_DIM; i++){
           for (j=0; j<HII_DIM; j++){
              for (k=0; k<HII_DIM; k++){
                 if(fwrite((float *)N_rec_unfiltered + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, F)!=1){
                     sprintf(error_message, "find_HII_bubbles.c: Write error occured while writting N_rec box.\n");
                     goto CLEANUP;
                 }
              }
           }
        }
        fclose(F);
        F = NULL;
     }

     // Fcoll_uni box -MG
     sprintf(filename, "../Boxes/Fcoll_z%06.2f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, HII_FILTER, MFP, HII_DIM, BOX_LEN);
     if (!(F = fopen(filename, "wb"))) {
       sprintf(error_message, "find_HII_bubbles.c: ERROR: unable to open file for writing Fcoll box.\n");
       goto CLEANUP;
     }
     else {
       for (x = 0; x < HII_DIM; x++) {
         for (y = 0; y < HII_DIM; y++) {
           for (z = 0; z < HII_DIM; z++) {
             if (fwrite(Fcoll_uni + HII_R_FFT_INDEX(x,y,z), sizeof(float), 1, F) != 1) {
                sprintf(error_message, "find_HII_bubbles.c: Write error occured while writing Fcoll box.\n");
                goto CLEANUP;
             }
           }
         }
       }
       fclose(F);
       F = NULL;
     }

     // Fcoll_fs box -MG
     sprintf(filename, "../Boxes/Fcollfs_z%06.2f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, HII_FILTER, MFP, HII_DIM, BOX_LEN);
     if (!(F = fopen(filename, "wb"))) {
       sprintf(error_message, "find_HII_bubbles.c: ERROR: unable to open file for writing Fcollfs box.\n");
       goto CLEANUP;
     }
     else {
       for (x = 0; x < HII_DIM; x++) {
         for (y = 0; y < HII_DIM; y++) {
           for (z = 0; z < HII_DIM; z++) {
             if (fwrite(Fcoll_fs + HII_R_FFT_INDEX(x,y,z), sizeof(float), 1, F) != 1) {
                sprintf(error_message, "find_HII_bubbles.c: Write error occured while writing Fcollfs box.\n");
                goto CLEANUP;
             }
           }
         }
       }
       fclose(F);
       F = NULL;
     }

     // mass box -MG
     sprintf(filename, "../Boxes/mass_z%06.2f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, HII_FILTER, MFP, HII_DIM, BOX_LEN);
     if (!(F = fopen(filename, "wb"))) {
       sprintf(error_message, "find_HII_bubbles.c: ERROR: unable to open file for writing mass box.\n");
       goto CLEANUP;
     }
     else {
       for (x = 0; x < HII_DIM; x++) {
         for (y = 0; y < HII_DIM; y++) {
           for (z = 0; z < HII_DIM; z++) {
             if (fwrite(mass_box + HII_R_FFT_INDEX(x,y,z), sizeof(float), 1, F) != 1) {
                sprintf(error_message, "find_HII_bubbles.c: Write error occured while writing mass box.\n");
                goto CLEANUP;
             }
           }
         }
       }
       fclose(F);
       F = NULL;
     }

     // collapsed mass box -MG
     sprintf(filename, "../Boxes/collmass_z%06.2f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, HII_FILTER, MFP, HII_DIM, BOX_LEN);
     if (!(F = fopen(filename, "wb"))) {
       sprintf(error_message, "find_HII_bubbles.c: ERROR: unable to open file for writing collapsed mass box.\n");
       goto CLEANUP;
     }
     else {
       for (x = 0; x < HII_DIM; x++) {
         for (y = 0; y < HII_DIM; y++) {
           for (z = 0; z < HII_DIM; z++) {
             if (fwrite(collapsed_mass + HII_R_FFT_INDEX(x,y,z), sizeof(float), 1, F) != 1) {
                sprintf(error_message, "find_HII_bubbles.c: Write error occured while writing collapsed mass box.\n");
                goto CLEANUP;
             }
           }
         }
       }
       fclose(F);
       F = NULL;
     }

     // stellar mass box 1 -MG
     sprintf(filename, "../Boxes/starmass_z%06.2f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, HII_FILTER, MFP, HII_DIM, BOX_LEN);
     if (!(F = fopen(filename, "wb"))) {
       sprintf(error_message, "find_HII_bubbles.c: ERROR: unable to open file for writing stellar mass box 1.\n");
       goto CLEANUP;
     }
     else {
       for (x = 0; x < HII_DIM; x++) {
         for (y = 0; y < HII_DIM; y++) {
           for (z = 0; z < HII_DIM; z++) {
             if (fwrite(star_mass + HII_R_FFT_INDEX(x,y,z), sizeof(float), 1, F) != 1) {
                sprintf(error_message, "find_HII_bubbles.c: Write error occured while writing stellar mass box 1.\n");
                goto CLEANUP;
             }
           }
         }
       }
       fclose(F);
       F = NULL;
     }

     // stellar mass box 2 -MG
     /*sprintf(filename, "../Boxes/stellarmass_z%06.2f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, HII_FILTER, MFP, HII_DIM, BOX_LEN);
     if (!(F = fopen(filename, "wb"))) {
       sprintf(error_message, "find_HII_bubbles.c: ERROR: unable to open file for writing stellar mass box 2.\n");
       goto CLEANUP;
     }
     else {
       for (x = 0; x < HII_DIM; x++) {
         for (y = 0; y < HII_DIM; y++) {
           for (z = 0; z < HII_DIM; z++) {
             if (fwrite(stellar_mass + HII_R_FFT_INDEX(x,y,z), sizeof(float), 1, F) != 1) {
                sprintf(error_message, "find_HII_bubbles.c: Write error occured while writing stellar mass box 2.\n");
                goto CLEANUP;
             }
           }
         }
       }
       fclose(F);
       F = NULL;
     }*/

     // SFR box -MG
     sprintf(filename, "../Boxes/SFR_z%06.2f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, HII_FILTER, MFP, HII_DIM, BOX_LEN);
     if (!(F = fopen(filename, "wb"))) {
       sprintf(error_message, "find_HII_bubbles.c: ERROR: unable to open file for writing SFR box.\n");
       goto CLEANUP;
     }
     else {
       for (x = 0; x < HII_DIM; x++) {
         for (y = 0; y < HII_DIM; y++) {
           for (z = 0; z < HII_DIM; z++) {
             if (fwrite(SFR + HII_R_FFT_INDEX(x,y,z), sizeof(float), 1, F) != 1) {
                sprintf(error_message, "find_HII_bubbles.c: Write error occured while writing SFR box.\n");
                goto CLEANUP;
             }
           }
         }
       }
       fclose(F);
       F = NULL;
     }

     // metallicity box -MG
     sprintf(filename, "../Boxes/metallicity_z%06.2f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, HII_FILTER, MFP, HII_DIM, BOX_LEN);
     if (!(F = fopen(filename, "wb"))) {
       sprintf(error_message, "find_HII_bubbles.c: ERROR: unable to open file for writing metallicity box.\n");
       goto CLEANUP;
     }
     else {
       for (x = 0; x < HII_DIM; x++) {
         for (y = 0; y < HII_DIM; y++) {
           for (z = 0; z < HII_DIM; z++) {
             if (fwrite(metallicity + HII_R_FFT_INDEX(x,y,z), sizeof(float), 1, F) != 1) {
                sprintf(error_message, "find_HII_bubbles.c: Write error occured while writing metallicity box.\n");
                goto CLEANUP;
             }
           }
         }
       }
       fclose(F);
       F = NULL;
     }

     // Lya_igm box -MG
     sprintf(filename, "../Boxes/Lyaigm_z%06.2f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, HII_FILTER, MFP, HII_DIM, BOX_LEN);
     if (!(F = fopen(filename, "wb"))) {
       sprintf(error_message, "find_HII_bubbles.c: ERROR: unable to open file for writing Lya_igm box.\n");
       goto CLEANUP;
     }
     else {
       for (x = 0; x < HII_DIM; x++) {
         for (y = 0; y < HII_DIM; y++) {
           for (z = 0; z < HII_DIM; z++) {
             if (fwrite(Lya_igm + HII_R_FFT_INDEX(x,y,z), sizeof(float), 1, F) != 1) {
                sprintf(error_message, "find_HII_bubbles.c: Write error occured while writing Lya_igm box.\n");
                goto CLEANUP;
             }
           }
         }
       }
       fclose(F);
       F = NULL;
     }

     if (USE_CUSTOM_SEDS && compute_seds) {
       // Lya_sf intensity box -MG
       sprintf(filename, "../Boxes/Lyasf_z%06.2f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, HII_FILTER, MFP, HII_DIM, BOX_LEN);
       if (!(F = fopen(filename, "wb"))) {
         sprintf(error_message, "find_HII_bubbles.c: ERROR: unable to open file for writing Lya_sf box.\n");
         goto CLEANUP;
       }
       else {
         for (x = 0; x < HII_DIM; x++) {
           for (y = 0; y < HII_DIM; y++) {
             for (z = 0; z < HII_DIM; z++) {
               if (fwrite(Lya_sf + HII_R_FFT_INDEX(x,y,z), sizeof(float), 1, F) != 1) {
                  sprintf(error_message, "find_HII_bubbles.c: Write error occured while writing Lya_sf box.\n");
                  goto CLEANUP;
               }
             }
           }
         }
         fclose(F);
         F = NULL;
       }

       // Halpha_sf box -MG
       sprintf(filename, "../Boxes/Halphasf_z%06.2f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, HII_FILTER, MFP, HII_DIM, BOX_LEN);
       if (!(F = fopen(filename, "wb"))) {
         sprintf(error_message, "find_HII_bubbles.c: ERROR: unable to open file for writing Halpha_sf box.\n");
         goto CLEANUP;
       }
       else {
         for (x = 0; x < HII_DIM; x++) {
           for (y = 0; y < HII_DIM; y++) {
             for (z = 0; z < HII_DIM; z++) {
               if (fwrite(Halpha_sf + HII_R_FFT_INDEX(x,y,z), sizeof(float), 1, F) != 1) {
                  sprintf(error_message, "find_HII_bubbles.c: Write error occured while writing Halpha_sf box.\n");
                  goto CLEANUP;
               }
             }
           }
         }
         fclose(F);
         F = NULL;
       }
     }

      // Write z_re in the box
      sprintf(filename, "../Boxes/z_first_ionization_z%06.2f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, HII_FILTER, MFP, HII_DIM, BOX_LEN);
      if (!(F = fopen(filename, "wb"))){
  sprintf(error_message, "find_HII_bubbles: ERROR: unable to open file for writting z_re box!\n");
  goto CLEANUP;
      }
      else{
  if (mod_fwrite(z_re, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
    sprintf(error_message, "find_HII_bubbles: ERROR: unable to open file for writting z_re box!\n");
    goto CLEANUP;
  }
  fclose(F);
  F = NULL;
}

      // Gamma12 box
      sprintf(filename, "../Boxes/Gamma12aveHII_z%06.2f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, HII_FILTER, MFP, HII_DIM, BOX_LEN);
      if (!(F = fopen(filename, "wb"))){
  sprintf(error_message, "find_HII_bubbles: ERROR: unable to open file for writting gamma box!\n");
  goto CLEANUP;
      }
      else{
  if (mod_fwrite(Gamma12, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
    sprintf(error_message, "find_HII_bubbles.c: Write error occured while writting gamma box.\n");
    goto CLEANUP;
  }
  fclose(F);
  F = NULL;
      }
    }


    // print out the xH box
    switch(FIND_BUBBLE_ALGORITHM){
    case 2:
    if (HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY != 0) {
        if (USE_HALO_FIELD)
      sprintf(filename, "../Boxes/xH_z%06.2f_nf%f_Fstar%.4%.4f_Fesc%.4f_escPL%.4f_Mturn%.2e_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, HII_FILTER, MFP, HII_DIM, BOX_LEN);
        else
      sprintf(filename, "../Boxes/xH_nohalos_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_Mturn%.2e_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, HII_FILTER, MFP, HII_DIM, BOX_LEN);
    }
    else {
        if (USE_HALO_FIELD)
        sprintf(filename, "../Boxes/xH_z%06.2f_nf%f_eff%.1f_effPLindex0_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, ION_EFF_FACTOR, HII_FILTER, M_MIN, MFP, HII_DIM, BOX_LEN);
        else
        sprintf(filename, "../Boxes/xH_nohalos_z%06.2f_nf%f_eff%.1f_effPLindex0_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, ION_EFF_FACTOR, HII_FILTER, M_MIN, MFP, HII_DIM, BOX_LEN);
      }
      break;
    default:
    if (HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY != 0) {
        if (USE_HALO_FIELD)
      sprintf(filename, "../Boxes/sphere_xH_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_Mturn%.2e_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, HII_FILTER, MFP, HII_DIM, BOX_LEN);
        else
      sprintf(filename, "../Boxes/sphere_xH_nohalos_z%06.2f_nf%f_Fstar%.4f_starPL%.4f_Fesc%.4f_escPL%.4f_Mturn%.2e_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, F_STAR10, ALPHA_STAR, F_ESC10, ALPHA_ESC, M_TURN, HII_FILTER, MFP, HII_DIM, BOX_LEN);
    }
    else {
        if (USE_HALO_FIELD)
        sprintf(filename, "../Boxes/sphere_xH_z%06.2f_nf%f_eff%.1f_effPLindex0_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, ION_EFF_FACTOR, HII_FILTER, M_MIN, MFP, HII_DIM, BOX_LEN);
        else
        sprintf(filename, "../Boxes/sphere_xH_nohalos_z%06.2f_nf%f_eff%.1f_effPLindex0_HIIfilter%i_Mmin%.1e_RHIImax%.0f_%i_%.0fMpc", REDSHIFT, global_xH, ION_EFF_FACTOR, HII_FILTER, M_MIN, MFP, HII_DIM, BOX_LEN);
    }
    }
    if (!(F = fopen(filename, "wb"))){
        fprintf(stderr, "find_HII_bubbles: ERROR: unable to open file %s for writting!\n", filename);
        fprintf(LOG, "find_HII_bubbles: ERROR: unable to open file %s for writting!\n", filename);
        global_xH = -1;
    }
    else {
        fprintf(LOG, "Neutral fraction is %f\nNow writting xH box at %s\n", global_xH, filename);
        fprintf(stderr, "Neutral fraction is %f\nNow writting xH box at %s\n", global_xH, filename);
        fflush(LOG);
        if (mod_fwrite(xH, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
            fprintf(stderr, "find_HII_bubbles.c: Write error occured while writting xH box.\n");
            fprintf(LOG, "find_HII_bubbles.c: Write error occured while writting xH box.\n");
            global_xH = -1;
        }
        fclose(F);
      F = NULL;
    }

    /* TEST computing time */
    end = clock();
    time_taken = (end - start)/(CLOCKS_PER_SEC);
    printf("\n");
    printf("time taken = %1.6e sec \n", time_taken);
    printf("\n");


    fprintf(stderr, error_message);
      fprintf(LOG, error_message);
      fftwf_free(xH);
      fftwf_free(z_re);
      fftwf_free(Gamma12);
      fclose(LOG);
      fftwf_free(deltax_unfiltered);
      fftwf_free(deltax_filtered);
      fftwf_free(M_coll_unfiltered);
      fftwf_free(M_coll_filtered);
      fftwf_cleanup_threads();
      if (F){ fclose(F);}
      free_ps();
      free_MHR();
      fftwf_free(xe_filtered);
      fftwf_free(xe_unfiltered);
      fftwf_free(N_rec_unfiltered);
      fftwf_free(N_rec_filtered);
      if (HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY) {destroy_21cmMC_arrays();}

      return 0;
}
