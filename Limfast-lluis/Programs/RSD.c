#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"
#include "../Parameter_files/HEAT_PARAMS.H"
#include "../Parameter_files/SOURCES.H"


/*
  USAGE: RSD [-p <NUM THREADS>] <redshift> <LIM filename> <output key>

  generates RSD-included redshift space boxes for LIM signals

  NOTE: the optional argument of thread number including the -p flag, MUST
  be the first two arguments.  If these are omitted, num_threads defaults
  to NUMCORES in INIT_PARAMS.H
*/


int main(int argc, char ** argv){
  fftwf_complex *deldel_T;
  fftwf_plan plan;
  char filename[1000], psoutputdir[1000], *token;
  float *deltax, REDSHIFT, growth_factor, dDdt, pixel_x_HI, pixel_deltax, *delta_T, *v, H, dummy;
  FILE *F, *LOG;
  int i,j,k, n_x, n_y, n_z, NUM_BINS, curr_Pop, arg_offset,num_th;
  double dvdx, ave, *p_box, *k_ave, max_v_deriv;
  unsigned long long ct, *in_bin_ct, nonlin_ct, temp_ct;
  float nf, max, maxi, maxj, maxk, maxdvdx, min, mini, minj, mink, mindvdx;
  float k_x, k_y, k_z, k_mag, k_sq, k_floor, k_ceil, k_max, k_first_bin_ceil, k_factor;
  float *xH, const_factor, *Ts, T_rad, pixel_Ts_factor, curr_alphaX, curr_MminX;
  float *I_lim, *I_lim_rsd;
  double ave_Ts, min_Ts, max_Ts, temp, curr_zetaX;
  double ave_I_lim, ave_I_lim_rsd;
  int ii;
  double checkpt;

  float d1_low, d1_high, d2_low, d2_high, gradient_component, min_gradient_component, subcell_width, x_val1, x_val2, subcell_displacement;
  float RSD_pos_new, RSD_pos_new_boundary_low,RSD_pos_new_boundary_high, fraction_within, fraction_outside, cell_distance;

  float *x_pos_offset, *x_pos, *I_lim_RSD_LOS;

  x_pos = calloc(N_RSD_STEPS,sizeof(float));
  x_pos_offset = calloc(N_RSD_STEPS,sizeof(float));
  I_lim_RSD_LOS = calloc(HII_DIM,sizeof(float));

  int HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY = 0;

  /************  BEGIN INITIALIZATION ****************************/
  if (SHARP_CUTOFF) HALO_MASS_DEPENDENT_IONIZING_EFFICIENCY = 1;
  max = -1e3;
  min = 1e3;
  ave = 0;
  nonlin_ct=0.0;

  // check arguments
  if (argc < 3){
    fprintf(stderr, "USAGE: RSD [-p <NUM THREADS>] <redshift> <LIM filename> <output key> \nAborting!\n");
    return -1;
  }
  if ( (argv[1][0]=='-') && ((argv[1][1]=='p') || (argv[1][1]=='P')) ){
    // user specified num proc
    num_th = atoi(argv[2]);
    fprintf(stderr, "RSD: threading with user-specified %i threads\n", num_th);
    arg_offset = 2;
  }
  else{
    num_th = NUMCORES;
    fprintf(stderr, "RSD: threading with default %i threads\n", num_th);
    arg_offset = 0;
  }
  if (fftwf_init_threads()==0){
    fprintf(stderr, "RSD: ERROR: problem initializing fftwf threads\nAborting\n.");
    return -1;
  }
  fftwf_plan_with_nthreads(num_th);

  // open LOG file
  REDSHIFT = atof(argv[1+arg_offset]);
  system("mkdir ../Log_files");
  sprintf(filename, "../Log_files/RSD_log_file_%d", getpid());
  LOG = fopen(filename, "w");
  if (!LOG){ fprintf(stderr, "RSD.c: Error opening log file %s\n", filename);}
  T_rad = T_background(T_cmb, RADIO_EXCESS_FRAC, REDSHIFT);
  H = hubble(REDSHIFT);
  //const_factor = 27 * (OMb*hlittle*hlittle/0.023) *
  //  sqrt( (0.15/OMm/hlittle/hlittle) * (1+REDSHIFT)/10.0 );
  //system("mkdir ../Output_files/");
  //system("mkdir ../Output_files/Deldel_T_power_spec");
  system("mkdir ../Log_files");

  // get the neutral fraction and HII filter from the filename
  //strcpy(filename, argv[2+arg_offset]);
  //strtok(filename, "f");
  //token = strtok(filename, "f");
  //nf = atof(strtok(NULL, "_"));


  // initialize power spectrum
  init_ps();
  growth_factor = dicke(REDSHIFT); // normalized to 1 at z=0

  // allocate memory for line intensity box and read it in
  I_lim = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
  if (!I_lim){
    fprintf(stderr, "RSD: Error allocating memory for I_lim box\nAborting...\n");
    fprintf(LOG, "RSD: Error allocating memory for I_lim box\nAborting...\n");
    fclose(LOG); fftwf_cleanup_threads(); return -1;
  }
  if (!(F = fopen(argv[2+arg_offset], "rb"))){
    fprintf(stderr, "RSD: unable to open I_lim box at %s\nAborting...\n", argv[2+arg_offset]);
    fprintf(LOG, "RSD: unable to open I_lim box at %s\nAborting...\n", argv[2+arg_offset]);
    free(I_lim);
    fclose(LOG); fftwf_cleanup_threads(); return -1;
  }
  fprintf(stderr, "Reading in I_lim box\n");
  fprintf(LOG, "Reading in I_lim box\n");
  if (mod_fread(I_lim, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
    fprintf(stderr, "RSD: Read error occured while reading line intensity box.\n");
    fprintf(LOG, "RSD: Read error occured while reading line intensity box.\n");
    fclose(F); free(I_lim);
    fclose(LOG); fftwf_cleanup_threads(); return -1;
  }
  fclose(F);

  // allocate memory for our output rsd-corrected line intensity box
  I_lim_rsd = (float *) malloc(sizeof(float)*HII_TOT_NUM_PIXELS);
  if (!I_lim_rsd){
    fprintf(stderr, "RSD: Error allocating memory for I_lim_rsd box\nAborting...\n");
    fprintf(LOG, "RSD: Error allocating memory for I_lim_rsd box\nAborting...\n");
    free(I_lim);
    fclose(LOG); fftwf_cleanup_threads(); return -1;
  }

  // allocate memory for the velocity box and read it in
  v = (float *) malloc(sizeof(float)*HII_TOT_FFT_NUM_PIXELS);
  if (!v){
    fprintf(stderr, "RSD: Error allocating memory for velocity box\nAborting...\n");
    fprintf(LOG, "RSD: Error allocating memory for velocity box\nAborting...\n");
    free(I_lim); free(I_lim_rsd);
    fclose(LOG); fftwf_cleanup_threads(); return -1;
  }
  switch(VELOCITY_COMPONENT){
  case 1:  sprintf(filename, "../Boxes/updated_vx_z%06.2f_%i_%.0fMpc", REDSHIFT, HII_DIM, BOX_LEN);
    break;
  case 3:  sprintf(filename, "../Boxes/updated_vz_z%06.2f_%i_%.0fMpc", REDSHIFT, HII_DIM, BOX_LEN);
    break;
  default: sprintf(filename, "../Boxes/updated_vy_z%06.2f_%i_%.0fMpc", REDSHIFT, HII_DIM, BOX_LEN);
  }
  if (T_USE_VELOCITIES){
    if (!(F=fopen(filename, "rb"))){
      fprintf(stderr, "RSD: Error opening velocity file at %s\n", filename);
      fprintf(LOG, "RSD: Error opening velocity file at %s\n", filename);
      free(I_lim); free(I_lim_rsd); free(v);
      fclose(LOG); fftwf_cleanup_threads(); return -1;
    }
    for (i=0; i<HII_DIM; i++){
      for (j=0; j<HII_DIM; j++){
  for (k=0; k<HII_DIM; k++){
    if (fread((float *)v + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, F)!=1){
      fprintf(stderr, "RSD: Read error occured while reading velocity box.\n");
      fprintf(LOG, "RSD: Read error occured while reading velocity box.\n");
      fclose(F); free(I_lim); free(I_lim_rsd); free(v);
      fclose(LOG); fftwf_cleanup_threads(); return -1;
    }
  }
      }
    }
    fclose(F);
  }

  /************  END INITIALIZATION ****************************/

  x_val1 = 0.;
  x_val2 = 1.;

  subcell_width = (BOX_LEN/(float)HII_DIM)/(float)N_RSD_STEPS;

  float max_cell_distance;

  max_cell_distance = 0.;

  // let's take the derivative in k-space
  plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)v, (fftwf_complex *)v, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
  fftwf_cleanup();
  for (n_x=0; n_x<HII_DIM; n_x++){
    if (n_x>HII_MIDDLE)
      k_x =(n_x-HII_DIM) * DELTA_K;  // wrap around for FFT convention
    else
      k_x = n_x * DELTA_K;

    for (n_y=0; n_y<HII_DIM; n_y++){
      if (n_y>HII_MIDDLE)
        k_y =(n_y-HII_DIM) * DELTA_K;
      else
        k_y = n_y * DELTA_K;

      for (n_z=0; n_z<=HII_MIDDLE; n_z++){
         k_z = n_z * DELTA_K;

  // take partial deriavative along the line of sight
  switch(VELOCITY_COMPONENT){
  case 1:
    *((fftwf_complex *) v + HII_C_INDEX(n_x,n_y,n_z)) *= k_x*I/(float)HII_TOT_NUM_PIXELS;
    break;
  case 3:
    *((fftwf_complex *) v + HII_C_INDEX(n_x,n_y,n_z)) *= k_z*I/(float)HII_TOT_NUM_PIXELS;
    break;
  default:
    *((fftwf_complex *) v + HII_C_INDEX(n_x,n_y,n_z)) *= k_y*I/(float)HII_TOT_NUM_PIXELS;
  }
      }
    }
  }
  plan = fftwf_plan_dft_c2r_3d(HII_DIM, HII_DIM, HII_DIM, (fftwf_complex *)v, (float *)v, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
  fftwf_cleanup();

  min_gradient_component = 1.0;

  ave_I_lim = 0.0;
  ave_I_lim_rsd = 0.0;



  // now apply velocity correction to the Intensity

  max_v_deriv = fabs(MAX_DVDR*H);
  for (i=0; i<HII_DIM; i++){
      for (j=0; j<HII_DIM; j++){
          for (k=0; k<HII_DIM; k++){

              //gradient_component = fabs(v[HII_R_FFT_INDEX(i,j,k)]/H + 1.0);
              //fprintf(stderr, "Check gradient: %.3e\n", gradient_component);

              // Calculate the gradient-corrected line intensity box
              //if(gradient_component < FRACT_FLOAT_ERR) {
                // simply set the gradient component to the tolerance
              //  fprintf(stderr, "BANG: gradient_component = %.3e\n", gradient_component);
              //  gradient_component = FRACT_FLOAT_ERR;
              //}
              //I_lim_rsd[HII_R_INDEX(i,j,k)] = I_lim[HII_R_INDEX(i,j,k)]/gradient_component;

              dvdx = v[HII_R_FFT_INDEX(i,j,k)];

              // set maximum allowed gradient for this linear approximation
              if (fabs(dvdx) > max_v_deriv){
                  if (dvdx < 0) dvdx = -max_v_deriv;
                  else dvdx = max_v_deriv;
                  nonlin_ct++;
              }

              I_lim_rsd[HII_R_INDEX(i,j,k)] = I_lim[HII_R_INDEX(i,j,k)] / (dvdx/H + 1.0);
              ave_I_lim += I_lim[HII_R_INDEX(i,j,k)];
              //ave_I_lim_rsd += I_lim_rsd[HII_R_INDEX(i,j,k)];
          }
      }
  }

  ave_I_lim /= (float)HII_TOT_NUM_PIXELS;
  //ave_I_lim_rsd /= (float)HII_TOT_NUM_PIXELS;
  //fprintf(stderr, "CHECK BEFORE: ave_I_lim_rsd = %.3e, ave_I_lim_rsd = %.3e\n", ave_I_lim, ave_I_lim_rsd);

  fprintf(stderr, "%llu out of %llu voxels (fraction=%e) exceeded max allowed velocity gradient\n", nonlin_ct, HII_TOT_NUM_PIXELS, nonlin_ct/(double)HII_TOT_NUM_PIXELS);

  // normalised units of cell length. 0 equals beginning of cell, 1 equals end of cell
  // These are the sub-cell central positions (x_pos_offset), and the corresponding normalised value (x_pos) between 0 and 1
  for(ii=0;ii<N_RSD_STEPS;ii++) {
      x_pos_offset[ii] = subcell_width*(float)ii + subcell_width/2.;
      x_pos[ii] = x_pos_offset[ii]/( BOX_LEN/(float)HII_DIM );
  }
  // Note to convert the velocity v, to a displacement in redshift space, convert from s -> r + (1+z)*v/H(z)
  // To convert the velocity within the array v to km/s, it is a*dD/dt*delta. Where the scale factor a comes from the continuity equation
  // The array v as defined in 21cmFAST is (ik/k^2)*dD/dt*delta, as it is defined as a comoving quantity (scale factor is implicit).
  // However, the conversion between real and redshift space also picks up a scale factor, therefore the scale factors drop out and therefore
  // the displacement of the sub-cells is purely determined from the array, v and the Hubble factor: v/H.

  ave_I_lim_rsd = 0.0;

  for (i=0; i<HII_DIM; i++){
      for (j=0; j<HII_DIM; j++){

          // Generate the optical-depth for the specific line-of-sight with R.S.D
          for(k=0;k<HII_DIM;k++) {
                I_lim_RSD_LOS[k] = 0.0;
            }

            for (k=0; k<HII_DIM; k++){

                if(fabs(I_lim_rsd[HII_R_INDEX(i,j,k)]) > 0.0) {

                    if(k==0) {
                        d1_low = v[HII_R_FFT_INDEX(i,j,HII_DIM-1)]/H;
                        d2_low = v[HII_R_FFT_INDEX(i,j,k)]/H;
                    }
                    else {
                        d1_low = v[HII_R_FFT_INDEX(i,j,k-1)]/H;
                        d2_low = v[HII_R_FFT_INDEX(i,j,k)]/H;
                    }
                    // Displacements (converted from velocity) for the original cell centres straddling half of the sub-cells (cell after)
                    if(k==(HII_DIM-1)) {
                        d1_high = v[HII_R_FFT_INDEX(i,j,k)]/H;
                        d2_high = v[HII_R_FFT_INDEX(i,j,0)]/H;
                    }
                    else {
                        d1_high = v[HII_R_FFT_INDEX(i,j,k)]/H;
                        d2_high = v[HII_R_FFT_INDEX(i,j,k+1)]/H;
                    }

                    for(ii=0;ii<N_RSD_STEPS;ii++) {

                        // linearly interpolate the displacements to determine the corresponding displacements of the sub-cells
                        // Checking of 0.5 is for determining if we are left or right of the mid-point of the original cell (for the linear interpolation of the displacement)
                        // to use the appropriate cell

                        if(x_pos[ii] <= 0.5) {
                            subcell_displacement = d1_low + ( (x_pos[ii] + 0.5 ) - x_val1)*( d2_low - d1_low )/( x_val2 - x_val1 );
                        }
                        else {
                            subcell_displacement = d1_high + ( (x_pos[ii] - 0.5 ) - x_val1)*( d2_high - d1_high )/( x_val2 - x_val1 );
                        }

                        // The new centre of the sub-cell post R.S.D displacement. Normalised to units of cell width for determining it's displacement
                        RSD_pos_new = (x_pos_offset[ii] + subcell_displacement)/( BOX_LEN/(float)HII_DIM );
                        // The sub-cell boundaries of the sub-cell, for determining the fractional contribution of the sub-cell to neighbouring cells when
                        // the sub-cell straddles two cell positions
                        RSD_pos_new_boundary_low = RSD_pos_new - (subcell_width/2.)/( BOX_LEN/(float)HII_DIM );
                        RSD_pos_new_boundary_high = RSD_pos_new + (subcell_width/2.)/( BOX_LEN/(float)HII_DIM );

                        if(RSD_pos_new_boundary_low >= 0.0 && RSD_pos_new_boundary_high < 1.0) {
                            // sub-cell has remained in the original cell (just add it back to the original cell)
                            //checkpt = I_lim_rsd[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                            //fprintf(stderr, "Check checkpt: %.3e\n", checkpt);
                            I_lim_RSD_LOS[k] += I_lim_rsd[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                        }
                        else if(RSD_pos_new_boundary_low < 0.0 && RSD_pos_new_boundary_high < 0.0) {
                            // sub-cell has moved completely into a new cell (toward the observer)

                            // determine how far the sub-cell has moved in units of original cell boundary
                            cell_distance = ceil(fabs(RSD_pos_new_boundary_low))-1.;

                            // Determine the location of the sub-cell relative to the original cell binning
                            if(fabs(RSD_pos_new_boundary_high) > cell_distance) {
                                // sub-cell is entirely contained within the new cell (just add it to the new cell)
                               // check if the new cell position is at the edge of the box. If so, periodic boundary conditions
                                if(k<((int)cell_distance+1)) {
                                    I_lim_RSD_LOS[k-((int)cell_distance+1) + HII_DIM] += I_lim_rsd[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                                }
                                else {
                                    I_lim_RSD_LOS[k-((int)cell_distance+1)] += I_lim_rsd[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                                }
                            }
                            else {
                                // sub-cell is partially contained within the cell

                                // Determine the fraction of the sub-cell which is in either of the two original cells
                                fraction_outside = (fabs(RSD_pos_new_boundary_low) - cell_distance)/(subcell_width/( BOX_LEN/(float)HII_DIM ));
                                fraction_within = 1. - fraction_outside;

                                // Check if the first part of the sub-cell is at the box edge
                                if(k<(((int)cell_distance))) {
                                    I_lim_RSD_LOS[k-((int)cell_distance) + HII_DIM] += fraction_within*I_lim_rsd[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                                }
                                else {
                                    I_lim_RSD_LOS[k-((int)cell_distance)] += fraction_within*I_lim_rsd[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                                }
                                // Check if the second part of the sub-cell is at the box edge
                                if(k<(((int)cell_distance + 1))) {
                                    I_lim_RSD_LOS[k-((int)cell_distance+1) + HII_DIM] += fraction_outside*I_lim_rsd[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                                }
                                else {
                                    I_lim_RSD_LOS[k-((int)cell_distance+1)] += fraction_outside*I_lim_rsd[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                                }
                            }
                        }
                       else if(RSD_pos_new_boundary_low < 0.0 && (RSD_pos_new_boundary_high > 0.0 && RSD_pos_new_boundary_high < 1.0)) {
                            // sub-cell has moved partially into a new cell (toward the observer)

                            // Determine the fraction of the sub-cell which is in either of the two original cells
                            fraction_within = RSD_pos_new_boundary_high/(subcell_width/( BOX_LEN/(float)HII_DIM ));
                            fraction_outside = 1. - fraction_within;

                            // Check the periodic boundaries conditions and move the fraction of each sub-cell to the appropriate new cell
                            if(k==0) {
                                I_lim_RSD_LOS[HII_DIM-1] += fraction_outside*I_lim_rsd[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                                I_lim_RSD_LOS[k] += fraction_within*I_lim_rsd[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                            }
                            else {
                                I_lim_RSD_LOS[k-1] += fraction_outside*I_lim_rsd[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                                I_lim_RSD_LOS[k] += fraction_within*I_lim_rsd[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                            }
                        }
                        else if((RSD_pos_new_boundary_low >= 0.0 && RSD_pos_new_boundary_low < 1.0) && (RSD_pos_new_boundary_high >= 1.0)) {
                            // sub-cell has moved partially into a new cell (away from the observer)

                            // Determine the fraction of the sub-cell which is in either of the two original cells
                            fraction_outside = (RSD_pos_new_boundary_high - 1.)/(subcell_width/( BOX_LEN/(float)HII_DIM ));
                            fraction_within = 1. - fraction_outside;

                            // Check the periodic boundaries conditions and move the fraction of each sub-cell to the appropriate new cell
                            if(k==(HII_DIM-1)) {
                                I_lim_RSD_LOS[k] += fraction_within*I_lim_rsd[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                                I_lim_RSD_LOS[0] += fraction_outside*I_lim_rsd[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                            }
                            else {
                                I_lim_RSD_LOS[k] += fraction_within*I_lim_rsd[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                                I_lim_RSD_LOS[k+1] += fraction_outside*I_lim_rsd[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                            }
                        }
                       else {
                            // sub-cell has moved completely into a new cell (away from the observer)

                            // determine how far the sub-cell has moved in units of original cell boundary
                            cell_distance = floor(fabs(RSD_pos_new_boundary_high));

                            if(RSD_pos_new_boundary_low >= cell_distance) {
                                // sub-cell is entirely contained within the new cell (just add it to the new cell)

                                // check if the new cell position is at the edge of the box. If so, periodic boundary conditions
                                if(k>(HII_DIM - 1 - (int)cell_distance)) {
                                    I_lim_RSD_LOS[k+(int)cell_distance - HII_DIM] += I_lim_rsd[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                                }
                                else {
                                    I_lim_RSD_LOS[k+(int)cell_distance] += I_lim_rsd[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                                }
                            }
                            else {
                                // sub-cell is partially contained within the cell

                                // Determine the fraction of the sub-cell which is in either of the two original cells
                                fraction_outside = (RSD_pos_new_boundary_high - cell_distance)/(subcell_width/( BOX_LEN/(float)HII_DIM ));
                                fraction_within = 1. - fraction_outside;

                                // Check if the first part of the sub-cell is at the box edge
                                if(k>(HII_DIM - 1 - ((int)cell_distance-1))) {
                                    I_lim_RSD_LOS[k+(int)cell_distance-1 - HII_DIM] += fraction_within*I_lim_rsd[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                                }
                                else {
                                    I_lim_RSD_LOS[k+(int)cell_distance-1] += fraction_within*I_lim_rsd[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                                }
                                // Check if the second part of the sub-cell is at the box edge
                                if(k>(HII_DIM - 1 - ((int)cell_distance))) {
                                    I_lim_RSD_LOS[k+(int)cell_distance - HII_DIM] += fraction_outside*I_lim_rsd[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                                }
                                else {
                                    I_lim_RSD_LOS[k+(int)cell_distance] += fraction_outside*I_lim_rsd[HII_R_INDEX(i,j,k)]/(float)N_RSD_STEPS;
                                }
                           }
                        }
                    }
                }
            }

          for(k=0;k<HII_DIM;k++) {
              I_lim_rsd[HII_R_INDEX(i,j,k)] = I_lim_RSD_LOS[k];
              ave_I_lim_rsd += I_lim_rsd[HII_R_INDEX(i,j,k)];
          }
      }
  }

  ave_I_lim_rsd /= (float)HII_TOT_NUM_PIXELS;
  fprintf(stderr, "CHECK AFTER: ave_I_lim = %.3e, ave_I_lim_rsd = %.3e\n", ave_I_lim, ave_I_lim_rsd);

  // now write out the line intensity box with velocity correction
  sprintf(filename, "../Boxes/I_%s_RSD_maxdvdr%06.2f_z%06.2f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", argv[3], MAX_DVDR, REDSHIFT, HII_FILTER, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
  F = fopen(filename, "wb");
  fprintf(stderr, "Writting output RSD corrected line intensity box: %s\n", filename);
  if (mod_fwrite(I_lim_rsd, sizeof(float)*HII_TOT_NUM_PIXELS, 1, F)!=1){
    fprintf(stderr, "RSD: Write error occured while writting line intensity box.\n");
  }
  fclose(F);

  // deallocate what we aren't using anymore
 free(v);


/******  PRINT OUT THE POWERSPECTRUM  *********/

 k_factor = 1.5;
  k_first_bin_ceil = DELTA_K;
  k_max = DELTA_K*HII_DIM;
  // initialize arrays
  // ghetto counting (lookup how to do logs of arbitrary bases in c...)
  NUM_BINS = 0;
  k_floor = 0;
  k_ceil = k_first_bin_ceil;
  while (k_ceil < k_max){
    NUM_BINS++;
    k_floor=k_ceil;
    k_ceil*=k_factor;
  }

  p_box =  (double *)malloc(sizeof(double)*NUM_BINS);
  k_ave =  (double *)malloc(sizeof(double)*NUM_BINS);
  in_bin_ct = (unsigned long long *)malloc(sizeof(unsigned long long)*NUM_BINS);
  if (!p_box || !in_bin_ct || !k_ave){ // a bit sloppy, but whatever..
    fprintf(stderr, "RSD.c: Error allocating memory.\nAborting...\n");
    fprintf(LOG, "RSD.c: Error allocating memory.\nAborting...\n");
    free(I_lim_rsd); fclose(LOG);
    fftwf_cleanup_threads(); return -1;
  }
  for (ct=0; ct<NUM_BINS; ct++){
    p_box[ct] = k_ave[ct] = 0;
    in_bin_ct[ct] = 0;
  }

  deldel_T = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*HII_KSPACE_NUM_PIXELS);
  if (!deldel_T){
    fprintf(stderr, "Unable to allocate memory for the deldel_T box!\n");
    fprintf(LOG, "Unable to allocate memory for the deldel_T box!\n");
    free(I_lim_rsd); fclose(LOG); free(p_box); free(k_ave); free(in_bin_ct);
    fftwf_cleanup_threads(); return -1;
  }

  // fill-up the real-space of the deldel box
  for (i=0; i<HII_DIM; i++){
    for (j=0; j<HII_DIM; j++){
      for (k=0; k<HII_DIM; k++){
  *((float *)deldel_T + HII_R_FFT_INDEX(i,j,k)) = (I_lim_rsd[HII_R_INDEX(i,j,k)]/ave_I_lim_rsd - 1)*VOLUME/(HII_TOT_NUM_PIXELS+0.0);
  if (DIMENSIONAL_T_POWER_SPEC){
    *((float *)deldel_T + HII_R_FFT_INDEX(i,j,k)) *= ave_I_lim_rsd;
  }
  // Note: we include the V/N factor for the scaling after the fft
      }
    }
  }

  // transform to k-space
  plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)deldel_T, (fftwf_complex *)deldel_T, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
  fftwf_cleanup();

  // now construct the power spectrum file
  for (n_x=0; n_x<HII_DIM; n_x++){
    if (n_x>HII_MIDDLE)
      k_x =(n_x-HII_DIM) * DELTA_K;  // wrap around for FFT convention
    else
      k_x = n_x * DELTA_K;

    for (n_y=0; n_y<HII_DIM; n_y++){
      if (n_y>HII_MIDDLE)
  k_y =(n_y-HII_DIM) * DELTA_K;
      else
  k_y = n_y * DELTA_K;

      for (n_z=0; n_z<=HII_MIDDLE; n_z++){
  k_z = n_z * DELTA_K;

  k_mag = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);

  // now go through the k bins and update
  ct = 0;
  k_floor = 0;
  k_ceil = k_first_bin_ceil;
  while (k_ceil < k_max){
    // check if we fal in this bin
    if ((k_mag>=k_floor) && (k_mag < k_ceil)){
      in_bin_ct[ct]++;
      p_box[ct] += pow(k_mag,3)*pow(cabs(deldel_T[HII_C_INDEX(n_x, n_y, n_z)]), 2)/(2.0*PI*PI*VOLUME);
      // note the 1/VOLUME factor, which turns this into a power density in k-space

      k_ave[ct] += k_mag;
      break;
    }

    ct++;
    k_floor=k_ceil;
    k_ceil*=k_factor;
  }
      }
    }
  } // end looping through k box


  if (DIMENSIONAL_T_POWER_SPEC)
    sprintf(psoutputdir, "../Output_files/Deldel_T_power_spec");
  else
    sprintf(psoutputdir, "../Output_files/Deldel_T_power_spec/Dimensionless");
  sprintf(filename, "mkdir %s", psoutputdir);
  system(filename);

  // now lets print out the k bins
  if (T_USE_VELOCITIES){
    sprintf(filename, "%s/ps_%s_z%06.2f_%i_%.0fMpc_maxdvdr%06.2f_RSD", psoutputdir, argv[3], REDSHIFT,   HII_DIM, BOX_LEN,  MAX_DVDR);
  }
  else{
    sprintf(filename, "%s/ps_%s_z%06.2f_%i_%.0fMpc_maxdvdr%06.2f_RSD", psoutputdir, argv[3], REDSHIFT,  HII_DIM, BOX_LEN, MAX_DVDR);
  }
  F = fopen(filename, "w");
  if (!F){
    fprintf(stderr, "RSD.c: Couldn't open file %s for writting!\n", filename);
    fprintf(LOG, "RSD.c: Couldn't open file %s for writting!\n", filename);
    free(I_lim_rsd); fclose(LOG); free(p_box); free(k_ave); free(in_bin_ct); fftwf_free(deldel_T);
  }
  for (ct=1; ct<NUM_BINS; ct++){
    if (in_bin_ct[ct]>0)
      fprintf(F, "%e\t%e\t%e\n", k_ave[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0)/sqrt(in_bin_ct[ct]+0.0));
  }
  fclose(F); free(p_box); free(k_ave); free(in_bin_ct); fftwf_free(deldel_T);

  /****** END POWER SPECTRUM STUFF   ************/


  // deallocate
  free(I_lim); free(I_lim_rsd);
  fclose(LOG);
  fftwf_cleanup_threads(); return 0;
}


