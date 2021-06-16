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
  Program compute_line_fields reads previously computed boxes for metallicity and stellar mass
  and uses a look up table of SEDs and emissivities to compute radiation fields of various
  emission lines for each cell in the box.
*/

//make sure these are the same high and low as the selected driver file

#define POP2_METALLICITY_SAMPLES (int) 24
#define POP2_ION_PARAM_SAMPLES (int) 4
#define POP3_ION_SAMPLES (int) 3
#define POP3_COLUMNS (int) 8

#define HALPHA_POP2_FILENAME (const char *) "../External_tables/Ha_PopII_table_bursty.dat"
#define LYA_POP2_FILENAME (const char *) "../External_tables/Lya_PopII_table_bursty.dat"
#define O2_POP2_FILENAME (const char *) "../External_tables/O2_PopII_table_bursty.dat"
#define O3_POP2_FILENAME (const char *) "../External_tables/O3_PopII_table_bursty.dat"
#define HeII_POP2_FILENAME (const char *) "../External_tables/HeII_PopII_table_bursty.dat"
#define POP3_FILENAME (const char *) "../External_tables/pop3_luminosity.dat"

//int init_pop2_tables(float Lya_table_2[POP2_METALLICITY_SAMPLES+1][POP2_ION_PARAM_SAMPLES+1], float Halpha_table_2[POP2_METALLICITY_SAMPLES+1][POP2_ION_PARAM_SAMPLES+1]);
//int init_pop3_table(double table[POP3_ION_SAMPLES][POP3_COLUMNS]);
//float lookup_pop2(float table[POP2_METALLICITY_SAMPLES+1][POP2_ION_PARAM_SAMPLES+1], float metallicity, int ion_param);
//float lookup_pop3(double table[POP3_ION_SAMPLES][POP3_COLUMNS], int ion_param, int col);
//float get_luminosity(float table_pop2[POP2_METALLICITY_SAMPLES+1][POP2_ION_PARAM_SAMPLES+1], double table_pop3[POP3_ION_SAMPLES][POP3_COLUMNS], float metallicity, int ion_param, double fPopIII, int col);

int init_pop2_tables(float Lya_table_2[POP2_METALLICITY_SAMPLES+1][POP2_ION_PARAM_SAMPLES+1], float Halpha_table_2[POP2_METALLICITY_SAMPLES+1][POP2_ION_PARAM_SAMPLES+1], float O2_table_2[POP2_METALLICITY_SAMPLES+1][POP2_ION_PARAM_SAMPLES+1], float O3_table_2[POP2_METALLICITY_SAMPLES+1][POP2_ION_PARAM_SAMPLES+1], float HeII_table_2[POP2_METALLICITY_SAMPLES+1][POP2_ION_PARAM_SAMPLES+1]) {
  int i;
  FILE *F;
  //initialize Lya table
  if (!(F = fopen(LYA_POP2_FILENAME, "r"))) {
    fprintf(stderr, "Unable to open file: Lya_table_pop2.dat for reading\nAborting\n");
    return -1;
  }
  for (i = 0; i < POP2_METALLICITY_SAMPLES + 1; i++) {
    fscanf(F, "%e %e %e %e %e", &Lya_table_2[i][0], &Lya_table_2[i][1], &Lya_table_2[i][2], &Lya_table_2[i][3], &Lya_table_2[i][4]);
  }
  fclose(F);
  //initialize Halpha table
  if (!(F = fopen(HALPHA_POP2_FILENAME, "r"))) {
    fprintf(stderr, "Unable to open file: Halpha_table_pop2.dat for reading\nAborting\n");
    return -1;
  }
  for (i = 0; i < POP2_METALLICITY_SAMPLES + 1; i++) {
    fscanf(F, "%e %e %e %e %e", &Halpha_table_2[i][0], &Halpha_table_2[i][1], &Halpha_table_2[i][2], &Halpha_table_2[i][3], &Halpha_table_2[i][4]);
  }
  fclose(F);
    //initialize O2 table
  if (!(F = fopen(O2_POP2_FILENAME, "r"))) {
    fprintf(stderr, "Unable to open file: Halpha_table_pop2.dat for reading\nAborting\n");
    return -1;
  }
  for (i = 0; i < POP2_METALLICITY_SAMPLES + 1; i++) {
    fscanf(F, "%e %e %e %e %e", &O2_table_2[i][0], &O2_table_2[i][1], &O2_table_2[i][2], &O2_table_2[i][3], &O2_table_2[i][4]);
  }
  fclose(F);
      //initialize O3 table
  if (!(F = fopen(O3_POP2_FILENAME, "r"))) {
    fprintf(stderr, "Unable to open file: Halpha_table_pop2.dat for reading\nAborting\n");
    return -1;
  }
  for (i = 0; i < POP2_METALLICITY_SAMPLES + 1; i++) {
    fscanf(F, "%e %e %e %e %e", &O3_table_2[i][0], &O3_table_2[i][1], &O3_table_2[i][2], &O3_table_2[i][3], &O3_table_2[i][4]);
  }
  fclose(F);

      //initialize HeII table
  if (!(F = fopen(HeII_POP2_FILENAME, "r"))) {
    fprintf(stderr, "Unable to open file: Halpha_table_pop2.dat for reading\nAborting\n");
    return -1;
  }
  for (i = 0; i < POP2_METALLICITY_SAMPLES + 1; i++) {
    fscanf(F, "%e %e %e %e %e", &HeII_table_2[i][0], &HeII_table_2[i][1], &HeII_table_2[i][2], &HeII_table_2[i][3], &HeII_table_2[i][4]);
  }
  fclose(F);

  return 0;
}

int init_pop3_table(double table[POP3_ION_SAMPLES][POP3_COLUMNS]) {
  int i;
  FILE *F;
  //initialize pop 3 table
  if (!(F = fopen(POP3_FILENAME, "r"))) {
    fprintf(stderr, "Unable to open file: pop3_luminosity.dat for reading\nAborting\n");
    return -1;
  }
  for (i = 0; i < POP3_ION_SAMPLES; i++) {
    fscanf(F, "%le %le %le %le %le %le %le %le", &table[i][0], &table[i][1], &table[i][2], &table[i][3], &table[i][4], &table[i][5], &table[i][6], &table[i][7]);
    //printf("%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\n", table[i][0], table[i][1], table[i][2], table[i][3], table[i][4], table[i][5], table[i][6], table[i][7]);
  }
  fclose(F);

  return 0;
}
/*
//read the contents of the files at each z to respective boxes
void read_boxes(float z) {
  int i, j, k;
  FILE *test;
  //metallicity box
  //sprintf(filename, "../Boxes/metallicity_z%06.2f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", z, HII_FILTER, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
  sprintf(filename, "../Boxes/metallicity_z016.83_HIIfilter1_RHIImax50_200_300Mpc");
  if (!(test = fopen(filename, "rb"))) {
    fprintf(stderr, "Error opening file %s for reading.\n", filename);
    error = 1;
    return;
  }
  else {
    for (i = 0; i < HII_DIM; i++){
      for (j = 0; j < HII_DIM; j++){
        for (k = 0; k < HII_DIM; k++){
          float val;
          int output = fread(&val, sizeof(float), 1, test);
          fprintf(stderr, "i,j,k: %d %d %d   val: %f   output: %d\n", i, j, k, val, output);
          metallicity[HII_R_FFT_INDEX(i,j,k)] = val;
          /*if (fread((float *)metallicity + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, test) != 1){
            fprintf(stderr, "Error reading file %s.\n", filename);
            error = 1;
            return;
          }
        }
      }
    }
    fclose(test);
  }
  file = NULL;
  //starmass box
  sprintf(filename, "../Boxes/starmass_z%06.2f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", z, HII_FILTER, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
  if (!(file = fopen(filename, "rb"))) {
    fprintf(stderr, "Error opening file %s for reading.\n", filename);
    error = 1;
    return;
  }
  else {
    for (i = 0; i < HII_DIM; i++){
      for (j = 0; j < HII_DIM; j++){
        for (k = 0; k < HII_DIM; k++){
          if (fread((float *)starmass + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, file) != 1){
            fprintf(stderr, "Error reading file %s.\n", filename);
            error = 1;
            return;
          }
        }
      }
    }
    fclose(file);
  }
  file = NULL;
}

//write the current contents of each emline box to a file
void write_boxes(float z) {
  int i, j, k;
  //metallicity box
  sprintf(filename, "../Boxes/Halpha_z%06.2f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", z, HII_FILTER, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
  if (!(file = fopen(filename, "wb"))) {
    fprintf(stderr, "Error opening file %s for writing.\n", filename);
    error = 1;
    return;
  }
  else {
    for (i = 0; i < HII_DIM; i++){
      for (j = 0; j < HII_DIM; j++){
        for (k = 0; k < HII_DIM; k++){
          if (fwrite((float *)Halpha_box + HII_R_FFT_INDEX(i,j,k), sizeof(float), 1, file) != 1){
            fprintf(stderr, "Error writing file %s.\n", filename);
            error = 1;
            return;
          }
        }
      }
    }
    fclose(file);
  }
  file = NULL;
}*/
float lookup_pop2(float table[POP2_METALLICITY_SAMPLES+1][POP2_ION_PARAM_SAMPLES+1], float metallicity, int ion_param) {
  float val1, val2;
  float metal1, metal2;
  if(metallicity <= table[1][0]) { //edge case where metallicity is lower than lowest value
    val1 = table[1][ion_param];
    val2 = table[2][ion_param];
    metal1 = table[1][0];
    metal2 = table[2][0];
  }
  else if(metallicity >= table[POP2_METALLICITY_SAMPLES][0]) { //edge case where metallicity is higher than highest value
    val1 = table[POP2_METALLICITY_SAMPLES-1][ion_param];
    val2 = table[POP2_METALLICITY_SAMPLES][ion_param];
    metal1 = table[POP2_METALLICITY_SAMPLES-1][0];
    metal2 = table[POP2_METALLICITY_SAMPLES][0];
  }
  else {//normal case, find correct rows to linearly interpolate between
    int row;
    for (row = 2; row <= POP2_METALLICITY_SAMPLES; row++) {
      if (table[row][0] >= metallicity) {//if metallicity is within range
        val1 = table[row-1][ion_param];
        val2 = table[row][ion_param];
        metal1 = table[row-1][0];
        metal2 = table[row][0];
        break;
      }
    }
  }
  //y = b + mx
  return val1 + (val2-val1) * ((metallicity-metal1) / (metal2-metal1));
}

float lookup_pop3(double table[POP3_ION_SAMPLES][POP3_COLUMNS], int ion_param, int col) {
  int row, ion;
  for (row = 0; row < POP3_ION_SAMPLES; row++) {
    ion = (int) table[row][1];
    if (ion_param == ion) { //matches exactly
      return table[row][col] * table[row][POP3_COLUMNS-1] / table[row][2];
    }
    else if (ion_param < ion) { //interpolate
      double x = ((double) ion_param-table[row-1][1]) / (table[row][1]-table[row-1][1]);
      double lum = table[row-1][col] + x * (table[row][col]-table[row-1][col]);
      double q_h_sed = table[row-1][POP3_COLUMNS-1] + x * (table[row][POP3_COLUMNS-1]-table[row-1][POP3_COLUMNS-1]);
      double q_h = table[row-1][2] + x * (table[row][2]-table[row-1][2]);
      return (float) (lum * q_h_sed / q_h);
    }
  }
}

float get_luminosity(float table_pop2[POP2_METALLICITY_SAMPLES+1][POP2_ION_PARAM_SAMPLES+1], double table_pop3[POP3_ION_SAMPLES][POP3_COLUMNS], float metallicity, int ion_param, double fPopIII, int col) {
  return lookup_pop3(table_pop3, -2, col) * (float) fPopIII + lookup_pop2(table_pop2, metallicity, ion_param) * (float) (1.0 - fPopIII);
}

/*void update_cell(float z, float z_prev, unsigned long long cell) {
  float cell_metallicity = log10(metallicity[cell] / 0.012); //switch to units of log(Z_Z_sun)
  //do lookup table stuff to determine J_tot
  float Halpha_lum = lookup(Halpha_table, cell_metallicity, -2) * starmass[cell]; //normalize for total stellar mass
  Halpha_raw[cell] += Halpha_lum * (double)(3e5) * (z-z_prev) / (nu0 * hubble(z) * CELL_VOLUME);
  Halpha_box[cell] = Halpha_raw[cell] * pow(1+z, 2) / (4 * PI * pow(CMperMPC, 2));
  //emline_box2[cell] += J_alpha2;
  //emline_box3[cell] += J_alpha3;
}

//use metallicity and stellar mass with lookup table to determine radiation at current cell
void update_boxes(float z, float z_prev) {
  int i, j, k;
  //use lookup table to find radiation (and emissivity?) at each cell
  for (i = 0; i < HII_DIM; i++){
    for (j = 0; j < HII_DIM; j++){
      for (k = 0; k < HII_DIM; k++){
        update_cell(z, z_prev, HII_R_FFT_INDEX(i, j, k));
      }
    }
  }
}

int main(int argc, char **argv) {
  float redshift, redshift_prev;
  //boxes are initialized before
  Halpha_raw = (float *) malloc(sizeof(float) * HII_TOT_FFT_NUM_PIXELS);
  Halpha_box = (float *) malloc(sizeof(float) * HII_TOT_FFT_NUM_PIXELS);
  //emline_box2 = malloc(sizeof(float) * HII_TOT_FFT_NUM_PIXELS);
  //emline_box3 = malloc(sizeof(float) * HII_TOT_FFT_NUM_PIXELS);
  metallicity = (float *) malloc(sizeof(float) * HII_TOT_FFT_NUM_PIXELS);
  starmass = (float *) malloc(sizeof(float) * HII_TOT_FFT_NUM_PIXELS);
 //initialize lookup tables for luminosity values
  if (init_tables() != 0)
    return 0;

  redshift = Z_LOW * 1.0001; //follow same z step convention as in driver files
  while (redshift < Z_HIGH)
    redshift = ((1 + redshift) * ZPRIME_STEP_FACTOR - 1);
  redshift_prev = redshift; //assign this as previous redshift
  redshift = ((1 + redshift) / ZPRIME_STEP_FACTOR - 1); //go forward one step
  //for each redshift step, compute radiation and add to previous amount per cell
  error = 0;
  while (redshift >= Z_LOW && error == 0) {
    read_boxes(redshift); //read metallicity and stellar mass box at current redshift
    update_boxes(redshift, redshift_prev); //use lookup table to update cells for each emission line box
    write_boxes(redshift); //write out contents of each box
    redshift_prev = redshift;
    redshift = ((1 + redshift) / ZPRIME_STEP_FACTOR - 1); //update redshift according to trend in driver files
  }
  //de-allocate boxes here?
  free(Halpha_raw);
  free(Halpha_box);
  //free(emline_box2);
  //free(emline_box3);
  free(metallicity);
  free(starmass);
  return 0;
}*/
