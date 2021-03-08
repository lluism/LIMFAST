#ifndef _SOURCES_
#define _SOURCES_

// Function definitions for source populations
// Rick Mebane 2018

#include "INIT_PARAMS.H"
#include "ANAL_PARAMS.H"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


// more general, allows for halos with multiple populations
typedef struct{
    double (*fesc)(double, double); //escape fraction as a function of halo mass and redshift
    double (*Nion)(double, double); //Nion, generally just set to 4e4 for Pop III, 5e3 for Pop II, and weighted average if a halo has both
    double (*fstar)(double, double); //star formation efficiency, defined as fraction of a halo's gas in stars. combine both Pop II and Pop III here, differences set in Nion and fPopIII
                                     //NOTE: fstar here is CUMULATIVE, which is very important for Pop III... make sure to add up all the stellar mass that has been in the halo, not that is there currently
    double (*fx)(double, double); //boost to xray emissivity, fiducially 1 for Pop II and 100 for Pop III in Yuxiang's 21cmFAST w/ minihalos
    double (*minMass)(double); //minimum halo mass for star formation (both Pop II and Pop III)
    double (*fPopIII)(double); //mass fraction of Pop III halos at redshift z. If no mixed populations, just set to 0 or 1. Needed for spectral emissivity calculation
} sources;

// atomic cooling threshold
// used for min mass comparisons
double minMassATC(double z)
{
    double mu;
    mu = 1.22;
    return (double) TtoM(z, 1.0e4, mu);
}

// initialize source struct
sources setSources(double (*fesc)(double, double),
                   double (*Nion)(double, double),
                   double (*fstar)(double, double),
                   double (*fx)(double, double),
                   double (*minMass)(double),
                   double (*fPopIII)(double))
{
    sources s;
    s.fesc = fesc;
    s.Nion = Nion;
    s.fstar = fstar;
    s.minMass = minMass;
    s.fx = fx;
    s.fPopIII = fPopIII;
    return s;
}

double zeta_at_M(sources s, double z, double M)
{
    return 1.22 * s.fstar(z, M) * s.fesc(z, M) * s.Nion(z, M);
}

// BEGIN SPECIFIC SOURCE FUNCTIONS

double fstar_jordan(double z, double M)
{
    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read = 0;
    double table[26*61];

    const char *s = " ";
    char *token = NULL;

    int i = 0;
    int j;

    fp = fopen("../Parameter_files/smhm_uvlf.txt", "r");
    if(fp == NULL)
    {
        printf("Error opening");
        exit(EXIT_FAILURE);
    }

    while ((read=getline(&line, &len, fp)) != -1)
    {
        token = strtok(line, s);
        while(token != NULL)
        {
            table[i] = atof(token);
            token=strtok(NULL, s);
            i++;
        }
    }

    fclose(fp);

    // get fstar from table
    // Not general, but just wanted to get this working quick
    int zind, Mind;
    zind = (int) round(z-5);
    Mind = (int) round(((log10(M) - 8) / 5) * 50 + 10);
    return table[Mind * 26 + zind] * .28 / .046 * 2.0;

}

// Eq. 1 from Feng & Holder 2018
// Not really a source pop function, but putting it here for now...
// Tcmb given at z=0
double T_background(double Tcmb, double excess_frac, double z)
{
    double beta, Tr;
    beta = -2.62;
    Tr = 1.19;
    return (Tcmb + excess_frac * Tr * pow(1.428 / (1.0 + z), beta)) * (1.0 + z);
}

double Mmin_jordan(double z)
{
    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read = 0;
    double table[26*61];

    const char *s = " ";
    char *token = NULL;

    int i = 0;
    int j;

    fp = fopen("../Parameter_files/smhm_uvlf.txt", "r");
    if(fp == NULL)
    {
        printf("Error opening");
        exit(EXIT_FAILURE);
    }

    while ((read=getline(&line, &len, fp)) != -1)
    {
        token = strtok(line, s);
        while(token != NULL)
        {
            table[i] = atof(token);
            token=strtok(NULL, s);
            i++;
        }
    }

    fclose(fp);

    double logM = 7;
    int zind, Mind;
    while(logM < 10)
    {
        Mind = (int) round(((logM - 8) / 5) * 50 + 10);
        zind = (int) round(z-5);

        if(table[Mind * 26 + zind] > 0)
        {
            return pow(10.0, logM);
        }
        logM += 0.1;
    }
    return 1.0e8;

}

// just returns 1 if res is > 1
// useful for fstar and fesc functions
double maxOne(double res)
{
    if(res > 1.0) return 1.0;
    return res;
}

// Weighted by the exponential cutoff
// Check to make sure fstar isn't above one BEFORE applying cutoff
double fstar_21cmfid(double z, double M)
{
    double res;
    res = maxOne(0.05 * pow(M / 1.0e10, 0.5));
    return res * exp(-5e8 / M);
}

double fesc_21cmfid(double z, double M)
{
    double res;
    res = 0.1 * pow(M / 1.0e10, -0.5);
    return maxOne(res);
}

double fxDef(double z, double M)
{
    if(M > minMassATC(z)) return 1.0;
    return 100.0;
}

double minMassDef(double z)
{
    //return 1.0e4;
    double res;
    res = pow(10.0, -3.0/44.0 * z + 8.4); // quick linear fit to fiducial min mass in log space for testing
                                           // will do this right once we make sure everything works
    if(res > minMassATC(z)) return minMassATC(z);
    return res;
}

double minMassDefMom(double z)
{
    double res;
    res = pow(10.0, -1.0/14.0 * z + 60.0/7.0);
    if(res > minMassATC(z)) return minMassATC(z);
    return res;
}

double NionDef(double z, double M)
{
    if(M > minMassATC(z)) return 5000.0;
    return 40000.0;
}

double fstarDef(double z, double M)
{
    if(M > minMassATC(z)) return maxOne(0.03 * pow(M / 1.0e10, 0.4));
    return 0.0005;
}

double fstarTest(double z, double M)
{
    if(M > 1.0e8) return 0.05;
    return 0.005;
}

double fescDef(double z, double M)
{
    if(M > minMassATC(z)) return 0.1;
    return 1.0;
}

double etaSteveEnergy(double z, double M)
{
    return pow(pow(10, 11.5) / M, 2.0/3.0) * 9.0 / (1.0 + z);
}

double fstarSteveEnergy(double z, double M)
{
    return maxOne(1.0 / (1.0 + etaSteveEnergy(z, M)));
}

double etaSteveMom(double z, double M)
{
    return 0.2 * pow(pow(10, 11.5) / M, 1.0/3.0) * pow(9.0 / (1.0 + z), 1.0/2.0);
}

double fstarSteveMom(double z, double M)
{
    return maxOne(1.0 / (1.0 + etaSteveMom(z, M)));
}

double fstarFid(double z, double M)
{
    if (M > minMassATC(z)) return fstarSteveEnergy(z, M);
    return 0.005;
}

double fstarFidMom(double z, double M)
{
    if (M > minMassATC(z)) return fstarSteveMom(z, M);
    return 0.005;
}

double noPopIII(double z)
{
    return 0;
}

double allPopIII(double z)
{
    return 1;
}

double halfPopIII(double z) {
  return 0.5;
}

double customPopIII(double z) {
  return 0.5 * (1.0 + erf((z - 14.0)/4.0));
}


double custom_Nion(double z, double M) {
    double frac_pop3 = customPopIII(z);
    z = ((1+z)*ZPRIME_STEP_FACTOR - 1);
    float ion_pop3, ion_pop2;
    float avg = 0.0;
    char filename[1000];
    FILE *F;

    //read ion_pop3 from pop3 SED
    if (!(F = fopen(POP3_SPECTRA_FILENAME, "r"))) {
        fprintf(stderr, "custom_Nion: Unable to open file: pop3_spectra.dat for reading\n");
        ion_pop3 = 0.0;
    }
    else {
        fscanf(F, "%e", &ion_pop3);
    }
    fclose(F);

    float vals[] = {1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.0,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1.0,-1.1,
                    -1.2,-1.3,-1.4,-1.5,-1.6,-1.7,-1.8,-1.9,-2.0,-2.1,-2.2,-2.3,-2.4,-2.5,-2.6,-2.7,-2.8,-2.9,-3.0};
    float closest_val = -3.0;

    //attempt to open metallicity box at current redshift, if does not exist then avg=0
    sprintf(filename, "../Boxes/metallicity_z%06.2f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", z, HII_FILTER, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
    if (!(F = fopen(filename, "rb"))){
        fprintf(stderr, "Error opening file %s for reading.\n", filename);
    }
    else { //otherwise take average of metallicity box
        float *curr = malloc(sizeof(float));
        for (long i = 0; i < HII_DIM * HII_DIM * HII_DIM; i++) {
            if (fread(curr, sizeof(float), 1, F) != 1) {
                fprintf(stderr, "Error reading file %s.\n", filename);
                avg = 0.0;
                break;
            }
            avg += *curr;
        }
        fclose(F);
        avg /= (HII_DIM * HII_DIM * HII_DIM);
        avg = (float) log10(avg / 0.014); //in units of log/msun
        free(curr);

        //use average to find pop2 SED with closest metallicity and read value from there
        float smallest_diff = fabs(avg - vals[0]);
        for (size_t i = 0; i < sizeof(vals) / sizeof(vals[0]); i++) {
            if (fabs(avg - vals[i]) < smallest_diff) {
                smallest_diff = fabs(avg - vals[i]);
                closest_val = vals[i];
            }
        }
    }

    //read value from this file as ion_pop2
    sprintf(filename, "../External_tables/pop2_SEDs/FSPS_SED_Salpeter_Age3Myr_logZ%.1f_per_Msun.dat", closest_val);
    if (!(F = fopen(filename, "r"))) {
        fprintf(stderr, "Error opening file %s for reading.\n", filename);
        ion_pop2 = 0.0;
    }
    else {
        fscanf(F, "%e", &ion_pop2);
        fprintf(stderr, "ion_pop2: %f\n", ion_pop2);
        fclose(F);
    }
    return frac_pop3*(double)ion_pop3 + (1.0-frac_pop3)*(double)ion_pop2;
}

double histogram_Nion(double z, double M) {
    FILE *F;
    char filename[1000];
    double frac_pop3 = customPopIII(z);
    double ion_pop3;

    //multiply z one redshift step up
    double newZ = ((1+z)*ZPRIME_STEP_FACTOR - 1);

    /***** Pop III Calculations *****/

    //read ion_pop3 from pop3 SED
    if (!(F = fopen(POP3_SPECTRA_FILENAME, "r"))) {
        fprintf(stderr, "custom_Nion: Unable to open file: pop3_spectra.dat for reading\n");
        ion_pop3 = 0.0;
    }
    else {
        fscanf(F, "%e", &ion_pop3);
        fclose(F);
    }

    /***** Pop II Calculations *****/

    //range of values in SEDs
    float vals[] = {1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.0,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1.0,-1.1,
                    -1.2,-1.3,-1.4,-1.5,-1.6,-1.7,-1.8,-1.9,-2.0,-2.1,-2.2,-2.3,-2.4,-2.5,-2.6,-2.7,-2.8,-2.9,-3.0};
    //number of cells corresponding to values above
    int num_cells[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double cell_met = 0.0;
    double ion_pop2 = 0.0;

    // Attempt to open metallicity box at z
    sprintf(filename, "../Boxes/metallicity_z%06.2f_HIIfilter%i_RHIImax%.0f_%i_%.0fMpc", newZ, HII_FILTER, R_BUBBLE_MAX, HII_DIM, BOX_LEN);
    if (!(F = fopen(filename, "rb"))) { //If unsuccessful, read ion_pop2 from logz=-3.0 SED only
        //fprintf(stderr, "Error opening file %s for reading.\n", filename);
        return NionDef(z, M);
        // Read value from this file as ion_pop2
        /*sprintf(filename, "../External_tables/pop2_SEDs/FSPS_SED_Salpeter_Age3Myr_logZ%.1f_per_Msun.dat", -3.0);
        if (!(F = fopen(filename, "rb"))) { // If unsuccessful, assume 0
            fprintf(stderr, "Error opening file %s for reading.\n", filename);
            ion_pop2 = 0.0;
        }
        else { // If successfully opened, load Nion for this SED
            fscanf(F, "%e", &ion_pop2);
            fclose(F);
        }*/
    }
    else { // Successfully opened box, loop through each cell
        float *cell = malloc(sizeof(float));
        for (long i = 0; i < HII_DIM*HII_DIM*HII_DIM; i++) {
            fread(&cell, sizeof(float), 1, F); //load current value into cell
            cell_met = log10((*cell) / 0.014); //relate to solar metallicity
            // Find closest metallicity from list
            for (int i = 0; i < 41; i++) {
                if (vals[i] - 0.05 < cell_met && cell_met <= vals[i] + 0.05) {
                    num_cells[i]++; //increment corresponding index
                    break;
                }
            }
        }
        fclose(F);
        free(cell);

        // Open each SED and calc ionizing photon contribution
        float *SED_nion = malloc(sizeof(float));

        for (int i = 0; i < 41; i++) {
            if (num_cells[i] > 0) {
                // Attempt to open SED
                sprintf(filename, "../External_tables/pop2_SEDs/FSPS_SED_Salpeter_Age3Myr_logZ%.1f_per_Msun.dat", vals[i]);
                if (!(F = fopen(filename, "r"))) { // Unsuccessful, assume Nion for this SED is 0
                    fprintf(stderr, "Error opening file %s for reading.\n", filename);
                    *SED_nion = 0.0;
                }
                else { // Successfully opened, load Nion for SED
                    fscanf(F, "%e", &SED_nion);
                    fclose(F);
                }
                // apply Nion = (Nion of SED) * (num cells in range) / (total cells)
                ion_pop2 += (double) (*SED_nion) * num_cells[i] / (HII_DIM*HII_DIM*HII_DIM);
            }
        }

        free(SED_nion);
    }

    return frac_pop3*ion_pop3 + (1.0-frac_pop3)*ion_pop2;
}

double PopIII_frac_test(double z)
{
    if(z < 20) return 0;
    if(z > 25) return 1.0;
    return 0.5;
}

double etaSteveMom_fig2(double z, double M)
{
    return 5.0 * pow(pow(10, 11.5) / M, 1.0/3.0) * pow(9.0 / (1.0 + z), 0.5);
}
double fshockSteveMom_fig2(double z, double M)
{
    return maxOne(0.47 * pow(pow(10, 12.0) / M, 0.25) * pow((1.0 + z) / 4.0, 0.38));
}
double fstarSteveMom_fig2(double z, double M)
{
    return maxOne(fshockSteveMom_fig2(z, M) / (5.0 + etaSteveMom_fig2(z, M)));
}
double fstarFidMom_fig2(double z, double M)
{
    if (M > minMassATC(z)) return fstarSteveMom_fig2(z, M);
    return 0.005;
}

double minMass_21cmfid(double z) {return 5e8 / 50.0;}

// source functions for comparing with sharp cutoff model
double minMassCutoff(double z) {return 5e8;}
double NionCutoff(double z, double M) {return 5000;}
double fescCutoff(double z, double M) {return 0.1;}
double fstarCutoff(double z, double M) {return 0.05;}
double fxCutoff(double z, double M) {return 1.0;}

sources MomPopIII()
{
    sources s;
    s = setSources(fescDef, NionDef, fstarFidMom, fxDef, minMassDefMom, noPopIII);
    return s;
}

sources EnergyPopIII()
{
    sources s;
    s = setSources(fescDef, NionDef, fstarFid, fxDef, minMassDef, noPopIII);
    return s;
}

sources MomPopII()
{
    sources s;
    s = setSources(fescDef, NionDef, fstarSteveMom, fxDef, minMassATC, noPopIII);
    return s;
}

sources EnergyPopII()
{
    sources s;
    s = setSources(fescDef, NionDef, fstarSteveEnergy, fxDef, minMassATC, noPopIII);
    return s;
}

sources fid_21cmFAST()
{
    sources s;
    s = setSources(fesc_21cmfid, NionCutoff, fstar_21cmfid, fxCutoff, minMass_21cmfid, noPopIII);
    return s;
}

sources fig2() {
  return setSources(fesc_21cmfid, custom_Nion, fstarFidMom_fig2, fxDef, minMassDefMom, customPopIII);
}

double fstar_test(double z, double M) {
    return 0.05;
}

double minmass_test(double z) {
    return 1e10;
}

sources test() {
    return setSources(fesc_21cmfid, histogram_Nion, fstarFidMom_fig2, fxDef, minMassDefMom, noPopIII);
}

// set the sources here
sources defaultSources()
{
    sources s;
    //change according to desired model
    //s = MomPopIII();
    s = fig2();
    //s = test();
    return s;
}

#endif
