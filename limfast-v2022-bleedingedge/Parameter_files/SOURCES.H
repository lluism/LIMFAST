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

// Include the bathtub model tables
#include "../Cosmo_c_files/compute_bathtub.c"
#include "../Cosmo_c_files/compute_nion.c"
#include "../Cosmo_c_files/compute_all_lines.c"


// more general, allows for halos with multiple populations
typedef struct{
    double (*fesc)(double, double); //escape fraction as a function of halo mass and redshift
    double (*Nion)(double, double); //Nion, generally just set to 4e4 for Pop III, 5e3 for Pop II, and weighted average if a halo has both
    double (*fstar)(double, double); //star formation efficiency, defined as fraction of a halo's gas in stars. combine both Pop II and Pop III here, differences set in Nion and fPopIII
                                     //NOTE: fstar here is CUMULATIVE, which is very important for Pop III... make sure to add up all the stellar mass that has been in the halo, not that is there currently
    double (*sfr)(double, double);
    double (*Mg)(double, double);
    double (*Ms)(double, double);
    double (*metal)(double, double);
    double (*linelum)(double, double, double, int);
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
                   double (*sfr)(double, double),
                   double (*Mg)(double, double),
                   double (*Ms)(double, double),
                   double (*metal)(double, double),
                   double (*linelum)(double, double, double, int),
                   double (*fx)(double, double),
                   double (*minMass)(double),
                   double (*fPopIII)(double))
{
    sources s;
    s.fesc = fesc;
    s.Nion = Nion;
    s.fstar = fstar;
    s.sfr = sfr;
    s.Mg = Mg;
    s.Ms = Ms;
    s.metal = metal;
    s.linelum = linelum;
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
    if(M > minMassATC(z)) return 0.100;
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

double PopIII_frac_test(double z)
{
    if(z < 20) return 0;
    if(z > 25) return 1.0;
    return 0.5;
}

double minMass_21cmfid(double z) {return 5e8 / 50.0;}

// source functions for comparing with sharp cutoff model
double minMassCutoff(double z) {return 5e8;}
double NionCutoff(double z, double M) {return 5000;}
double fescCutoff(double z, double M) {return 0.1;}
double fstarCutoff(double z, double M) {return 0.05;}
double fxCutoff(double z, double M) {return 1.0;}




/*
Halo properties as a function of redshift and halo mass from bathtub model
*/
// Star formation rate
double sfr_from_bathtub(double z, double M) {
    if (z >= 30.0) {return 1.0e-9;}
    else{
      if (M >= 1e16) {return get_sfr(fmax(z, 5.001), 9.999e15);}
      if (M <= 1e7) {return get_sfr(fmax(z, 5.001), 1.001e7);}
      else {return get_sfr(fmax(z, 5.001), M);}
    }
}

// Gas mass
double Mg_from_bathtub(double z, double M) {
    if (z >= 30.0) {return 1.0e-9;}
    else{
      if (M >= 1e16) {return get_Mg(fmax(z, 5.001), 9.999e15);}
      if (M <= 1e7) {return get_Mg(fmax(z, 5.001), 1.001e7);}
      else {return get_Mg(fmax(z, 5.001), M);}
    }
}

// Stellar mass
double Ms_from_bathtub(double z, double M) {
    if (z >= 30.0) {return 1.0e-9;}
    else{
      if (M >= 1e16) {return get_Ms(fmax(z, 5.001), 9.999e15);}
      if (M <= 1e7) {return get_Ms(fmax(z, 5.001), 1.001e7);}
      else {return get_Ms(fmax(z, 5.001), M);}
    }
}

// Metallicity (absolute)
double Z_from_bathtub(double z, double M) {
    if (z >= 30.0) {return 1.0e-9;}
    else{
      if (M >= 1e16) {return get_Z(fmax(z, 5.001), 9.999e15);}
      if (M <= 1e7) {return get_Z(fmax(z, 5.001), 1.001e7);}
      else {return get_Z(fmax(z, 5.001), M);}
    }
}

// SFE
double fstar_from_bathtub(double z, double M) {
    if (z >= 30.0) {return 1.0e-9;}
    else{
      if (M >= 1e16) {return get_fstar(fmax(z, 5.001), 9.999e15);}
      if (M <= 1e7) {return get_fstar(fmax(z, 5.001), 1.001e7);}
      else {return get_fstar(fmax(z, 5.001), M);}
    }
}

// Line luminosity (from Cloudy)
double line_luminosity(double Z, double logU, int flag) {
    if (flag == 0) {
      // Don't go beyond the [1e-5, 4e-2] interpolation range
      return get_lya_luminosity(fmax(fmin(Z, 0.04), 0.00001), logU);
    }
    else if (flag == 1) {
      return get_ha_luminosity(fmax(fmin(Z, 0.04), 0.00001), logU);
    }
    else if (flag == 2) {
      return get_hb_luminosity(fmax(fmin(Z, 0.04), 0.00001), logU);
    }
    else if (flag == 3) {
      return get_he2_luminosity(fmax(fmin(Z, 0.04), 0.00001), logU);
    }
    else if (flag == 4) {
      return get_o3_luminosity(fmax(fmin(Z, 0.04), 0.00001), logU);
    }
    else if (flag == 5) {
      return get_o2s_luminosity(fmax(fmin(Z, 0.04), 0.00001), logU);
    }
    else if (flag == 6) {
      return get_o2l_luminosity(fmax(fmin(Z, 0.04), 0.00001), logU);
    }
    else if (flag == 7) {
      return get_c2_luminosity(fmax(fmin(Z, 0.04), 0.00001), logU);
    }
    else if (flag == 8) {
      return get_co10_luminosity(fmax(fmin(Z, 0.04), 0.00001), logU);
    }
    else if (flag == 9) {
      return get_o3fir_luminosity(fmax(fmin(Z, 0.04), 0.00001), logU);
    }
    else {
      fprintf(stderr, "No such line, please check!\n");
      return 0;
    }
}

// Luminosity-halo mass relation
double LzM(double z, double M, double logU, int flag) {
    // if flag < 7, optical/UV lines are scaled from the SFR
    if (flag < 7) {
      return line_luminosity(Z_from_bathtub(z, M), logU, flag) * sfr_from_bathtub(z, M);
      //return Z_from_bathtub(z, M) * sfr_from_bathtub(z, M);
    }
    else {
      return line_luminosity(Z_from_bathtub(z, M), logU, flag) * Mg_from_bathtub(z, M)/4.7e4;
    }
}

// Nion
double Nion_LIMFAST(double z, double M) {
    float metal_res;
    metal_res = fmax(fmin(Z_from_bathtub(z, M), 0.04), 0.00001);
    return get_nion(metal_res);
}


sources test() {
    return setSources(fescDef, Nion_LIMFAST,
                      fstar_from_bathtub, sfr_from_bathtub, Mg_from_bathtub, Ms_from_bathtub, Z_from_bathtub, LzM,
                      fxDef, minMassATC, noPopIII);
}

// set the sources here
sources defaultSources()
{
    sources s;
    s = test();
    return s;
}

#endif
