#ifndef _PS_
#define _PS_

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#include "../Parameter_files/COSMOLOGY.H"
#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/SOURCES.H"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include "cosmo_progs.c"
#include "misc.c"

/* New in v1.1 */
#define ERFC_NPTS (int) 75
#define ERFC_PARAM_DELTA (float) 0.1
static double log_erfc_table[ERFC_NPTS], erfc_params[ERFC_NPTS];
static gsl_interp_accel *erfc_acc;
static gsl_spline *erfc_spline;

#define NR_END 1
#define FREE_ARG char*

#define MM 7
#define NSTACK 50

#define FUNC(x,y,z,xx,yy) ((*func)(x,y,z,xx,yy))
#define FUNC2(x1,x2,x3,x4,x5,x6) ((*func)(x1,x2,x3,x4,x5,x6))
#define EPS2 3.0e-11

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

#define NGaussLegendre 40  //defines the number of points in the Gauss-Legendre quadrature integration

#define SPLINE_NPTS (int) 250
#define NGLhigh 100
#define NGLlow 100

#define Nhigh 200
#define Nlow 100
#define NMass 2000

/* New in v1.4 - part 1 of 4 start */
#define NSFR_high 20 //for quicker run, set to 20, fiducial 200
#define NSFR_low 20 //for quicker run, set to 20, fiducial 250
#define NGL_SFR 100
/* Number of interpolation points for the interpolation table for z''
                This is the same parameter in 21CMMC */
#define zpp_interp_points (int) (30) //for quicker run, set to 30, fiducial 300
/* New in v1.4 - part 1 of 4 end */

static gsl_interp_accel *Fcoll_spline_acc;
static gsl_spline *Fcoll_spline;

/* New in v1.4 - part 2 of 4: begin */
static double log10_overdense_spline_SFR[NSFR_low], log10_Fcoll_spline_SFR[NSFR_low], log10_zeta_spline_SFR[NSFR_low];
static double log10_overdense_spline_SFRPOP2[NSFR_low], log10_Fcoll_spline_SFRPOP2[NSFR_low]; //-GS
static double log10_overdense_spline_gas[NSFR_low], log10_Fcoll_spline_gas[NSFR_low]; //-GS
static double log10_overdense_spline_met[NSFR_low], log10_Fcoll_spline_met[NSFR_low]; //-GS
static double log10_overdense_spline_LINELUMPOP2A[NSFR_low], log10_Fcoll_spline_LINELUMPOP2A[NSFR_low]; //-GS
static double log10_overdense_spline_LINELUMPOP2B[NSFR_low], log10_Fcoll_spline_LINELUMPOP2B[NSFR_low]; //-GS
static double log10_overdense_spline_LINELUMPOP2C[NSFR_low], log10_Fcoll_spline_LINELUMPOP2C[NSFR_low]; //-GS
static double log10_overdense_spline_LINELUMPOP2D[NSFR_low], log10_Fcoll_spline_LINELUMPOP2D[NSFR_low]; //-GS
static double log10_overdense_spline_LINELUMPOP2E[NSFR_low], log10_Fcoll_spline_LINELUMPOP2E[NSFR_low]; //-GS
static double log10_overdense_spline_LINELUMPOP2F[NSFR_low], log10_Fcoll_spline_LINELUMPOP2F[NSFR_low]; //-GS
static double log10_overdense_spline_LINELUMPOP2G[NSFR_low], log10_Fcoll_spline_LINELUMPOP2G[NSFR_low]; //-GS
static gsl_interp_accel *FcollLow_spline_acc;
static gsl_interp_accel *FcollLow_spline_acc_SFRPOP2; //-GS
static gsl_interp_accel *FcollLow_spline_acc_gas; //-GS
static gsl_interp_accel *FcollLow_spline_acc_met; //-GS
static gsl_interp_accel *FcollLow_spline_acc_LINELUMPOP2A; //-GS
static gsl_interp_accel *FcollLow_spline_acc_LINELUMPOP2B; //-GS
static gsl_interp_accel *FcollLow_spline_acc_LINELUMPOP2C; //-GS
static gsl_interp_accel *FcollLow_spline_acc_LINELUMPOP2D; //-GS
static gsl_interp_accel *FcollLow_spline_acc_LINELUMPOP2E; //-GS
static gsl_interp_accel *FcollLow_spline_acc_LINELUMPOP2F; //-GS
static gsl_interp_accel *FcollLow_spline_acc_LINELUMPOP2G; //-GS
static gsl_interp_accel *zetaLow_spline_acc;
static gsl_spline *FcollLow_spline;
static gsl_spline *FcollLow_spline_SFRPOP2; //-GS
static gsl_spline *FcollLow_spline_gas; //-GS
static gsl_spline *FcollLow_spline_met; //-GS
static gsl_spline *FcollLow_spline_LINELUMPOP2A; //-GS
static gsl_spline *FcollLow_spline_LINELUMPOP2B; //-GS
static gsl_spline *FcollLow_spline_LINELUMPOP2C; //-GS
static gsl_spline *FcollLow_spline_LINELUMPOP2D; //-GS
static gsl_spline *FcollLow_spline_LINELUMPOP2E; //-GS
static gsl_spline *FcollLow_spline_LINELUMPOP2F; //-GS
static gsl_spline *FcollLow_spline_LINELUMPOP2G; //-GS
static gsl_spline *zetaLow_spline;

void initialiseGL_FcollSFR(int n, float M_TURN, float M_Max, float z);
void FcollSpline_SFR(float Overdensity, float *splined_value);
void FcollSpline_SFRPOP2(float Overdensity, float *splined_value);
void FcollSpline_gas(float Overdensity, float *splined_value);
void FcollSpline_met(float Overdensity, float *splined_value);
void FcollSpline_LINELUMPOP2A(float Overdensity, float *splined_value);
void FcollSpline_LINELUMPOP2B(float Overdensity, float *splined_value);
void FcollSpline_LINELUMPOP2C(float Overdensity, float *splined_value);
void FcollSpline_LINELUMPOP2D(float Overdensity, float *splined_value);
void FcollSpline_LINELUMPOP2E(float Overdensity, float *splined_value);
void FcollSpline_LINELUMPOP2F(float Overdensity, float *splined_value);
void FcollSpline_LINELUMPOP2G(float Overdensity, float *splined_value);
void initialiseFcollSFR_spline(float z, float Mmax, float MassTurnover, float Alpha_star, float Alpha_esc, float Fstar10, float Fesc10, float Mlim_Fstar, float Mlim_Fesc);
void initialiseFcollSFRPOP2_spline(float z, float Mmax, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar);
void initialiseFcollgas_spline(float z, float Mmax, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar);
void initialiseFcollmet_spline(float z, float Mmax, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar);
void initialiseFcollLINELUMPOP2A_spline(float z, float Mmax, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar, int LineID);
void initialiseFcollLINELUMPOP2B_spline(float z, float Mmax, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar, int LineID);
void initialiseFcollLINELUMPOP2C_spline(float z, float Mmax, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar, int LineID);
void initialiseFcollLINELUMPOP2D_spline(float z, float Mmax, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar, int LineID);
void initialiseFcollLINELUMPOP2E_spline(float z, float Mmax, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar, int LineID);
void initialiseFcollLINELUMPOP2F_spline(float z, float Mmax, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar, int LineID);
void initialiseFcollLINELUMPOP2G_spline(float z, float Mmax, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar, int LineID);

void initialiseGL_zetaSFR(int n, float M_Min, float M_Max);
void zetaSpline_SFR(float Overdensity, float *splined_value);
void initialisezetaSFR_spline(float z, float Mmax, float MassTurnover, float Alpha_star, float Alpha_esc, float Fstar10, float Fesc10, float Mlim_Fstar, float Mlim_Fesc);

float Mass_limit (float logM, float PL, float FRAC);
void bisection(float *x, float xlow, float xup, int *iter);
float Mass_limit_bisection(float Mmin, float Mmax, float PL, float FRAC);

    /* For Ts.c: being under development below this line */
static double z_val[zpp_interp_points],Fcollz_val[zpp_interp_points]; // For Ts.c
static double z_X_val[zpp_interp_points], FcollzX_val[zpp_interp_points], FcollzLya_val[zpp_interp_points];
static gsl_interp_accel *Fcollz_spline_acc;
static gsl_spline *Fcollz_spline;
static gsl_interp_accel *FcollzX_spline_acc;
static gsl_spline *FcollzX_spline;
static gsl_interp_accel *FcollzLya_spline_acc;
static gsl_spline *FcollzLya_spline;
void initialise_FgtrM_st_SFR_spline(int Nbin, float zmin, float zmax, float MassTurn, float Alpha_star, float Alpha_esc, float Fstar10, float Fesc10);
void FgtrM_st_SFR_z(float z, float *splined_value);
void initialise_Xray_FgtrM_st_SFR_spline(int Nbin, float zmin, float zmax, float MassTurn, float Alpha_star, float Fstar10);
void initialise_Lya_FgtrM_st_SFR_spline(int Nbin, float zmin, float zmax, float MassTurn, float Alpha_star, float Fstar10);
void FgtrM_st_SFR_Xray_z(float z, float *splined_value);
void FgtrM_st_SFR_Lya_z(float z, float *splined_value);
float *zpp_table;
double *log10_overdense_low_table, **log10_FcollzXray_SFR_low_table, **log10_FcollzLya_SFR_low_table;
float *Overdense_high_table, **FcollzXray_SFR_high_table, **FcollzLya_SFR_high_table;
void initialise_Xray_Fcollz_SFR_Conditional_table(int Nsteps_zp, int Nfilter, float z[], double R[], float MassTurnover, float Alpha_star, float Fstar10);
void initialise_Lya_Fcollz_SFR_Conditional_table(int Nsteps_zp, int Nfilter, float z[], double R[], float MassTurnover, float Alpha_star, float Fstar10);
//void initialise_Xray_Fcollz_SFR_Conditional(int R_ct, int zp_int1, int zp_int2);
void free_interpolation();
    /* For Ts.c: being under development above this line */


struct parameters_gsl_SFR_int_{
  double z_obs;
  double Mdrop;
  double pl_star;
  double pl_esc;
	double frac_star;
	double frac_esc;
	double LimitMass_Fstar;
	double LimitMass_Fesc;
  int lid;
};

struct parameters_gsl_SFR_con_int_{
  double z_obs;
  double Mval;
  double delta1;
  double delta2;
  double Mdrop;
  double pl_star;
  double pl_esc;
	double frac_star;
	double frac_esc;
	double LimitMass_Fstar;
	double LimitMass_Fesc;
  int lid;
};
/* New in v1.4 - part 2 of 4: end */


unsigned long *lvector(long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);

float *vector(long nl, long nh);
void free_vector(float *v, long nl, long nh);

void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);

void gauleg(float x1, float x2, float x[], float w[], int n);

double FgtrlnM_general(double lnM, void *params);
double FgtrM_general(float z, float M1, float M_Max, float M2, float MFeedback, float alpha, float delta1, float delta2);

float FgtrConditionalM_second(float z, float M1, float M2, float MFeedback, float alpha, float delta1, float delta2);
float dNdM_conditional_second(float z, float M1, float M2, float delta1, float delta2);

float *Mass_Spline, *Sigma_Spline, *dSigmadm_Spline, *second_derivs_sigma, *second_derivs_dsigma;

void initialiseSplinedSigmaM(float M_Min, float M_Max);
void initialiseGL_Fcoll(int n_low, int n_high, float M_Min, float M_Max);
void initialiseGL_FcollDblPl(int n_low, int n_high, float M_Min, float M_feedback, float M_Max);
void initialiseFcoll_spline(float z, float Mmin, float Mmax, float Mval, float MFeedback, float alphapl);

double dFdlnM_st_PL (double lnM, void *params);
double FgtrM_st_PL(double z, double Mmin, double MFeedback, double alpha_pl);

/* New in v1.4 - part 3 of 4: start */
float *Overdense_spline_SFR,*Fcoll_spline_SFR,*zeta_spline_SFR,*second_derivs_SFR,*zeta_second_derivs_SFR;
float *Overdense_spline_SFRPOP2,*Fcoll_spline_SFRPOP2,*second_derivs_SFRPOP2; //-GS
float *Overdense_spline_gas,*Fcoll_spline_gas,*second_derivs_gas; //-GS
float *Overdense_spline_met,*Fcoll_spline_met,*second_derivs_met; //-GS
float *Overdense_spline_LINELUMPOP2A,*Fcoll_spline_LINELUMPOP2A,*second_derivs_LINELUMPOP2A; //-GS
float *Overdense_spline_LINELUMPOP2B,*Fcoll_spline_LINELUMPOP2B,*second_derivs_LINELUMPOP2B; //-GS
float *Overdense_spline_LINELUMPOP2C,*Fcoll_spline_LINELUMPOP2C,*second_derivs_LINELUMPOP2C; //-GS
float *Overdense_spline_LINELUMPOP2D,*Fcoll_spline_LINELUMPOP2D,*second_derivs_LINELUMPOP2D; //-GS
float *Overdense_spline_LINELUMPOP2E,*Fcoll_spline_LINELUMPOP2E,*second_derivs_LINELUMPOP2E; //-GS
float *Overdense_spline_LINELUMPOP2F,*Fcoll_spline_LINELUMPOP2F,*second_derivs_LINELUMPOP2F; //-GS
float *Overdense_spline_LINELUMPOP2G,*Fcoll_spline_LINELUMPOP2G,*second_derivs_LINELUMPOP2G; //-GS
float *xi_SFR,*wi_SFR;
float *xi_SFR_zeta,*wi_SFR_zeta;

float FgtrConditionallnM_GL_SFR(float M, struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con);
float FgtrConditionallnM_GL_SFRPOP2(float lnM, struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con); //-GS
float FgtrConditionallnM_GL_gas(float lnM, struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con); //-GS
float FgtrConditionallnM_GL_met(float lnM, struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con); //-GS
float FgtrConditionallnM_GL_LINELUMPOP2A(float lnM, struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con); //-GS
float FgtrConditionallnM_GL_LINELUMPOP2B(float lnM, struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con); //-GS
float FgtrConditionallnM_GL_LINELUMPOP2C(float lnM, struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con); //-GS
float FgtrConditionallnM_GL_LINELUMPOP2D(float lnM, struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con); //-GS
float FgtrConditionallnM_GL_LINELUMPOP2E(float lnM, struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con); //-GS
float FgtrConditionallnM_GL_LINELUMPOP2F(float lnM, struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con); //-GS
float FgtrConditionallnM_GL_LINELUMPOP2G(float lnM, struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con); //-GS
float GaussLegendreQuad_FcollSFR(int n, float z, float M2, float delta1, float delta2, float MassTurnover, float Alpha_star, float Alpha_esc, float Fstar10, float Fesc10, float Mlim_Fstar, float Mlim_Fesc);
float GaussLegendreQuad_FcollSFRPOP2(int n, float z, float M2, float delta1, float delta2, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar); // -GS
float GaussLegendreQuad_Fcollgas(int n, float z, float M2, float delta1, float delta2, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar); // -GS
float GaussLegendreQuad_Fcollmet(int n, float z, float M2, float delta1, float delta2, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar); // -GS
float GaussLegendreQuad_FcollLINELUMPOP2A(int n, float z, float M2, float delta1, float delta2, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar, int LineID); // -GS
float GaussLegendreQuad_FcollLINELUMPOP2B(int n, float z, float M2, float delta1, float delta2, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar, int LineID); // -GS
float GaussLegendreQuad_FcollLINELUMPOP2C(int n, float z, float M2, float delta1, float delta2, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar, int LineID); // -GS
float GaussLegendreQuad_FcollLINELUMPOP2D(int n, float z, float M2, float delta1, float delta2, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar, int LineID); // -GS
float GaussLegendreQuad_FcollLINELUMPOP2E(int n, float z, float M2, float delta1, float delta2, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar, int LineID); // -GS
float GaussLegendreQuad_FcollLINELUMPOP2F(int n, float z, float M2, float delta1, float delta2, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar, int LineID); // -GS
float GaussLegendreQuad_FcollLINELUMPOP2G(int n, float z, float M2, float delta1, float delta2, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar, int LineID); // -GS

double dFgtrConditionallnM_SFR(double lnM, void *params);
double dFgtrConditionallnM_SFRPOP2(double lnM, void *params);
double dFgtrConditionallnM_gas(double lnM, void *params);
double dFgtrConditionallnM_met(double lnM, void *params);
double dFgtrConditionallnM_LINELUMPOP2A(double lnM, void *params);
double dFgtrConditionallnM_LINELUMPOP2B(double lnM, void *params);
double dFgtrConditionallnM_LINELUMPOP2C(double lnM, void *params);
double dFgtrConditionallnM_LINELUMPOP2D(double lnM, void *params);
double dFgtrConditionallnM_LINELUMPOP2E(double lnM, void *params);
double dFgtrConditionallnM_LINELUMPOP2F(double lnM, void *params);
double dFgtrConditionallnM_LINELUMPOP2G(double lnM, void *params);
double FgtrConditionalM_SFR(double z, double M1, double M2, double delta1, double delta2, double MassTurnover, double Alpha_star, double Alpha_esc, double Fstar10, double Fesc10, double Mlim_Fstar, double Mlim_Fesc);
double FgtrConditionalM_SFRPOP2(double z, double M1, double M2, double delta1, double delta2, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar); // -GS
double FgtrConditionalM_gas(double z, double M1, double M2, double delta1, double delta2, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar); // -GS
double FgtrConditionalM_met(double z, double M1, double M2, double delta1, double delta2, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar); // -GS
double FgtrConditionalM_LINELUMPOP2A(double z, double M1, double M2, double delta1, double delta2, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar, int LineID); // -GS
double FgtrConditionalM_LINELUMPOP2B(double z, double M1, double M2, double delta1, double delta2, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar, int LineID); // -GS
double FgtrConditionalM_LINELUMPOP2C(double z, double M1, double M2, double delta1, double delta2, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar, int LineID); // -GS
double FgtrConditionalM_LINELUMPOP2D(double z, double M1, double M2, double delta1, double delta2, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar, int LineID); // -GS
double FgtrConditionalM_LINELUMPOP2E(double z, double M1, double M2, double delta1, double delta2, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar, int LineID); // -GS
double FgtrConditionalM_LINELUMPOP2F(double z, double M1, double M2, double delta1, double delta2, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar, int LineID); // -GS
double FgtrConditionalM_LINELUMPOP2G(double z, double M1, double M2, double delta1, double delta2, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar, int LineID); // -GS

double dFdlnM_st_SFR (double lnM, void *params);
double dFdlnM_st_SFRPOP2 (double lnM, void *params); //-GS
double dFdlnM_st_gas (double lnM, void *params); //-GS
double dFdlnM_st_met (double lnM, void *params); //-GS
double dFdlnM_st_ZNUMEPOP2 (double lnM, void *params); //-GS, added for tracking SED evolution
double dFdlnM_st_ZDENOPOP2 (double lnM, void *params); //-GS, added for tracking SED evolution
double dFdlnM_st_LINELUMPOP2A (double lnM, void *params); //-GS
double dFdlnM_st_LINELUMPOP2B (double lnM, void *params); //-GS
double dFdlnM_st_LINELUMPOP2C (double lnM, void *params); //-GS
double dFdlnM_st_LINELUMPOP2D (double lnM, void *params); //-GS
double dFdlnM_st_LINELUMPOP2E (double lnM, void *params); //-GS
double dFdlnM_st_LINELUMPOP2F (double lnM, void *params); //-GS
double dFdlnM_st_LINELUMPOP2G (double lnM, void *params); //-GS
double FgtrM_st_SFR (double z, double MassTurnover, double Alpha_star, double Alpha_esc, double Fstar10, double Fesc10, double Mlim_Fstar, double Mlim_Fesc);
double FgtrM_st_SFRPOP2(double z, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar); // -GS
double FgtrM_st_gas(double z, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar); // -GS
double FgtrM_st_met(double z, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar); // -GS
double FgtrM_st_ZNUMEPOP2(double z, double MassTurnover, double Alpha_star, double Fstar10); // -GS, added for tracking SED evolution
double FgtrM_st_ZDENOPOP2(double z, double MassTurnover, double Alpha_star, double Fstar10); // -GS, added for tracking SED evolution
double FgtrM_st_LINELUMPOP2A(double z, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar, int LineID); // -GS
double FgtrM_st_LINELUMPOP2B(double z, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar, int LineID); // -GS
double FgtrM_st_LINELUMPOP2C(double z, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar, int LineID); // -GS
double FgtrM_st_LINELUMPOP2D(double z, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar, int LineID); // -GS
double FgtrM_st_LINELUMPOP2E(double z, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar, int LineID); // -GS
double FgtrM_st_LINELUMPOP2F(double z, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar, int LineID); // -GS
double FgtrM_st_LINELUMPOP2G(double z, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar, int LineID); // -GS
/* New in v1.4 - part 3 of 4: end */

double sigma_norm, R, theta_cmb, omhh, z_equality, y_d, sound_horizon, alpha_nu, f_nu, f_baryon, beta_c, d2fact, R_CUTOFF, DEL_CURR, SIG_CURR;


/*****     FUNCTION PROTOTYPES     *****/
double init_ps(); /* initialize global variables, MUST CALL THIS FIRST!!! returns R_CUTOFF */
void free_ps(); /* deallocates the gsl structures from init_ps */
float mean_SFRD(double z); // returns the mean star formation rate density at z in M_sun yr^-1 Mpc^-3
double splined_erfc(double); /* returns erfc for x>=0, using cubic spline in logy-x space */
double deltolindel(float del, float z); /* converts a non-linear overdensity, del, at z to a linear overdensity at z=0 */
double lindeltodel(float lindel, float z); /* converts a linear overdensity, del, at z=0 to a non-linear overdensity at redshift z */
double power_in_k(double k); /* Returns the value of the linear power spectrum density (i.e. <|delta_k|^2>/V) at a given k mode at z=0 */
double RtoM(double); /* R in Mpc, M in Msun */
double MtoR(double); /* R in Mpc, M in Msun */
double M_J_WDM(); /* returns the "effective Jeans mass" corresponding to the gas analog of WDM ; eq. 10 in BHO 2001 */
double sheth_delc(double del, double sig);
double dNdM_st(double z, double M);
double dNdM(double z, double M);
double dnbiasdM(double M, float z, double M_o, float del_o); /* dnbiasdM */
double FgtrM(double z, double M);  //calculates the fraction of mass contained in haloes with mass > M at redshift z
double FgtrM_st(double z, double M);  //calculates the fraction of mass contained in haloes with mass > M at redshift z, with Sheth-Tormen correction
double FgtrM_bias(double z, double M, double del_bias, double sig_bias);  //calculates the fraction of mass contained in haloes with mass > M at redshift z, in regions with a linear overdensity of del_bias, and standard deviation sig_bias
double sigmaparam_FgtrM_bias(float z, float sigsmallR, float del_bias, float sig_bias);/* Uses sigma parameters instead of Mass for scale */

double FgtrM_bias_BL08(double z, double M, double del_bias, double sig_bias); // as above, but this version uses the hybrid perscription of Barkana & Loeb 2004 (specifically the separate integral version of eq. 2 in Barkana & Loeb 2008)
double dicke(double z); //calculates the dicke growth function at redshift z
double ddickedz(double z); /* Redshift derivative of the growth function at z */
double ddickedt(double z); /* Time derivative of the growth function at z */
double sigma_z0(double M); //calculates sigma at z=0 (no dicke)
double dsigmasqdm_z0(double M); //calculates d(sigma^2)/dm at z=0 (i.e. does not include dicke growth)
double TFmdm(double k); //Eisenstien & Hu power spectrum transfer function
void TFset_parameters();
float get_R_c();  // returns R_CUTOFF
double get_M_min_ion(float z);
/***************************************/



/* Returns the minimum source mass for ionizing sources, according to user specifications */
double get_M_min_ion(float z){
  double MMIN;

  MMIN = M_TURNOVER;

  // check for WDM
  if (P_CUTOFF && ( MMIN < M_J_WDM()))
    MMIN = M_J_WDM();
  return MMIN;
}

/* Returns the minimum source mass for x-ray sources, according to user specifications */
double get_M_min_xray(float z){
  return get_M_min_ion(z);
}


/* returns the "effective Jeans mass" in Msun
   corresponding to the gas analog of WDM ; eq. 10 in Barkana+ 2001 */
double M_J_WDM(){
  double z_eq, fudge=60;
  if (!P_CUTOFF)
    return 0;
  z_eq = 3600*(OMm-OMb)*hlittle*hlittle/0.15;
  return fudge*3.06e8 * (1.5/g_x) * sqrt((OMm-OMb)*hlittle*hlittle/0.15) * pow(M_WDM, -4) * pow(z_eq/3000.0, 1.5);
}


/* converts a non-linear overdensity, del, at z to a linear overdensity at z=0 */
double deltolindel(float del, float z){
  float onepdel = 1.0+del;
  return ( 1.68647 - 1.35*pow(onepdel,-2/3.0) + 0.78785*pow(onepdel,-0.58661) - 1.12431*pow(onepdel,-0.5) )/dicke(z);
}

/* converts a linear overdensity, del, at z=0 to a non-linear overdensity at redshift z */
double lindeltodel(float lindel, float z){
  float prev_lindelguess, delcrit, delguess;
  float lindelguess, delmin, delmax, epsilon = 1.0e-7;

  // set the critical density corresponding to virialization
  // this will be maximum allowed del
  delcrit = Deltac_nonlinear(z)*rho_critz(z)/(OMm*RHOcrit*pow(1+z, 3)) - 1;

  delmin = -1;
  delmax = 500;
  prev_lindelguess = -1e10;
  while (1){
    delguess = 0.5*(delmax+delmin);
    lindelguess = deltolindel(delguess, z);
    //fprintf(stderr, "%e\t%e\n", delmin, delmax);
	      //    fprintf(stderr, "%e\t%e\t%e\n\n", delguess, lindelguess, lindel);
    if ((fabs((lindelguess-lindel)/lindel) < epsilon ) ||
	(fabs(lindelguess-lindel) < epsilon ) ||
	(fabs(prev_lindelguess - lindelguess) < TINY ))// close enough, or resolution loop
      return delguess;

    if (lindelguess > lindel)
      delmax = delguess;
    else
      delmin = delguess;

    // check if we are above delcrit (see above)
    if (delmin > delcrit){
      //      printf("exced max at lindel=%e\n", lindel);
      return delcrit;
    }

    prev_lindelguess = lindelguess;
  }
}


/* R in Mpc, M in Msun */
double RtoM(double R){
  // set M according to M<->R conversion defined by the filter type in ../Parameter_files/COSMOLOGY.H
  if (FILTER == 0) //top hat M = (4/3) PI <rho> R^3
    return (4.0/3.0)*PI*pow(R,3)*(OMm*RHOcrit);
  else if (FILTER == 1) //gaussian: M = (2PI)^1.5 <rho> R^3
    return pow(2*PI, 1.5) * OMm*RHOcrit * pow(R, 3);
  else // filter not defined
    fprintf(stderr, "No such filter = %i.\nResults are bogus.\n", FILTER);
  return -1;
}

/* R in Mpc, M in Msun */
double MtoR(double M){
  // set R according to M<->R conversion defined by the filter type in ../Parameter_files/COSMOLOGY.H
  if (FILTER == 0) //top hat M = (4/3) PI <rho> R^3
    return pow(3*M/(4*PI*OMm*RHOcrit), 1.0/3.0);
  else if (FILTER == 1) //gaussian: M = (2PI)^1.5 <rho> R^3
    return pow( M/(pow(2*PI, 1.5) * OMm * RHOcrit), 1.0/3.0 );
  else // filter not defined
    fprintf(stderr, "No such filter = %i.\nResults are bogus.\n", FILTER);
  return -1;
}


/* equation (5) from jenkis et al. (2001) */
double f_jenkins(float del, double sigsq){
  if (del < 0){  fprintf(stderr, "ERROR:  In function f_jenkins del_o must be less than del_1 = del_crit/dicke(z)!\nAborting...\n"); return 0; }

  //  fprintf(stderr, "%f\t%f\n", del, sqrt(sigsq));
  return sqrt(2/PI) * del/sqrt(sigsq) * pow(E, -0.5*del*del/sigsq);
}


float get_R_c(){
  return R_CUTOFF;
}

/* sheth correction to delta crit */
double sheth_delc(double del, double sig){
  return sqrt(SHETH_a)*del*(1 + SHETH_b*pow(sig*sig/(SHETH_a*del*del), SHETH_c));
}


/* dnbiasdM */
double dnbiasdM(double M, float z, double M_o, float del_o){
  double sigsq, del, sig_one, sig_o;

  if ((M_o-M) < TINY){
    fprintf(stderr, "WARNING:  In function dnbiasdM: M must be less than M_o!\nAborting...\n");
    return -1;
  }
  del = Deltac/dicke(z) - del_o;
  if (del < 0){  fprintf(stderr, "ERROR:  In function dnbiasdM: del_o must be less than del_1 = del_crit/dicke(z)!\nAborting...\n"); return 0; }
  sig_o = sigma_z0(M_o);
  sig_one = sigma_z0(M);
  sigsq = sig_one*sig_one - sig_o*sig_o;
  return -(RHOcrit*OMm)/M /sqrt(2*PI) *del*pow(sigsq,-1.5)*pow(E, -0.5*del*del/sigsq)*dsigmasqdm_z0(M);
}


/*
  FUNCTION dNdM(z, M)
  Computes the Press_schechter mass function with Sheth-Torman correction for ellipsoidal collapse at
  redshift z, and dark matter halo mass M (in solar masses).

  The return value is the number density per unit mass of halos in the mass range M to M+dM in units of:
  comoving Mpc^-3 Msun^-1

  Reference: Sheth, Mo, Torman 2001
*/
double dNdM_st(double z, double M){
  double sigma, dsigmadm, nuhat, dicke_growth;

  dicke_growth = dicke(z);
  sigma = sigma_z0(M) * dicke_growth;
  dsigmadm = dsigmasqdm_z0(M) * dicke_growth*dicke_growth/(2.0*sigma);
  nuhat = sqrt(SHETH_a) * Deltac / sigma;

  return (-OMm*RHOcrit/M) * (dsigmadm/sigma) * sqrt(2/PI)*SHETH_A * (1+ pow(nuhat, -2*SHETH_p)) * nuhat * pow(E, -nuhat*nuhat/2.0);
}



/*
  FUNCTION dNdM(z, M)
  Computes the Press_schechter mass function at
  redshift z, and dark matter halo mass M (in solar masses).

  The return value is the number density per unit mass of halos in the mass range M to M+dM in units of:
  comoving Mpc^-3 Msun^-1

  Reference: Padmanabhan, pg. 214
*/
double dNdM(double z, double M){
  double sigma, dsigmadm, dicke_growth;

  dicke_growth = dicke(z);
  sigma = sigma_z0(M) * dicke_growth;
  dsigmadm = dsigmasqdm_z0(M) * (dicke_growth*dicke_growth/(2*sigma));

  return (-OMm*RHOcrit/M) * sqrt(2/PI) * (Deltac/(sigma*sigma)) * dsigmadm * pow(E, -(Deltac*Deltac)/(2*sigma*sigma));
}


/*
  FUNCTION FgtrM_st(z, M)
  Computes the fraction of mass contained in haloes with mass > M at redshift z
  Uses Sheth-Torman correction
*/
double dFdlnM_st (double lnM, void *params){
  double z = *(double *)params;
  double M = exp(lnM);
  return dNdM_st(z, M) * M * M;
}
double FgtrM_st(double z, double M){
  double result, error, lower_limit, upper_limit;
  gsl_function F;
  double rel_tol  = 0.001; //<- relative tolerance
  gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

  F.function = &dFdlnM_st;
  F.params = &z;
  lower_limit = log(M);
  upper_limit = log(FMAX(1e16, M*100));

  gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
		       1000, GSL_INTEG_GAUSS61, w, &result, &error);
  gsl_integration_workspace_free (w);

  return result / (OMm*RHOcrit);
}


/*
  FUNCTION FgtrM(z, M)
  Computes the fraction of mass contained in haloes with mass > M at redshift z
*/
double FgtrM(double z, double M){
  double del, sig;

  del = Deltac/dicke(z); //regular spherical collapse delta
  sig = sigma_z0(M);

  return splined_erfc(del / (sqrt(2)*sig));
}


/*
  calculates the fraction of mass contained in haloes with mass > M at redshift z, in regions with a linear overdensity of del_bias, and standard deviation sig_bias
 */
double FgtrM_bias(double z, double M, double del_bias, double sig_bias){
  double del, sig, sigsmallR;

  sigsmallR = sigma_z0(M);

  if (!(sig_bias < sigsmallR)){ // biased region is smaller that halo!
    fprintf(stderr, "FgtrM_bias: Biased region is smaller than halo!\nResult is bogus.\n");
    return 0;
  }

  del = Deltac/dicke(z) - del_bias;
  sig = sqrt(sigsmallR*sigsmallR - sig_bias*sig_bias);

   return splined_erfc(del / (sqrt(2)*sig));
}


/* Uses sigma parameters instead of Mass for scale */
double sigmaparam_FgtrM_bias(float z, float sigsmallR, float del_bias, float sig_bias){
  double del, sig;

  if (!(sig_bias < sigsmallR)){ // biased region is smaller that halo!
    fprintf(stderr, "local_FgtrM_bias: Biased region is smaller than halo!\nResult is bogus.\n");
printf("In function: zpp = %.4f, sigma_Tmin = %.4e, sigma_atR = %.4e, delta = %.5f\n",z,sigsmallR,sig_bias,del_bias);
    return 0;
  }

  del = Deltac/dicke(z) - del_bias;
  sig = sqrt(sigsmallR*sigsmallR - sig_bias*sig_bias);

  return splined_erfc(del / (sqrt(2)*sig));
}


/*
  Calculates the fraction of mass contained in haloes with mass > M at redshift z, in regions with a linear overdensity of del_bias, and standard deviation sig_bias.
  This version uses the hybrid perscription of Barkana & Loeb 2004 (specifically the separate
  integral version of eq. 2 in Barkana & Loeb 2008)
 */
double FgtrM_bias_BL08(double z, double M, double del_bias, double sig_bias){
  return FgtrM_st(z, M) / FgtrM(z, M) * FgtrM_bias(z, M, del_bias, sig_bias);
}


/*
  FUNCTION dicke(z)
  Computes the dicke growth function at redshift z, i.e. the z dependance part of sigma

  References: Peebles, "Large-Scale...", pg.53 (eq. 11.16). Includes omega<=1
  Nonzero Lambda case from Liddle et al, astro-ph/9512102, eqs. 6-8.
  and quintessence case from Wang et al, astro-ph/9804015

  Normalized to dicke(z=0)=1
*/
double dicke(double z){
  double omegaM_z, dick_z, dick_0, x, x_0;
  double tiny = 1e-4;

  if (fabs(OMm-1.0) < tiny){ //OMm = 1 (Einstein de-Sitter)
    return 1.0/(1.0+z);
  }
  else if ( (OMl > (-tiny)) && (fabs(OMl+OMm+OMr-1.0) < 0.01) && (fabs(wl+1.0) < tiny) ){
    //this is a flat, cosmological CONSTANT universe, with only lambda, matter and radiation
    //it is taken from liddle et al.
    omegaM_z = OMm*pow(1+z,3) / ( OMl + OMm*pow(1+z,3) + OMr*pow(1+z,4) );
    dick_z = 2.5*omegaM_z / ( 1.0/70.0 + omegaM_z*(209-omegaM_z)/140.0 + pow(omegaM_z, 4.0/7.0) );
    dick_0 = 2.5*OMm / ( 1.0/70.0 + OMm*(209-OMm)/140.0 + pow(OMm, 4.0/7.0) );
    return dick_z / (dick_0 * (1.0+z));
  }
  else if ( (OMtot < (1+tiny)) && (fabs(OMl) < tiny) ){ //open, zero lambda case (peebles, pg. 53)
    x_0 = 1.0/(OMm+0.0) - 1.0;
    dick_0 = 1 + 3.0/x_0 + 3*log(sqrt(1+x_0)-sqrt(x_0))*sqrt(1+x_0)/pow(x_0,1.5);
    x = fabs(1.0/(OMm+0.0) - 1.0) / (1+z);
    dick_z = 1 + 3.0/x + 3*log(sqrt(1+x)-sqrt(x))*sqrt(1+x)/pow(x,1.5);
    return dick_z/dick_0;
  }
  else if ( (OMl > (-tiny)) && (fabs(OMtot-1.0) < tiny) && (fabs(wl+1) > tiny) ){
    fprintf(stderr, "IN WANG\n");
    return -1;
  }

  fprintf(stderr, "No growth function!!! Output will be fucked up.");
  return -1;
}


/* redshift derivative of the growth function at z */
double ddicke_dz(double z){
  float dz = 1e-10;
  double omegaM_z, ddickdz, dick_0, x, x_0, domegaMdz;

   return (dicke(z+dz)-dicke(z))/dz;
}

/* Time derivative of the growth function at z */
double ddickedt(double z){
  float dz = 1e-10;
  double omegaM_z, ddickdz, dick_0, x, x_0, domegaMdz;
  double tiny = 1e-4;

   return (dicke(z+dz)-dicke(z))/dz/dtdz(z); // lazy non-analytic form getting

  if (fabs(OMm-1.0) < tiny){ //OMm = 1 (Einstein de-Sitter)
    return -pow(1+z,-2)/dtdz(z);
  }
  else if ( (OMl > (-tiny)) && (fabs(OMl+OMm+OMr-1.0) < 0.01) && (fabs(wl+1.0) < tiny) ){
    //this is a flat, cosmological CONSTANT universe, with only lambda, matter and radiation
    //it is taken from liddle et al.
    omegaM_z = OMm*pow(1+z,3) / ( OMl + OMm*pow(1+z,3) + OMr*pow(1+z,4) );
    domegaMdz = omegaM_z*3/(1+z) - OMm*pow(1+z,3)*pow(OMl + OMm*pow(1+z,3) + OMr*pow(1+z,4), -2) * (3*OMm*(1+z)*(1+z) + 4*OMr*pow(1+z,3));
    dick_0 = OMm / ( 1.0/70.0 + OMm*(209-OMm)/140.0 + pow(OMm, 4.0/7.0) );

    ddickdz = (domegaMdz/(1+z)) * (1.0/70.0*pow(omegaM_z,-2) + 1.0/140.0 + 3.0/7.0*pow(omegaM_z, -10.0/3.0)) * pow(1.0/70.0/omegaM_z + (209.0-omegaM_z)/140.0 + pow(omegaM_z, -3.0/7.0) , -2);
    ddickdz -= pow(1+z,-2)/(1.0/70.0/omegaM_z + (209.0-omegaM_z)/140.0 + pow(omegaM_z, -3.0/7.0));

    return ddickdz / dick_0 / dtdz(z);
  }

  fprintf(stderr, "No growth function!!! Output will be fucked up.");
  return -1;
}

/*
  FUNCTION sigma_z0(M)
  Returns the standard deviation of the normalized, density excess (delta(x)) field,
  smoothed on the comoving scale of M (see filter definitions for M<->R conversion).
  The sigma is evaluated at z=0, with the time evolution contained in the dicke(z) factor,
  i.e. sigma(M,z) = sigma_z0(m) * dicke(z)

  normalized so that sigma_z0(M->8/h Mpc) = SIGMA8 in ../Parameter_files/COSMOLOGY.H

  NOTE: volume is normalized to = 1, so this is equvalent to the mass standard deviation

  M is in solar masses

  References: Padmanabhan, pg. 210, eq. 5.107
*/
double dsigma_dk(double k, void *params){
  double p, w, T, gamma, q, aa, bb, cc, kR;

  // get the power spectrum.. choice of 5:
  if (POWER_SPECTRUM == 0){ // Eisenstein & Hu
    T = TFmdm(k);
    // check if we should cuttoff power spectrum according to Bode et al. 2000 transfer function
    if (P_CUTOFF) T *= pow(1 + pow(BODE_e*k*R_CUTOFF, 2*BODE_v), -BODE_n/BODE_v);
    p = pow(k, POWER_INDEX) * T * T;
  }
  else if (POWER_SPECTRUM == 1){ // BBKS
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    q = k / (hlittle*gamma);
    T = (log(1.0+2.34*q)/(2.34*q)) *
      pow( 1.0+3.89*q + pow(16.1*q, 2) + pow( 5.46*q, 3) + pow(6.71*q, 4), -0.25);
     p = pow(k, POWER_INDEX) * T * T;
  }
  else if (POWER_SPECTRUM == 2){ // Efstathiou,G., Bond,J.R., and White,S.D.M., MNRAS,258,1P (1992)
    gamma = 0.25;
    aa = 6.4/(hlittle*gamma);
    bb = 3.0/(hlittle*gamma);
    cc = 1.7/(hlittle*gamma);
    p = pow(k, POWER_INDEX) / pow( 1+pow( aa*k + pow(bb*k, 1.5) + pow(cc*k,2), 1.13), 2.0/1.13 );
  }
  else if (POWER_SPECTRUM == 3){ // Peebles, pg. 626
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    aa = 8.0 / (hlittle*gamma);
    bb = 4.7 / pow(hlittle*gamma, 2);
    p = pow(k, POWER_INDEX) / pow(1 + aa*k + bb*k*k, 2);
  }
  else if (POWER_SPECTRUM == 4){ // White, SDM and Frenk, CS, 1991, 379, 52
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    aa = 1.7/(hlittle*gamma);
    bb = 9.0/pow(hlittle*gamma, 1.5);
    cc = 1.0/pow(hlittle*gamma, 2);
    p = pow(k, POWER_INDEX) * 19400.0 / pow(1 + aa*k + bb*pow(k, 1.5) + cc*k*k, 2);
  }
  else{
    fprintf(stderr, "No such power spectrum defined: %i\nOutput is bogus.\n", POWER_SPECTRUM);
    p = 0;
  }


  // now get the value of the window function
  // NOTE: only use top hat for SIGMA8 normalization
  kR = k*R;
  if ( (FILTER == 0) || (sigma_norm < 0) ){ // top hat
    if ( (kR) < 1.0e-4 ){ w = 1.0;} // w converges to 1 as (kR) -> 0
    else { w = 3.0 * (sin(kR)/pow(kR, 3) - cos(kR)/pow(kR, 2));}
  }
  else if (FILTER == 1){ // gaussian of width 1/R
    w = pow(E, -kR*kR/2.0);
  }
  else {
    fprintf(stderr, "No such filter: %i\nOutput is bogus.\n", FILTER);
    w=0;
  }

  return k*k*p*w*w;
}
double sigma_z0(double M){
  double result, error, lower_limit, upper_limit;
  gsl_function F;
  double rel_tol  = FRACT_FLOAT_ERR*10; //<- relative tolerance
  gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);
  double kstart, kend;

  R = MtoR(M);

  // now lets do the integral for sigma and scale it with sigma_norm
  kstart = 1.0e-99/R;
  kend = 350.0/R;
  lower_limit = kstart;//log(kstart);
  upper_limit = kend;//log(kend);

  F.function = &dsigma_dk;
  gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
		       1000, GSL_INTEG_GAUSS61, w, &result, &error);
  gsl_integration_workspace_free (w);

  return sigma_norm * sqrt(result);
}


/*
  Returns the value of the linear power spectrum DENSITY (i.e. <|delta_k|^2>/V)
  at a given k mode linearly extrapolated to z=0
*/
double power_in_k(double k){
  double p, T, gamma, q, aa, bb, cc;

  // get the power spectrum.. choice of 5:
  if (POWER_SPECTRUM == 0){ // Eisenstein & Hu
    T = TFmdm(k);
    // check if we should cuttoff power spectrum according to Bode et al. 2000 transfer function
    if (P_CUTOFF) T *= pow(1 + pow(BODE_e*k*R_CUTOFF, 2*BODE_v), -BODE_n/BODE_v);
    p = pow(k, POWER_INDEX) * T * T;
    //p = pow(k, POWER_INDEX - 0.05*log(k/0.05)) * T * T; //running, alpha=0.05
  }
  else if (POWER_SPECTRUM == 1){ // BBKS
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    q = k / (hlittle*gamma);
    T = (log(1.0+2.34*q)/(2.34*q)) *
      pow( 1.0+3.89*q + pow(16.1*q, 2) + pow( 5.46*q, 3) + pow(6.71*q, 4), -0.25);
     p = pow(k, POWER_INDEX) * T * T;
  }
  else if (POWER_SPECTRUM == 2){ // Efstathiou,G., Bond,J.R., and White,S.D.M., MNRAS,258,1P (1992)
    gamma = 0.25;
    aa = 6.4/(hlittle*gamma);
    bb = 3.0/(hlittle*gamma);
    cc = 1.7/(hlittle*gamma);
    p = pow(k, POWER_INDEX) / pow( 1+pow( aa*k + pow(bb*k, 1.5) + pow(cc*k,2), 1.13), 2.0/1.13 );
  }
  else if (POWER_SPECTRUM == 3){ // Peebles, pg. 626
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    aa = 8.0 / (hlittle*gamma);
    bb = 4.7 / pow(hlittle*gamma, 2);
    p = pow(k, POWER_INDEX) / pow(1 + aa*k + bb*k*k, 2);
  }
  else if (POWER_SPECTRUM == 4){ // White, SDM and Frenk, CS, 1991, 379, 52
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    aa = 1.7/(hlittle*gamma);
    bb = 9.0/pow(hlittle*gamma, 1.5);
    cc = 1.0/pow(hlittle*gamma, 2);
    p = pow(k, POWER_INDEX) * 19400.0 / pow(1 + aa*k + bb*pow(k, 1.5) + cc*k*k, 2);
  }
  else{
    fprintf(stderr, "No such power spectrum defined: %i\nOutput is bogus.\n", POWER_SPECTRUM);
    p = 0;
  }


  return p*TWOPI*PI*sigma_norm*sigma_norm;
}


/*
  FUNCTION dsigmasqdm_z0(M)
  returns  d/dm (sigma^2) (see function sigma), in units of Msun^-1
*/
double dsigmasq_dm(double k, void *params){
  double p, w, T, gamma, q, aa, bb, cc, dwdr, drdm, kR;

  // get the power spectrum.. choice of 5:
  if (POWER_SPECTRUM == 0){ // Eisenstein & Hu ApJ, 1999, 511, 5
    T = TFmdm(k);
    // check if we should cuttoff power spectrum according to Bode et al. 2000 transfer function
    if (P_CUTOFF) T *= pow(1 + pow(BODE_e*k*R_CUTOFF, 2*BODE_v), -BODE_n/BODE_v);
    p = pow(k, POWER_INDEX) * T * T;
    //p = pow(k, POWER_INDEX - 0.05*log(k/0.05)) * T * T; //running, alpha=0.05
  }
  else if (POWER_SPECTRUM == 1){ // BBKS
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    q = k / (hlittle*gamma);
    T = (log(1.0+2.34*q)/(2.34*q)) *
      pow( 1.0+3.89*q + pow(16.1*q, 2) + pow( 5.46*q, 3) + pow(6.71*q, 4), -0.25);
    p = pow(k, POWER_INDEX) * T * T;
  }
  else if (POWER_SPECTRUM == 2){ // Efstathiou,G., Bond,J.R., and White,S.D.M., MNRAS,258,1P (1992)
    gamma = 0.25;
    aa = 6.4/(hlittle*gamma);
    bb = 3.0/(hlittle*gamma);
    cc = 1.7/(hlittle*gamma);
    p = pow(k, POWER_INDEX) / pow( 1+pow( aa*k + pow(bb*k, 1.5) + pow(cc*k,2), 1.13), 2.0/1.13 );
  }
  else if (POWER_SPECTRUM == 3){ // Peebles, pg. 626
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    aa = 8.0 / (hlittle*gamma);
    bb = 4.7 / (hlittle*gamma);
    p = pow(k, POWER_INDEX) / pow(1 + aa*k + bb*k*k, 2);
  }
  else if (POWER_SPECTRUM == 4){ // White, SDM and Frenk, CS, 1991, 379, 52
    gamma = OMm * hlittle * pow(E, -OMb - OMb/OMm);
    aa = 1.7/(hlittle*gamma);
    bb = 9.0/pow(hlittle*gamma, 1.5);
    cc = 1.0/pow(hlittle*gamma, 2);
    p = pow(k, POWER_INDEX) * 19400.0 / pow(1 + aa*k + pow(bb*k, 1.5) + cc*k*k, 2);
  }
  else{
    fprintf(stderr, "No such power spectrum defined: %i\nOutput is bogus.\n", POWER_SPECTRUM);
    p = 0;
  }


  // now get the value of the window function
  kR = k * R;
  if (FILTER == 0){ // top hat
    if ( (kR) < 1.0e-4 ){ w = 1.0; }// w converges to 1 as (kR) -> 0
    else { w = 3.0 * (sin(kR)/pow(kR, 3) - cos(kR)/pow(kR, 2));}

    // now do d(w^2)/dm = 2 w dw/dr dr/dm
    if ( (kR) < 1.0e-10 ){  dwdr = 0;}
    else{ dwdr = 9*cos(kR)*k/pow(kR,3) + 3*sin(kR)*(1 - 3/(kR*kR))/(kR*R);}
	    //3*k*( 3*cos(kR)/pow(kR,3) + sin(kR)*(-3*pow(kR, -4) + 1/(kR*kR)) );}
      //     dwdr = -1e8 * k / (R*1e3);
    drdm = 1.0 / (4.0*PI * OMm*RHOcrit * R*R);
  }
  else if (FILTER == 1){ // gaussian of width 1/R
    w = pow(E, -kR*kR/2.0);
    dwdr = - k*kR * w;
    drdm = 1.0 / (pow(2*PI, 1.5) * OMm*RHOcrit * 3*R*R);
  }
  else {
    fprintf(stderr, "No such filter: %i\nOutput is bogus.\n", FILTER);
    w=0;
  }

  //  printf("%e\t%e\t%e\t%e\t%e\t%e\t%e\n", k, R, p, w, dwdr, drdm, dsigmadk[1]);
  return k*k*p*2*w*dwdr*drdm * d2fact;
}
double dsigmasqdm_z0(double M){
  double result, error, lower_limit, upper_limit;
  gsl_function F;
  double rel_tol  = FRACT_FLOAT_ERR*10; //<- relative tolerance
  gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);
  double kstart, kend;

  R = MtoR(M);

  // now lets do the integral for sigma and scale it with sigma_norm
  kstart = 1.0e-99/R;
  kend = 350.0/R;
  lower_limit = kstart;//log(kstart);
  upper_limit = kend;//log(kend);
  d2fact = M*10000/sigma_z0(M);

  F.function = &dsigmasq_dm;
  gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
		       1000, GSL_INTEG_GAUSS61, w, &result, &error);
  gsl_integration_workspace_free (w);

  return sigma_norm * sigma_norm * result /d2fact;
}



/*
  FUNCTION TFmdm is the power spectrum transfer function from Eisenstein & Hu ApJ, 1999, 511, 5
*/
double TFmdm(double k){
  double q, gamma_eff, q_eff, TF_m, q_nu;

  q = k*pow(theta_cmb,2)/omhh;
  gamma_eff=sqrt(alpha_nu) + (1.0-sqrt(alpha_nu))/(1.0+pow(0.43*k*sound_horizon, 4));
  q_eff = q/gamma_eff;
  TF_m= log(E+1.84*beta_c*sqrt(alpha_nu)*q_eff);
  TF_m /= TF_m + pow(q_eff,2) * (14.4 + 325.0/(1.0+60.5*pow(q_eff,1.11)));
  q_nu = 3.92*q/sqrt(f_nu/N_nu);
  TF_m *= 1.0 + (1.2*pow(f_nu,0.64)*pow(N_nu,0.3+0.6*f_nu)) /
    (pow(q_nu,-1.6)+pow(q_nu,0.8));

  //   printf("%f  %e  %f  %f  %f  %f\n",omhh,f_nu,f_baryon,N_nu,y_d,alpha_nu);
  // printf("%f  %f  %f  %f\n", beta_c,sound_horizon,theta_cmb,z_equality);
  //printf("%f  %e  %f  %f  %f\n\n",q, k, gamma_eff, q_nu, TF_m);


  return TF_m;
}


void TFset_parameters(){
  double z_drag, R_drag, R_equality, p_c, p_cb, f_c, f_cb, f_nub, k_equality;

  z_equality = 25000*omhh*pow(theta_cmb, -4) - 1.0;
  k_equality = 0.0746*omhh/(theta_cmb*theta_cmb);

  z_drag = 0.313*pow(omhh,-0.419) * (1 + 0.607*pow(omhh, 0.674));
  z_drag = 1 + z_drag*pow(OMb*hlittle*hlittle, 0.238*pow(omhh, 0.223));
  z_drag *= 1291 * pow(omhh, 0.251) / (1 + 0.659*pow(omhh, 0.828));

  y_d = (1 + z_equality) / (1.0 + z_drag);

  R_drag = 31.5 * OMb*hlittle*hlittle * pow(theta_cmb, -4) * 1000 / (1.0 + z_drag);
  R_equality = 31.5 * OMb*hlittle*hlittle * pow(theta_cmb, -4) * 1000 / (1.0 + z_equality);

  sound_horizon = 2.0/3.0/k_equality * sqrt(6.0/R_equality) *
    log( (sqrt(1+R_drag) + sqrt(R_drag+R_equality)) / (1.0 + sqrt(R_equality)) );

  p_c = -(5 - sqrt(1 + 24*(1 - f_nu-f_baryon)))/4.0;
  p_cb = -(5 - sqrt(1 + 24*(1 - f_nu)))/4.0;
  f_c = 1 - f_nu - f_baryon;
  f_cb = 1 - f_nu;
  f_nub = f_nu+f_baryon;

  alpha_nu = (f_c/f_cb) * (2*(p_c+p_cb)+5)/(4*p_cb+5.0);
  alpha_nu *= 1 - 0.553*f_nub+0.126*pow(f_nub,3);
  alpha_nu /= 1-0.193*sqrt(f_nu)+0.169*f_nu;
  alpha_nu *= pow(1+y_d, p_c-p_cb);
  alpha_nu *= 1+ (p_cb-p_c)/2.0 * (1.0+1.0/(4.0*p_c+3.0)/(4.0*p_cb+7.0))/(1.0+y_d);
  beta_c = 1.0/(1.0-0.949*f_nub);
}




double init_ps(){
  double result, error, lower_limit, upper_limit;
  gsl_function F;
  double rel_tol  = FRACT_FLOAT_ERR*10; //<- relative tolerance
  gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);
  double kstart, kend;
  int i;
  double x;

  // Set cuttoff scale for WDM (eq. 4 in Barkana et al. 2001) in comoving Mpc
  R_CUTOFF = 0.201*pow((OMm-OMb)*hlittle*hlittle/0.15, 0.15)*pow(g_x/1.5, -0.29)*pow(M_WDM, -1.15);


  if (P_CUTOFF){
    fprintf(stderr, "For M_DM = %.2e keV, R_CUTOFF is: %.2e comoving Mpc\n", M_WDM, R_CUTOFF);
  }

  omhh = OMm*hlittle*hlittle;
  theta_cmb = T_cmb / 2.7;

  // Translate Parameters into forms GLOBALVARIABLES form
  f_nu = OMn/OMm;
  f_baryon = OMb/OMm;
  if (f_nu < TINY) f_nu = 1e-10;
  if (f_baryon < TINY) f_baryon = 1e-10;


  TFset_parameters();

  sigma_norm = -1;

  R = 8.0/hlittle;
  kstart = 1.0e-99/R;
  kend = 350.0/R;
  lower_limit = kstart;//log(kstart);
  upper_limit = kend;//log(kend);

  F.function = &dsigma_dk;

  gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
		       1000, GSL_INTEG_GAUSS61, w, &result, &error);
  gsl_integration_workspace_free (w);

  sigma_norm = SIGMA8/sqrt(result); //takes care of volume factor


  /* initialize the lookup table for erfc */
  /*
  for (i=0; i<=ERFC_NPTS; i++){
    erfc_params[i] = i*ERFC_PARAM_DELTA;
    log_erfc_table[i] = log(erfcc(erfc_params[i]));
  }
  // Set up spline table
  erfc_acc   = gsl_interp_accel_alloc ();
  erfc_spline  = gsl_spline_alloc (gsl_interp_cspline, ERFC_NPTS);
  gsl_spline_init(erfc_spline, erfc_params, log_erfc_table, ERFC_NPTS);
  */

  return R_CUTOFF;
}

void free_ps(){
  /*    gsl_spline_free (erfc_spline);
    gsl_interp_accel_free(erfc_acc);
  */
  return;
}

double splined_erfc(double x){
  if (x < 0){
    //    fprintf(stderr, "WARNING: Negative value %e passed to splined_erfc. Returning 1\n", x);
    return 1;
  }
  return erfcc(x); // the interpolation below doesn't seem to be stable in Ts.c
  if (x > ERFC_PARAM_DELTA*(ERFC_NPTS-1))
    return erfcc(x);
  else
    return exp(gsl_spline_eval(erfc_spline, x, erfc_acc));
}

float FgtrConditionalM_second(float z, float M1, float M2, float MFeedback, float alpha, float delta1, float delta2) {

    return exp(M1)*pow(exp(M1)/MFeedback,alpha)*dNdM_conditional_second(z,M1,M2,delta1,delta2)/sqrt(2.*PI);
}

float dNdM_conditional_second(float z, float M1, float M2, float delta1, float delta2){

    float sigma1, sigma2, dsigmadm, dicke_growth,dsigma_val;

    M1 = exp(M1);
    M2 = exp(M2);

    dicke_growth = dicke(z);

    splint(Mass_Spline-1,Sigma_Spline-1,second_derivs_sigma-1,(int)NMass,M1,&(sigma1));
    splint(Mass_Spline-1,Sigma_Spline-1,second_derivs_sigma-1,(int)NMass,M2,&(sigma2));

    sigma1 = sigma1*sigma1;
    sigma2 = sigma2*sigma2;

	//printf("sigma_min^2 = %.4e sigma^2 = %.4e\n",sigma1,sigma2);

    splint(Mass_Spline-1,dSigmadm_Spline-1,second_derivs_dsigma-1,(int)NMass,M1,&(dsigma_val));

    dsigmadm = -pow(10.,dsigma_val)/(2.0*sigma1); // This is actually sigma1^{2} as calculated above, however, it should just be sigma1. It cancels with the same factor below. Why I have decided to write it like that I don't know!
    //if(z > 33 && z < 35)
	//printf("dNdM_conditional: z = %.4f, sigma_min^2 = %.4e sigma^2 = %.4e, delta1 = %.4e, delta2 = %.4e \n",z,sigma1,sigma2,delta1,delta2);
    if((sigma1 > sigma2)) {

        return -(( delta1 - delta2 )/dicke_growth)*( 2.*sigma1*dsigmadm )*( exp( - ( delta1 - delta2 )*( delta1 - delta2 )/( 2.*dicke_growth*dicke_growth*( sigma1 - sigma2 ) ) ) )/(pow( sigma1 - sigma2, 1.5));
    }
    else if(sigma1==sigma2) {

        return -(( delta1 - delta2 )/dicke_growth)*( 2.*sigma1*dsigmadm )*( exp( - ( delta1 - delta2 )*( delta1 - delta2 )/( 2.*dicke_growth*dicke_growth*( 1.e-6 ) ) ) )/(pow( 1.e-6, 1.5));

    }
    else {
        return 0.;
    }
}

void gauleg(float x1, float x2, float x[], float w[], int n)
//Given the lower and upper limits of integration x1 and x2, and given n, this routine returns arrays x[1..n] and w[1..n] of length n,
//containing the abscissas and weights of the Gauss- Legendre n-point quadrature formula.
{

    int m,j,i;
    double z1,z,xm,xl,pp,p3,p2,p1;

    m=(n+1)/2;
    xm=0.5*(x2+x1);
    xl=0.5*(x2-x1);
    for (i=1;i<=m;i++) {
        //High precision is a good idea for this routine.
        //The roots are symmetric in the interval, so we only have to find half of them.
        //Loop over the desired roots.

        z=cos(3.141592654*(i-0.25)/(n+0.5));

        //Starting with the above approximation to the ith root, we enter the main loop of refinement by Newton’s method.
        do {
            p1=1.0;
            p2=0.0;
            for (j=1;j<=n;j++) {
                //Loop up the recurrence relation to get the Legendre polynomial evaluated at z.
                p3=p2;
                p2=p1;
                p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
            }
            //p1 is now the desired Legendre polynomial. We next compute pp, its derivative, by a standard relation involving also p2,
            //the polynomial of one lower order.
            pp=n*(z*p1-p2)/(z*z-1.0);
            z1=z;
            z=z1-p1/pp;
        } while (fabs(z-z1) > EPS2);
        x[i]=xm-xl*z;
        x[n+1-i]=xm+xl*z;
        w[i]=2.0*xl/((1.0-z*z)*pp*pp);
        w[n+1-i]=w[i];
    }
}

void nrerror(char error_text[])
{
    fprintf(stderr,"Numerical Recipes run-time error...\n");
    fprintf(stderr,"%s\n",error_text);
    fprintf(stderr,"...now exiting to system...\n");
    exit(1);
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
    float *v;
    v = (float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
    if(!v) nrerror("allocation failure in vector()");
    return v - nl + NR_END;
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
    free((FREE_ARG) (v+nl-NR_END));
}

void spline(float x[], float y[], int n, float yp1, float ypn, float y2[])
/*Given arrays x[1..n] and y[1..n] containing a tabulated function, i.e., yi = f(xi), with
 x1 <x2 < :: : < xN, and given values yp1 and ypn for the first derivative of the interpolating
 function at points 1 and n, respectively, this routine returns an array y2[1..n] that contains
 the second derivatives of the interpolating function at the tabulated points xi. If yp1 and/or
 ypn are equal to 1e30 or larger, the routine is signaled to set the corresponding boundary
 condition for a natural spline, with zero second derivative on that boundary.*/
{
    int i,k;
    float p,qn,sig,un,*u;
    int na,nb,check;
    u=vector(1,n-1);
    if (yp1 > 0.99e30)                     // The lower boundary condition is set either to be "natural"
        y2[1]=u[1]=0.0;
    else {                                 // or else to have a specified first derivative.
        y2[1] = -0.5;
        u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
    }
    for (i=2;i<=n-1;i++) {                              //This is the decomposition loop of the tridiagonal algorithm.
        sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);                //y2 and u are used for temporary
        na = 1;
        nb = 1;
        check = 0;
        while(((float)(x[i+na*1]-x[i-nb*1])==(float)0.0)) {
            check = check + 1;
            if(check%2==0) {
                na = na + 1;
            }
            else {
                nb = nb + 1;
            }
            sig=(x[i]-x[i-1])/(x[i+na*1]-x[i-nb*1]);
        }
        p=sig*y2[i-1]+2.0;                                //storage of the decomposed
        y2[i]=(sig-1.0)/p;                                //  factors.
        u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
        u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;

        if(((float)(x[i+1]-x[i])==(float)0.0) || ((float)(x[i]-x[i-1])==(float)0.0)) {
            na = 0;
            nb = 0;
            check = 0;
            while((float)(x[i+na*1]-x[i-nb])==(float)(0.0) || ((float)(x[i+na]-x[i-nb*1])==(float)0.0)) {
                check = check + 1;
                if(check%2==0) {
                    na = na + 1;
                }
                else {
                    nb = nb + 1;
                }
            }
            u[i]=(y[i+1]-y[i])/(x[i+na*1]-x[i-nb]) - (y[i]-y[i-1])/(x[i+na]-x[i-nb*1]);

            u[i]=(6.0*u[i]/(x[i+na*1]-x[i-nb*1])-sig*u[i-1])/p;

        }
    }
    if (ypn > 0.99e30)                        //The upper boundary condition is set either to be "natural"
        qn=un=0.0;
    else {                                    //or else to have a specified first derivative.
        qn=0.5;
        un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
    }
    y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);

    for (k=n-1;k>=1;k--) {                      //This is the backsubstitution loop of the tridiagonal
        y2[k]=y2[k]*y2[k+1]+u[k];               //algorithm.
    }
    free_vector(u,1,n-1);
}


void splint(float xa[], float ya[], float y2a[], int n, float x, float *y)
/*Given the arrays xa[1..n] and ya[1..n], which tabulate a function (with the xai's in order),
 and given the array y2a[1..n], which is the output from spline above, and given a value of
 x, this routine returns a cubic-spline interpolated value y.*/
{
    void nrerror(char error_text[]);
    int klo,khi,k;
    float h,b,a;
    klo=1;                                                  // We will find the right place in the table by means of
    khi=n;                                                  //bisection. This is optimal if sequential calls to this
    while (khi-klo > 1) {                                   //routine are at random values of x. If sequential calls
        k=(khi+klo) >> 1;                                     //are in order, and closely spaced, one would do better
        if (xa[k] > x) khi=k;                                 //to store previous values of klo and khi and test if
        else klo=k;                                           //they remain appropriate on the next call.
    }                                                           // klo and khi now bracket the input value of x.
    h=xa[khi]-xa[klo];
    if (h == 0.0) nrerror("Bad xa input to routine splint");    //The xa's must be distinct.
    a=(xa[khi]-x)/h;
    b=(x-xa[klo])/h;                                            //Cubic spline polynomial is now evaluated.
    *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
    unsigned long *v;
    v = (unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
    if(!v) nrerror("allocation failure in lvector()");
    return v - nl + NR_END;
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
    free((FREE_ARG) (v+nl-NR_END));
}


void initialiseSplinedSigmaM(float M_Min, float M_Max)
{
    int i;
    float Mass;

    Mass_Spline = calloc(NMass,sizeof(float));
    Sigma_Spline = calloc(NMass,sizeof(float));
    dSigmadm_Spline = calloc(NMass,sizeof(float));
    second_derivs_sigma = calloc(NMass,sizeof(float));
    second_derivs_dsigma = calloc(NMass,sizeof(float));

    for(i=0;i<NMass;i++) {
        Mass_Spline[i] = pow(10., log10(M_Min) + (float)i/(NMass-1)*( log10(M_Max) - log10(M_Min) ) );
        Sigma_Spline[i] = sigma_z0(Mass_Spline[i]);
        dSigmadm_Spline[i] = log10(-dsigmasqdm_z0(Mass_Spline[i]));
    }
    spline(Mass_Spline-1,Sigma_Spline-1,NMass,0,0,second_derivs_sigma-1);
    spline(Mass_Spline-1,dSigmadm_Spline-1,NMass,0,0,second_derivs_dsigma-1);
}


float Mass_limit (float logM, float PL, float FRAC) {
	return FRAC*pow(pow(10.,logM)/1e10,PL);
}
void bisection(float *x, float xlow, float xup, int *iter){
    *x=(xlow + xup)/2.;
    ++(*iter);
}

float Mass_limit_bisection(float Mmin, float Mmax, float PL, float FRAC){
  int i, iter, max_iter=200;
  float rel_tol=0.001;
  float logMlow, logMupper, x, x1;
  iter = 0;
  logMlow = log10(Mmin);
  logMupper = log10(Mmax);

  if (PL < 0.) {
    if (Mass_limit(logMlow,PL,FRAC) < 1.) {
      return Mmin;
    }
  }
  else if (PL > 0.) {
    if (Mass_limit(logMupper,PL,FRAC) < 1.) {
      return Mmax;
    }
  }
  else
	return 0;
  bisection(&x, logMlow, logMupper, &iter);
  do {
    if((Mass_limit(logMlow,PL,FRAC)-1.)*(Mass_limit(x,PL,FRAC)-1.) < 0.)
      logMupper = x;
    else
      logMlow = x;
    bisection(&x1, logMlow, logMupper, &iter);
    if(fabs(x1-x) < rel_tol) {
      return pow(10.,x1);
    }
    x = x1;
  }
  while(iter < max_iter);
  printf("\n Failed to find a mass limit to regulate stellar fraction/escape fraction is between 0 and 1.\n");
  printf(" The solution does not converge or iterations are not sufficient\n");
  return -1;
}

/*
 FUNCTION FgtrM_st_SFR(z, Mturn)
 Computes the fraction of mass contained in haloes with mass > M at redshift z
 Uses Sheth-Torman correction
 */
double dFdlnM_st_SFR(double lnM, void *params){
    struct parameters_gsl_SFR_int_ vals = *(struct parameters_gsl_SFR_int_ *)params;

    double M = exp(lnM);
    float z = vals.z_obs;
    double MassTurnover = vals.Mdrop;
    double Alpha_star = vals.pl_star;
    double Alpha_esc = vals.pl_esc;
	double Fstar10 = vals.frac_star;
	double Fesc10 = vals.frac_esc;
	double Mlim_Fstar = vals.LimitMass_Fstar;
	double Mlim_Fesc = vals.LimitMass_Fesc;

	double Fstar, Fesc;


	if (Alpha_star > 0. && M > Mlim_Fstar)
		Fstar = 1./Fstar10;
	else if (Alpha_star < 0. && M < Mlim_Fstar)
		Fstar = 1/Fstar10;
	else
		Fstar = pow(M/1e10,Alpha_star);

	if (Alpha_esc > 0. && M > Mlim_Fesc)
		Fesc = 1./Fesc10;
	else if (Alpha_esc < 0. && M < Mlim_Fesc)
		Fesc = 1./Fesc10;
	else
		Fesc = pow(M/1e10,Alpha_esc);
  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    return dNdM_st(z,M) * M * M * src.fstar(z, M) * src.fesc(z, M) * src.Nion(z, M);
  }
    return dNdM_st(z,M) * M * M * exp(-MassTurnover/M) * Fstar * Fesc;
}

double FgtrM_st_SFR(double z, double MassTurnover, double Alpha_star, double Alpha_esc, double Fstar10, double Fesc10, double Mlim_Fstar, double Mlim_Fesc){
	double M_Min = MassTurnover/50.;
	//double M_Min = 1e8;
    double result, error, lower_limit, upper_limit;
    gsl_function F;
    double rel_tol = 0.001; //<- relative tolerance
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

    struct parameters_gsl_SFR_int_ parameters_gsl_SFR = {
        .z_obs = z,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
        .pl_esc = Alpha_esc,
		.frac_star = Fstar10,
		.frac_esc = Fesc10,
		.LimitMass_Fstar = Mlim_Fstar,
		.LimitMass_Fesc = Mlim_Fesc,
    };

    if(USE_GENERAL_SOURCES)
    {
      sources s;
      s = defaultSources();
      M_Min = s.minMass(z);
    }

    F.function = &dFdlnM_st_SFR;
    F.params = &parameters_gsl_SFR;
    lower_limit = log(M_Min);
    upper_limit = log(FMAX(1e16, M_Min*100));

    gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                        1000, GSL_INTEG_GAUSS61, w, &result, &error);
    gsl_integration_workspace_free (w);

    return result / (OMm*RHOcrit);

}

/*
 FUNCTION FgtrM_st_SFRPOP2(z, Mturn)
 */
double dFdlnM_st_SFRPOP2(double lnM, void *params){
  struct parameters_gsl_SFR_int_ vals = *(struct parameters_gsl_SFR_int_ *)params;
  double M = exp(lnM);
  float z = vals.z_obs;
  double MassTurnover = vals.Mdrop;
  double Alpha_star = vals.pl_star;
	double Fstar10 = vals.frac_star;
	double Mlim_Fstar = vals.LimitMass_Fstar;
	double Fstar;

  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    return dNdM_st(z, M) * M * src.sfr(z, M);
  }
    return 0;
}

double dFdlnM_st_gas(double lnM, void *params){
  struct parameters_gsl_SFR_int_ vals = *(struct parameters_gsl_SFR_int_ *)params;
  double M = exp(lnM);
  float z = vals.z_obs;
  double MassTurnover = vals.Mdrop;
  double Alpha_star = vals.pl_star;
  double Fstar10 = vals.frac_star;
  double Mlim_Fstar = vals.LimitMass_Fstar;
  double Fstar;

  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    return dNdM_st(z, M) * M * src.Mg(z, M);
  }
    return 0;
}
double dFdlnM_st_met(double lnM, void *params){
  struct parameters_gsl_SFR_int_ vals = *(struct parameters_gsl_SFR_int_ *)params;
  double M = exp(lnM);
  float z = vals.z_obs;
  double MassTurnover = vals.Mdrop;
  double Alpha_star = vals.pl_star;
  double Fstar10 = vals.frac_star;
  double Mlim_Fstar = vals.LimitMass_Fstar;
  double Fstar;

  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    return dNdM_st(z, M) * M * src.Mz(z, M);
  }
    return 0;
}

double FgtrM_st_SFRPOP2(double z, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar){
  double M_Min = MassTurnover/50.;
  double result, error, lower_limit, upper_limit;
  gsl_function F;
  double rel_tol = 0.001; //<- relative tolerance
  gsl_integration_workspace * w
  = gsl_integration_workspace_alloc (1000);

  struct parameters_gsl_SFR_int_ parameters_gsl_SFR = {
      .z_obs = z,
      .Mdrop = MassTurnover,
      .pl_star = Alpha_star,
	    .frac_star = Fstar10,
	    .LimitMass_Fstar = Mlim_Fstar,
  };

  if(USE_GENERAL_SOURCES)
  {
    sources s;
    s = defaultSources();
    M_Min = s.minMass(z);
  }

  F.function = &dFdlnM_st_SFRPOP2;
  F.params = &parameters_gsl_SFR;
  lower_limit = log(M_Min);
  upper_limit = log(FMAX(1e15, M_Min*100));

  gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                      1000, GSL_INTEG_GAUSS61, w, &result, &error);
  gsl_integration_workspace_free (w);

  return result / (OMm*RHOcrit);

}

double FgtrM_st_gas(double z, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar){
  double M_Min = MassTurnover/50.;
  double result, error, lower_limit, upper_limit;
  gsl_function F;
  double rel_tol = 0.001; //<- relative tolerance
  gsl_integration_workspace * w
  = gsl_integration_workspace_alloc (1000);

  struct parameters_gsl_SFR_int_ parameters_gsl_SFR = {
      .z_obs = z,
      .Mdrop = MassTurnover,
      .pl_star = Alpha_star,
      .frac_star = Fstar10,
      .LimitMass_Fstar = Mlim_Fstar,
  };

  if(USE_GENERAL_SOURCES)
  {
    sources s;
    s = defaultSources();
    M_Min = s.minMass(z);
  }

  F.function = &dFdlnM_st_gas;
  F.params = &parameters_gsl_SFR;
  lower_limit = log(M_Min);
  upper_limit = log(FMAX(1e15, M_Min*100));

  gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                      1000, GSL_INTEG_GAUSS61, w, &result, &error);
  gsl_integration_workspace_free (w);

  return result / (OMm*RHOcrit);

}

double FgtrM_st_met(double z, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar){
  double M_Min = MassTurnover/50.;
  double result, error, lower_limit, upper_limit;
  gsl_function F;
  double rel_tol = 0.001; //<- relative tolerance
  gsl_integration_workspace * w
  = gsl_integration_workspace_alloc (1000);

  struct parameters_gsl_SFR_int_ parameters_gsl_SFR = {
      .z_obs = z,
      .Mdrop = MassTurnover,
      .pl_star = Alpha_star,
      .frac_star = Fstar10,
      .LimitMass_Fstar = Mlim_Fstar,
  };

  if(USE_GENERAL_SOURCES)
  {
    sources s;
    s = defaultSources();
    M_Min = s.minMass(z);
  }

  F.function = &dFdlnM_st_met;
  F.params = &parameters_gsl_SFR;
  lower_limit = log(M_Min);
  upper_limit = log(FMAX(1e15, M_Min*100));

  gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                      1000, GSL_INTEG_GAUSS61, w, &result, &error);
  gsl_integration_workspace_free (w);

  return result / (OMm*RHOcrit);

}





/*
 FUNCTION FgtrM_st_ZNUMPOP2(z, Mturn)
 */
double dFdlnM_st_ZNUMEPOP2(double lnM, void *params){
  struct parameters_gsl_SFR_int_ vals = *(struct parameters_gsl_SFR_int_ *)params;
  double M = exp(lnM);
  float z = vals.z_obs;
  double MassTurnover = vals.Mdrop;
  double Alpha_star = vals.pl_star;
  double Fstar10 = vals.frac_star;
  double Fstar;
  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    return dNdM_st(z, M) * M * src.metal(z, M) * src.Mg(z, M);
  }
    return 0;
}
double FgtrM_st_ZNUMEPOP2(double z, double MassTurnover, double Alpha_star, double Fstar10){
  double M_Min = MassTurnover/50.;
  double result, error, lower_limit, upper_limit;
  gsl_function F;
  double rel_tol = 0.001; //<- relative tolerance
  gsl_integration_workspace * w
  = gsl_integration_workspace_alloc (1000);
  struct parameters_gsl_SFR_int_ parameters_gsl_SFR = {
      .z_obs = z,
      .Mdrop = MassTurnover,
      .pl_star = Alpha_star,
      .frac_star = Fstar10,
  };
  if(USE_GENERAL_SOURCES)
  {
    sources s;
    s = defaultSources();
    M_Min = s.minMass(z);
  }
  F.function = &dFdlnM_st_ZNUMEPOP2;
  F.params = &parameters_gsl_SFR;
  lower_limit = log(M_Min);
  upper_limit = log(FMAX(1e15, M_Min*100));
  gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                      1000, GSL_INTEG_GAUSS61, w, &result, &error);
  gsl_integration_workspace_free (w);
  return result / (OMm*RHOcrit);
}
/*
 FUNCTION FgtrM_st_ZDENOPOP2(z, Mturn)
 */
double dFdlnM_st_ZDENOPOP2(double lnM, void *params){
  struct parameters_gsl_SFR_int_ vals = *(struct parameters_gsl_SFR_int_ *)params;
  double M = exp(lnM);
  float z = vals.z_obs;
  double MassTurnover = vals.Mdrop;
  double Alpha_star = vals.pl_star;
  double Fstar10 = vals.frac_star;
  double Fstar;
  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    return dNdM_st(z, M) * M * src.Mg(z, M);
  }
    return 0;
}
double FgtrM_st_ZDENOPOP2(double z, double MassTurnover, double Alpha_star, double Fstar10){
  double M_Min = MassTurnover/50.;
  double result, error, lower_limit, upper_limit;
  gsl_function F;
  double rel_tol = 0.001; //<- relative tolerance
  gsl_integration_workspace * w
  = gsl_integration_workspace_alloc (1000);
  struct parameters_gsl_SFR_int_ parameters_gsl_SFR = {
      .z_obs = z,
      .Mdrop = MassTurnover,
      .pl_star = Alpha_star,
      .frac_star = Fstar10,
  };
  if(USE_GENERAL_SOURCES)
  {
    sources s;
    s = defaultSources();
    M_Min = s.minMass(z);
  }
  F.function = &dFdlnM_st_ZDENOPOP2;
  F.params = &parameters_gsl_SFR;
  lower_limit = log(M_Min);
  upper_limit = log(FMAX(1e15, M_Min*100));
  gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                      1000, GSL_INTEG_GAUSS61, w, &result, &error);
  gsl_integration_workspace_free (w);
  return result / (OMm*RHOcrit);
}






/*
 FUNCTION FgtrM_st_LINELUMPOP2A(z, Mturn)
 */
double dFdlnM_st_LINELUMPOP2A(double lnM, void *params){
  struct parameters_gsl_SFR_int_ vals = *(struct parameters_gsl_SFR_int_ *)params;
  double M = exp(lnM);
  float z = vals.z_obs;
  double MassTurnover = vals.Mdrop;
  double Alpha_star = vals.pl_star;
	double Fstar10 = vals.frac_star;
	double Mlim_Fstar = vals.LimitMass_Fstar;
  double LineID = vals.lid;
	double Fstar;

  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    return dNdM_st(z, M) * M * src.linelum(z, M, -2., LineID);
  }
    return 0;
}

double FgtrM_st_LINELUMPOP2A(double z, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar, int LineID){
  double M_Min = MassTurnover/50.;
  double result, error, lower_limit, upper_limit;
  gsl_function F;
  double rel_tol = 0.001; //<- relative tolerance
  gsl_integration_workspace * w
  = gsl_integration_workspace_alloc (1000);

  struct parameters_gsl_SFR_int_ parameters_gsl_SFR = {
      .z_obs = z,
      .Mdrop = MassTurnover,
      .pl_star = Alpha_star,
	    .frac_star = Fstar10,
	    .LimitMass_Fstar = Mlim_Fstar,
      .lid = LineID,
  };

  if(USE_GENERAL_SOURCES)
  {
    sources s;
    s = defaultSources();
    M_Min = s.minMass(z);
  }

  F.function = &dFdlnM_st_LINELUMPOP2A;
  F.params = &parameters_gsl_SFR;
  lower_limit = log(M_Min);
  upper_limit = log(FMAX(1e15, M_Min*100));

  gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                      1000, GSL_INTEG_GAUSS61, w, &result, &error);
  gsl_integration_workspace_free (w);

  return result / (OMm*RHOcrit);

}

/*
 FUNCTION FgtrM_st_LINELUMPOP2B(z, Mturn)
 */
double dFdlnM_st_LINELUMPOP2B(double lnM, void *params){
  struct parameters_gsl_SFR_int_ vals = *(struct parameters_gsl_SFR_int_ *)params;
  double M = exp(lnM);
  float z = vals.z_obs;
  double MassTurnover = vals.Mdrop;
  double Alpha_star = vals.pl_star;
	double Fstar10 = vals.frac_star;
	double Mlim_Fstar = vals.LimitMass_Fstar;
  double LineID = vals.lid;
	double Fstar;

  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    return dNdM_st(z, M) * M * src.linelum(z, M, -2., LineID);
  }
    return 0;
}

double FgtrM_st_LINELUMPOP2B(double z, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar, int LineID){
  double M_Min = MassTurnover/50.;
  double result, error, lower_limit, upper_limit;
  gsl_function F;
  double rel_tol = 0.001; //<- relative tolerance
  gsl_integration_workspace * w
  = gsl_integration_workspace_alloc (1000);

  struct parameters_gsl_SFR_int_ parameters_gsl_SFR = {
      .z_obs = z,
      .Mdrop = MassTurnover,
      .pl_star = Alpha_star,
	    .frac_star = Fstar10,
	    .LimitMass_Fstar = Mlim_Fstar,
      .lid = LineID,
  };

  if(USE_GENERAL_SOURCES)
  {
    sources s;
    s = defaultSources();
    M_Min = s.minMass(z);
  }

  F.function = &dFdlnM_st_LINELUMPOP2B;
  F.params = &parameters_gsl_SFR;
  lower_limit = log(M_Min);
  upper_limit = log(FMAX(1e15, M_Min*100));

  gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                      1000, GSL_INTEG_GAUSS61, w, &result, &error);
  gsl_integration_workspace_free (w);

  return result / (OMm*RHOcrit);

}

/*
 FUNCTION FgtrM_st_LINELUMPOP2C(z, Mturn)
 */
double dFdlnM_st_LINELUMPOP2C(double lnM, void *params){
  struct parameters_gsl_SFR_int_ vals = *(struct parameters_gsl_SFR_int_ *)params;
  double M = exp(lnM);
  float z = vals.z_obs;
  double MassTurnover = vals.Mdrop;
  double Alpha_star = vals.pl_star;
	double Fstar10 = vals.frac_star;
	double Mlim_Fstar = vals.LimitMass_Fstar;
  double LineID = vals.lid;
	double Fstar;

  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    return dNdM_st(z, M) * M * src.linelum(z, M, -2., LineID);
  }
    return 0;
}

double FgtrM_st_LINELUMPOP2C(double z, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar, int LineID){
  double M_Min = MassTurnover/50.;
  double result, error, lower_limit, upper_limit;
  gsl_function F;
  double rel_tol = 0.001; //<- relative tolerance
  gsl_integration_workspace * w
  = gsl_integration_workspace_alloc (1000);

  struct parameters_gsl_SFR_int_ parameters_gsl_SFR = {
      .z_obs = z,
      .Mdrop = MassTurnover,
      .pl_star = Alpha_star,
	    .frac_star = Fstar10,
	    .LimitMass_Fstar = Mlim_Fstar,
      .lid = LineID,
  };

  if(USE_GENERAL_SOURCES)
  {
    sources s;
    s = defaultSources();
    M_Min = s.minMass(z);
  }

  F.function = &dFdlnM_st_LINELUMPOP2C;
  F.params = &parameters_gsl_SFR;
  lower_limit = log(M_Min);
  upper_limit = log(FMAX(1e15, M_Min*100));

  gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                      1000, GSL_INTEG_GAUSS61, w, &result, &error);
  gsl_integration_workspace_free (w);

  return result / (OMm*RHOcrit);

}

/*
 FUNCTION FgtrM_st_LINELUMPOP2D(z, Mturn)
 */
double dFdlnM_st_LINELUMPOP2D(double lnM, void *params){
  struct parameters_gsl_SFR_int_ vals = *(struct parameters_gsl_SFR_int_ *)params;
  double M = exp(lnM);
  float z = vals.z_obs;
  double MassTurnover = vals.Mdrop;
  double Alpha_star = vals.pl_star;
	double Fstar10 = vals.frac_star;
	double Mlim_Fstar = vals.LimitMass_Fstar;
  double LineID = vals.lid;
	double Fstar;

  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    return dNdM_st(z, M) * M * src.linelum(z, M, -2., LineID);
  }
    return 0;
}

double FgtrM_st_LINELUMPOP2D(double z, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar, int LineID){
  double M_Min = MassTurnover/50.;
  double result, error, lower_limit, upper_limit;
  gsl_function F;
  double rel_tol = 0.001; //<- relative tolerance
  gsl_integration_workspace * w
  = gsl_integration_workspace_alloc (1000);

  struct parameters_gsl_SFR_int_ parameters_gsl_SFR = {
      .z_obs = z,
      .Mdrop = MassTurnover,
      .pl_star = Alpha_star,
	    .frac_star = Fstar10,
	    .LimitMass_Fstar = Mlim_Fstar,
      .lid = LineID,
  };

  if(USE_GENERAL_SOURCES)
  {
    sources s;
    s = defaultSources();
    M_Min = s.minMass(z);
  }

  F.function = &dFdlnM_st_LINELUMPOP2D;
  F.params = &parameters_gsl_SFR;
  lower_limit = log(M_Min);
  upper_limit = log(FMAX(1e15, M_Min*100));

  gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                      1000, GSL_INTEG_GAUSS61, w, &result, &error);
  gsl_integration_workspace_free (w);

  return result / (OMm*RHOcrit);

}

/*
 FUNCTION FgtrM_st_LINELUMPOP2E(z, Mturn)
 */
double dFdlnM_st_LINELUMPOP2E(double lnM, void *params){
  struct parameters_gsl_SFR_int_ vals = *(struct parameters_gsl_SFR_int_ *)params;
  double M = exp(lnM);
  float z = vals.z_obs;
  double MassTurnover = vals.Mdrop;
  double Alpha_star = vals.pl_star;
	double Fstar10 = vals.frac_star;
	double Mlim_Fstar = vals.LimitMass_Fstar;
  double LineID = vals.lid;
	double Fstar;

  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    return dNdM_st(z, M) * M * src.linelum(z, M, -2., LineID);
  }
    return 0;
}

double FgtrM_st_LINELUMPOP2E(double z, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar, int LineID){
  double M_Min = MassTurnover/50.;
  double result, error, lower_limit, upper_limit;
  gsl_function F;
  double rel_tol = 0.001; //<- relative tolerance
  gsl_integration_workspace * w
  = gsl_integration_workspace_alloc (1000);

  struct parameters_gsl_SFR_int_ parameters_gsl_SFR = {
      .z_obs = z,
      .Mdrop = MassTurnover,
      .pl_star = Alpha_star,
	    .frac_star = Fstar10,
	    .LimitMass_Fstar = Mlim_Fstar,
      .lid = LineID,
  };

  if(USE_GENERAL_SOURCES)
  {
    sources s;
    s = defaultSources();
    M_Min = s.minMass(z);
  }

  F.function = &dFdlnM_st_LINELUMPOP2E;
  F.params = &parameters_gsl_SFR;
  lower_limit = log(M_Min);
  upper_limit = log(FMAX(1e15, M_Min*100));

  gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                      1000, GSL_INTEG_GAUSS61, w, &result, &error);
  gsl_integration_workspace_free (w);

  return result / (OMm*RHOcrit);

}

/*
 FUNCTION FgtrM_st_LINELUMPOP2F(z, Mturn)
 */
double dFdlnM_st_LINELUMPOP2F(double lnM, void *params){
  struct parameters_gsl_SFR_int_ vals = *(struct parameters_gsl_SFR_int_ *)params;
  double M = exp(lnM);
  float z = vals.z_obs;
  double MassTurnover = vals.Mdrop;
  double Alpha_star = vals.pl_star;
	double Fstar10 = vals.frac_star;
	double Mlim_Fstar = vals.LimitMass_Fstar;
  double LineID = vals.lid;
	double Fstar;

  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    return dNdM_st(z, M) * M * src.linelum(z, M, -2., LineID);
  }
    return 0;
}

double FgtrM_st_LINELUMPOP2F(double z, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar, int LineID){
  double M_Min = MassTurnover/50.;
  double result, error, lower_limit, upper_limit;
  gsl_function F;
  double rel_tol = 0.001; //<- relative tolerance
  gsl_integration_workspace * w
  = gsl_integration_workspace_alloc (1000);

  struct parameters_gsl_SFR_int_ parameters_gsl_SFR = {
      .z_obs = z,
      .Mdrop = MassTurnover,
      .pl_star = Alpha_star,
	    .frac_star = Fstar10,
	    .LimitMass_Fstar = Mlim_Fstar,
      .lid = LineID,
  };

  if(USE_GENERAL_SOURCES)
  {
    sources s;
    s = defaultSources();
    M_Min = s.minMass(z);
  }

  F.function = &dFdlnM_st_LINELUMPOP2F;
  F.params = &parameters_gsl_SFR;
  lower_limit = log(M_Min);
  upper_limit = log(FMAX(1e15, M_Min*100));

  gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                      1000, GSL_INTEG_GAUSS61, w, &result, &error);
  gsl_integration_workspace_free (w);

  return result / (OMm*RHOcrit);

}

/*
 FUNCTION FgtrM_st_LINELUMPOP2G(z, Mturn)
 */
double dFdlnM_st_LINELUMPOP2G(double lnM, void *params){
  struct parameters_gsl_SFR_int_ vals = *(struct parameters_gsl_SFR_int_ *)params;
  double M = exp(lnM);
  float z = vals.z_obs;
  double MassTurnover = vals.Mdrop;
  double Alpha_star = vals.pl_star;
	double Fstar10 = vals.frac_star;
	double Mlim_Fstar = vals.LimitMass_Fstar;
  double LineID = vals.lid;
	double Fstar;

  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    return dNdM_st(z, M) * M * src.linelum(z, M, -2., LineID);
  }
    return 0;
}

double FgtrM_st_LINELUMPOP2G(double z, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar, int LineID){
  double M_Min = MassTurnover/50.;
  double result, error, lower_limit, upper_limit;
  gsl_function F;
  double rel_tol = 0.001; //<- relative tolerance
  gsl_integration_workspace * w
  = gsl_integration_workspace_alloc (1000);

  struct parameters_gsl_SFR_int_ parameters_gsl_SFR = {
      .z_obs = z,
      .Mdrop = MassTurnover,
      .pl_star = Alpha_star,
	    .frac_star = Fstar10,
	    .LimitMass_Fstar = Mlim_Fstar,
      .lid = LineID,
  };

  if(USE_GENERAL_SOURCES)
  {
    sources s;
    s = defaultSources();
    M_Min = s.minMass(z);
  }

  F.function = &dFdlnM_st_LINELUMPOP2G;
  F.params = &parameters_gsl_SFR;
  lower_limit = log(M_Min);
  upper_limit = log(FMAX(1e15, M_Min*100));

  gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                      1000, GSL_INTEG_GAUSS61, w, &result, &error);
  gsl_integration_workspace_free (w);

  return result / (OMm*RHOcrit);

}

// Same as above, but with fx and fesc = 1 and no Nion
double dFdlnM_st_SFR_Xray(double lnM, void *params){
    struct parameters_gsl_SFR_int_ vals = *(struct parameters_gsl_SFR_int_ *)params;

    double M = exp(lnM);
    float z = vals.z_obs;
    double MassTurnover = vals.Mdrop;
    double Alpha_star = vals.pl_star;
    double Alpha_esc = vals.pl_esc;
  double Fstar10 = vals.frac_star;
  double Fesc10 = vals.frac_esc;
  double Mlim_Fstar = vals.LimitMass_Fstar;
  double Mlim_Fesc = vals.LimitMass_Fesc;

  double Fstar, Fesc;


  if (Alpha_star > 0. && M > Mlim_Fstar)
    Fstar = 1./Fstar10;
  else if (Alpha_star < 0. && M < Mlim_Fstar)
    Fstar = 1/Fstar10;
  else
    Fstar = pow(M/1e10,Alpha_star);

  if (Alpha_esc > 0. && M > Mlim_Fesc)
    Fesc = 1./Fesc10;
  else if (Alpha_esc < 0. && M < Mlim_Fesc)
    Fesc = 1./Fesc10;
  else
    Fesc = pow(M/1e10,Alpha_esc);
  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    /*
    int n_ct_fcoll
    for (n_ct_fcoll=NSPEC_MAX; n_ct>=2; n_ct--){
    if (zpp > zmax(zp, n_ct))
    continue;

  nuprime = nu_n(n_ct)*(1+zpp)/(1.0+zp);
  // remove spectral_emissivity calculation, as it is added into one of the fcoll integrals if USE_GENERAL_SOURCES
  if(USE_GENERAL_SOURCES) sum_lyn[R_ct] += frecycle(n_ct);
  else sum_lyn[R_ct] += frecycle(n_ct) * spectral_emissivity(nuprime, 0, Pop, 0);
      }*/
    //return dNdM_st(z,M) * M * M * src.fstar(z, M) * src.fx(z, M);
    return dNdM_st(z,M) * M * (src.sfr(z, M)/f_baryon)/SperYR * src.fx(z, M);
    //return dNdM_st(z,M) * M * M * src.fstar(z, M) * hubble(z) / 0.5 * src.fx(z, M);
  }
    return dNdM_st(z,M) * M * M * exp(-MassTurnover/M) * Fstar * Fesc;
}

double FgtrM_st_SFR_Xray(double z, double MassTurnover, double Alpha_star, double Alpha_esc, double Fstar10, double Fesc10, double Mlim_Fstar, double Mlim_Fesc){
  double M_Min = MassTurnover/50.;
  //double M_Min = 1e8;
    double result, error, lower_limit, upper_limit;
    gsl_function F;
    double rel_tol = 0.001; //<- relative tolerance
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

    struct parameters_gsl_SFR_int_ parameters_gsl_SFR = {
        .z_obs = z,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
        .pl_esc = Alpha_esc,
    .frac_star = Fstar10,
    .frac_esc = Fesc10,
    .LimitMass_Fstar = Mlim_Fstar,
    .LimitMass_Fesc = Mlim_Fesc,
    };

    if(USE_GENERAL_SOURCES)
    {
      sources s;
      s = defaultSources();
      M_Min = s.minMass(z);
    }

    F.function = &dFdlnM_st_SFR_Xray;
    F.params = &parameters_gsl_SFR;
    lower_limit = log(M_Min);
    upper_limit = log(FMAX(1e16, M_Min*100));

    gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                        1000, GSL_INTEG_GAUSS61, w, &result, &error);
    gsl_integration_workspace_free (w);

    return result / (OMm*RHOcrit);

}

/* Normalization factor for Xray */
double FgtrM_st_SFR_Xray_Norm(double z, double MassTurnover, double Alpha_star, double Alpha_esc, double Fstar10, double Fesc10, double Mlim_Fstar, double Mlim_Fesc){
    double result, error, lower_limit, upper_limit;
    gsl_function F;
    double rel_tol = 0.001; //<- relative tolerance
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);
    struct parameters_gsl_SFR_int_ parameters_gsl_SFR = {
    .z_obs = z,
    .Mdrop = MassTurnover,
    .pl_star = Alpha_star,
    .pl_esc = Alpha_esc,
    .frac_star = Fstar10,
    .frac_esc = Fesc10,
    .LimitMass_Fstar = Mlim_Fstar,
    .LimitMass_Fesc = Mlim_Fesc,
    };
    F.function = &dFdlnM_st_SFR_Xray;
    F.params = &parameters_gsl_SFR;
    //lower_limit = log(M_Min);
    //upper_limit = log(FMAX(1e16, M_Min*100));
    gsl_integration_qag (&F, log(1.0e1), log(1.0e16), 0, rel_tol, 1000, GSL_INTEG_GAUSS61, w, &result, &error);
    gsl_integration_workspace_free (w);
    return result / (OMm*RHOcrit);
}

// Same as above, but with fx and fesc = 1 and no Nion and no fx
double dFdlnM_st_SFR_Lya(double lnM, void *params){
    struct parameters_gsl_SFR_int_ vals = *(struct parameters_gsl_SFR_int_ *)params;

    double M = exp(lnM);
    float z = vals.z_obs;
    double MassTurnover = vals.Mdrop;
    double Alpha_star = vals.pl_star;
    double Alpha_esc = vals.pl_esc;
  double Fstar10 = vals.frac_star;
  double Fesc10 = vals.frac_esc;
  double Mlim_Fstar = vals.LimitMass_Fstar;
  double Mlim_Fesc = vals.LimitMass_Fesc;

  double Fstar, Fesc;


  if (Alpha_star > 0. && M > Mlim_Fstar)
    Fstar = 1./Fstar10;
  else if (Alpha_star < 0. && M < Mlim_Fstar)
    Fstar = 1/Fstar10;
  else
    Fstar = pow(M/1e10,Alpha_star);

  if (Alpha_esc > 0. && M > Mlim_Fesc)
    Fesc = 1./Fesc10;
  else if (Alpha_esc < 0. && M < Mlim_Fesc)
    Fesc = 1./Fesc10;
  else
    Fesc = pow(M/1e10,Alpha_esc);
  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    /*
    int n_ct_fcoll
    for (n_ct_fcoll=NSPEC_MAX; n_ct>=2; n_ct--){
    if (zpp > zmax(zp, n_ct))
    continue;

  nuprime = nu_n(n_ct)*(1+zpp)/(1.0+zp);
  // remove spectral_emissivity calculation, as it is added into one of the fcoll integrals if USE_GENERAL_SOURCES
  if(USE_GENERAL_SOURCES) sum_lyn[R_ct] += frecycle(n_ct);
  else sum_lyn[R_ct] += frecycle(n_ct) * spectral_emissivity(nuprime, 0, Pop, 0);
      }*/
    //return dNdM_st(z,M) * M * M * src.fstar(z, M);
    return dNdM_st(z,M) * M * (src.sfr(z, M)/f_baryon)/SperYR;
    //return dNdM_st(z,M) * M * M * src.fstar(z, M) * hubble(z) / 0.5;
  }
    return dNdM_st(z,M) * M * M * exp(-MassTurnover/M) * Fstar * Fesc;
}

double FgtrM_st_SFR_Lya(double z, double MassTurnover, double Alpha_star, double Alpha_esc, double Fstar10, double Fesc10, double Mlim_Fstar, double Mlim_Fesc){
  double M_Min = MassTurnover/50.;
  //double M_Min = 1e8;
    double result, error, lower_limit, upper_limit;
    gsl_function F;
    double rel_tol = 0.001; //<- relative tolerance
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

    struct parameters_gsl_SFR_int_ parameters_gsl_SFR = {
        .z_obs = z,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
        .pl_esc = Alpha_esc,
    .frac_star = Fstar10,
    .frac_esc = Fesc10,
    .LimitMass_Fstar = Mlim_Fstar,
    .LimitMass_Fesc = Mlim_Fesc,
    };

    if(USE_GENERAL_SOURCES)
    {
      sources s;
      s = defaultSources();
      M_Min = s.minMass(z);
    }

    F.function = &dFdlnM_st_SFR_Lya;
    F.params = &parameters_gsl_SFR;
    lower_limit = log(M_Min);
    upper_limit = log(FMAX(1e16, M_Min*100));

    gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                        1000, GSL_INTEG_GAUSS61, w, &result, &error);
    gsl_integration_workspace_free (w);

    return result / (OMm*RHOcrit);

}

/* Normalization factor for Lya */
double FgtrM_st_SFR_Lya_Norm(double z, double MassTurnover, double Alpha_star, double Alpha_esc, double Fstar10, double Fesc10, double Mlim_Fstar, double Mlim_Fesc){
    double result, error, lower_limit, upper_limit;
    gsl_function F;
    double rel_tol = 0.001; //<- relative tolerance
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);
    struct parameters_gsl_SFR_int_ parameters_gsl_SFR = {
    .z_obs = z,
    .Mdrop = MassTurnover,
    .pl_star = Alpha_star,
    .pl_esc = Alpha_esc,
    .frac_star = Fstar10,
    .frac_esc = Fesc10,
    .LimitMass_Fstar = Mlim_Fstar,
    .LimitMass_Fesc = Mlim_Fesc,
    };
    F.function = &dFdlnM_st_SFR_Lya;
    F.params = &parameters_gsl_SFR;
    //lower_limit = log(M_Min);
    //upper_limit = log(FMAX(1e16, M_Min*100));
    gsl_integration_qag (&F, log(1.0e1), log(1.0e16), 0, rel_tol, 1000, GSL_INTEG_GAUSS61, w, &result, &error);
    gsl_integration_workspace_free (w);
    return result / (OMm*RHOcrit);
}

// - RM

//void initialiseGL_zetaSFR(int n, float M_Min, float M_Max)
//{
//  gauleg(log(M_Min),log(M_Max),xi_SFR_zeta,wi_SFR_zeta,n);
//}

/*
    double z_obs;
    double Mval;
    double delta1;
    double delta2;
    double Mdrop;
    double pl_star;
    double pl_esc;
  double frac_star;
  double frac_esc;
  double LimitMass_Fstar;
  double LimitMass_Fesc;
*/

double zeta_denom_integrand(double lnM, void *params)
{
  struct parameters_gsl_SFR_con_int_ vals = *(struct parameters_gsl_SFR_con_int_ *)params;
      double M = exp(lnM); // linear scale
    double z = vals.z_obs;
    double M2 = vals.Mval; // natural log scale
    double del1 = vals.delta1;
    double del2 = vals.delta2;
    double MassTurnover = vals.Mdrop;
    double Alpha_star = vals.pl_star;
    double Alpha_esc = vals.pl_esc;
  double Fstar10 = vals.frac_star;
  double Fesc10 = vals.frac_esc;

  return M*dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI);
}

double zeta_denom(double z, double M1, double M2, double delta1, double delta2, int flag)
{
  double result, error, lower_limit, upper_limit, Mmin;
  gsl_function F;
  double rel_tol = 0.005;
  gsl_integration_workspace * w
  = gsl_integration_workspace_alloc (1000);

  // jsut setting things we don't care about to 0 so we don't have to make another struct
  struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
      .z_obs = z,
      .Mval = M2,
      .delta1 = delta1,
      .delta2 = delta2,
      .Mdrop = 0,
      .pl_star = 0,
      .pl_esc = 0,
      .frac_star = 0,
      .frac_esc = 0,
      .LimitMass_Fstar = 0,
      .LimitMass_Fesc = 0
    };

    sources src;
    src = defaultSources();
    Mmin = src.minMass(z);


    F.function = &zeta_denom_integrand;
    F.params = &parameters_gsl_SFR_con;
    lower_limit = log(Mmin);
    upper_limit = M2;

    gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                        1000, GSL_INTEG_GAUSS61, w, &result, &error);
    //testing...
    gsl_integration_workspace_free (w);
    //if(flag==1) printf("%le %le %le %le %le\n", z, log(M1), M2, delta1, delta2);
    //if(flag==1) printf("res = %le\n", result);
    return result;
}

float zeta_GL_SFR(float lnM, struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con)
{
      float M = exp(lnM);
    float z = parameters_gsl_SFR_con.z_obs;
    float M2 = parameters_gsl_SFR_con.Mval;
    float del1 = parameters_gsl_SFR_con.delta1;
    float del2 = parameters_gsl_SFR_con.delta2;
    float MassTurnover = parameters_gsl_SFR_con.Mdrop;
    float Alpha_star = parameters_gsl_SFR_con.pl_star;
    float Alpha_esc = parameters_gsl_SFR_con.pl_esc;
  float Fstar10 = parameters_gsl_SFR_con.frac_star;
  float Fesc10 = parameters_gsl_SFR_con.frac_esc;

  float Fstar = pow(M/1e10,Alpha_star);
  float Fesc = pow(M/1e10,Alpha_esc);

  if (Fstar > 1.) Fstar = 1./Fstar10;
  else if (Fstar < 0.) Fstar = FRACT_FLOAT_ERR;

  if (Fesc > 1.) Fesc = 1./Fesc10;
  else if (Fstar < 0.) Fesc = FRACT_FLOAT_ERR;

  sources s;
  s = defaultSources();

  //printf("%e %e %e %e %e", z, M, M2, del1, del2);
  return zeta_at_M(s, z, M) * M*dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI);
}


float GaussLegendreQuad_zetaSFR(int n, float z, float M2, float delta1, float delta2, float MassTurnover, float Alpha_star, float Alpha_esc, float Fstar10, float Fesc10)
{
      //Performs the Gauss-Legendre quadrature.
    int i;
    float integrand, x, denom_integrand;
    integrand = 0.;
    denom_integrand = 0.;
    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
        .z_obs = z,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
        .pl_esc = Alpha_esc,
    .frac_star = Fstar10,
    .frac_esc = Fesc10
    };

    if(delta2 > delta1){
        return 1.;
    }
    else{
        for(i=1; i<(n+1); i++){
            x = xi_SFR[i];
            x = xi_SFR[i];
            integrand += wi_SFR[i]*zeta_GL_SFR(x,parameters_gsl_SFR_con);
        }
        // MassTurnover here is really Mmin
        return integrand / ((float) zeta_denom((double) z, (double) MassTurnover, (double) M2, (double) delta1, (double) delta2, 0));
    }

}

double dzetalnM_SFR(double lnM, void *params)
{
  struct parameters_gsl_SFR_con_int_ vals = *(struct parameters_gsl_SFR_con_int_ *)params;
    double M = exp(lnM); // linear scale
    double z = vals.z_obs;
    double M2 = vals.Mval; // natural log scale
    double del1 = vals.delta1;
    double del2 = vals.delta2;
    double MassTurnover = vals.Mdrop;
    double Alpha_star = vals.pl_star;
    double Alpha_esc = vals.pl_esc;
  double Fstar10 = vals.frac_star;
  double Fesc10 = vals.frac_esc;

  double Fstar = pow(M/1e10,Alpha_star);
  double Fesc = pow(M/1e10,Alpha_esc);

  if (Fstar > 1./Fstar10) Fstar = 1./Fstar10;
  else if (Fstar < 0.) Fstar = FRACT_FLOAT_ERR;

  if (Fesc > 1./Fstar10) Fesc = 1./Fesc10;
  else if (Fstar < 0.) Fesc = FRACT_FLOAT_ERR;

  sources s;
  s = defaultSources();

  return zeta_at_M(s, z, M) * M*dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI);

}

double zeta_SFR(double z, double M1, double M2, double delta1, double delta2, double MassTurnover, double Alpha_star, double Alpha_esc, double Fstar10, double Fesc10)
{
    double result, error, lower_limit, upper_limit, res;
    gsl_function F;
    double rel_tol = 0.005; //<- relative tolerance
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);
    //printf("FgtrConditionalM_SFR start\n");

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
        .z_obs = z,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
        .pl_esc = Alpha_esc,
    .frac_star = Fstar10,
    .frac_esc = Fesc10
    };

    F.function = &dzetalnM_SFR;
    F.params = &parameters_gsl_SFR_con;
    lower_limit = M1;
    upper_limit = M2;

    gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                        1000, GSL_INTEG_GAUSS61, w, &result, &error);
    gsl_integration_workspace_free (w);

    //printf("FgtrConditionalM_SFR end\n");
    if(delta2 >= delta1) {
        result = 1.0;
        return result / zeta_denom(z, exp(M1), M2, delta1, delta2, 1);
    }
    else {
        res = result / zeta_denom(z, exp(M1), M2, delta1, delta2, 1);
        return res; //bug is here (in zeta_denom?)
    }
}

// params struct for integration
struct parameters_gsl_ionEff_int_{
    sources s;
    double z;
};

//integrand for ionEff
// returns dn/dlnM * fesc * fstar * Nion * A_He
// average case
double ionEffIntegrand(double lnM, void *params)
{
    sources s;
    double A_He, z, M;
    struct parameters_gsl_ionEff_int_ vals =
        *(struct parameters_gsl_ionEff_int_ *)params;
    s = vals.s;
    z = vals.z;
    M = exp(lnM);
    A_He = 1.22;
    // this is a Pop II star forming halo
    if(M >= s.minMass(z))
        return A_He * s.fstar(z, M) * s.fesc(z, M) * s.Nion(z, M)
               * dNdM_st(z, M) * M;
    // this halo is below the minimum mass
    return 0.0;
}

// integrand for integral over mass function
double massFuncInt(double lnM, void *params)
{
    double z, M;
    z = *(double *)params;
    M = exp(lnM);
    return dNdM_st(z, M) * M;
}


// integrate over the mass function to calculate average ionizing efficiency
// for the given source popultion
double ionEff(double z, sources s)
{
    double A_He, rel_tol, lower_limit, upper_limit;
    double resultEff, resultMf, error, test;
    gsl_function F;

    rel_tol = 0.001;

    lower_limit = log(s.minMass(z));
    upper_limit = log(FMAX(1.0e16, 100*exp(lower_limit)));

    gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);

    struct parameters_gsl_ionEff_int_ parameters_gsl_ionEff = {
        .s = s,
        .z = z,
    };

    F.function = &ionEffIntegrand;
    F.params = &parameters_gsl_ionEff;

    gsl_integration_qag(&F, lower_limit, upper_limit, 0, rel_tol, 1000,
                        GSL_INTEG_GAUSS61, w, &resultEff, &error);

    F.function = &massFuncInt;
    F.params = &z;


    gsl_integration_qag(&F, lower_limit, upper_limit, 0, rel_tol, 1000,
                        GSL_INTEG_GAUSS61, w, &resultMf, &error);

    gsl_integration_workspace_free(w);

    return resultEff / resultMf;
}

void print_arr(float* arr, int N)
{
  int i;
  for(i=0; i<N; i++) printf("%e\n", arr[i]);
}

void initialisezetaSFR_spline(float z, float Mmax, float MassTurnover, float Alpha_star, float Alpha_esc, float Fstar10, float Fesc10, float Mlim_Fstar, float Mlim_Fesc)
{
    //double overdense_val,overdense_small_low,overdense_small_high,overdense_large_low,overdense_large_high;
    double overdense_val;
    double overdense_large_high = Deltac, overdense_large_low = 1.5;
    double overdense_small_high = 1.5, overdense_small_low = -1. + 9e-8;
    int i;
    zetaLow_spline_acc = gsl_interp_accel_alloc ();
    zetaLow_spline = gsl_spline_alloc (gsl_interp_cspline, NSFR_low);

    for (i=0; i<NSFR_low; i++){

        overdense_val = log10(1. + overdense_small_low) + (double)i/((double)NSFR_low-1.)*(log10(1.+overdense_small_high)-log10(1.+overdense_small_low));

        log10_overdense_spline_SFR[i] = overdense_val;
        log10_zeta_spline_SFR[i] = log10(GaussLegendreQuad_zetaSFR(NGL_SFR,z,log(Mmax),Deltac,pow(10.,overdense_val)-1.,MassTurnover,Alpha_star,Alpha_esc,Fstar10,Fesc10));


        if(log10_zeta_spline_SFR[i] < -40.){
            log10_zeta_spline_SFR[i] = -40.;
        }
    }
    //printf("LOW ZETA\n");
    //for(i=0; i<NSFR_low; i++) printf("%le\n", log10_zeta_spline_SFR[i]);
    gsl_spline_init(zetaLow_spline, log10_overdense_spline_SFR, log10_zeta_spline_SFR, NSFR_low);
    for(i=0;i<NSFR_high;i++) {

        Overdense_spline_SFR[i] = overdense_large_low + (float)i/((float)NSFR_high-1.)*(overdense_large_high - overdense_large_low);
        zeta_spline_SFR[i] = zeta_SFR(z,log(MassTurnover),log(Mmax),Deltac,Overdense_spline_SFR[i],MassTurnover,Alpha_star,Alpha_esc,Fstar10,Fesc10);
        if(i == NSFR_high - 1) zeta_spline_SFR[i] = zeta_spline_SFR[i-1]; // this is wrong, but get nan otherwise, need to fix this
        if(zeta_spline_SFR[i]<0.) {
            zeta_spline_SFR[i]=pow(10.,-40.0);
        }
    }
    //printf("HIGH ZETA\n");
    //print_arr(zeta_spline_SFR, NSFR_high);
    spline(Overdense_spline_SFR-1,zeta_spline_SFR-1,NSFR_high,0,0,zeta_second_derivs_SFR-1);
}

void zetaSpline_SFR(float Overdensity, float *splined_value)
{
    int i;
    float returned_value;

    if(Overdensity<1.5) {
        if(Overdensity<-1.) {
            returned_value = 0;
        }
        else {
            returned_value = gsl_spline_eval(zetaLow_spline, log10(Overdensity+1.), zetaLow_spline_acc);
            returned_value = pow(10.,returned_value);
        }
    }
    else {
        if(Overdensity<0.99*Deltac) {
            splint(Overdense_spline_SFR-1,zeta_spline_SFR-1,zeta_second_derivs_SFR-1,NSFR_high,Overdensity,&(returned_value));
        }
        else {
            returned_value = 1.;
        }
    }
  //if(returned_value > 1.) returned_value = 1.;
    *splined_value = returned_value;

}

void initialiseGL_FcollSFR(int n, float M_TURN, float M_Max, float z){
	float M_Min = M_TURN/50.;
  if(USE_GENERAL_SOURCES)
  {
    sources s;
    s = defaultSources();
    M_Min = s.minMass(z);
  }
	//float M_Min = 1e8;
    //calculates the weightings and the positions for Gauss-Legendre quadrature.
    gauleg(log(M_Min),log(M_Max),xi_SFR,wi_SFR,n);

}

float FgtrConditionallnM_GL_SFR(float lnM, struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con){
    float M = exp(lnM);
    float z = parameters_gsl_SFR_con.z_obs;
    float M2 = parameters_gsl_SFR_con.Mval;
    float del1 = parameters_gsl_SFR_con.delta1;
    float del2 = parameters_gsl_SFR_con.delta2;
    float MassTurnover = parameters_gsl_SFR_con.Mdrop;
    float Alpha_star = parameters_gsl_SFR_con.pl_star;
    float Alpha_esc = parameters_gsl_SFR_con.pl_esc;
	float Fstar10 = parameters_gsl_SFR_con.frac_star;
	float Fesc10 = parameters_gsl_SFR_con.frac_esc;
	float Mlim_Fstar = parameters_gsl_SFR_con.LimitMass_Fstar;
	float Mlim_Fesc = parameters_gsl_SFR_con.LimitMass_Fesc;

	float Fstar,Fesc;

	if (Alpha_star > 0. && M > Mlim_Fstar)
		Fstar = 1./Fstar10;
	else if (Alpha_star < 0. && M < Mlim_Fstar)
		Fstar = 1./Fstar10;
	else
		Fstar = pow(M/1e10,Alpha_star);

	if (Alpha_esc > 0. && M > Mlim_Fesc)
		Fesc = 1./Fesc10;
	else if (Alpha_esc < 0. && M < Mlim_Fesc)
		Fesc = 1./Fesc10;
	else
		Fesc = pow(M/1e10,Alpha_esc);

  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    return M * dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI) * src.fstar(z, M) * src.fesc(z, M) * src.Nion(z, M);
  }

    return M*exp(-MassTurnover/M)*Fstar*Fesc*dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI);

}

float GaussLegendreQuad_FcollSFR(int n, float z, float M2, float delta1, float delta2, float MassTurnover, float Alpha_star, float Alpha_esc, float Fstar10, float Fesc10, float Mlim_Fstar, float Mlim_Fesc) {
    //Performs the Gauss-Legendre quadrature.
    int i;

    float integrand, x;
    integrand = 0.;

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
        .z_obs = z,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
        .pl_esc = Alpha_esc,
		.frac_star = Fstar10,
		.frac_esc = Fesc10,
		.LimitMass_Fstar = Mlim_Fstar,
		.LimitMass_Fesc = Mlim_Fesc
    };

    if(delta2 > delta1){
        return 1.;
    }
    else{
        for(i=1; i<(n+1); i++){
            x = xi_SFR[i];
            x = xi_SFR[i];
            integrand += wi_SFR[i]*FgtrConditionallnM_GL_SFR(x,parameters_gsl_SFR_con);
        }
        return integrand;
    }

}

float FgtrConditionallnM_GL_SFRPOP2(float lnM, struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con){
    float M = exp(lnM);
    float z = parameters_gsl_SFR_con.z_obs;
    float M2 = parameters_gsl_SFR_con.Mval;
    float del1 = parameters_gsl_SFR_con.delta1;
    float del2 = parameters_gsl_SFR_con.delta2;
    float MassTurnover = parameters_gsl_SFR_con.Mdrop;
    float Alpha_star = parameters_gsl_SFR_con.pl_star;
    float Alpha_esc = parameters_gsl_SFR_con.pl_esc;
  	float Fstar10 = parameters_gsl_SFR_con.frac_star;
  	float Fesc10 = parameters_gsl_SFR_con.frac_esc;
  	float Mlim_Fstar = parameters_gsl_SFR_con.LimitMass_Fstar;
  	float Mlim_Fesc = parameters_gsl_SFR_con.LimitMass_Fesc;

  	float Fstar,Fesc;

    if(USE_GENERAL_SOURCES)
    {
      sources src;
      src = defaultSources();
      return dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI) * src.sfr(z, M);
    }

    return 0;

}

float FgtrConditionallnM_GL_gas(float lnM, struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con){
    float M = exp(lnM);
    float z = parameters_gsl_SFR_con.z_obs;
    float M2 = parameters_gsl_SFR_con.Mval;
    float del1 = parameters_gsl_SFR_con.delta1;
    float del2 = parameters_gsl_SFR_con.delta2;
    float MassTurnover = parameters_gsl_SFR_con.Mdrop;
    float Alpha_star = parameters_gsl_SFR_con.pl_star;
    float Alpha_esc = parameters_gsl_SFR_con.pl_esc;
    float Fstar10 = parameters_gsl_SFR_con.frac_star;
    float Fesc10 = parameters_gsl_SFR_con.frac_esc;
    float Mlim_Fstar = parameters_gsl_SFR_con.LimitMass_Fstar;
    float Mlim_Fesc = parameters_gsl_SFR_con.LimitMass_Fesc;

    float Fstar,Fesc;

    if(USE_GENERAL_SOURCES)
    {
      sources src;
      src = defaultSources();
      return dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI) * src.Mg(z, M);
    }

    return 0;

}

float FgtrConditionallnM_GL_met(float lnM, struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con){
    float M = exp(lnM);
    float z = parameters_gsl_SFR_con.z_obs;
    float M2 = parameters_gsl_SFR_con.Mval;
    float del1 = parameters_gsl_SFR_con.delta1;
    float del2 = parameters_gsl_SFR_con.delta2;
    float MassTurnover = parameters_gsl_SFR_con.Mdrop;
    float Alpha_star = parameters_gsl_SFR_con.pl_star;
    float Alpha_esc = parameters_gsl_SFR_con.pl_esc;
    float Fstar10 = parameters_gsl_SFR_con.frac_star;
    float Fesc10 = parameters_gsl_SFR_con.frac_esc;
    float Mlim_Fstar = parameters_gsl_SFR_con.LimitMass_Fstar;
    float Mlim_Fesc = parameters_gsl_SFR_con.LimitMass_Fesc;

    float Fstar,Fesc;

    if(USE_GENERAL_SOURCES)
    {
      sources src;
      src = defaultSources();
      return dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI) * src.Mz(z, M);
    }

    return 0;

}

float GaussLegendreQuad_FcollSFRPOP2(int n, float z, float M2, float delta1, float delta2, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar) {
    //Performs the Gauss-Legendre quadrature.
    int i;

    float integrand, x;
    integrand = 0.;

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
      .z_obs = z,
      .Mval = M2,
      .delta1 = delta1,
      .delta2 = delta2,
      .Mdrop = MassTurnover,
      .pl_star = Alpha_star,
  		.frac_star = Fstar10,
  		.LimitMass_Fstar = Mlim_Fstar,
    };

    if(delta2 > delta1){
        return 1.;
    }
    else{
        for(i=1; i<(n+1); i++){
            x = xi_SFR[i];
            x = xi_SFR[i];
            integrand += wi_SFR[i]*FgtrConditionallnM_GL_SFRPOP2(x,parameters_gsl_SFR_con);
        }
        return integrand;
    }

}

float GaussLegendreQuad_Fcollgas(int n, float z, float M2, float delta1, float delta2, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar) {
    //Performs the Gauss-Legendre quadrature.
    int i;

    float integrand, x;
    integrand = 0.;

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
      .z_obs = z,
      .Mval = M2,
      .delta1 = delta1,
      .delta2 = delta2,
      .Mdrop = MassTurnover,
      .pl_star = Alpha_star,
      .frac_star = Fstar10,
      .LimitMass_Fstar = Mlim_Fstar,
    };

    if(delta2 > delta1){
        return 1.;
    }
    else{
        for(i=1; i<(n+1); i++){
            x = xi_SFR[i];
            x = xi_SFR[i];
            integrand += wi_SFR[i]*FgtrConditionallnM_GL_gas(x,parameters_gsl_SFR_con);
        }
        return integrand;
    }

}

float GaussLegendreQuad_Fcollmet(int n, float z, float M2, float delta1, float delta2, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar) {
    //Performs the Gauss-Legendre quadrature.
    int i;

    float integrand, x;
    integrand = 0.;

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
      .z_obs = z,
      .Mval = M2,
      .delta1 = delta1,
      .delta2 = delta2,
      .Mdrop = MassTurnover,
      .pl_star = Alpha_star,
      .frac_star = Fstar10,
      .LimitMass_Fstar = Mlim_Fstar,
    };

    if(delta2 > delta1){
        return 1.;
    }
    else{
        for(i=1; i<(n+1); i++){
            x = xi_SFR[i];
            x = xi_SFR[i];
            integrand += wi_SFR[i]*FgtrConditionallnM_GL_met(x,parameters_gsl_SFR_con);
        }
        return integrand;
    }

}

float FgtrConditionallnM_GL_LINELUMPOP2A(float lnM, struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con){
  float M = exp(lnM);
  float z = parameters_gsl_SFR_con.z_obs;
  float M2 = parameters_gsl_SFR_con.Mval;
  float del1 = parameters_gsl_SFR_con.delta1;
  float del2 = parameters_gsl_SFR_con.delta2;
  float MassTurnover = parameters_gsl_SFR_con.Mdrop;
  float Alpha_star = parameters_gsl_SFR_con.pl_star;
	float Fstar10 = parameters_gsl_SFR_con.frac_star;
	float Mlim_Fstar = parameters_gsl_SFR_con.LimitMass_Fstar;
  int LineID = parameters_gsl_SFR_con.lid;
	float Fstar;

  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    return src.linelum(z, M, -2., LineID) * dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI);
  }
    return 0;
}

float GaussLegendreQuad_FcollLINELUMPOP2A(int n, float z, float M2, float delta1, float delta2, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar, int LineID) {
    //Performs the Gauss-Legendre quadrature.
    int i;

    float integrand, x;
    integrand = 0.;

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
        .z_obs = z,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
		    .frac_star = Fstar10,
        .LimitMass_Fstar = Mlim_Fstar,
        .lid = LineID,
    };

    if(delta2 > delta1){
        return 1.;
    }
    else{
        for(i=1; i<(n+1); i++){
            x = xi_SFR[i];
            x = xi_SFR[i];
            integrand += wi_SFR[i]*FgtrConditionallnM_GL_LINELUMPOP2A(x,parameters_gsl_SFR_con);
        }
        return integrand;
    }

}

float FgtrConditionallnM_GL_LINELUMPOP2B(float lnM, struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con){
  float M = exp(lnM);
  float z = parameters_gsl_SFR_con.z_obs;
  float M2 = parameters_gsl_SFR_con.Mval;
  float del1 = parameters_gsl_SFR_con.delta1;
  float del2 = parameters_gsl_SFR_con.delta2;
  float MassTurnover = parameters_gsl_SFR_con.Mdrop;
  float Alpha_star = parameters_gsl_SFR_con.pl_star;
	float Fstar10 = parameters_gsl_SFR_con.frac_star;
	float Mlim_Fstar = parameters_gsl_SFR_con.LimitMass_Fstar;
  int LineID = parameters_gsl_SFR_con.lid;
	float Fstar;

  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    return src.linelum(z, M, -2., LineID) * dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI);
  }
    return 0;
}

float GaussLegendreQuad_FcollLINELUMPOP2B(int n, float z, float M2, float delta1, float delta2, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar, int LineID) {
    //Performs the Gauss-Legendre quadrature.
    int i;

    float integrand, x;
    integrand = 0.;

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
        .z_obs = z,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
		    .frac_star = Fstar10,
        .LimitMass_Fstar = Mlim_Fstar,
        .lid = LineID,
    };

    if(delta2 > delta1){
        return 1.;
    }
    else{
        for(i=1; i<(n+1); i++){
            x = xi_SFR[i];
            x = xi_SFR[i];
            integrand += wi_SFR[i]*FgtrConditionallnM_GL_LINELUMPOP2B(x,parameters_gsl_SFR_con);
        }
        return integrand;
    }

}

float FgtrConditionallnM_GL_LINELUMPOP2C(float lnM, struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con){
  float M = exp(lnM);
  float z = parameters_gsl_SFR_con.z_obs;
  float M2 = parameters_gsl_SFR_con.Mval;
  float del1 = parameters_gsl_SFR_con.delta1;
  float del2 = parameters_gsl_SFR_con.delta2;
  float MassTurnover = parameters_gsl_SFR_con.Mdrop;
  float Alpha_star = parameters_gsl_SFR_con.pl_star;
	float Fstar10 = parameters_gsl_SFR_con.frac_star;
	float Mlim_Fstar = parameters_gsl_SFR_con.LimitMass_Fstar;
  int LineID = parameters_gsl_SFR_con.lid;
	float Fstar;

  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    return src.linelum(z, M, -2., LineID) * dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI);
  }
    return 0;
}

float GaussLegendreQuad_FcollLINELUMPOP2C(int n, float z, float M2, float delta1, float delta2, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar, int LineID) {
    //Performs the Gauss-Legendre quadrature.
    int i;

    float integrand, x;
    integrand = 0.;

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
        .z_obs = z,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
		    .frac_star = Fstar10,
        .LimitMass_Fstar = Mlim_Fstar,
        .lid = LineID,
    };

    if(delta2 > delta1){
        return 1.;
    }
    else{
        for(i=1; i<(n+1); i++){
            x = xi_SFR[i];
            x = xi_SFR[i];
            integrand += wi_SFR[i]*FgtrConditionallnM_GL_LINELUMPOP2C(x,parameters_gsl_SFR_con);
        }
        return integrand;
    }

}

float FgtrConditionallnM_GL_LINELUMPOP2D(float lnM, struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con){
  float M = exp(lnM);
  float z = parameters_gsl_SFR_con.z_obs;
  float M2 = parameters_gsl_SFR_con.Mval;
  float del1 = parameters_gsl_SFR_con.delta1;
  float del2 = parameters_gsl_SFR_con.delta2;
  float MassTurnover = parameters_gsl_SFR_con.Mdrop;
  float Alpha_star = parameters_gsl_SFR_con.pl_star;
	float Fstar10 = parameters_gsl_SFR_con.frac_star;
	float Mlim_Fstar = parameters_gsl_SFR_con.LimitMass_Fstar;
  int LineID = parameters_gsl_SFR_con.lid;
	float Fstar;

  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    return src.linelum(z, M, -2., LineID) * dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI);
  }
    return 0;
}

float GaussLegendreQuad_FcollLINELUMPOP2D(int n, float z, float M2, float delta1, float delta2, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar, int LineID) {
    //Performs the Gauss-Legendre quadrature.
    int i;

    float integrand, x;
    integrand = 0.;

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
        .z_obs = z,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
		    .frac_star = Fstar10,
        .LimitMass_Fstar = Mlim_Fstar,
        .lid = LineID,
    };

    if(delta2 > delta1){
        return 1.;
    }
    else{
        for(i=1; i<(n+1); i++){
            x = xi_SFR[i];
            x = xi_SFR[i];
            integrand += wi_SFR[i]*FgtrConditionallnM_GL_LINELUMPOP2D(x,parameters_gsl_SFR_con);
        }
        return integrand;
    }

}

float FgtrConditionallnM_GL_LINELUMPOP2E(float lnM, struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con){
  float M = exp(lnM);
  float z = parameters_gsl_SFR_con.z_obs;
  float M2 = parameters_gsl_SFR_con.Mval;
  float del1 = parameters_gsl_SFR_con.delta1;
  float del2 = parameters_gsl_SFR_con.delta2;
  float MassTurnover = parameters_gsl_SFR_con.Mdrop;
  float Alpha_star = parameters_gsl_SFR_con.pl_star;
	float Fstar10 = parameters_gsl_SFR_con.frac_star;
	float Mlim_Fstar = parameters_gsl_SFR_con.LimitMass_Fstar;
  int LineID = parameters_gsl_SFR_con.lid;
	float Fstar;

  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    return src.linelum(z, M, -2., LineID) * dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI);
  }
    return 0;
}

float GaussLegendreQuad_FcollLINELUMPOP2E(int n, float z, float M2, float delta1, float delta2, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar, int LineID) {
    //Performs the Gauss-Legendre quadrature.
    int i;

    float integrand, x;
    integrand = 0.;

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
        .z_obs = z,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
		    .frac_star = Fstar10,
        .LimitMass_Fstar = Mlim_Fstar,
        .lid = LineID,
    };

    if(delta2 > delta1){
        return 1.;
    }
    else{
        for(i=1; i<(n+1); i++){
            x = xi_SFR[i];
            x = xi_SFR[i];
            integrand += wi_SFR[i]*FgtrConditionallnM_GL_LINELUMPOP2E(x,parameters_gsl_SFR_con);
        }
        return integrand;
    }

}

float FgtrConditionallnM_GL_LINELUMPOP2F(float lnM, struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con){
  float M = exp(lnM);
  float z = parameters_gsl_SFR_con.z_obs;
  float M2 = parameters_gsl_SFR_con.Mval;
  float del1 = parameters_gsl_SFR_con.delta1;
  float del2 = parameters_gsl_SFR_con.delta2;
  float MassTurnover = parameters_gsl_SFR_con.Mdrop;
  float Alpha_star = parameters_gsl_SFR_con.pl_star;
	float Fstar10 = parameters_gsl_SFR_con.frac_star;
	float Mlim_Fstar = parameters_gsl_SFR_con.LimitMass_Fstar;
  int LineID = parameters_gsl_SFR_con.lid;
	float Fstar;

  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    return src.linelum(z, M, -2., LineID) * dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI);
  }
    return 0;
}

float GaussLegendreQuad_FcollLINELUMPOP2F(int n, float z, float M2, float delta1, float delta2, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar, int LineID) {
    //Performs the Gauss-Legendre quadrature.
    int i;

    float integrand, x;
    integrand = 0.;

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
        .z_obs = z,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
		    .frac_star = Fstar10,
        .LimitMass_Fstar = Mlim_Fstar,
        .lid = LineID,
    };

    if(delta2 > delta1){
        return 1.;
    }
    else{
        for(i=1; i<(n+1); i++){
            x = xi_SFR[i];
            x = xi_SFR[i];
            integrand += wi_SFR[i]*FgtrConditionallnM_GL_LINELUMPOP2F(x,parameters_gsl_SFR_con);
        }
        return integrand;
    }

}

float FgtrConditionallnM_GL_LINELUMPOP2G(float lnM, struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con){
  float M = exp(lnM);
  float z = parameters_gsl_SFR_con.z_obs;
  float M2 = parameters_gsl_SFR_con.Mval;
  float del1 = parameters_gsl_SFR_con.delta1;
  float del2 = parameters_gsl_SFR_con.delta2;
  float MassTurnover = parameters_gsl_SFR_con.Mdrop;
  float Alpha_star = parameters_gsl_SFR_con.pl_star;
	float Fstar10 = parameters_gsl_SFR_con.frac_star;
	float Mlim_Fstar = parameters_gsl_SFR_con.LimitMass_Fstar;
  int LineID = parameters_gsl_SFR_con.lid;
	float Fstar;

  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    return src.linelum(z, M, -2., LineID) * dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI);
  }
    return 0;
}

float GaussLegendreQuad_FcollLINELUMPOP2G(int n, float z, float M2, float delta1, float delta2, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar, int LineID) {
    //Performs the Gauss-Legendre quadrature.
    int i;

    float integrand, x;
    integrand = 0.;

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
        .z_obs = z,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
		    .frac_star = Fstar10,
        .LimitMass_Fstar = Mlim_Fstar,
        .lid = LineID,
    };

    if(delta2 > delta1){
        return 1.;
    }
    else{
        for(i=1; i<(n+1); i++){
            x = xi_SFR[i];
            x = xi_SFR[i];
            integrand += wi_SFR[i]*FgtrConditionallnM_GL_LINELUMPOP2G(x,parameters_gsl_SFR_con);
        }
        return integrand;
    }

}

// same as above, but with fx and fesc = 1 and no Nion
float FgtrConditionallnM_GL_SFR_Xray(float lnM, struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con){
    float M = exp(lnM);
    float z = parameters_gsl_SFR_con.z_obs;
    float M2 = parameters_gsl_SFR_con.Mval;
    float del1 = parameters_gsl_SFR_con.delta1;
    float del2 = parameters_gsl_SFR_con.delta2;
    float MassTurnover = parameters_gsl_SFR_con.Mdrop;
    float Alpha_star = parameters_gsl_SFR_con.pl_star;
    float Alpha_esc = parameters_gsl_SFR_con.pl_esc;
  float Fstar10 = parameters_gsl_SFR_con.frac_star;
  float Fesc10 = parameters_gsl_SFR_con.frac_esc;
  float Mlim_Fstar = parameters_gsl_SFR_con.LimitMass_Fstar;
  float Mlim_Fesc = parameters_gsl_SFR_con.LimitMass_Fesc;

  float Fstar,Fesc;

  if (Alpha_star > 0. && M > Mlim_Fstar)
    Fstar = 1./Fstar10;
  else if (Alpha_star < 0. && M < Mlim_Fstar)
    Fstar = 1./Fstar10;
  else
    Fstar = pow(M/1e10,Alpha_star);

  if (Alpha_esc > 0. && M > Mlim_Fesc)
    Fesc = 1./Fesc10;
  else if (Alpha_esc < 0. && M < Mlim_Fesc)
    Fesc = 1./Fesc10;
  else
    Fesc = pow(M/1e10,Alpha_esc);

  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    //return M * dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI) * src.fstar(z, M) * src.fx(z, M);
    return dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI) * (src.sfr(z, M)/f_baryon)/SperYR * src.fx(z, M);
    //return M * dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI) * src.fstar(z, M) * hubble(z) / 0.5 * src.fx(z, M);
  }

    return M*exp(-MassTurnover/M)*Fstar*Fesc*dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI);

}

float GaussLegendreQuad_FcollSFR_Xray(int n, float z, float M2, float delta1, float delta2, float MassTurnover, float Alpha_star, float Alpha_esc, float Fstar10, float Fesc10, float Mlim_Fstar, float Mlim_Fesc) {
    //Performs the Gauss-Legendre quadrature.
    int i;

    float integrand, x;
    integrand = 0.;

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
        .z_obs = z,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
        .pl_esc = Alpha_esc,
    .frac_star = Fstar10,
    .frac_esc = Fesc10,
    .LimitMass_Fstar = Mlim_Fstar,
    .LimitMass_Fesc = Mlim_Fesc
    };

    if(delta2 > delta1){
        return 1.;
    }
    else{
        for(i=1; i<(n+1); i++){
            x = xi_SFR[i];
            x = xi_SFR[i];
            integrand += wi_SFR[i]*FgtrConditionallnM_GL_SFR_Xray(x,parameters_gsl_SFR_con);
        }
        return integrand;
    }

}

// same as above, but with fx and fesc = 1 and no Nion and no fx
float FgtrConditionallnM_GL_SFR_Lya(float lnM, struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con){
    float M = exp(lnM);
    float z = parameters_gsl_SFR_con.z_obs;
    float M2 = parameters_gsl_SFR_con.Mval;
    float del1 = parameters_gsl_SFR_con.delta1;
    float del2 = parameters_gsl_SFR_con.delta2;
    float MassTurnover = parameters_gsl_SFR_con.Mdrop;
    float Alpha_star = parameters_gsl_SFR_con.pl_star;
    float Alpha_esc = parameters_gsl_SFR_con.pl_esc;
  float Fstar10 = parameters_gsl_SFR_con.frac_star;
  float Fesc10 = parameters_gsl_SFR_con.frac_esc;
  float Mlim_Fstar = parameters_gsl_SFR_con.LimitMass_Fstar;
  float Mlim_Fesc = parameters_gsl_SFR_con.LimitMass_Fesc;

  float Fstar,Fesc;

  if (Alpha_star > 0. && M > Mlim_Fstar)
    Fstar = 1./Fstar10;
  else if (Alpha_star < 0. && M < Mlim_Fstar)
    Fstar = 1./Fstar10;
  else
    Fstar = pow(M/1e10,Alpha_star);

  if (Alpha_esc > 0. && M > Mlim_Fesc)
    Fesc = 1./Fesc10;
  else if (Alpha_esc < 0. && M < Mlim_Fesc)
    Fesc = 1./Fesc10;
  else
    Fesc = pow(M/1e10,Alpha_esc);

  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    //return M * dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI) * src.fstar(z, M);
    return dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI) * (src.sfr(z, M)/f_baryon)/SperYR;
    //return M * dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI) * src.fstar(z, M) * hubble(z) / 0.5;
  }

    return M*exp(-MassTurnover/M)*Fstar*Fesc*dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI);

}

float GaussLegendreQuad_FcollSFR_Lya(int n, float z, float M2, float delta1, float delta2, float MassTurnover, float Alpha_star, float Alpha_esc, float Fstar10, float Fesc10, float Mlim_Fstar, float Mlim_Fesc) {
    //Performs the Gauss-Legendre quadrature.
    int i;

    float integrand, x;
    integrand = 0.;

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
        .z_obs = z,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
        .pl_esc = Alpha_esc,
    .frac_star = Fstar10,
    .frac_esc = Fesc10,
    .LimitMass_Fstar = Mlim_Fstar,
    .LimitMass_Fesc = Mlim_Fesc
    };

    if(delta2 > delta1){
        return 1.;
    }
    else{
        for(i=1; i<(n+1); i++){
            x = xi_SFR[i];
            x = xi_SFR[i];
            integrand += wi_SFR[i]*FgtrConditionallnM_GL_SFR_Lya(x,parameters_gsl_SFR_con);
        }
        return integrand;
    }

}

double dFgtrConditionallnM_SFR(double lnM, void *params) {
    struct parameters_gsl_SFR_con_int_ vals = *(struct parameters_gsl_SFR_con_int_ *)params;
    double M = exp(lnM); // linear scale
    double z = vals.z_obs;
    double M2 = vals.Mval; // natural log scale
    double del1 = vals.delta1;
    double del2 = vals.delta2;
    double MassTurnover = vals.Mdrop;
    double Alpha_star = vals.pl_star;
    double Alpha_esc = vals.pl_esc;
	double Fstar10 = vals.frac_star;
	double Fesc10 = vals.frac_esc;
	double Mlim_Fstar = vals.LimitMass_Fstar;
	double Mlim_Fesc = vals.LimitMass_Fesc;

	double Fstar,Fesc;

	if (Alpha_star > 0. && M > Mlim_Fstar)
		Fstar = 1./Fstar10;
	else if (Alpha_star < 0. && M < Mlim_Fstar)
		Fstar = 1./Fstar10;
	else
		Fstar = pow(M/1e10,Alpha_star);

	if (Alpha_esc > 0. && M > Mlim_Fesc)
		Fesc = 1./Fesc10;
	else if (Alpha_esc < 0. && M < Mlim_Fesc)
		Fesc = 1./Fesc10;
	else
		Fesc = pow(M/1e10,Alpha_esc);

  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    return M * dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI) * src.fstar(z, M) * src.fesc(z, M) * src.Nion(z, M);
  }

    return M*exp(-MassTurnover/M)*Fstar*Fesc*dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI);
}

double FgtrConditionalM_SFR(double z, double M1, double M2, double delta1, double delta2, double MassTurnover, double Alpha_star, double Alpha_esc, double Fstar10, double Fesc10, double Mlim_Fstar, double Mlim_Fesc) {
    double result, error, lower_limit, upper_limit;
    gsl_function F;
    double rel_tol = 0.005; //<- relative tolerance
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
        .z_obs = z,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
        .pl_esc = Alpha_esc,
		.frac_star = Fstar10,
		.frac_esc = Fesc10,
		.LimitMass_Fstar = Mlim_Fstar,
		.LimitMass_Fesc = Mlim_Fesc
    };

    F.function = &dFgtrConditionallnM_SFR;
    F.params = &parameters_gsl_SFR_con;
    lower_limit = M1;
    upper_limit = M2;

    if(USE_GENERAL_SOURCES)
    {
      sources s;
      s = defaultSources();
      lower_limit = log(s.minMass(z));
    }

    gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                        1000, GSL_INTEG_GAUSS61, w, &result, &error);
    gsl_integration_workspace_free (w);

    if(delta2 >= delta1) {
        result = 1.;
        return result;
    }
    else {
        return result;
    }

}

double dFgtrConditionallnM_SFRPOP2(double lnM, void *params) {
    struct parameters_gsl_SFR_con_int_ vals = *(struct parameters_gsl_SFR_con_int_ *)params;
    double M = exp(lnM); // linear scale
    double z = vals.z_obs;
    double M2 = vals.Mval; // natural log scale
    double del1 = vals.delta1;
    double del2 = vals.delta2;
    double MassTurnover = vals.Mdrop;
    double Alpha_star = vals.pl_star;
    double Alpha_esc = vals.pl_esc;
  	double Fstar10 = vals.frac_star;
  	double Fesc10 = vals.frac_esc;
  	double Mlim_Fstar = vals.LimitMass_Fstar;
  	double Mlim_Fesc = vals.LimitMass_Fesc;

	  double Fstar,Fesc;

    if(USE_GENERAL_SOURCES)
    {
      sources src;
      src = defaultSources();
      return dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI) * src.sfr(z, M);
    }

      return 0;
}

double dFgtrConditionallnM_gas(double lnM, void *params) {
    struct parameters_gsl_SFR_con_int_ vals = *(struct parameters_gsl_SFR_con_int_ *)params;
    double M = exp(lnM); // linear scale
    double z = vals.z_obs;
    double M2 = vals.Mval; // natural log scale
    double del1 = vals.delta1;
    double del2 = vals.delta2;
    double MassTurnover = vals.Mdrop;
    double Alpha_star = vals.pl_star;
    double Alpha_esc = vals.pl_esc;
    double Fstar10 = vals.frac_star;
    double Fesc10 = vals.frac_esc;
    double Mlim_Fstar = vals.LimitMass_Fstar;
    double Mlim_Fesc = vals.LimitMass_Fesc;

    double Fstar,Fesc;

    if(USE_GENERAL_SOURCES)
    {
      sources src;
      src = defaultSources();
      return dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI) * src.Mg(z, M);
    }

      return 0;
}

double dFgtrConditionallnM_met(double lnM, void *params) {
    struct parameters_gsl_SFR_con_int_ vals = *(struct parameters_gsl_SFR_con_int_ *)params;
    double M = exp(lnM); // linear scale
    double z = vals.z_obs;
    double M2 = vals.Mval; // natural log scale
    double del1 = vals.delta1;
    double del2 = vals.delta2;
    double MassTurnover = vals.Mdrop;
    double Alpha_star = vals.pl_star;
    double Alpha_esc = vals.pl_esc;
    double Fstar10 = vals.frac_star;
    double Fesc10 = vals.frac_esc;
    double Mlim_Fstar = vals.LimitMass_Fstar;
    double Mlim_Fesc = vals.LimitMass_Fesc;

    double Fstar,Fesc;

    if(USE_GENERAL_SOURCES)
    {
      sources src;
      src = defaultSources();
      return dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI) * src.Mz(z, M);
    }

      return 0;
}

double FgtrConditionalM_SFRPOP2(double z, double M1, double M2, double delta1, double delta2, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar) {
    double result, error, lower_limit, upper_limit;
    gsl_function F;
    double rel_tol = 0.005; //<- relative tolerance
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
      .z_obs = z,
      .Mval = M2,
      .delta1 = delta1,
      .delta2 = delta2,
      .Mdrop = MassTurnover,
      .pl_star = Alpha_star,
  		.frac_star = Fstar10,
  		.LimitMass_Fstar = Mlim_Fstar,
    };

    F.function = &dFgtrConditionallnM_SFRPOP2;
    F.params = &parameters_gsl_SFR_con;
    lower_limit = M1;
    upper_limit = M2;

    if(USE_GENERAL_SOURCES)
    {
      sources s;
      s = defaultSources();
      lower_limit = log(s.minMass(z));
    }

    gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                        1000, GSL_INTEG_GAUSS61, w, &result, &error);
    gsl_integration_workspace_free (w);

    if(delta2 >= delta1) {
        result = 1.;
        return result;
    }
    else {
        return result;
    }

}

double FgtrConditionalM_gas(double z, double M1, double M2, double delta1, double delta2, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar) {
    double result, error, lower_limit, upper_limit;
    gsl_function F;
    double rel_tol = 0.005; //<- relative tolerance
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
      .z_obs = z,
      .Mval = M2,
      .delta1 = delta1,
      .delta2 = delta2,
      .Mdrop = MassTurnover,
      .pl_star = Alpha_star,
      .frac_star = Fstar10,
      .LimitMass_Fstar = Mlim_Fstar,
    };

    F.function = &dFgtrConditionallnM_gas;
    F.params = &parameters_gsl_SFR_con;
    lower_limit = M1;
    upper_limit = M2;

    if(USE_GENERAL_SOURCES)
    {
      sources s;
      s = defaultSources();
      lower_limit = log(s.minMass(z));
    }

    gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                        1000, GSL_INTEG_GAUSS61, w, &result, &error);
    gsl_integration_workspace_free (w);

    if(delta2 >= delta1) {
        result = 1.;
        return result;
    }
    else {
        return result;
    }

}

double FgtrConditionalM_met(double z, double M1, double M2, double delta1, double delta2, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar) {
    double result, error, lower_limit, upper_limit;
    gsl_function F;
    double rel_tol = 0.005; //<- relative tolerance
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
      .z_obs = z,
      .Mval = M2,
      .delta1 = delta1,
      .delta2 = delta2,
      .Mdrop = MassTurnover,
      .pl_star = Alpha_star,
      .frac_star = Fstar10,
      .LimitMass_Fstar = Mlim_Fstar,
    };

    F.function = &dFgtrConditionallnM_met;
    F.params = &parameters_gsl_SFR_con;
    lower_limit = M1;
    upper_limit = M2;

    if(USE_GENERAL_SOURCES)
    {
      sources s;
      s = defaultSources();
      lower_limit = log(s.minMass(z));
    }

    gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                        1000, GSL_INTEG_GAUSS61, w, &result, &error);
    gsl_integration_workspace_free (w);

    if(delta2 >= delta1) {
        result = 1.;
        return result;
    }
    else {
        return result;
    }

}

// used below for LINELUMPOP2A -GS
double dFgtrConditionallnM_LINELUMPOP2A(double lnM, void *params) {
  struct parameters_gsl_SFR_con_int_ vals = *(struct parameters_gsl_SFR_con_int_ *)params;
  double M = exp(lnM); // linear scale
  double z = vals.z_obs;
  double M2 = vals.Mval; // natural log scale
  double del1 = vals.delta1;
  double del2 = vals.delta2;
  double MassTurnover = vals.Mdrop;
  double Alpha_star = vals.pl_star;
	double Fstar10 = vals.frac_star;
	double Mlim_Fstar = vals.LimitMass_Fstar;
  int LineID = vals.lid;
	double Fstar;

  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    return src.linelum(z, M, -2., LineID) * dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI);
  }
    return 0;
}

// used when only taking into account Fstar -GS
double FgtrConditionalM_LINELUMPOP2A(double z, double M1, double M2, double delta1, double delta2, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar, int LineID) {
    double result, error, lower_limit, upper_limit;
    gsl_function F;
    double rel_tol = 0.005; //<- relative tolerance
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
        .z_obs = z,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
		    .frac_star = Fstar10,
		    .LimitMass_Fstar = Mlim_Fstar,
        .lid = LineID,
    };

    F.function = &dFgtrConditionallnM_LINELUMPOP2A;
    F.params = &parameters_gsl_SFR_con;
    lower_limit = M1;
    upper_limit = M2;

    if(USE_GENERAL_SOURCES)
    {
      sources s;
      s = defaultSources();
      lower_limit = log(s.minMass(z));
    }

    gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                        1000, GSL_INTEG_GAUSS61, w, &result, &error);
    gsl_integration_workspace_free (w);

    if(delta2 >= delta1) {
        result = 1.;
        return result;
    }
    else {
        return result;
    }

}

// used below for LINELUMPOP2B -GS
double dFgtrConditionallnM_LINELUMPOP2B(double lnM, void *params) {
  struct parameters_gsl_SFR_con_int_ vals = *(struct parameters_gsl_SFR_con_int_ *)params;
  double M = exp(lnM); // linear scale
  double z = vals.z_obs;
  double M2 = vals.Mval; // natural log scale
  double del1 = vals.delta1;
  double del2 = vals.delta2;
  double MassTurnover = vals.Mdrop;
  double Alpha_star = vals.pl_star;
	double Fstar10 = vals.frac_star;
	double Mlim_Fstar = vals.LimitMass_Fstar;
  int LineID = vals.lid;
	double Fstar;

  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    return src.linelum(z, M, -2., LineID) * dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI);
  }
    return 0;
}

// used when only taking into account Fstar -GS
double FgtrConditionalM_LINELUMPOP2B(double z, double M1, double M2, double delta1, double delta2, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar, int LineID) {
    double result, error, lower_limit, upper_limit;
    gsl_function F;
    double rel_tol = 0.005; //<- relative tolerance
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
        .z_obs = z,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
		    .frac_star = Fstar10,
		    .LimitMass_Fstar = Mlim_Fstar,
        .lid = LineID,
    };

    F.function = &dFgtrConditionallnM_LINELUMPOP2B;
    F.params = &parameters_gsl_SFR_con;
    lower_limit = M1;
    upper_limit = M2;

    if(USE_GENERAL_SOURCES)
    {
      sources s;
      s = defaultSources();
      lower_limit = log(s.minMass(z));
    }

    gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                        1000, GSL_INTEG_GAUSS61, w, &result, &error);
    gsl_integration_workspace_free (w);

    if(delta2 >= delta1) {
        result = 1.;
        return result;
    }
    else {
        return result;
    }

}

// used below for LINELUMPOP2C -GS
double dFgtrConditionallnM_LINELUMPOP2C(double lnM, void *params) {
  struct parameters_gsl_SFR_con_int_ vals = *(struct parameters_gsl_SFR_con_int_ *)params;
  double M = exp(lnM); // linear scale
  double z = vals.z_obs;
  double M2 = vals.Mval; // natural log scale
  double del1 = vals.delta1;
  double del2 = vals.delta2;
  double MassTurnover = vals.Mdrop;
  double Alpha_star = vals.pl_star;
	double Fstar10 = vals.frac_star;
	double Mlim_Fstar = vals.LimitMass_Fstar;
  int LineID = vals.lid;
	double Fstar;

  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    return src.linelum(z, M, -2., LineID) * dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI);
  }
    return 0;
}

// used when only taking into account Fstar -GS
double FgtrConditionalM_LINELUMPOP2C(double z, double M1, double M2, double delta1, double delta2, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar, int LineID) {
    double result, error, lower_limit, upper_limit;
    gsl_function F;
    double rel_tol = 0.005; //<- relative tolerance
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
        .z_obs = z,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
		    .frac_star = Fstar10,
		    .LimitMass_Fstar = Mlim_Fstar,
        .lid = LineID,
    };

    F.function = &dFgtrConditionallnM_LINELUMPOP2C;
    F.params = &parameters_gsl_SFR_con;
    lower_limit = M1;
    upper_limit = M2;

    if(USE_GENERAL_SOURCES)
    {
      sources s;
      s = defaultSources();
      lower_limit = log(s.minMass(z));
    }

    gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                        1000, GSL_INTEG_GAUSS61, w, &result, &error);
    gsl_integration_workspace_free (w);

    if(delta2 >= delta1) {
        result = 1.;
        return result;
    }
    else {
        return result;
    }

}

// used below for LINELUMPOP2D -GS
double dFgtrConditionallnM_LINELUMPOP2D(double lnM, void *params) {
  struct parameters_gsl_SFR_con_int_ vals = *(struct parameters_gsl_SFR_con_int_ *)params;
  double M = exp(lnM); // linear scale
  double z = vals.z_obs;
  double M2 = vals.Mval; // natural log scale
  double del1 = vals.delta1;
  double del2 = vals.delta2;
  double MassTurnover = vals.Mdrop;
  double Alpha_star = vals.pl_star;
	double Fstar10 = vals.frac_star;
	double Mlim_Fstar = vals.LimitMass_Fstar;
  int LineID = vals.lid;
	double Fstar;

  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    return src.linelum(z, M, -2., LineID) * dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI);
  }
    return 0;
}

// used when only taking into account Fstar -GS
double FgtrConditionalM_LINELUMPOP2D(double z, double M1, double M2, double delta1, double delta2, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar, int LineID) {
    double result, error, lower_limit, upper_limit;
    gsl_function F;
    double rel_tol = 0.005; //<- relative tolerance
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
        .z_obs = z,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
		    .frac_star = Fstar10,
		    .LimitMass_Fstar = Mlim_Fstar,
        .lid = LineID,
    };

    F.function = &dFgtrConditionallnM_LINELUMPOP2D;
    F.params = &parameters_gsl_SFR_con;
    lower_limit = M1;
    upper_limit = M2;

    if(USE_GENERAL_SOURCES)
    {
      sources s;
      s = defaultSources();
      lower_limit = log(s.minMass(z));
    }

    gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                        1000, GSL_INTEG_GAUSS61, w, &result, &error);
    gsl_integration_workspace_free (w);

    if(delta2 >= delta1) {
        result = 1.;
        return result;
    }
    else {
        return result;
    }

}

// used below for LINELUMPOP2E -GS
double dFgtrConditionallnM_LINELUMPOP2E(double lnM, void *params) {
  struct parameters_gsl_SFR_con_int_ vals = *(struct parameters_gsl_SFR_con_int_ *)params;
  double M = exp(lnM); // linear scale
  double z = vals.z_obs;
  double M2 = vals.Mval; // natural log scale
  double del1 = vals.delta1;
  double del2 = vals.delta2;
  double MassTurnover = vals.Mdrop;
  double Alpha_star = vals.pl_star;
	double Fstar10 = vals.frac_star;
	double Mlim_Fstar = vals.LimitMass_Fstar;
  int LineID = vals.lid;
	double Fstar;

  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    return src.linelum(z, M, -2., LineID) * dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI);
  }
    return 0;
}

// used when only taking into account Fstar -GS
double FgtrConditionalM_LINELUMPOP2E(double z, double M1, double M2, double delta1, double delta2, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar, int LineID) {
    double result, error, lower_limit, upper_limit;
    gsl_function F;
    double rel_tol = 0.005; //<- relative tolerance
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
        .z_obs = z,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
		    .frac_star = Fstar10,
		    .LimitMass_Fstar = Mlim_Fstar,
        .lid = LineID,
    };

    F.function = &dFgtrConditionallnM_LINELUMPOP2E;
    F.params = &parameters_gsl_SFR_con;
    lower_limit = M1;
    upper_limit = M2;

    if(USE_GENERAL_SOURCES)
    {
      sources s;
      s = defaultSources();
      lower_limit = log(s.minMass(z));
    }

    gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                        1000, GSL_INTEG_GAUSS61, w, &result, &error);
    gsl_integration_workspace_free (w);

    if(delta2 >= delta1) {
        result = 1.;
        return result;
    }
    else {
        return result;
    }

}

// used below for LINELUMPOP2F -GS
double dFgtrConditionallnM_LINELUMPOP2F(double lnM, void *params) {
  struct parameters_gsl_SFR_con_int_ vals = *(struct parameters_gsl_SFR_con_int_ *)params;
  double M = exp(lnM); // linear scale
  double z = vals.z_obs;
  double M2 = vals.Mval; // natural log scale
  double del1 = vals.delta1;
  double del2 = vals.delta2;
  double MassTurnover = vals.Mdrop;
  double Alpha_star = vals.pl_star;
	double Fstar10 = vals.frac_star;
	double Mlim_Fstar = vals.LimitMass_Fstar;
  int LineID = vals.lid;
	double Fstar;

  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    return src.linelum(z, M, -2., LineID) * dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI);
  }
    return 0;
}

// used when only taking into account Fstar -GS
double FgtrConditionalM_LINELUMPOP2F(double z, double M1, double M2, double delta1, double delta2, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar, int LineID) {
    double result, error, lower_limit, upper_limit;
    gsl_function F;
    double rel_tol = 0.005; //<- relative tolerance
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
        .z_obs = z,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
		    .frac_star = Fstar10,
		    .LimitMass_Fstar = Mlim_Fstar,
        .lid = LineID,
    };

    F.function = &dFgtrConditionallnM_LINELUMPOP2F;
    F.params = &parameters_gsl_SFR_con;
    lower_limit = M1;
    upper_limit = M2;

    if(USE_GENERAL_SOURCES)
    {
      sources s;
      s = defaultSources();
      lower_limit = log(s.minMass(z));
    }

    gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                        1000, GSL_INTEG_GAUSS61, w, &result, &error);
    gsl_integration_workspace_free (w);

    if(delta2 >= delta1) {
        result = 1.;
        return result;
    }
    else {
        return result;
    }

}

// used below for LINELUMPOP2G -GS
double dFgtrConditionallnM_LINELUMPOP2G(double lnM, void *params) {
  struct parameters_gsl_SFR_con_int_ vals = *(struct parameters_gsl_SFR_con_int_ *)params;
  double M = exp(lnM); // linear scale
  double z = vals.z_obs;
  double M2 = vals.Mval; // natural log scale
  double del1 = vals.delta1;
  double del2 = vals.delta2;
  double MassTurnover = vals.Mdrop;
  double Alpha_star = vals.pl_star;
	double Fstar10 = vals.frac_star;
	double Mlim_Fstar = vals.LimitMass_Fstar;
  int LineID = vals.lid;
	double Fstar;

  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    return src.linelum(z, M, -2., LineID) * dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI);
  }
    return 0;
}

// used when only taking into account Fstar -GS
double FgtrConditionalM_LINELUMPOP2G(double z, double M1, double M2, double delta1, double delta2, double MassTurnover, double Alpha_star, double Fstar10, double Mlim_Fstar, int LineID) {
    double result, error, lower_limit, upper_limit;
    gsl_function F;
    double rel_tol = 0.005; //<- relative tolerance
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
        .z_obs = z,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
		    .frac_star = Fstar10,
		    .LimitMass_Fstar = Mlim_Fstar,
        .lid = LineID,
    };

    F.function = &dFgtrConditionallnM_LINELUMPOP2G;
    F.params = &parameters_gsl_SFR_con;
    lower_limit = M1;
    upper_limit = M2;

    if(USE_GENERAL_SOURCES)
    {
      sources s;
      s = defaultSources();
      lower_limit = log(s.minMass(z));
    }

    gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                        1000, GSL_INTEG_GAUSS61, w, &result, &error);
    gsl_integration_workspace_free (w);

    if(delta2 >= delta1) {
        result = 1.;
        return result;
    }
    else {
        return result;
    }

}

/*
// same as above, but with fx and fesc = 1 and no Nion
double dFgtrConditionallnM_SFR_Xray(double lnM, void *params) {
    struct parameters_gsl_SFR_con_int_ vals = *(struct parameters_gsl_SFR_con_int_ *)params;
    double M = exp(lnM); // linear scale
    double z = vals.z_obs;
    double M2 = vals.Mval; // natural log scale
    double del1 = vals.delta1;
    double del2 = vals.delta2;
    double MassTurnover = vals.Mdrop;
    double Alpha_star = vals.pl_star;
    double Alpha_esc = vals.pl_esc;
  double Fstar10 = vals.frac_star;
  double Fesc10 = vals.frac_esc;
  double Mlim_Fstar = vals.LimitMass_Fstar;
  double Mlim_Fesc = vals.LimitMass_Fesc;

  double Fstar,Fesc;

  if (Alpha_star > 0. && M > Mlim_Fstar)
    Fstar = 1./Fstar10;
  else if (Alpha_star < 0. && M < Mlim_Fstar)
    Fstar = 1./Fstar10;
  else
    Fstar = pow(M/1e10,Alpha_star);

  if (Alpha_esc > 0. && M > Mlim_Fesc)
    Fesc = 1./Fesc10;
  else if (Alpha_esc < 0. && M < Mlim_Fesc)
    Fesc = 1./Fesc10;
  else
    Fesc = pow(M/1e10,Alpha_esc);

  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    return M * dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI) * src.fstar(z, M) * src.fx(z, M);
  }

    return M*exp(-MassTurnover/M)*Fstar*Fesc*dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI);
}
*/

// same as above, but with fx and fesc = 1 and no Nion
double dFgtrConditionallnM_SFR_Xray(double lnM, void *params) {
    struct parameters_gsl_SFR_con_int_ vals = *(struct parameters_gsl_SFR_con_int_ *)params;
    double M = exp(lnM); // linear scale
    double z = vals.z_obs;
    double M2 = vals.Mval; // natural log scale
    double del1 = vals.delta1;
    double del2 = vals.delta2;
    double MassTurnover = vals.Mdrop;
    double Alpha_star = vals.pl_star;
    double Alpha_esc = vals.pl_esc;
    double Fstar10 = vals.frac_star;
    double Fesc10 = vals.frac_esc;
    double Mlim_Fstar = vals.LimitMass_Fstar;
    double Mlim_Fesc = vals.LimitMass_Fesc;
    double Fstar,Fesc;

    if(USE_GENERAL_SOURCES)
    {
      sources src;
      src = defaultSources();
      //return M * dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI) * src.fstar(z, M) * src.fx(z, M);
      return dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI) * (src.sfr(z, M)/f_baryon)/SperYR * src.fx(z, M);
      //return M * dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI) * src.fstar(z, M) * hubble(z) / 0.5 * src.fx(z, M);
    }

    return 0;
}


double FgtrConditionalM_SFR_Xray(double z, double M1, double M2, double delta1, double delta2, double MassTurnover, double Alpha_star, double Alpha_esc, double Fstar10, double Fesc10, double Mlim_Fstar, double Mlim_Fesc) {
    double result, error, lower_limit, upper_limit;
    gsl_function F;
    double rel_tol = 0.005; //<- relative tolerance
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
        .z_obs = z,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
        .pl_esc = Alpha_esc,
    .frac_star = Fstar10,
    .frac_esc = Fesc10,
    .LimitMass_Fstar = Mlim_Fstar,
    .LimitMass_Fesc = Mlim_Fesc
    };

    F.function = &dFgtrConditionallnM_SFR_Xray;
    F.params = &parameters_gsl_SFR_con;
    lower_limit = M1;
    upper_limit = M2;

    if(USE_GENERAL_SOURCES)
    {
      sources s;
      s = defaultSources();
      lower_limit = log(s.minMass(z));
    }

    gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                        1000, GSL_INTEG_GAUSS61, w, &result, &error);
    gsl_integration_workspace_free (w);

    if(delta2 >= delta1) {
        result = 1.;
        return result;
    }
    else {
        return result;
    }

}

/*
// same as above, but with fx and fesc = 1 and no Nion and no fx
double dFgtrConditionallnM_SFR_Lya(double lnM, void *params) {
    struct parameters_gsl_SFR_con_int_ vals = *(struct parameters_gsl_SFR_con_int_ *)params;
    double M = exp(lnM); // linear scale
    double z = vals.z_obs;
    double M2 = vals.Mval; // natural log scale
    double del1 = vals.delta1;
    double del2 = vals.delta2;
    double MassTurnover = vals.Mdrop;
    double Alpha_star = vals.pl_star;
    double Alpha_esc = vals.pl_esc;
  double Fstar10 = vals.frac_star;
  double Fesc10 = vals.frac_esc;
  double Mlim_Fstar = vals.LimitMass_Fstar;
  double Mlim_Fesc = vals.LimitMass_Fesc;

  double Fstar,Fesc;

  if (Alpha_star > 0. && M > Mlim_Fstar)
    Fstar = 1./Fstar10;
  else if (Alpha_star < 0. && M < Mlim_Fstar)
    Fstar = 1./Fstar10;
  else
    Fstar = pow(M/1e10,Alpha_star);

  if (Alpha_esc > 0. && M > Mlim_Fesc)
    Fesc = 1./Fesc10;
  else if (Alpha_esc < 0. && M < Mlim_Fesc)
    Fesc = 1./Fesc10;
  else
    Fesc = pow(M/1e10,Alpha_esc);

  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    return M * dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI) * src.fstar(z, M);
  }

    return M*exp(-MassTurnover/M)*Fstar*Fesc*dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI);
}
*/

// same as above, but with fx and fesc = 1 and no Nion and no fx
double dFgtrConditionallnM_SFR_Lya(double lnM, void *params) {
    struct parameters_gsl_SFR_con_int_ vals = *(struct parameters_gsl_SFR_con_int_ *)params;
    double M = exp(lnM); // linear scale
    double z = vals.z_obs;
    double M2 = vals.Mval; // natural log scale
    double del1 = vals.delta1;
    double del2 = vals.delta2;
    double MassTurnover = vals.Mdrop;
    double Alpha_star = vals.pl_star;
    double Alpha_esc = vals.pl_esc;
    double Fstar10 = vals.frac_star;
    double Fesc10 = vals.frac_esc;
    double Mlim_Fstar = vals.LimitMass_Fstar;
    double Mlim_Fesc = vals.LimitMass_Fesc;
    double Fstar,Fesc;

    if(USE_GENERAL_SOURCES)
    {
      sources src;
      src = defaultSources();
      //return M * dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI) * src.fstar(z, M);
      return dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI) * (src.sfr(z, M)/f_baryon)/SperYR;
      //return M * dNdM_conditional_second(z,log(M),M2,del1,del2)/sqrt(2.*PI) * src.fstar(z, M) * hubble(z) / 0.5;
    }

    return 0;
}


double FgtrConditionalM_SFR_Lya(double z, double M1, double M2, double delta1, double delta2, double MassTurnover, double Alpha_star, double Alpha_esc, double Fstar10, double Fesc10, double Mlim_Fstar, double Mlim_Fesc) {
    double result, error, lower_limit, upper_limit;
    gsl_function F;
    double rel_tol = 0.005; //<- relative tolerance
    gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);

    struct parameters_gsl_SFR_con_int_ parameters_gsl_SFR_con = {
        .z_obs = z,
        .Mval = M2,
        .delta1 = delta1,
        .delta2 = delta2,
        .Mdrop = MassTurnover,
        .pl_star = Alpha_star,
        .pl_esc = Alpha_esc,
    .frac_star = Fstar10,
    .frac_esc = Fesc10,
    .LimitMass_Fstar = Mlim_Fstar,
    .LimitMass_Fesc = Mlim_Fesc
    };

    F.function = &dFgtrConditionallnM_SFR_Lya;
    F.params = &parameters_gsl_SFR_con;
    lower_limit = M1;
    upper_limit = M2;

    if(USE_GENERAL_SOURCES)
    {
      sources s;
      s = defaultSources();
      lower_limit = log(s.minMass(z));
    }

    gsl_integration_qag (&F, lower_limit, upper_limit, 0, rel_tol,
                        1000, GSL_INTEG_GAUSS61, w, &result, &error);
    gsl_integration_workspace_free (w);

    if(delta2 >= delta1) {
        result = 1.;
        return result;
    }
    else {
        return result;
    }

}

// Set up interpolation table for collapse fraction and initialise interpolation.
// At a redshift compute conditional mass function over the density field.
// Split density range into low and high density, since to increase accuray in high density region.
void initialiseFcollSFR_spline(float z, float Mmax, float MassTurnover, float Alpha_star, float Alpha_esc, float Fstar10, float Fesc10, float Mlim_Fstar, float Mlim_Fesc){
	float Mmin = MassTurnover/50.;
  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    Mmin = src.minMass(z);
  }
	//float Mmin = 1e8; // TEST
    double overdense_val;
    double overdense_large_high = Deltac, overdense_large_low = 1.5;
    double overdense_small_high = 1.5, overdense_small_low = -1. + 9e-8;
    int i;
    FcollLow_spline_acc = gsl_interp_accel_alloc ();
    FcollLow_spline = gsl_spline_alloc (gsl_interp_cspline, NSFR_low);

    for (i=0; i<NSFR_low; i++){
        overdense_val = log10(1. + overdense_small_low) + (double)i/((double)NSFR_low-1.)*(log10(1.+overdense_small_high)-log10(1.+overdense_small_low));

        log10_overdense_spline_SFR[i] = overdense_val;
        log10_Fcoll_spline_SFR[i] = log10(GaussLegendreQuad_FcollSFR(NGL_SFR,z,log(Mmax),Deltac,pow(10.,overdense_val)-1.,MassTurnover,Alpha_star,Alpha_esc,Fstar10,Fesc10,Mlim_Fstar,Mlim_Fesc));

        if(log10_Fcoll_spline_SFR[i] < -40.){
            log10_Fcoll_spline_SFR[i] = -40.;
        }
    }
    gsl_spline_init(FcollLow_spline, log10_overdense_spline_SFR, log10_Fcoll_spline_SFR, NSFR_low);


    for(i=0;i<NSFR_high;i++) {
        Overdense_spline_SFR[i] = overdense_large_low + (float)i/((float)NSFR_high-1.)*(overdense_large_high - overdense_large_low);
        Fcoll_spline_SFR[i] = FgtrConditionalM_SFR(z,log(Mmin),log(Mmax),Deltac,Overdense_spline_SFR[i],MassTurnover,Alpha_star,Alpha_esc,Fstar10,Fesc10,Mlim_Fstar,Mlim_Fesc);

        if(Fcoll_spline_SFR[i]<0.) {
            Fcoll_spline_SFR[i]=pow(10.,-40.0);
        }
    }
    spline(Overdense_spline_SFR-1,Fcoll_spline_SFR-1,NSFR_high,0,0,second_derivs_SFR-1);
}

// Find collapse fraction given overdensity using interpolation.
void FcollSpline_SFR(float Overdensity, float *splined_value){
    int i;
    float returned_value;

    if(Overdensity<1.5) {
        if(Overdensity<-1.) {
            returned_value = 0;
        }
        else {
            returned_value = gsl_spline_eval(FcollLow_spline, log10(Overdensity+1.), FcollLow_spline_acc);
            returned_value = pow(10.,returned_value);
        }
    }
    else {
        if(Overdensity<0.99*Deltac) {
		// Usage of 0.99*Deltac arises due to the fact that close to the critical density, the collapsed fraction becomes a little unstable
		// However, such densities should always be collapsed, so just set f_coll to unity.
		// Additionally, the fraction of points in this regime relative to the entire simulation volume is extremely small.
            splint(Overdense_spline_SFR-1,Fcoll_spline_SFR-1,second_derivs_SFR-1,NSFR_high,Overdensity,&(returned_value));
        }
        else {
            returned_value = 1.;
        }
    }
	if(returned_value > 1.) returned_value = 1.;
    *splined_value = returned_value;
}

void initialiseFcollSFRPOP2_spline(float z, float Mmax, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar){
	float Mmin = MassTurnover/50.;
  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    Mmin = src.minMass(z);
  }
	//float Mmin = 1e8; // TEST
    double overdense_val;
    double overdense_large_high = Deltac, overdense_large_low = 1.5;
    double overdense_small_high = 1.5, overdense_small_low = -1. + 9e-8;
    int i;
    FcollLow_spline_acc_SFRPOP2 = gsl_interp_accel_alloc ();
    FcollLow_spline_SFRPOP2 = gsl_spline_alloc (gsl_interp_cspline, NSFR_low);
    for (i=0; i<NSFR_low; i++){
        overdense_val = log10(1. + overdense_small_low) + (double)i/((double)NSFR_low-1.)*(log10(1.+overdense_small_high)-log10(1.+overdense_small_low));

        log10_overdense_spline_SFRPOP2[i] = overdense_val;
        log10_Fcoll_spline_SFRPOP2[i] = log10(GaussLegendreQuad_FcollSFRPOP2(NGL_SFR, z, log(Mmax), Deltac, pow(10.,overdense_val)-1., MassTurnover, Alpha_star, Fstar10, Mlim_Fstar));

        if(log10_Fcoll_spline_SFRPOP2[i] < -40.){
            log10_Fcoll_spline_SFRPOP2[i] = -40.;
        }
    }
    gsl_spline_init(FcollLow_spline_SFRPOP2, log10_overdense_spline_SFRPOP2, log10_Fcoll_spline_SFRPOP2, NSFR_low);

    for(i=0;i<NSFR_high;i++) {
        Overdense_spline_SFRPOP2[i] = overdense_large_low + (float)i/((float)NSFR_high-1.)*(overdense_large_high - overdense_large_low);
        Fcoll_spline_SFRPOP2[i] = FgtrConditionalM_SFRPOP2(z, log(Mmin), log(Mmax), Deltac, Overdense_spline_SFRPOP2[i], MassTurnover, Alpha_star, Fstar10, Mlim_Fstar);
        if(Fcoll_spline_SFRPOP2[i]<0.) {
            Fcoll_spline_SFRPOP2[i]=pow(10.,-40.0);
        }
    }
    spline(Overdense_spline_SFRPOP2-1, Fcoll_spline_SFRPOP2-1, NSFR_high, 0, 0, second_derivs_SFRPOP2-1);
}

void initialiseFcollgas_spline(float z, float Mmax, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar){
  float Mmin = MassTurnover/50.;
  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    Mmin = src.minMass(z);
  }
  //float Mmin = 1e8; // TEST
    double overdense_val;
    double overdense_large_high = Deltac, overdense_large_low = 1.5;
    double overdense_small_high = 1.5, overdense_small_low = -1. + 9e-8;
    int i;
    FcollLow_spline_acc_gas = gsl_interp_accel_alloc ();
    FcollLow_spline_gas = gsl_spline_alloc (gsl_interp_cspline, NSFR_low);
    for (i=0; i<NSFR_low; i++){
        overdense_val = log10(1. + overdense_small_low) + (double)i/((double)NSFR_low-1.)*(log10(1.+overdense_small_high)-log10(1.+overdense_small_low));

        log10_overdense_spline_gas[i] = overdense_val;
        log10_Fcoll_spline_gas[i] = log10(GaussLegendreQuad_Fcollgas(NGL_SFR, z, log(Mmax), Deltac, pow(10.,overdense_val)-1., MassTurnover, Alpha_star, Fstar10, Mlim_Fstar));

        if(log10_Fcoll_spline_gas[i] < -40.){
            log10_Fcoll_spline_gas[i] = -40.;
        }
    }
    gsl_spline_init(FcollLow_spline_gas, log10_overdense_spline_gas, log10_Fcoll_spline_gas, NSFR_low);

    for(i=0;i<NSFR_high;i++) {
        Overdense_spline_gas[i] = overdense_large_low + (float)i/((float)NSFR_high-1.)*(overdense_large_high - overdense_large_low);
        Fcoll_spline_gas[i] = FgtrConditionalM_gas(z, log(Mmin), log(Mmax), Deltac, Overdense_spline_gas[i], MassTurnover, Alpha_star, Fstar10, Mlim_Fstar);
        if(Fcoll_spline_gas[i]<0.) {
            Fcoll_spline_gas[i]=pow(10.,-40.0);
        }
    }
    spline(Overdense_spline_gas-1, Fcoll_spline_gas-1, NSFR_high, 0, 0, second_derivs_gas-1);
}

void initialiseFcollmet_spline(float z, float Mmax, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar){
  float Mmin = MassTurnover/50.;
  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    Mmin = src.minMass(z);
  }
  //float Mmin = 1e8; // TEST
    double overdense_val;
    double overdense_large_high = Deltac, overdense_large_low = 1.5;
    double overdense_small_high = 1.5, overdense_small_low = -1. + 9e-8;
    int i;
    FcollLow_spline_acc_met = gsl_interp_accel_alloc ();
    FcollLow_spline_met = gsl_spline_alloc (gsl_interp_cspline, NSFR_low);
    for (i=0; i<NSFR_low; i++){
        overdense_val = log10(1. + overdense_small_low) + (double)i/((double)NSFR_low-1.)*(log10(1.+overdense_small_high)-log10(1.+overdense_small_low));

        log10_overdense_spline_met[i] = overdense_val;
        log10_Fcoll_spline_met[i] = log10(GaussLegendreQuad_Fcollmet(NGL_SFR, z, log(Mmax), Deltac, pow(10.,overdense_val)-1., MassTurnover, Alpha_star, Fstar10, Mlim_Fstar));

        if(log10_Fcoll_spline_met[i] < -40.){
            log10_Fcoll_spline_met[i] = -40.;
        }
    }
    gsl_spline_init(FcollLow_spline_met, log10_overdense_spline_met, log10_Fcoll_spline_met, NSFR_low);

    for(i=0;i<NSFR_high;i++) {
        Overdense_spline_met[i] = overdense_large_low + (float)i/((float)NSFR_high-1.)*(overdense_large_high - overdense_large_low);
        Fcoll_spline_met[i] = FgtrConditionalM_met(z, log(Mmin), log(Mmax), Deltac, Overdense_spline_met[i], MassTurnover, Alpha_star, Fstar10, Mlim_Fstar);
        if(Fcoll_spline_met[i]<0.) {
            Fcoll_spline_met[i]=pow(10.,-40.0);
        }
    }
    spline(Overdense_spline_met-1, Fcoll_spline_met-1, NSFR_high, 0, 0, second_derivs_met-1);
}

void FcollSpline_SFRPOP2(float Overdensity, float *splined_value){
    int i;
    float returned_value;

    if(Overdensity<1.5) {
        if(Overdensity<-1.) {
            returned_value = 0;
        }
        else {
            returned_value = gsl_spline_eval(FcollLow_spline_SFRPOP2, log10(Overdensity+1.), FcollLow_spline_acc_SFRPOP2);
            returned_value = pow(10.,returned_value);
        }
    }
    else {
        if(Overdensity<0.99*Deltac) {
		// Usage of 0.99*Deltac arises due to the fact that close to the critical density, the collapsed fraction becomes a little unstable
		// However, such densities should always be collapsed, so just set f_coll to unity.
		// Additionally, the fraction of points in this regime relative to the entire simulation volume is extremely small.
            splint(Overdense_spline_SFRPOP2-1, Fcoll_spline_SFRPOP2-1, second_derivs_SFRPOP2-1, NSFR_high, Overdensity, &(returned_value));
        }
        else {
            // We set overdensity of regions close to critical density to critical density to avoid instability.
            //returned_value = 1.;
            splint(Overdense_spline_SFRPOP2-1, Fcoll_spline_SFRPOP2-1, second_derivs_SFRPOP2-1, NSFR_high, Overdensity, &(returned_value));
        }
    }
    // We don't need the following ceiling condition any more... -GS
	  //if(returned_value > 1.) returned_value = 1.;
    *splined_value = returned_value;
}

void FcollSpline_gas(float Overdensity, float *splined_value){
    int i;
    float returned_value;

    if(Overdensity<0) {
        if(Overdensity<-1.) {
            returned_value = 0;
        }
        else {
            returned_value = gsl_spline_eval(FcollLow_spline_gas, log10(Overdensity+1.), FcollLow_spline_acc_gas);
            returned_value = pow(10.,returned_value);
        }
    }
    else {
        if(Overdensity<0.99*Deltac) {
    // Usage of 0.99*Deltac arises due to the fact that close to the critical density, the collapsed fraction becomes a little unstable
    // However, such densities should always be collapsed, so just set f_coll to unity.
    // Additionally, the fraction of points in this regime relative to the entire simulation volume is extremely small.
            splint(Overdense_spline_gas-1, Fcoll_spline_gas-1, second_derivs_gas-1, NSFR_high, Overdensity, &(returned_value));
        }
        else {
            // We set overdensity of regions close to critical density to critical density to avoid instability.
            //returned_value = 1.;
            splint(Overdense_spline_gas-1, Fcoll_spline_gas-1, second_derivs_gas-1, NSFR_high, Overdensity, &(returned_value));
        }
    }
    // We don't need the following ceiling condition any more... -GS
    //if(returned_value > 1.) returned_value = 1.;
    *splined_value = returned_value;
}

void FcollSpline_met(float Overdensity, float *splined_value){
    int i;
    float returned_value;

    if(Overdensity<0) {
        if(Overdensity<-1.) {
            returned_value = 0;
        }
        else {
            returned_value = gsl_spline_eval(FcollLow_spline_met, log10(Overdensity+1.), FcollLow_spline_acc_met);
            returned_value = pow(10.,returned_value);
        }
    }
    else {
        if(Overdensity<0.99*Deltac) {
    // Usage of 0.99*Deltac arises due to the fact that close to the critical density, the collapsed fraction becomes a little unstable
    // However, such densities should always be collapsed, so just set f_coll to unity.
    // Additionally, the fraction of points in this regime relative to the entire simulation volume is extremely small.
            splint(Overdense_spline_met-1, Fcoll_spline_met-1, second_derivs_met-1, NSFR_high, Overdensity, &(returned_value));
        }
        else {
            // We set overdensity of regions close to critical density to critical density to avoid instability.
            //returned_value = 1.;
            splint(Overdense_spline_met-1, Fcoll_spline_met-1, second_derivs_met-1, NSFR_high, Overdensity, &(returned_value));
        }
    }
    // We don't need the following ceiling condition any more... -GS
    //if(returned_value > 1.) returned_value = 1.;
    *splined_value = returned_value;
}

// Interpolation for LINELUMPOP2A -GS
void initialiseFcollLINELUMPOP2A_spline(float z, float Mmax, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar, int LineID){
	float Mmin = MassTurnover/50.;
  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    Mmin = src.minMass(z);
  }
	//float Mmin = 1e8; // TEST
    double overdense_val;
    double overdense_large_high = Deltac, overdense_large_low = 1.5;
    double overdense_small_high = 1.5, overdense_small_low = -1. + 9e-8;
    int i;
    FcollLow_spline_acc_LINELUMPOP2A = gsl_interp_accel_alloc ();
    FcollLow_spline_LINELUMPOP2A = gsl_spline_alloc (gsl_interp_cspline, NSFR_low);
    for (i=0; i<NSFR_low; i++){
        overdense_val = log10(1. + overdense_small_low) + (double)i/((double)NSFR_low-1.)*(log10(1.+overdense_small_high)-log10(1.+overdense_small_low));

        log10_overdense_spline_LINELUMPOP2A[i] = overdense_val;
        log10_Fcoll_spline_LINELUMPOP2A[i] = log10(GaussLegendreQuad_FcollLINELUMPOP2A(NGL_SFR, z, log(Mmax), Deltac, pow(10.,overdense_val)-1., MassTurnover, Alpha_star, Fstar10, Mlim_Fstar, LineID));

        if(log10_Fcoll_spline_LINELUMPOP2A[i] < -40.){
            log10_Fcoll_spline_LINELUMPOP2A[i] = -40.;
        }
    }
    gsl_spline_init(FcollLow_spline_LINELUMPOP2A, log10_overdense_spline_LINELUMPOP2A, log10_Fcoll_spline_LINELUMPOP2A, NSFR_low);

    for(i=0;i<NSFR_high;i++) {
        Overdense_spline_LINELUMPOP2A[i] = overdense_large_low + (float)i/((float)NSFR_high-1.)*(overdense_large_high - overdense_large_low);
        Fcoll_spline_LINELUMPOP2A[i] = FgtrConditionalM_LINELUMPOP2A(z, log(Mmin), log(Mmax), Deltac, Overdense_spline_LINELUMPOP2A[i], MassTurnover, Alpha_star, Fstar10, Mlim_Fstar, LineID);
        if(Fcoll_spline_LINELUMPOP2A[i]<0.) {
            Fcoll_spline_LINELUMPOP2A[i]=pow(10.,-40.0);
        }
    }
    spline(Overdense_spline_LINELUMPOP2A-1, Fcoll_spline_LINELUMPOP2A-1, NSFR_high, 0, 0, second_derivs_LINELUMPOP2A-1);
}

void FcollSpline_LINELUMPOP2A(float Overdensity, float *splined_value){
    int i;
    float returned_value;

    if(Overdensity<1.5) {
        if(Overdensity<-1.) {
            returned_value = 0;
        }
        else {
            returned_value = gsl_spline_eval(FcollLow_spline_LINELUMPOP2A, log10(Overdensity+1.), FcollLow_spline_acc_LINELUMPOP2A);
            returned_value = pow(10.,returned_value);
        }
    }
    else {
        if(Overdensity<0.99*Deltac) {
		// Usage of 0.99*Deltac arises due to the fact that close to the critical density, the collapsed fraction becomes a little unstable
		// However, such densities should always be collapsed, so just set f_coll to unity.
		// Additionally, the fraction of points in this regime relative to the entire simulation volume is extremely small.
            splint(Overdense_spline_LINELUMPOP2A-1, Fcoll_spline_LINELUMPOP2A-1, second_derivs_LINELUMPOP2A-1, NSFR_high, Overdensity, &(returned_value));
        }
        else {
            // We set overdensity of regions close to critical density to critical density to avoid instability.
            //returned_value = 1.;
            splint(Overdense_spline_LINELUMPOP2A-1, Fcoll_spline_LINELUMPOP2A-1, second_derivs_LINELUMPOP2A-1, NSFR_high, Overdensity, &(returned_value));
        }
    }
    // We don't need the following ceiling condition any more... -GS
	  //if(returned_value > 1.) returned_value = 1.;
    *splined_value = returned_value;
}

// Interpolation for LINELUMPOP2B -GS
void initialiseFcollLINELUMPOP2B_spline(float z, float Mmax, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar, int LineID){
	float Mmin = MassTurnover/50.;
  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    Mmin = src.minMass(z);
  }
	//float Mmin = 1e8; // TEST
    double overdense_val;
    double overdense_large_high = Deltac, overdense_large_low = 1.5;
    double overdense_small_high = 1.5, overdense_small_low = -1. + 9e-8;
    int i;
    FcollLow_spline_acc_LINELUMPOP2B = gsl_interp_accel_alloc ();
    FcollLow_spline_LINELUMPOP2B = gsl_spline_alloc (gsl_interp_cspline, NSFR_low);
    for (i=0; i<NSFR_low; i++){
        overdense_val = log10(1. + overdense_small_low) + (double)i/((double)NSFR_low-1.)*(log10(1.+overdense_small_high)-log10(1.+overdense_small_low));

        log10_overdense_spline_LINELUMPOP2B[i] = overdense_val;
        log10_Fcoll_spline_LINELUMPOP2B[i] = log10(GaussLegendreQuad_FcollLINELUMPOP2B(NGL_SFR, z, log(Mmax), Deltac, pow(10.,overdense_val)-1., MassTurnover, Alpha_star, Fstar10, Mlim_Fstar, LineID));

        if(log10_Fcoll_spline_LINELUMPOP2B[i] < -40.){
            log10_Fcoll_spline_LINELUMPOP2B[i] = -40.;
        }
    }
    gsl_spline_init(FcollLow_spline_LINELUMPOP2B, log10_overdense_spline_LINELUMPOP2B, log10_Fcoll_spline_LINELUMPOP2B, NSFR_low);

    for(i=0;i<NSFR_high;i++) {
        Overdense_spline_LINELUMPOP2B[i] = overdense_large_low + (float)i/((float)NSFR_high-1.)*(overdense_large_high - overdense_large_low);
        Fcoll_spline_LINELUMPOP2B[i] = FgtrConditionalM_LINELUMPOP2B(z, log(Mmin), log(Mmax), Deltac, Overdense_spline_LINELUMPOP2B[i], MassTurnover, Alpha_star, Fstar10, Mlim_Fstar, LineID);
        if(Fcoll_spline_LINELUMPOP2B[i]<0.) {
            Fcoll_spline_LINELUMPOP2B[i]=pow(10.,-40.0);
        }
    }
    spline(Overdense_spline_LINELUMPOP2B-1, Fcoll_spline_LINELUMPOP2B-1, NSFR_high, 0, 0, second_derivs_LINELUMPOP2B-1);
}

void FcollSpline_LINELUMPOP2B(float Overdensity, float *splined_value){
    int i;
    float returned_value;

    if(Overdensity<1.5) {
        if(Overdensity<-1.) {
            returned_value = 0;
        }
        else {
            returned_value = gsl_spline_eval(FcollLow_spline_LINELUMPOP2B, log10(Overdensity+1.), FcollLow_spline_acc_LINELUMPOP2B);
            returned_value = pow(10.,returned_value);
        }
    }
    else {
        if(Overdensity<0.99*Deltac) {
		// Usage of 0.99*Deltac arises due to the fact that close to the critical density, the collapsed fraction becomes a little unstable
		// However, such densities should always be collapsed, so just set f_coll to unity.
		// Additionally, the fraction of points in this regime relative to the entire simulation volume is extremely small.
            splint(Overdense_spline_LINELUMPOP2B-1, Fcoll_spline_LINELUMPOP2B-1, second_derivs_LINELUMPOP2B-1, NSFR_high, Overdensity, &(returned_value));
        }
        else {
            // We set overdensity of regions close to critical density to critical density to avoid instability.
            //returned_value = 1.;
            splint(Overdense_spline_LINELUMPOP2B-1, Fcoll_spline_LINELUMPOP2B-1, second_derivs_LINELUMPOP2B-1, NSFR_high, Overdensity, &(returned_value));
        }
    }
    // We don't need the following ceiling condition any more... -GS
	  //if(returned_value > 1.) returned_value = 1.;
    *splined_value = returned_value;
}

// Interpolation for LINELUMPOP2C -GS
void initialiseFcollLINELUMPOP2C_spline(float z, float Mmax, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar, int LineID){
	float Mmin = MassTurnover/50.;
  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    Mmin = src.minMass(z);
  }
	//float Mmin = 1e8; // TEST
    double overdense_val;
    double overdense_large_high = Deltac, overdense_large_low = 1.5;
    double overdense_small_high = 1.5, overdense_small_low = -1. + 9e-8;
    int i;
    FcollLow_spline_acc_LINELUMPOP2C = gsl_interp_accel_alloc ();
    FcollLow_spline_LINELUMPOP2C = gsl_spline_alloc (gsl_interp_cspline, NSFR_low);
    for (i=0; i<NSFR_low; i++){
        overdense_val = log10(1. + overdense_small_low) + (double)i/((double)NSFR_low-1.)*(log10(1.+overdense_small_high)-log10(1.+overdense_small_low));

        log10_overdense_spline_LINELUMPOP2C[i] = overdense_val;
        log10_Fcoll_spline_LINELUMPOP2C[i] = log10(GaussLegendreQuad_FcollLINELUMPOP2C(NGL_SFR, z, log(Mmax), Deltac, pow(10.,overdense_val)-1., MassTurnover, Alpha_star, Fstar10, Mlim_Fstar, LineID));

        if(log10_Fcoll_spline_LINELUMPOP2C[i] < -40.){
            log10_Fcoll_spline_LINELUMPOP2C[i] = -40.;
        }
    }
    gsl_spline_init(FcollLow_spline_LINELUMPOP2C, log10_overdense_spline_LINELUMPOP2C, log10_Fcoll_spline_LINELUMPOP2C, NSFR_low);

    for(i=0;i<NSFR_high;i++) {
        Overdense_spline_LINELUMPOP2C[i] = overdense_large_low + (float)i/((float)NSFR_high-1.)*(overdense_large_high - overdense_large_low);
        Fcoll_spline_LINELUMPOP2C[i] = FgtrConditionalM_LINELUMPOP2C(z, log(Mmin), log(Mmax), Deltac, Overdense_spline_LINELUMPOP2C[i], MassTurnover, Alpha_star, Fstar10, Mlim_Fstar, LineID);
        if(Fcoll_spline_LINELUMPOP2C[i]<0.) {
            Fcoll_spline_LINELUMPOP2C[i]=pow(10.,-40.0);
        }
    }
    spline(Overdense_spline_LINELUMPOP2C-1, Fcoll_spline_LINELUMPOP2C-1, NSFR_high, 0, 0, second_derivs_LINELUMPOP2C-1);
}

void FcollSpline_LINELUMPOP2C(float Overdensity, float *splined_value){
    int i;
    float returned_value;

    if(Overdensity<1.5) {
        if(Overdensity<-1.) {
            returned_value = 0;
        }
        else {
            returned_value = gsl_spline_eval(FcollLow_spline_LINELUMPOP2C, log10(Overdensity+1.), FcollLow_spline_acc_LINELUMPOP2C);
            returned_value = pow(10.,returned_value);
        }
    }
    else {
        if(Overdensity<0.99*Deltac) {
		// Usage of 0.99*Deltac arises due to the fact that close to the critical density, the collapsed fraction becomes a little unstable
		// However, such densities should always be collapsed, so just set f_coll to unity.
		// Additionally, the fraction of points in this regime relative to the entire simulation volume is extremely small.
            splint(Overdense_spline_LINELUMPOP2C-1, Fcoll_spline_LINELUMPOP2C-1, second_derivs_LINELUMPOP2C-1, NSFR_high, Overdensity, &(returned_value));
        }
        else {
            // We set overdensity of regions close to critical density to critical density to avoid instability.
            //returned_value = 1.;
            splint(Overdense_spline_LINELUMPOP2C-1, Fcoll_spline_LINELUMPOP2C-1, second_derivs_LINELUMPOP2C-1, NSFR_high, Overdensity, &(returned_value));
        }
    }
    // We don't need the following ceiling condition any more... -GS
	  //if(returned_value > 1.) returned_value = 1.;
    *splined_value = returned_value;
}

// Interpolation for LINELUMPOP2D -GS
void initialiseFcollLINELUMPOP2D_spline(float z, float Mmax, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar, int LineID){
	float Mmin = MassTurnover/50.;
  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    Mmin = src.minMass(z);
  }
	//float Mmin = 1e8; // TEST
    double overdense_val;
    double overdense_large_high = Deltac, overdense_large_low = 1.5;
    double overdense_small_high = 1.5, overdense_small_low = -1. + 9e-8;
    int i;
    FcollLow_spline_acc_LINELUMPOP2D = gsl_interp_accel_alloc ();
    FcollLow_spline_LINELUMPOP2D = gsl_spline_alloc (gsl_interp_cspline, NSFR_low);
    for (i=0; i<NSFR_low; i++){
        overdense_val = log10(1. + overdense_small_low) + (double)i/((double)NSFR_low-1.)*(log10(1.+overdense_small_high)-log10(1.+overdense_small_low));

        log10_overdense_spline_LINELUMPOP2D[i] = overdense_val;
        log10_Fcoll_spline_LINELUMPOP2D[i] = log10(GaussLegendreQuad_FcollLINELUMPOP2D(NGL_SFR, z, log(Mmax), Deltac, pow(10.,overdense_val)-1., MassTurnover, Alpha_star, Fstar10, Mlim_Fstar, LineID));

        if(log10_Fcoll_spline_LINELUMPOP2D[i] < -40.){
            log10_Fcoll_spline_LINELUMPOP2D[i] = -40.;
        }
    }
    gsl_spline_init(FcollLow_spline_LINELUMPOP2D, log10_overdense_spline_LINELUMPOP2D, log10_Fcoll_spline_LINELUMPOP2D, NSFR_low);

    for(i=0;i<NSFR_high;i++) {
        Overdense_spline_LINELUMPOP2D[i] = overdense_large_low + (float)i/((float)NSFR_high-1.)*(overdense_large_high - overdense_large_low);
        Fcoll_spline_LINELUMPOP2D[i] = FgtrConditionalM_LINELUMPOP2D(z, log(Mmin), log(Mmax), Deltac, Overdense_spline_LINELUMPOP2D[i], MassTurnover, Alpha_star, Fstar10, Mlim_Fstar, LineID);
        if(Fcoll_spline_LINELUMPOP2D[i]<0.) {
            Fcoll_spline_LINELUMPOP2D[i]=pow(10.,-40.0);
        }
    }
    spline(Overdense_spline_LINELUMPOP2D-1, Fcoll_spline_LINELUMPOP2D-1, NSFR_high, 0, 0, second_derivs_LINELUMPOP2D-1);
}

void FcollSpline_LINELUMPOP2D(float Overdensity, float *splined_value){
    int i;
    float returned_value;

    if(Overdensity<1.5) {
        if(Overdensity<-1.) {
            returned_value = 0;
        }
        else {
            returned_value = gsl_spline_eval(FcollLow_spline_LINELUMPOP2D, log10(Overdensity+1.), FcollLow_spline_acc_LINELUMPOP2D);
            returned_value = pow(10.,returned_value);
        }
    }
    else {
        if(Overdensity<0.99*Deltac) {
		// Usage of 0.99*Deltac arises due to the fact that close to the critical density, the collapsed fraction becomes a little unstable
		// However, such densities should always be collapsed, so just set f_coll to unity.
		// Additionally, the fraction of points in this regime relative to the entire simulation volume is extremely small.
            splint(Overdense_spline_LINELUMPOP2D-1, Fcoll_spline_LINELUMPOP2D-1, second_derivs_LINELUMPOP2D-1, NSFR_high, Overdensity, &(returned_value));
        }
        else {
            // We set overdensity of regions close to critical density to critical density to avoid instability.
            //returned_value = 1.;
            splint(Overdense_spline_LINELUMPOP2D-1, Fcoll_spline_LINELUMPOP2D-1, second_derivs_LINELUMPOP2D-1, NSFR_high, Overdensity, &(returned_value));
        }
    }
    // We don't need the following ceiling condition any more... -GS
	  //if(returned_value > 1.) returned_value = 1.;
    *splined_value = returned_value;
}

// Interpolation for LINELUMPOP2E -GS
void initialiseFcollLINELUMPOP2E_spline(float z, float Mmax, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar, int LineID){
	float Mmin = MassTurnover/50.;
  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    Mmin = src.minMass(z);
  }
	//float Mmin = 1e8; // TEST
    double overdense_val;
    double overdense_large_high = Deltac, overdense_large_low = 1.5;
    double overdense_small_high = 1.5, overdense_small_low = -1. + 9e-8;
    int i;
    FcollLow_spline_acc_LINELUMPOP2E = gsl_interp_accel_alloc ();
    FcollLow_spline_LINELUMPOP2E = gsl_spline_alloc (gsl_interp_cspline, NSFR_low);
    for (i=0; i<NSFR_low; i++){
        overdense_val = log10(1. + overdense_small_low) + (double)i/((double)NSFR_low-1.)*(log10(1.+overdense_small_high)-log10(1.+overdense_small_low));

        log10_overdense_spline_LINELUMPOP2E[i] = overdense_val;
        log10_Fcoll_spline_LINELUMPOP2E[i] = log10(GaussLegendreQuad_FcollLINELUMPOP2E(NGL_SFR, z, log(Mmax), Deltac, pow(10.,overdense_val)-1., MassTurnover, Alpha_star, Fstar10, Mlim_Fstar, LineID));

        if(log10_Fcoll_spline_LINELUMPOP2E[i] < -40.){
            log10_Fcoll_spline_LINELUMPOP2E[i] = -40.;
        }
    }
    gsl_spline_init(FcollLow_spline_LINELUMPOP2E, log10_overdense_spline_LINELUMPOP2E, log10_Fcoll_spline_LINELUMPOP2E, NSFR_low);

    for(i=0;i<NSFR_high;i++) {
        Overdense_spline_LINELUMPOP2E[i] = overdense_large_low + (float)i/((float)NSFR_high-1.)*(overdense_large_high - overdense_large_low);
        Fcoll_spline_LINELUMPOP2E[i] = FgtrConditionalM_LINELUMPOP2E(z, log(Mmin), log(Mmax), Deltac, Overdense_spline_LINELUMPOP2E[i], MassTurnover, Alpha_star, Fstar10, Mlim_Fstar, LineID);
        if(Fcoll_spline_LINELUMPOP2E[i]<0.) {
            Fcoll_spline_LINELUMPOP2E[i]=pow(10.,-40.0);
        }
    }
    spline(Overdense_spline_LINELUMPOP2E-1, Fcoll_spline_LINELUMPOP2E-1, NSFR_high, 0, 0, second_derivs_LINELUMPOP2E-1);
}

void FcollSpline_LINELUMPOP2E(float Overdensity, float *splined_value){
    int i;
    float returned_value;

    if(Overdensity<1.5) {
        if(Overdensity<-1.) {
            returned_value = 0;
        }
        else {
            returned_value = gsl_spline_eval(FcollLow_spline_LINELUMPOP2E, log10(Overdensity+1.), FcollLow_spline_acc_LINELUMPOP2E);
            returned_value = pow(10.,returned_value);
        }
    }
    else {
        if(Overdensity<0.99*Deltac) {
		// Usage of 0.99*Deltac arises due to the fact that close to the critical density, the collapsed fraction becomes a little unstable
		// However, such densities should always be collapsed, so just set f_coll to unity.
		// Additionally, the fraction of points in this regime relative to the entire simulation volume is extremely small.
            splint(Overdense_spline_LINELUMPOP2E-1, Fcoll_spline_LINELUMPOP2E-1, second_derivs_LINELUMPOP2E-1, NSFR_high, Overdensity, &(returned_value));
        }
        else {
            // We set overdensity of regions close to critical density to critical density to avoid instability.
            //returned_value = 1.;
            splint(Overdense_spline_LINELUMPOP2E-1, Fcoll_spline_LINELUMPOP2E-1, second_derivs_LINELUMPOP2E-1, NSFR_high, Overdensity, &(returned_value));
        }
    }
    // We don't need the following ceiling condition any more... -GS
	  //if(returned_value > 1.) returned_value = 1.;
    *splined_value = returned_value;
}

// Interpolation for LINELUMPOP2F -GS
void initialiseFcollLINELUMPOP2F_spline(float z, float Mmax, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar, int LineID){
	float Mmin = MassTurnover/50.;
  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    Mmin = src.minMass(z);
  }
	//float Mmin = 1e8; // TEST
    double overdense_val;
    double overdense_large_high = Deltac, overdense_large_low = 1.5;
    double overdense_small_high = 1.5, overdense_small_low = -1. + 9e-8;
    int i;
    FcollLow_spline_acc_LINELUMPOP2F = gsl_interp_accel_alloc ();
    FcollLow_spline_LINELUMPOP2F = gsl_spline_alloc (gsl_interp_cspline, NSFR_low);
    for (i=0; i<NSFR_low; i++){
        overdense_val = log10(1. + overdense_small_low) + (double)i/((double)NSFR_low-1.)*(log10(1.+overdense_small_high)-log10(1.+overdense_small_low));

        log10_overdense_spline_LINELUMPOP2F[i] = overdense_val;
        log10_Fcoll_spline_LINELUMPOP2F[i] = log10(GaussLegendreQuad_FcollLINELUMPOP2F(NGL_SFR, z, log(Mmax), Deltac, pow(10.,overdense_val)-1., MassTurnover, Alpha_star, Fstar10, Mlim_Fstar, LineID));

        if(log10_Fcoll_spline_LINELUMPOP2F[i] < -40.){
            log10_Fcoll_spline_LINELUMPOP2F[i] = -40.;
        }
    }
    gsl_spline_init(FcollLow_spline_LINELUMPOP2F, log10_overdense_spline_LINELUMPOP2F, log10_Fcoll_spline_LINELUMPOP2F, NSFR_low);

    for(i=0;i<NSFR_high;i++) {
        Overdense_spline_LINELUMPOP2F[i] = overdense_large_low + (float)i/((float)NSFR_high-1.)*(overdense_large_high - overdense_large_low);
        Fcoll_spline_LINELUMPOP2F[i] = FgtrConditionalM_LINELUMPOP2F(z, log(Mmin), log(Mmax), Deltac, Overdense_spline_LINELUMPOP2F[i], MassTurnover, Alpha_star, Fstar10, Mlim_Fstar, LineID);
        if(Fcoll_spline_LINELUMPOP2F[i]<0.) {
            Fcoll_spline_LINELUMPOP2F[i]=pow(10.,-40.0);
        }
    }
    spline(Overdense_spline_LINELUMPOP2F-1, Fcoll_spline_LINELUMPOP2F-1, NSFR_high, 0, 0, second_derivs_LINELUMPOP2F-1);
}

void FcollSpline_LINELUMPOP2F(float Overdensity, float *splined_value){
    int i;
    float returned_value;

    if(Overdensity<1.5) {
        if(Overdensity<-1.) {
            returned_value = 0;
        }
        else {
            returned_value = gsl_spline_eval(FcollLow_spline_LINELUMPOP2F, log10(Overdensity+1.), FcollLow_spline_acc_LINELUMPOP2F);
            returned_value = pow(10.,returned_value);
        }
    }
    else {
        if(Overdensity<0.99*Deltac) {
		// Usage of 0.99*Deltac arises due to the fact that close to the critical density, the collapsed fraction becomes a little unstable
		// However, such densities should always be collapsed, so just set f_coll to unity.
		// Additionally, the fraction of points in this regime relative to the entire simulation volume is extremely small.
            splint(Overdense_spline_LINELUMPOP2F-1, Fcoll_spline_LINELUMPOP2F-1, second_derivs_LINELUMPOP2F-1, NSFR_high, Overdensity, &(returned_value));
        }
        else {
            // We set overdensity of regions close to critical density to critical density to avoid instability.
            //returned_value = 1.;
            splint(Overdense_spline_LINELUMPOP2F-1, Fcoll_spline_LINELUMPOP2F-1, second_derivs_LINELUMPOP2F-1, NSFR_high, Overdensity, &(returned_value));
        }
    }
    // We don't need the following ceiling condition any more... -GS
	  //if(returned_value > 1.) returned_value = 1.;
    *splined_value = returned_value;
}

// Interpolation for LINELUMPOP2G -GS
void initialiseFcollLINELUMPOP2G_spline(float z, float Mmax, float MassTurnover, float Alpha_star, float Fstar10, float Mlim_Fstar, int LineID){
	float Mmin = MassTurnover/50.;
  if(USE_GENERAL_SOURCES)
  {
    sources src;
    src = defaultSources();
    Mmin = src.minMass(z);
  }
	//float Mmin = 1e8; // TEST
    double overdense_val;
    double overdense_large_high = Deltac, overdense_large_low = 1.5;
    double overdense_small_high = 1.5, overdense_small_low = -1. + 9e-8;
    int i;
    FcollLow_spline_acc_LINELUMPOP2G = gsl_interp_accel_alloc ();
    FcollLow_spline_LINELUMPOP2G = gsl_spline_alloc (gsl_interp_cspline, NSFR_low);
    for (i=0; i<NSFR_low; i++){
        overdense_val = log10(1. + overdense_small_low) + (double)i/((double)NSFR_low-1.)*(log10(1.+overdense_small_high)-log10(1.+overdense_small_low));

        log10_overdense_spline_LINELUMPOP2G[i] = overdense_val;
        log10_Fcoll_spline_LINELUMPOP2G[i] = log10(GaussLegendreQuad_FcollLINELUMPOP2G(NGL_SFR, z, log(Mmax), Deltac, pow(10.,overdense_val)-1., MassTurnover, Alpha_star, Fstar10, Mlim_Fstar, LineID));

        if(log10_Fcoll_spline_LINELUMPOP2G[i] < -40.){
            log10_Fcoll_spline_LINELUMPOP2G[i] = -40.;
        }
    }
    gsl_spline_init(FcollLow_spline_LINELUMPOP2G, log10_overdense_spline_LINELUMPOP2G, log10_Fcoll_spline_LINELUMPOP2G, NSFR_low);

    for(i=0;i<NSFR_high;i++) {
        Overdense_spline_LINELUMPOP2G[i] = overdense_large_low + (float)i/((float)NSFR_high-1.)*(overdense_large_high - overdense_large_low);
        Fcoll_spline_LINELUMPOP2G[i] = FgtrConditionalM_LINELUMPOP2G(z, log(Mmin), log(Mmax), Deltac, Overdense_spline_LINELUMPOP2G[i], MassTurnover, Alpha_star, Fstar10, Mlim_Fstar, LineID);
        if(Fcoll_spline_LINELUMPOP2G[i]<0.) {
            Fcoll_spline_LINELUMPOP2G[i]=pow(10.,-40.0);
        }
    }
    spline(Overdense_spline_LINELUMPOP2G-1, Fcoll_spline_LINELUMPOP2G-1, NSFR_high, 0, 0, second_derivs_LINELUMPOP2G-1);
}

void FcollSpline_LINELUMPOP2G(float Overdensity, float *splined_value){
    int i;
    float returned_value;

    if(Overdensity<1.5) {
        if(Overdensity<-1.) {
            returned_value = 0;
        }
        else {
            returned_value = gsl_spline_eval(FcollLow_spline_LINELUMPOP2G, log10(Overdensity+1.), FcollLow_spline_acc_LINELUMPOP2G);
            returned_value = pow(10.,returned_value);
        }
    }
    else {
        if(Overdensity<0.99*Deltac) {
		// Usage of 0.99*Deltac arises due to the fact that close to the critical density, the collapsed fraction becomes a little unstable
		// However, such densities should always be collapsed, so just set f_coll to unity.
		// Additionally, the fraction of points in this regime relative to the entire simulation volume is extremely small.
            splint(Overdense_spline_LINELUMPOP2G-1, Fcoll_spline_LINELUMPOP2G-1, second_derivs_LINELUMPOP2G-1, NSFR_high, Overdensity, &(returned_value));
        }
        else {
            // We set overdensity of regions close to critical density to critical density to avoid instability.
            //returned_value = 1.;
            splint(Overdense_spline_LINELUMPOP2G-1, Fcoll_spline_LINELUMPOP2G-1, second_derivs_LINELUMPOP2G-1, NSFR_high, Overdensity, &(returned_value));
        }
    }
    // We don't need the following ceiling condition any more... -GS
	  //if(returned_value > 1.) returned_value = 1.;
    *splined_value = returned_value;
}

  // For Ts.c : The functions below this line are under development.
void initialise_FgtrM_st_SFR_spline(int Nbin, float zmin, float zmax, float MassTurn, float Alpha_star, float Alpha_esc, float Fstar10, float Fesc10){
	int i;
	float Mmin = MassTurn/50., Mmax = 1e16;
	float Mlim_Fstar, Mlim_Fesc;

	Mlim_Fstar = Mass_limit_bisection(Mmin, Mmax, Alpha_star, Fstar10);
	Mlim_Fesc = Mass_limit_bisection(Mmin, Mmax, Alpha_esc, Fesc10);

	Fcollz_spline_acc = gsl_interp_accel_alloc ();
	Fcollz_spline = gsl_spline_alloc (gsl_interp_cspline, Nbin);
	for (i=0; i<Nbin; i++){
		z_val[i] = zmin + (double)i/((double)Nbin-1.)*(zmax - zmin);
		Fcollz_val[i] = FgtrM_st_SFR(z_val[i], MassTurn, Alpha_star, Alpha_esc, Fstar10, Fesc10, Mlim_Fstar, Mlim_Fesc);
	}
	gsl_spline_init(Fcollz_spline, z_val, Fcollz_val, Nbin);
}

void FgtrM_st_SFR_z(float z, float *splined_value){
	float returned_value;

	returned_value = gsl_spline_eval(Fcollz_spline, z, Fcollz_spline_acc);
	*splined_value = returned_value;
}

void initialise_Xray_FgtrM_st_SFR_spline(int Nbin, float zmin, float zmax, float MassTurn, float Alpha_star, float Fstar10){
	int i;
	float Mmin = MassTurn/50., Mmax = 1e16;
	float Mlim_Fstar;

	Mlim_Fstar = Mass_limit_bisection(Mmin, Mmax, Alpha_star, Fstar10);

	FcollzX_spline_acc = gsl_interp_accel_alloc ();
	FcollzX_spline = gsl_spline_alloc (gsl_interp_cspline, Nbin);
	for (i=0; i<Nbin; i++){
		z_X_val[i] = zmin + (double)i/((double)Nbin-1.)*(zmax - zmin);
		FcollzX_val[i] = FgtrM_st_SFR_Xray(z_val[i], MassTurn, Alpha_star, 0., Fstar10, 1.,Mlim_Fstar,0.);
	}
	gsl_spline_init(FcollzX_spline, z_X_val, FcollzX_val, Nbin);
}

void initialise_Lya_FgtrM_st_SFR_spline(int Nbin, float zmin, float zmax, float MassTurn, float Alpha_star, float Fstar10){
  int i;
  float Mmin = MassTurn/50., Mmax = 1e16;
  float Mlim_Fstar;

  Mlim_Fstar = Mass_limit_bisection(Mmin, Mmax, Alpha_star, Fstar10);

  FcollzLya_spline_acc = gsl_interp_accel_alloc ();
  FcollzLya_spline = gsl_spline_alloc (gsl_interp_cspline, Nbin);
  for (i=0; i<Nbin; i++){
    z_X_val[i] = zmin + (double)i/((double)Nbin-1.)*(zmax - zmin);
    FcollzLya_val[i] = FgtrM_st_SFR_Lya(z_val[i], MassTurn, Alpha_star, 0., Fstar10, 1.,Mlim_Fstar,0.);
  }
  gsl_spline_init(FcollzLya_spline, z_X_val, FcollzLya_val, Nbin);
}

void FgtrM_st_SFR_Xray_z(float z, float *splined_value){
	float returned_value;

	returned_value = gsl_spline_eval(FcollzX_spline, z, FcollzX_spline_acc);
	*splined_value = returned_value;
}

void FgtrM_st_SFR_Lya_z(float z, float *splined_value){
  float returned_value;

  returned_value = gsl_spline_eval(FcollzLya_spline, z, FcollzLya_spline_acc);
  *splined_value = returned_value;
}


void initialise_Xray_Fcollz_SFR_Conditional_table(int Nsteps_zp, int Nfilter, float z[], double R[], float MassTurnover, float Alpha_star, float Fstar10){
    double overdense_val;
    double overdense_large_high = Deltac, overdense_large_low = 1.5;
    double overdense_small_high = 1.5, overdense_small_low = -1. + 9e-8;
	double overdense_low_table[NSFR_low];
	float Mmin,Mmax,Mlim_Fstar;
    int i,j,k,i_tot;
    //int Nfilter = n;
    sources src;
    src = defaultSources();
    Mmin = MassTurnover/50;
    if(USE_GENERAL_SOURCES)
    {
      Mmin = 1e5; // just set this to something really low to make sure we get the full range of possible masses in initialiseSplinedSigmaM
    }
	Mmax = RtoM(R[Nfilter-1]);
	Mlim_Fstar = Mass_limit_bisection(Mmin, Mmax, Alpha_star, Fstar10);
    initialiseSplinedSigmaM(Mmin,Mmax);
    fprintf(stderr, "In initialise_Fcollz_SFR_Conditional_table: Rmin = %6.4f, Rmax = %6.4f, Mmin = %.4e, Mmax = %.4e\n",
																R[0],R[Nfilter-1],RtoM(R[0]),RtoM(R[Nfilter-1]));
    for (i=0; i<NSFR_low; i++) {
      overdense_val = log10(1. + overdense_small_low) + (double)i/((double)NSFR_low-1.)*(log10(1.+overdense_small_high)-log10(1.+overdense_small_low));
      log10_overdense_low_table[i] = overdense_val;
	  overdense_low_table[i] = pow(10.,log10_overdense_low_table[i]);
    }
    for (i=0; i<NSFR_high;i++) {
      Overdense_high_table[i] = overdense_large_low + (float)i/((float)NSFR_high-1.)*(overdense_large_high - overdense_large_low);
    }
    //    printf("Nfilter = %d\n",Nfilter);
    for (k=0; k < Nsteps_zp; k++) {
	  i_tot = Nfilter*k;
      //printf("R filter = %d\n",k);
      //Mmax = RtoM(R[k]);
      for (j=0; j < Nfilter; j++) {
        Mmax = RtoM(R[j]);
        if(USE_GENERAL_SOURCES) Mmin = src.minMass(z[i_tot + j]);
        initialiseGL_FcollSFR(NGL_SFR, Mmin, Mmax, 0);
        for (i=0; i<NSFR_low; i++){
            log10_FcollzXray_SFR_low_table[i_tot+j][i] = log10(GaussLegendreQuad_FcollSFR_Xray(NGL_SFR,z[i_tot+j],log(Mmax),Deltac,overdense_low_table[i]-1.,MassTurnover,Alpha_star,0.,Fstar10,1.,Mlim_Fstar,0.));
            if(log10_FcollzXray_SFR_low_table[i_tot+j][i] < -40.) log10_FcollzXray_SFR_low_table[i_tot+j][i] = -40.;
		    //if(i==0)printf("zp_step=%d, filter_step=%d, arr=%d\n",k, j,i_tot+j); //TEST
        }

        for(i=0;i<NSFR_high;i++) {
            FcollzXray_SFR_high_table[i_tot+j][i] = FgtrConditionalM_SFR_Xray(z[i_tot+j],log(Mmin),log(Mmax),Deltac,Overdense_high_table[i],MassTurnover,Alpha_star,0.,Fstar10,1.,Mlim_Fstar,0.);
            if(FcollzXray_SFR_high_table[i_tot+j][i]<0.) FcollzXray_SFR_high_table[i_tot+j][i]=pow(10.,-40.0);
        }
      }
    }
}

void initialise_Lya_Fcollz_SFR_Conditional_table(int Nsteps_zp, int Nfilter, float z[], double R[], float MassTurnover, float Alpha_star, float Fstar10){
    double overdense_val;
    double overdense_large_high = Deltac, overdense_large_low = 1.5;
    double overdense_small_high = 1.5, overdense_small_low = -1. + 9e-8;
  double overdense_low_table[NSFR_low];
  float Mmin,Mmax,Mlim_Fstar;
    int i,j,k,i_tot;
    //int Nfilter = n;
    sources src;
    src = defaultSources();
    Mmin = MassTurnover/50;
    if(USE_GENERAL_SOURCES)
    {
      Mmin = 1e5; // just set this to something really low to make sure we get the full range of possible masses in initialiseSplinedSigmaM
    }
  Mmax = RtoM(R[Nfilter-1]);
  Mlim_Fstar = Mass_limit_bisection(Mmin, Mmax, Alpha_star, Fstar10);
    initialiseSplinedSigmaM(Mmin,Mmax);
    fprintf(stderr, "In initialise_Fcollz_SFR_Conditional_table: Rmin = %6.4f, Rmax = %6.4f, Mmin = %.4e, Mmax = %.4e\n",
                                R[0],R[Nfilter-1],RtoM(R[0]),RtoM(R[Nfilter-1]));
    for (i=0; i<NSFR_low; i++) {
      overdense_val = log10(1. + overdense_small_low) + (double)i/((double)NSFR_low-1.)*(log10(1.+overdense_small_high)-log10(1.+overdense_small_low));
      log10_overdense_low_table[i] = overdense_val;
    overdense_low_table[i] = pow(10.,log10_overdense_low_table[i]);
    }
    for (i=0; i<NSFR_high;i++) {
      Overdense_high_table[i] = overdense_large_low + (float)i/((float)NSFR_high-1.)*(overdense_large_high - overdense_large_low);
    }
    //    printf("Nfilter = %d\n",Nfilter);
    for (k=0; k < Nsteps_zp; k++) {
    i_tot = Nfilter*k;
      //printf("R filter = %d\n",k);
      //Mmax = RtoM(R[k]);
      for (j=0; j < Nfilter; j++) {
        Mmax = RtoM(R[j]);
        if(USE_GENERAL_SOURCES) Mmin = src.minMass(z[i_tot + j]);
        initialiseGL_FcollSFR(NGL_SFR, Mmin, Mmax, 0);
        for (i=0; i<NSFR_low; i++){
            log10_FcollzLya_SFR_low_table[i_tot+j][i] = log10(GaussLegendreQuad_FcollSFR_Lya(NGL_SFR,z[i_tot+j],log(Mmax),Deltac,overdense_low_table[i]-1.,MassTurnover,Alpha_star,0.,Fstar10,1.,Mlim_Fstar,0.));
            if(log10_FcollzLya_SFR_low_table[i_tot+j][i] < -40.) log10_FcollzLya_SFR_low_table[i_tot+j][i] = -40.;
        //if(i==0)printf("zp_step=%d, filter_step=%d, arr=%d\n",k, j,i_tot+j); //TEST
        }

        for(i=0;i<NSFR_high;i++) {
            FcollzLya_SFR_high_table[i_tot+j][i] = FgtrConditionalM_SFR_Lya(z[i_tot+j],log(Mmin),log(Mmax),Deltac,Overdense_high_table[i],MassTurnover,Alpha_star,0.,Fstar10,1.,Mlim_Fstar,0.);
            if(FcollzLya_SFR_high_table[i_tot+j][i]<0.) FcollzLya_SFR_high_table[i_tot+j][i]=pow(10.,-40.0);
        }
      }
    }
}

void free_interpolation() {
    gsl_spline_free (FcollzX_spline);
    gsl_spline_free (FcollzLya_spline);
    gsl_interp_accel_free (FcollzX_spline_acc);
    gsl_interp_accel_free (FcollzLya_spline_acc);
    gsl_spline_free (Fcollz_spline);
    gsl_interp_accel_free (Fcollz_spline_acc);
}
// For Ts.c : The functions above this line are under development.

/* New in v1.4: end */

#endif
