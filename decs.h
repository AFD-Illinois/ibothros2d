
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include "constants.h"

#define NX	256
#define NY	256

#define NDIM	4
#define NPRIM	8	

/* mnemonics for primitive vars; conserved vars */
#define KRHO     0
#define UU      1
#define U1      2
#define U2      3
#define U3      4
#define B1      5
#define B2      6
#define B3      7

/* numerical convenience */
#define SMALL	1.e-40

extern int pflag ;

/* some coordinate parameters */
extern double a;
extern double R0 ;
extern double Rin ;
extern double Rout ;
extern double Rh ;
extern double hslope ;
extern double th_len,th_beg;
extern double startx[NDIM], stopx[NDIM], dx[NDIM];
extern double gam ;

/* HARM model globals */
extern double M_unit;
extern double L_unit;
extern double T_unit;
extern double RHO_unit;
extern double U_unit;
extern double B_unit;
extern double Ne_unit;
extern double Thetae_unit;

extern int N1, N2, N3;

/* model parameter globals */
extern double MBH, M_unit, Dsource, TP_OVER_TE, THETAE_MAX;
extern double pin, FRAC, Gmin, Gmax;


/** model-independent subroutines **/
/* core routines */
void init(int i, int j, double Xcam[4], double Ucam[4], double fovx, double fovy,
	double X[4], double Kcon[4]) ;
void null_normalize(double Kcon[NDIM], double fnorm) ;
void normalize(double *vcon, double gcov[][NDIM]);
//double approximate_solve(double Ii, double ji, double ki, double jf, double kf, double dl) ;
double approximate_solve(double Ii, double S, double dtau);
void get_jkinv(double X[NDIM], double Kcon[NDIM], double *jnuinv, double *knuinv) ;
int get_skinv(double X[NDIM], double Kcon[NDIM], double *snuinv, double *knuinv);
int    stop_forward_integration(double X[NDIM], double Kcon[NDIM],
	double Xcam[NDIM]) ;
int    stop_backward_integration(double X[NDIM], double Kcon[NDIM],
	double Xcam[NDIM]) ;
void dump(double image[NX][NY], char *fname) ;

/* geodesic integration */
void   push_photon(double X[NDIM], double Kcon[NDIM], double dl);
void   push_photon_gsl(double X[NDIM], double Kcon[NDIM], double dl);

/* model */
void   gcov_func(double *X, double gcov[][NDIM]);
void   gcon_func(double gcov[][NDIM], double gcon[][NDIM]);
double gdet_func(double gcov[][NDIM]);
void   get_connection(double *X, double lconn[][NDIM][NDIM]);
void   lower(double *ucon, double Gcov[NDIM][NDIM], double *ucov);
double stepsize(double X[NDIM], double K[NDIM]);
void init_model(char *args[]) ;
void get_fluid_params(double X[NDIM], double *Ne, double *Thetae, double *B, double Ucov[NDIM], double Bcov[NDIM]);
double get_model_thetae(double X[NDIM]) ;
double get_model_b(double X[NDIM]) ;
double get_model_ne(double X[NDIM]) ;
void get_model_bcov(double X[NDIM], double Bcov[NDIM]) ;
void get_model_bcon(double X[NDIM], double Bcon[NDIM]) ;
void get_model_ucov(double X[NDIM], double Ucov[NDIM]) ;
void get_model_ucon(double X[NDIM], double Ucon[NDIM]) ;

/* harm2d model */
void init_harm_data(char *fname) ;

/* harm utilities */
/*
void interp_fourv(double X[NDIM], double ***fourv, double Fourv[NDIM]) ;
double interp_scalar(double X[NDIM], double **var) ;
void Xtoij(double X[NDIM], int *i, int *j, double del[NDIM]) ;
*/
void   bl_coord(double *X, double *r, double *th);
//void coord(int i, int j, double *X) ;
void set_units(char *instr) ;
void init_physical_quantities(void) ;

/*
void *malloc_rank1(int n1, int size) ;
void **malloc_rank2(int n1, int n2, int size) ;
void ***malloc_rank3(int n1, int n2, int n3, int size) ;
void ****malloc_rank4(int n1, int n2, int n3, int n4, int size) ;
*/

void init_storage(void) ;

/* tetrad related */
void   coordinate_to_tetrad(double Ecov[NDIM][NDIM], double K[NDIM], double K_tetrad[NDIM]);
void   tetrad_to_coordinate(double Ecov[NDIM][NDIM], double K_tetrad[NDIM], double K[NDIM]);
double delta(int i, int j);
void   make_tetrad(double t0[NDIM], double t1[NDIM], double t2[NDIM],
	double Gcov[NDIM][NDIM], double Econ[NDIM][NDIM], 
	double Ecov[NDIM][NDIM]) ;


/* imaging */
void make_ppm(double p[NX][NY], double freq, char filename[]) ;
void john_pal(double data, double min, double max,
	      int *pRed, int *pGreen, int *pBlue) ;

/* radiation */
double Bnu_inv(double nu, double Thetae) ;
double jnu_inv(double nu, double Thetae, double Ne, double B, double theta) ;
double alpha_inv_abs(double nu, double Thetae, double Ne, double B, double theta) ;
double get_fluid_nu(double Kcon[NDIM], double Ucov[NDIM]) ;
//double get_bk_angle(double X[NDIM], double Kcon[NDIM], double Ucov[NDIM]) ;
double get_bk_angle(double X[NDIM], double Kcon[NDIM], double Ucov[NDIM], double Bcov[NDIM], double B);

/* emissivity */
double jnu_synch(double nu, double Ne, double Thetae, double B, double theta) ;

/* set_model_param.c */
void set_model_param(void);

#define DLOOP  for(k=0;k<NDIM;k++)for(l=0;l<NDIM;l++)

