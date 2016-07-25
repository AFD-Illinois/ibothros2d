
#include "decs.h"

/*

get iharm3d_v3 simulation data from fname

checks for consistency of coordinates in data file with
values of coordinate parameters 

Uses standard iharm3d_v3 data file format

CFG 29 apr 16

*/

extern double ***bcon;
extern double ***bcov;
extern double ***ucon;
extern double ***ucov;
extern double ***p;
extern double **ne;
extern double **thetae;
extern double **b;

void init_harm_data(char *fname)
{
	FILE *fp;
	double x[4];
	int i, j, k ;

	/* header variables not used except locally */
	double t,tf,cour,DTd,DTl,DTi,dt ;
	int nstep,DTr,dump_cnt,image_cnt,rdump_cnt;
	double divb,vmin,vmax,gdet ;
	double J ;
	double dMact, Ladv ;

	fp = fopen(fname, "r");

	if (fp == NULL) {
		fprintf(stderr, "can't open sim data file\n");
		exit(1);
	} else {
		fprintf(stderr,
			"successfully opened %s\n",fname);
	}

	/* get standard HARM header */
	fscanf(fp, "%lf ", &t) ;
	fscanf(fp, "%d ",  &N1) ;
	fscanf(fp, "%d ",  &N2) ;
	fscanf(fp, "%lf ", &startx[1]) ;
	fscanf(fp, "%lf ", &startx[2]) ;
	fscanf(fp, "%lf ", &dx[1]) ;
	fscanf(fp, "%lf ", &dx[2]) ;
	fscanf(fp, "%lf ", &tf) ;
	fscanf(fp, "%d ",  &nstep) ;
	fscanf(fp, "%lf ", &gam) ;
	fscanf(fp, "%lf ", &cour) ;
	fscanf(fp, "%lf ", &DTd) ;
	fscanf(fp, "%lf ", &DTl) ;
	fscanf(fp, "%lf ", &DTi) ;
	fscanf(fp, "%d ",  &DTr) ;
	fscanf(fp, "%d ",  &dump_cnt) ;
	fscanf(fp, "%d ",  &image_cnt) ;
	fscanf(fp, "%d ",  &rdump_cnt) ;
	fscanf(fp, "%lf ", &dt) ;
	
	double dum;
        for(i=0;i<21;i++) fscanf(fp, "%lf",&dum);

	/* not set automatically */
	a = 0.9375;
	Rin = 0.98 * (1. + sqrt(1. - a*a)) ;
	Rout = 40.;
	hslope = 0.3;
	R0 = 0.0;

	fprintf(stderr,"coordinate parameters: Rin,Rout,hslope,R0,dx[1],dx[2]: %g %g %g %g %g %g\n",
		Rin,Rout,hslope,R0,dx[1],dx[2]) ;

	/* nominal non-zero values for axisymmetric simulations */
	startx[0] = 0.;
	startx[3] = 0.;

	stopx[0] = 1. ;
	stopx[1] = startx[1] + N1*dx[1] ;
	stopx[2] = startx[2] + N2*dx[2] ;
	stopx[3] = 2.*M_PI ;

	dx[0] = 1. ;
	dx[3] = 2.*M_PI;

	/* Allocate storage for all model size dependent variables */
	init_storage();

	dMact = 0. ;
	Ladv = 0. ;
	for(k = 0 ; k < N1*N2 ; k++) {
			j = k%N2 ;
			i = (k - j)/N2 ;
			fscanf(fp, "%lf %lf ", &x[1], &x[2]);
			fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf",
			       &p[i][j][KRHO], 
			       &p[i][j][UU],
			       &p[i][j][U1], 
			       &p[i][j][U2], 
			       &p[i][j][U3],
			       &p[i][j][B1], 
			       &p[i][j][B2], 
			       &p[i][j][B3]);

			fscanf(fp, "%lf", &divb);

			fscanf(fp, "%lf %lf %lf %lf",
			       &ucon[i][j][0], &ucon[i][j][1],
			       &ucon[i][j][2], &ucon[i][j][3]);
			fscanf(fp, "%lf %lf %lf %lf", &ucov[i][j][0],
			       &ucov[i][j][1], &ucov[i][j][2],
			       &ucov[i][j][3]);
			fscanf(fp, "%lf %lf %lf %lf", &bcon[i][j][0],
			       &bcon[i][j][1], &bcon[i][j][2],
			       &bcon[i][j][3]);
			fscanf(fp, "%lf %lf %lf %lf", &bcov[i][j][0],
			       &bcov[i][j][1], &bcov[i][j][2],
			       &bcov[i][j][3]);
			fscanf(fp, "%lf ", &vmin);
			fscanf(fp, "%lf ", &vmax);
			fscanf(fp, "%lf ", &vmin);
			fscanf(fp, "%lf ", &vmax);
			fscanf(fp, "%lf ", &gdet);

			/* additional stuff: current */
			fscanf(fp, "%lf ", &J) ;
			fscanf(fp, "%lf ", &J) ;
			fscanf(fp, "%lf ", &J) ;
			fscanf(fp, "%lf ", &J) ;

			fscanf(fp, "%lf ", &J) ;
			fscanf(fp, "%lf ", &J) ;
			fscanf(fp, "%lf ", &J) ;
			fscanf(fp, "%lf\n", &J) ;

			/* check accretion rate */
			if(i <= 20) dMact += gdet * p[i][j][KRHO] * ucon[i][j][1] ;
			if(i >= 20 && i < 40) Ladv += gdet * p[i][j][UU] * ucon[i][j][1] * ucov[i][j][0] ;

		}

	dMact *= dx[3]*dx[2] ;
	dMact /= 21. ;
	Ladv *= dx[3]*dx[2] ;
	Ladv /= 21. ;

	fprintf(stderr,"dMact: %g, Ladv: %g\n",dMact,Ladv) ;

	/* done! */
}
