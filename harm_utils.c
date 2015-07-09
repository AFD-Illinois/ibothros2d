
#include "decs.h"

/** HARM utilities **/

extern double ***bcon;
extern double ***bcov;
extern double ***ucon;
extern double ***ucov;
extern double ***p;
extern double **ne;
extern double **thetae;
extern double **b;

void Xtoij(double X[NDIM], int *i, int *j, double del[NDIM]) ;

/********************************************************************

				Interpolation routines

 ********************************************************************/

/* return fluid four-vector in simulation units */
void interp_fourv(double X[NDIM], double ***fourv, double Fourv[NDIM])
{
	double del[NDIM],d1,d2,d3,d4;
	int i, j, ip1, jp1;

	/* find the current zone location and offsets del[0], del[1] */
	Xtoij(X, &i, &j, del);

	ip1 = i + 1;
	jp1 = j + 1;

	d1 = (1. - del[1]) * (1. - del[2]);
	d3 = del[1] * (1. - del[2]);
	d2 = del[2] * (1. - del[1]);
	d4 = del[1] * del[2];

	Fourv[0] = d1*fourv[i][j][0];
	Fourv[1] = d1*fourv[i][j][1];
	Fourv[2] = d1*fourv[i][j][2];
	Fourv[3] = d1*fourv[i][j][3];

	Fourv[0] += d2*fourv[i][jp1][0];
	Fourv[1] += d2*fourv[i][jp1][1];
	Fourv[2] += d2*fourv[i][jp1][2];
	Fourv[3] += d2*fourv[i][jp1][3];

	Fourv[0] += d3*fourv[ip1][j][0];
	Fourv[1] += d3*fourv[ip1][j][1];
	Fourv[2] += d3*fourv[ip1][j][2];
	Fourv[3] += d3*fourv[ip1][j][3];

	Fourv[0] += d4*fourv[ip1][jp1][0];
	Fourv[1] += d4*fourv[ip1][jp1][1];
	Fourv[2] += d4*fourv[ip1][jp1][2];
	Fourv[3] += d4*fourv[ip1][jp1][3];

}

/* return scalar in cgs units */
double interp_scalar(double X[NDIM], double **var)
{
	double del[NDIM];
	int i, j;

	/* find the current zone location and offsets del[0], del[1] */
	Xtoij(X, &i, &j, del);

	/* use bilinear interpolation to find rho; piecewise constant
	   near the boundaries */
	return(
		(
			  (1. - del[1]) * (1. - del[2]) * var[i][j]
			+ (1. - del[1]) *       del[2]  * var[i][j + 1] 
			+ del[1] *        (1. - del[2]) * var[i + 1][j]  
			+ del[1] *              del[2]  * var[i + 1][j + 1]
			)
	) ;

}

/***********************************************************************************

					End interpolation routines

 ***********************************************************************************/


void Xtoij(double X[NDIM], int *i, int *j, double del[NDIM])
{

	*i = (int) ((X[1] - startx[1]) / dx[1] - 0.5 + 1000) - 1000;
	*j = (int) ((X[2] - startx[2]) / dx[2] - 0.5 + 1000) - 1000;

	if(*i < 0) {
		*i = 0 ;
		del[1] = 0. ;
	}
	else if(*i > N1-2) {
		*i = N1-2 ;
		del[1] = 1. ;
	}
	else {
		del[1] = (X[1] - ((*i + 0.5) * dx[1] + startx[1])) / dx[1];
	}

	if(*j < 0) {
		*j = 0 ;
		del[2] = 0. ;
	}
	else if(*j > N2-2) {
		*j = N2-2 ;
		del[2] = 1. ;
	}
	else {
		del[2] = (X[2] - ((*j + 0.5) * dx[2] + startx[2])) / dx[2];
	}

	return;
}

/* return boyer-lindquist coordinate of point */
void bl_coord(double *X, double *r, double *th)
{

	*r = exp(X[1]) + R0 ;
	*th = M_PI * X[2] + ((1. - hslope) / 2.) * sin(2. * M_PI * X[2]);

	return;
}

void coord(int i, int j, double *X)
{

	/* returns zone-centered values for coordinates */
	X[0] = startx[0];
	X[1] = startx[1] + (i + 0.5) * dx[1];
	X[2] = startx[2] + (j + 0.5) * dx[2];
	X[3] = startx[3];

	return;
}

void set_units(char *instr) {
    
    double two_temp_gam ;

	MBH *= MSUN ; /* for Sgr A* */

	/** from this, calculate units of length, time, mass,
	    and derivative units **/
	L_unit = GNEWT * MBH / (CL * CL);
	T_unit = L_unit / CL;

	fprintf(stderr,"\nUNITS\n") ;
	fprintf(stderr,"L,T,M: %g %g %g\n",L_unit,T_unit,M_unit) ;

	RHO_unit = M_unit / pow(L_unit, 3);
	U_unit = RHO_unit * CL * CL;
	B_unit = CL * sqrt(4.*M_PI*RHO_unit);
	Ne_unit = RHO_unit / (MP + ME) ;

	two_temp_gam = 0.5 * ((1. + 2. / 3. * (TP_OVER_TE + 1.) / 
		(TP_OVER_TE + 2.)) + gam);
        Thetae_unit = (two_temp_gam - 1.) * (MP / ME) / (1. + TP_OVER_TE);

	fprintf(stderr,"rho,u,B: %g %g %g\n",RHO_unit,U_unit,B_unit) ;
}

void init_physical_quantities(void)
{
	int i, j;

	for (i = 0; i < N1; i++) 
	for (j = 0; j < N2; j++) {

		ne[i][j] = p[i][j][KRHO] * Ne_unit ;

		/* NOTE WELL: first 0.5 is a fudge factor */
		thetae[i][j] = (p[i][j][UU]/p[i][j][KRHO])* Thetae_unit ;
		if(thetae[i][j] > THETAE_MAX) thetae[i][j] = THETAE_MAX ;
		b[i][j] = sqrt( bcon[i][j][0] * bcov[i][j][0] + 
			        bcon[i][j][1] * bcov[i][j][1] + 
			        bcon[i][j][2] * bcov[i][j][2] + 
			        bcon[i][j][3] * bcov[i][j][3] ) * B_unit ;

	}

	return ;
}

void *malloc_rank1(int n1, int size)
{
	void *A;

	if((A = malloc(n1*size)) == NULL){
		fprintf(stderr,"malloc failure in malloc_rank1\n");
		exit(123);
	}

	return A;
}

void **malloc_rank2(int n1, int n2, int size)
{

	void **A;
	int i;

	if((A=malloc(n1*sizeof(void *))) == NULL) {
		fprintf(stderr,"malloc failure in malloc_rank2\n");
		exit(124);
	}

	for(i=0;i<n1;i++){
		A[i] = malloc_rank1(n2, size);
	}

	return A;
}

void ***malloc_rank3(int n1, int n2, int n3, int size)
{

	void ***A;
	int i;

	if((A=malloc(n1*sizeof(void *))) == NULL) {
		fprintf(stderr,"malloc failure in malloc_rank3\n");
		exit(125);
	}

	for(i=0;i<n1;i++){
		A[i] = malloc_rank2(n2, n3, size);
	}

	return A;
}

void ****malloc_rank4(int n1, int n2, int n3, int n4, int size)
{

	void ****A;
	int i;

	if((A=malloc(n1*sizeof(void *))) == NULL) {
		fprintf(stderr,"malloc failure in malloc_rank4\n");
		exit(126);
	}

	for(i=0;i<n1;i++){
		A[i] = malloc_rank3(n2, n3, n4, size);
	}

	return A;
}

void init_storage(void)
{

	bcon = (double ***)malloc_rank3(N1,N2,NDIM,sizeof(double));
	bcov = (double ***)malloc_rank3(N1,N2,NDIM,sizeof(double));
	ucon = (double ***)malloc_rank3(N1,N2,NDIM,sizeof(double));
	ucov = (double ***)malloc_rank3(N1,N2,NDIM,sizeof(double));
	p = (double ***)malloc_rank3(N1,N2,NPRIM,sizeof(double));
	ne = (double **)malloc_rank2(N1,N2,sizeof(double));
	thetae = (double **)malloc_rank2(N1,N2,sizeof(double));
	b = (double **)malloc_rank2(N1,N2,sizeof(double));

	return;
}

