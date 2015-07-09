
#include "decs.h"

void push_photon_ham(double X[NDIM], double Kcon[][NDIM], double dl[])
{
	double dgcon[NDIM][NDIM][NDIM];
	double dKcov[NDIM], Kcov[NDIM];
	double Xh[NDIM], Kconh[NDIM], Kcovh[NDIM];
	double gcon[NDIM][NDIM], gcov[NDIM][NDIM] ;
	int i, j, k;

	/* advance X */

	/* 2nd order: scheme
	   take half-step and then evaluate derivatives 
	   at half-step. */

	//fprintf(stderr,"---> %g %g %g %g %g %g %g %g\n", X[0], X[1], X[2], X[3], Kcon[0], Kcon[1], Kcon[2], Kcon[3]) ;
	/** half-step **/
	dgcon_calc(X, dgcon) ;

	/* advance K */
	gcov_func(X, gcov) ;
	lower(Kcon[0], gcov, Kcov) ;

	for (k = 0; k < 4; k++)
		dKcov[k] = 0.;

	for (i = 1; i < 3; i++)
		for (j = 0; j < 4; j++)
			for (k = 0; k < 4; k++)
				dKcov[i] -= 0.5 * dl[0] * 0.5 * dgcon[i][j][k] * Kcov[j] * Kcov[k];

	for (k = 0; k < 4; k++)
		Kcovh[k] = Kcov[k] + dKcov[k];

	/* advance X */
	for (i = 0; i < 4; i++)
		Xh[i] = X[i] + 0.5 * dl[0] * Kcon[0][i] ;

	/** full step **/
	dgcon_calc(Xh, dgcon) ;

	/* advance K */
	for (k = 0; k < 4; k++)
		dKcov[k] = 0.;
	for (i = 1; i < 3; i++)
		for (j = 0; j < 4; j++)
			for (k = 0; k < 4; k++)
				dKcov[i] -= dl[0] * 0.5 * dgcon[i][j][k] * Kcovh[j] * Kcovh[k];

	for (k = 0; k < 4; k++)
		Kcov[k] += dKcov[k];

	gcon_func(Xh, gcon) ;

	/* lower does the same thing as raise... */
	lower(Kcovh, gcon, Kconh) ;

	/* advance X */
	for (k = 0; k < 4; k++)
		X[k] += dl[0] * Kconh[k];

	gcon_func(X, gcon) ;
	lower(Kcov, gcon, Kcon[0]) ;

	//fprintf(stderr,"---> %g %g %g %g %g %g %g %g\n", X[0], X[1], X[2], X[3], Kcon[0], Kcon[1], Kcon[2], Kcon[3]) ;

	/* done! */
}

