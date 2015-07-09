
#include "decs.h"


/*

	HARM model specification routines 

*/

/*
	HARM 2d grid functions
*/
double ***bcon;
double ***bcov;
double ***ucon;
double ***ucov;
double ***p;
double **ne;
double **thetae;
double **b;

void interp_fourv(double X[NDIM], double ***fourv, double Fourv[NDIM]) ;
double interp_scalar(double X[NDIM], double **var) ;

/* 
   Current metric: modified Kerr-Schild, squashed in theta
   to give higher resolution at the equator 
*/

/* mnemonics for dimensional indices */
#define TT      0
#define RR      1
#define TH      2
#define PH      3

void gcov_func(double *X, double gcov[][NDIM])
{
	int k, l;
	double sth, cth, s2, rho2;
	double r, th;
	double tfac, rfac, hfac, pfac;

	DLOOP gcov[k][l] = 0.;

	bl_coord(X, &r, &th);

	cth = cos(th);
	sth = fabs(sin(th));
	if (sth < SMALL)
		sth = SMALL;
	s2 = sth * sth;
	rho2 = r * r + a * a * cth * cth;

	/* transformation for Kerr-Schild -> modified Kerr-Schild */
	tfac = 1.;
	rfac = r - R0;
	hfac = M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X[2]);
	pfac = 1.;

	gcov[TT][TT] = (-1. + 2. * r / rho2) * tfac * tfac;
	gcov[TT][1] = (2. * r / rho2) * tfac * rfac;
	gcov[TT][3] = (-2. * a * r * s2 / rho2) * tfac * pfac;

	gcov[1][TT] = gcov[TT][1];
	gcov[1][1] = (1. + 2. * r / rho2) * rfac * rfac;
	gcov[1][3] = (-a * s2 * (1. + 2. * r / rho2)) * rfac * pfac;

	gcov[2][2] = rho2 * hfac * hfac;

	gcov[3][TT] = gcov[TT][3];
	gcov[3][1] = gcov[1][3];
	gcov[3][3] =
	    s2 * (rho2 + a * a * s2 * (1. + 2. * r / rho2)) * pfac * pfac;
}

#undef TT
#undef RR
#undef TH
#undef PH

/* 

   connection calculated analytically for modified Kerr-Schild
   	coordinates 


   this gives the connection coefficient
	\Gamma^{i}_{j,k} = conn[..][i][j][k]
   where i = {1,2,3,4} corresponds to, e.g., {t,ln(r),theta,phi}
*/
void get_connection(double X[4], double lconn[4][4][4])
{
	double r1,r2,r3,r4,x2,sx,cx ;
	double th,sth,cth,sth2,cth2,sth4,cth4,s2th,c2th,c4th ;
	double a2,a3,a4,rho2,rho22,rho23,fac1,fac3,fac4,fac5,fac6,fac52 ;
	double r1api,a2pr1r1m2,r1r1p1,r1fac1_rho23,r1_rho23,fac1_rho23,asth2,fac42_rho2 ;

	r1 = exp(X[1]) ;
	r2 = r1*r1 ;
	r3 = r2*r1 ;
	r4 = r3*r1 ;

	x2 = X[2] ;
	sx = sin(2.*M_PI*x2) ;
	cx = cos(2.*M_PI*x2) ;

	th = M_PI*X[2] + (1. - hslope)*sx/2. ;
	sth = sin(th) ;
	cth = cos(th) ;
	sth2 = sth*sth ;
	cth2 = cth*cth ;
	sth4 = sth2*sth2 ;
	cth4 = cth2*cth2 ;
	s2th = sin(2.*th) ;
	c2th = cos(2.*th) ;
	c4th = cth4 - 6.*cth2*sth2 + sth4 ;

	a2 = a*a ;
	a3 = a2*a ;
	a4 = a2*a2 ;

	rho2 = r2 + a2*cth2 ;
	rho22 = rho2*rho2 ;
	rho23 = rho22*rho2 ;

	fac1 = r2 - a2*cth2 ;
	fac3 = r1*(2.+r1) + a2*cth2 ;
	fac4 = 1. + (1. - hslope)*cx ;
	fac5 = a2 + 2.*r2 + a2*c2th ;
	fac52 = fac5*fac5 ;
	fac6 = a*r1*s2th/(M_PI*rho23*fac4) ;

	r1api = r1*a*M_PI ;
	a2pr1r1m2 = a2 + r1*(r1 - 2.) ;
	r1r1p1 = r1*(r1 + 1.) ;
	r1fac1_rho23 = r1*fac1/rho23 ;
	r1_rho23 = r1/rho23;
	fac1_rho23 = fac1/rho23 ;
	asth2 = a*sth2 ;
	fac42_rho2 = fac4*fac4/rho2;

	lconn[0][0][0] = 2.*r1fac1_rho23 ;
	lconn[0][0][1] = fac3*r1fac1_rho23 ;
	lconn[0][0][2] = -r1api*a*fac4*s2th/rho22 ;
	lconn[0][0][3] = -2.*asth2*r1fac1_rho23 ;

	lconn[0][1][0] = lconn[0][0][1] ;
	lconn[0][1][1] = 2.*r2*(r3 + r4 - a2*r1*cth2 - a4*cth4)/rho23 ;
	lconn[0][1][2] = -a2*M_PI*r2*fac4*s2th/rho22 ;
	lconn[0][1][3] = -asth2*r1fac1_rho23*fac3 ;

	lconn[0][2][0] = lconn[0][0][2] ;
	lconn[0][2][1] = lconn[0][1][2] ;
	lconn[0][2][2] = -2.*r2*M_PI*M_PI*fac42_rho2 ;
	lconn[0][2][3] = 2.*a2*r1api*cth*fac4*sth2*sth/rho22 ;

	lconn[0][3][0] = lconn[0][0][3] ;
	lconn[0][3][1] = lconn[0][1][3] ;
	lconn[0][3][2] = lconn[0][2][3] ;
	lconn[0][3][3] = 2.*r1_rho23*sth2*(-r1*rho22 + a2*fac1*sth2) ;

	lconn[1][0][0] = fac1_rho23*a2pr1r1m2/r1 ;
	lconn[1][0][1] = fac1_rho23*(-2.*r1 + a2*sth2) ;
	lconn[1][0][2] = 0. ;
	lconn[1][0][3] = -asth2*a2pr1r1m2*fac1_rho23/r1 ;

	lconn[1][1][0] = lconn[1][0][1] ;
	lconn[1][1][1] = 1. + 4.*r1*(a2 - 2.*r2 + a2*c2th)*(r1*(2.+r1)+a2*c2th)/(fac52*fac5) ;
	lconn[1][1][2] = -a2*M_PI*fac4*s2th/fac5 ;
	lconn[1][1][3] = asth2*(a4*r1*cth4 + r2*(2.*r1 + r3 - a2*sth2) + 
				a2*cth2*(2.*r1*(r2 - 1.) + a2*sth2))/rho23 ;

	lconn[1][2][0] = lconn[1][0][2] ;
	lconn[1][2][1] = lconn[1][1][2] ;
	lconn[1][2][2] = -M_PI*M_PI*a2pr1r1m2*fac42_rho2 ;
	lconn[1][2][3] = 0. ;

	lconn[1][3][0] = lconn[1][0][3] ;
	lconn[1][3][1] = lconn[1][1][3] ;
	lconn[1][3][2] = lconn[1][2][3] ;
	lconn[1][3][3] = -(a2pr1r1m2*sth2*(r1*rho22 - a2*fac1*sth2))/(r1*rho23) ;

	lconn[2][0][0] = -a*fac6 ;
	lconn[2][0][1] = r1*lconn[2][0][0] ;
	lconn[2][0][2] = 0. ;
	lconn[2][0][3] = (a2 + r2)*fac6 ;

	lconn[2][1][0] = lconn[2][0][1] ;
	lconn[2][1][1] = -r2*a*fac6 ;
	lconn[2][1][2] = r2/rho2 ;
	lconn[2][1][3] = a*r1_rho23*cth*sth*(r3*(2.+r1) + a2*(2.*r1r1p1*cth2 + a2*cth4 + 2.*r1*sth2))/(M_PI*fac4) ;

	lconn[2][2][0] = lconn[2][0][2] ;
	lconn[2][2][1] = lconn[2][1][2] ;
	lconn[2][2][2] = -4.*M_PI*(th - M_PI*x2)/fac4 - a2*M_PI*fac4*s2th/(2.*rho2) ;
	lconn[2][2][3] = 0. ;

	lconn[2][3][0] = lconn[2][0][3] ;
	lconn[2][3][1] = lconn[2][1][3] ;
	lconn[2][3][2] = lconn[2][2][3] ;
	lconn[2][3][3] = -(cth*sth/(M_PI*fac4))*(1. + 2.*a4*r1_rho23*sth4 + 
			sth2*(4.*a2*r1 + a2*r2 + a4*cth2)/rho22) ;

	lconn[3][0][0] = a*fac1_rho23 ;
	lconn[3][0][1] = a*r1fac1_rho23 ;
	lconn[3][0][2] = -2.*r1api*fac4*cth/(sth*rho22) ;
	lconn[3][0][3] = -a2*fac1_rho23*sth2 ;

	lconn[3][1][0] = lconn[3][0][1] ;
	lconn[3][1][1] = a*r2*fac1_rho23 ;
	lconn[3][1][2] = -2.*r1api*(a2 + 2.*r1*(2.+r1) + a2*c2th)*fac4*cth/(sth*fac52) ;
	lconn[3][1][3] = r1_rho23*(r1*rho22 - a2*fac1*sth2) ;

	lconn[3][2][0] = lconn[3][0][2] ;
	lconn[3][2][1] = lconn[3][1][2] ;
	lconn[3][2][2] = -r1api*M_PI*fac42_rho2 ;
	lconn[3][2][3] = M_PI*(3.*a4 + 8.*r4 + 8.*a2*r1r1p1 +
			4.*a2*(a2 + 2.*r1*(r1-1.))*c2th + a4*c4th)*fac4*cth/(2.*sth*fac52) ;

	lconn[3][3][0] = lconn[3][0][3] ;
	lconn[3][3][1] = lconn[3][1][3] ;
	lconn[3][3][2] = lconn[3][2][3] ;
	lconn[3][3][3] = -asth2*r1_rho23*rho22 + a3*fac1_rho23*sth4 ;

}

#define EPS	0.03

double stepsize(double X[NDIM], double Kcon[NDIM])
{
	double dl, dlx1, dlx2, dlx3;
	double idlx1,idlx2,idlx3 ;

	dlx1 = EPS / (fabs(Kcon[1]) + SMALL*SMALL) ;
	dlx2 = EPS * GSL_MIN(X[2], 1. - X[2]) / (fabs(Kcon[2]) + SMALL*SMALL) ;
	dlx3 = EPS / (fabs(Kcon[3]) + SMALL*SMALL) ;

	idlx1 = 1./(fabs(dlx1) + SMALL*SMALL) ;
	idlx2 = 1./(fabs(dlx2) + SMALL*SMALL) ;
	idlx3 = 1./(fabs(dlx3) + SMALL*SMALL) ;

	dl = 1. / (idlx1 + idlx2 + idlx3) ;

	return (dl);
}

void init_model(char *argv[])
{
	char filename[512] ;

	/* Reads in model parameters; put frequently changed
	   parameters into command line */
    	set_model_param(); 

	sscanf(argv[3],"%s",filename) ;

        fprintf(stderr, "getting simulation data...\n");
        init_harm_data(filename);  /* read in HARM simulation data */
        fprintf(stderr, "done.\n\n");

        Rh = 1. + sqrt(1. - a * a) ;

        /* find dimensional quantities from black hole
                mass and its accretion rate */
        set_units(argv[4]);

        /* pre-compute densities, field strengths, etc. */
        init_physical_quantities() ;
}

double get_model_thetae(double X[NDIM])
{
        if(X[1] < startx[1] || 
           X[1] > stopx[1]  || 
           X[2] < startx[2] || 
           X[2] > stopx[2]) {
                return(0.) ;
        }
        
        return(interp_scalar(X, thetae)) ;
}

double get_model_b(double X[NDIM])
{
        if(X[1] < startx[1] || 
           X[1] > stopx[1]  || 
           X[2] < startx[2] || 
           X[2] > stopx[2]) {
                return(0.) ;
        }

        return(interp_scalar(X, b)) ;
}

double get_model_ne(double X[NDIM])
{
        if(X[1] < startx[1] || 
           X[1] > stopx[1]  || 
           X[2] < startx[2] || 
           X[2] > stopx[2]) {
                return(0.) ;
        }

        return(interp_scalar(X, ne)) ;
}

void get_model_bcov(double X[NDIM], double Bcov[NDIM])
{
        if(X[1] < startx[1] || 
           X[1] > stopx[1]  || 
           X[2] < startx[2] || 
           X[2] > stopx[2]) {

                Bcov[0] = 0. ;
                Bcov[1] = 0. ;
                Bcov[2] = 0. ;
                Bcov[3] = 0. ;

                return ;
        }
        interp_fourv(X, bcov, Bcov) ;
}

void get_model_bcon(double X[NDIM], double Bcon[NDIM])
{
        if(X[1] < startx[1] || 
           X[1] > stopx[1]  || 
           X[2] < startx[2] || 
           X[2] > stopx[2]) {

                Bcon[0] = 0. ;
                Bcon[1] = 0. ;
                Bcon[2] = 0. ;
                Bcon[3] = 0. ;

                return ;
        }
        interp_fourv(X, bcon, Bcon) ;
}

void get_model_ucon(double X[NDIM], double Ucon[NDIM])
{

        double gcov[NDIM][NDIM] ;
        double gcon[NDIM][NDIM] ;
        double tmp[NDIM] ;

        if(X[1] < startx[1] || 
           X[1] > stopx[1]  || 
           X[2] < startx[2] || 
           X[2] > stopx[2]) {
                /* sensible default value */
                gcov_func(X, gcov) ;

                tmp[0] = -1./sqrt(-gcov[0][0]) ;
                tmp[1] = 0. ;
                tmp[2] = 0. ;
                tmp[3] = 0. ;

                gcon_func(gcov, gcon) ;
                Ucon[0] = 
                        tmp[0]*gcon[0][0] +
                        tmp[1]*gcon[0][1] +
                        tmp[2]*gcon[0][2] +
                        tmp[3]*gcon[0][3] ;
                Ucon[1] = 
                        tmp[0]*gcon[1][0] +
                        tmp[1]*gcon[1][1] +
                        tmp[2]*gcon[1][2] +
                        tmp[3]*gcon[1][3] ;
                Ucon[2] = 
                        tmp[0]*gcon[2][0] +
                        tmp[1]*gcon[2][1] +
                        tmp[2]*gcon[2][2] +
                        tmp[3]*gcon[2][3] ;
                Ucon[3] = 
                        tmp[0]*gcon[3][0] +
                        tmp[1]*gcon[3][1] +
                        tmp[2]*gcon[3][2] +
                        tmp[3]*gcon[3][3] ;
        
                return ;
        }
           
        interp_fourv(X, ucon, Ucon) ;
}

void get_fluid_params(double X[NDIM], double *Ne, double *Thetae, double *B, double Ucov[NDIM], double Bcov[NDIM])
{

	/* super-slow way to do this.  Could be rewritten */
	*Ne = get_model_ne(X) ;
	*Thetae = get_model_thetae(X) ;
	*B = get_model_b(X) ;

	get_model_ucov(X,Ucov) ;
	get_model_bcov(X,Bcov) ;
}

/* 

	these supply basic model data to grmonty

*/

void get_model_ucov(double X[NDIM], double Ucov[NDIM])
{
	double gcov[NDIM][NDIM], Ucon[NDIM] ;

	gcov_func(X, gcov);

	if(X[1] < startx[1] || 
	   X[1] > stopx[1]  || 
	   X[2] < startx[2] || 
	   X[2] > stopx[2]) {
	   
	   	/* sensible default value */
		Ucov[0] = -1./sqrt(-gcov[0][0]) ;
		Ucov[1] = 0. ;
		Ucov[2] = 0. ;
		Ucov[3] = 0. ;

		return ;
	}

	get_model_ucon(X, Ucon);
	lower(Ucon, gcov, Ucov);

}

