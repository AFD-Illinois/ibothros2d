
#include "decs.h"
#include "defs.h"

#define MAXNSTEP 	10000
#define NSUBSTEP 	5

struct of_traj {
	double j ;
	double k ;
	double dl ;
	double X[NDIM] ;
	double Kcon[NDIM] ;
} traj[MAXNSTEP] ;

int main(int argc, char *argv[]) 
{
	double X[NDIM],Kcon[NDIM] ;
	double dl,I,dlsubstep ;
	double DX,DY,fovx,fovy ;
	double x2cam,phicam,rcam,Xcam[NDIM],Ucam[NDIM] ;
	double Gcov[NDIM][NDIM] ;
	double image[NX][NY] ;
	double freq,freqcgs ;
	double Ftot ;
	int i,j,k,l,nstep ;
	double Xi[NDIM],Xf[NDIM],Kconi[NDIM],Kconf[NDIM],ki,kf,si,sf ;

	if(argc < 3) {
		fprintf(stderr,"usage: ibothros2d x2cam freq filename\n") ;
		exit(0) ;
	}

	sscanf(argv[1],"%lf",&x2cam) ;
	sscanf(argv[2],"%lf",&freqcgs) ;
    
	init_model(argv) ;

	/* normalize frequency to electron rest-mass energy */
	freq = freqcgs*HPL/(ME*CL*CL) ;

	/* initialize local parameters */
	/* fix camera worldline */
	rcam = 1.e3 ;
	phicam = 0. ;
	Xcam[0] = 0. ;
	Xcam[1] = log(rcam) ;
	Xcam[2] = x2cam;
	Xcam[3] = phicam ;

	Ucam[0] = 1. ;
	Ucam[1] = 0. ;
	Ucam[2] = 0. ;
	Ucam[3] = 0. ;
	gcov_func(Xcam,Gcov) ;
	normalize(Ucam,Gcov) ;

	/* fix camera field of view */
	DX = 50. ;	/* size of field of view in the plane of the black hole, in units of GM/c^2 */
	DY = 50. ;
	fovx = DX/rcam ;
	fovy = DY/rcam ;

	/* loop over pixels */
	for(i=0;i<NX;i++) {
		fprintf(stderr,"%d ",i) ;
	for(j=0;j<NY;j++) {

		/* initialize wavevectors, positions */
		init(i,j,Xcam,Ucam,fovx,fovy,X,Kcon) ;

		/* normalize Kcon to desired frequency */
		for(k=0;k<NDIM;k++){
			Kcon[k] *= freq ;
		}

		/* integrate backwards along trajectory */
		nstep = 0 ;
		while(!stop_backward_integration(X,Kcon,Xcam)) {

			/* get emissivity at current location */
			/* This stepsize function can be troublesome inside of R = 2M,
			   and should be used cautiously in this region. */
			dl = stepsize(X,Kcon) ;

			/* record data */
			traj[nstep].dl = dl*L_unit /( ME*CL*CL/HPL ) ;
			for(k=0;k<NDIM;k++) traj[nstep].X[k] = X[k] ;
			for(k=0;k<NDIM;k++) traj[nstep].Kcon[k] = Kcon[k] ;

			/* move photon */
			push_photon(X,Kcon,-dl) ;

			nstep++ ;
			if(nstep > MAXNSTEP-2) {
				fprintf(stderr,"MAXNSTEP exceeded on j=%d, %d\n",j,nstep) ;
				break ;
			}
		}
		nstep-- ;	/* final step violated the "stop" condition,
				   so don't record it */

		/* integrate forwards along trajectory, including
		   radiative transfer equation */
		I = 0. ;
		while(nstep > 0) {

			/* find start point emissivity */
			for(l=0;l<NDIM;l++) {
				Xi[l] = traj[nstep].X[l] ;
				Kconi[l] = traj[nstep].Kcon[l] ;
			}
			//get_jkinv(Xi, Kconi,&ji,&ki) ;
			get_skinv(Xi, Kconi, &si, &ki);

			/* loop over substeps.  stepsize really should
			   be set adaptively */
			dlsubstep = traj[nstep-1].dl/NSUBSTEP ;
			for(k=0;k<NSUBSTEP;k++) {

				/* end point */
				for(l=0;l<NDIM;l++) {
					Xf[l] = 
					(1. - (((double) (k+1))/NSUBSTEP)) * traj[nstep].X[l] +
					(((double)(k+1))/NSUBSTEP) * traj[nstep-1].X[l] ;
					Kconf[l] = 
					(1. - (((double)(k+1))/NSUBSTEP)) * traj[nstep].Kcon[l] +
					(((double) (k+1))/NSUBSTEP) * traj[nstep-1].Kcon[l] ;
				}
				//get_jkinv(Xf,Kconf,&jf,&kf) ;
				get_skinv(Xf, Kconf, &sf, &kf);

				//I = approximate_solve(I,ji,ki,jf,kf,dlsubstep) ;
				I = approximate_solve(I, 0.5*(si+sf), 0.5*dlsubstep*(ki+kf));

				/* swap start and finish */
				//ji = jf ;
				si = sf;
				ki = kf ;
			}
			nstep-- ;
		}

		/* deposit intensity in pixel */
		image[i][j] = I*pow(freqcgs,3) ;
	}
	}

	Ftot = 0. ;
	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) Ftot += image[i][j] ;

	Dsource *= PC ;
	Ftot *= (DX*L_unit/NX)*(DY*L_unit/NY)/(Dsource*Dsource) ;

	Ftot = Ftot/JY ;
	double nLn ;
	nLn = Ftot * JY * freqcgs * 4.*M_PI*Dsource*Dsource ;
	fprintf(stderr,"\nfreq, Ftot, nuLnu: %g %g %g\n",freqcgs,Ftot,nLn) ;


	/* image, dump result */
	make_ppm(image, freq, "ibothros2d_fnu.ppm") ;
	dump(image, "ibothros2d.dat") ;
	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) image[i][j] = log(image[i][j] + 1.e-50) ;
	make_ppm(image, freq, "ibothros2d_lfnu.ppm") ;

	/* done! */
	return(0) ;
}

void dump(double image[NX][NY], char *fname) 
{
	FILE *fp ;
	int i,j ;

	fp = fopen(fname,"w") ;
	if(fp == NULL) {
		fprintf(stderr,"unable to open %s\n",fname) ;
		exit(1) ;
	}

	for(i=0;i<NX;i++)
	for(j=0;j<NY;j++) fprintf(fp,"%d %d %15.10g\n",i,j,image[i][j]) ;

	fclose(fp) ;

}

void init(int i, int j, 
	double Xcam[4], double Ucam[4], 
	double fovx, double fovy,	/* field of view, in radians */
	double X[4], double Kcon[4]	/* position, wavevector */
) 
{
	double Gcov[NDIM][NDIM] ;
	double Econ[NDIM][NDIM] ;
	double Ecov[NDIM][NDIM] ;
	double Kcon_tetrad[NDIM] ;
	double trial1[NDIM] ;
	double trial2[NDIM] ;
	int k ;

	/* construct orthonormal tetrad.
		e^0 along Ucam
		e^1 inward along radius vector
		e^2 toward north pole of coordinate system
			("y" for the image plane)
		e^3 in the remaining direction
			("x" for the image plane)
		note: this is *modified* from the make_tetrad 
			routine used in grmonty

		this could easily be modified for a fly-through.
	*/
	/* set up trial vectors */
	trial1[0] = 0. ; trial1[1] = -1. ; trial1[2] = 0. ; trial1[3] = 0. ;
	trial2[0] = 0. ; trial2[1] = 0. ; trial2[2] = -1. ; trial2[3] = 0. ;
	gcov_func(Xcam, Gcov) ;
	make_tetrad(Ucam, trial1, trial2, Gcov, Econ, Ecov) ;

	/* construct *outgoing* wavevectors */
	Kcon_tetrad[0] = 0. ;
	Kcon_tetrad[1] = -1. ;
	Kcon_tetrad[2] = -(j/((double)NY) - 0.5)*fovy ;
	Kcon_tetrad[3] = -(i/((double)NX) - 0.5)*fovx ;

	/* normalize */
	null_normalize(Kcon_tetrad, 1.) ;

	/* translate into coordinate frame */
	tetrad_to_coordinate(Econ, Kcon_tetrad, Kcon) ;
	
	/* set position */
	for(k=0;k<NDIM;k++) X[k] = Xcam[k] ;

	/* done! */
}

/* normalize null vector in a tetrad frame */
void null_normalize(double Kcon[NDIM], double fnorm) 
{
	double inorm ;

	inorm = sqrt( Kcon[1]*Kcon[1] + Kcon[2]*Kcon[2] + Kcon[3]*Kcon[3]) ;

	Kcon[0] = fnorm ;
	Kcon[1] *= fnorm/inorm ;
	Kcon[2] *= fnorm/inorm ;
	Kcon[3] *= fnorm/inorm ;

}

/* 
   must be a stable, approximate solution to radiative transfer
   that runs between points w/ initial intensity I, emissivity
   ji, opacity ki, and ends with emissivity jf, opacity kf.

   Return final intensity
*/
/*
double approximate_solve(double Ii,
	double ji,
	double ki,
	double jf,
	double kf,
	double dl)
{
	double efac,If,javg,kavg,dtau ;

	javg = (ji + jf)/2. ;
	kavg = (ki + kf)/2. ;

	dtau = dl*kavg ;

	if(dtau < 1.e-5) {
		If = Ii + (javg - Ii*kavg)*dl*(1. - (dtau/2.)*(1. - dtau/3.)) ;
	}
	else {
		efac = exp(-dtau) ;
		If = Ii*efac + (javg/kavg) * (1. - efac) ;
	}


	return(If) ;
}
*/
double approximate_solve(double Ii, double S, double dtau)
{
	double If,efac;

	if(dtau < 1.e-5)
		If = Ii - (Ii - S) * ( 0.166666667*dtau * (6. - dtau * (3. - dtau)));
	else{
		efac = exp(-dtau);
		If = Ii*efac + S*(1. - efac);
	}

	return If;
}

/* condition for stopping the backwards-in-lambda
   integration of the photon geodesic */

#define LRMAX (log(1.1*Rout))
#define LRMIN (log(1.05*Rh))

int stop_backward_integration(
	double X[NDIM],
	double Kcon[NDIM],
	double Xcam[NDIM])
{

	if(
		(X[1] > LRMAX && Kcon[1] < 0.) || 	/* out far */
		 X[1] < LRMIN				/* in deep */
	) return(1) ;
	else return(0) ;                                /* neither out far nor in deep */

}

int get_skinv(double X[NDIM], double Kcon[NDIM], double *snuinv, double *knuinv)
{
	double nu,theta,B,Thetae,Ne,Bnuinv,Ucov[NDIM],Bcov[NDIM],jnuinv ;

	/* get fluid parameters */
	get_fluid_params(X, &Ne, &Thetae, &B, Ucov, Bcov);

	if(Ne == 0.) {
		*snuinv = 0. ;
		*knuinv = 0. ;
		return 0;
	}

	/* get covariant four-velocity of fluid for use in get_bk_angle and get_fluid_nu */

	theta = get_bk_angle(X,Kcon,Ucov,Bcov,B) ;	/* angle between k & b */
	if(theta <= 0. || theta >= M_PI) {	/* no emission along field */
		*snuinv = 0. ;
		*knuinv = 0. ;
		return 0;
	}

	nu = get_fluid_nu(Kcon,Ucov);


	/* assume emission is thermal */
	Bnuinv = Bnu_inv(nu,Thetae) ;

	jnuinv = jnu_inv(nu,Thetae,Ne,B,theta) ;

	if(Bnuinv < SMALL)
		*knuinv = 0.;
	else
		*knuinv = jnuinv/Bnuinv ;

	*snuinv = Bnuinv;

	/* check for errors */
	if(isnan(snuinv[0]) || isnan(knuinv[0])) {
		fprintf(stderr,"isnan get_jkinv\n") ;
		fprintf(stderr,">> %g %g %g %g %g %g %g %g\n",snuinv[0],knuinv[0],
			Ne,theta,nu,B,Thetae,Bnuinv) ;
	}

	return 0;
}

