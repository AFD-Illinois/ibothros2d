
/*

make_ppm

modified 6 June 2013 CFG
	- adapted this from iHARM2d image.c routine for ibothros2d.
	- simpler, shorter!

modified 17 June 2012 CFG
        - eliminate r8 option (never used).
        - eliminate read-in of color map (map.ppm no longer required)
        - insert analytic function for john.pal color map
        - removed failimage capability (probably unwise, but eliminates
                global arrays for failimage) 

*/

#include "decs.h"
#include <ctype.h>

void make_ppm(double p[NX][NY], double freq, char *fname)
{
	int i,j,k ;
	double max, min ;
	static double q[NX * NY];
	static double f[NX * NY];
        int red, green, blue ;
	int compare_doubles(const void *a, const void *b);
        void john_pal(double data, double min, double max, int *pRed, int *pGreen, int *pBlue) ;
	FILE *fp;

	if ((fp = fopen(fname, "w")) == NULL) {
		fflush(stderr);
		fprintf(stderr, "image(): Cannot open %s !! \n", fname);
		fflush(stderr);
		exit(1);
	}

	/* put the data into a 1D array */
	k = 0 ;
        for (j = NY-1; j >= 0; j--) 
        for (i = 0; i < NX; i++) {
		f[k] = p[i][j] ;
		k++ ;
	}

	/*  mapping is in 255 steps lmax and lmin */
	for (i = 0; i < NX * NY; i++) q[i] = f[i];
                // this bit sorts the values and lops off the top and bottom
                // end of the distribution to create max and min.
	qsort(q, NX * NY, sizeof(double), compare_doubles);
	min = q[NX * NY / 128];
	max = q[NX * NY - NY * NY / 128];

        /* write out header information */
        fprintf(fp, "P6\n#  min=%g  , max=%g \n#  frequency=%g \n%d %d\n%d\n",
                min, max, freq, NX, NY, 255);
        fflush(fp);

	for (i = 0; i < NX * NY; i++) {
                john_pal(f[i], min, max ,&red,&green,&blue) ;
		fputc((char) red, fp);
		fputc((char) green, fp);
		fputc((char) blue, fp);
	}

	fclose(fp);

	return;
}

/* this is needed for qsort, which is used to set colormap min and max */
int compare_doubles(const void *a, const void *b)
{
	const double *da = (const double *) a;
	const double *db = (const double *) b;

	return (*da > *db) - (*da < *db);
}

/* 
color palette, based on john.pal

        input: integer 0-255 
        output: red, green, blue integers 

        author: Bryan M. Johnson
*/

void john_pal(double data, double min, double max, int *pRed, int *pGreen, int *pBlue)
{
  double a, b, c, d, e, f;
  double x, y;
  double max_min = max - min ;

  if(max_min > 0.0) { // trust no one
    x = (data - min)/(max_min) ;
    
    /* ========== Red ============ */
    a = 4.0*x - 1.52549019607844;
    b = 4.52941176470589 - 4.0*x;
    y = a < b ? a : b;
    *pRed = (int)(255.0*y);
    *pRed = *pRed >   0 ? *pRed :   0;
    *pRed = *pRed < 255 ? *pRed : 255;

    /* ========== Green ========== */
    a = 4.0*x - 0.521568627450979;
    b = 2.52549019607844 - 4.0*x;
    c = a < b ? a : b;
    d = 4.0*x - 1.53725490196073;
    e = 3.52941176470581 - 4.0*x;
    f = d < e ? d : e;
    y = c > f ? c : f;
    *pGreen = (int)(255.0*y);
    *pGreen = *pGreen >   0 ? *pGreen :   0;
    *pGreen = *pGreen < 255 ? *pGreen : 255;

    /* ========== Blue =========== */
    a = 4.0*x + 0.498039215686276;
    b = 2.50980392156862 - 4.0*x;
    y = a < b ? a : b;
    *pBlue = (int)(255.0*y);
    *pBlue = *pBlue >   0 ? *pBlue :   0;
    *pBlue = *pBlue < 255 ? *pBlue : 255;

  }
  else {
    *pRed = *pGreen = *pBlue = (data > max ? 255: 0) ;
  }

  return;
}

