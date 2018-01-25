/*
 * perform least-square fit of a superellipse
 *
 * use Powell's method for non-linear optimisation
 *
 * Paul Rosin
 * April 1993
 *
 * =========================================================
 *
 * added extra overflow check for DEC Alpha
 * Paul Rosin
 * October 1994
 *
 * added alternatives: LS or LMedS error
 * Paul Rosin
 * January 1995
 *
 * =========================================================
 *
 * added alternative distance measures
 * also tidied up func()
 * Paul Rosin
 * March 1998
 *
 * =========================================================
 *
 * added bias to avoid star shapes
 * uses ellipse_specific fit for initialisation
 * Paul Rosin
 * July 2005
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "nr.h"
#include "nrutil.h"

#ifndef FALSE
# define FALSE 0
# define TRUE (!FALSE)
#endif

#define ALGEBRAIC   1
#define GRADIENT    2
#define DIRECTIONAL 3
#define GROSS       4
#define RAY         5
#define SAFAEE      6
#define EXPAND      7
#define CONSTANT    8
#define QUADRARC    9
#define RECTANGLE   10

#define NO_PIXELS 10000

#define LMEDS FALSE

#define TINY 1.0e-20
#define NDIM 6
#define FTOL 1.0e-6

#define PI 3.1415926

#define SQR(x)   ((x)*(x))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define ABS(a)   ((a)>=0 ? a : -(a))
#define R2D(x)   ((x)*180.0/PI)
#define D2R(x)   ((x)*PI/180.0)

extern float x_trans3[],y_trans3[];
extern int nseg2;
extern int representation_ok;
extern int measure;

extern int bias;
extern double mean_error;

fit_superellipse(a,b,e,xc,yc,rot)
float *a,*b,*e,*xc,*yc,*rot;
{
    double aT,bT,eT,xcT,ycT,rotT;
    float radius2;
    float compute_orientation();

    /* fit ellipse to get initial estimate for superellipse */
    ellipse_fit(&xcT,&ycT,&aT,&bT,&rotT);
    //printf("initialisation:  ctr: %.2f %.2f  ab: %.2f %.2f  r: %.2f\n",xcT,ycT,aT,bT,rotT);

    rotT = R2D(rotT);
    eT = 1;
    testfit(&aT,&bT,&eT,&xcT,&ycT,&rotT);
    rotT = D2R(rotT);

    *a = aT;
    *b = bT;
    *e = eT;
    *xc = xcT;
    *yc = ycT;
    *rot = rotT;
}

/* avoid error messages from math library */
/* HAD TO ADD EXTRA OVERFLOW CHECK FOR DEC - PLR Oct 1994 */
float POW(a,b)
float a,b;
{
    float r;
    float LN_MAX_FLOAT = 88;
    float LN_MIN_FLOAT = -88;

    if (a == 0)
        r = 0;
    else {
        /* trap overflow */
        if (b*log((float)a) > LN_MAX_FLOAT)
            r = 9e9;
        /* trap underflow */
        else if (b*log((float)a) < LN_MIN_FLOAT)
            r = 0;
        else
            r = pow((double)a,(double)b);
    }

    return(r);
}

float pow2(a,b)
float a,b;
{
    float r;

    if (a < 0) a = -a;

    if (a == 0)
        r = 0;
    else
        r = pow((double)a,(double)b);

    r = ABS(r);
    return(r);
}

/*********************************************
 definition of function to optimize
**********************************************/
float func(p)
float p[];
{
    int i;
    float fi,fu;
    double xi,yi;
    float dx,dy;
    float tx,ty;
    float x,y;
    float sine2,cosine2;
    float list[NO_PIXELS];
    float maj_axis = p[1];
    float min_axis = p[2];
    float e = p[3];

    /* abort if a, b, or e have impossible values */
    /* ### set lower bound on e to larger value = 0.1; helps convergence ###*/
    if ((maj_axis <= 0) || (min_axis <= 0) || (e <= 0.1)) {
        return(9e9);
    }

    cosine2 = cos(D2R(-p[6]));
    sine2 = sin(D2R(-p[6]));

    fu = 0;
    for (i = 1; i <= nseg2; i++) {
        tx = x_trans3[i] - p[4];
        ty = y_trans3[i] - p[5];

        x = ABS(tx * cosine2 - ty * sine2);
        y = ABS(tx * sine2 + ty * cosine2);

        /* skip points on axes - too much trouble with divide by zero errors */
        if ((x == 0) || (y == 0))
            continue;

        /* algebraic distance */
        if (measure == ALGEBRAIC) {
            fi = pow2(x/maj_axis,2.0/e) + pow2(y/min_axis,2.0/e) - 1;
        }

        /* Gross & Boult */
        else if (measure == GROSS) {
            fi = pow2(x/maj_axis,2.0/e) + pow2(y/min_axis,2.0/e) - 1;
            fi = ABS(fi);
            fi = pow2(fi,e);
        }

        /* algebraic distance divided by gradient */
        else if (measure == GRADIENT) {
            double gx,gy,grad;

            gx = 2.0 / (e * x) * pow2(x/maj_axis,2.0/e);
            gy = 2.0 / (e * y) * pow2(y/min_axis,2.0/e);
            grad = sqrt(SQR(gx) + SQR(gy));
            fi = pow2(x/maj_axis,2.0/e) + pow2(y/min_axis,2.0/e) - 1;
            fi /= grad;
        }

        /* weight algebraic distance by directional derivative */
        else if (measure == DIRECTIONAL) {
            double gx,gy,theta,deriv;

            theta = atan2((double)y,(double)x);
            gx = 2.0 / (e * x) * pow2(x/maj_axis,2.0/e);
            gy = 2.0 / (e * y) * pow2(y/min_axis,2.0/e);
            deriv = cos(theta) * gx + sin(theta) * gy;
            fi = pow2(x/maj_axis,2.0/e) + pow2(y/min_axis,2.0/e) - 1;
            fi /= deriv;
        }

        /* expand the axes by a factor */
        else if (measure == EXPAND) {
            double mean,da,db,a_new,b_new;
            double s,te,t1,t2;

            te = 2.0 / e;
            t1 = pow2((double)x/(double)maj_axis,te);
            t2 = pow2((double)y/(double)min_axis,te);
            s  = pow2(t1+t2,e/2.0);
            a_new = maj_axis * s;
            b_new = min_axis * s;

            fi = maj_axis - a_new;
        }

        /* expand the axes by adding a constant to their length
         */
        else if (measure == CONSTANT) {
            double s1,s2;
            double A,B,C;
            double te,ae,be,ae1,be1,xe,ye;
            double a,b,t;

            te = 2.0 / e;
            ae = pow2((double)maj_axis,te);
            be = pow2((double)min_axis,te);
            ae1 = pow2((double)maj_axis,te-1.0);
            be1 = pow2((double)min_axis,te-1.0);
            xe = pow2((double)x,te);
            ye = pow2((double)y,te);

            A = 4*ae1*be1/SQR(e);
            B = 2*(ae1*be + ae*be1 - xe*be1 - ye*ae1)/e;
            C = ae*be - xe*be - ye*ae;

            s1 = (-B + sqrt(SQR(B) - 4*A*C)) / (2*A);
            s2 = (-B - sqrt(SQR(B) - 4*A*C)) / (2*A);

            fi = pow2(s1,e/2.0);
        }

        /* Knowlton's quadrarc */
        else if (measure == QUADRARC) {
            double xi,yi,t,h,k;

            xi = maj_axis/pow2(2.0,e/2.0);
            yi = min_axis/pow2(2.0,e/2.0);

            h = SQR(maj_axis)-SQR(xi)-SQR(yi);
            h /= 2*(maj_axis-xi);
            k = SQR(min_axis)-SQR(xi)-SQR(yi);
            k /= 2*(min_axis-yi);
            k = -k;

            t = y*(maj_axis-min_axis) +
                x*( (min_axis*yi/(xi-h)) - (maj_axis*(k+yi)/xi) ) +
                (maj_axis*k-min_axis*h);

            if (t < 0) {
                fi = sqrt(SQR(h-x)+SQR(y))-(maj_axis-h);
            }
            else {
                fi = sqrt(SQR(x)+SQR(y+k))-(min_axis+k);
            }
        }

        /* combine ellipse & rectangle distance */
        else if (measure == RECTANGLE) {
            double Aq,Bq,Cq;
            double A,A1,A2,F,X,Y;
            double ae,be,ah,bh; /* ellipse & hyperbola axes */
            double xi,yi,euc1,euc2,euc3,euc4;
            double edist,rdist;

            /* ####### ellipse distance ####### */

            ae = MAX(maj_axis,min_axis);
            be = MIN(maj_axis,min_axis);

            /* solve to find confocal hyperbola through point (x,y) */

            X = SQR(x);
            Y = SQR(y);
            F = SQR(ae) - SQR(be);

            Aq = 1;
            Bq = -(X + Y + F);
            Cq = X * F;
            A1 = (-Bq - sqrt(SQR(Bq) - 4 * Aq * Cq)) / (2 * Aq);
            A2 = (-Bq + sqrt(SQR(Bq) - 4 * Aq * Cq)) / (2 * Aq);

            if ((A1 < 0) && (A2 < 0)) {
                fprintf(stderr,"ERROR: only negative solutions to ah\n");
                exit(-1);
            }

            if (A1 < 0)
                A = A2;
            else if (A2 < 0)
                A = A1;
            else {
                if (A1 < F)
                    A = A1;
                else
                    A = A2;
            }

            if ((F - A) < 0) {
                /*
                fprintf(stderr,"ERROR: only negative solutions to bh\n");
                exit(-1);
                */
                ah = bh = 0;
            }
            else {
                ah = sqrt(A);
                bh = sqrt(ABS(F - A));
            }

            /* find intersection of ellipse & hyperbola */

            xi = ah * sqrt(
                        (SQR(ae) * (SQR(be) + SQR(bh))) /
                        (SQR(ah) * SQR(be) + SQR(ae) * SQR(bh))
                      );
            yi = (be * bh * sqrt(SQR(ae) - SQR(ah))) /
                  sqrt(SQR(ah) * SQR(be) + SQR(ae) * SQR(bh));

            euc1 = SQR(x-xi) + SQR(y-yi);
            euc2 = SQR(x+xi) + SQR(y-yi);
            euc3 = SQR(x-xi) + SQR(y+yi);
            euc4 = SQR(x+xi) + SQR(y+yi);
            edist = sqrt(
                        MIN(
                            MIN(euc1,euc2),
                            MIN(euc3,euc4)
                        )
                   );

            /* ####### rectangle distance ####### */

            if ((x > ae) && (y > be))
                rdist = sqrt(SQR(x-ae)+SQR(y-be));
            else if ((x < ae) && (y < be))
                rdist = MIN(ABS(y-be),ABS(x-ae));
            else if (x <= ae)
                rdist = ABS(y-be);
            else if (y <= be)
                rdist = ABS(x-ae);

            /* ####### combined ellipse & rectangle distance ####### */

            fi = edist*e + rdist*(1.0-e);
        }

        /* ray */
        else if (measure == RAY) {
            xi = pow2(p[1],2.0/p[3]);
            if (xi != 0)
                xi = 1.0 / xi + pow2(y/(x*p[2]),2.0/p[3]);
            else
                xi = 1e10;

            xi = 1.0 / xi;
            xi = ABS(xi);
            xi = pow2(xi,p[3]/2.0);

            yi = xi * y/x;

            dx = x - xi;
            dy = y - yi;

            fi = sqrt(dx*dx + dy*dy);
        }

        /* Safaee-Rad */
        else if (measure == SAFAEE) {
            double centre_to_point,centre_to_ellipse,ellipse_to_point;
            double w,m,t1,t2,t3,xi,yi;
            double te = 2.0 / e;

            /* compute intersection of ray & ellipse from centre to (x,y) */
            if (x == 0) {
                xi = 0;
                yi = min_axis;
            }
            else {
                m = (double)y / (double)x;
                t1 = pow2(1.0/(double)maj_axis, te);
                t2 = pow2(m/(double)min_axis, te);
                t3 = 1.0 / (t1 + t2);
                xi = pow2(t3,e/2.0);
                yi = m * xi;
            }
            /* distance from the centre to the point */
            centre_to_point = sqrt((double)SQR(x) + (double)SQR(y));

            /* distance from centre to the ellipse */
            centre_to_ellipse = sqrt(SQR(xi) + SQR(yi));

            /* difference between the two */
            ellipse_to_point = ABS(centre_to_point - centre_to_ellipse);

            /* full Safaee-Rad */
            w = centre_to_ellipse;
            w = w * (1 + ellipse_to_point / (2.0 * maj_axis));
            w = w / (1 + ellipse_to_point / (2.0 * centre_to_ellipse));

            fi = pow2(x/maj_axis,2.0/e) + pow2(y/min_axis,2.0/e) - 1;
            fi *= w;
        }

        else {
            fprintf(stderr,"ERROR: %d unknown distance measure\n",measure);
            exit(-1);
        }

        /* bias the superellipse to have 0 <= e <= 1 */
        if (bias)
            if (p[3] > 1)
                //fi *= 10 * SQR(p[3]);
                fi *= SQR(p[3]);

#if LMEDS
        list[i] = ABS(fi);
#else
        fu += ABS(fi);
#endif
    }

#if LMEDS
    sort(nseg2,list);
    return(list[nseg2/2]);
#else
    return(fu);
#endif
}

/**********************************************
 Fitting data to an ellipse
***********************************************/
testfit(a,b,e,xo,yo,rot)
double *a,*b,*e,*xo,*yo,*rot;
{
    int i,iter,j;
    float fret,**xi;
    float p[7];

    /* starting point */
    p[1] = *a;
    p[2] = *b;
    p[3] = *e;
    p[4] = *xo;
    p[5] = *yo;
    p[6] = *rot;

    /* set of initial directions */
    xi = matrix(1,NDIM,1,NDIM);
    for (i = 1; i <= NDIM; i++)
        for (j = 1; j <= NDIM; j++)
            xi[i][j] = (i == j ? 1.0 : 0.0);

    powell(p,xi,NDIM,FTOL,&iter,&fret,func);

    /***
    fprintf(stderr,"Iterations: %d\n",iter);
    fprintf(stderr,"Minimum found at:\n");
    for (i = 1; i <= NDIM; i++)
        fprintf(stderr,"%.2f  ",p[i]);
    ***/

    *a = p[1];
    *b = p[2];
    *e = p[3];
    *xo = p[4];
    *yo = p[5];
    *rot = p[6];

    /***
    fprintf(stderr,"\nMinimum function value = %.2f\n\n",fret);
    ***/
    mean_error = fret / nseg2;

    free_matrix(xi,1,NDIM,1,NDIM);
}

/* from numerical recipes in C */
sort(n,ra)
int n;
float ra[];
{
   int l,j,ir,i;
   float rra;

   l=(n >> 1)+1;
   ir=n;
   for (;;) {
      if (l > 1)
         rra=ra[--l];
      else {
         rra=ra[ir];
         ra[ir]=ra[1];
         if (--ir == 1) {
            ra[1]=rra;
            return;
         }
      }
      i=l;
      j=l << 1;
      while (j <= ir) {
         if (j < ir && ra[j] < ra[j+1]) ++j;
         if (rra < ra[j]) {
            ra[i]=ra[j];
            j += (i=j);
         }
         else j=ir+1;
      }
      ra[i]=rra;
   }
}
