/*
 * perform least-square fit of an ellipse
 *
 * Paul Rosin
 * April 1993
 */

#include <stdio.h>
#include <math.h>

#ifndef FALSE
# define FALSE 0
# define TRUE (!FALSE)
#endif

#define MAX_POINTS 10000

#define TINY 1.0e-20

#define HALF_PI 1.57079

#define SQR(x) ((x)*(x))

extern float x_trans3[],y_trans3[];
extern int nseg2;
extern int representation_ok;

float tempx[MAX_POINTS],tempy[MAX_POINTS];

/* lms conic fit; normalisation: f=1 */
lms_fit(xc,yc,maj,mina,rot)
float *xc,*yc,*maj,*mina,*rot;
{
   int i;
   float m[7];   /* coefficients of best fit ellipse */
   float x_org,y_org;

   representation_ok = TRUE;

   /* shift data to the origin to minimise numerical errors */
   x_org = y_org = 0;
   for (i = 1; i <= nseg2; i++) {
      x_org += x_trans3[i];
      y_org += y_trans3[i];
   }
   x_org /= nseg2; y_org /= nseg2;
   for (i = 1; i <= nseg2; i++) {
      tempx[i] = x_trans3[i] - x_org;
      tempy[i] = y_trans3[i] - y_org;
   }
   /* rescale from 1->511 to 0-1 to prevent enormous numbers */
   for (i = 1; i <= nseg2; i++) {
      tempx[i] /= 512.0;
      tempy[i] /= 512.0;
   }

   determine_ellipse_lms(m);

   if (representation_ok && ellipse_type(m[1],m[2],m[3])) {
      determine_parameters(m[1],m[2],m[3],m[4],m[5],m[6],
                           xc,yc,maj,mina,rot);
      *xc *= 512; *yc *= 512;
      *maj *= 512; *mina *= 512;
      *xc += x_org;
      *yc += y_org;
   }
   else {
      representation_ok = FALSE;
   }
}

determine_parameters(a,b,c,d,e,f,x_cent,y_cent,major_axis,minor_axis,rot_angle)
double a,b,c,d,e,f;
double *x_cent,*y_cent,*major_axis,*minor_axis,*rot_angle;
{
    float ca,sa;
    float t,u,v,w;
    float mina,maja;

    if (a == c)
       *rot_angle = 0;
    else
       *rot_angle = atan(b/(a-c)) / 2.0;
    ca = cos(*rot_angle);
    sa = sin(*rot_angle);
    t = a*SQR(ca)+b*ca*sa+c*SQR(sa);
    u = a*SQR(sa)-b*ca*sa+c*SQR(ca);
    v = d*ca+e*sa;
    w = -d*sa+e*ca;
    *x_cent = (-d*c/2.0 + e*b/4.0) / (a*c - b*b/4.0);
    *y_cent = (-a*e/2.0 + d*b/4.0) / (a*c - b*b/4.0);

    maja = sqrt(((SQR(v)/(4*t))+(SQR(w)/(4*u))-f)/t);
    mina = sqrt(((SQR(v)/(4*t))+(SQR(w)/(4*u))-f)/u);
    if (maja < mina) {
        *minor_axis = maja;
        *major_axis = mina;
        *rot_angle += HALF_PI;
    }
    else {
        *minor_axis = mina;
        *major_axis = maja;
    }
}

/* is conic an ellipse? */
int ellipse_type(c1,c2,c3)
float c1,c2,c3;
{
   float f;

   f = c1 * c3 - c2 * c2 / 4;
   return(f > 0);
}

/* least square fit; f = 1 */
determine_ellipse_lms(m)
float m[7];
{
    int i,j;
    float sx,sx2,sx3,sx4,sy,sy2,sy3,sy4,sxy,sx2y,sx3y,sx2y2,sxy2,sxy3;
    float x2,y2,x3,y3;
    float m1[6][6];
    float m2[6];
    int indx[6];
    float d;

    sx = 0;
    sx2 = 0;
    sx3 = 0;
    sx4 = 0;
    sy = 0;
    sy2 = 0;
    sy3 = 0;
    sy4 = 0;
    sxy = 0;
    sx2y = 0;
    sx3y = 0;
    sx2y2 = 0;
    sxy2 = 0;
    sxy3 = 0;
 
    for (j = 1; j <= nseg2; j++) {
       x2 = tempx[j] * tempx[j];
       x3 = x2 * tempx[j];
       y2 = tempy[j] * tempy[j];
       y3 = y2 * tempy[j];
 
       sx += tempx[j];
       sx2 += x2;
       sx3 += x3;
       sx4 += (x3 * tempx[j]);
       sy += tempy[j];
       sy2 += y2;
       sy3 += y3;
       sy4 += (y3 * tempy[j]);
       sxy += (tempx[j] * tempy[j]);
       sx2y = sx2y + x2 * tempy[j];
       sx3y = sx3y + x3 * tempy[j];
       sx2y2 = sx2y2 + x2 * y2;
       sxy2 = sxy2 + tempx[j] * y2;
       sxy3 = sxy3 + tempx[j] * y3;
    }
 
    m1[1][1] = sx4;
    m1[2][1] = sx3y;
    m1[3][1] = sx2y2;
    m1[4][1] = sx3;
    m1[5][1] = sx2y;
    m1[1][2] = sx3y;
    m1[2][2] = sx2y2;
    m1[3][2] = sxy3;
    m1[4][2] = sx2y;
    m1[5][2] = sxy2;
    m1[1][3] = sx2y2;
    m1[2][3] = sxy3;
    m1[3][3] = sy4;
    m1[4][3] = sxy2;
    m1[5][3] = sy3;
    m1[1][4] = sx3;
    m1[2][4] = sx2y;
    m1[3][4] = sxy2;
    m1[4][4] = sx2;
    m1[5][4] = sxy;
    m1[1][5] = sx2y;
    m1[2][5] = sxy2;
    m1[3][5] = sy3;
    m1[4][5] = sxy;
    m1[5][5] = sy2;
 
    m2[1] = -sx2;
    m2[2] = -sxy;
    m2[3] = -sy2;
    m2[4] = -sx;
    m2[5] = -sy;
 
    /* solve simultaneous equations */
    ludcmp(m1,5,indx,&d);
	if (representation_ok)
    	lubksb(m1,5,indx,m2);
	else
		return;
 
    m[1] = m2[1];
    m[2] = m2[2];
    m[3] = m2[3];
    m[4] = m2[4];
    m[5] = m2[5];
    m[6] = 1;
}

/* LU Decomposition - from Numerical Recipes in C */
/* not-dynamic - assumes n=5; i.e. a[1..5][1..5], indx[1..5] */
ludcmp(a,n,indx,d)
float a[6][6];
int n,indx[6];
float *d;
{
        int i,imax,j,k;
        float big,dum,sum,temp;
        float vv[6];

        *d = 1.0;
        for (i = 1; i <= n; i++) {
                big = 0.0;
                for (j = i; j <= n; j++)
                        if ((temp = fabs(a[i][j])) > big) big = temp;
                if (big == 0.0) {
					representation_ok = FALSE;
					return;
				}
                vv[i] = 1.0/big;
        }
        for (j = 1; j <= n; j++) {
                for (i = 1; i < j; i++) {
                        sum = a[i][j];
                        for (k = 1; k < i; k++) sum -= a[i][k] * a[k][j];
                        a[i][j] = sum;
                }
                big = 0.0;
                for (i = j; i <= n;i++) {
                        sum = a[i][j];
                        for (k = 1; k < j; k++) sum -= a[i][k] * a[k][j];
                        a[i][j] = sum;
                        if ( (dum = vv[i]*fabs(sum)) >= big) {
                                big = dum;
                                imax = i;
                        }
                }
                if (j != imax) {
                        for (k = 1; k <= n; k++) {
                                dum = a[imax][k];
                                a[imax][k] = a[j][k];
                                a[j][k] = dum;
                        }
                        *d = -(*d);
                        vv[imax] = vv[j];
                }
                indx[j] = imax;
                if (a[j][j] == 0.0) a[j][j] = TINY;
                if (j != n) {
                        dum = 1.0/(a[j][j]);
                        for (i = j+1; i <= n ;i++) a[i][j] *= dum;
                }
        }
}

/* forward substitution & backsubstitution */
/* use with LU Decomposition - from Numerical Recipes in C */
/* not-dynamic - assumes n=5; i.e. a[1..5][1..5], indx[1..5] */
lubksb(a,n,indx,b)
float a[6][6],b[6];
int n,indx[6];
{
        int i,ii=0,ip,j;
        float sum;

        for (i=1;i<=n;i++) {
                ip = indx[i];
                sum = b[ip];
                b[ip] = b[i];
                if (ii)
                        for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
                else if (sum) ii = i;
                b[i] = sum;
        }
        for (i=n;i>=1;i--) {
                sum = b[i];
                for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
                b[i] = sum/a[i][i];
        }
}
