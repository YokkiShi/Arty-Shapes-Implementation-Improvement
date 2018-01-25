/*
 * Implentation of lms circle fitting algorithm from:
 * 'A simple approach for the estimation of circular arc centre and
 * its *radius'
 * by S M Thomas and Y T Chan, CVGIP Vol 45, pp 362-370, 1989
 */

#define FALSE 0
#define TRUE (!FALSE)

#define SQR(x) ((x) * (x))

extern int representation_ok;

compute_circle(x_cent,y_cent,radius2,nseg2,X,Y)
float *x_cent,*y_cent,*radius2;
int nseg2;
float X[],Y[];
{
   float sx,sy,sx2,sy2,sxy,sx3,sy3,sx2y,sxy2;
   float a1,a2,b1,b2,c1,c2;
   int loop1;
   float temp;

   representation_ok = TRUE;
   loop1 = 0;
   /* initialise variables */;
   sx = 0.0;
   sy = 0.0;
   sx2 = 0.0;
   sy2 = 0.0;
   sxy = 0.0;
   sx3 = 0.0;
   sy3 = 0.0;
   sx2y = 0.0;
   sxy2 = 0.0;
   /* compute summations */
   for (loop1=1; loop1<=nseg2; loop1++) {
      sx = sx + X[loop1];
      sy = sy + Y[loop1];
      sx2 = sx2 + SQR(X[loop1]);
      sy2 = sy2 + SQR(Y[loop1]);
      sxy = sxy + X[loop1] * Y[loop1];
      sx3 = sx3 + X[loop1] * X[loop1] * X[loop1];
      sy3 = sy3 + Y[loop1] * Y[loop1] * Y[loop1];
      sx2y = sx2y + SQR(X[loop1]) * Y[loop1];
      sxy2 = sxy2 + SQR(Y[loop1]) * X[loop1];
   }
   /* compute a's,b's,c's */
   a1 = 2.0 * (SQR(sx) - sx2 * nseg2);
   a2 = 2.0 * (sx * sy - sxy * nseg2);
   b1 = a2;
   b2 = 2.0 * (SQR(sy) - sy2 * nseg2);
   c1 = sx2 * sx - sx3 * nseg2 + sx * sy2 - sxy2 * nseg2;
   c2 = sx2 * sy - sy3 * nseg2 + sy * sy2 - sx2y * nseg2;
   /* compute centre */
   temp = a1*b2 - a2*b1;
   if (temp == 0) {
        representation_ok = FALSE;
        return;
   }
   *x_cent = (c1 * b2 - c2 * b1) / (a1 * b2 - a2 * b1);
   *y_cent = (a1 * c2 - a2 * c1) / ( a1 * b2 - a2 * b1);
   /* compute radius squared */
   *radius2 = (sx2 - 2.0 * sx * *x_cent + SQR(*x_cent) * nseg2
             + sy2 - 2.0 * sy * *y_cent + SQR(*y_cent) * nseg2) / nseg2;
}
