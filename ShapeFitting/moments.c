/*
 * compute principle axis of data using central moments
 *
 * Paul Rosin
 * May 1993
 *
 * update to convert atan to atan2 - avoids divide-by-zero errors
 * Paul Rosin
 * January 1995
 */

#include <math.h>

float compute_orientation(no_points,X,Y)
int no_points;
float X[],Y[];
{
   int loop1;
   float x_avg,y_avg;
   float compute_moment();
   float m11,m20,m02;
   float orient;

   /* compute centroid of data */
   x_avg = y_avg = 0;
   for (loop1 = 1; loop1 <= no_points; loop1++) {
      x_avg += X[loop1];
      y_avg += Y[loop1];
   }
   x_avg /= (float) no_points;
   y_avg /= (float) no_points;

   m11 = compute_moment(X,Y,no_points,x_avg,y_avg,1,1);
   m20 = compute_moment(X,Y,no_points,x_avg,y_avg,2,0);
   m02 = compute_moment(X,Y,no_points,x_avg,y_avg,0,2);

   /***
   orient = 0.5 * atan((2.0 * m11) / (m20 - m02));
   ***/
   orient = 0.5 * atan2(2.0*m11,m20-m02);

   return(orient);
}

/* compute p,q'th central moment */
float compute_moment(X,Y,no_points,x_avg,y_avg,p,q)
float X[],Y[];
float x_avg,y_avg;
int no_points;
int p,q;
{
    int i;
    float m = 0;
    float POW();

    for (i = 1; i <= no_points; i++)
        m += pow2((X[i] - x_avg),(float)p) * pow2((Y[i] - y_avg),(float)q);

    return(m);
}
