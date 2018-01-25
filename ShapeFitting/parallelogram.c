/*
 *    calculate minimal enclosing parallelogram of CONVEX polygon
 *    inefficient approach!
 *
 *    garbage output if input not convex!
 *
 *    could use for shape measure
 *
 *    Paul Rosin
 *    October 2005
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.1415927

#define MAX_PIXELS 10000

#ifndef FALSE
#define FALSE 0
#define TRUE (!FALSE)
#endif

#define ABS(x)       (((x)<0.0)? (-(x)): (x))
#define SQR(x)       ((x)*(x))

int no_pixels;
int x[MAX_PIXELS],y[MAX_PIXELS];

double orientation[MAX_PIXELS];
double max_dist[MAX_PIXELS];
int max_loc[MAX_PIXELS];

double ise();

main(argc,argv)
int argc;
char *argv[];
{
    FILE *fp1;
    char file_type[50];
    int i,j,k,ii,kk;
    int endoffile;
    double dist;
    double dx,dy;
    double ang,area,min_area;

    if (argc != 2) {
        printf("usage: %s input_file \n",argv[0]);
        exit(-1);
    }

    if ((fp1=fopen(argv[1],"r")) == NULL) {
        printf("cant open %s\n",argv[1]);
        exit(-1);
    }
    /* read magic word for format of file */
    fscanf(fp1,"%s\n",file_type);
    j = strcmp(file_type,"pixel");
    if (j != 0) {
        printf("not link data file - aborting\n");
        exit(-1);
    }

    printf("pixel\n");
    do {
        read_link_data(fp1,&endoffile);

        /* duplicate end points */
        if ((x[0] != x[no_pixels-1]) || (y[0] != y[no_pixels-1])) {
            x[no_pixels] = x[0];
            y[no_pixels] = y[0];
        }
        else
            no_pixels--;

        /* for each edge i->i+1 */
        for (i = 0; i < no_pixels; i++) {
            dx = x[i] - x[i+1];
            dy = y[i] - y[i+1];
            orientation[i] = atan2(dy,dx);

            /* calculate distance to each point j and save max */
            max_dist[i] = 0;
            for (j = 0; j < no_pixels; j++) {
                dist = ise(x[i],y[i],x[i+1],y[i+1],x[j],y[j]);
                if (dist > max_dist[i]) {
                    max_dist[i] = dist;
                    max_loc[i] = j;
                }
            }
        }

        min_area = 9e9;
        for (i = 0; i < no_pixels-1; i++) {
            for (k = i+1; k < no_pixels; k++) {
                //printf("%2d/  (%3d, %3d)  -> %d = %.1f\n",i,x[i],y[i],max_loc[i],max_dist[i]);

                ang = orientation[i] - orientation[k];
                area = (max_dist[i] * max_dist[k]) / ABS(sin(ang));
                //printf("%d %d  AREA %f\n",i,k,area);

                if (area < min_area) {
                    min_area = area;
                    ii = i;
                    kk = k;
                }
            }
        }

        //plot(ii,kk);
        vertices(ii,kk);

        printf("AREA %f\n",min_area);

    } while (!endoffile);
    fclose(fp1);
}

/* plot bounding lines of parallelogram */
plot(i,k)
int i,k;
{
    int j;
    int xx,yy;

    printf("list: 0\n");
    xx = -500;
    yy = tan(orientation[i])*(xx-x[i])+y[i];
    printf("%d %d\n",xx,yy);
    xx = 700;
    yy = tan(orientation[i])*(xx-x[i])+y[i];
    printf("%d %d\n",xx,yy);
    printf("-1 0\n");

    printf("list: 0\n");
    j = max_loc[i];
    xx = -500;
    yy = tan(orientation[i])*(xx-x[j])+y[j];
    printf("%d %d\n",xx,yy);
    xx = 700;
    yy = tan(orientation[i])*(xx-x[j])+y[j];
    printf("%d %d\n",xx,yy);
    printf("-1 0\n");

    /* -----------------------------------------*/

    printf("list: 0\n");

    xx = -500;
    yy = tan(orientation[k])*(xx-x[k])+y[k];
    printf("%d %d\n",xx,yy);
    xx = 700;
    yy = tan(orientation[k])*(xx-x[k])+y[k];
    printf("%d %d\n",xx,yy);
    printf("-1 0\n");

    printf("list: 0\n");
    j = max_loc[k];
    xx = -500;
    yy = tan(orientation[k])*(xx-x[j])+y[j];
    printf("%d %d\n",xx,yy);
    xx = 700;
    yy = tan(orientation[k])*(xx-x[j])+y[j];
    printf("%d %d\n",xx,yy);
    printf("-1 -1\n");
}

/* plot proper lines between vertices of parallelogram */
vertices(i,k)
int i,k;
{
    int xi,yi;
    double m,M,a,A,b,B;
    int ii,kk;
    int xx,yy;

    printf("list: 0\n");

    ii = max_loc[i];
    kk = max_loc[k];

    m = tan(orientation[i]);
    M = tan(orientation[k]);

    a = x[i]; b = y[i];
    A = x[k]; B = y[k];
    xi = (-b + B + a*m - A*M)/(m - M);
    yi = (B*m - (b + (-a + A)*m)*M)/(m - M);
    printf("%d %d\n",xi,yi);
    xx = xi; yy = yi;

    kk = max_loc[k];
    a = x[i]; b = y[i];
    A = x[kk]; B = y[kk];
    xi = (-b + B + a*m - A*M)/(m - M);
    yi = (B*m - (b + (-a + A)*m)*M)/(m - M);
    printf("%d %d\n",xi,yi);

    a = x[ii]; b = y[ii];
    A = x[kk]; B = y[kk];
    xi = (-b + B + a*m - A*M)/(m - M);
    yi = (B*m - (b + (-a + A)*m)*M)/(m - M);
    printf("%d %d\n",xi,yi);

    a = x[ii]; b = y[ii];
    A = x[k]; B = y[k];
    xi = (-b + B + a*m - A*M)/(m - M);
    yi = (B*m - (b + (-a + A)*m)*M)/(m - M);
    printf("%d %d\n",xi,yi);
    printf("%d %d\n",xx,yy);
    printf("-1 -1\n");
}

read_link_data(fp,endoffile)
FILE *fp;
int *endoffile;
{
    char dumstring[50];
    int j;

    fscanf(fp,"%s %d\n",dumstring,&j);
    j = -1;
    do{
       j++;
       fscanf(fp,"%d %d\n",&x[j],&y[j]);
    } while(x[j] != -1);
    *endoffile = (y[j] == -1);
    if (feof(fp) && !(*endoffile)) {
        fprintf(stderr,"Incorrectly terminated file - assuming end of list\n");
        *endoffile = TRUE;
    }
    no_pixels = j;
    if (no_pixels >= MAX_PIXELS) {
        fprintf(stderr,"ERROR: Too many pixels\n");
        exit(-1);
    }
}

double ise(xst,yst,xfi,yfi,xt,yt)
int xst,yst,xfi,yfi,xt,yt;
{
    double dist;

    dist = (yst - yfi) * xt - (xst - xfi) * yt -
            xfi * yst + xst * yfi;
    dist = SQR(dist) /
           (double)(SQR(xst - xfi) + SQR(yst - yfi));
    dist = sqrt(dist);

    return dist;
}
