/*
 * find convex hull of a set of points
 *
 * input & output format: pixel lists (integer, float, or curvature)
 *
 * program converted from code provided by (I think) Michael Swain
 *
 * doesn't work well on some obscure examples - fixed with hack
 *
 * Paul Rosin
 * April 1994
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef FALSE
#define FALSE 0
#define TRUE (!FALSE)
#endif

#define MAX_POINTS 100000

/* format used for data in pixel lists */
#define INTEGER 0
#define FLOAT 1
#define CURVATURE 2

/* temporary arrays for reading in pixel data */
int no_pixels;
float tmpx[MAX_POINTS];
float tmpy[MAX_POINTS];

/* file i/o stuff */
char file_type[50];
int file_data_type;
FILE *fp1, *fp2;

main(argc, argv)
int argc;
char *argv[];
{
    double p[MAX_POINTS * 2];
    int n;
    int ch[MAX_POINTS];
    int nch;
    static char directions[] = "norm";
    double epsilon = 0.001;
    int i;

    if (argc != 3) {
        fprintf(stderr, "%s infile outfile\n", argv[0]);
        exit(-1);
    }

    if ((fp1 = fopen(argv[1], "r")) == NULL) {
        fprintf(stderr, "cant open %s\n", argv[1]);
        exit(-1);
    }
    if ((fp2 = fopen(argv[2], "w")) == NULL) {
        fprintf(stderr, "cant open %s\n", argv[2]);
        exit(-1);
    }

    /* read pixel data */
    read_data(fp1);

    for (i = 0; i < no_pixels; i++) {
        p[i * 2] = tmpx[i];
        p[i * 2 + 1] = tmpy[i];
    }
    n = no_pixels;

    convex_hull(p, n, ch, &nch, directions, epsilon);

    fprintf(fp2, "pixel\nlist: 0\n");
    for (i = 0; i < nch; i++) {
        n = ch[i];
        fprintf(fp2, "%.0f %.0f\n", tmpx[n], tmpy[n]);
    }
    /* first point again to form closed boundary */
    n = ch[0];
    fprintf(fp2, "%.0f %.0f\n", tmpx[n], tmpy[n]);

    fprintf(fp2, "-1 -1\n");
    fclose(fp2);
}

read_data(fp)
FILE *fp;
{
    int endoffile, list_no;

    /* read magic word for format of file */
    fscanf(fp, "%s\n", file_type);
    if (strcmp(file_type, "pixel") == 0) {
        fprintf(stderr, "integer pixel data file\n");
        file_data_type = INTEGER;
    } else if (strcmp(file_type, "pixel_float") == 0) {
        fprintf(stderr, "float pixel data file\n");
        file_data_type = FLOAT;
    } else if (strcmp(file_type, "pixel_curvature") == 0) {
        fprintf(stderr, "curvature pixel data file\n");
        file_data_type = CURVATURE;
    } else {
        fprintf(stderr,
            "input file not pixel, pixel_float or pixel_curvature type - aborting\n");
        exit(-1);
    }

    switch (file_data_type) {
    case INTEGER:
        read_integer_data(fp, &endoffile, &list_no);
        break;
    case FLOAT:
        read_float_data(fp, &endoffile, &list_no);
        break;
    case CURVATURE:
        read_curvature_data(fp, &endoffile, &list_no);
        break;
    }

    if (no_pixels >= MAX_POINTS) {
        fprintf(stderr, "ERROR: too many pixels\n");
        exit(-1);
    }

    fclose(fp);
}

read_integer_data(fp, endoffile, list_no)
FILE *fp;
int *endoffile, *list_no;
{
    char dumstring[50];
    int j;
    int tx, ty;

    fscanf(fp, "%s %d\n", dumstring, list_no);
    j = -1;
    do {
        j++;
        fscanf(fp, "%d %d\n", &tx, &ty);
        tmpx[j] = tx;
        tmpy[j] = ty;
    } while (tmpx[j] != -1);
    *endoffile = (tmpy[j] == -1);
    if (feof(fp) && !(*endoffile)) {
        fprintf(stderr, "Incorrectly terminated file - assuming end of list\n");
        *endoffile = TRUE;
    }
    no_pixels = j;
}

read_float_data(fp, endoffile, list_no)
FILE *fp;
int *endoffile, *list_no;
{
    char dumstring[50];
    int j;
    float tx, ty;

    fscanf(fp, "%s %d\n", dumstring, list_no);
    j = -1;
    do {
        j++;
        fscanf(fp, "%f %f\n", &tx, &ty);
        tmpx[j] = tx;
        tmpy[j] = ty;
    } while (tmpx[j] != -1);
    *endoffile = (tmpy[j] == -1);
    if (feof(fp) && !(*endoffile)) {
        fprintf(stderr, "Incorrectly terminated file - assuming end of list\n");
        *endoffile = TRUE;
    }
    no_pixels = j;
}

read_curvature_data(fp, endoffile, list_no)
FILE *fp;
int *endoffile, *list_no;
{
    char dumstring[50];
    int j;
    float sigma;
    float tx, ty, K;
    int label;

    fscanf(fp, "%s %d\n", dumstring, list_no);
    fscanf(fp, "%s %f\n", dumstring, &sigma);
    j = -1;
    do {
        j++;
        fscanf(fp, "%f %f %f %d\n", &tx, &ty, &K, &label);
        tmpx[j] = tx;
        tmpy[j] = ty;
    } while (tmpx[j] != -1);
    *endoffile = (tmpy[j] == -1);
    if (feof(fp) && !(*endoffile)) {
        fprintf(stderr, "Incorrectly terminated file - assuming end of list\n");
        *endoffile = TRUE;
    }
    no_pixels = j;
}

/* ############################################################ */

/* 
convexhull.c

Find the convex hull of a set of points.

The algorithm used is as describe in 

Shamos, Michael,  "Problems
in Computational Geometry"  May, 1975 (a set of photocopies -- 
QA 447.S52 1983 in Carlson).  

It originally appeared in a more complicated form in 

Jarvis, R.A., "On the Identification of the Convex Hull of a 
Finite Set of Points in the Plane", Info. Proc. Letters 2(1973),
pp. 18-21.

The algorithm is of complexity k*n, where n is the number of points in the
input and k the number of points in the convex hull.

usage:
convex_hull(p,n,ch,&nch,directions);
where p is an n*2 array of doubles containing the set of points,
n is the number of points,
ch is an array of size n integers to hold the list of points in the 
    convex hull, numbered 0 to nch-1;
In nch the number of points in the convex hull is returned.
directions is either "full" or "norm".  If directions="full" all the
possible points on the convex hull are returned.  If directions="norm"
a minimal set of points to describe the convex hull is returned.

epsilon is the angle tolerance in radians.  When two angles are closer than
PI radians they are considered equal.
*/

/* ### CONTENTS OF INCLUDE FILE INSERTED HERE ### */

/* incl.h */

#define boolean int
#define true 1
#define false 0

#define ERROR -1
#define OK 0

#define EPSILON 1.e-6
#define BIGNUM 1.e20
#define PI 3.1415927

#define sqr(x) ((x)*(x))

/* ### -------------------------------------- ### */

static int debug = 0;

convex_hull(p, n, ch, p_nch, directions, epsilon)
double *p;
int n, *p_nch, *ch;
char *directions;
double epsilon;
{
    double phi, max_phi, dist, max_cen_dist, min_dist, max_dist,
        x, y, cen_x, cen_y, xp, yp, xx, yy, m[2][2];
    int i, i_keep, vertex, furthest, ch_vert;
    boolean full;

    if (!strcmp(directions, "full")) {
        full = true;
    } else if (!strcmp(directions, "norm")) {
        full = false;
    } else {
        fprintf(stderr, "convex_hull: invalid argument \"%s\"\n", directions);
        exit(ERROR);
    }
    centroid(p, n, &cen_x, &cen_y);
    /* find the point furthest from the centroid */
    max_cen_dist = 0.;
    for (i = 0; i < n; i++) {
        x = p[2 * i];
        y = p[2 * i + 1];
        dist = sqrt(sqr(x - cen_x) + sqr(y - cen_y));
        if (dist > max_cen_dist) {
            max_cen_dist = dist;
            furthest = i;
        }
    }

    /* 
       Determine rotation matrix so that coordinate system for determining
       the orientations of line segments is wrt the line from the point
       under consideration to the centroid.  Then all angles will be 
       < 90 degrees.  To maintain the strict inequality the furthest 
       point along the extreme ray will be used as the next point on 
       the convex hull.
     */

    make_rotation_matrix((cen_x - p[furthest * 2]) / max_cen_dist,
                 (cen_y - p[furthest * 2 + 1]) / max_cen_dist, (double *) m);
    ch_vert = 0;
    vertex = furthest;

    /* HACK: to avoid loop with problem data - PLR */
    ch[1] = -999;

    do {
        ch[ch_vert] = vertex;
        /* Find the ray with the maximum and minimum angle in the new 
           coordinate system */
        if (debug)
            printf("vertex %d\n", vertex);
        max_phi = -BIGNUM;
        min_dist = BIGNUM;
        max_dist = 0.;
        for (i = 0; i < n; i++) {
            /* calculate phi, the angle in the new coordinate system */
            x = p[2 * i] - p[2 * vertex];
            y = p[2 * i + 1] - p[2 * vertex + 1];
            xp = m[0][0] * x + m[0][1] * y;
            yp = m[1][0] * x + m[1][1] * y;
            if (debug)
                printf("\ti %d x %f y %f xp %f yp %f\n", i, x, y, xp, yp);
            if ((dist = sqrt(sqr(x) + sqr(y))) > EPSILON) {
                phi = atan2(yp, xp);
                if (debug)
                    printf("\tphi %f\n", phi);
                if (phi > max_phi - epsilon) {
                    if (full) {
                        /* use the closest point */
                        if (phi > max_phi + epsilon || dist < min_dist) {
                            min_dist = dist;
                            max_phi = phi;
                            i_keep = i;
                        }
                    } else {
                        /* use the furthest point */
                        if (phi > max_phi + epsilon || dist > max_dist) {
                            max_dist = dist;
                            max_phi = phi;
                            i_keep = i;
                        }
                    }
                    if (debug)
                        printf("\t\tmax_phi %f i_keep %d\n", max_phi,
                               i_keep);
                }
            }
        }
        vertex = i_keep;
        xx = cen_x - p[vertex * 2];
        yy = cen_y - p[vertex * 2 + 1];
        dist = sqrt(sqr(xx) + sqr(yy));
        make_rotation_matrix(xx / dist, yy / dist, (double *) m);
        ch_vert++;
        /*
           } while (vertex != ch[0]);
         */
        /* HACK: to avoid loop with problem data - PLR */
    } while ((vertex != ch[0]) && (vertex != ch[1]));
    *p_nch = ch_vert;
}

/* 
Determines the centroid of the set of n points in the n*2 array of
doubles pts and returns its x-coordinate and y-coordinate in p_cen_x
and p_cen_y respectively.
*/
centroid(pts, n, p_cen_x, p_cen_y)
double *pts, *p_cen_x, *p_cen_y;
int n;
{
    double sumx, sumy;
    int i;

    sumx = sumy = 0;
    for (i = 0; i < n; i++) {
        sumx += pts[2 * i];
        sumy += pts[2 * i + 1];
    }
    *p_cen_x = sumx / n;
    *p_cen_y = sumy / n;
}


/* 
Given a unit vector (vx,vy) system, finds the matrix m that
corresponds to rotating the coordinate system so that its x-axis 
lies in the direction defined by (vx,vy).
*/
make_rotation_matrix(cos_theta, sin_theta, m)
double cos_theta, sin_theta, *m;
{
    *m = cos_theta;
    *(m + 1) = sin_theta;
    *(m + 2) = -sin_theta;
    *(m + 3) = cos_theta;
}
