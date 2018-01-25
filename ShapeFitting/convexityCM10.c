/* GENETIC ALGORITHM 2, ONE-POINT CROSSOVER AND MUTATION
 * REPAIRS CHROMOSOME BY REPLACING DERIVED SHAPE WITH ITS CONVEX HULL
 * WILL NEED CHANGING IF REAL NUMBERS USED FOR COORDINATES
 *
 * does not write repaired chromosome back to the population
 *
 * measure convexity of a region
 * also measures intrusiveness & protrusiveness
 *
 * outputs data for perturb program
 *
 * ############### NEED TO CHANGE CP NORMALISATION ###############
 *
 * ratio of area of region to area of convex hull
 * 1 = perfect convexity; 0 = terrible convexity
 *
 * input is a pixel list of boundary coordinates
 *
 * Paul Rosin
 * December 1998
 *
 * ammended by Christine Mumford, August, September, October, November 2003
 * to generate all subsets of pixels of size >= 3 pixels
 * to hunt for a convex shape with minimum XOR-ed area with original polygon
 *
 * slight mods by Paul Rosin, November 2003
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gpc.h"

#ifndef FALSE
#define FALSE 0
#define TRUE (!FALSE)
#endif

#define MAX_POINTS 30000
#define VERY_BIG 1e9
#define POPSIZE 1000
#define NUMBER_GENERATIONS 100
#define NO_MUTATIONS 2

#define ABS(x)    (((x)<0.0)? (-(x)): (x))

FILE *results_file = NULL;
FILE *results_file2 = NULL;

double poly_area();

gpc_polygon poly1;
gpc_polygon poly4, result2, result3, result4;
gpc_polygon result_intrusion,result_protrusion;

/* hack added to trap gpc errors - PLR */
int gpc_error_flag;

/* input pixel data */
int no_pixels;
int x[MAX_POINTS];
int y[MAX_POINTS];

double convexity;
double area_original,area_ch,area_ch2,area_xor,area_int;
double area_intrusion,area_protrusion;
double best_error;

int bits_set;
int best_points[MAX_POINTS];

/* random number seeds */
long rs1, rs2, rs3;

double ran3();

/* population of bits to identify vertices included in the derived shape */
int population[POPSIZE][MAX_POINTS];

/* array storing XORed error values between original and derived shapes */
float xor_error[POPSIZE];

/* array storing number of bits set to 1 */
int bit_set_count[POPSIZE];

main(argc,argv)
int argc;
char *argv[];
{

    char filename[20];
    gpc_polygon poly2, result;
    int first;
    int flag = 0;

    /* convex hull data for original shape*/
    int no_pixels2;
    int x2[MAX_POINTS];
    int y2[MAX_POINTS];

    int i, j;

    if (argc == 5) {
        flag = atoi(argv[4]);
    }
    else if (argc != 4) {
       fprintf(stderr,"%s infile outfile1 outfile2 [0/1 = initial randomisation]\n",argv[0]);
       exit(-1);
    }

    /* read pixel data */
    read_link_data(argv[1],x,y,&no_pixels);

    if (no_pixels < 4) {
       fprintf(stderr,"ERROR: %d is too few pixels\n",no_pixels);
       exit(-1);
    }

    /* prepare output file */

    if ((results_file=fopen(argv[2],"w")) == NULL) {
        fprintf(stderr,"cant open %s\n",argv[2]);
        exit(-1);
    }

    if ((results_file2=fopen(argv[3],"w")) == NULL) {
        fprintf(stderr,"cant open %s\n",argv[3]);
        exit(-1);
    }

    if (flag == 1) {
        rs1 = ran3() * 100000;
        rs2 = ran3() * 100000;
        rs3 = ran3() * 100000;
        printf("randomising seeds: %d %d %d\n",rs1,rs2,rs3);
    }
    else {
        rs1 = 1000; rs2 = 2000; rs3 = 3000;
    }

    /*
    printf("type in 3 random number seeds separated by return\n");
    scanf("%d %d %d", &rs1, &rs2, &rs3);
    */

    /* create a polygon structure with the original data */
    gpc_create_polygon(x,y,no_pixels,&poly1);

    /* compute convex hull */
    convex_hull(x,y,no_pixels,x2,y2,&no_pixels2);

    /* create a polygon structure with the convex hull */
    gpc_create_polygon(x2,y2,no_pixels2,&poly2);

    /* XOR polygons */
    gpc_polygon_clip(GPC_XOR, &poly1, &poly2, &result);

    /* calculate areas */
    area_original =  poly_area(&poly1);
    area_ch =        poly_area(&poly2);
    area_xor =       poly_area(&result);

    printf("areas: original %g, convex hull %g, xor %g\n",
        area_original,area_ch,area_xor);

    convexity = area_original / area_ch;

    printf("convexity %f\n",convexity);

    best_error = VERY_BIG;

    generate_population();
    printf("population generated\n");

    ga();

    /* print out vertices of best shape, plus error */
    printf("best error between original and new shape is %g\n",best_error);
    printf("original weighted error %f\n",1-best_error/area_original);
    printf("scaled ORIGINAL WEIGHTED error %f\n",
        255*(1-best_error/area_original));

	// C_P
    printf("renormalised original weighted error %f\n",1.0/(1.0+best_error/area_original));

    //printf("intersection weighted error %f\n",best_error/area_int);
    //printf("CH weighted error %f\n",best_error/area_ch2);

	// C_Q
    printf("CH normalised weighted error %f\n",1.0/(1.0+best_error/area_ch2));

    printf("intrusiveness %f\n",area_intrusion/area_ch2);
    printf("protrusiveness %f\n",area_protrusion/(area_protrusion+area_ch2));

    print_subset(best_points,no_pixels, bits_set);

    /*
    fprintf(results_file,"%d\n", no_pixels);

    for (j=0;j<no_pixels;j++) {
        fprintf(results_file,"%d %d\n", x[j], y[j]);
    }

    fprintf(results_file,"%d\n", bits_set);
    */

    fprintf(results_file,"pixel\nlist: 0\n");
    first = -1;
    for (j=0;j<no_pixels;j++)  {
        if (best_points[j] == 1) {
            fprintf(results_file,"%d %d\n", x[j], y[j]);
            if (first == -1) first = j;
        }
    }
    fprintf(results_file,"%d %d\n-1 -1\n",x[first],y[first]);

    /* ---------- output for perturb program ---------- */
    fprintf(results_file2,"%d\n", no_pixels);

    for (j=0;j<no_pixels;j++)
        fprintf(results_file2,"%d %d\n", x[j], y[j]);

    fprintf(results_file2,"%d\n", bits_set);
            
    for (j=0;j<no_pixels;j++) 
        if (best_points[j] == 1)
            fprintf(results_file2,"%d %d\n", x[j], y[j]);

    fclose(results_file);
    fclose(results_file2);
}

/* ---------------------------------------------------------------------- */

convex_hull(int x[], int y[], int no_pixels, int x2[], int y2[], int *no_pixels2)
{
    double p[MAX_POINTS*2];
    int ch[MAX_POINTS];
    int nch;
    static char directions[] = "norm";
    double epsilon = 0.001;
    int i,n;

    for (i = 0; i < no_pixels; i++) {
        p[i*2] = x[i];
        p[i*2+1] = y[i];
    }

    convex_hull0(p,no_pixels,ch,&nch,directions,epsilon);

    for (i = 0; i < nch; i++) {
        n = ch[i];
        x2[i] = x[n];
        y2[i] = y[n];
    }
    *no_pixels2 = nch;
}

read_link_data(char *fname, int x[], int y[], int *no_pixels)
{
    char dumstring[50];
    int j;
    int endoffile;
    FILE *fp;
    char file_type[50];

    if ((fp=fopen(fname,"r")) == NULL) {
        fprintf(stderr,"cant open %s\n",fname);
        exit(-1);
    }

    /* read magic word for format of file */
    fscanf(fp,"%s\n",file_type);
    if (strcmp(file_type,"pixel") != 0) {
        fprintf(stderr,"input file not pixel type - aborting\n");
        exit(-1);
    }

    fscanf(fp,"%s %d\n",dumstring,&j);
    j = -1;
    do{
       j++;
       fscanf(fp,"%d %d\n",&x[j],&y[j]);
    } while(x[j] != -1);
    endoffile = (y[j] == -1);
    if (feof(fp) && !(endoffile)) {
        fprintf(stderr,"Incorrectly terminated file - assuming end of list\n");
        endoffile = TRUE;
    }
    *no_pixels = j;
    fclose(fp);
}

/* calculate ABSOLUTE summed area of polygons */
double poly_area(gpc_polygon *p)
{
  int c, v;
  double area,sum = 0;
  int flag;
  int x[MAX_POINTS],y[MAX_POINTS];
  double xc,yc;

  for (c = 0; c < p->num_contours; c++) {
      for (v = 0; v < p->contour[c].num_vertices; v++) {
          x[v] = p->contour[c].vertex[v].x;
          y[v] = p->contour[c].vertex[v].y;
      }
      flag = polyCentroid(x,y,p->contour[c].num_vertices,&xc,&yc,&area);
      if (flag == 0)
          sum += ABS(area);
  }

  return(sum);
}

/* ############################################################ */

/*********************************************************************
polyCentroid: Calculates the centroid (xCentroid, yCentroid) and area
of a polygon, given its vertices (x[0], y[0]) ... (x[n-1], y[n-1]). It
is assumed that the contour is closed, i.e., that the vertex following
(x[n-1], y[n-1]) is (x[0], y[0]).  The algebraic sign of the area is
positive for counterclockwise ordering of vertices in x-y plane;
otherwise negative.

Returned values:  0 for normal execution;  1 if the polygon is
degenerate (number of vertices < 3);  and 2 if area = 0 (and the
centroid is undefined).
**********************************************************************/
int polyCentroid(int x[], int y[], int n,
         double *xCentroid, double *yCentroid, double *area)
{
     register int i, j;
     double ai, atmp = 0, xtmp = 0, ytmp = 0;

     if (n < 3) return 1;
     for (i = n-1, j = 0; j < n; i = j, j++) {
         ai = x[i] * y[j] - x[j] * y[i];
         atmp += ai;
         xtmp += (x[j] + x[i]) * ai;
         ytmp += (y[j] + y[i]) * ai;
     }
     *area = atmp / 2;

     /* !!!!! make area unsigned !!!!! */
     *area = ABS(*area);

     if (atmp != 0) {
         *xCentroid = xtmp / (3 * atmp);
         *yCentroid = ytmp / (3 * atmp);
         return 0;
     }
     return 2;
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
convex_hull0(p,n,ch,&nch,directions);
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

static int debug=0;

convex_hull0(p,n,ch,p_nch,directions,epsilon)
double *p;
int n, *p_nch, *ch;
char *directions;
double epsilon;
{
    double phi, max_phi, dist, max_cen_dist, min_dist, max_dist,
        x, y, cen_x, cen_y, xp, yp,
    xx, yy, m[2][2];
    int i, i_keep, vertex, furthest, ch_vert;
    boolean full;

    if (!strcmp(directions,"full")) {
    full = true;
    }
    else if (!strcmp(directions,"norm")) {
    full = false;
    }
    else {
    fprintf(stderr,"convex_hull: invalid argument \"%s\"\n",directions);
    exit(ERROR);
    }
    centroid(p,n,&cen_x,&cen_y);
    /* find the point furthest from the centroid */
    max_cen_dist = 0.;
    for (i=0; i<n; i++) {
    x = p[2*i];
    y = p[2*i+1];
    dist = sqrt(sqr(x-cen_x) + sqr(y-cen_y));
    if (dist>max_cen_dist) {
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

    make_rotation_matrix((cen_x-p[furthest*2])/max_cen_dist,
             (cen_y-p[furthest*2+1])/max_cen_dist,(double *)m);
    ch_vert = 0;
    vertex = furthest;

    /* HACK: to avoid loop with problem data - PLR */
    ch[1] = -999;

    do {
        ch[ch_vert] = vertex;
        /* Find the ray with the maximum and minimum angle in the new
           coordinate system */
        if (debug) printf("vertex %d\n",vertex);
        max_phi = - BIGNUM;
        min_dist = BIGNUM;
        max_dist = 0.;
        for (i=0; i<n; i++) {
            /* calculate phi, the angle in the new coordinate system */
            x = p[2*i] - p[2*vertex];
            y = p[2*i+1] - p[2*vertex+1];
            xp = m[0][0]*x + m[0][1]*y;
            yp = m[1][0]*x + m[1][1]*y;
            if (debug) printf("\ti %d x %f y %f xp %f yp %f\n", i,x,y,xp,yp);
            if ((dist=sqrt(sqr(x)+sqr(y))) > EPSILON) {
                phi = atan2(yp,xp);
                if (debug) printf("\tphi %f\n", phi);
                if (phi > max_phi-epsilon) {
                    if (full) {
                        /* use the closest point */
                        if (phi > max_phi + epsilon || dist < min_dist) {
                            min_dist = dist;
                            max_phi = phi;
                            i_keep = i;
                        }
                    }
                    else {
                        /* use the furthest point */
                        if (phi > max_phi + epsilon || dist > max_dist) {
                            max_dist = dist;
                            max_phi = phi;
                            i_keep = i;
                        }
                    }
                    if (debug) printf("\t\tmax_phi %f i_keep %d\n", max_phi,i_keep);
                }
            }
        }
        vertex = i_keep;
        xx = cen_x - p[vertex*2];
        yy = cen_y - p[vertex*2+1];
        dist = sqrt(sqr(xx) + sqr(yy));
        make_rotation_matrix(xx/dist,yy/dist,(double *)m);
        ch_vert++;
    /*
    } while (vertex != ch[0]);
    */
    /* HACK: to avoid loop with problem data - PLR */
    } while ((vertex != ch[0]) && (vertex != ch[1]));

    *p_nch = ch_vert;
}

centroid(pts,n,p_cen_x, p_cen_y)
/*
Determines the centroid of the set of n points in the n*2 array of
doubles pts and returns its x-coordinate and y-coordinate in p_cen_x
and p_cen_y respectively.
*/
double *pts, *p_cen_x, *p_cen_y;
int n;
{
    double sumx, sumy;
    int i;

    sumx = sumy = 0;
    for (i=0; i<n; i++) {
        sumx += pts[2*i];
        sumy += pts[2*i+1];
    }
    *p_cen_x = sumx / n;
    *p_cen_y = sumy / n;
}


make_rotation_matrix(cos_theta,sin_theta,m)
double cos_theta, sin_theta, *m;
/*
Given a unit vector (vx,vy) system, finds the matrix m that
corresponds to rotating the coordinate system so that its x-axis
lies in the direction defined by (vx,vy).
*/
{
    *m = cos_theta;
    *(m+1) = sin_theta;
    *(m+2) = - sin_theta;
    *(m+3) = cos_theta;
}



print_subset(int flag_array[MAX_POINTS], int no_pixels, int count)
{
    int i;

    if (count > 2) {
        /*for (i=0;i<no_pixels;i++) {
        printf(" %d", flag_array[i]);
        } */
        printf("%d pixels in final derived shape \n", count);
    }
}

double rn()
/*define constants for random number generator*/
#define p           30269
#define q           30307
#define r           30323

{
  /*Returns a random sample drawn from the interval 0.0 - 1.0*/
  /*rs1, rs2, rs3 are the 3 seeds which must be declared as  */
  /*integers, and given initial values.                      */
  /*Generate the next number in the sequence.                */

  static double z;

  rs1 = abs(rs1 * 171 % p);
  rs2 = abs(rs2 * 172 % q);
  rs3 = abs(rs3 * 170 % r);
  /*reduce these to a sample on 0..3*/
  z = (double)rs1 / p + (double)rs2 / q + (double)rs3 / r;
  /*take the fractional part - this is uniformly distributed on 0..1*/
  return (z - (long)z);
}/*rn*/

generate_subset(int no_pixels, int flag_array[MAX_POINTS], int *count)


{
    int i = 0;


    *count = 0;

    for (i=0;i<no_pixels;i++) {
        if (rn() < 0.5) {
            flag_array[i] =1;
            *count = *count + 1;
        }
        else flag_array[i] = 0;
    }

}


int repair(int flag_array[MAX_POINTS], float *error, int *count) {

/* works out whether the derived shape is convex or not*/
/* if it is not convex, it is replaced with its convex hull */
/* any shape will be rejected out of hand if it is not a true polygon */
/* i.e. if it has less than 3 vertices, or 3 vertices in a line */

    /* subset of input pixel data */
    int no_pixels3;
    int x3[MAX_POINTS];
    int y3[MAX_POINTS];

    /* convex hull data for new polygon */
    int no_pixels4;
    int x4[MAX_POINTS];
    int y4[MAX_POINTS];

    double area_new, area_xor2, area_xor3;
    int j, k;
    int legal = 1;

    /*printf("best error at i = %d is %g\n", i, best_error);*/

    /* for valid polygons only (vertices >= 3*/
    /* copy count cordinates from x, y indicated by '1'*/
    /* in flag-array to x3,y3 */
    if (*count > 2) {
        /*print_subset(flag_array,no_pixels, count);*/
        k=0;
        for (j=0;j<no_pixels;j++) {
            if (flag_array[j] == 1) {
                x3[k] = x[j];
                y3[k] = y[j];
                k++;
            }
        }
        no_pixels3 = k;

        /* compute convex hull of new polygon */
        convex_hull(x3,y3,no_pixels3,x4,y4,&no_pixels4);

        if (no_pixels4 > 2) {
            /* create a polygon structure with the convex hull */
            gpc_create_polygon(x4,y4,no_pixels4,&poly4);

            /* write back convex hull as new poly */
            for (j=0;j<no_pixels;j++) {
                flag_array[j] = 0;
            }

            /*THIS COULD POSSIBLY BE RECODED MORE EFFICIENTLY*/
            /*IT SEARCHES THE WHOLE LIST OF PIXELS IN THE ORIGINAL POLY */
            /*TO MATCH THOSE IN THE DERIVED POLY*/
            *count = no_pixels4;
            for (j=0;j<no_pixels4;j++) {
                k = -1;
                do {
                    k++;
                    if (k >= no_pixels) exit(1);
                } while ((x[k] != x4[j])||(y[k] != y4[j]));
                flag_array[k] = 1;
            }

            /* XOR the new poly area with original poly area*/
            gpc_polygon_clip(GPC_XOR, &poly1, &poly4, &result3);
            area_xor3 = poly_area(&result3);
            *error = (float) area_xor3;
            /*printf("xored area of new and orignal %f \n", area_xor3);*/

            /* minimize this XORed area, provided perfect convexity (= 1) */
            /* keep running total of best so far in a flag array, plus error */

            if (area_xor3 < best_error) {
                gpc_polygon_clip(GPC_INT, &poly1, &poly4, &result4);
                area_int = poly_area(&result4);
                gpc_free_polygon(&result4);

                gpc_polygon_clip(GPC_DIFF, &poly4, &poly1, &result_intrusion);
                area_intrusion = poly_area(&result_intrusion);
                /*
                fpp=fopen("iiii","w");
                gpc_write_polygon2(fpp, &result_intrusion);
                fclose(fpp);
                */
                gpc_free_polygon(&result_intrusion);

                gpc_polygon_clip(GPC_DIFF, &poly1, &poly4, &result_protrusion);
                area_protrusion = poly_area(&result_protrusion);
                /*
                fpp=fopen("pppp","w");
                gpc_write_polygon2(fpp, &result_protrusion);
                fclose(fpp);
                */
                gpc_free_polygon(&result_protrusion);

                /* calculates area of fitted convex hull */
                area_ch2 = poly_area(&poly4);

                best_error = area_xor3;
                /* printf("improved best error is %g\n", best_error); */
                for (j=0;j<no_pixels;j++) best_points[j] = flag_array[j];
                bits_set = *count;

            }
            gpc_free_polygon(&result3);
            gpc_free_polygon(&poly4);
        }
        else legal = 0;
    }
    else legal = 0;
    return legal;
}

crossover(int parent2, int subset_points[MAX_POINTS], int *count)
{

    /*simple one-point crossover */
    int i, position;

    position = (int)(rn()*no_pixels);

    *count = 0;

    for (i=position;i<no_pixels;i++) {
        subset_points[i] = population[parent2][i];
    }

    for (i=0;i<no_pixels;i++) {
        if (subset_points[i] == 1) *count = *count + 1;
    }

}

mutation(int subset_points[MAX_POINTS], int *count)
{
    int position;
    position = (int)(rn()*no_pixels);

    if (subset_points[position] == 0)
    {
        subset_points[position] = 1;
        *count = *count + 1;
    }
    else {
        subset_points[position] = 0;
        *count = *count - 1;
    }
}

int noduplicates(error,count)
{
    int i =0;
    int flag = 1;

    do {
        if ((error == xor_error[i])&&(count == bit_set_count[i]))
            flag = 0;
        i++;
    }while ((flag == 1)&&(i<POPSIZE));
    return flag;
}

ga()
{
    int j, k, count1, legal, weaker;
    int count;
    int subset_points[MAX_POINTS], temp_points[MAX_POINTS];
    float error;
    int generation;
    int attempts;
    int parent2;

    for (generation=0;generation<NUMBER_GENERATIONS;generation++)
    {
        printf(" generation %d\n", generation);
        for (j=0;j<POPSIZE;j++) {

            /* select second parent */
            do {
                parent2 = (int)(rn()*POPSIZE);
            } while (parent2 == j);


            attempts = 0;

            /* crossover and mutation */
            do {
                attempts++;
                count = bit_set_count[j];
                for (k=0;k<no_pixels;k++) {
                    subset_points[k] = population[j][k];
                }/*endFOR*/

                crossover(parent2, subset_points, &count);

                for (count1=0;count1<NO_MUTATIONS;count1++)
                    mutation(subset_points, &count);




                for (k=0;k<no_pixels;k++) {
                     temp_points[k] = subset_points[k];
                }

                legal = repair(subset_points, &error, &count);
            } while ((!legal)&&(attempts<100));
            if (xor_error[parent2] > xor_error[j]) 
                weaker = parent2; 
            else weaker = j; 
            if ((error < xor_error[weaker])&&(attempts<100) 
                &&(noduplicates(error,count)))
            {
                for (k=0;k<no_pixels;k++) {
                    population[weaker][k] = subset_points[k]; 
                }/*endFOR*/
                xor_error[weaker] = error; 
                bit_set_count[weaker] = count; 
            }/*endIF*/
        }/*endFOR*/
    }/*endFOR*/
}/*endGA*/


generate_population()
{
    long i;
    int j, k, legal;
    int subset_points[MAX_POINTS], temp_points[MAX_POINTS];
    float error;
    int count;

    for (j=0;j<POPSIZE;j++) {
        do {
            generate_subset(no_pixels, subset_points,&count);
            for (k=0;k<no_pixels;k++) {
                    temp_points[k] = subset_points[k];
            }
            legal = repair(subset_points, &error, &count);
        } while (!legal);
        xor_error[j] = error;
        bit_set_count[j] = count;
    for (k=0;k<no_pixels;k++) {
            population[j][k] = temp_points[k];
        }

    }
}

/* checks if the polygon is simple or not */

double AreaTriangle2();
typedef struct point { double x, y; } POINT;
POINT polygon[MAX_POINTS];

test_simple(x4,y4,no_pixels4)
int x4[],y4[],no_pixels4;
{
    int i,j,jj;
    int flag = FALSE;

    for (i=0; i<no_pixels4; i++) {
        polygon[i].x = x4[i];
        polygon[i].y = y4[i];
    }

    for (i=0; (i<no_pixels4-1) && !flag; i++)
        for (j=i+1; (j<no_pixels4) && !flag; j++) {

            jj = (j+1)%no_pixels4;

            flag = doIntersectProper(polygon[i],polygon[i+1],polygon[j],polygon[jj]);

            if (flag) {
                printf("pixel\n");
                printf("list: 0\n");
                printf("%.0f %.0f  %.0f %.0f -1 0\n",
                    polygon[i].x,polygon[i].y,polygon[i+1].x,polygon[i+1].y);
                printf("list: 0\n");
                printf("%.0f %.0f  %.0f %.0f -1 -1\n",
                    polygon[j].x,polygon[j].y,polygon[jj].x,polygon[jj].y);
                printf("indices: %d %d  %d %d\n",i,i+1,j,jj);
                printf("SSS %d\n",no_pixels4);
            }
        }

    if (flag) {
        printf("not-simple\n");
        printf("pixel\nlist: 1\n");
        for (i=0; i<no_pixels4; i++)
            printf("%d %d\n",x4[i],y4[i]);
        printf("-1 -1\n");
        exit(-1);
    }

    return(!flag);
}

// Test whether lines intersect properly (i.e. lines cross completely).
// (from 'Computational Geometry in C', J. O'Rourke, ode 1.6)
doIntersectProper(POINT a, POINT b, POINT c, POINT d)
{
    if ((isColinear(a,b,c))||(isColinear(a,b,d))||(isColinear(c,d,a))||(isColinear(c,d,b)))
        return FALSE;
    else
        return (isLeft(a,b,c)^isLeft(a,b,d)) && (isLeft(c,d,a)^isLeft(c,d,b));
}

isBetween(POINT a, POINT b, POINT c)
{
    if (isColinear(a,b,c)==FALSE)
        return FALSE;
    // If line a-b not vertical check betweenness on x else on y
    if (a.x != b.x)
        return ((a.x<=c.x)&&(c.x<=b.x)) || ((a.x>=c.x)&&(c.x>=b.x));
    else
        return ((a.y<=c.y)&&(c.y<=b.y)) || ((a.y>=c.y)&&(c.y>=b.y));
}

isColinear(POINT p1, POINT p2, POINT p3)
{
    if (AreaTriangle2(p1,p2,p3)==0)
        return TRUE;
    else
        return FALSE;
}

isLeft(POINT p1, POINT p2, POINT p3)
{
    if (AreaTriangle2(p1,p2,p3)>0)
        return TRUE;
    else
        return FALSE;
}

double AreaTriangle2(POINT a, POINT b, POINT c)
{
    // returns twice the signed area of the triangle,m positive is ccw from a to b to c
    return  a.x*b.y-a.y*b.x +
            a.y*c.x-a.x*c.y +
            b.x*c.y-c.x*b.y;
}

/*#####################################*/

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

/* generates a uniform random value between 0 & 1 */
double ran3()
{
    static int inext,inextp;
    static long ma[56];
    static int iff=0;
    long mj,mk;
    int i,ii,k;
    static idum = -1;

    if (idum == -1)
        idum = random_seed();

    if (idum < 0 || iff == 0) {
        iff=1;
        mj=MSEED-(idum < 0 ? -idum : idum);
        mj %= MBIG;
        ma[55]=mj;
        mk=1;
        for (i=1;i<=54;i++) {
            ii=(21*i) % 55;
            ma[ii]=mk;
            mk=mj-mk;
            if (mk < MZ) mk += MBIG;
            mj=ma[ii];
        }
        for (k=1;k<=4;k++)
            for (i=1;i<=55;i++) {
                ma[i] -= ma[1+(i+30) % 55];
                if (ma[i] < MZ) ma[i] += MBIG;
            }
        inext=0;
        inextp=31;
        idum=1;
    }
    if (++inext == 56) inext=1;
    if (++inextp == 56) inextp=1;
    mj=ma[inext]-ma[inextp];
    if (mj < MZ) mj += MBIG;
    ma[inext]=mj;
    return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

/* ######## ADDED Paul Rosin 5/94 ######## */
/* hack routine to get initial random number to seed the random number generator
 * uses the date function and combines the second & minute digits
 *
 * returns a random integer
 */
int random_seed()
{
    char line[1000];
    FILE *fp;
    int c,sum;

    system("date > /tmp/tmpdate");
    if ((fp=fopen("/tmp/tmpdate","r")) == NULL) {
        fprintf(stderr,"can't open %s\n","/tmp/tmpdate");
        exit(-1);
    }
    fgets(line,1000,fp);
    sum = 0;
    /* seconds */
    c = (unsigned char)line[18] - '0';
    sum = sum * 10 + c;
    c = (unsigned char)line[17] - '0';
    sum = sum * 10 + c;
    /* minutes */
    c = (unsigned char)line[15] - '0';
    sum = sum * 10 + c;
    c = (unsigned char)line[14] - '0';
    sum = sum * 10 + c;

    fclose(fp);

    return(sum);
}

