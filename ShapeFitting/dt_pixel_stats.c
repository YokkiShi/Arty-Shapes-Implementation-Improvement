/*    input 2 (dense) pixel lists
 *    compute distance transform of 1st list
 *    print statistics of distances at positions of 2nd list
 *
 *    Paul Rosin
 *    Cardiff University
 *    March 2008
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef FALSE
# define FALSE 0
# define TRUE (!FALSE)
#endif

#define MAX_PIXELS 10000
#define MAX_SIZE 1024

#include "pgmio.h"

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define ABS(x)   (((x)<0.0) ? (-(x)): (x))
#define SQR(x)   ((x)*(x))

unsigned char image[MAX_SIZE][MAX_SIZE];
double dt[MAX_SIZE][MAX_SIZE];

int height,width,depth;

int no_pixels1;
int xdata1[MAX_PIXELS],ydata1[MAX_PIXELS];
int no_pixels2;
int xdata2[MAX_PIXELS],ydata2[MAX_PIXELS];

double errors[MAX_PIXELS];

main(argc,argv)
int argc;
char **argv;
{
    int i,j,x,y;
    FILE *fp_pix;
    char file_type[50];
    int minX1,maxX1,minX2,maxX2;
    int minY1,maxY1,minY2,maxY2;
    double mean,sd,skew,kurtosis,adev,var;
    double area;

    if (argc != 3) {
        printf("usage: %s model-pixel-list1 region-pixel-list2\n",argv[0]);
        exit(-1);
    }

    // ------------- read in first shape -------------

    if ((fp_pix=fopen(argv[1],"r")) == NULL){
        printf("cant open %s\n",argv[1]);
        exit(-1);
    }
    /* read magic word for format of file */
    fscanf(fp_pix,"%s\n",file_type);
    j = strcmp(file_type,"pixel");
    if (j != 0){
        printf("not link data file - aborting\n");
        exit(-1);
    }
    read_link_data(fp_pix,xdata1,ydata1,&no_pixels1);
    fclose(fp_pix);

    // find range of data
    minX1 = minY1 = +99999;
    maxX1 = maxY1 = -99999;
    for (i = 0; i < no_pixels1; i++) {
        minX1 = MIN(xdata1[i],minX1);
        minY1 = MIN(ydata1[i],minY1);
        maxX1 = MAX(xdata1[i],maxX1);
        maxY1 = MAX(ydata1[i],maxY1);
    }

    // ------------- read in second shape -------------

    if ((fp_pix=fopen(argv[2],"r")) == NULL){
        printf("cant open %s\n",argv[2]);
        exit(-1);
    }
    /* read magic word for format of file */
    fscanf(fp_pix,"%s\n",file_type);
    j = strcmp(file_type,"pixel");
    if (j != 0){
        printf("not link data file - aborting\n");
        exit(-1);
    }
    read_link_data(fp_pix,xdata2,ydata2,&no_pixels2);
    fclose(fp_pix);

    // find range of data
    minX2 = minY2 = +99999;
    maxX2 = maxY2 = -99999;
    for (i = 0; i < no_pixels2; i++) {
        minX2 = MIN(xdata2[i],minX2);
        minY2 = MIN(ydata2[i],minY2);
        maxX2 = MAX(xdata2[i],maxX2);
        maxY2 = MAX(ydata2[i],maxY2);
    }

    minX1 = MIN(minX1,minX2);
    minY1 = MIN(minY1,minY2);
    maxX1 = MAX(maxX1,maxX2);
    maxY1 = MAX(maxY1,maxY2);

    // hack
    minX1--; minY1--;
    maxX1 += 2; maxY1 += 2;
    // translate both curves (by same amount)
    for (i = 0; i < no_pixels1; i++) {
        xdata1[i] -= minX1;
        ydata1[i] -= minY1;
    }
    for (i = 0; i < no_pixels2; i++) {
        xdata2[i] -= minX1;
        ydata2[i] -= minY1;
    }

    // ------------- compute distances -------------

    width = maxX1 - minX1;
    height = maxY1 - minY1;

    for (y = 0; y < height; y++)
        for (x = 0; x < width; x++)
            image[x][y] = 255;
    for (i = 0; i < no_pixels1; i++)
        image[xdata1[i]][ydata1[i]] = 0;
    edt(image,dt);

    // distances at 2nd curve (index from 1)
    for (i = 0; i < no_pixels2; i++)
        errors[i+1] = dt[xdata2[i]][ydata2[i]];

    // rescale to normalise by square root of region area
    polyCentroid(xdata2,ydata2,no_pixels2,&area);
    for (i = 0; i < no_pixels2; i++)
        errors[i+1] /= sqrt(area);

    moment(errors,no_pixels2,&mean,&adev,&sd,&var,&skew,&kurtosis);
    printf("%f, %f, %f, %f, ",mean,sd,skew,kurtosis);

    /***
    printf("mean %f\n",mean);
    printf("sd %f\n",sd);
    printf("skew %f\n",skew);
    printf("kurtosis %f\n",kurtosis);
    ***/

    /*****************************
    for (y = 0; y < height; y++)
        for (x = 0; x < width; x++)
            image[x][y] = dt[x][y];
    write_pgm(image,"dt-result",width,height);
    *****************************/
}

// ------------------------ EDT ------------------------

#define MASK_SIZE1 5
int x_offset1[MASK_SIZE1] = { -1,  0,  1, -1,  0 };
int y_offset1[MASK_SIZE1] = { -1, -1, -1,  0,  0 };

#define MASK_SIZE2 4
int x_offset2[MASK_SIZE2] = {  1,  1,  1,  0 };
int y_offset2[MASK_SIZE2] = { -1,  0,  1,  0 };

/* perform Euclidean distance transform */
edt(image1,image2)
unsigned char image1[MAX_SIZE][MAX_SIZE];
double image2[MAX_SIZE][MAX_SIZE];
{
    int i,x,y;
    double d[MASK_SIZE1];
    double min_measure;
    double dist;
    int min_loc;
    int xn,yn;
    int dx,dy;
    int new_x,new_y;
    static int image_dist_x[MAX_SIZE][MAX_SIZE],image_dist_y[MAX_SIZE][MAX_SIZE];

    for (y = 0; y < height; y++)
        for (x = 0; x < width; x++)
            image_dist_x[x][y] = image_dist_y[x][y] = image1[x][y];

    /* propagate distance etc */
    /* forward pass */
    for (y = 1; y < height-1; y++) {
        for (x = 1; x < width-1; x++) {
            for (i = 0; i < MASK_SIZE1; i++) {
                xn = x + x_offset1[i];
                yn = y + y_offset1[i];
                dx = image_dist_x[xn][yn] + ABS(x_offset1[i]);
                dy = image_dist_y[xn][yn] + ABS(y_offset1[i]);
                d[i] = SQR(dx) + SQR(dy);
                d[i] = sqrt(d[i]);
            }

            min_loc = 0;
            min_measure = d[0];
            for (i = 1; i < MASK_SIZE1; i++) {
                if (d[i] < min_measure) {
                    min_measure = d[i];
                    min_loc = i;
                }
            }
            xn = x + x_offset1[min_loc];
            yn = y + y_offset1[min_loc];
            new_x = image_dist_x[xn][yn] + ABS(x_offset1[min_loc]);
            new_y = image_dist_y[xn][yn] + ABS(y_offset1[min_loc]);

            image_dist_x[x][y] = new_x;
            image_dist_y[x][y] = new_y;
        }
    }

    /* backward pass */
    for (y = height-2; y >= 1; y--) {
        for (x = 1; x < width-1; x++) {
            for (i = 0; i < MASK_SIZE1; i++) {
                xn = x - x_offset1[i];
                yn = y - y_offset1[i];
                dx = image_dist_x[xn][yn] + ABS(x_offset1[i]);
                dy = image_dist_y[xn][yn] + ABS(y_offset1[i]);
                d[i] = SQR(dx) + SQR(dy);
                d[i] = sqrt(d[i]);
            }

            min_loc = 0;
            min_measure = d[0];
            for (i = 1; i < MASK_SIZE1; i++) {
                if (d[i] < min_measure) {
                    min_measure = d[i];
                    min_loc = i;
                }
            }
            xn = x - x_offset1[min_loc];
            yn = y - y_offset1[min_loc];
            new_x = image_dist_x[xn][yn] + ABS(x_offset1[min_loc]);
            new_y = image_dist_y[xn][yn] + ABS(y_offset1[min_loc]);

            image_dist_x[x][y] = new_x;
            image_dist_y[x][y] = new_y;
        }
    }

    /* sideways pass */
    for (x = width-2; x >= 1; x--) {
        for (y = 1; y < height-1; y++) {
            for (i = 0; i < MASK_SIZE2; i++) {
                xn = x + x_offset2[i];
                yn = y + y_offset2[i];
                dx = image_dist_x[xn][yn] + ABS(x_offset2[i]);
                dy = image_dist_y[xn][yn] + ABS(y_offset2[i]);
                d[i] = SQR(dx) + SQR(dy);
                d[i] = sqrt(d[i]);
            }

            min_loc = 0;
            min_measure = d[0];
            for (i = 1; i < MASK_SIZE2; i++) {
                if (d[i] < min_measure) {
                    min_measure = d[i];
                    min_loc = i;
                }
            }
            xn = x + x_offset2[min_loc];
            yn = y + y_offset2[min_loc];
            new_x = image_dist_x[xn][yn] + ABS(x_offset2[min_loc]);
            new_y = image_dist_y[xn][yn] + ABS(y_offset2[min_loc]);

            image_dist_x[x][y] = new_x;
            image_dist_y[x][y] = new_y;
        }
    }

    /* calculate measure values */
    for (y = 1; y < height-1; y++) {
        for (x = 1; x < width-1; x++) {
            dist = SQR(image_dist_x[x][y]) + SQR(image_dist_y[x][y]);
            dist = sqrt(dist);
            image2[x][y] = (double)dist;
        }
    }

    /* set borders by copying adjacent values */
    for (x = 0; x < width; x++) {
        image2[x][0] = image2[x][1];
        image2[x][height-1] = image2[x][height-2];
    }
    for (y = 0; y < height; y++) {
        image2[0][y] = image2[1][y];
        image2[width-1][y] = image2[width-2][y];
    }
}

read_link_data(fp,xdata,ydata,no_pixels)
FILE *fp;
int xdata[],ydata[];
int *no_pixels;
{
    char dumstring[50];
    int j;
    int endoffile;

    fscanf(fp,"%s %d\n",dumstring,&j);
    j = -1;
    do{
       j++;
       fscanf(fp,"%d %d\n",&xdata[j],&ydata[j]);
    } while(xdata[j] != -1);
    endoffile = (ydata[j] == -1);
    if (feof(fp) && !(endoffile)) {
        fprintf(stderr,"Incorrectly terminated file - assuming end of list\n");
        endoffile = TRUE;
    }
    *no_pixels = j;
    if (*no_pixels >= MAX_PIXELS) {
        fprintf(stderr,"ERROR: Too many pixels\n");
        exit(-1);
    }
}

moment(data,n,ave,adev,sdev,svar,skew,curt)
int n;
double data[],*ave,*adev,*sdev,*svar,*skew,*curt;
{
    int j;
    double s,p;
    void nrerror();

    if (n <= 1) nrerror("n must be at least 2 in MOMENT");
    s=0.0;
    for (j=1;j<=n;j++) s += data[j];
    *ave=s/n;
    *adev=(*svar)=(*skew)=(*curt)=0.0;
    for (j=1;j<=n;j++) {
        *adev += ABS(s=data[j]-(*ave));
        *svar += (p=s*s);
        *skew += (p *= s);
        *curt += (p *= s);
    }
    *adev /= n;
    *svar /= (n-1);
    *sdev=sqrt(*svar);
    if (*svar) {
        *skew /= (n*(*svar)*(*sdev));
        *curt=(*curt)/(n*(*svar)*(*svar))-3.0;
    } else nrerror("No skew/kurtosis when variance = 0 (in MOMENT)");
}

void nrerror(error_text)
char error_text[];
{
    fprintf(stderr,"Numerical Recipes run-time error...\n");
    fprintf(stderr,"%s\n",error_text);
    fprintf(stderr,"...now exiting to system...\n");
    exit(-1);
}

/*********************************************************************
polyCentroid: Calculates                                     and area
of a polygon, given its vertices (x[0], y[0]) ... (x[n-1], y[n-1]). It
is assumed that the contour is closed, i.e., that the vertex following
(x[n-1], y[n-1]) is (x[0], y[0]).  The algebraic sign of the area is
positive for counterclockwise ordering of vertices in x-y plane;
otherwise negative.

Returned values:  0 for normal execution;  1 if the polygon is
degenerate (number of vertices < 3);  and 2 if area = 0 (and the
centroid is undefined).
**********************************************************************/
int polyCentroid(int x[], int y[], int n, double *area)
{
    register int i, j;
    double ai, atmp = 0, xtmp = 0, ytmp = 0;
    if (n < 3)
        return 1;
    for (i = n - 1, j = 0; j < n; i = j, j++) {
        ai = x[i] * y[j] - x[j] * y[i];
        atmp += ai;
        xtmp += (x[j] + x[i]) * ai;
        ytmp += (y[j] + y[i]) * ai;
    }
    *area = atmp / 2;

    /* !!!!! make area unsigned !!!!! */
    *area = ABS(*area);
}
