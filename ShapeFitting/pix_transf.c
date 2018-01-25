/*  ---------------------------------------------------------------
 *    transform pixel data (integer, float, or curvature)
 *        scale
 *        translate
 *        rotate
 *  ---------------------------------------------------------------
 *    select pixel lists
 *  ---------------------------------------------------------------
 *    convert (& round) pixel_float -> pixel
 *  ---------------------------------------------------------------
 *    modify pixel data (only integer)
 *        subsample
 *        4-way to 8-way connectivity
 *        8-way to 4-way connectivity
 *        convert to 1 pixel separation (8-way connectivity)
 *        average co-ordinate values
 *  ---------------------------------------------------------------
 *    convert pixel data -> super data
 *  ---------------------------------------------------------------
 *
 *    Paul Rosin
 *    Joint Research Centre
 *    Ispra, Italy
 *    March 1995
 *    paul.rosin@jrc.it
 *
 *    converted float -> double
 *    PLR September 2006
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef FALSE
#define FALSE 0
#define TRUE (!FALSE)
#endif

#define ABS(x)       (((x)<0.0)? (-(x)): (x))
#define MIN(a,b)     (((a) < (b)) ? (a) : (b))
#define MAX(a,b)     (((a) > (b)) ? (a) : (b))
#define SQR(x)       ((x)*(x))

/* format used for data in pixel lists */
#define INTEGER 0
#define FLOAT 1
#define CURVATURE 2
#define SUPER 999

#define PI 3.1415926

#define MAX_PIXELS 1000000

int no_pixels,no_new_pixels;
double x[MAX_PIXELS],y[MAX_PIXELS];
double xo[MAX_PIXELS],yo[MAX_PIXELS];
double k[MAX_PIXELS];
int l[MAX_PIXELS];
double sigma;
int xpos,ypos;

double xscale=1,yscale=1;
double xyscale;
double xshift=0,yshift=0;
double theta;
int reverse = FALSE;
int do_xyscale = FALSE;
int clip = FALSE;
int centre_list = TRUE;  // centre each list - turn off for multiple list shapes

int file_data_type;
int file_data_type_output;

int lists[10000];
int no_lists;

main(argc,argv)
int argc;
char *argv[];
{
    FILE *fp_in,*fp_out;
    char file_type[50];
    int a,i,s;
    int endoffile,closed_loop,list_no;
    char *infile=NULL,*outfile=NULL;
    char *list_file=NULL;
    int subsample = FALSE;
    int c4t8 = FALSE;
    int c8t4 = FALSE;
    int space1 = FALSE;
    int average = FALSE;
    int gauss = FALSE;
    int do_transform = FALSE;
    int float2int = FALSE;
    int int2float = FALSE;
    int do_top = FALSE;
    int do_ctr = FALSE;
    double var;
    int orig_file_data_type;
    int make_super = FALSE;
    int Christine = FALSE;
    int qhull = FALSE;
    int make2 = FALSE;
    double g;
    int do_gaps = FALSE;

    for (i = 1; i < argc; i++) {
       if (argv[i][0] == '-')
          switch (argv[i][1]) {
              case 'X':
                   i++;
                   sscanf(argv[i],"%lf",&xscale);
                   printf("X scaling = %f\n",xscale);
                   do_transform = TRUE;
                   break;
              case 'Y':
                   i++;
                   sscanf(argv[i],"%lf",&yscale);
                   printf("Y scaling = %f\n",yscale);
                   do_transform = TRUE;
                   break;
              case 'Z':
                   i++;
                   sscanf(argv[i],"%lf",&xyscale);
                   printf("scale to max = %f\n",xyscale);
                   do_transform = TRUE;
                   do_xyscale = TRUE;
                   break;
              case 'z':
                   clip = TRUE;
                   break;
              case 'D':
                   centre_list = FALSE;
                   break;
              case 'l':
                   i++;
                   list_file = argv[i];
                   break;
              case 'x':
                   i++;
                   sscanf(argv[i],"%lf",&xshift);
                   printf("X shift = %f\n",xshift);
                   do_transform = TRUE;
                   break;
              case 'y':
                   i++;
                   sscanf(argv[i],"%lf",&yshift);
                   printf("Y shift = %f\n",yshift);
                   do_transform = TRUE;
                   break;
              case 'r':
                   i++;
                   sscanf(argv[i],"%lf",&theta);
                   printf("rotation angle = %f\n",theta);
                   do_transform = TRUE;
                   break;
              case 'a':
                   i++;
                   sscanf(argv[i],"%d",&a);
                   printf("averaging in %d window\n",a);
                   average = TRUE;
                   break;
              case 'G':
                   i++;
                   sscanf(argv[i],"%lf",&g);
                   printf("filling gaps > %f\n",g);
                   do_gaps = TRUE;
                   break;
              case 's':
                   i++;
                   sscanf(argv[i],"%d",&s);
                   printf("subsampling 1 in %d point\n",s);
                   subsample = TRUE;
                   break;
              case 'g':
                   i++;
                   sscanf(argv[i],"%lf",&var);
                   printf("adding Gaussian noise (var = %.2f)\n",var);
                   gauss = TRUE;
                   break;
              case 'd':
                   make2 = TRUE;
                   break;
              case 'q':
                   qhull = TRUE;
                   break;
              case 'C':
                   Christine = TRUE;
                   break;
              case 'S':
                   make_super = TRUE;
                   break;
              case 'R':
                   reverse = TRUE;
                   break;
              case 't':
                   do_top = TRUE;
                   break;
              case 'c':
                   do_ctr = TRUE;
                   break;
              case 'I':
                   float2int = TRUE;
                   break;
              case 'F':
                   int2float = TRUE;
                   break;
              case '1':
                   space1 = TRUE;
                   break;
              case '4':
                   c8t4 = TRUE;
                   break;
              case '8':
                   c4t8 = TRUE;
                   break;
              case 'i':
                   i++;
                   infile = argv[i];
                   break;
              case 'o':
                   i++;
                   outfile = argv[i];
                   break;
              default:
                   printf("illegal option -%c\n",argv[i][1]);
                   print_usage(argv[0]);
          }
        else {
           printf("illegal option %c\n",argv[i][0]);
           print_usage(argv[0]);
        }
    }

    theta = theta * PI / 180.0;

    if ((infile == NULL) || (outfile == NULL)) {
       print_usage(argv[0]);
    }

    if ((fp_in=fopen(infile,"r")) == NULL) {
        printf("cant open %s\n",infile);
        exit(-1);
    }

    /* read magic word for format of file */
    fscanf(fp_in,"%s\n",file_type);
    if (strcmp(file_type,"pixel") == 0) {
        /***
        fprintf(stderr,"integer pixel data file\n");
        ***/
        file_data_type = INTEGER;
    }
    else if (strcmp(file_type,"pixel_float") == 0) {
        /***
        fprintf(stderr,"float pixel data file\n");
        ***/
        file_data_type = FLOAT;
    }
    else if (strcmp(file_type,"pixel_curvature") == 0) {
        /***
        fprintf(stderr,"curvature pixel data file\n");
        ***/
        file_data_type = CURVATURE;
    }
    else{
        fprintf(stderr,"input file not pixel, pixel_float or pixel_curvature type - aborting\n");
        exit(-1);
    }

    if (float2int) {
        orig_file_data_type = file_data_type;
        /***
        if (file_data_type != FLOAT) {
            fprintf(stderr,"ERROR: input must be pixel_float\n");
            exit(-1);
        }
        ***/
        file_data_type = INTEGER;
    }
    else if (int2float) {
        orig_file_data_type = file_data_type;
        /***
        if (file_data_type != INTEGER) {
            fprintf(stderr,"ERROR: input must be pixel\n");
            exit(-1);
        }
        ***/
        file_data_type = FLOAT;
    }
    /*** PLR IS THIS NECESSARY? July 2004
    else if (!do_transform) {
        if (file_data_type != INTEGER) {
            fprintf(stderr,"ONLY INT!!!\n");
            exit(-1);
        }
    }
    ***/

    if (list_file != NULL)
        read_list_ids(list_file);

    file_data_type_output = file_data_type;
    if (make_super)
        file_data_type_output = SUPER;

    if ((fp_out=fopen(outfile,"w")) == NULL) {
        printf("cant open %s\n",outfile);
        exit(-1);
    }

    if (!Christine && !qhull)
    switch (file_data_type_output) {
        case INTEGER:
            fprintf(fp_out,"pixel\n");
            break;
        case FLOAT:
            fprintf(fp_out,"pixel_float\n");
            break;
        case CURVATURE:
            fprintf(fp_out,"pixel_curvature\n");
            break;
        case SUPER:
            fprintf(fp_out,"super\n");
            break;
    }

    if (float2int) {
        do {
            file_data_type = FLOAT;
            read_link_data(fp_in,&endoffile,&closed_loop,&list_no);
            if (orig_file_data_type == INTEGER)
                for (i = 0; i < no_pixels; i++) {
                    xo[i] = x[i];
                    yo[i] = y[i];
                }
            /* round */
            else
                for (i = 0; i < no_pixels; i++) {
                    /*
                    xo[i] = x[i] + 0.5;
                    yo[i] = y[i] + 0.5;
                    */
                    xo[i] = x[i] + 0.49;
                    yo[i] = y[i] + 0.49;
                }
            no_new_pixels = no_pixels;
            file_data_type = INTEGER;
            store_link_data(fp_out,endoffile,list_no);
        } while (!endoffile);

        /* may need to add dummy list to terminate file */
        fprintf(fp_out,"list: -1\n-1 -1\n");
    }
    else if (Christine) {
        do {
            read_link_data(fp_in,&endoffile,&closed_loop,&list_no);

            fprintf(fp_out,"%d\n",no_pixels);
            switch (file_data_type) {
                case INTEGER:
                    for (i = 0; i < no_pixels; i++)
                        fprintf(fp_out,"%.0f %.0f\n",x[i],y[i]);
                    break;
                case FLOAT:
                    for (i = 0; i < no_pixels; i++)
                        fprintf(fp_out,"%f %f\n",x[i],y[i]);
                    break;
            }

            fclose(fp_out);
        } while (!endoffile);
        exit(0);
    }
    else if (qhull) {
        do {
            read_link_data(fp_in,&endoffile,&closed_loop,&list_no);
            fprintf(fp_out,"2 voronoi vertices\n",no_pixels);
            fprintf(fp_out,"%d\n",no_pixels);
            for (i = 0; i < no_pixels; i++)
                fprintf(fp_out,"%.0f %.0f\n",x[i],y[i]);
            fclose(fp_out);
        } while (!endoffile);
        exit(0);
    }
    else if (make2) {
        do {
            read_link_data(fp_in,&endoffile,&closed_loop,&list_no);
            fprintf(fp_out,"list: 0\n");
            for (i = 1; i < no_pixels; i++) {
                fprintf(fp_out,"%.0f %.0f\n",x[i-1],y[i-1]);
                fprintf(fp_out,"%.0f %.0f\n",x[i],y[i]);
            }
            fprintf(fp_out,"-1 -1\n");
            fclose(fp_out);
        } while (!endoffile);
        exit(0);
    }
    else if (subsample) {
        do {
            read_link_data(fp_in,&endoffile,&closed_loop,&list_no);

            if (no_pixels > s) {
                subs(s);
                store_link_data(fp_out,endoffile,list_no);
            }
            /*
            else {
                subs(1);
                store_link_data(fp_out,endoffile,list_no);
            }
            */
        } while (!endoffile);

        /* may need to add dummy list to terminate file */
        fprintf(fp_out,"list: -1\n-1 -1\n");
    }
    else {
        do {
            read_link_data(fp_in,&endoffile,&closed_loop,&list_no);
    
            if (list_file != NULL)
                if (!list_member(list_no))
                    continue;

            if (c4t8)
                convert4to8();
            else if (c8t4)
                convert8to4();
            else if (c8t4)
                convert8to4();
            else if (space1)
                do_space1();
            else if (average)
                avg(a,closed_loop);
            else if (gauss)
                add_gauss(var);
            else if (do_transform)
                transform();
            else if (do_top) {
                find_min(&xshift,&yshift);
                printf("calculated shift: %f %f\n",xshift,yshift);
                xshift = -xshift; yshift = -yshift;
                transform();
            }
            else if (do_ctr) {
                centroid(&xshift,&yshift);
                printf("centroid: %f %f\n",xshift,yshift);
                /* arbitrarily put centroid at (300,300) */
                xshift = 300-xshift; yshift = 300-yshift;
                transform();
            }
            else if (do_gaps)
                fillgaps(g,closed_loop);
            /* default: copy input to output */
            else {
                for (i = 0; i < no_pixels; i++) {
                    xo[i] = x[i];
                    yo[i] = y[i];
                }
                no_new_pixels = no_pixels;
printf("copied %d items\n",no_pixels);
            }
    
            store_link_data(fp_out,endoffile,list_no);
        } while (!endoffile);
    }

    fclose(fp_in);
    fclose(fp_out);
}

/* check if given ID is in list of IDs */
int list_member(id)
int id;
{
    int i;
    int found = FALSE;

    for (i = 0; i < no_lists && !found; i++) {
        if (id == lists[i])
            found = TRUE;
    }
    return(found);
}

/* read in a list of pixel list ID's */
read_list_ids(infile)
char *infile;
{
    FILE *fp;
    int j = 0;

    if ((fp=fopen(infile,"r")) == NULL) {
        printf("cant open %s\n",infile);
        exit(-1);
    }

    while (fscanf(fp,"%d",&lists[j]) != EOF) {
       j++;
    }
    no_lists = j;

    fclose(fp);
}

transform()
{
    int i;
    double temp;
    double xt,yt;
    double sine,cosine;
    double xc,yc;

    /* rescale to make max(x,y) = xyscale after shifting data to XY axes
     * or clip: rescale only if max > xyscale
     */
    if (do_xyscale) {
        double x1,y1,x2,y2;

        find_min(&x1,&y1);
        find_max(&x2,&y2);
        temp = MAX(x2-x1,y2-y1);

        if (clip) {
            if (temp > xyscale)
                xyscale /= temp;
            else
                xyscale = 1;
        }
        else
            xyscale /= temp;

        for (i = 0; i < no_pixels; i++) {
            xo[i] = (x[i]-x1) * xyscale;
            yo[i] = (y[i]-y1) * xyscale;
        }
        no_new_pixels = no_pixels;
        return;
    }
    /* rescale to make max(x,y) = xyscale */
    /*
    if (do_xyscale) {
        find_max(&xt,&yt);
        temp = MAX(xt,yt);
        xyscale /= temp;
        for (i = 0; i < no_pixels; i++) {
            xo[i] = x[i] * xyscale;
            yo[i] = y[i] * xyscale;
        }
        no_new_pixels = no_pixels;
        return;
    }
    */
    
    /* translate */
    xc = yc = 0;
    for (i = 0; i < no_pixels; i++) {
        xo[i] = x[i] + xshift;
        yo[i] = y[i] + yshift;
        xc += x[i]; yc += y[i];
    }
    xc /= no_pixels; yc /= no_pixels;

    /* scale */
    for (i = 0; i < no_pixels; i++) {
        xo[i] *= xscale;
        yo[i] *= yscale;
    }

    /* rotate */
    if (theta != 0) {
        sine = sin(-theta);
        cosine = cos(-theta);
        for (i = 0; i < no_pixels; i++) {
            xt = xo[i];
            yt = yo[i];

            if (centre_list) {
                xt -= xc;
                yt -= yc;
            }

            temp = xt*cosine - yt*sine;
            xo[i] = temp;
            temp = xt*sine + yt*cosine;
            yo[i] = temp;

            if (centre_list) {
                xo[i] += xc;
                yo[i] += yc;
            }

        }
    }
    no_new_pixels = no_pixels;
}

find_min(xt,yt)
double *xt,*yt;
{
    int i;
    double xtt,ytt;

    xtt = x[0]; ytt = y[0];

    for (i = 1; i < no_pixels; i++) {
        xtt = MIN(xtt,x[i]);
        ytt = MIN(ytt,y[i]);
    }
    *xt = xtt; *yt = ytt;
}

find_max(xt,yt)
double *xt,*yt;
{
    int i;
    double xtt,ytt;

    xtt = x[0]; ytt = y[0];

    for (i = 1; i < no_pixels; i++) {
        xtt = MAX(xtt,x[i]);
        ytt = MAX(ytt,y[i]);
    }
    *xt = xtt; *yt = ytt;
}

centroid(xt,yt)
double *xt,*yt;
{
    int i;
    double sx,sy;

    sx = sy = 0;

    for (i = 0; i < no_pixels; i++) {
        sx += x[i];
        sy += y[i];
    }
    *xt = sx/no_pixels; *yt = sy/no_pixels;
}

/* add Gaussian noise to each point */
add_gauss(var)
double var;
{
    double random_gauss();
    int loop1;

    for (loop1 = 0; loop1 < no_pixels; loop1++) {
        xo[loop1] = x[loop1] + random_gauss(0.0,var);
        yo[loop1] = y[loop1] + random_gauss(0.0,var);
    }
    no_new_pixels = no_pixels;
}


/* convert 4-way connected curve to 8-way connected */
convert8to4()
{
    int loop1;
    int xp,yp;
    int xn,yn;
    int xd,yd;

    xp = x[0]; yp = y[0];
    xo[0] = xp; yo[0] = yp;
    no_new_pixels = 1;
    for (loop1 = 1; loop1 < no_pixels; loop1++) {
        xn = x[loop1]; yn = y[loop1];
        xd = xn - xp; yd = yn - yp;
        if ((xd != 0) && (yd != 0)) {
            xo[no_new_pixels] = xp+xd; yo[no_new_pixels] = yp;
            no_new_pixels++;
        }
        xo[no_new_pixels] = xn; yo[no_new_pixels] = yn;
        no_new_pixels++;
        xp = xn; yp = yn;
    }
}

/* convert 4-way connected curve to 8-way connected */
convert4to8()
{
    int loop1;
    int x_prev,y_prev;
    int x_cur,y_cur;
    int x_next,y_next;
    int xd1,yd1;
    int xd2,yd2;
    int kept_prev;

    x_prev = x[0]; y_prev = y[0];
    xo[0] = x_prev; yo[0] = y_prev;
    no_new_pixels = 1;
    kept_prev = TRUE;
    for (loop1 = 1; loop1 < no_pixels-1; loop1++) {
        x_cur = x[loop1]; y_cur = y[loop1];
        x_next = x[loop1+1]; y_next = y[loop1+1];
        xd1 = x_cur - x_prev; yd1 = y_cur - y_prev;
        xd2 = x_next - x_cur; yd2 = y_next - y_cur;
        if (kept_prev) {
            if (((xd1 != 0) && (yd1 == 0)) &&
                ((xd2 == 0) && (yd2 != 0)))
            {
                /* delete pixel */
                kept_prev = FALSE;
            }
            else
            if (((xd1 == 0) && (yd1 != 0)) &&
                ((xd2 != 0) && (yd2 == 0)))
            {
                /* delete pixel */
                kept_prev = FALSE;
            }
            else {
                xo[no_new_pixels] = x_cur; yo[no_new_pixels] = y_cur;
                no_new_pixels++;
                kept_prev = TRUE;
            }
        }
        else {
            xo[no_new_pixels] = x_cur; yo[no_new_pixels] = y_cur;
            no_new_pixels++;
            kept_prev = TRUE;
        }
        x_prev = x_cur; y_prev = y_cur;
    }
    xo[no_new_pixels] = x_next; yo[no_new_pixels] = y_next;
    no_new_pixels++;

}

/* subsample curve */
subs(s)
int s;
{
    int j,k=0;

    for (j = 0; j < no_pixels-1; j += s) {
        xo[k] = x[j]; yo[k] = y[j];
        k++;
    }
    xo[k] = x[no_pixels-1]; yo[k] = y[no_pixels-1];
    k++;
    no_new_pixels = k;
}

/* average curve */
avg(w,closed)
int w,closed;
{
    int i,j,ii;
    int ax,ay;

    if (closed) {
        for (j = 0; j < no_pixels; j++) {
            ax = ay = 0;
            for (i = 0; i < w; i++) {
                ii = (j+i) % no_pixels;
                ax += x[ii];
                ay += y[ii];
            }
            xo[j] = ax / w;
            yo[j] = ay / w;
        }
        no_new_pixels = no_pixels;
    }
    else {
        for (j = 0; j < no_pixels-w; j++) {
            ax = ay = 0;
            for (i = 0; i < w; i++) {
                ax += x[j+i];
                ay += y[j+i];
            }
            xo[j] = ax / w;
            yo[j] = ay / w;
        }
        no_new_pixels = no_pixels-w;
    }
}

/* fill gaps in curve by interpolation */
fillgaps(g,closed)
double g;
int closed;
{
    int i,j,ii;
    double dist,dx,dy;
    int tmp_no_pixels;

    x[no_pixels] = x[0];
    y[no_pixels] = y[0];
    if (closed)
        tmp_no_pixels = no_pixels + 1;
    else
        tmp_no_pixels = no_pixels;

    i = 0;
    for (j = 0; j < tmp_no_pixels-1; j++) {
        dx = x[j] - x[j+1];
        dy = y[j] - y[j+1];
        dist = sqrt(SQR(dx)+SQR(dy));
        xo[i] = x[j];
        yo[i] = y[j];
        i++;

        // interpolate extra point in middle
        if (dist > g) {
            xo[i] = (x[j] + x[j+1]) / 2.0;
            yo[i] = (y[j] + y[j+1]) / 2.0;
            i++;
        }
    }
    if (!closed) {
        xo[i] = x[tmp_no_pixels-1];
        yo[i] = y[tmp_no_pixels-1];
        i++;
    }

    no_new_pixels = i;

    printf("%d points added\n",no_new_pixels-no_pixels);
}

do_space1()
{
    int loop1;

    no_new_pixels = 0;
    moveto_image((int)x[0],(int)y[0]);
    for (loop1 = 1; loop1 < no_pixels; loop1++)
        lineto_image((int)x[loop1],(int)y[loop1]);
}

lineto_image(xd,yd)
int xd,yd;
{
     int xnew,ynew;
     
     xnew = xd;
     ynew = yd;
     lineby_image(xnew-xpos,ynew-ypos);
}

moveto_image(xd,yd)
int xd,yd;
{
    xpos = xd;
    ypos = yd;
}

lineby_image(xd,yd)
int xd,yd;
{
    int loop1,loop2,i,xstep,ystep;

    if (xd == 0) {
        xstep = 0;
        xd = 1;
    }
    else if (xd < 0) {
        xstep = -1;
        xd = -xd;
    }
    else
       xstep = 1;
    if (yd == 0) {
        ystep = 0;
        yd = 1;
    }
    else if (yd < 0) {
        ystep = -1;
        yd = -yd;
    }
    else
        ystep = 1;
    draw_point(xpos,ypos);
    i = 0;
    for (loop1=1; loop1<=xd; loop1++) {
        for (loop2=1; loop2<=yd; loop2++) {
            i++;
            if (i >= xd) {
                i = 0;
                ypos = ypos + ystep;
                draw_point(xpos,ypos);
            }
        }
        xpos = xpos + xstep;
        draw_point(xpos,ypos);
    }
}

draw_point(xpos,ypos)
int xpos,ypos;
{
    xo[no_new_pixels] = xpos;
    yo[no_new_pixels] = ypos;
    no_new_pixels++;
}

read_link_data(fp,endoffile,closed_loop,list_no)
FILE *fp;
int *endoffile,*closed_loop,*list_no;
{
    switch (file_data_type) {
        case INTEGER:
            read_integer_data(fp,endoffile,closed_loop,list_no);
            break;
        case FLOAT:
            read_float_data(fp,endoffile,closed_loop,list_no);
            break;
        case CURVATURE:
            read_curvature_data(fp,endoffile,closed_loop,list_no);
            break;
    }
}

read_integer_data(fp,endoffile,closed_loop,list_no)
FILE *fp;
int *endoffile,*closed_loop,*list_no;
{
    char dumstring[50];
    int j;
    int tx,ty;
    int xd,yd;

    fscanf(fp,"%s %d\n",dumstring,list_no);
    j = -1;
    do {
       j++;
       fscanf(fp,"%d %d",&tx,&ty);
       x[j] = tx; y[j] = ty;
    } while(x[j] != -1);
    *endoffile = (y[j] == -1);
    if (feof(fp) && !(*endoffile)) {
        fprintf(stderr,"Incorrectly terminated file - assuming end of list\n");
        *endoffile = TRUE;
    }
    xd = ABS(x[0] - x[j-1]);
    yd = ABS(y[0] - y[j-1]);
    *closed_loop = ((xd <= 1) && (yd <= 1));
    no_pixels = j;
}

read_float_data(fp,endoffile,closed_loop,list_no)
FILE *fp;
int *endoffile,*closed_loop,*list_no;
{
    char dumstring[50];
    int j;
    double tx,ty;
    int xd,yd;
	double tmp;

    fscanf(fp,"%s %d\n",dumstring,list_no);
    j = -1;
    do {
       j++;
       fscanf(fp,"%lf %lf",&tx,&ty);
       x[j] = tx; y[j] = ty;
    } while(x[j] != -1);
    *endoffile = (y[j] == -1);
    if (feof(fp) && !(*endoffile)) {
        fprintf(stderr,"Incorrectly terminated file - assuming end of list\n");
        *endoffile = TRUE;
    }
    xd = ABS(x[0] - x[j-1]);
    yd = ABS(y[0] - y[j-1]);
    *closed_loop = ((xd <= 1) && (yd <= 1));
    no_pixels = j;
}

read_curvature_data(fp,endoffile,closed_loop,list_no)
FILE *fp;
int *endoffile,*closed_loop,*list_no;
{
    char dumstring[50];
    int j;
    double tx,ty,K;
    int label;
    int xd,yd;

    fscanf(fp,"%s %d\n",dumstring,list_no);
    fscanf(fp,"%s %lf\n",dumstring,&sigma);
    j = -1;
    do {
       j++;
       fscanf(fp,"%lf %lf %lf %d",&tx,&ty,&K,&label);
       x[j] = tx; y[j] = ty; k[j] = K; l[j] = label;
    } while(x[j] != -1);
    *endoffile = (y[j] == -1);
    if (feof(fp) && !(*endoffile)) {
        fprintf(stderr,"Incorrectly terminated file - assuming end of list\n");
        *endoffile = TRUE;
    }
    xd = ABS(x[0] - x[j-1]);
    yd = ABS(y[0] - y[j-1]);
    *closed_loop = ((xd <= 1) && (yd <= 1));
    no_pixels = j;
}

store_link_data(fp,endoffile,list_no)
FILE *fp;
int endoffile,list_no;
{
    int j;

    switch (file_data_type_output) {
        case INTEGER:
            store_integer_data(fp,endoffile,list_no);
            break;
        case FLOAT:
            store_float_data(fp,endoffile,list_no);
            break;
        case CURVATURE:
            store_curvature_data(fp,endoffile,list_no);
            break;
        case SUPER:
            store_super_data(fp,endoffile,list_no);
            break;
    }
}

store_super_data(fp,endoffile,list_no)
FILE *fp;
int endoffile,list_no;
{
    int j;

    fprintf(fp,"list: %d\n",list_no);
    for (j = 0; j < no_new_pixels-1; j++)
        fprintf(fp,"line: 0.0 %.0f %.0f   %.0f %.0f\n",xo[j],yo[j],xo[j+1],yo[j+1]);
    fprintf(fp,"endl:\n");

    if (endoffile)
        fprintf(fp,"endf:\n");
}

store_integer_data(fp,endoffile,list_no)
FILE *fp;
int endoffile,list_no;
{
    int j;

    fprintf(fp,"list: %d\n",list_no);
    if (reverse)
        for (j = no_new_pixels-1; j >= 0; j--)
            fprintf(fp,"%.0f %.0f\n",xo[j],yo[j]);
    else
        for (j = 0; j < no_new_pixels; j++)
            fprintf(fp,"%.0f %.0f\n",xo[j],yo[j]);
    if (endoffile)
        fprintf(fp,"-1 -1\n");
    else
        fprintf(fp,"-1 0\n");
}

store_float_data(fp,endoffile,list_no)
FILE *fp;
int endoffile,list_no;
{
    int j;

    fprintf(fp,"list: %d\n",list_no);
    if (reverse)
        for (j = no_new_pixels-1; j >= 0; j--)
            fprintf(fp,"%f %f\n",xo[j],yo[j]);
    else
        for (j = 0; j < no_new_pixels; j++)
            fprintf(fp,"%f %f\n",xo[j],yo[j]);
    if (endoffile)
        fprintf(fp,"-1 -1\n");
    else
        fprintf(fp,"-1 0\n");
}

store_curvature_data(fp,endoffile,list_no)
FILE *fp;
int endoffile,list_no;
{
    int j;

    fprintf(fp,"list: %d\n",list_no);
    fprintf(fp,"sigma: %f\n",sigma);
    if (reverse)
        for (j = no_new_pixels-1; j >= 0; j--)
            fprintf(fp,"%f %f %f %d\n",xo[j],yo[j],k[j],l[j]);
    else
        for (j = 0; j < no_new_pixels; j++)
            fprintf(fp,"%f %f %f %d\n",xo[j],yo[j],k[j],l[j]);
    if (endoffile)
        fprintf(fp,"-1 -1\n");
    else
        fprintf(fp,"-1 0\n");
}

double random_gauss(mean,var)
double mean,var;
{
    int i;
    double sum = 0;
    double ran3();

    for (i = 0; i < 12; i++)
        sum += (ran3() - 0.5);
    sum *= var;
    sum += mean;
    return(sum);
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

print_usage(pname)
char *pname;
{
    printf("usage: %s\n",pname);
    printf("     -i file    name of input file\n");
    printf("     -o file    name of output file\n");
    printf("     -Y float   Y scaling factor\n");
    printf("     -X float   X scaling factor\n");
    printf("     -Z float   scale data so max X or Y = float\n");
    printf("     -z         use with -Z, but only apply if necessary\n");
    printf("     -x int     X translation\n");
    printf("     -y int     Y translation\n");
    printf("     -c         translate so that centroid is at (300,300)\n");
    printf("     -t         translate so that min X & Y = 0\n");
    printf("     -r float   rotation (degrees)\n");
    printf("     -l file    list of pixel list IDs to select\n");
    printf("     -a int     average over window\n");
    printf("     -s int     subsample\n");
    printf("     -g float   add Gaussian random noise - specify variance\n");
    printf("     -D         do NOT centre each list on rotation\n");
    printf("     -I         convert from pixel_float to pixel (real->integer)\n");
    printf("     -F         convert from pixel to pixel_float (integer->real-)\n");
    printf("     -1         convert to 1 pixel separation\n");
    printf("     -4         convert to 4-way connectivity\n");
    printf("     -8         convert to 8-way connectivity\n");
    printf("     -R         reverse list\n");
    printf("     -S         convert pixel data to super data\n");
    printf("     -C         output in Christine's format\n");
    printf("     -q         output in qhull's format\n");
    printf("     -d         output double (i.e. duplicated pairs) pixel format\n");
    printf("     -G float   interpolate to fill gaps bigger than int\n");
    exit(-1);
}
