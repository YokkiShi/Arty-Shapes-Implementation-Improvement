/*
    Basterdised version of superellipses.c
    fits one superellipse to complete pixel list without recursing or segmenting

    Input: pixel data format.
    Output: super data format (only arcs).

    -------

    modified: Paul Rosin, July 2005
    uses ellipse_specific fit for initialisation

	computes an ovalness measure
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TESTING 0         /* output in non-super format */
#define DO_ARC  0         /* fit arc, otherwise fit complete superellipse */

#define NO_LINE_SEGS 10000 /* max number of lines per list */
#define NO_ARCS 200       /* max number of arcs per list  */
#define LINE 1
#define ARC 2
#define PI 3.141591
#define CLOCKWISE 1       /* direction of arc from start to finish */
#define ANTICLOCKWISE 2

#define DEBUG 0           /* compile debug output if 1 */
#define STOP_EARLY 1      /* just fit to all data in a list and stop if 1 */
                          /* normal recursive algorithm if = 0 */
#define NO_PIXELS 10000
#define NO_ARCS 200       /* max number of arcs per list  */

#define PI 3.141591
#define NINETY 1.570796327
#define KEEP 0         /* KEEP and REJECT used to tell if arcs to be kept */
#define REJECT 1
#define LARGE_SIG 99999
#define MIN_LENGTH 7
#define FALSE 0
#define TRUE (!FALSE)

#define ELLIPSE 0
#define PARABOLA 1
#define HYPERBOLA 2

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

#define SQR(x)    ((x) * (x))
#define ABS(x)    (((x)<0.0)? (-(x)): (x))

/* temp array for manipulatable data from x_c,y_c */
float x_trans3[NO_LINE_SEGS],y_trans3[NO_LINE_SEGS];

/* ellipse parameterisation arrays */
float arc_centre_x[NO_ARCS],arc_centre_y[NO_ARCS];
int arc_start_x[NO_ARCS],arc_start_y[NO_ARCS];
int arc_end_x[NO_ARCS],arc_end_y[NO_ARCS];
float arc_squareness[NO_ARCS];
float maj_axis[NO_ARCS],min_axis[NO_ARCS];
float rot_ang[NO_ARCS];
int arc_dirs[NO_ARCS];
float arc_sig[NO_ARCS];

int number_arcs;
int nseg2;   /* number of points in x_trans3, y_trans3 */

int arc_status[NO_ARCS];
int x_c[NO_PIXELS],y_c[NO_PIXELS];
int arc_start[NO_ARCS],arc_finish[NO_ARCS];

int number_pixels;

int global_pos;
int representation_ok;
int measure = ALGEBRAIC;

FILE *fp,*fp_out;
float distance();
float first_root();

int bias = FALSE;
double mean_error;

main(argc,argv)
int argc;
char *argv[];
{
    int i,j;
    float sig;    /* not really used - just dummy param for segment */
    int list_no,end_of_file;
    float temp;
    char file_type[50];
    char *file_name_in;
    char *file_name_out;
    char *ch;
    float delete = -1;

    int flag;
    double xc, yc, area;

    file_name_in = file_name_out = NULL;

    /* parse command line */
    for(i = 1; i < argc; i++) {
        if (argv[i][0] == '-') {
            switch(argv[i][1]) {
                case 'i':
                    i++;
                    file_name_in = argv[i];
                    break;
                case 'o':
                    i++;
                    file_name_out = argv[i];
                    break;
                case 'd':
                    i++;
                    delete = atof(argv[i]);
                    break;
                case 'm':
                    i++;
                    measure = atoi(argv[i]);
                    break;
                case 'b':
                    bias = TRUE;
                    break;
                default:
                    printf("unknown option %s\n",argv[i]);
                    options(argv[0]);
            }
        }
        else {
            printf("unknown option %s\n",argv[i]);
            options(argv[0]);
        }
    }

    if (file_name_in == NULL) {
        printf("no input file specified - aborting\n");
        options(argv[0]);
    }
    if (file_name_out == NULL) {
        printf("no output file specified - aborting\n");
        options(argv[0]);
    }

    switch(measure) {
        case ALGEBRAIC: printf("using ALGEBRAIC measure\n"); break;
        case GROSS: printf("using GROSS measure\n"); break;
        case GRADIENT: printf("using GRADIENT measure\n"); break;
        case DIRECTIONAL: printf("using DIRECTIONAL measure\n"); break;
        case RAY: printf("using RAY measure\n"); break;
        case SAFAEE: printf("using SAFAEE measure\n"); break;
        case EXPAND: printf("using EXPAND measure\n"); break;
        case CONSTANT: printf("using CONSTANT measure\n"); break;
        case QUADRARC: printf("using QUADRARC measure\n"); break;
        case RECTANGLE: printf("using RECTANGLE measure\n"); break;

        default:
            fprintf(stderr,"ERROR: %d unknown distance measure\n",measure);
            exit(-1);
    }

    if ((fp = fopen(file_name_in,"r")) == NULL) {
        printf("cant open %s\n",file_name_in);
        exit(-1);
    }
    /* read magic word for format of file */
    fscanf(fp,"%s\n",file_type);
    j = strcmp(file_type,"pixel");
    if (j != 0) {
        printf("not pixel data file - aborting\n");
        exit(-1);
    }

    if ((fp_out = fopen(file_name_out,"w")) == NULL) {
        printf("cant open %s\n",file_name_out);
        exit(-1);
    }
    /* write magic word for format of file */
#if !TESTING
    fprintf(fp_out,"super\n");
#endif

    do{
        read_list(&number_pixels,&end_of_file,&list_no);

        /* add these values for compatability with polyCentroid */
        x_c[0] = x_c[number_pixels];
        y_c[0] = y_c[number_pixels];
        flag = polyCentroid(x_c, y_c, number_pixels, &xc, &yc, &area);

        /* delete some data to simulate occlusion */
        if (delete != -1)
            number_pixels *= delete;

        number_arcs = 0;
        segment(1,number_pixels,&sig);

        printf("ovalness: %f\n",1.0 / (1 + mean_error/sqrt(area)));

        /* super data output format */
        for (j = 1;j <= number_arcs;j++)
           if (arc_status[j] == KEEP) {
              if (arc_sig[j] < (LARGE_SIG - 1)) {
                 /* !! HACK !! */
                 if (maj_axis[j] < min_axis[j]) {
                    temp = maj_axis[j];
                    maj_axis[j] = min_axis[j];
                    min_axis[j] = temp;
                    rot_ang[j] += PI/2;
                 }
#if TESTING
                 fprintf(fp_out,"ctr %f %f\n",arc_centre_x[j],arc_centre_y[j]);
                 fprintf(fp_out,"a b %f %f\n",maj_axis[j],min_axis[j]);
                 fprintf(fp_out,"theta %f\n",rot_ang[j]);
                 fprintf(fp_out,"epsilon %f\n",arc_squareness[j]);
#else
                 fprintf(fp_out,"superellipse: %f %.0f %.0f %d %d %d %d %f %f %f %f %d\n",
                    arc_sig[j],
                    arc_centre_x[j],arc_centre_y[j],
                    arc_start_x[j],arc_start_y[j],
                    arc_start_x[j],arc_start_y[j],  // we want a complete SE
                    //arc_end_x[j],arc_end_y[j],
                    maj_axis[j],min_axis[j],arc_squareness[j],rot_ang[j],
                    arc_dirs[j]);
#endif
              }
              else
                 fprintf(fp_out,"line: 0.0 %d %d %d %d\n",
                    arc_start_x[j],arc_start_y[j],
                    arc_end_x[j],arc_end_y[j]);
           }
#if !TESTING
        fprintf(fp_out,"endl:\n");
        if (end_of_file == 1)
            fprintf(fp_out,"endf:\n");
#endif
    }while(end_of_file == 0);
}

read_list(number_pixels,end_of_file,list_no)
int *end_of_file,*number_pixels,*list_no;
{
    char dumstring[50];
    int j;
    int tx,ty;

    fscanf(fp,"%s %d\n",dumstring,list_no);
    j = 0;
    do{
        j++;
        fscanf(fp,"%d %d\n",&tx,&ty);
        x_c[j] = tx;
        y_c[j] = ty;
    }while(x_c[j] != -1);
    if (y_c[j] == -1)
       *end_of_file = 1;
    else
       *end_of_file = 0;
    *number_pixels = --j;
}

segment(start_in,finish_in,sig_out)
int start_in,finish_in;
float *sig_out;
{
    int i;
    int pos;
    float max_dev;
    float sig1,sig2,sig3,max_sig;
    float centre_x, centre_y;
    int arc_length,arc_dir;
    float major_axis,minor_axis,squareness,rot_angle;

    /* compute significance at this level */

    if ((finish_in - start_in) < MIN_LENGTH)
        representation_ok = FALSE;
    else
        determine_superellipse_fit(start_in,finish_in,&major_axis,&minor_axis,
                &squareness,&rot_angle,&centre_x,&centre_y,
                &max_dev,&arc_length,&arc_dir,&sig1);

    if (representation_ok == FALSE)
       sig1 = LARGE_SIG;

    pos = global_pos + start_in - 1;

    number_arcs++;
    arc_sig[number_arcs] = sig1;
    arc_centre_x[number_arcs] = centre_x;
    arc_centre_y[number_arcs] = centre_y;
    arc_start_x[number_arcs] = x_c[start_in];   /* to change? */
    arc_start_y[number_arcs] = y_c[start_in];   /*      "     */
    arc_end_x[number_arcs] = x_c[finish_in];    /*      "     */
    arc_end_y[number_arcs] = y_c[finish_in];    /*      "     */
    arc_dirs[number_arcs] = arc_dir;
    maj_axis[number_arcs] = major_axis;
    min_axis[number_arcs] = minor_axis;
    arc_squareness[number_arcs] = squareness;
    rot_ang[number_arcs] = rot_angle;
    arc_start[number_arcs] = start_in;
    arc_finish[number_arcs] = finish_in;
    arc_status[number_arcs] = KEEP;
    if (((finish_in - start_in) < MIN_LENGTH) ||
         (max_dev < 3) ||
         (STOP_EARLY == 1))
    {
        /* at lowest level of recursion */
        /* save arc match at this lowest of levels */
        /* and significance */

        *sig_out = sig1;
    }
    else{
        /* recurse to next level down */
        segment(start_in,pos,&sig2);
        segment(pos,finish_in,&sig3);
        /* get best significance from lower level */
        if (sig2 < sig3)
            max_sig = sig2;
        else
            max_sig = sig3;
        if (max_sig < sig1) {
            /* return best significance, keep lower level description */
            *sig_out = max_sig;
        /* remove the single arc from start_in to finish_in */
        for (i = 0;i <= number_arcs;i++)
            if ((arc_start[i] == start_in) &&
                (arc_finish[i] == finish_in))
            arc_status[i] = REJECT;
        }
        else{
            /* line at this level is more significant */
            /* remove arcs at lower levels */
            /* i.e between start_in and finish_in */
            /* DOES IT EVER GET HERE ???? */
            printf("JUST REPLACED LOWER LEVEL ARCS\n");
            *sig_out = sig1;
            for (i = 0;i <= number_arcs;i++) {
                if ((arc_start[i] == start_in) &&
                    (arc_finish[i] == finish_in))
                {
                    /* do nothing - leave it */
                }
                else if ((arc_start[i] >= start_in) &&
                         (arc_finish[i] <= finish_in))
                {
                    arc_status[i] = REJECT;
                }
            }
        }
    }
}

/*int round(value)
float value;
{
    int i;

    i = floor(value + 0.5);
    return(i);
}*/

float angle(x1,y1,x2,y2)
float x1,y1,x2,y2;
{

    float angle_temp;
    float xx,yy;

    xx = x2 - x1;
    yy = y2 - y1;
    if (xx == 0.0)
        angle_temp = PI/2.0;
    else
        angle_temp = atan(fabs(yy/xx));
    if ((xx < 0.0) && (yy >= 0.0))
       angle_temp = PI - angle_temp;
    else if ((xx < 0.0) && (yy < 0.0))
       angle_temp = PI + angle_temp;
    else if ((xx >= 0.0) && ( yy < 0.0))
       angle_temp = PI*2.0 - angle_temp;
    return(angle_temp);
}

float euclid(x1,y1,x2,y2)
float x1,y1,x2,y2;
{
    float temp1,temp2;
    float dist;

    temp1 = fabs(x1-x2);
    temp2 = fabs(y1-y2);
    dist = sqrt(SQR(temp1)+SQR(temp2));
    return(dist);
}

/* rotate polygon segment (x1,y1) -> (x2,y2) to lie along x axis,
 * compute the minimum distance between the line segment and
 * a point on the circle defined as (x_cir,y_cir).
 * take perpendicular distance if x_cir between end points else
 * take euclidian distance to nearest end point.
 */
compute_minimum_diff(deviation,x_cir,y_cir,x1,y1,x2,y2)
float x_cir,y_cir;
float x1,y1,x2,y2;
float *deviation;
{
    float angle1;
    float x_off,y_off,cosine,sine;
    float temp;
    float min1,min2;

/*
    printf("entered compute_minimum_diff\n");
    printf("x_cir: %f y_cir %f\n",x_cir,y_cir);
    printf("x1: %f y1 %f x2 %f y2 %f\n",x1,y1,x2,y2);
*/
    angle1 = angle(x1,y1,x2,y2);
    cosine = cos(-angle1);
    sine = sin(-angle1);
    x_off = x1;
    y_off = y1;
    /* offset points so x1,y1 at origin */
    x1 = 0;
    y1 = 0;
    x2 = x2 - x_off;
    y2 = y2 - y_off;
    x_cir = x_cir - x_off;
    y_cir = y_cir - y_off;
    /* rotate points with x2,y2 on x axis */
    temp = x2*cosine - y2*sine;
    y2 = x2*sine + y2*cosine;
    x2 = temp;
    temp = x_cir*cosine - y_cir*sine;
    y_cir = x_cir*sine + y_cir*cosine;
    x_cir = temp;
/*
    printf("x_cir: %f y_cir %f\n",x_cir,y_cir);
    printf("x1: %f y1 %f x2 %f y2 %f\n",x1,y1,x2,y2);
*/
    min1 = euclid(x_cir,y_cir,x1,y1);
    min2 = euclid(x_cir,y_cir,x2,y2);
    if (x_cir < x1) {
        *deviation = min1;
    }else if (x_cir > x2) {
        *deviation = min2;
    }else{
        *deviation = fabs(y_cir);
    }
/*    printf("deviation: %f\n",*deviation,); */
}

/*
 * New version that walks around hypothesised superellipse
 * finding minimum distance to polygon at each step.
 * Nearest point is either the perpendicular distance,
 * or the euclidian distance to the nearest end point.
 * For all cases the nearest end point is taken as the
 * possible break point. max_dev is the deviation,
 * max_pos is the vertex for max_dev,
 * major_axis,minor_axis define the ellipse
 */
compute_dev(max_dev,max_pos,major_axis,minor_axis,squareness,arc_dir)
float major_axis,minor_axis,squareness,*max_dev;
int arc_dir,*max_pos;
{
    int l2;
    float l1,step;
    float s_ang,f_ang,temp;
    float x_cir,y_cir;
    float deviation;
    float min_dev;
    int min_pos;
    float t;

    *max_dev = 0.0;

    /* determine angles */
    s_ang = angle(0.0,0.0,
            first_root(x_trans3[1]/major_axis,squareness),
            first_root(y_trans3[1]/minor_axis,squareness));
    f_ang = angle(0.0,0.0,
                  first_root(x_trans3[nseg2]/major_axis,squareness),
                  first_root(y_trans3[nseg2]/minor_axis,squareness));
#if DEBUG
    printf("start coords: %f %f\n",x_trans3[1],y_trans3[1]);
    printf("finish coords: %f %f\n",x_trans3[nseg2],y_trans3[nseg2]);
    printf("angles: %f %f\n",s_ang,f_ang);
#endif
    /* default is anticlockwise, swap if clockwise */
    if (arc_dir == CLOCKWISE) {
        temp = s_ang;
        s_ang = f_ang;
        f_ang = temp;
    }
/* ang4 must be bigger than ang3 */
    if (f_ang < s_ang) {
        f_ang += PI*2.0;
    }
#if DEBUG
    printf("start angle: %f finish angle %f\n",s_ang,f_ang);
/* walk around circle from angle1 to angle2 */
    printf("walking around circle\n");
#endif
    step = (f_ang-s_ang)/10.0;    /* use 11 steps presently */
    l1 = s_ang;
    do{
        min_dev = 1000000;
        x_cir = major_axis * first_root(cos(l1),squareness);
        y_cir = minor_axis * first_root(sin(l1),squareness);
        l2 = 1;
        do{
            /* printf("for line %d\n",l2);   */
            compute_minimum_diff(&deviation,x_cir,y_cir,x_trans3[l2],
                                 y_trans3[l2],x_trans3[l2+1],y_trans3[l2+1]);
            if (deviation < min_dev) {
                min_dev = deviation;
                min_pos = l2;
            }
            l2++;
        }while(l2 <= nseg2-1);

        if (min_dev > *max_dev) {
            *max_dev = min_dev;
            *max_pos = min_pos;
        }
        l1 += step;
    }while (l1 <= f_ang);
#if DEBUG
    printf("maximum deviation: %f\n",*max_dev);
#endif
}

/* computes length of arc
 * works on temporary data from determine_superellipse_fit
 * i.e. data: (x_trans3,y_trans3),major_axis,minor_axis.
 * new version computes start and finish angles, then
 * arc length from angle and radius
 */
compute_lgt(lgt,major_axis,minor_axis,squareness,arc_dir)
float *lgt,major_axis,minor_axis,squareness;
int arc_dir;
{
    float s_ang,f_ang,temp;
    float xold,yold,xnew,ynew;
    float step;
    float t;

    s_ang = angle(0.0,0.0,
            first_root(x_trans3[1]/major_axis,squareness),
            first_root(y_trans3[1]/minor_axis,squareness));
    f_ang = angle(0.0,0.0,
            first_root(x_trans3[nseg2]/major_axis,squareness),
            first_root(y_trans3[nseg2]/minor_axis,squareness));
    /* default is anticlockwise, so swap angles if clockwise */
    if (arc_dir == CLOCKWISE) {
       temp = s_ang;
       s_ang = f_ang;
       f_ang = temp;
    }
    /* start angle must be less that f_ang */
    if (f_ang < s_ang)
       f_ang += PI*2.0;
    *lgt = 0.0;

    xold = major_axis * first_root(cos(s_ang),squareness);
    yold = minor_axis * first_root(sin(s_ang),squareness);

    step = (f_ang - s_ang) / 20.0;
    temp = s_ang;
    do{
        temp += step;
        xnew = major_axis * first_root(cos(temp),squareness);
        ynew = minor_axis * first_root(sin(temp),squareness);
        *lgt = *lgt + distance(xold,yold,xnew,ynew);
        xold = xnew;
        yold = ynew;
    }while(temp <= f_ang);
}

compute_poly_lgt(lgt)
float *lgt;
{
    int loop1;
    float dx,dy;

    *lgt = 0;
    for (loop1 = 1;loop1 < nseg2;loop1++) {
        dx = x_trans3[loop1] - x_trans3[loop1+1];
        dy = y_trans3[loop1] - y_trans3[loop1+1];
        *lgt += sqrt(SQR(dx) + SQR(dy));
    }
}

/*
 * determine the best fit superellipse
 */
determine_superellipse_fit(st,fi,final_major_axis,final_minor_axis,
                           final_squareness,final_rot_angle,
                           final_xc,final_yc,
                           final_dev,final_lgt,arc_dir,sig)
int st,fi;
float *final_major_axis,*final_minor_axis;
float *final_squareness,*final_rot_angle;
float *final_xc,*final_yc;
int *final_lgt,*arc_dir;
float *final_dev,*sig;
{
    /* new variables for lms fitting */
    float c_ang;

    /* original variables */
    float xt,yt;   /* temp variables */
    float x_cent_t,y_cent_t; /* temp centre coords - after transformation */
    float sine,cosine;
    float max_dev,arc_length,poly_length;
    int loop1;
    int max_pos;
    float sum;
    float ratio = 1;
    float temp;
    float x_off,y_off;

    float a,b,e,xc,yc,rot;

    float angle_s,angle_f,diff_angle;

    *sig = 0;

    /*
       set representation_ok to TRUE - set to FALSE by any code that
       cannot fit the representation to the data
    */
    representation_ok = TRUE;

    /* get data into temp array */
    nseg2 = 0;
    for (loop1 = st; loop1 <= fi; loop1++) {
       nseg2++;
       x_trans3[nseg2] = x_c[loop1];
       y_trans3[nseg2] = y_c[loop1];
    }

    fit_superellipse(&a,&b,&e,&xc,&yc,&rot);

    if (representation_ok == FALSE)
       return;

#if DO_ARC
   /* determine arc_size
    * transform temp to lie along X axis, then compute average displacement
    * of y coords. If +ve then circle on +ve side of axis else on -ve. This
    * in combination with the position of the centre will indicate if the are
    * is large (>180 degrees) or small (<=180 degrees). Also can determine if
    * arc is clockwise or anticlockwise
    * rotate and translate temp data.
    * take chord between points, rotate so that chord is
    * along x axis - rotate and translate all other points the same.
    */
    c_ang = angle(x_trans3[1],y_trans3[1],x_trans3[nseg2],y_trans3[nseg2]);
    sine = sin(-c_ang);
    cosine = cos(-c_ang);

    /* translate so first point is at origin */
    x_off = x_trans3[1];
    y_off = y_trans3[1];
    for (loop1 = 1;loop1 <= nseg2;loop1++) {
        x_trans3[loop1] = x_trans3[loop1] - x_off;
        y_trans3[loop1] = y_trans3[loop1] - y_off;
    }
    /* rotate to align chord with x axis */
    for (loop1 = 1;loop1 <= nseg2;loop1++) {
        xt = x_trans3[loop1];
        yt = y_trans3[loop1];
        temp = xt * cosine - yt * sine;
        x_trans3[loop1] = temp;
        temp = xt * sine + yt * cosine;
        y_trans3[loop1] = temp;
    }

    /* do same for centre */
    xt = xc - x_off;
    yt = yc - y_off;
    temp = xt*cosine - yt*sine;
    x_cent_t = temp;
    temp = xt*sine + yt*cosine;
    y_cent_t = temp;

    /* compute average Y coord of all points */
    sum = 0.0;
    for (loop1=1;loop1<=nseg2;loop1++)
        sum += y_trans3[loop1];

    /* determine size and sense of arc */
    if ((sum >= 0.0) && (y_cent_t >= 0.0)) {
        *arc_dir = CLOCKWISE;
    }
    else if ((sum < 0.0) && (y_cent_t < 0.0)) {
        *arc_dir = ANTICLOCKWISE;
    }
    else if ((sum >= 0.0) && (y_cent_t < 0.0)) {
        *arc_dir = CLOCKWISE;
    }
    /* CHANGED Y_CENT TO Y_CENT_T - PLR */
    else if ((sum < 0.0) && (y_cent_t >= 0.0)) {
        *arc_dir = ANTICLOCKWISE;
    }
    /* !!! FOR VERY FLAT SECTIONS !!! */
    if (fabs(sum/nseg2) < 3.0) {
       angle_s = angle(x_cent_t,y_cent_t,x_trans3[1],y_trans3[1]);
       angle_f = angle(x_cent_t,y_cent_t,x_trans3[nseg2],y_trans3[nseg2]);
       diff_angle = angle_f - angle_s;
       if (diff_angle < 0) diff_angle += (PI*2.0);
       if (diff_angle < PI)
          *arc_dir = ANTICLOCKWISE;
       else
          *arc_dir = CLOCKWISE;
    }

    /*
    for (loop2=1;loop2<=nseg2;loop2++)
        printf("x_trans3: %f y_trans3: %f\n",x_trans3[loop2],y_trans3[loop2]);
    */

    /* transform data so major_axis axis of superellipse is along X axis */
    /* get data into temp arrays */
    nseg2 = 0;
    for (loop1 = st; loop1 <= fi; loop1++) {
        nseg2++;
        x_trans3[nseg2] = x_c[loop1];
        y_trans3[nseg2] = y_c[loop1];
    }
    /* first translate */
    for (loop1=1; loop1<=nseg2; loop1++) {
        x_trans3[loop1] = x_trans3[loop1] - xc;
        y_trans3[loop1] = y_trans3[loop1] - yc;
    }
    /* now rotate */
    sine = sin(-rot);
    cosine = cos(-rot);
    for (loop1=1; loop1<=nseg2; loop1++) {
        xt = x_trans3[loop1];
        yt = y_trans3[loop1];
        temp = xt*cosine - yt*sine;
        x_trans3[loop1] = temp;
        temp = xt*sine + yt*cosine;
        y_trans3[loop1] = temp;
    }
    x_cent_t = 0.0;
    y_cent_t = 0.0;

    /*
     * determine length of ellipse for this circle approx - use data
     * transformed to put major axis along x axis
     */
    compute_lgt(&arc_length,a,b,e,*arc_dir);

    /*
     * determine actual length of data - a polygon here
     */
    compute_poly_lgt(&poly_length);

    /*
     * determine position and value of the maximum deviation - use data
     * transformed to put major axis along x axis
     */
    compute_dev(&max_dev,&max_pos,a,b,e,*arc_dir);

    /*
     * hack to stop breakpoint near ends - to be modified
     * break at middle of data if too close to end points
     */
    if ((max_pos <= MIN_LENGTH) || (max_pos > nseg2-MIN_LENGTH))
        max_pos = nseg2 / 2;

#if DEBUG
    printf("arc length: %f\n",arc_length);
    printf("poly length: %f\n",poly_length);
    printf("maximum deviation: %f\n",max_dev);
    printf("point of maximum deviation: %d\n",max_pos);
#endif

    /*
     * compute significance based on max_dev, arc_length and
     * some heuristics!
     */
    /* polygonal length should be similar to arc length */
    ratio = poly_length / arc_length;
    if (ratio < 1)
        ratio = arc_length / poly_length;

    *final_lgt = arc_length;
    max_dev = max_dev * ratio;        /* temp to modify sig */
    *final_dev = max_dev;

    if ((max_dev != 0.0) && (arc_length != 0))
        *sig = max_dev / (float)arc_length;
    else
        *sig = LARGE_SIG;

    /* save position of maximum deviation */
    global_pos = max_pos;
#else
    *arc_dir = 1;
#endif

    /* save parameters for best ellipse as integers */
    *final_major_axis = a;
    *final_minor_axis = b;
    *final_squareness = e;
    *final_rot_angle = rot;
    *final_xc = xc;
    *final_yc = yc;
}

float distance(x1,y1,x2,y2)
float x1,y1,x2,y2;
{
    return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}

int curve_type(array)
char array[50];
{
    int i,j;

    if ((j = strcmp(array,"list:")) == 0)
        i = 1;
    else if ((j = strcmp(array,"line:")) == 0)
        i = 2;
    else if ((j = strcmp(array,"arc:")) == 0)
        i = 3;
    else if ((j = strcmp(array,"endl:")) == 0)
        i = 4;
    else if ((j = strcmp(array,"endf:")) == 0)
        i = 5;
    return(i);
}

int strcmp(s,t)
char s[],t[];
{
    int i;

    i = 0;
    while(s[i] == t[i])


        if (s[i++] == '\0')
            return(0);
    return(s[i] - t[i]);
}

/*
 *  code below used to be in determ3.c
 *  fits ellipses to data using a kalman filter
 *  abort when any error condition occurs - that is cannot fit an ellipse
 *  Paul Rosin
 *  January 1990
 */

int find_conic_type(c1,c2,c3)
float c1,c2,c3;
{
    float fx;
    int temp;

    fx = c1 * c3 - c2 * c2 / 4;
        if (fx > 0)
           temp = ELLIPSE;
    else if (fx < 0)
           temp = HYPERBOLA;
    else
           temp = PARABOLA;
/*      printf("fx: %f\n",fx); */
    return(temp);
}

error(s)
char *s;
{
    printf("ERROR: %s\n",s);
    exit(-1);
}

float first_root(t,e)
float t,e;
/* get first root of the value t for the power e */
{
    if (t > 0)
        return pow((double)t,(double)e);
    else if (t < 0)
        return -pow(-(double)t,(double)e);
    else
        return 0.0;
}

/*
 * ANSI C code from the article
 * "Centroid of a Polygon"
 * by Gerard Bashein and Paul R. Detmer,
 *  (gb@locke.hs.washington.edu, pdetmer@u.washington.edu)
 * in "Graphics Gems IV", Academic Press, 1994
 */

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
int polyCentroid(int x[], int y[], int n, double *xCentroid, double *yCentroid, double *area)
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

    if (atmp != 0) {
        *xCentroid = xtmp / (3 * atmp);
        *yCentroid = ytmp / (3 * atmp);
        return 0;
    }
    return 2;
}

options(progname)
char *progname;
{
    printf("usage: %s [options]\n",progname);
    printf("     -i file    input pixel list\n");
    printf("     -o file    output superellipse\n");
    printf("     -d float   delete - only keep proportion\n");
    printf("     -m int     distance measure\n");
    printf("     -b         bias fit to avoid stars\n");
    exit(-1);
}
