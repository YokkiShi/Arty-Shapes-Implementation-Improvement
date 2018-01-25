/*
 * convert superdata line to pixel format

 not finished - partially done lines and superellipse
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define INTEGER 0
#define FLOAT 1

#define CLOCKWISE 1
#define ANTICLOCKWISE 2
 
#define PI 3.141591

#define MIN(a,b)  ((a)<(b)?(a):(b))

/* types of data in files */
#define LISTS 1
#define LINES 2
#define ARCS 3
#define ENDL 4
#define ENDF 5
#define ELLIPSES 6
#define SBREAK 7
#define CORNERS 8
#define PIXELS 9
#define CIRCLES 10
#define CODONS 1
#define CROSSES 12
#define VERTEX 13
#define SUPERELLIPSES 14
#define DTRIANGLE 15
#define UTRIANGLE 16
#define BOX 17
#define DIAMOND 18

/* types of files */
#define SUPER 1
#define LIST 4
#define CORNER 5

#ifndef FALSE
# define FALSE 0
# define TRUE (!FALSE)
#endif

#define ONLY2 FALSE

char *infile,*outfile;

float first_root(),angle();

int xpos,ypos;

FILE *fp_in,*fp_out;

int closed = FALSE;
int pairs = FALSE;
int interpolate = FALSE;

main(argc,argv)
int argc;
char *argv[];
{
    int i;

    for (i = 1; i < argc; i++) {
       if (argv[i][0] == '-')
          switch (argv[i][1]) {
              case 'p':
                   pairs = TRUE;
                   break;
              case 'I':
                   interpolate = TRUE;
                   break;
              case 'c':
                   closed = TRUE;
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

    if ((infile == NULL) || (outfile == NULL)) {
       print_usage(argv[0]);
    }

#if ONLY2
/* ??????
    printf("pixel\nlist: 0\n");
*/
#endif

    translate(infile,outfile);

#if ONLY2
/* ??????
    printf("-1 -1\n");
*/
#endif
}

translate(filename1,filename2)
char *filename1,*filename2;
{
    char data_type[40];
    char ch;
    int i,j;
    float xx_s,yy_s,xx_f,yy_f;
    float as_x,as_y,af_x,af_y;
    int old_list_no,list_no;
    float signif;
    char file_type[50];
    int file_data_type;
    int count = 0;
    float x_cent,y_cent;
    int int_radius;
    float temp,ang1,ang2,step,arc_length,l1;
    float x,y,x1,y1,x2,y2,x3,y3;
    int x1d,y1d,x2d,y2d,x3d,y3d;
    float maj_axis,min_axis,squareness,rot_angle;
    int x_s,y_s,x_f,y_f,arc_dir;
    int prevx,prevy;

    prevx = prevy = -9999;

    if ((fp_in = fopen(filename1,"r")) == NULL) {
        printf("cannot open file: %s\n",filename1);
        exit(-1);
    }
    if ((fp_out = fopen(filename2,"w")) == NULL) {
        printf("cannot open file: %s\n",filename2);
        exit(-1);
    }

    fscanf(fp_in,"%s\n",file_type);
    if (strcmp(file_type,"super") == 0) {
        file_data_type = INTEGER;
        fprintf(fp_out,"pixel\n");
    }
    else if (strcmp(file_type,"super_float") == 0) {
        file_data_type = FLOAT;
        fprintf(fp_out,"pixel_float\n");
    }
    else {
        printf("input file not super or super_float type - aborting\n");
        exit(-1);
    }

    list_no = -1;
    do {
        fscanf(fp_in,"%s",data_type);

        j = curve_type(data_type);
        if (j == LISTS) {
            old_list_no = list_no;
            if (old_list_no != -1)
                fprintf(fp_out,"-1 0\n");

            fscanf(fp_in,"%d\n",&list_no);
            fprintf(fp_out,"list: %d\n",list_no);
        }
        else if (j == ENDL) {
            fscanf(fp_in,"\n");
        }
        else if (j == LINES) {
            if (file_data_type == INTEGER) {
                fscanf(fp_in,"%f %d %d %d %d\n",&signif,&x_s,&y_s,&x_f,&y_f);
                x1d = x_s; y1d = y_s; x2d = x_f; y2d = y_f;
            }
            else{
                fscanf(fp_in,"%f %f %f %f %f\n",&signif,&x1,&y1,&x,&y);
                x1d = x1; y1d = y1; x2d = x2; y2d = y2;
            }
            if (interpolate) {
                move(x1d,y1d);
                draw(x2d,y2d);
            }
            if (pairs) {
                fprintf(fp_out,"%d %d\n",x1d,y1d);
                fprintf(fp_out,"%d %d\n",x2d,y2d);
            }
            else {
                if ((prevx != x1d) || (prevy != y1d))
                    fprintf(fp_out,"%d %d\n",x1d,y1d);
                fprintf(fp_out,"%d %d\n",x2d,y2d);
                prevx = x2d;
                prevy = y2d;
            }
        }
        else if (j == ELLIPSES) {
            if (file_data_type == INTEGER) {
                fscanf(fp_in,"%f %d %d %d %d %d %d %f %f %f %d\n",&signif,
                    &x1d,&y1d,&x2d,&y2d,&x3d,&y3d,&maj_axis,
                    &min_axis,&rot_angle,&arc_dir);
                x_cent = x1d;
                y_cent = y1d;
                as_x = x2d;
                as_y = y2d;
                af_x = x3d;
                af_y = y3d;
            }
            else{
                fscanf(fp_in,"%f %f %f %f %f %f %f %f %f %f %d\n",&signif,
                    &x_cent,&y_cent,&as_x,&as_y,&af_x,&af_y,&maj_axis,
                    &min_axis,&rot_angle,&arc_dir);
            }
            plot_ellipse(as_x,as_y,af_x,af_y,x_cent,y_cent,
                             maj_axis,min_axis,rot_angle,arc_dir,
                             1,1);
        }
        else if (j == CIRCLES) {
            if (file_data_type == INTEGER) {
                fscanf(fp_in,"%f %d %d %d\n",&signif,
                    &x1d,&y1d,&int_radius);
					x_cent = x1d;
					y_cent = y1d;
					maj_axis = int_radius;
            }
            else{
                fscanf(fp_in,"%f %f %f %f\n",&signif,
                    &x_cent,&y_cent,&maj_axis);
            }
			plot_circle(x_cent,y_cent,maj_axis,0);
        }
        else if (j == SUPERELLIPSES) {
            if (file_data_type == INTEGER) {
                fscanf(fp_in,"%f %d %d %d %d %d %d %f %f %f %f %d\n",&signif,
                    &x1d,&y1d,&x2d,&y2d,&x3d,&y3d,&maj_axis,
                    &min_axis,&squareness,&rot_angle,&arc_dir);
                    x_cent=x1d; y_cent=y1d;
                    as_x=x2d; as_y=y2d; af_x=x3d; af_y=y3d;
            }
            else{
                fscanf(fp_in,"%f %f %f %f %f %f %f %f %f %f %f %d\n",&signif,
                    &x_cent,&y_cent,&as_x,&as_y,&af_x,&af_y,&maj_axis,
                    &min_axis,&squareness,&rot_angle,&arc_dir);
            }
            plot_superellipse(as_x,as_y,af_x,af_y,x_cent,y_cent,
                              maj_axis,min_axis,squareness,rot_angle,arc_dir);
        }
    } while (j != ENDF);
    fclose(fp_in);

    fprintf(fp_out,"-1 -1\n");
    fclose(fp_out);
}

int curve_type(array)
char array[];
{
    int i,j;

    if ((j = strcmp(array,"list:")) == 0)
        i = LISTS;
    else if ((j = strcmp(array,"line:")) == 0)
        i = LINES;
    else if ((j = strcmp(array,"arc:")) == 0)
        i = ARCS;
    else if ((j = strcmp(array,"endl:")) == 0)
        i = ENDL;
    else if ((j = strcmp(array,"endf:")) == 0)
        i = ENDF;
    else if ((j = strcmp(array,"ellipse:")) == 0)
        i = ELLIPSES;
    else if ((j = strcmp(array,"sbreak:")) == 0)
        i = SBREAK;
    else if ((j = strcmp(array,"corner:")) == 0)
        i = CORNERS;
    else if ((j = strcmp(array,"pixel:")) == 0)
        i = PIXELS;
    else if ((j = strcmp(array,"circle:")) == 0)
        i = CIRCLES;
    else if ((j = strcmp(array,"codon:")) == 0)
        i = CODONS;
    else if ((j = strcmp(array,"cross:")) == 0)
        i = CROSSES;
    else if ((j = strcmp(array,"vertex:")) == 0)
        i = VERTEX;
    else if ((j = strcmp(array,"superellipse:")) == 0)
        i = SUPERELLIPSES;
    else if ((j = strcmp(array,"dtriangle:")) == 0)
        i = DTRIANGLE;
    else if ((j = strcmp(array,"utriangle:")) == 0)
        i = UTRIANGLE;
    else if ((j = strcmp(array,"box:")) == 0)
        i = BOX;
    else if ((j = strcmp(array,"diamond:")) == 0)
        i = DIAMOND;
    return(i);
}

int strcmp(s,t)
char s[],t[];
{
    int i;

    i = 0;
    while (s[i] == t[i])
        if (s[i++] == '\0')
            return(0);
    return(s[i] - t[i]);
}

print_usage(pname)
char *pname;
{
    printf("usage: %s\n",pname);
    printf("     -i file    input superdata file\n");
    printf("     -o file    output pixel file\n");
    printf("     -I         interpolate dense pixel list\n");
    printf("     -c         close list\n");
    printf("     -p         output data as pairs of points per list\n");
    exit(-1);
}

point(xs,ys)
float xs,ys;
{
    unsigned char col1,col2;
    float x,y;

    x = xs;
    y = ys;
    move((int)x,(int)y);
    draw((int)x,(int)y);
}

plot_circle(xc,yc,r,style)
float xc,yc,r;
int style;
{
    int i;
    float tx,ty;
    float angle;
    float factor;
    float x,y;

    factor = r * 0.5;
    /* do a few more points - especially for small circles */
    factor += 100 / r;
    tx = r * cos(0.0);
    ty = r * sin(0.0);
    x = xc + tx;
    y = yc + ty;
    move((int)x,(int)y);
    for (i=0; i < (factor+1); i++) {
        angle = i * 2 * PI / factor;
        tx = r * cos(angle);
        ty = r * sin(angle);
        x = xc + tx;
        y = yc + ty;
        if (style == 0)
            draw((int)x,(int)y);
        else if (style == 1)
            draw((int)x,(int)y);
        else if (style == 2)
            draw((int)x,(int)y);
    }

    /* reset position to original co-oridinates */
    x = xc;
    y = yc;
    move((int)x,(int)y);
}



/*
   plots ellipse using all parameters by generating ellipse around
   origin and transforming it to image space
   parameters: xs,ys to xf,yf - start to finish coords
                xc,yc - centre of ellipse
                major_axis,minor_axis - of the ellipse
                rot_angle - rotation of major axis clockwise about
                origin starting at +ve x axis
   if full is set then plots remainder of ellipse dotted
*/
plot_ellipse(xs,ys,xf,yf,xc,yc,major_axis,minor_axis,rot_angle,arc_dir,bold,full)
float xs,ys,xf,yf,xc,yc;
float major_axis,minor_axis,rot_angle;
int arc_dir,bold,full;
{
    double x1,y1,x2,y2;
    double sine,cosine,sine2,cosine2;
    double temp;
    float step,temp_angle;
    float xold,yold,xe,ye;
    float st_ang,fi_ang;

    if (full) {
        st_ang = 0;
        fi_ang = PI * 2.0;
    }
    else {
        /* determine start and finish angles wrt to centre at the origin */
        /* first transform start and finish to origin offset */
        x1 = xs - xc;
        y1 = ys - yc;
        x2 = xf - xc;
        y2 = yf - yc;
        /* rotate */
        sine = sin(-rot_angle);
        cosine = cos(-rot_angle);
        temp = x1 * cosine - y1 * sine;
        y1 = x1 * sine + y1 * cosine;
        x1 = temp;
        temp = x2 * cosine - y2 * sine;
        y2 = x2 * sine + y2 * cosine;
        x2 = temp;
        st_ang = angle(0.0,0.0,x1/major_axis,y1/minor_axis);
        fi_ang = angle(0.0,0.0,x2/major_axis,y2/minor_axis);
        if(arc_dir == CLOCKWISE) {
            temp = st_ang;
            st_ang = fi_ang;
            fi_ang = temp;
        }
        /* fi_ang must be bigger than st_ang */
        if (fi_ang <= st_ang) {
            fi_ang += PI * 2.0;
        }
    }

    /* generate transform from origin centred ellipse to image space */
    cosine2 = cos(rot_angle);
    sine2 = sin(rot_angle);

    /* moveto first point on ellipse in image */
    xe = major_axis * cos(st_ang);
    ye = minor_axis * sin(st_ang);
    temp = xe * cosine2 - ye * sine2;
    ye = xe * sine2 + ye * cosine2;
    xe = temp;
    xe = xe + xc;
    ye = ye + yc;
    move((int)xe,(int)ye);

    /* successively plot lines around ellipse */
    if (major_axis < 100.0) {
        step = (fi_ang - st_ang) / major_axis;
        step = MIN(step,0.1);
    }
    else
        step = (fi_ang - st_ang) / 100.0;
    temp_angle = st_ang;
    for (temp_angle=st_ang;temp_angle<=fi_ang;temp_angle += step) {
        xe = major_axis * cos(temp_angle);
        ye = minor_axis * sin(temp_angle);
        temp = xe * cosine2 - ye * sine2;
        ye = xe * sine2 + ye * cosine2;
        xe = temp;
        xe = xe + xc;
        ye = ye + yc;
        draw((int)xe,(int)ye);
    }
    xe = major_axis * cos(fi_ang);
    ye = minor_axis * sin(fi_ang);
    temp = xe * cosine2 - ye * sine2;
    ye = xe * sine2 + ye * cosine2;
    xe = temp;
    xe = xe + xc;
    ye = ye + yc;
    draw((int)xe,(int)ye);
}


/*
   plots superellipse using all parameters by generating superellipse
   around origin and transforming it to image space
   parameters: xs,ys to xf,yf - start to finish coords
               xc,yc - centre
               a,b - axes lengths
               e - squareness
               rot_angle - rotation of major axis clockwise about
                           origin starting at +ve X axis
   if full is set then plots remainder of superellipse dotted
*/
plot_superellipse(xs,ys,xf,yf,xc,yc,a,b,e,rot,arc_dir)
float xs,ys,xf,yf,xc,yc;
float a,b,e,rot;
int arc_dir;
{
    float sine,cosine;
    float sine2,cosine2;
    float temp,step,temp_angle;
    float xe,ye,xst,yst;
    float tx,ty;
    double x1,y1,x2,y2;
    float st_ang,fi_ang;
    float t;

    /* determine start and finish angles wrt to centre at the origin */
    /* first transform start and finish to origin offset */
    x1 = xs - xc;
    y1 = ys - yc;
    x2 = xf - xc;
    y2 = yf - yc;
    /* rotate */
    sine = sin(-rot);
    cosine = cos(-rot);
    temp = x1 * cosine - y1 * sine;
    y1 = x1 * sine + y1 * cosine;
    x1 = temp;
    temp = x2 * cosine - y2 * sine;
    y2 = x2 * sine + y2 * cosine;
    x2 = temp;
    st_ang = angle(0.0,0.0,first_root(x1/a,1/e),first_root(y1/b,1/e));
    fi_ang = angle(0.0,0.0,first_root(x2/a,1/e),first_root(y2/b,1/e));
    if(arc_dir == CLOCKWISE) {
        temp = st_ang;
        st_ang = fi_ang;
        fi_ang = temp;
    }
    /* fi_ang must be bigger than st_ang */
    if (fi_ang <= st_ang) {
        fi_ang += PI * 2.0;
    }

    cosine2 = cos(rot);
    sine2 = sin(rot);

    /* move to first point on superellipse in image */
    xe = a * first_root(cos(st_ang),e);
    ye = b * first_root(sin(st_ang),e);
    /* rotate */
    tx = xe * cosine2 - ye * sine2;
    ty = xe * sine2 + ye * cosine2;
    /* shift */
    xe = tx + xc;
    ye = ty + yc;
    move((int)xe,(int)ye);
    if (closed) {
        xst = xe;
        yst = ye;
    }

    /* successively plot lines around ellipse */
    step = 2 * PI / 100;

    for (temp_angle = st_ang; temp_angle <= fi_ang; temp_angle += step) {
        xe = a * first_root(cos(temp_angle),e);
        ye = b * first_root(sin(temp_angle),e);
        /* rotate */
        tx = xe * cosine2 - ye * sine2;
        ty = xe * sine2 + ye * cosine2;
        /* shift */
        xe = tx + xc;
        ye = ty + yc;

        draw((int)xe,(int)ye);
    }

    /* draw to last point on superellipse in image */
    xe = a * first_root(cos(fi_ang),e);
    ye = b * first_root(sin(fi_ang),e);
    /* rotate */
    tx = xe * cosine2 - ye * sine2;
    ty = xe * sine2 + ye * cosine2;
    /* shift */
    xe = tx + xc;
    ye = ty + yc;
    if (closed)
        draw((int)xst,(int)yst);
}


plot_line(xs,ys,xf,yf)
float xs,ys,xf,yf;
{
    move((int)xs,(int)ys);
    draw((int)xf,(int)yf);
}

float first_root(t,e)
float t,e;
/* get first root of the value t for the power e */
{
    if (t > 0)
        return pow(t,e);
    else if (t < 0)
        return -pow(-t,e);
    else
        return 0.0;
}

float angle(x1, y1, x2, y2)
float x1, y1, x2, y2;
{
    float angle_temp;
    float xx, yy;
 
    xx = x2 - x1;
    yy = y2 - y1;
    if (xx == 0.0)
        angle_temp = PI / 2.0;
    else
        angle_temp = atan(fabs(yy / xx));
    if ((xx < 0.0) && (yy >= 0.0))
        angle_temp = PI - angle_temp;
    else if ((xx < 0.0) && (yy < 0.0))
        angle_temp = PI + angle_temp;
    else if ((xx >= 0.0) && (yy < 0.0))
        angle_temp = PI * 2.0 - angle_temp;
    return (angle_temp);
}

/* +++++++++++++++++++++++++ IMAGE DRAWING +++++++++++++++++++++++++ */

draw_point(x,y)
int x,y;
{
    fprintf(fp_out,"%d %d\n",x,y);
}

lineby_image(xd,yd)
int xd,yd;
{
    int loop1, loop2, i, xstep, ystep;

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
    for (loop1=1;loop1<=xd;loop1++) {
        for (loop2=1;loop2<=yd;loop2++) {
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

draw(xd,yd)
int xd,yd;
{
    int xnew,ynew;

    xnew = xd;
    ynew = yd;
#if ONLY2
    printf("%d %d\n",xd,yd);
#endif
    lineby_image(xnew-xpos,ynew-ypos);
}

move(xd,yd)
int xd,yd;
{
    xpos = xd;
    ypos = yd;
#if ONLY2
    printf("%d %d\n",xd,yd);
#endif
}

