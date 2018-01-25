/******************************************
 *
 * pgmio.h:
 *
 * read/write binary portable gray map    
 * images.
 *
 *******************************************/



#ifndef PGMIO_FRANCOIS_H
#define PGMIO_FRANCOIS_H


typedef struct {
    int height;
    int width;
    unsigned char *storage;
} GreyImage;




#define PIXVAL(img,i,j)  ((img).storage[(i)*(img).width + (j)])

/* 
 * all these function will return 1 
 * upon successful completion and 0
 * otherwise.
 */

int createGreyImage(GreyImage *gi, int width, int height);
int readGreyImage(GreyImage *gi, const char *filename);
int writeGreyImage(GreyImage *gi, const char *filename);
int destroyGreyImage(GreyImage *gi);



#endif /* ifndef PGMIO_FRANCOIS_H */
