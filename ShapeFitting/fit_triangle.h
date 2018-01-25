/*************************************************************************/
/* defs.h                                                                */
/* Funktion: Definition der Fehlermeldungen Konstanten                   */
/*************************************************************************/

/* Kontrolle von Speicheranforderungen und Freigaben */

/************************************************
 *                                              *
 *       Fehlernachrichten                      *
 *                                              *
 ************************************************/

#define M_0 			""
#define M_INTERN		"Softwarefehler in ICE"
#define M_OS			"Betriebssystemfehler "
#define M_NO_MEM		"Nicht genug Speicher vorhanden"
#define M_SYSTEM_1		"Systemfehler: Code 1"
#define M_SYSTEM_2		"Systemfehler: Code 2"
#define M_SYSTEM_3		"Systemfehler: Code 3"
#define M_SYSTEM_4		"Systemfehler: Code 4"
#define M_SYSTEM_5		"Systemfehler: Code 5"
#define M_SYSTEM_6		"Systemfehler: Code 6"
#define M_SYSTEM_7		"Systemfehler: Code 7"
#define M_SYSTEM_8		"Systemfehler: Code 8"
#define M_SYSTEM_9		"Systemfehler: Code 9"
#define M_SYSTEM_10		"Systemfehler: Code 10"

/* Allgemeine Parameterfehler */
#define M_NOT_FOUND		"Nicht gefunden"
#define M_WRONG_PARAM		"Falsche Parameterbelegung"
#define M_WRONG_MODE		"Unzulaessiger Modus"
#define M_WRONG_LEN		"Unzulaessige Laengenangabe"
#define M_WRONG_VAL		"Unzulaessiger Wert-Parameter"
#define M_WRONG_MAGNITUDE	"Unzulaessige Groessenangabe"
#define M_WRONG_KOORD		"Unzulaessige Koordinatenangaben"
#define M_WRONG_PARA2		"Negative Radiusangabe"
#define M_INV_STRUCT   		"Fehlende/falsche Parameter in Datenstruktur"
#define M_WRONG_POINTLIST 	"Unzulaessige Punktliste"
#define M_TO_MUCH_POINTS        "Zu wenig Punkte in Punktliste"
#define M_POINTLIST_ADR		"Unzulaessige Listenadressierung"
#define M_WRONG_PTR		"Unzulaessiger Pointer"
#define M_XTOOSMALL             "x-Koordinate zu klein"
#define M_XTOOLARGE             "x-Koordinate zu groﬂ"
#define M_YTOOSMALL             "y-Koordinate zu klein"
#define M_YTOOLARGE             "y-Koordinate zu groﬂ"
#define M_VALTOOSMALL           "Pixelwert zu klein"
#define M_VALTOOLARGE           "Pixelwert zu groﬂ"

/* String */
#define M_WRONG_STRPTR 		"Unzulaessiger Stringpointer"
#define M_STR_TO_LONG		"String zu lang"

/* Dateiarbeit */
#define M_WRONG_FILE		"Dateifehler"
#define M_FILE_OPEN		"Fehler beim File Oeffnen"
#define M_WRONG_WRITE		"Fehler beim File Schreiben"

/* Fehler bei Visualisierung */
#define M_NOT_VIS		"Bild nicht zur Visualisierung angemeldet"

/* Arbeit mit Bildern */
#define M_WRONG_IMGPTR 		"Unzulaessiger Bildpointer"
#define M_WRONG_IMGDPTR		"Unzulaessiger Double-Bildpointer"
#define M_WRONG_IMGSIZE 	"Falsche Bildgroesse"

/* Arbeit mit Fenstern */
#define M_WRONG_WINDOW 		"Passt nicht in aktuellen Fensterbereich"
#define M_WRONG_WINDOW2		"Gemeinsamer Fensterbereich zu klein"
#define M_WRONG_WNDW3		"Unzulaessiger Fensterbereich"

/* Konturen*/
#define M_WRONG_STARTPOINT 	"Unzulaessiger Startpunkt"
#define M_WRONG_STARTPOINT2 	"Isolierter Startpunkt"
#define M_WRONG_STARTPOINT3 	"Startpunkt liegt im Objekt"
#define M_WRONG_CPTR		"Unzulaessiger Konturpointer"
#define M_NO_SEGMENT 		"Kein Segment in Kontur gefunden"
#define M_CLINE			"Systemfehler bei Aufruf von \"AddLineContur\""
#define M_CONOCLOSED		"Kontur nicht geschlossen !"
#define M_NO_ELLIPLSE

/* Numerik */
#define M_DIVISION_ZERO  	"Division durch Null"
#define M_NO_SOLUTION		"Keine Loesung vorhanden"
#define M_SOL_MANIFOLD          "Loesungsvielfalt"
#define M_NUM_INSTABILITY	"Numerische Instabilitaet"
#define M_ZERO_VECTOR           "Nullvektor"
#define M_NO_ORTHO		"Keine Orthonormalmatrix"
#define M_WRONG_VECTOR          "Unzulaessiger Vektorpointer"
#define M_EQUSYS_ERROR          "Fehler in Prozedur Equsys"
#define M_ZERO_DET	  	"Nicht loesbar (Determinante 0)"
#define M_NO_QUAD		"Nichtquadratischen Matrix"
#define M_WRONG_MATRIX          "Unzulaessige Matrix"
#define M_MAT_NO_COMPAT         "Matrizen nicht kompatibel"
#define M_MATRIX_SINGULAR       "Singularmatrix"
#define M_WRONG_MATRIXTYPE      "Falscher Matrix-Typ"
#define M_NO_REGULAR            "Matrix nicht regulaer"
#define M_NO_SYMM               "Matrix nicht symmetrisch"
#define M_VEC_DEPEND		"Vektoren nicht unabhaengig"
#define M_NO_INVERS		"Inverse existiert nicht"
#define M_WRONG_START		"Ungeeigneter Startwert"
#define M_POINT_IDENTIC		"Punkte identisch"
#define M_WRONG_POINTS		"Unzulaessige Punktauswahl"
#define M_WRONG_TRANS		"Ungzulaessige Transformationsmatrix"
#define M_NO_PROJ               "Keine projektive Matrix (dritte Zeile 0)"


/* Histogramme */
#define M_DESTR_HISTCHAIN	"verkettete Liste der Histogramme zerstoert"
#define M_WRONG_HISTPTR  	"unzulaessiger Histogramm-Pointer"

/* Klassifikatoren */
#define M_CLASSPTR      	"Nicht initialisierte Klassifikatorstruktur"
#define M_TEACH		      	"Unvollstaendig angelernter Klassifikator"
#define M_CLASS      		"Falsche Klassenangabe"

/* Statistik */
#define M_STATIST_PTR      	"Nicht initialisierter Statistikbereich"
#define M_STAT_NOENTRY		"Keine Eintraege im Statistikbereich"
#define M_ZERO_VARIANZ		"Eine Zufallskomponente besitzt Varianz 0"
#define M_EMPTY_STAT		"Leerer Statistikbereich"

/* Kalibrierung - Rekonstruktion - multikulare BV (Triangulation) */
#define M_CALIBERR      	"Kalibrierung aus Punktkorrespondenzen nicht moeglich"
#define M_PRJZ          	"Projektionszentrum nicht berechenbar"
#define M_MISS_P        	"Anzahl Punktkorrespondenzen zu klein"
#define M_UNCALIBRATED_CAM      "Kamera ist nicht kalibriert"
#define M_TOO_FEW_CAM           "unzureichende Anzahl an Kamerastrukturen"
#define M_NO_VIEWRAY            "Sehstrahl nicht verfuegbar"
#define M_NO_DISTANCE           "Distanz Raumpunkt-Sehstrahl nicht berechenbar"
#define M_WRONG_CAMERA_POS	"Kameraposition ueberpruefen!"
#define M_DIFF_PRJCENTERS	"Unterschiedliche Projektionszentren"
#define M_NO_ELLIPSE		"keine Ellipse"
#define M_INFINIT_POINT		"Abbildung in unendlich fernen Punkt"
/************************************************
 *                                              *
 *    Konstanten, Fehler- und Rueckkehr-Kodes   *
 *                                              *
 ************************************************/

#ifndef ERROR
  #define ERROR           -255
#endif

#define WRONG_PARAM		-1
#define NO_GRAPHIC		-2
#define TOO_LARGE		-3
#define NOT_FOUND		-4
#define NO_MEM			-5
#define FILE_NOT_FOUND		-6
#define WRONG_FILE		-7
#define WRONG_WINDOW    	-8
#define	NO_SOLUTION     	-9
#define	VARIOUS_SOLUTION   	-10
#define	NUM_INSTABILITY    	-11
#define INVALID_STRUCTURE  	-12
#define INTERN_ERROR 		-13
#define NO_UNIQUE_SOLUTION 	-14
#define WRONG_POINTER	 	-15
#define WRONG_STARTPOINT	-16
#define CHAIN_DESTROYED         -17
#define DIVISION_ZERO           -18
#define ZERO_VECTOR             -19 
#define MISSING_DATA            -20
#define WRONG_MODE		-21
#define CONOCLOSED		-22
#define WRONG_PRECISION         -23
#define WRONG_VECTOR		-24
#define POINT_IDENTIC		-25
#define WRONG_TRANS             -26
#define INFINIT_POINT		-27	
#define NO_ELLIPSE		-28
#define WRONG_MATRIX            -29
#define MAT_NO_COMPAT           -30
	
#define OK		0
#define INVALID		1

#ifndef TRUE
  #define TRUE		(1==1)	/* Boolsche Konstanten */
#endif
#ifndef FALSE
  #define FALSE		(1==0)
#endif

#define OFF		0       /* Schalter fuer Visualisierung .. */
#define ON		1
#define ENABLE		2
#define DISABLE		3
#define REORG		4
#define RESET		5
#define UPDATE		6
#define OVERLAY		7
#define SET             8
#ifndef RGB
 #define RGB		9
#endif
#define Rgb		10
#define rGb		11
#define rgB		12

#define OVERLAY0        13
#define OVERLAY1        14
#define OVERLAY2        15

/* Informations-Konstante */

#define DEFAULT		0
#define MAXX		1
#define MAXY		2
#define TABLE_LEN	3
#define REAL_TIME_COLOR 4
#define REAL_TIME_ZOOM  5
#define CHANNEL		6
#define VIRTUAL_X	7
#define VIRTUAL_Y	8
#define IMAGES		9
#define TRUECOLOR	10
#define ALL		11
#define WNDW		12

#define EQUAL		13
#define NOT_EQUAL	14
#define NORMAL          15 /* Hintransformation */
#define INVERS          16 /* Ruecktransformation */
#define REAL		17 /* reelle Loesung */
#define COMPLEX		18 /* komplexe Loesung */
#define MAXVAL          19

#define NO_REGULAR	30 /* keine reguleare Matrix */
#define ZERO            31 /* Nullmatrix */ 
#define UNIT		32 /* Einheitsmatrix */
#define ORTHONORMAL	33 /* Orthonormalmatrix */
#define DIAGONAL	34 /* Diagonalmatrix */
#define SYMMETRIC	35 /* Symmetr. Matrix */

#define EXIST		40
#define NO_EXIST	41
#define DEGENERATE      42 /* entartete Kurve */
#define LINE		43 /* eine Gerade */
#define CUTTING_LINES   44 /* sich schneidende Geraden */
#define PARALLEL_LINES  45 /* parallele Geraden */
#define PARABEL		46
#define HYPERBEL        47
#define ELLIPSE         48
#define INTERPOL	49 /* Interpolation */
#define IDENTICAL	50
#define PARALLEL	51
#define ROW		52
#define NOFILL		53

/*************************************************************************
complex.h
*************************************************************************/

#ifndef _COMPLEX_H
#define _COMPLEX_H

/***** Strukturen *****/
typedef struct Complex_
 {
 double re,im;
 } Complex;

int PrintCom
  (char *str,Complex c);
Complex	*SetCom
  (double re, double im, Complex *c);
Complex *MoveCom
  (Complex *c1, Complex *c2);
double DistCom
  (Complex *c1, Complex *c2);
Complex *AddCom
  (Complex *c1, Complex *c2, Complex *c3);
Complex *SubCom
  (Complex *c1, Complex *c2, Complex *c3);
Complex *MulCom
  (Complex *c1, Complex *c2, Complex *c3);
Complex *ScaleCom
  (Complex *c1, double fac , Complex *c2);
Complex *DivCom
  (Complex *c1, Complex *c2, Complex *c3);
Complex *InvertCom
  (Complex *c1, Complex *c2);
Complex	*KonjCom	
  (Complex *c1, Complex *c2);
int ConvComPolar
  (Complex *c, double *rad, double *arcus);
Complex *ConvPolarCom
  (double rad, double arcus, Complex *c);
#endif

/*************************************************************************
numbase.h
*************************************************************************/

double	Arcus		(double degree);
double	Degree		(double arcus);
int	Max		(int val1, int val2);
int	Min		(int val1, int val2);
double	MaxD		(double val1, double val2);
double	MinD		(double val1, double val2);
double	Sqr 		(double val);
double	Cub 		(double val);
double	CubRoot		(double val);
int	Sign		(int val);
double	SignD		(double val);
void	ChangeI		(int *i1,int *i2);
void 	ChangeD		(double *d1,double *d2);
double  Round		(double val);
void	Randomize	(void);
int	Random		(int val);
double	RandomD		(void);
double	GaussRandom	(double sigma);

double QuadrFunc(double *par,double x,double y);


/* Numerische Konstanten */

#define M_E         2.7182818284590452354
#define M_LOG2E     1.4426950408889634074
#define M_LOG10E    0.43429448190325182765
#define M_LN2       0.69314718055994530942
#define M_LN10      2.30258509299404568402
#define M_PI        3.14159265358979323846
#define M_PI_2      1.57079632679489661923
#define M_1_PI      0.31830988618379067154
#define M_PI_4      0.78539816339744830962
#define M_2_PI      0.63661977236758134308
#define M_2_SQRTPI  1.12837916709551257390
#define M_SQRT2     1.41421356237309504880
#define M_SQRT1_2   0.70710678118654752440

#define PI  M_PI
#define PI2  M_PI_2

/*************************************************************************
matrix.h
*************************************************************************/

/********************************************************************
  header zu "matrix.c" 
  funktionen zur matrixrechnung  

  neubauer 10/95     
********************************************************************/


/*** Matrizenalgebra ************************************************/
int PrintMatrix
  (char *str,double *h,int row, int col);

double *DefMatrix
  (double *m,int row,int col,int option);

double *SetMatrix
  (double *m,int row,int col,double val);

double *MoveMatrix
  (double *m1, int row, int col, double *m2);

double *AddMatrix
  (double *m1, double *m2, int row, int col, double *m3);

double *SubMatrix
  (double *m1, double *m2, int row, int col, double *m3);

double *TranspMatrix
  (double *m1, int row, int col, double *m2);

double *ScaleMatrix
  (double *m1,int row,int col,double fac,double *m2);

double *MulMatrix
  (double *m1, double *m2, int row1, int col1, int col2, double *m3);

double *InvertMatrix
  (double *m1, int dim, double *m2);

int NormMatrix
  (double *a,int row, int col, int mode, double *norm);

int DetMatrix
  (double *a,int row,int col,double *det,int *rank);

double *ChangeMatrixRow
  (double *m,int row1,int row2,int col);

double *ChangeMatrixCol
  (double *m,int col1,int col2,int row,int col);

double *OrthoMatrix
  (double *m1, int row, int col, double *m2);

int GenOrthoSystem
  (double n[3],double m[3][3]);

int IsMatrixZero
  (double *a,int row,int col,double *eps);
int IsMatrixUnit
  (double *a,int dim,double *eps);
int IsMatrixOrtho
  (double *a,int dim,double *eps);
int IsMatrixDiag
  (double *a,int dim,double *eps);
int IsMatrixSymm
  (double *a,int dim,double *eps);
int IsMatrixAnti
  (double *a,int dim,double *eps);
int IsMatrixRegular
  (double *a,int dim,double *eps);
int IsMatrixPosDef
  (double *a,int dim);

/*************************************************************************
fit.h
*************************************************************************/

/*******************************************************************
     -	Arbeit mit Punklisten		
	Approximation von Punktlisten an elementare Funktionen 
	  (Gerad, Kreis, Ellipse, quadrat. Formen)
	Approximation von Punktlisten an Segmente 
	  (Geradensegment,Kreissegment,Ellipsensegment)
********************************************************************/


/*******Datenstrukturen********************************************/
typedef struct PointList_
{
  int lng;                         /*Anzahl der Punkte*/
  double *xptr;                    /*x-Koordinaten*/
  double *yptr;                    /*y-Koordinaten*/
  double *wptr;                    /*Gewichte*/
} *PointList;

/******************************************************************/

typedef struct Segment_
{
  double          p0[2],p1[2];     /*Anfangs- und Endpunkt*/
  int             typ;             /*Segmenttyp: 1 - Geradensegment*/
				   /*            2 - Kreissegment*/
				   /*            3 - Ellipsensegment*/
  double          par[7];          /*Geradensegment: p,phi*/
				   /*Kreissegment:   xm,ym,r,psi1,psi2*/
				   /*Ellipsensegment:xm,ym,a,b,phi,psi1,psi2*/
  struct Segment_ *prev,*next;     /*Listenverkettung*/
}*Segment;
/******************************************************************/


#if 0

/*******Funktionen*************************************************/
PointList NewPointList(int len);
int PutPoint(PointList pl, int adr, double x, double y, double weight);
int FreePointList(PointList pl);
PointList ConturPointList(Contur c, int diff);
/******************************************************************/
int FitLine(PointList pl,int ad1,int ad2,int step,double *par,double *mdist,int *madr);
int FitCircle(PointList pl,int ad1,int ad2,int step,double *par,double *mdist,int *madr);
int FitEllipse(PointList pl,int ad1,int ad2,int step,double *par,double *mdist,int *madr);
/******************************************************************/
int FreeSegmentList(Segment sl);
Segment AddSegment(Segment sl,Segment s);
Segment FirstSegment(Segment sl);
Segment NextSegment(Segment sl);
/******************************************************************/
Segment FitLineSegment(PointList pl,int ad1,int ad2,int step,double *mdist,int *madr);
Segment FitCircleSegment(PointList pl,int ad1,int ad2,int step,double *mdist,int *madr);
Segment FitEllipseSegment(PointList pl,int ad1,int ad2,int step,double *mdist,int *madr);

#endif

/*************************************************************************
matdef.h
*************************************************************************/
typedef struct Matrix_
{
  int type;
  int rsize,csize;
  double **data;
  int **datai;
  unsigned char **datac;
  
  struct Matrix_ *prev,*next;
}  *Matrix;

#define MAT_DOUBLE 0
#define MAT_INT 1
#define MAT_CHAR 2
