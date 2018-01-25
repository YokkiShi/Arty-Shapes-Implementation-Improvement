/****************************************************************************
 * a severly bastardised version of ICE from Universitat Jena
 * developed by Prof. K. Voss, Herbert Suesse, et.al.  
 * plus some separate code from Herbert Suesse
 *
 * this fragment fits shapes using moments:
 *      circle 
 *      ellipse 
 *      triangle 
 *      rectangle
 *      affine superellipse
 *
 * there are further shapes that could be fit, they need to be converted from the
 * code in herbert-extra1.c and herbert-extra2.c
 *
 * probably best if the data is made into a dense pixel list first
 *
 * Paul Rosin
 * Cardiff University
 * March 2000
 ****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#if 0
#include "defs.h"
#include "complex.h"
#include "numbase.h"
#include "matrix.h"
#include "fit.h"
#endif
#include "fit_triangle.h"

void Message();

double *TransPoint();

int e_code;        /* error code */
int e_message=1;     /* Flag for displaying the messages */
                     /* set to 1 to stop program halting - PLR 6/2005*/
char e_msg[80];

void SetOk(void)
{
  e_code=OK;
}

/*******************************************************************/
#define FNAME "NewPointList"
PointList NewPointList(int lng)
{
  PointList pl;

  pl=(PointList) malloc(sizeof(struct PointList_));
  if (pl==0)
  {
    Message(FNAME,M_NO_MEM,ERROR);
    return NULL;
  }
  pl->lng=lng;
  pl->xptr=(double*)malloc(lng*sizeof(double));
  pl->yptr=(double*)malloc(lng*sizeof(double));
  pl->wptr=(double*)malloc(lng*sizeof(double));
  SetOk();
  return pl;
}
#undef FNAME

/************************************************************/
/*  Fitting a circle to an object                           */
/************************************************************/
int FitCircleMoments(double m[],double *x0,double *y0,double *radius)
{
    double sx,sy;
    double mtrans[15];
    double mtrop[15];
    double alpha;
    sx=m[1]/m[0];
    sy=m[2]/m[0];

    TranslateMoments(m,-sx,-sy,mtrans);
    alpha=1.0/sqrt(mtrans[0]);
    ScaleMoments(mtrans,alpha,alpha,mtrop);
/************************************************************************/
/* On mtrop the moments of the object are in standard position          */
/* m00=1, m10=m01=0                                                     */
/************************************************************************/
/*  Now to it a set in standard position is to be adapted */
/************************************************************************/
    *x0=sx;*y0=sy;
    *radius=1.0/sqrt(PI)/alpha;
    return(0);
}

/*************************************************************/
/*  Fitting an ellipse to an object                          */
/*************************************************************/
int FitEllipseMoments(double m[],double *ell_par, double *ellipticity)
{
    double sx,sy;
    double ell_koef[6];
    double mtrans[15],det;
    int type;
    /* added for ellipticity - PLR */
    int i;
    double sum;
    double maf[21];
    double atr[3][3];
    double me[15];

    sx=m[1]/m[0];
    sy=m[2]/m[0];
    TranslateMoments(m,-sx,-sy,mtrans);

    det=mtrans[5]*mtrans[3]-mtrans[4]*mtrans[4];
    ell_koef[0]=mtrans[5];
    ell_koef[1]=mtrans[3];
    ell_koef[2]=0;  
    ell_koef[3]=0;  
    ell_koef[4]=-2*mtrans[4];
    ell_koef[5]=-4*det/mtrans[0];         

    FeatureQuadrFunc(ell_koef,ell_par,&type);
    if (type!=ELLIPSE) {printf("\n Error - non ellipse fitted");return(1);}
    ell_par[0]=sx;
    ell_par[1]=sy;

    /* calculate ellipticity - PLR */

    /* affine normalised moments of data */
    AffinIterateMoments(m,maf,atr);

    /* normalised moments of ellipse */
    /* moments in order: m00,m10,m01,m20,m11,m02,m30,m21,m12,m03,m40,m31 */
    for (i=0;i<15;++i)
        me[i] = 0;
    me[0] = 2.0*sqrt(PI);
    me[3] = me[5] = 1;
    me[12] = 1.0 / (3.0 * sqrt(PI));
    me[10] = me[14] = 1.0 / sqrt(PI);

    /* RMS difference in moments */
    sum = 0;
    for (i=6;i<15;i++) {
       sum += (me[i]-maf[i])*(me[i]-maf[i]);
       /***
       moments m00 and m11 doesn't seem to match exactly - PLR
       printm(i);
       printf("%f\n",(me[i]-maf[i])*(me[i]-maf[i]));
       ***/
   }
    sum = sqrt(sum);
    *ellipticity = sum;

    return(0);
}

/******************************************************************/
/***********************************************************************/
/****************************************************************************/
/* Es wird an das Objekt mit den Momenten m ein Dreieck mit Hilfe           */
/* von Flaechenmomenten und darauf basierenden Standardlagen angepaﬂt       */
/****************************************************************************/
#define FNAME "FitTriangle"
int FitTriangle(double m[15],double p[3][2],double *triangularity)
{
  double tr[3][3],mh[15];
  double d[3][2]=
  {
    {-1.3642616, 0.7876568},
    { 0.0000000,-1.5753136},
    { 1.3642616, 0.7876568}
  };
  int i;
  PointList plm;
  double mm[15],pc0[2],sum;
  int flag;

  *triangularity = -9999;
  flag = PolyNormMoments(m,mh,tr);
  if (flag == ERROR)
      return;
  InvertTrans(tr);
  for(i=0;i<3;i++)
  {
    p[i][0]=d[i][0]*tr[0][0]+d[i][1]*tr[0][1]+tr[0][2];
    p[i][1]=d[i][0]*tr[1][0]+d[i][1]*tr[1][1]+tr[1][2];
  }

  /* calculate triangularity - PLR */
  plm = NewPointList(3);
  plm->xptr[0] = -1.3642616; plm->yptr[0] =  0.7876568;
  plm->xptr[1] =  0.0000000; plm->yptr[1] = -1.5753136;
  plm->xptr[2] =  1.3642616; plm->yptr[2] =  0.7876568;
  MomentPolygon(plm,mm,pc0);
  sum = 0;
  for (i=6;i<15;++i) {
     sum += (mm[i]-mh[i])*(mm[i]-mh[i]);
     /***
     moment m00 doesn't seem to match exactly - PLR
     printm(i);
     printf("%f\n",(mm[i]-mh[i])*(mm[i]-mh[i]));
     ***/
  }
  sum = sqrt(sum);
  *triangularity = sum;

  return(OK);
}
#undef FNAME

int FitRectangle(double m[15],double p[4][2],double *rectangularity)
{
/************************************************************/
/* Fitting of a rectangle to an object */
/* Shift, rotation, isotropic scaling                */
/* Canonical frame wird anistrop skaliert und ¸ber dieses    */
/* Verh‰ltnis optimiert, da f¸r Rechtecke keine Gruppe       */
/* existiert                                                 */
/*************************************************************/
double sx,sy,b,b_opt,s[2];
double min1,sum1;
double mtrans[15];
double mrot[15];
double mtrop[15],mp_h[15],mp_opt_n[15];
double mp_h_n[15];
double sinphi,cosphi,tanphi,phi;
double alpha,beta;
double norm;
double tr[3][3],tr_r[3][3],p1[2],p2[2];
char diskret;
PointList pl;
int isx,isy,i;

sx=m[1]/m[0];
sy=m[2]/m[0];
isx=(int)sx;
isy=(int)sy;
TranslateMoments(m,-sx,-sy,mtrans);
if (fabs(mtrans[11])<=1.0E-13 || fabs(mtrans[13])<=1.0E-13)
  {
    sinphi=0;cosphi=1;phi=0;
  }
else
  {
    b=(mtrans[5]-mtrans[3])/mtrans[4];
    tanphi=-b/2.0+sqrt(b*b/4.0+1.0);
    phi=atan(tanphi);
    norm=sqrt(tanphi*tanphi+1.0);
    sinphi=tanphi/norm;
    cosphi=1.0/norm;
  }  
RotateMoments(mtrans,cosphi,sinphi,mrot);
alpha=1.0/sqrt(mrot[0]);
ScaleMoments(mrot,alpha,alpha,mtrop);

/************************************************************************/
/* On mtrop the moments of the object are in standard position             */
/* m00=1, m10=m01=0,m11=0                                               */
/************************************************************************/
/*
 * Now to it a rectangle in standard position is to be
 * adapted to Wehlen a rectangle in
 * axle-parallel position with sidelong 1 and b in addition
 * is for a b the normalized moments to be
 * calculated and the b with the best Fit to the normalized
 * position of the object to be determined
*/

    diskret='h';
    b=0.0;min1=1.0E20;
    pl=NewPointList(4);
    b_opt=1.0;
do
    {
      if (b>1.5)
        {b+=0.1;}
      else
        {b+=0.01;}
    pl->xptr[0]=b/2.0;
    pl->yptr[0]=0.5;
    pl->xptr[1]=-b/2.0;
    pl->yptr[1]=0.5;    
    pl->xptr[2]=-b/2.0;
    pl->yptr[2]=-0.5;    
    pl->xptr[3]=b/2.0;
    pl->yptr[3]=-0.5;    
    MomentPolygon(pl,mp_h,s);
    beta=1.0/sqrt(b);
    ScaleMoments(mp_h,beta,beta,mp_h_n);
  sum1=0.0;
for (i=3;i<15;++i) {
   sum1+=(mtrop[i]-mp_h_n[i])*(mtrop[i]-mp_h_n[i]);
   /***
   printm(i);
   printf("%f\n",(mtrop[i]-mp_h_n[i])*(mtrop[i]-mp_h_n[i]));
   ***/
  }
sum1=sqrt(sum1);
  if (sum1<min1)
    {
      min1=sum1;b_opt=b;
      for (i=3;i<15;++i) mp_opt_n[i]=mp_h_n[i];
    }  
    } 
while (b<=100.0);

    /* return rectangularity measure - PLR */
    *rectangularity = min1;

/********************************************************************/
    InitTrans(tr);
    ShiftTrans(-sx,-sy,tr);
    RotTrans(0.0,0.0,phi,tr);
    ScaleTrans(0.0,0.0,alpha,alpha,tr);
    InvertTrans(tr);
    beta=1.0/sqrt(b_opt);
    InitTrans(tr_r);
    ScaleTrans(0.0,0.0,beta,beta,tr_r);
/****************************************************************************/
    p1[0]=b_opt/2.0;p1[1]=0.5;
    TransPoint(p1,tr_r,p2);
    TransPoint(p2,tr,p1);
    p[0][0]=p1[0];
    p[0][1]=p1[1];
/***************************/
    p1[0]=-b_opt/2.0;p1[1]=0.5;
    TransPoint(p1,tr_r,p2);
    TransPoint(p2,tr,p1);
    p[1][0]=p1[0];
    p[1][1]=p1[1];
/*****************************/    
    p1[0]=-b_opt/2.0;p1[1]=-0.5;
    TransPoint(p1,tr_r,p2);
    TransPoint(p2,tr,p1);
    p[2][0]=p1[0];
    p[2][1]=p1[1];
/*****************************/                    
    p1[0]=b_opt/2.0;p1[1]=-0.5;
    TransPoint(p1,tr_r,p2);
    TransPoint(p2,tr,p1);
    p[3][0]=p1[0];
    p[3][1]=p1[1];
    return(0);
} 

/****************************************************************/
/* Transformation of the moments with translation               */
/****************************************************************/
#define FNAME "TranslateMoments"
int TranslateMoments(double m1[15],double x,double y,double m2[15])
{
  int i;
  double m[15],*ms,*md;
  double x2,y2,x3,y3,x4,y4;
  if (m1==m2)
  {
    ms=m;
    md=m1;
    for(i=0;i<15;i++) ms[i]=md[i];
  }
  else
  {
    ms=m1;
    md=m2;
    for(i=0;i<15;i++) md[i]=ms[i];
  }
  x2=x*x;x3=x2*x;x4=x3*x;
  y2=y*y;y3=y2*y;y4=y3*y;
  md[1]+=x*ms[0];
  md[2]+=y*ms[0];
  md[3]+=2*x*ms[1]+x2*ms[0];
  md[4]+=y*ms[1]+x*ms[2]+x*y*ms[0];
  md[5]+=2*y*ms[2]+y2*ms[0];
  md[6]+=3*x*ms[3]+3*x2*ms[1]+x3*ms[0];
  md[7]+=2*x*ms[4]+y*ms[3]+x2*ms[2]+2*x*y*ms[1]+x2*y*ms[0];
  md[8]+=x*ms[5]+2*y*ms[4]+2*x*y*ms[2]+y2*ms[1]+x*y2*ms[0];
  md[9]+=3*y*ms[5]+3*y2*ms[2]+y3*ms[0];
  md[10]+=4*x*ms[6]+6*x2*ms[3]+4*x3*ms[1]+x4*ms[0];
  md[11]+=3*x*ms[7]+y*ms[6]+3*x2*ms[4]+3*x*y*ms[3]+x3*ms[2]+3*x2*y*ms[1]+x3*y*ms[0];
  md[12]+=2*x*ms[8]+2*y*ms[7]+x2*ms[5]+4*x*y*ms[4]+y2*ms[3]+\
    2*x2*y*ms[2]+2*x*y2*ms[1]+x2*y2*ms[0];
  md[13]+=x*ms[9]+3*y*ms[8]+3*x*y*ms[5]+3*y2*ms[4]+3*x*y2*ms[2]+y3*ms[1]+x*y3*ms[0];
  md[14]+=4*y*ms[9]+6*y2*ms[5]+4*y3*ms[2]+y4*ms[0];
  return(OK);
}
#undef FNAME

/****************************************************************/
/* Transformation der Momente bei Rotation                      */
/****************************************************************/
#define FNAME "RotateMoments"
int RotateMoments(double m1[15],double c,double s, double m2[15])
{
  int i;
  double c2,s2;
  double mx[15],*m,*mr;
  if (m1==m2)
  {
    m=mx; mr=m2;
    for(i=0;i<15;i++) m[i]=mr[i];
  }
  else
  {
    m=m1; mr=m2;
    for(i=0;i<15;i++) mr[i]=m[i];
  }
  c2=c*c;
  s2=s*s;
  mr[3]=c2*m[3]+s*(-2*c*m[4]+s*m[5]);
  mr[4]=c*s*(m[3]-m[5])+(c2-s2)*m[4];
  mr[5]=s*(s*m[3]+2*c*m[4])+c2*m[5];
  mr[6]=c2*(c*m[6]-3*s*m[7])+s2*(3*c*m[8]-s*m[9]);
  mr[7]=c2*(c*m[7]+s*(m[6]-2*m[8]))+s2*(s*m[8]+c*(m[9]-2*m[7]));
  mr[8]=c2*(c*m[8]-s*(m[9]-2*m[7]))-s2*(s*m[7]-c*(m[6]-2*m[8]));
  mr[9]=s2*(s*m[6]+3*c*m[7])+c2*(3*s*m[8]+c*m[9]);
  mr[10]=c2*(c*(c*m[10]-4*s*m[11])+6*s2*m[12])+s2*s*(-4*c*m[13]+s*m[14]);
  mr[11]=c2*(c*(c*m[11]+s*(m[10]-3*m[12]))+s2*3*(-m[11]+m[13]))+\
    s2*(s*(c*(3*m[12]-m[14])-s*m[13]));
  mr[12]=c2*(s2*(m[10]-4*m[12]+m[14])+c*(2*s*(m[11]-m[13])+c*m[12]))+\
    s2*s*(s*m[12]-2*c*(m[11]-m[13]));
  mr[13]=s2*(s*(c*(m[10]-3*m[12])-s*m[11])+3*c2*(m[11]-m[13]))+\
    c2*(c*(s*(3*m[12]-m[14])+c*m[13]));
  mr[14]=s2*(s*(s*m[10]+4*c*m[11])+6*c2*m[12])+c2*c*(4*s*m[13]+c*m[14]);
  return(OK);
}

#undef FNAME
/***********************************************************/
/* Standardisation of the moments with the polynomial method (Voss) */
/* X-shearing, y-shearing, scaling      */
/***********************************************************/
#define FNAME "PolyNormMoments"
int PolyNormMoments(double m1[15],double m2[15],double atr[3][3])
{
  double mx[3][15];
  Complex sol[3];
  double a[3],b[3],c[3],d[3],s,smin,xs,ys,h;
  int i,imin,k,l;
  double sum;

  xs=0;ys=0;
  if (m1[0]<0)
    for(i=0;i<15;i++) m2[i]=-m1[i];
  else
    for(i=0;i<15;i++) m2[i]=m1[i];
  /* translation standardisation */
  if (m2[1]!=0 || m2[2]!=0)
  {
    xs=-m2[1]/m2[0]; ys=-m2[2]/m2[0];
    TranslateMoments(m2,xs,ys,m2);

    /* added later by Herbert - PLR */
    sum = 0.0;
    for (i=6;i<10;++i)
        sum += fabs(m2[i]);
    
    /* this is a normalization to the area = relative error */
    sum = sum/sqrt(m2[0]*m2[0]*m2[0]*m2[0]*m2[0]);
    if (sum<1e-6) {
        fprintf(stderr,"triangle fitting failed\n");
        return(ERROR);
    }

  }
  if (fabs(m2[9])<1e-15)
  {
    if (m2[8]!=0)
    {
      h=Sqr(m2[7]/(2*m2[8]))-m2[6]/(3*m2[8]);
      if (h<0) return(ERROR);
      sol[0].re=-m2[7]/(2*m2[8])+sqrt(h);sol[0].im=0;
      sol[1].re=-m2[7]/(2*m2[8])-sqrt(h);sol[1].im=0;
      k=0;l=2;
    }
    else
    {
      if (m2[7]==0) return(ERROR);
      else
      {
    sol[0].re=-m2[6]/(3*m2[7]);
    sol[0].im=0;
    k=0;l=1;
      }
    }
  }
  /*m03!=0 --> Equation 3. Degree */
  else
  {
    Root3(3*m2[8]/m2[9],3*m2[7]/m2[9],m2[6]/m2[9],sol);
    /*Lˆsungen mit gleichem Realteil aussortieren*/
    k=0;l=3;
    if (sol[0].re==sol[1].re) k=1;
    if (sol[1].re==sol[2].re) l=2;
  }
  /*Normallagen fur alle unterschiedlichen reellen Lˆsungen berechnen*/
  smin=1e32;
  imin=-1;
  for(i=k;i<l;i++)
  {
    if (sol[i].im==0)
    {
      a[i]=sol[i].re;
      XShearMoments(m2,a[i],mx[i]);
      b[i]=-mx[i][4]/mx[i][3];
      YShearMoments(mx[i],b[i],mx[i]);
      if (mx[i][5]*mx[i][3] > 1e-10)
      {
    c[i]=pow(mx[i][5]/(mx[i][3]*mx[i][3]*mx[i][3]),0.125);
    d[i]=pow(mx[i][3]/(mx[i][5]*mx[i][5]*mx[i][5]),0.125);
    ScaleMoments(mx[i],c[i],d[i],mx[i]);
    /*Normallage mit m40+m04 --> min. ausw‰hlen*/
    s=mx[i][10]+mx[i][14];
    if (s<smin)
    {
      smin=s;
      imin=i;
    }
      }
    }
  }
  if (imin<0)
  {
    return(ERROR);
  }
  for(i=0;i<15;i++) m2[i]=mx[imin][i];
  /*affine Transformation fur augew‰hlte Normallage*/
  InitTrans(atr);
  atr[0][0]=c[imin];
  atr[0][1]=a[imin]*c[imin];
  atr[1][0]=b[imin]*d[imin];
  atr[1][1]=(1+a[imin]*b[imin])*d[imin];
  if (m2[7]<0)
  {
    atr[0][0]=-atr[0][0];atr[0][1]=-atr[0][1];
    atr[1][0]=-atr[1][0];atr[1][1]=-atr[1][1];
    m2[7]=-m2[7];m2[8]=-m2[8];m2[9]=-m2[9];
  }
/*   NormMoments(m2,atr); */
  atr[0][2]=atr[0][0]*xs+atr[0][1]*ys;
  atr[1][2]=atr[1][0]*xs+atr[1][1]*ys;
  return(OK);
}
#undef FNAME

/* Initialize a transformation matrix (unit matrix) */
int InitTrans(double Trans[3][3])
{
  DefMatrix((double *)Trans,3,3,UNIT);
  return OK;
}  

/*******************************************************************/
#define FNAME "ShiftTrans"
/* Premultiplikation von "Trans" mit Verschiebung (x0,y0)*/
int ShiftTrans
  (double x0, double y0, double Trans[3][3])
  {
  double m[3][3];
  double eps=1e-10;
  double fac;

  if (Trans==NULL )
  {
    Message(FNAME,M_WRONG_PTR,WRONG_POINTER);
    return ERROR;
  }
  DefMatrix((double *)m,3,3,UNIT);
  m[0][2]=x0;m[1][2]=y0;
  MulMatrix((double *)m,(double *)Trans,3,3,3,(double *)Trans);
  /* Normierung nach Summe der Quadrate der unteren Zeile */
  fac=Sqr(Trans[2][0])+Sqr(Trans[2][1])+Sqr(Trans[2][2]);
  if (fac<eps)
    {
    Message(FNAME,M_NO_PROJ,ERROR);
    return ERROR;
    }
  ScaleMatrix((double *)Trans,3,3,fac,(double *)Trans);
  SetOk();return OK;
  }  

#undef FNAME

/*******************************************************************/
#define FNAME "RotTrans"
/* Premultiplikation von "Trans" mit Rotation um Punkt (x0,y0) */
int RotTrans
  (double x0, double y0, double phi, double Trans[3][3])
  {
  double rot[3][3];
  double eps=1e-10;
  double fac;

  if (Trans==NULL)
  {
    Message(FNAME,M_WRONG_PTR,WRONG_POINTER);
    return ERROR;
  }

  InitTrans(rot);
  rot[0][0]= cos(phi); rot[0][1]=-sin(phi);
  rot[1][0]= sin(phi); rot[1][1]= cos(phi);

  rot[0][2]= -cos(phi)*x0+sin(phi)*y0+x0; 
  rot[1][2]= -sin(phi)*x0-cos(phi)*y0+y0;

  MulMatrix((double *)rot,(double *)Trans,3,3,3,(double *)Trans);
  /* Normierung nach Summe der Quadrate der unteren Zeile */
  fac=Sqr(Trans[2][0])+Sqr(Trans[2][1])+Sqr(Trans[2][2]);
  if (fac<eps)
    {
    Message(FNAME,M_NO_PROJ,ERROR);
    return ERROR;
    }
  ScaleMatrix((double *)Trans,3,3,fac,(double *)Trans);
  SetOk();return OK;
  }
#undef FNAME

/*******************************************************************/
#define FNAME "ScaleTrans"
/* Premultiplikation von "Trans" mit Skalierung um Punkt (x0,y0) */
int ScaleTrans
  (double x0, double y0, double facx, double facy, double Trans[3][3])
  {
  double scal[3][3];
  double eps=1e-10;
  double fac;

  if (Trans==NULL)
  {
    Message(FNAME,M_WRONG_PTR,WRONG_POINTER);
    return ERROR;
  }

  InitTrans(scal);
  scal[0][0]= facx; scal[0][2]=-facx*x0+x0; 
  scal[1][1]= facy; scal[1][2]=-facy*y0+y0;
  MulMatrix((double *)scal,(double *)Trans,3,3,3,(double *)Trans);
  /* Normierung nach Summe der Quadrate der unteren Zeile */
  fac=Sqr(Trans[2][0])+Sqr(Trans[2][1])+Sqr(Trans[2][2]);
  if (fac<eps)
    {
    Message(FNAME,M_NO_PROJ,ERROR);
    return ERROR;
    }
  ScaleMatrix((double *)Trans,3,3,fac,(double *)Trans);
  SetOk();return OK;
  }
#undef FNAME

#define FNAME "InvertTrans"
/* Invertieren einer Transformationsmatrix */
int InvertTrans(double Trans[3][3])
  {
  double m[3][3];
  double eps=1e-10;
  double fac;

  MoveMatrix((double *)Trans,3,3,(double *)m);
  InvertMatrix((double *)m,3,(double *)Trans);
  /* Normierung nach Summe der Quadrate der unteren Zeile */
  fac=Sqr(Trans[2][0])+Sqr(Trans[2][1])+Sqr(Trans[2][2]);
  if (fac<eps)
    {
    Message(FNAME,M_NO_PROJ,ERROR);
    return ERROR;
    }
  ScaleMatrix((double *)Trans,3,3,1/fac,(double *)Trans);
  SetOk();return OK;
  }
#undef FNAME

#define FNAME "MulMatrix"
double *MulMatrix
  (double *m1, double *m2, int row1, int col1, int col2, double *m3)
{
int i,j,k;
double *mh,*mptr;
if (m1==NULL || m2==NULL)
  {
  Message("MulMatrix", M_WRONG_VECTOR, WRONG_VECTOR);
  return NULL;
  }
if (row1<1 || col1<1 || col2<1)
  {
  Message("MulMatrix", M_WRONG_PARAM, WRONG_PARAM);
  return NULL;
  }
if (m3==NULL)
  {
  mptr=(double *)malloc(row1*col2*sizeof(double));
  for (i=0;i<row1;i++)
    {
    for (j=0;j<col2;j++)
      {
      mptr[i*col2+j]=0;
      for (k=0;k<col1;k++)
        mptr[i*col2+j]+=m1[i*col1+k] * m2[k*col2+j];
      }
    }
  return mptr;
  }
mh=(double *) malloc( row1*col2*sizeof(double) );
if (mh==NULL)
  {  Message("MulMatrix", M_NO_MEM, NO_MEM);
     return NULL;
  }
for (i=0;i<row1;i++)
  {
  for (j=0;j<col2;j++)
    {
    mh[i*col2+j]=0;
    for (k=0;k<col1;k++)
      mh[i*col2+j]+=m1[i*col1+k] * m2[k*col2+j];
    }
  }
MoveMatrix( (double *)mh, row1, col2, (double *)m3);
free(mh);
return m3;
}
#undef FNAME 

#define FNAME "TransPoint"
double *TransPoint(double p1[2],double t[3][3],double p2[2])
{
  double ph[3];
  double *ptr;

  if ((p1==NULL) || (t==NULL))
  {
    Message(FNAME,M_WRONG_PTR,WRONG_POINTER);
    return NULL;
  }
  ph[0]=p1[0];ph[1]=p1[1];ph[2]=1;
  MulMatrix((double *)t,(double *)ph,3,3,1,(double *)ph);
  if (fabs(ph[2])<1e-10)
  {
    Message(FNAME,M_WRONG_TRANS,WRONG_TRANS);
    return NULL;
  }
  if (p2==NULL)
  {
    ptr=(double *)malloc(2*sizeof(double));
    ptr[0]=ph[0]/ph[2];  
    ptr[1]=ph[1]/ph[2];  
    SetOk();
    return ptr;
  }
  p2[0]=ph[0]/ph[2];  
  p2[1]=ph[1]/ph[2];  
  SetOk();
  return p2;
}
#undef FNAME

/*******************************************************************/
void ConvCartesPolar (double p[2], double *rad, double *arcus)
{
  *rad=sqrt (Sqr(p[0]) + Sqr(p[1]) );
  *arcus=atan2(p[1],p[0]);
}

/*******************************************************************/
#define FNAME "Root3"
#define DEBUG 0

int Root3 (double p2, double p1, double p0, Complex c[3])
{
  double dis,h,b2,b3;
  double b3sign,r,x,y,phi,u,v;
  double hv[2];
  double eps=1e-10;

  #if DEBUG
    printf("p2: %f p1: %f p0: %f\n",p2,p1,p0);
  #endif

  b2=-(p2*p2/3.0)+p1;
  h =-p1*p2/3.0;
  b3=(2.0*p2*p2*p2/27.0) + p0 + h;
  h =b2*b2*b2/27.0;
  dis=(b3*b3/4.0) + h;

  #if DEBUG
    printf("b2 %f b3  %f dis %f\n",b2, b3, dis);
  #endif
  if (dis<0)
  {
    if (fabs(b3)<eps) b3sign=1;
    else b3sign=SignD(b3);
    r=sqrt(fabs(b2/3.0))*b3sign;
    x=b3/2.0/r/r/r;
    y=sqrt(MaxD(-x*x+1,0));
    hv[0]=x; hv[1]=y;
    ConvCartesPolar(hv,&h,&phi);
    #if DEBUG 
      printf("x: %f y: %f phi: %f ",x,y,phi);
    #endif
    h=p2/3.0;
    c[0].re = -2*r*cos(phi/3.0)-h;
    c[0].im =  0.0;
    c[1].re =  2*r*cos((-phi+M_PI)/3.0)-h;
    c[1].im =  0.0;
    c[2].re =  2*r*cos((phi+M_PI)/3.0)-h;
    c[2].im =  0.0;
    SetOk();
    return OK;
  }
  dis=sqrt(dis); 
  h=-b3/2.0+dis; u=CubRoot(h);
  h=-b3/2.0-dis; v=CubRoot(h);
  h=p2/3.0;
  c[0].re =  u+v-h;
  c[0].im =  0.0;
  c[1].re = -(u+v)/2.0-h;
  c[1].im =  (u-v)/2.0*sqrt(3.0);
  c[2].re =  c[1].re;
  c[2].im = -c[1].im;
  SetOk();
  return OK;
}
#undef FNAME 
#undef DEBUG

/*******************************************************************
 numbase.c - Numerische Grundfunktionen:
*******************************************************************/

/*******************************************************************/
/*** Bogenma· eines Winkels ***/
double Arcus(double val){return val*M_PI/180;}

/*******************************************************************/
/*** Winkel aus Bogenma·  ***/
double Degree(double arcus){return arcus*180/M_PI;}

/*******************************************************************/
/*** Maximum von "int"-Werten ***/
int Max(int val1,int val2){return val1>val2?val1:val2;}

/*******************************************************************/
/*** Minimum von "int"-Werten ***/
int Min(int val1,int val2){return val1<val2?val1:val2;}

/*******************************************************************/
/*** Maximum von "double"-Werten ***/
double MaxD(double val1,double val2){return val1>val2?val1:val2;}

/*******************************************************************/
/*** Minimum von "double"-Werten ***/
double MinD(double val1,double val2){return val1<val2?val1:val2;}

/*******************************************************************/
/*** Quadrat ***/
double Sqr(double val){return(val*val);}

/*******************************************************************/
/*** Kubikwert ***/
double Cub(double val) {return(val*val*val);}

/*******************************************************************/
/*** Kubikwurzel ***/
double CubRoot(double val)
{
if (val==0.0)
  return 0.0;
return exp(log(fabs(val))/3.0)*SignD(val);
}

/*******************************************************************/
/*** Vorzeichen eines "int"-Wertes ***/
int Sign(int val)
{
if (val<0) return -1;
if (val>0) return 1;
return 0;
}

/*******************************************************************/
/*** Vorzeichen eines "double"-Wertes ***/
double SignD(double val)
{
if (val<0.0) return -1.0;
if (val>0.0) return 1.0;
return 0.0;
}

/*******************************************************************/
/* Tausch von Integerwerten */
void ChangeI(int *i1,int *i2)
{int h;h=*i1;*i1=*i2;*i2=h;}

/*******************************************************************/
/* Tausch von Doublewerten */
void ChangeD(double *d1,double *d2)
{double h;h=*d1;*d1=*d2;*d2=h;}

/*******************************************************************/
/* Runden von Doublewerten */
double Round(double val)
{
  return (floor(val+0.5));
}

/*******************************************************************/
double QuadrFunc(double *par,double x,double y)
/*
    Es wird der Funktionswert der quadratischen Funktion
        f(x,y)=axx+byy+cx+dy+exy+f
    berechnet
*/
{

  return par[0]*x*x + par[1]*y*y + par[2]*x + par[3]*y + par[4]*x*y + par[5];
}

void Message(char *loc, char *msg, int code)
{
  char ch=' ';
  if (strcmp(msg,M_0)!=0)
    {
    e_code=code;
    strcpy(e_msg,msg);
    }
  if (e_message==0)
  {
    printf("\n%s() - %s\n",loc,e_msg);
    printf("Type 'I' to ignore or any other key to abort\n");
    ch=getchar();
    if ((ch=='i')||(ch=='I')) return;
    exit(code);     
  }
}

/*******************************************************************/
int PrintMatrix(char *str,double *h,int row, int col)
{
int i,j;
double *dptr;

printf("*** %s\n",str);
dptr=h;
for (i=0;i<row;i++)
  {
  for (j=0;j<col;j++)
    printf("%f   ",*dptr++);
  printf("\n");
  }
printf("\n\n");
return OK;
}
#undef FNAME 

/*******************************************************************/
#define FNAME "MoveMatrix"
double *MoveMatrix(double *m1, int row, int col, double *m2)
{
int i;
double *mptr;
if (m1==NULL)
  {
  Message(FNAME, M_WRONG_VECTOR, WRONG_VECTOR);
  return NULL;
  }
if (row<1 || col<1)
  {
  Message(FNAME, M_WRONG_PARAM, WRONG_PARAM);
  return NULL;
  }
if (m2==NULL)
  {
  mptr=(double *)malloc(row*col*sizeof(double));
  for (i=0;i<row*col;i++) *(mptr+i) = m1[i];
  return mptr;
  }
for (i=0;i<row*col;i++)
  m2[i]=m1[i];
return m2;
}
#undef FNAME

/*******************************************************************/
#define FNAME "DefMatrix"
double *DefMatrix
  (double *m,int row,int col,int option)
  {
  int i,j;
  if (row<1 || col<1)
    {
    Message(FNAME, M_WRONG_PARAM, WRONG_PARAM);
    return NULL;
    }
  if (m==NULL) m=(double *)malloc(row*col*sizeof(double));
  switch (option)
    {
    case UNIT:
      for (j=0;j<row;j++) 
        for (i=0;i<col;i++)
          if (i==j) *(m+(j*col)+i)=1;
          else *(m+(i*col)+j)=0;
      return m;
    case ZERO:
      for (i=0;i<row*col;i++) *(m+i)=0;
      return m;
    default:
      Message(FNAME, M_WRONG_PARAM, WRONG_PARAM);
      return NULL;
    }
  }
#undef FNAME 

/*******************************************************************/
#define FNAME "ScaleMatrix"
double *ScaleMatrix(double *m1, int row, int col, double fac,double *m2)
{
int i;
double *mptr;
if (m1==NULL)
  {
  Message(FNAME , M_WRONG_VECTOR, WRONG_VECTOR);
  return NULL;
  }
if (row<1 || col<1)
  {
  Message(FNAME , M_WRONG_PARAM, WRONG_PARAM);
  return NULL;
  }
if (m2==NULL)
  {
  mptr=(double *)malloc(row*col*sizeof(double));
  for (i=0;i<row*col;i++) *(mptr+i) = fac*m1[i];
  return mptr;
  }
for (i=0;i<row*col;i++)
  m2[i]=fac*m1[i];
return m2;
}
#undef FNAME

/********************************************************************/
#define FNAME "InvertMatrix"
double *InvertMatrix (double *m1, int dim, double *m2)
  {
  double *mptr=NULL;
  int c,i,j,k,col,colh,row,rowh,offs,ret;
  double maxa,max,fh,fha;
  double epsinst=1e-40;  /* Grenze der InstabilitÑt */
  double epsnull=1e-200; /* Definition 0 */
  double *dpa,*dpc;
  int size_a;
  col=0;row=0;maxa=0;
  if (m1==NULL)
    {
    Message("InvertMatrix", M_WRONG_VECTOR, WRONG_VECTOR);
    return NULL;
    }
  if (m2==NULL) mptr=(double *)malloc(row*col*sizeof(double));
  size_a=dim*dim*sizeof(double);
  /* Anforderung des dynamischen Speichers */
  dpa=(double *) malloc(size_a);
  if (dpa==NULL)
    {
    Message("InvertMatrix",M_NO_MEM,NO_MEM);
    if (mptr!=NULL) free(mptr);
    return NULL;
    }
  dpc=(double *)malloc(size_a);
  if (dpc==NULL)
    {
    free(dpa);
    if (mptr!=NULL) free(mptr);
    Message("InvertMatrix",M_NO_MEM,NO_MEM);
    return NULL;
    }
  /* kopieren in dynamischen speicherbereich */
  memcpy(dpa,m1,size_a);
  /* fÅllen der Einheitsmatrix */
  for(i=0;i<dim*dim;i++) dpc[i]=0;
  for(i=0;i<dim;i++) dpc[i*dim+i]=1;

  #if defined(TESTM)
  printf("dynamischer speicher:\n");
  printmat(dim,dim,dpa);
  printmat(dim,dim,dpc);
  getc(stdin);
  #endif
  /* transformation des dynamischen A|C    */
  offs=0;col=0;ret=OK;
  while (col<dim&&!offs)

    /* setzen der spalte "col" zu 0 */
    {
    #if defined(TESTM)
    printf("beginn der bearbeitung spalte: %d\n",col);
    #endif
    max=0;colh=col-offs;rowh=colh;
    /* bestimmen des maximums der spalte "col" */
    row=colh;
    while (row<dim)
      {
      fha=dpa[dim*row+col];
      fh= fabs(fha);
      if  (fh>max)
    {maxa=fha;max=fh;rowh=row;}
      row++;
      }
      #if defined(TESTM)
      printf("maximum der spalte %d: %.3e\n",col,(double) max);
      #endif
      /* transformation der zeilen bezÅglich Spalte "col" */
      if (rowh!=colh)
    /* tauschen der zeilen "rowh" und "colh" */
    {
    for (c=col;c<dim;c++)
      {
      fh=dpa[dim*rowh+c];
      dpa[dim*rowh+c]=dpa[dim*colh+c];
      dpa[dim*colh+c]=fh;
      }
    for (c=0;c<dim;c++)
      {
      fh=dpc[dim*rowh+c];
      dpc[dim*rowh+c]=dpc[dim*colh+c];
      dpc[dim*colh+c]=fh;
      }
    #if defined(TESTM)
    printf("ende tauschen zeilen %d  %d\n",rowh,colh);
    printmat(dim,dim,dpa);
    printmat(dim,dim,dpc);
    getc(stdin);
    #endif
    }
      /* setzen der spalte "col" zu 0 ab zeile "colh+1" */
      if (fabs(dpa[dim*colh+col]/(*dpa))<epsinst) ret=NUM_INSTABILITY;
      if (fabs(dpa[dim*colh+col])>epsnull)
      {  i=colh+1;
     while(i<dim)
     {
     fh=dpa[dim*i+col]/maxa;
     #if defined(TESTM)
     dpa[dim*i+col]=0;
     #endif
     for (j=col+1;j<dim;j++)
       dpa[dim*i+j]-=dpa[dim*colh+j]*fh;
     for(k=0;k<dim;k++) dpc[dim*i+k]-=dpc[colh*dim+k]*fh;
     i++;
     }

       }
      else
    offs=1;

    col++;
    #if defined(TESTM)
    printf("ende transformation:\n");
    printmat(dim,dim,dpa);
    printmat(dim,dim,dpc);
    getc(stdin);
    #endif
    };
  if (offs==0)
    {
    /* berechnen der Lîsungsvektoren */
    if (m2==NULL)
      {
      for(k=0;k<dim;k++)
    {
    for (i=dim-1;i>=0;i--)
      {
      fh=0;j=dim-1;
      while (j>i)
        {
        fh+=dpa[dim*i+j]*mptr[j*dim+k];
        j--;
        }
      mptr[dim*i+k]=(dpc[i*dim+k]-fh)/dpa[dim*i+i];
      }
    }
      #if defined(TESTM)
      printmat(dim,dim,m2);
      #endif
      }
    else
      {
      for(k=0;k<dim;k++)
    {
    for (i=dim-1;i>=0;i--)
      {
      fh=0;j=dim-1;
      while (j>i)
        {
        fh+=dpa[dim*i+j]*m2[j*dim+k];
        j--;
        }
      m2[dim*i+k]=(dpc[i*dim+k]-fh)/dpa[dim*i+i];
      }
    }
      #if defined(TESTM)
      printmat(dim,dim,m2);
      #endif
      }
    }
  else
    {
    Message("InvertMatrix",M_NO_INVERS,NO_SOLUTION);
    return NULL;
    }

  free(dpa);free(dpc);
  if (ret==NUM_INSTABILITY)
    {
    Message("InvertMatrix",M_NUM_INSTABILITY,NUM_INSTABILITY);
    return NULL;
    }
  if (m2==NULL) return mptr;
  SetOk();
  return m2;
  }
#undef FNAME 

/****************************************************************/
/* Transformation der Momente bei X-Scherung                    */
/****************************************************************/
#define FNAME "XShearMoments"
int XShearMoments(double m1[15],double a,double m2[15])
{
  int i;
  double mx[15],*m,*ms;
  if (m1==m2)
  {
    m=mx; ms=m2;
    for(i=0;i<15;i++) m[i]=ms[i];
  }
  else
  {
    m=m1; ms=m2;
    for(i=0;i<15;i++) ms[i]=m[i];
  }
  ms[3]+=(m[5]*a+2*m[4])*a;
  ms[4]+=m[5]*a;
  ms[6]+=((m[9]*a+3*m[8])*a+3*m[7])*a;
  ms[7]+=(m[9]*a+2*m[8])*a;
  ms[8]+=m[9]*a;
  ms[10]+=(((m[14]*a+4*m[13])*a+6*m[12])*a+4*m[11])*a;
  ms[11]+=((m[14]*a+3*m[13])*a+3*m[12])*a;
  ms[12]+=(m[14]*a+2*m[13])*a;
  ms[13]+=m[14]*a;
  return(OK);
}
#undef FNAME
/****************************************************************/
/* Transformation der Momente bei Y-Scherung                    */
/****************************************************************/
#define FNAME "YShearMoments"
int YShearMoments(double m1[15],double b,double m2[15])
{
  int i;
  double mx[15],*m,*ms;
  if (m1==m2)
  {
    m=mx; ms=m2;
    for(i=0;i<15;i++) m[i]=ms[i];
  }
  else
  {
    m=m1; ms=m2;
    for(i=0;i<15;i++) ms[i]=m[i];
  }
  ms[4]+=m[3]*b;
  ms[5]+=(m[3]*b+2*m[4])*b;
  ms[7]+=m[6]*b;
  ms[8]+=(m[6]*b+2*m[7])*b;
  ms[9]+=((m[6]*b+3*m[7])*b+3*m[8])*b;
  ms[11]+=m[10]*b;
  ms[12]+=(m[10]*b+2*m[11])*b;
  ms[13]+=((m[10]*b+3*m[11])*b+3*m[12])*b;
  ms[14]+=(((m[10]*b+4*m[11])*b+6*m[12])*b+4*m[13])*b;
  return(OK);
}
#undef FNAME

/****************************************************************/
/*  Transformation of the moments with anisotropic scaling */
/****************************************************************/
#define FNAME "ScaleMoments"
int ScaleMoments(double m1[15],double a,double b, double m2[15])
{
  int i;
  double *m;
  double a2,a3,a4,a5;
  double b2,b3,b4,b5;
  if (m2!=m1) for(i=0;i<15;i++) m2[i]=m1[i];
  m=m2;
  a2=a*a;a3=a2*a;a4=a3*a;a5=a4*a;
  b2=b*b;b3=b2*b;b4=b3*b;b5=b4*b;
  m[0]*=a*b;
  m[1]*=a2*b;
  m[2]*=b2*a;
  m[3]*=a3*b;
  m[4]*=a2*b2;
  m[5]*=a*b3;
  m[6]*=a4*b;
  m[7]*=a3*b2;
  m[8]*=a2*b3;
  m[9]*=a*b4;
  m[10]*=a5*b;
  m[11]*=a4*b2;
  m[12]*=a3*b3;
  m[13]*=a2*b4;
  m[14]*=a*b5;
  return(OK);
}
#undef FNAME

/************************************************/
/* Calculation of the moments of a polygon      */
/************************************************/
#define FNAME "MomentPolygon"
int MomentPolygon(PointList pl,double m[15],double centre[2])
{
  int i,j;
  double x1,y1,x2,y2,x12,x22,y12,y22,a;
  for(i=0;i<15;i++) m[i]=0;
  for(i=0,j=1;i<pl->lng;i++,j++)
  {
    if (j==pl->lng) j=0;
    x1=pl->xptr[i]; y1=pl->yptr[i];
    x2=pl->xptr[j]; y2=pl->yptr[j];
    a=(x1*y2-x2*y1)/2;
    x12=x1*x1; x22=x2*x2;
    y12=y1*y1; y22=y2*y2;
    m[0]+=a;
    m[1]+=(x1+x2)*a;
    m[2]+=(y1+y2)*a;
    m[3]+=(x1*(x1+x2)+x22)*a;
    m[4]+=(y1*(x1+x1+x2)+y2*(x1+x2+x2))*a;
    m[5]+=(y1*(y1+y2)+y22)*a;
    m[6]+=(x12+x22)*(x1+x2)*a;
    m[7]+=(y1*(x1*(3*x1+2*x2)+x22)+y2*(x12+x2*(2*x1+3*x2)))*a;
    m[8]+=(x1*(y1*(3*y1+2*y2)+y22)+x2*(y12+y2*(2*y1+3*y2)))*a;
    m[9]+=(y12+y22)*(y1+y2)*a;
    m[10]+=((((x1+x2)*x1+x22)*x1+x22*x2)*x1+x22*x22)*a;
    m[11]+=((y1+4*y2)*x22*x2+(2*y1+3*y2)*x1*x22+\
        (3*y1+2*y2)*x12*x2+(4*y1+y2)*x12*x1)*a;
    m[12]+=((2*y12+6*y1*y2+12*y22)*x22+(6*(y12+y22)+8*y1*y2)*x1*x2+\
        (12*y12+6*y1*y2+2*y22)*x12)*a;
    m[13]+=((x1+4*x2)*y22*y2+(2*x1+3*x2)*y1*y22+\
        (3*x1+2*x2)*y12*y2+(4*x1+x2)*y12*y1)*a;
    m[14]+=((((y1+y2)*y1+y22)*y1+y22*y2)*y1+y22*y22)*a;
/*    m[15]+=((((x1+x2)*x1+x22)*x1+x22*x2)*x1+x22*x22)*x1+x22*x22*x2;*/
/*    m[16]+=(y1+5*y1)*x22*x22+(2*y1+4*y2)*x1*x2*x22+3*(y1+y2)*x12+x22+\*/
/*      (4*y1+2*y2)*x12*x1*x2+(5*y1+y2)*x12*x12;*/
/*    m[17]+=(2*y12+8*y1*y2+20*y22)*x22*x2+(6*y12+12*y2*(y1+y2))*x1*x22\*/
/*      (12*y1(y1+y2)+6*y22)*x12*x2+(20*y12+8*y1*y2+2*y22)*x12*x1;*/
/*    m[18]+=(2*x12+8*x1*x2+20*x22)*y22*y2+(6*x12+12*x2*(x1+x2))*y1*y22\*/
/*      (12*x1(x1+x2)+6*x22)*y12*y2+(20*x12+8*x1*x2+2*x22)*y12*y1;*/
/*    m[19]+=(x1+5*x1)*y22*y22+(2*x1+4*x2)*y1*y2*y22+3*(x1+x2)*y12+y22+\*/
/*      (4*x1+2*x2)*y12*y1*y2+(5*x1+x2)*y12*y12;*/
/*    m[20]+=((((y1+y2)*y1+y22)*y1+y22*y2)*y1+y22*y22)*y1+y22*y22*y2;*/
  }
  m[1]/=3;  m[2]/=3;
  m[3]/=6; m[4]/=12; m[5]/=6;
  m[6]/=10; m[7]/=30; m[8]/=30; m[9]/=10;
  m[10]/=15; m[11]/=60; m[12]/=180;
  m[13]/=60; m[14]/=15;
/*  m[15]/=21; m[16]/=105; m[17]/=420;*/
/*  m[18]/=420; m[19]/=105; m[20]/=21;*/
  if (fabs(m[0])<1e-20) return(ERROR);
  centre[0]=m[1]/m[0]; centre[1]=m[2]/m[0];
  return(OK);
}
#undef FNAME
/************************************************/
/* Calculation of the moments for a polygon as window of a point list  */
/************************************************/
#define FNAME "PointListMoment"
int PointListMoment(PointList pl,int a1,int a2,double m[15],double centre[2])
{
  int i,j,a3;
  double x1,y1,x2,y2;
  double x12,x22,y12,y22;
  double a;
  if (a1<0 || a2<0 || a1>=pl->lng || a2>=pl->lng) 
  {
    Message(FNAME,M_WRONG_PARAM,WRONG_PARAM);
    return(WRONG_PARAM);
  }
  for(i=0;i<15;i++) m[i]=0;
  a3=a2+1;
  i=a1; j=a1+1;
  do
  {
    if (i==pl->lng) i=0;
    if (j==a3) j=a1;
    if (j==pl->lng) j=0;
    x1=pl->xptr[i];
    y1=pl->yptr[i];
    x2=pl->xptr[j];
    y2=pl->yptr[j];
    a=(x1*y2-x2*y1)/2;
    x12=x1*x1; x22=x2*x2;
    y12=y1*y1; y22=y2*y2;
    
    m[0]+=a;
    m[1]+=(x1+x2)*a;
    m[2]+=(y1+y2)*a;    
    m[3]+=(x1*(x1+x2)+x22)*a;
    m[4]+=(y1*(x1+x1+x2)+y2*(x1+x2+x2))*a;
    m[5]+=(y1*(y1+y2)+y22)*a;
    m[6]+=(x12+x22)*(x1+x2)*a;
    m[7]+=(y1*(x1*(3*x1+2*x2)+x22)+y2*(x12+x2*(2*x1+3*x2)))*a;
    m[8]+=(x1*(y1*(3*y1+2*y2)+y22)+x2*(y12+y2*(2*y1+3*y2)))*a;
    m[9]+=(y12+y22)*(y1+y2)*a;
    m[10]+=((((x1+x2)*x1+x22)*x1+x22*x2)*x1+x22*x22)*a;
    m[11]+=((y1+4*y2)*x22*x2+(2*y1+3*y2)*x1*x22+\
        (3*y1+2*y2)*x12*x2+(4*y1+y2)*x12*x1)*a;
    m[12]+=((2*y12+6*y1*y2+12*y22)*x22+(6*(y12+y22)+8*y1*y2)*x1*x2+\
        (12*y12+6*y1*y2+2*y22)*x12)*a;
    m[13]+=((x1+4*x2)*y22*y2+(2*x1+3*x2)*y1*y22+\
        (3*x1+2*x2)*y12*y2+(4*x1+x2)*y12*y1)*a;
    m[14]+=((((y1+y2)*y1+y22)*y1+y22*y2)*y1+y22*y22)*a;
    i++;j++;
  }while(i!=a3);
  m[1]/=3; m[2]/=3;
  m[3]/=6; m[4]/=12; m[5]/=6;
  m[6]/=10; m[7]/=30; m[8]/=30; m[9]/=10;
  m[10]/=15; m[11]/=60; m[12]/=180;
  m[13]/=60; m[14]/=15;
  centre[0]=m[1]/m[0]; centre[1]=m[2]/m[0];
  return(OK);
}
#undef FNAME

/*******************************************************************/
int FitSuper_Ell_Moments(double moments_object[],double *sup_c1,double *sup_f1,double sup_tr1[][3],double *sup_c2,double *sup_f2,double sup_tr2[][3])
{
int i,j;
double gamma_fkt( double x);
double c,s,h,s1_opt,s2_opt;
double m00,m40,m22,m04;
double qm00,qm20,qm40,qm22,qm04;
double m1s00,m1s40,m1s22,m1s04;
double m2s00,m2s40,m2s22,m2s04;
double f1c,f2c,f1c_min,f2c_min,c1_min,c2_min;
double maf[21];
double atr[3][3];
f1c_min=f2c_min=1.0e+60;
c1_min=c2_min=0.05;
AffinIterateMoments(moments_object,maf,atr);
m00=maf[0];
m40=maf[10];
m22=maf[12];
m04=maf[14];
s1_opt=s2_opt=0.01;
/* On maf[ ] are after the iteration method normalized the moments of the given object */
for (c=0.01;c<100;c+=0.01)
  {
    qm00=4.0/(c*2.0)*gamma_fkt(1.0/c)*gamma_fkt(1.0/c)/gamma_fkt(2.0/c);
    qm20=4.0/(c*4.0)*gamma_fkt(3.0/c)*gamma_fkt(1.0/c)/gamma_fkt(4.0/c);
    qm40=4.0/(c*6.0)*gamma_fkt(5.0/c)*gamma_fkt(1.0/c)/gamma_fkt(6.0/c);
    qm22=4.0/(c*6.0)*gamma_fkt(3.0/c)*gamma_fkt(3.0/c)/gamma_fkt(6.0/c);
    qm04=4.0/(c*6.0)*gamma_fkt(1.0/c)*gamma_fkt(5.0/c)/gamma_fkt(6.0/c);
     s=1.0/pow(qm20,0.25);
    m1s00=s*s*qm00;
    m1s40=s*s*s*s*s*s*qm40;
    m1s22=s*s*s*s*s*s*qm22;
    m1s04=s*s*s*s*s*s*qm04;
    h=(m00-m1s00);
    f1c=h*h;
    h=(m40-m1s40);
    f1c+=h*h;
    h=(m22-m1s22);
    f1c+=h*h;
    h=(m40-m1s40);
    f1c+=h*h;
    if (f1c<f1c_min) {f1c_min=f1c;c1_min=c;s1_opt=s;}
    
    /* m1spq are normalized moments, now still rotation around 45 degrees of m2spq */
    m2s00=m1s00;
    m2s40=(m1s40+6*m1s22+m1s04)/4.0;
    m2s22=(m1s40-2*m1s22+m1s04)/4.0;
    m2s04=m2s40; 
    h=(m00-m2s00);
    f2c=h*h;
    h=(m40-m2s40);
    f2c+=h*h;
    h=(m22-m2s22);
    f2c+=h*h;
    h=(m40-m2s40);
    f2c+=h*h;
    if (f2c<f2c_min) {f2c_min=f2c;c2_min=c;s2_opt=s;}
}
*sup_c1=c1_min;
*sup_f1=f1c_min;
*sup_c2=c2_min;
*sup_f2=f2c_min;
for (i=0;i<3;++i)
for (j=0;j<3;++j)
  {sup_tr1[i][j]=atr[i][j]; }
ScaleTrans(0.0,0.0,1.0/s1_opt,1.0/s1_opt,sup_tr1);
RotTrans(0.0,0.0,M_PI/4.0,atr);
for (i=0;i<3;++i)
for (j=0;j<3;++j)
  {sup_tr2[i][j]=atr[i][j]; }
ScaleTrans(0.0,0.0,1.0/s2_opt,1.0/s2_opt,sup_tr2);
return(0);
}
/*******************************************************************/
/* Berechnung der Gamma-Funktion mit                               */
/* a) Reihenentwicklung                                            */
/* b) Horner-Schema                                                */
/* c) Rekursionsformel                                             */
/*******************************************************************/
double gamma_fkt( double x)
{
  double h,h1,h2;
  double z;
  double c0,c1,c2,c3,c4,c5,c6;
  if (x<=0) {printf("\n Error in the calculation of the gamma function\n"); exit(1);}
  c0=1.0;
  c1=1.0/12.0;
  c2=1.0/288.0;
  c3=-139.0/51840.0;
  c4=-134.0/583949.0;
  c5=101.0/128820.0;
  c6=1.0/14514.0;
  /* Use of G(1+x+6)=(x+6)*(x+5)*(x+4)*(x+3)*(x+2)*(x+1)*x*G(x) */
  z=x+6.0; 
  h1=pow(z/M_E,z);
  h1*=sqrt(2*M_PI*z);
  /* Horner-Schema */
  h2=c6/z+c5;
  h2=h2/z+c4;
  h2=h2/z+c3;
  h2=h2/z+c2;
  h2=h2/z+c1;
  h2=h2/z+c0;
  h=(h1*h2)/(x*(x+1.0)*(x+2.0)*(x+3.0)*(x+4.0)*(x+5.0)*(x+6.0));
  return(h);
  }
/**************************************************************/
/* Affine moments (iteration method Voss) */
/**************************************************************/
#define FNAME "AffinIterateMoments"
int AffinIterateMoments(double m[15],double maf[15],double atr[3][3])
{
  int i;
  double alpha,beta,gamma,delta,h1,h2,xs,ys;
/* Modification of the direction of rotation, if m00 < 0 */
  alpha=0;beta=0;xs=0;ys=0;
  if(m[0]<0)
    for(i=0;i<15;i++) maf[i]=-m[i];
  else
    for(i=0;i<15;i++) maf[i]=m[i];
  /*Translationsnormierung*/
  InitTrans(atr);
  if(maf[1]!=0 || maf[2]!=0)
  {
    xs=-maf[1]/maf[0]; ys=-maf[2]/maf[0];
    TranslateMoments(maf,xs,ys,maf);
  }
/*evtl. instabile Standardlage bei achsenparallelen Figuren, desh. vorher Drehung um bel. Winkel*/
  if(fabs(maf[13])<1e-10 || fabs(maf[11])<1e-10)
  {
    h1=cos(0.5); h2=sin(0.5);
    RotateMoments(maf,h1,h2,maf);
    atr[0][0]=h1; atr[0][1]=-h2;
    atr[1][0]=h2; atr[1][1]=h1;    
  }
/* iterative Scherungsnormierung maf13=maf31=0*/
  if((maf[14]<1e-20)&&(maf[10]<1e-20))
  {
    Message(FNAME,M_NUM_INSTABILITY,ERROR);
    return(ERROR);
  }
  i=0;
  do
  {
    if(fabs(maf[14])>1e-20)
    {
      alpha=-maf[13]/maf[14];
      XShearMoments(maf,alpha,maf);
      atr[0][0]+=alpha*atr[1][0];
      atr[0][1]+=alpha*atr[1][1];
    }
    if(fabs(maf[10])>1e-20)
    {
      beta=-maf[11]/maf[10];
      YShearMoments(maf,beta,maf);
      atr[1][0]+=beta*atr[0][0];
      atr[1][1]+=beta*atr[0][1];
    }
    i++;
  }while(((fabs(alpha)>1e-5)||(fabs(beta)>1e-5))&&(i<200));
/* anisotrope Skalierungsnormierung */
  if((fabs(maf[3])<1e-7)||(fabs(maf[5])<1e-7))
  {
    Message(FNAME,M_NUM_INSTABILITY,ERROR);
    return(ERROR);
  }
  h1=maf[5]/(maf[3]*maf[3]*maf[3]);
  h2=maf[3]/(maf[5]*maf[5]*maf[5]);
  if(h1<0||h2<0)
  {
    Message(FNAME,M_NUM_INSTABILITY,ERROR);
    return(ERROR);
  }
  gamma=pow(h1,0.125);
  delta=pow(h2,0.125);
  ScaleMoments(maf,gamma,delta,maf);
  atr[0][0]*=gamma;
  atr[0][1]*=gamma;
  atr[1][0]*=delta;
  atr[1][1]*=delta;
/*   NormMoments(maf,atr); */
  atr[0][2]=atr[0][0]*xs+atr[0][1]*ys;
  atr[1][2]=atr[1][0]*xs+atr[1][1]*ys;
  return(i);
}
#undef FNAME

/*************************************************************************/
void Draw_Sup_Ellipse(int nx, int ny,double c,double tr[][3])
{
int i,j;
double xs,ys,xs1,ys1;
double z1u,z1o;
for (j=0;j<ny;++j)
for (i=0;i<nx;++i)
  {
    xs=tr[0][0]*(double)i+tr[0][1]*(double)j+tr[0][2];
    ys=tr[1][0]*(double)i+tr[1][1]*(double)j+tr[1][2];
    xs1=tr[0][0]*((double)i+1.0)+tr[0][1]*(double)j+tr[0][2];
    ys1=tr[1][0]*((double)i+1.0)+tr[1][1]*(double)j+tr[1][2];
    xs=fabs(xs);
    ys=fabs(ys);
    xs1=fabs(xs1);
    ys1=fabs(ys1);
    z1u=pow(xs,c)+pow(ys,c)-1.0;
    z1o=pow(xs1,c)+pow(ys1,c)-1.0;
    if (z1u<=0 && z1o >= 0) printf("%d %d\n",i,j);
    if (z1u>=0 && z1o <=0 ) printf("%d %d\n",i,j);
  }

for (i=0;i<nx;++i)
for (j=0;j<ny;++j)
  {
    xs=tr[0][0]*(double)i+tr[0][1]*(double)j+tr[0][2];
    ys=tr[1][0]*(double)i+tr[1][1]*(double)j+tr[1][2];
    xs1=tr[0][0]*(double)i+tr[0][1]*((double)j+1.0)+tr[0][2];
    ys1=tr[1][0]*(double)i+tr[1][1]*((double)j+1.0)+tr[1][2];
    xs=fabs(xs);
    ys=fabs(ys);
    xs1=fabs(xs1);
    ys1=fabs(ys1);
    z1u=pow(xs,c)+pow(ys,c)-1.0;
    z1o=pow(xs1,c)+pow(ys1,c)-1.0;
    if (z1u<=0 && z1o >= 0) printf("%d %d\n",i,j);
    if (z1u>=0 && z1o <=0 ) printf("%d %d\n",i,j);
  }
}

/*** FreePointList() *******************************/
#define FNAME "FreePointList"
int FreePointList(PointList pl)
{
  if (pl==NULL)
  {
    Message(FNAME,M_WRONG_PTR, WRONG_PARAM);
    return WRONG_PARAM;
  }

  if (pl->lng != 0)
  {
    free(pl->xptr);
    free(pl->yptr);
    free(pl->wptr);
  }
  free(pl);
  SetOk();
  return(OK);
}
#undef FNAME

#define FNAME "FeatureQuadrFunc"
/* #define DEBUG */ 
int FeatureQuadrFunc(double par[6],double feat[5],int *type)
{
  double A,A33;
  double Ah,Bh;
  double a,b,c,d,e,f;
  double eps=1e-30;
  double xm,ym,phi,ha,hb,hc,hd,hf,hha,hhb,hhc,hhd,hhe;
  a=par[0];b=par[1];c=par[2];d=par[3];e=par[4];f=par[5];
  A=a*(b*f-d*d/4)-
    e/2*(e*f/2-d*c/4)+
    c/2*(e*d/4-b*c/2);
  #ifdef DEBUG
    printf("A %f\n",A);
    getchar();
  #endif
  if (fabs(A)<eps)
  {
    *type=DEGENERATE;
    return OK;
  }
  A33=a*b-e*e/4;
  #ifdef DEBUG
    printf("A33 %f\n",A33);
    getchar();
  #endif
  if (fabs(A33)<eps) /* Parabel*/
  {
    if (a<0) {a=-a;b=-b;c=-c;d=-d;e=-e;f=-f;}
    /* Bestimmung des Scheitelpunktes */
    #ifdef DEBUG
      printf("Beginn Parabel\n");
      getchar();
      printf("a %f b %f c %f d %f e %f f %f\n",a,b,c,d,e,f);
      printf("fabs(2*b*sqrt(b)-e*sqrt(a) %g\n",fabs(2*b*sqrt(b)-e*sqrt(a)));
    #endif
    if ( fabs(2*b*sqrt(b)-e*sqrt(a)) < eps) /* Parabelachse parallel zur y-Achse */
    {
      #ifdef DEBUG
    printf("parallel y-Achse\n");
    getchar();
      #endif
      xm = (sqrt(a)*c - sqrt(b)*d) / (sqrt(b)*e - 2*sqrt(a)*a);
      if (fabs(b)<eps) ym=-(a*xm*xm + c*xm + f) / (d+e*xm);
      else ym=(e*xm+d)/(-2*b);
      goto label;
    }
    if ( fabs(sqrt(b)*e-2*a*sqrt(a)) < eps) /* Parabelachse parallel zur x-Achse */
    {
      #ifdef DEBUG
    printf("parallel x-Achse\n");
    printf("a %f b %f c %f d %f e %f f %f\n",a,b,c,d,e,f);
    getchar();
      #endif
      ym = (sqrt(a)*c - sqrt(b)*d) / (2*sqrt(b)*b - sqrt(a)*e);
      if (fabs(a)<eps) xm=-(b*ym*ym + d*ym + f) / (c+e*ym);
      else xm=(c+e*ym)/(-2*a);
      goto label;
    }
    /* Allgemeiner Fall */
    #ifdef DEBUG
      printf("Allgemeiner Fall\n");
      getchar();
    #endif
    Ah=(2*a*sqrt(a)-sqrt(b)*e) / (2*b*sqrt(b)-e*sqrt(a));
    Bh=(sqrt(a)*c-sqrt(b)*d) / (2*b*sqrt(b)-e*sqrt(a));
    if (fabs(a+b*Ah*Ah+e*Ah)<eps)
      xm=(-b*Bh*Bh-d*Bh-f)/(2*b*Ah*Bh+c+d*Ah+e*Bh);
    else
      xm=(2*b*Ah*Bh+c+d*Ah+e*Bh) / (-2*(a+b*Ah*Ah+e*Ah));
    ym=Ah*xm+Bh;
    label:
    feat[0]=xm;feat[1]=ym;
    #ifdef DEBUG
      printf("xm %f ym %f\n",xm,ym);
      getchar();
    #endif

    /* Bestimmung der neuen Koeffizienten nach der Verschiebung */
    hc=2*a*xm+c+e*ym;
    hd=2*b*ym+d+e*xm;
    phi=atan2(hd,hc);
    #ifdef DEBUG
      printf("phi %f\n",phi);
      getchar();
    #endif

    hf=a*xm*xm+b*ym*ym+c*xm+d*ym+e*xm*ym+f;
      par[0]=hha= a*Sqr(cos(phi))+b*Sqr(sin(phi))+e*sin(phi)*cos(phi);
      par[1]=hhb= a*Sqr(sin(phi))+b*Sqr(cos(phi))-e*sin(phi)*cos(phi);
      par[2]=hhc= hc*cos(phi)+hd*sin(phi);
      par[3]=hhd=-hc*sin(phi)+hd*cos(phi);
      par[4]=hhe= 2*( a*sin(phi)*cos(phi) - b*sin(phi)*cos(phi) )+
          e*( Sqr(cos(phi))-Sqr(sin(phi)) );
      par[5]=hf;
    #ifdef DEBUG
      PrintVecRn("par nach RÅcknahme Verdrehung",par,6);
      getchar();
    #endif
    feat[2]=-hhc/(hhb*2);
    #ifdef DEBUG
      printf("p: %f\n",feat[2]);
      getchar();
    #endif
    if (feat[2]<0)
    {
      feat[2]=-feat[2];
      phi+=M_PI;
    }
    feat[3]=fmod(phi,2*M_PI);
    if (phi<0) feat[3]+=2*M_PI;
    *type=PARABEL;
    return OK;
  }
  else /* Ellipse/Hyperbel */
  {
    feat[0]=xm=(e*d/4 - c*b/2)/A33;
    feat[1]=ym=(c*e/4 - a*d/2)/A33;

    #ifdef DEBUG
      printf("ELLIPSE\n");
      printf("xm %f ym %f\n",xm,ym);
      getchar();
    #endif

    hf=a*xm*xm + b*ym*ym + c*xm + d*ym + e*xm*ym + f;
    if (fabs(e)<eps && fabs(a-b)<eps) phi=0;
    else phi=atan2(e,a-b)/2;
    #ifdef DEBUG
      printf("phi %f cos(phi) %f\n",phi,cos(phi));
      PrintVecRn("par",par,6);
      getchar();
    #endif
    ha=(a*Sqr(cos(phi))+b*Sqr(sin(phi))+e*sin(phi)*cos(phi))/(-hf);
    hb=(a*Sqr(sin(phi))+b*Sqr(cos(phi))-e*sin(phi)*cos(phi))/(-hf);
    #ifdef DEBUG
      printf("ha %f hb %f\n",ha,hb);
      getchar();
    #endif
    if (SignD(ha*hb)<0) /* Hyperbel */
    {
      #ifdef DEBUG
    printf("Hyperbel\n");
    getchar();
      #endif
      *type=HYPERBEL;
      if (ha<0)
      {
    phi+=M_PI/2;
    feat[2]=1/sqrt(hb);
    feat[3]=1/sqrt(-ha);
      }
      else
      {
    feat[2]=1/sqrt(ha);
    feat[3]=1/sqrt(-hb);
    if (phi<0) feat[2]+=M_PI;
      }
      feat[4]=phi;
      if (phi>M_PI)
    feat[4]-=M_PI;
      if (phi<0)
    feat[4]+=M_PI;
      return OK;
    }
    else /* Ellipse */
    {
      if (ha<0)
      {
    *type=DEGENERATE;
    return OK;
      }
      else
      {
    #ifdef DEBUG
      printf("Ellipse\n");
      getchar();
    #endif
    *type=ELLIPSE;
    if (ha>hb) /* Vertauschen der Achsen */
    {
      phi+=M_PI/2;
      feat[2]=1/sqrt(hb);
      feat[3]=1/sqrt(ha);
    }
    else
    {
      feat[2]=1/sqrt(ha);
      feat[3]=1/sqrt(hb);
    }
    #ifdef DEBUG
      printf("vor Auswertung phi: %f\n",phi);
      getchar();
    #endif
    feat[4]=phi;
    if (phi>M_PI)
      feat[4]-=M_PI;
    if (phi<0)
      feat[4]+=M_PI;
    return OK;
      }
    }
  }
}
#undef FNAME
/* ============================================================================ */

int no_pixels;
#define MAX_PIXELS 10000
int x[MAX_PIXELS],y[MAX_PIXELS];

main(argc,argv)
int argc;
char *argv[];
{
    FILE *fp_in,*fp_out;
    char *infile,*outfile,file_type[50];
    int i,j;
    int endoffile;
    int circle,ellipse,rectangle,superellipse,triangle;
    double ellipticity,rectangularity,triangularity;
    double sup_c1,sup_c2,sup_f1,sup_f2,sup_tr1[3][3],sup_tr2[3][3];
    double xc,yc,radius,ell_par[5];

    double mm[15];
    double pc0[2],t[3][2],r[4][2];
    PointList plm;

    circle = ellipse = rectangle = superellipse = triangle = FALSE;
    infile = outfile = NULL;

    for (i = 1; i < argc; i++) {
        if (argv[i][0] == '-') {
            switch(argv[i][1]) {
                case 'i':
                    i++;
                    infile = argv[i];
                    break;
                case 'o':
                    i++;
                    outfile = argv[i];
                    break;
                case 'c':
                    circle = TRUE;
                    break;
                case 'e':
                    ellipse = TRUE;
                    break;
                case 'r':
                    rectangle = TRUE;
                    break;
                case 's':
                    superellipse = TRUE;
                    break;
                case 't':
                    triangle = TRUE;
                    break;
                default:
                    fprintf(stderr,"ERROR: unknown option %s\n",argv[i]);
                    usage(argv[0]);
                    exit(-1);
            }
        }
        else {
            fprintf(stderr,"ERROR: unknown option %s\n",argv[i]);
            usage(argv[0]);
            exit(-1);
        }
    }

    if ((infile == NULL) || (outfile == NULL)) {
        fprintf(stderr,"need input & output files\n");
        usage(argv[0]);
        exit(-1);
    }

    if ((fp_in=fopen(infile,"r")) == NULL){
        printf("cant open %s\n",infile);
        exit(-1);
    }
    /* read magic word for format of file */
    fscanf(fp_in,"%s\n",file_type);
    j = strcmp(file_type,"pixel");
    if (j != 0){
        fprintf(stderr,"not link data file - aborting\n");
        exit(-1);
    }

    if ((fp_out = fopen(outfile,"w")) == NULL) {
        fprintf(stderr,"cant open %s\n",outfile);
        exit(-1);
    }

    if (circle || ellipse)
        fprintf(fp_out,"super_float\n");
    else
        fprintf(fp_out,"pixel_float\n");

    do {
        read_link_data(fp_in,&endoffile);
    
        plm = NewPointList(no_pixels);
        for (i=0; i < no_pixels; i++) {
            plm->xptr[i] = x[i];
            plm->yptr[i] = y[i];
        }
        MomentPolygon(plm,mm,pc0);
        if (mm[0]<0)
            for (i=0;i<15;++i)
                mm[i] = -mm[i];

        FreePointList(plm);
      
        fprintf(fp_out,"list: 0\n");
        if (circle) {
            FitCircleMoments(mm,&xc,&yc,&radius);
            fprintf(fp_out,"circle: 0 %f %f %f\n",xc,yc,radius);
        }
        else if (ellipse) {
            FitEllipseMoments(mm,ell_par,&ellipticity);
            fprintf(fp_out,"ellipse: 0 %f %f %f %f %f %f %f %f %f %d\n",
                ell_par[0],ell_par[1],
                0.0,0.0,0.0,0.0,
                ell_par[2],ell_par[3],ell_par[4],1);
            ellipticity = 1.0/(1+ellipticity);
            printf("ellipticity = %f\n",ellipticity);
        }
        else if (triangle) {
            FitTriangle(mm,t,&triangularity);
            for (i=0; i < 3; i++)
                fprintf(fp_out,"%f %f\n",t[i][0],t[i][1]);
            fprintf(fp_out,"%f %f\n",t[0][0],t[0][1]);
            if (endoffile)
                fprintf(fp_out,"-1 -1\n");
            else
                fprintf(fp_out,"-1 0\n");
            triangularity = 1.0/(1+triangularity);
            printf("triangularity = %f\n",triangularity);
        }
        else if (rectangle) {
            FitRectangle(mm,r,&rectangularity);
            for (i=0; i < 4; i++)
                fprintf(fp_out,"%f %f\n",r[i][0],r[i][1]);
            fprintf(fp_out,"%f %f\n",r[0][0],r[0][1]);
            if (endoffile)
                fprintf(fp_out,"-1 -1\n");
            else
                fprintf(fp_out,"-1 0\n");
            /*
            printf("rectangularity = %f\n",rectangularity);
            */
            rectangularity = 1.0/(1+rectangularity);
            printf("rectangularity = %f\n",rectangularity);
        }
        /* actually fits a supercircle with affine deformation */
        else if (superellipse) {
            FitSuper_Ell_Moments(mm,&sup_c1,&sup_f1,sup_tr1,&sup_c2,&sup_f2,sup_tr2);
            printf("\nSuperquadrik c1 f1 %f %f \n",sup_c1,sup_f1);
            printf("Superquadrik c2 f2 %f %f \n",sup_c2,sup_f2);
            /*
            Draw_Sup_Ellipse(512,512,sup_c2,sup_tr2);
            */
            Draw_Sup_Ellipse(512,512,sup_c1,sup_tr1);
            if (endoffile)
                fprintf(fp_out,"-1 -1\n");
            else
                fprintf(fp_out,"-1 0\n");
        }
    }  while (!endoffile);
    if (circle || ellipse)
        fprintf(fp_out,"endl:\nendf:\n");

    fclose(fp_out);
}

/* ----------------------------- my I/O stuff ----------------------------- */

/* store in [0..no_pixels-1] */
read_link_data(fp,endoffile)
FILE *fp;
int *endoffile;
{
    char dumstring[50];
    int j;

    fscanf(fp,"%s %d\n",dumstring,&j);
    printf("list: %d ",j);
    j = -1;
    do{
       j++;
       fscanf(fp,"%d %d\n",&x[j],&y[j]);
    } while(x[j] != -1);
    *endoffile = (y[j] == -1);
    no_pixels = j;
    if ((x[no_pixels-1] == x[0]) && (y[no_pixels-1] == y[0]))
        no_pixels--;
    if (no_pixels >= MAX_PIXELS) {
        fprintf(stderr,"ERROR: Too many pixels\n");
        exit(-1);
    }
}

usage(progname)
char *progname;
{
    printf("usage: %s -i infile -o outfile\n",progname);
    printf("options:\n");
    printf("         -c  fit circle\n");
    printf("         -e  fit ellipse\n");
    printf("         -r  fit rectangle\n");
    printf("         -s  fit superellipse\n");
    printf("         -t  fit triangle\n");
    exit(-1);
}

printm(i)
{
switch(i) {
case 0: printf("00 ");break;
case 1: printf("10 ");break;
case 2: printf("01 ");break;
case 3: printf("20 ");break;
case 4: printf("11 ");break;
case 5: printf("02 ");break;
case 6: printf("30 ");break;
case 7: printf("21 ");break;
case 8: printf("12 ");break;
case 9: printf("03 ");break;
case 10: printf("40 ");break;
case 11: printf("31 ");break;
}
}