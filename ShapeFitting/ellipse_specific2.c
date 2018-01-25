/* Bastardized C version of Java ellipse fitting code by Pilu
 * reads in pixel lists and fits ellipses to them, and outputs superdata
 *
 * Paul Rosin
 * March 1998
 *
 * extracted from stand-alone program
 *
 * Paul Rosin
 * July 2005
 */

#include <stdio.h>
#include <math.h>

#ifndef FALSE
# define FALSE 0
# define TRUE (!FALSE)
#endif

#define ABS(x)       (((x)<0.0)? (-(x)): (x))
#define SQR(x)       ((x)*(x))

#define MAX_POINTS 10000

extern float x_trans3[],y_trans3[];
extern int nseg2;

double zero=10e-20;
int solind=0;

void ROTATE(double a[7][7], int i, int j, int k, int l,
                    double tau, double s)
{
    double g,h;
    g=a[i][j]; h=a[k][l]; a[i][j]=g-s*(h+g*tau);
    a[k][l]=h+s*(g-h*tau);
}

void jacobi(double a[7][7], int n, double d[7] , double v[7][7], int nrot)
{
    int j,iq,ip,i;
    double tresh,theta,tau,t,sm,s,h,g,c;

    double b[7];
    double z[7];

    for (ip=1;ip<=n;ip++) {
      for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
      v[ip][ip]=1.0;
    }
    for (ip=1;ip<=n;ip++) {
      b[ip]=d[ip]=a[ip][ip];
      z[ip]=0.0;
    }
    nrot=0;
    for (i=1;i<=50;i++) {
      sm=0.0;
      for (ip=1;ip<=n-1;ip++) {
        for (iq=ip+1;iq<=n;iq++)
          sm += ABS(a[ip][iq]);
      }
      if (sm == 0.0) {
        /*    free_vector(z,1,n);
          free_vector(b,1,n);  */
        return;
      }
      if (i < 4)
        tresh=0.2*sm/(n*n);
      else
        tresh=0.0;
      for (ip=1;ip<=n-1;ip++) {
        for (iq=ip+1;iq<=n;iq++) {
          g=100.0*ABS(a[ip][iq]);
          if (i > 4 && (ABS(d[ip])+g == ABS(d[ip]))
          && (ABS(d[iq])+g == ABS(d[iq])))
        a[ip][iq]=0.0;
          else if (ABS(a[ip][iq]) > tresh) {
        h=d[iq]-d[ip];
        if (ABS(h)+g == ABS(h))
          t=(a[ip][iq])/h;
        else {
          theta=0.5*h/(a[ip][iq]);
          t=1.0/(ABS(theta)+sqrt(1.0+theta*theta));
          if (theta < 0.0) t = -t;
        }
        c=1.0/sqrt(1+t*t);
        s=t*c;
        tau=s/(1.0+c);
        h=t*a[ip][iq];
        z[ip] -= h;
        z[iq] += h;
        d[ip] -= h;
        d[iq] += h;
        a[ip][iq]=0.0;
        for (j=1;j<=ip-1;j++)
            ROTATE(a,j,ip,j,iq,tau,s);
        for (j=ip+1;j<=iq-1;j++)
            ROTATE(a,ip,j,j,iq,tau,s);
        for (j=iq+1;j<=n;j++)
            ROTATE(a,ip,j,iq,j,tau,s);
        for (j=1;j<=n;j++)
            ROTATE(v,j,ip,j,iq,tau,s);
        ++nrot;
          }
        }
      }
      for (ip=1;ip<=n;ip++) {
        b[ip] += z[ip];
        d[ip]=b[ip];
        z[ip]=0.0;
      }
    }
}


/*  Perform the Cholesky decomposition     */
/* Return the lower triangular L  such that L*L'=A   */
void choldc(double a[7][7], int n, double l[7][7])
{
    int i,j,k;
    double sum;
    double p[7];

    for (i=1; i<=n; i++)  {
      for (j=i; j<=n; j++)  {
        for (sum=a[i][j],k=i-1;k>=1;k--) sum -= a[i][k]*a[j][k];
        if (i == j) {
          if (sum<=0.0)
        {}
          else
        p[i]=sqrt(sum); }
        else
          {
        a[j][i]=sum/p[i];
          }
      }
    }
    for (i=1; i<=n; i++)
      for (j=i; j<=n; j++)
        if (i==j)
          l[i][i] = p[i];
        else
          {
        l[j][i]=a[j][i];
        l[i][j]=0.0;
          }
}

/********************************************************************/
/**    Calcola la inversa della matrice  B mettendo il risultato   **/
/**    in InvB . Il metodo usato per l'inversione e' quello di     **/
/**    Gauss-Jordan.   N e' l'ordine della matrice .               **/
/**    ritorna 0 se l'inversione  corretta altrimenti ritorna     **/
/**    SINGULAR .                                                  **/
/********************************************************************/
int inverse6(double TB[7][7], double InvB[7][7], int N)
{
  int k,i,j,p,q;
  double mult;
  double D,temp;
  double maxpivot;
  int npivot;
  double B[7][8];
  double A[7][14];
  double eps = 10e-20;

  for (k=1;k<=N;k++)
    for (j=1;j<=N;j++)
        B[k][j]=TB[k][j];

  for (k=1;k<=N;k++)
{
  for (j=1;j<=N+1;j++)
    A[k][j]=B[k][j];
  for (j=N+2;j<=2*N+1;j++)
    A[k][j]=(float)0;
  A[k][k-1+N+2]=(float)1;
}
  for (k=1;k<=N;k++)
{
  maxpivot=ABS((double)A[k][k]);
  npivot=k;
  for (i=k;i<=N;i++)
    if (maxpivot<ABS((double)A[i][k]))
      {
    maxpivot=ABS((double)A[i][k]);
    npivot=i;
      }
  if (maxpivot>=eps)
    {      if (npivot!=k)
         for (j=k;j<=2*N+1;j++)
           {
         temp=A[npivot][j];
         A[npivot][j]=A[k][j];
         A[k][j]=temp;
           } ;
       D=A[k][k];
       for (j=2*N+1;j>=k;j--)
         A[k][j]=A[k][j]/D;
       for (i=1;i<=N;i++)
         {
           if (i!=k)
         {
           mult=A[i][k];
           for (j=2*N+1;j>=k;j--)
             A[i][j]=A[i][j]-mult*A[k][j] ;
         }
         }
     }
  else
    {  /* printf("\n The matrix may be singular !!") ; */
       return(-1);
     };
}
  /**   Copia il risultato nella matrice InvB  ***/
  for (k=1,p=1;k<=N;k++,p++)
    for (j=N+2,q=1;j<=2*N+1;j++,q++)
          InvB[p][q]=A[k][j];
  return(0);
}            /*  End of INVERSE   */

void AperB6(double _A[7][7], double _B[7][7], double _res[7][7],
        int _righA, int _colA, int _righB, int _colB)
{
    int p,q,l;

    for (p=1;p<=_righA;p++)
    for (q=1;q<=_colB;q++)
      { _res[p][q]=0.0;
        for (l=1;l<=_colA;l++)
          _res[p][q]=_res[p][q]+_A[p][l]*_B[l][q];
      }
}

void A_TperB(double _A[7][7], double  _B[7][7], double _res[7][7],
             int _righA, int _colA, int _righB, int _colB)
{
    int p,q,l;
    for (p=1;p<=_colA;p++)
    for (q=1;q<=_colB;q++) {
        _res[p][q]=0.0;
        for (l=1;l<=_righA;l++)
          _res[p][q]=_res[p][q]+_A[l][p]*_B[l][q];
    }
}

void AperB_T(double _A[7][7], double _B[7][7], double _res[7][7],
             int _righA, int _colA, int _righB, int _colB)
{
    int p,q,l;
    for (p=1;p<=_colA;p++)
    for (q=1;q<=_colB;q++)
      { _res[p][q]=0.0;
        for (l=1;l<=_righA;l++)
          _res[p][q]=_res[p][q]+_A[p][l]*_B[q][l];
      }
}

ellipse_fit(x_cent,y_cent,major_axis,minor_axis,rot_angle)
double *x_cent,*y_cent,*major_axis,*minor_axis,*rot_angle;
{
    int np;           /* number of points */
    double D[MAX_POINTS+1][7];
    double S[7][7];
    double Const[7][7] ;
    double temp[7][7];
    double L[7][7];
    double C[7][7];

    double invL[7][7];
    double d[7];
    double V[7][7];
    double sol[7][7];
    double tx,ty;
    int nrot=0;
    int i,j;
    double pvec[7];

    double x_centT,y_centT,major_axisT,minor_axisT,rot_angleT;

    np = nseg2;

    if (np<6) {
        fprintf(stderr,"ERROR: too few points\n");
        return;
    }

    /* ### hack for compatibility ### */
    x_trans3[0] = x_trans3[nseg2];
    y_trans3[0] = y_trans3[nseg2];
    for (i=0;i<=6;i++) {
       d[i] = 0.0;
       pvec[i] = 0.0;

       for (j=0;j<=6;j++) {
          S[i][j] = 0.0;
          Const[i][j] = 0.0;
          temp[i][j] = 0.0;
          L[i][j] = 0.0;
          C[i][j] = 0.0;
          invL[i][j] = 0.0;
          V[i][j] = 0.0;
          sol[i][j] = 0.0;
       }

       for (j=0;j<=MAX_POINTS;j++)
          D[j][i] = 0.0;
    }

    Const[1][3] = -2;
    Const[2][2] = 1;
    Const[3][1] = -2;

    /* Now first fill design matrix */
    for (i=1; i <= np; i++) {
        tx = x_trans3[i-1];
        ty = y_trans3[i-1];
        D[i][1] = tx*tx;
        D[i][2] = tx*ty;
        D[i][3] = ty*ty;
        D[i][4] = tx;
        D[i][5] = ty;
        D[i][6] = 1.0;
    }

    /*pm(Const,"Constraint"); */
    /* Now compute scatter matrix  S */
    /*
    A_TperB(D,D,S,np,6,np,6);
    */

    {
    int p,q,l;
    for (p=1;p<=6;p++)
    for (q=1;q<=6;q++) {
        S[p][q]=0.0;
        for (l=1;l<=np;l++)
          S[p][q]=S[p][q]+D[l][p]*D[l][q];
    }
    }

    choldc(S,6,L);

    inverse6(L,invL,6);

    AperB_T(Const,invL,temp,6,6,6,6);
    AperB6(invL,temp,C,6,6,6,6);

    jacobi(C,6,d,V,nrot);

    A_TperB(invL,V,sol,6,6,6,6);

    /* Now normalize them  */
    for (j=1;j<=6;j++)  /* Scan columns */
      {
        double mod = 0.0;
        for (i=1;i<=6;i++)
          mod += sol[i][j]*sol[i][j];
        for (i=1;i<=6;i++)
          sol[i][j] /=  sqrt(mod);
      }

    for (i=1; i<=6; i++)
        if (d[i]<0 && ABS(d[i])>zero)
            solind = i;

    /* Now fetch the right solution */
    for (j=1;j<=6;j++)
      pvec[j] = sol[j][solind];

    determine_parameters(pvec[1],pvec[2],pvec[3],pvec[4],pvec[5],pvec[6],
                         &x_centT,&y_centT,&major_axisT,&minor_axisT,&rot_angleT);

printf("%f %f      %f %f    %f\n",x_centT,y_centT,major_axisT,minor_axisT,rot_angleT);
    *x_cent = x_centT;
    *y_cent = y_centT;
    *major_axis = major_axisT;
    *minor_axis = minor_axisT;
    *rot_angle = rot_angleT;
}

