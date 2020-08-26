#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

int main( int argc, char* argv[] )
{

    double data[19][2] = {{46,40},{49,50},{51,55},{52,63},{54,72},{56,70},{57,77},
        {58,73},{59,90},{60,93},{61,96},{62,88},{63,99},{64,110},
        {66,113},{67,120},{68,127},{71,137},{72,132}};
    double lslope = 0.; double rslope = -100.;

    // INIT
    int n = 19;
    
    double f[19];
    double t[19];
    double a[19],b[19],c[19],d[19];
    double alpha[19],mu[19],z[19],L[19],h[19];
    
    // memset
    memset(f,0,19*sizeof(double));
    memset(t,0,19*sizeof(double));
    memset(a,0,19*sizeof(double));
    memset(b,0,19*sizeof(double));
    memset(c,0,19*sizeof(double));
    memset(d,0,19*sizeof(double));
    memset(alpha,0,19*sizeof(double));
    memset(mu,0,19*sizeof(double));
    memset(z,0,19*sizeof(double));
    memset(L,0,19*sizeof(double));
    memset(h,0,19*sizeof(double));
    //
    for(int i=0;i<19;i++)
    {   f[i] = data[i][1];
        t[i] = data[i][0];      }

    h[0] = t[1] - t[0];
    alpha[0] = 3.*(f[1]-f[0])/h[0]-3.*lslope;
    L[0] = 2.*h[0];
    mu[0] = 0.5;
    z[0] = alpha[0]/L[0];
    b[0] = lslope;
    // End Left Init

    for(int k=1;k<n-1;k++)
    {
        h[k] = t[k+1] - t[k];
        
        alpha[k] = (3./(h[k]*h[k-1]))*(f[k+1]*h[k-1]-f[k]*(h[k]+h[k-1])+f[k-1]*h[k]);
        
        L[k] = 2.*(h[k]+h[k-1])-h[k-1]*mu[k-1];

        mu[k] = h[k]/L[k];

        z[k] = (alpha[k]-h[k-1]*z[k-1])/L[k];
    }

    // RE init
    alpha[n-1] = 3.*rslope - 3.*(f[n-1]-f[n-2])/h[n-2];
    L[n-1] = h[n-2]*(2.-mu[n-2]);
    z[n-1] = (alpha[n-1]-h[n-2]*z[n-2])/L[n-1];
    c[n-1] = z[n-1];

    for(int j=n-2;j>=0;j--)
    {
        c[j] = z[j]-mu[j]*c[j+1];

        b[j] = (f[j+1]-f[j])/h[j] - h[j]*(c[j+1]+2.*c[j])/3.;

        d[j] = (c[j+1]-c[j])/(3.*h[j]);

        a[j] = f[j];
    }
    // RE END
    //
    //
    //
    


    for(int i=0;i<n-1;i++)
    {
        printf("knot sta/end : %lf\t%lf    \t\t   %lf\t%lf\t%lf\t%lf\n",
                t[i],t[i+1], a[i], b[i], c[i], d[i]);
    }
    printf("\n\n");

    for(int i=0;i<n-1;i++)
    {
        // order : high -> low
        printf("knot sta/end : %lf//%lf    \t\t %lf*x**3 + %lf*x**2 + %lf*x  + %lf \n",
                t[i],t[i+1],
                d[i], c[i] - 3.*t[i]*d[i], b[i] -2.*c[i]*t[i]+3.*t[i]*t[i]*d[i],
                    a[i] - b[i]*t[i]+c[i]*t[i]*t[i]-d[i]*t[i]*t[i]*t[i]);
    }










    printf("%d\n",10/3);





    return 0;
}
