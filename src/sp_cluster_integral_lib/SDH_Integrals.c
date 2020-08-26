#include<stdio.h>
#include<gsl/gsl_sf_expint.h>
#include<gsl/gsl_math.h>

// gsl_sf_expint_M_Ei -> gsl format at the end

/* < S | dxH | X >  ==  < S | dyH | Y > */
double SDH_Integral_x_12_case_1( double ASP, double RSP, double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d )
{
double res =(ASP*(d + RSP)*pow(3,0.5)*pow(d,-3)*pow(M_E,-((d + k1 + k2)*pow(RSP,-1)))*pow(RSP,2)*
(-4*b1*d2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,2) - 4*b1*d2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2) - b1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) - 
5*b1*c2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) - 5*a1*d2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + b1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3) - 
5*b1*c2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3) - 5*a1*d2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3) - b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) - 
a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) - 6*b1*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) - 6*a1*c2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + 
b1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) + a1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) - 6*b1*b2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) - 
6*a1*c2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) - b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) - a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) - 
7*a2*b1*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) - 7*a1*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + b1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) + 
a1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) - 7*a2*b1*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) - 7*a1*b2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) - 
a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) - a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) - 8*a1*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) + 
a2*b1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,6) + a1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,6) - 8*a1*a2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,6) - 
a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,7) + a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,7) + 4*b1*d2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,2) + 
4*b1*d2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2) + b1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) + 5*b1*c2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) + 
5*a1*d2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - b1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3) + 5*b1*c2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3) + 
5*a1*d2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3) + b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) + a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) + 
6*b1*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) + 6*a1*c2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - b1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) - 
a1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) + 6*b1*b2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) + 6*a1*c2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) + 
b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) + a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) + 7*a2*b1*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) + 
7*a1*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - b1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) - a1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) + 
7*a2*b1*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) + 7*a1*b2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) + a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) + 
a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) + 8*a1*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) - a2*b1*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,6) - 
a1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,6) + 8*a1*a2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,6) + a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,7) - 
a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,7) + 8*b1*d2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 8*b1*d2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 
8*b1*d2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,2) - 8*b1*d2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,2) - 
15*b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 15*a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 
15*b1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 15*a1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 
24*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 24*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 
24*b1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 24*a1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 
35*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) - 35*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 
35*a2*b1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 35*a1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) - 
48*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) + 48*a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) + 
15*b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 15*a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 
15*b1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 15*a1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 
24*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 24*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 
24*b1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 24*a1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 
35*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) + 35*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) - 
35*a2*b1*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) - 35*a1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) + 
48*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) - 48*a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) + 8*b1*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) + 
30*b1*c2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) + 30*a1*d2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 8*b1*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 
30*b1*c2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 30*a1*d2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 8*b1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,3) - 
30*b1*c2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,3) - 30*a1*d2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,3) + 
8*b1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,3) + 30*b1*c2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,3) + 
30*a1*d2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,3) - 72*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) - 
72*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) - 72*b1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) - 
72*a1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) - 140*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) - 
140*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) - 140*a2*b1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) - 
140*a1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) - 240*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) - 
240*a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) + 72*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 
72*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 72*b1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 
72*a1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 140*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) + 
140*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) + 140*a2*b1*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) + 
140*a1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) + 240*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) + 
240*a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) + 30*b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 30*a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 
144*b1*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 144*a1*c2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 30*b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 
30*a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 144*b1*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 144*a1*c2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 
30*b1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) + 30*a1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) + 
144*b1*b2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) + 144*a1*c2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) - 
30*b1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) - 30*a1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) - 
144*b1*b2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) - 144*a1*c2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) - 
420*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) - 420*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 
420*a2*b1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 420*a1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) - 
960*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) + 960*a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) + 
420*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 420*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) - 
420*a2*b1*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) - 420*a1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 
960*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) - 960*a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) + 
d1*(d2*((k2 + 2*RSP)*pow(M_E,k1*pow(RSP,-1)) - (k1 + 2*RSP)*pow(M_E,k2*pow(RSP,-1)) + (k1 - 2*RSP)*pow(M_E,(2*k1 + k2)*pow(RSP,-1)) + 
(-k2 + 2*RSP)*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))) - 4*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,2) - 4*b2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2) - 
b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) - 5*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3) - 
5*a2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3) - a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) + 
4*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,2) + 4*b2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2) + b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) + 
5*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3) + 5*a2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3) + 
a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) + 8*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 
8*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 8*b2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,2) - 8*b2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,2) - 
15*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 15*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 
15*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 15*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 
c2*(pow(M_E,(2*k1 + k2)*pow(RSP,-1))*(-3*k1*RSP + pow(k1,2) + 3*pow(RSP,2)) - pow(M_E,k2*pow(RSP,-1))*(3*k1*RSP + pow(k1,2) + 3*pow(RSP,2)) - 
pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*(-3*k2*RSP + pow(k2,2) + 3*pow(RSP,2)) + pow(M_E,k1*pow(RSP,-1))*(3*k2*RSP + pow(k2,2) + 3*pow(RSP,2))) + 
8*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) + 30*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 8*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 
30*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 8*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,3) - 30*a2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,3) + 
8*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,3) + 30*a2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,3) + 30*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 
30*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 30*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) - 30*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4)) + 
144*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 144*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 840*a2*b1*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 
840*a1*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 144*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 144*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 
840*a2*b1*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 840*a1*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 144*b1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,5) - 
144*a1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,5) - 840*a2*b1*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,5) - 
840*a1*b2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,5) + 144*b1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,5) + 
144*a1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,5) + 840*a2*b1*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,5) + 
840*a1*b2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,5) - 2880*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) - 
2880*a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) + 2880*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) + 
2880*a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) + 
c1*(-5*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) - 5*b2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3) - b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) - 
6*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) - 6*a2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) - 
a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) + 5*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) + 
5*b2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3) + b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) + 6*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 
b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) + 6*a2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) + a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - 
a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) - 15*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 15*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 
24*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 24*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 
15*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 15*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 
24*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 24*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 
d2*(pow(M_E,(2*k1 + k2)*pow(RSP,-1))*(-3*k1*RSP + pow(k1,2) + 3*pow(RSP,2)) - pow(M_E,k2*pow(RSP,-1))*(3*k1*RSP + pow(k1,2) + 3*pow(RSP,2)) - 
pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*(-3*k2*RSP + pow(k2,2) + 3*pow(RSP,2)) + pow(M_E,k1*pow(RSP,-1))*(3*k2*RSP + pow(k2,2) + 3*pow(RSP,2))) + 
30*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 30*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 30*b2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,3) + 
30*b2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,3) - 72*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) - 
72*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 72*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 
72*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 
c2*(pow(M_E,(2*k1 + k2)*pow(RSP,-1))*(-4*RSP*pow(k1,2) + pow(k1,3) + 8*k1*pow(RSP,2) - 8*pow(RSP,3)) - 
pow(M_E,k2*pow(RSP,-1))*(4*RSP*pow(k1,2) + pow(k1,3) + 8*k1*pow(RSP,2) + 8*pow(RSP,3)) + 
pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*(4*RSP*pow(k2,2) - pow(k2,3) - 8*k2*pow(RSP,2) + 8*pow(RSP,3)) + 
pow(M_E,k1*pow(RSP,-1))*(4*RSP*pow(k2,2) + pow(k2,3) + 8*k2*pow(RSP,2) + 8*pow(RSP,3))) + 30*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 
144*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 30*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 144*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 
30*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) + 144*a2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) - 30*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) - 
144*a2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) + 144*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 144*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 
144*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,5) + 144*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,5)) + 840*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) + 
840*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) + 5760*a1*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) - 840*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) - 
840*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) - 5760*a1*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) + 840*a2*b1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,6) + 
840*a1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,6) + 5760*a1*a2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,6) - 
840*a2*b1*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,6) - 840*a1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,6) - 
5760*a1*a2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,6) + 5760*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,7) - 5760*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,7) - 
5760*a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,7) + 5760*a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,7)))/2.;

return res;
}

double SDH_Integral_x_12_case_2_sub_1( double ASP, double RSP, double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d )
{
double res =-(ASP*pow(3,0.5)*pow(d,-3)*(d + RSP + d*pow(M_E,2*d*pow(RSP,-1)) - RSP*pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,-((2*d + k2)*pow(RSP,-1)))*pow(RSP,2)*
(-(d1*d2*k2*pow(M_E,d*pow(RSP,-1))) - 2*d1*d2*RSP*pow(M_E,d*pow(RSP,-1)) - 3*c2*d1*k2*RSP*pow(M_E,d*pow(RSP,-1)) + d*d1*d2*pow(M_E,k2*pow(RSP,-1)) + 
3*c2*d*d1*RSP*pow(M_E,k2*pow(RSP,-1)) + 2*d1*d2*RSP*pow(M_E,k2*pow(RSP,-1)) + c2*d1*pow(d,2)*pow(M_E,k2*pow(RSP,-1)) + 4*b2*d1*RSP*pow(d,2)*pow(M_E,k2*pow(RSP,-1)) + 
b2*d1*pow(d,3)*pow(M_E,k2*pow(RSP,-1)) + 5*a2*d1*RSP*pow(d,3)*pow(M_E,k2*pow(RSP,-1)) + 5*a1*d2*RSP*pow(d,3)*pow(M_E,k2*pow(RSP,-1)) + 
a2*d1*pow(d,4)*pow(M_E,k2*pow(RSP,-1)) + a1*d2*pow(d,4)*pow(M_E,k2*pow(RSP,-1)) + 6*a1*c2*RSP*pow(d,4)*pow(M_E,k2*pow(RSP,-1)) + 
a1*c2*pow(d,5)*pow(M_E,k2*pow(RSP,-1)) + 7*a1*b2*RSP*pow(d,5)*pow(M_E,k2*pow(RSP,-1)) + a1*b2*pow(d,6)*pow(M_E,k2*pow(RSP,-1)) + 
8*a1*a2*RSP*pow(d,6)*pow(M_E,k2*pow(RSP,-1)) + a1*a2*pow(d,7)*pow(M_E,k2*pow(RSP,-1)) - c2*d1*pow(M_E,d*pow(RSP,-1))*pow(k2,2) - 
4*b2*d1*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,2) - b2*d1*pow(M_E,d*pow(RSP,-1))*pow(k2,3) - 5*a2*d1*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,3) - 
5*a1*d2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,3) - a2*d1*pow(M_E,d*pow(RSP,-1))*pow(k2,4) - a1*d2*pow(M_E,d*pow(RSP,-1))*pow(k2,4) - 
6*a1*c2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,4) - a1*c2*pow(M_E,d*pow(RSP,-1))*pow(k2,5) - 7*a1*b2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,5) - 
a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k2,6) - 8*a1*a2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,6) - a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,7) - 
3*c2*d1*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) - 8*b2*d1*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) + 3*c2*d1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 
8*b2*d*d1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 15*a2*d1*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 15*a1*d2*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 
24*a1*c2*pow(d,3)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 35*a1*b2*pow(d,4)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 48*a1*a2*pow(d,5)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) - 
15*a2*d1*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 15*a1*d2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 24*a1*c2*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 
35*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) - 48*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) - 8*b2*d1*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) - 
30*a2*d1*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) - 30*a1*d2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 8*b2*d1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 
30*a2*d*d1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 30*a1*d*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 72*a1*c2*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 
140*a1*b2*pow(d,3)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 240*a1*a2*pow(d,4)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 
72*a1*c2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 140*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) - 
240*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) - 30*a2*d1*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) - 30*a1*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) - 
144*a1*c2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 144*a1*c2*d*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 30*a2*d1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 
30*a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 420*a1*b2*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 960*a1*a2*pow(d,3)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 
420*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) - 960*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) - 144*a1*c2*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) - 
840*a1*b2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 144*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 840*a1*b2*d*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 
2880*a1*a2*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 2880*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) + 
c1*(-3*d2*k2*RSP*pow(M_E,d*pow(RSP,-1)) + 3*d*d2*RSP*pow(M_E,k2*pow(RSP,-1)) + d2*pow(d,2)*pow(M_E,k2*pow(RSP,-1)) + 5*b2*RSP*pow(d,3)*pow(M_E,k2*pow(RSP,-1)) + 
b2*pow(d,4)*pow(M_E,k2*pow(RSP,-1)) + 6*a2*RSP*pow(d,4)*pow(M_E,k2*pow(RSP,-1)) + a2*pow(d,5)*pow(M_E,k2*pow(RSP,-1)) - d2*pow(M_E,d*pow(RSP,-1))*pow(k2,2) - 
5*b2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,3) - b2*pow(M_E,d*pow(RSP,-1))*pow(k2,4) - 6*a2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,4) - a2*pow(M_E,d*pow(RSP,-1))*pow(k2,5) - 
3*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) + 3*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 15*b2*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 
24*a2*pow(d,3)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) - 15*b2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 24*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 
30*b2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 30*b2*d*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 72*a2*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 
72*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + c2*(4*RSP*pow(d,2)*pow(M_E,k2*pow(RSP,-1)) + pow(d,3)*pow(M_E,k2*pow(RSP,-1)) + 
8*d*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 8*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 
pow(M_E,d*pow(RSP,-1))*(4*RSP*pow(k2,2) + pow(k2,3) + 8*k2*pow(RSP,2) + 8*pow(RSP,3))) - 30*b2*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) - 
144*a2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 30*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 144*a2*d*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 
144*a2*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 144*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5)) - 840*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) - 
5760*a1*a2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) + 840*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) + 5760*a1*a2*d*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) + 
b1*(4*d2*RSP*pow(d,2)*pow(M_E,k2*pow(RSP,-1)) + d2*pow(d,3)*pow(M_E,k2*pow(RSP,-1)) + 6*b2*RSP*pow(d,4)*pow(M_E,k2*pow(RSP,-1)) + 
b2*pow(d,5)*pow(M_E,k2*pow(RSP,-1)) + 7*a2*RSP*pow(d,5)*pow(M_E,k2*pow(RSP,-1)) + a2*pow(d,6)*pow(M_E,k2*pow(RSP,-1)) - 4*d2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,2) - 
d2*pow(M_E,d*pow(RSP,-1))*pow(k2,3) - 6*b2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,4) - b2*pow(M_E,d*pow(RSP,-1))*pow(k2,5) - 7*a2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,5) - 
a2*pow(M_E,d*pow(RSP,-1))*pow(k2,6) - 8*d2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) + 8*d*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 
24*b2*pow(d,3)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 35*a2*pow(d,4)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) - 24*b2*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 
35*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) - 8*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 8*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 
72*b2*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 140*a2*pow(d,3)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 72*b2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 
140*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) - 144*b2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 144*b2*d*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 
420*a2*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 420*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 
c2*(5*RSP*pow(d,3)*pow(M_E,k2*pow(RSP,-1)) + pow(d,4)*pow(M_E,k2*pow(RSP,-1)) + 15*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 
30*d*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 30*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 
pow(M_E,d*pow(RSP,-1))*(5*RSP*pow(k2,3) + pow(k2,4) + 15*pow(k2,2)*pow(RSP,2) + 30*k2*pow(RSP,3) + 30*pow(RSP,4))) - 
144*b2*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) - 840*a2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 144*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 
840*a2*d*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 840*a2*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) + 840*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6)) - 
5760*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(RSP,7) + 5760*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,7)))/2.;

return res;
}

double SDH_Integral_x_12_case_2_sub_2( double ASP, double RSP, double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d )
{ 
double res =-(ASP*(d + RSP)*pow(3,0.5)*pow(d,-3)*pow(M_E,-((2*d + k1)*pow(RSP,-1)))*pow(RSP,2)*
(d1*d2*k1*pow(M_E,d*pow(RSP,-1)) + 2*d1*d2*RSP*pow(M_E,d*pow(RSP,-1)) + 3*c2*d1*k1*RSP*pow(M_E,d*pow(RSP,-1)) - d*d1*d2*pow(M_E,k1*pow(RSP,-1)) - 
3*c2*d*d1*RSP*pow(M_E,k1*pow(RSP,-1)) - 2*d1*d2*RSP*pow(M_E,k1*pow(RSP,-1)) - c2*d1*pow(d,2)*pow(M_E,k1*pow(RSP,-1)) - 4*b2*d1*RSP*pow(d,2)*pow(M_E,k1*pow(RSP,-1)) - 
b2*d1*pow(d,3)*pow(M_E,k1*pow(RSP,-1)) - 5*a2*d1*RSP*pow(d,3)*pow(M_E,k1*pow(RSP,-1)) - 5*a1*d2*RSP*pow(d,3)*pow(M_E,k1*pow(RSP,-1)) - 
a2*d1*pow(d,4)*pow(M_E,k1*pow(RSP,-1)) - a1*d2*pow(d,4)*pow(M_E,k1*pow(RSP,-1)) - 6*a1*c2*RSP*pow(d,4)*pow(M_E,k1*pow(RSP,-1)) - 
a1*c2*pow(d,5)*pow(M_E,k1*pow(RSP,-1)) - 7*a1*b2*RSP*pow(d,5)*pow(M_E,k1*pow(RSP,-1)) - a1*b2*pow(d,6)*pow(M_E,k1*pow(RSP,-1)) - 
8*a1*a2*RSP*pow(d,6)*pow(M_E,k1*pow(RSP,-1)) - a1*a2*pow(d,7)*pow(M_E,k1*pow(RSP,-1)) + d*d1*d2*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 
3*c2*d*d1*RSP*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 2*d1*d2*RSP*pow(M_E,(2*d + k1)*pow(RSP,-1)) + c2*d1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 
4*b2*d1*RSP*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + b2*d1*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 5*a2*d1*RSP*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 
5*a1*d2*RSP*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + a2*d1*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + a1*d2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 
6*a1*c2*RSP*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + a1*c2*pow(d,5)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 7*a1*b2*RSP*pow(d,5)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + 
a1*b2*pow(d,6)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 8*a1*a2*RSP*pow(d,6)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + a1*a2*pow(d,7)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 
d1*d2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1)) + 2*d1*d2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1)) + 3*c2*d1*k1*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1)) + 
c2*d1*pow(M_E,d*pow(RSP,-1))*pow(k1,2) + 4*b2*d1*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,2) - c2*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2) + 
4*b2*d1*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2) + b2*d1*pow(M_E,d*pow(RSP,-1))*pow(k1,3) + 5*a2*d1*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,3) + 
5*a1*d2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,3) - b2*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3) + 5*a2*d1*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3) + 
5*a1*d2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3) + a2*d1*pow(M_E,d*pow(RSP,-1))*pow(k1,4) + a1*d2*pow(M_E,d*pow(RSP,-1))*pow(k1,4) + 
6*a1*c2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,4) - a2*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4) - a1*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4) + 
6*a1*c2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4) + a1*c2*pow(M_E,d*pow(RSP,-1))*pow(k1,5) + 7*a1*b2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,5) - 
a1*c2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5) + 7*a1*b2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5) + a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k1,6) + 
8*a1*a2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,6) - a1*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,6) + 8*a1*a2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,6) + 
a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,7) - a1*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,7) + 3*c2*d1*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) + 
8*b2*d1*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) - 3*c2*d1*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 8*b2*d*d1*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 
15*a2*d1*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 15*a1*d2*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 24*a1*c2*pow(d,3)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 
35*a1*b2*pow(d,4)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 48*a1*a2*pow(d,5)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) + 3*c2*d1*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) + 
8*b2*d*d1*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) + 15*a2*d1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) + 
15*a1*d2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) + 24*a1*c2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) + 
35*a1*b2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) + 48*a1*a2*pow(d,5)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) - 
3*c2*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,2) - 8*b2*d1*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,2) + 15*a2*d1*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 
15*a1*d2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 15*a2*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 
15*a1*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 24*a1*c2*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 
24*a1*c2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 35*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) - 
35*a1*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 48*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) - 
48*a1*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) + 8*b2*d1*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 30*a2*d1*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 
30*a1*d2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) - 8*b2*d1*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 30*a2*d*d1*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
30*a1*d*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 72*a1*c2*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 140*a1*b2*pow(d,3)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
240*a1*a2*pow(d,4)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 8*b2*d1*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) - 30*a2*d*d1*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) - 
30*a1*d*d2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) - 72*a1*c2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) - 
140*a1*b2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) - 240*a1*a2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) + 
8*b2*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,3) + 30*a2*d1*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,3) + 
30*a1*d2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,3) + 72*a1*c2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
72*a1*c2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 140*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 
140*a1*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 240*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) + 
240*a1*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) + 30*a2*d1*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 30*a1*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 
144*a1*c2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) - 144*a1*c2*d*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 30*a2*d1*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 
30*a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 420*a1*b2*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 960*a1*a2*pow(d,3)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 
144*a1*c2*d*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) + 30*a2*d1*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) + 30*a1*d2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) + 
420*a1*b2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) + 960*a1*a2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) - 
30*a2*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,4) - 30*a1*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,4) - 
144*a1*c2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,4) + 420*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) - 
420*a1*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 960*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) - 
960*a1*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) + 144*a1*c2*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 840*a1*b2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) - 
144*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 840*a1*b2*d*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 2880*a1*a2*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 
144*a1*c2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,5) - 840*a1*b2*d*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,5) - 
2880*a1*a2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,5) + 144*a1*c2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,5) + 
840*a1*b2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,5) + 2880*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) + 
2880*a1*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) + 
c1*(3*d2*k1*RSP*pow(M_E,d*pow(RSP,-1)) - 3*d*d2*RSP*pow(M_E,k1*pow(RSP,-1)) - d2*pow(d,2)*pow(M_E,k1*pow(RSP,-1)) - 5*b2*RSP*pow(d,3)*pow(M_E,k1*pow(RSP,-1)) - 
b2*pow(d,4)*pow(M_E,k1*pow(RSP,-1)) - 6*a2*RSP*pow(d,4)*pow(M_E,k1*pow(RSP,-1)) - a2*pow(d,5)*pow(M_E,k1*pow(RSP,-1)) - 
3*d*d2*RSP*pow(M_E,(2*d + k1)*pow(RSP,-1)) + d2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 5*b2*RSP*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + 
b2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 6*a2*RSP*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + a2*pow(d,5)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + 
3*d2*k1*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1)) + d2*pow(M_E,d*pow(RSP,-1))*pow(k1,2) - d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2) + 
5*b2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,3) + 5*b2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3) + b2*pow(M_E,d*pow(RSP,-1))*pow(k1,4) + 
6*a2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,4) - b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4) + 6*a2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4) + 
a2*pow(M_E,d*pow(RSP,-1))*pow(k1,5) - a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5) + 3*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) - 
3*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 15*b2*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 24*a2*pow(d,3)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) + 
3*d2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) + 15*b2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) + 
24*a2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) - 3*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,2) + 15*b2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 
15*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 24*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 
24*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 30*b2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) - 30*b2*d*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
72*a2*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 30*b2*d*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) - 
72*a2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) + 30*b2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,3) + 
72*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 72*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
c2*(pow(d,3)*(-1 + pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1)) - 4*RSP*pow(d,2)*(1 + pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1)) + 
8*d*(-1 + pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 8*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 8*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) + 
pow(M_E,(d + 2*k1)*pow(RSP,-1))*(4*RSP*pow(k1,2) - pow(k1,3) - 8*k1*pow(RSP,2) + 8*pow(RSP,3)) + 
pow(M_E,d*pow(RSP,-1))*(4*RSP*pow(k1,2) + pow(k1,3) + 8*k1*pow(RSP,2) + 8*pow(RSP,3))) + 30*b2*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 
144*a2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) - 30*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 144*a2*d*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 
30*b2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) + 144*a2*d*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) - 30*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,4) - 
144*a2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,4) + 144*a2*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) - 144*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 
144*a2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,5) + 144*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,5)) + 840*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) + 
5760*a1*a2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) - 840*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) - 5760*a1*a2*d*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) + 
840*a1*b2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,6) + 5760*a1*a2*d*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,6) - 
840*a1*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,6) - 5760*a1*a2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,6) + 
b1*(-4*d2*RSP*pow(d,2)*pow(M_E,k1*pow(RSP,-1)) - d2*pow(d,3)*pow(M_E,k1*pow(RSP,-1)) - 6*b2*RSP*pow(d,4)*pow(M_E,k1*pow(RSP,-1)) - 
b2*pow(d,5)*pow(M_E,k1*pow(RSP,-1)) - 7*a2*RSP*pow(d,5)*pow(M_E,k1*pow(RSP,-1)) - a2*pow(d,6)*pow(M_E,k1*pow(RSP,-1)) - 
4*d2*RSP*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + d2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 6*b2*RSP*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + 
b2*pow(d,5)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 7*a2*RSP*pow(d,5)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + a2*pow(d,6)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + 
4*d2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,2) + 4*d2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2) + d2*pow(M_E,d*pow(RSP,-1))*pow(k1,3) - 
d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3) + 6*b2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,4) + 6*b2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4) + 
b2*pow(M_E,d*pow(RSP,-1))*pow(k1,5) + 7*a2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,5) - b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5) + 
7*a2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5) + a2*pow(M_E,d*pow(RSP,-1))*pow(k1,6) - a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,6) + 
8*d2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) - 8*d*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 24*b2*pow(d,3)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 
35*a2*pow(d,4)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) + 8*d*d2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) + 
24*b2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) + 35*a2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) - 
8*d2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,2) + 24*b2*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 
24*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 35*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) - 
35*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 8*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) - 8*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
72*b2*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 140*a2*pow(d,3)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 8*d2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) - 
72*b2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) - 140*a2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) + 
8*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,3) + 72*b2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
72*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 140*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 
140*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 144*b2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) - 144*b2*d*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 
420*a2*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 144*b2*d*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) + 
420*a2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) - 144*b2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,4) + 
420*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) - 420*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 
c2*(pow(d,4)*(-1 + pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1)) - 5*RSP*pow(d,3)*(1 + pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1)) + 
15*pow(d,2)*(-1 + pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 30*d*(1 + pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
30*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 30*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) - 
pow(M_E,(d + 2*k1)*pow(RSP,-1))*(-5*RSP*pow(k1,3) + pow(k1,4) + 15*pow(k1,2)*pow(RSP,2) - 30*k1*pow(RSP,3) + 30*pow(RSP,4)) + 
pow(M_E,d*pow(RSP,-1))*(5*RSP*pow(k1,3) + pow(k1,4) + 15*pow(k1,2)*pow(RSP,2) + 30*k1*pow(RSP,3) + 30*pow(RSP,4))) + 
144*b2*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 840*a2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) - 144*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 
840*a2*d*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 144*b2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,5) - 840*a2*d*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,5) + 
144*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,5) + 840*a2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,5) + 840*a2*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) - 
840*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) + 840*a2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,6) - 840*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,6)) + 
5760*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(RSP,7) - 5760*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,7) - 5760*a1*a2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,7) + 
5760*a1*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,7)))/2.;

return res;
}

double SDH_Integral_x_12_case_3( double ASP, double RSP, double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d )
{
double res =-(ASP*pow(3,0.5)*pow(d,-3)*(d + RSP + d*pow(M_E,2*d*pow(RSP,-1)) - RSP*pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,-((d + k1 + k2)*pow(RSP,-1)))*pow(RSP,2)*
(4*b1*d2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,2) + b1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + 5*b1*c2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + 
5*a1*d2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + 
6*b1*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + 6*a1*c2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + 
a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + 7*a2*b1*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + 7*a1*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + 
a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) + a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) + 8*a1*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) + 
a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,7) - 4*b1*d2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,2) - b1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - 
5*b1*c2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - 5*a1*d2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 
a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 6*b1*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 6*a1*c2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 
b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - 7*a2*b1*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - 
7*a1*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) - a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) - 
8*a1*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) - a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,7) - 8*b1*d2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) + 
8*b1*d2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 15*b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 15*a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 
24*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 24*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 
35*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 35*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 
48*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) - 15*b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 
15*a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 24*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 
24*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 35*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) - 
35*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) - 48*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) - 8*b1*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
30*b1*c2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 30*a1*d2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) + 8*b1*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 
30*b1*c2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 30*a1*d2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 72*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
72*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 140*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 
140*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 240*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) - 
72*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 72*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 
140*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) - 140*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) - 
240*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) - 30*b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 30*a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 
144*b1*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 144*a1*c2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 30*b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 
30*a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 144*b1*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 144*a1*c2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 
420*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 420*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 
960*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) - 420*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) - 
420*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) - 960*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) + 
d1*(-(d2*(k2 + 2*RSP)*pow(M_E,k1*pow(RSP,-1))) + d2*(k1 + 2*RSP)*pow(M_E,k2*pow(RSP,-1)) + 4*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,2) + 
b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + 5*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) - 
4*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,2) - b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - 5*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - 
a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 8*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) + 8*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 
15*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 15*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 
c2*pow(M_E,k2*pow(RSP,-1))*(3*k1*RSP + pow(k1,2) + 3*pow(RSP,2)) - c2*pow(M_E,k1*pow(RSP,-1))*(3*k2*RSP + pow(k2,2) + 3*pow(RSP,2)) - 
8*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 30*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) + 8*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 
30*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 30*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 30*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4)) - 
144*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 144*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 840*a2*b1*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 
840*a1*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 144*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 144*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 
840*a2*b1*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 840*a1*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 2880*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) - 
2880*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) + c1*
(5*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + 6*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + 
a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) - 5*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 
6*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) + 15*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 
24*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 15*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 24*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 
d2*pow(M_E,k2*pow(RSP,-1))*(3*k1*RSP + pow(k1,2) + 3*pow(RSP,2)) - d2*pow(M_E,k1*pow(RSP,-1))*(3*k2*RSP + pow(k2,2) + 3*pow(RSP,2)) - 
30*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) + 30*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 72*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) - 
72*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + c2*pow(M_E,k2*pow(RSP,-1))*(4*RSP*pow(k1,2) + pow(k1,3) + 8*k1*pow(RSP,2) + 8*pow(RSP,3)) - 
c2*pow(M_E,k1*pow(RSP,-1))*(4*RSP*pow(k2,2) + pow(k2,3) + 8*k2*pow(RSP,2) + 8*pow(RSP,3)) - 30*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 
144*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 30*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 144*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 
144*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 144*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5)) - 840*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) - 
840*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) - 5760*a1*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) + 840*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) + 
840*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) + 5760*a1*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) - 5760*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,7) + 
5760*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,7)))/2.;

return res;
}

/* < S | dxH | X >  M_M_END */









/* < X | dxH | Z >  ==  < Y | dyH | Z > */
double SDH_Integral_x_24_case_1( double AP, double RP, double k1, double k2, double a2, double b2, double c2, double d2, double d )
{
double res =(3*AP*pow(d,-4)*pow(M_E,-((d + k1 + k2)*pow(RP,-1)))*pow(RP,2)*(3*d*RP + pow(d,2) + 3*pow(RP,2))*
(-10*c2*d2*k2*RP*pow(M_E,k1*pow(RP,-1)) - k2*pow(d2,2)*pow(M_E,k1*pow(RP,-1)) - 4*RP*pow(d2,2)*pow(M_E,k1*pow(RP,-1)) + 10*c2*d2*k1*RP*pow(M_E,k2*pow(RP,-1)) + 
k1*pow(d2,2)*pow(M_E,k2*pow(RP,-1)) + 4*RP*pow(d2,2)*pow(M_E,k2*pow(RP,-1)) - 3*RP*gsl_sf_expint_Ei(-(k1*pow(RP,-1)))*pow(d2,2)*pow(M_E,(k1 + k2)*pow(RP,-1)) + 
3*RP*gsl_sf_expint_Ei(k1*pow(RP,-1))*pow(d2,2)*pow(M_E,(k1 + k2)*pow(RP,-1)) + 3*RP*gsl_sf_expint_Ei(-(k2*pow(RP,-1)))*pow(d2,2)*pow(M_E,(k1 + k2)*pow(RP,-1)) - 
3*RP*gsl_sf_expint_Ei(k2*pow(RP,-1))*pow(d2,2)*pow(M_E,(k1 + k2)*pow(RP,-1)) - 10*c2*d2*k1*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1)) + 
k1*pow(d2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) - 4*RP*pow(d2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) + 10*c2*d2*k2*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) - 
k2*pow(d2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) + 4*RP*pow(d2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) + 2*c2*d2*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 
12*b2*d2*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 6*RP*pow(c2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 2*c2*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) - 
12*b2*d2*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) - 6*RP*pow(c2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) + 2*b2*d2*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 
14*b2*c2*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 14*a2*d2*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + pow(c2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 
2*b2*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) - 14*b2*c2*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) - 14*a2*d2*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + 
pow(c2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + 2*b2*c2*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 2*a2*d2*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 
16*a2*c2*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 8*RP*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 2*b2*c2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) + 
2*a2*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) - 16*a2*c2*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) - 
8*RP*pow(b2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) + 2*a2*c2*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 18*a2*b2*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 
pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 2*a2*c2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) - 18*a2*b2*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) + 
pow(b2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) + 2*a2*b2*pow(M_E,k2*pow(RP,-1))*pow(k1,6) + 10*RP*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,6) + 
2*a2*b2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6) - 10*RP*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6) + pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,7) + 
pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,7) - 2*c2*d2*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 12*b2*d2*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 
6*RP*pow(c2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 2*c2*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) + 12*b2*d2*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) + 
6*RP*pow(c2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) - 2*b2*d2*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 14*b2*c2*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 
14*a2*d2*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - pow(c2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 2*b2*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) + 
14*b2*c2*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) + 14*a2*d2*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - pow(c2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - 
2*b2*c2*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 2*a2*d2*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 16*a2*c2*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 
8*RP*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 2*b2*c2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) - 2*a2*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) + 
16*a2*c2*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) + 8*RP*pow(b2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) - 2*a2*c2*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 
18*a2*b2*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 2*a2*c2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) + 
18*a2*b2*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) - pow(b2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) - 2*a2*b2*pow(M_E,k1*pow(RP,-1))*pow(k2,6) - 
10*RP*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,6) - 2*a2*b2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6) + 10*RP*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6) - 
pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,7) - pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,7) - 16*c2*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 
30*b2*d2*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 15*k2*pow(c2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 16*c2*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 
30*b2*d2*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 15*k1*pow(c2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 16*c2*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) + 
30*b2*d2*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) + 15*k1*pow(c2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) - 
16*c2*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) - 30*b2*d2*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) - 
15*k2*pow(c2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) + 48*b2*c2*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
48*a2*d2*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 48*b2*c2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
48*a2*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 70*a2*c2*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
35*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 70*a2*c2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
35*pow(b2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 96*a2*b2*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 
96*a2*b2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 63*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 
63*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 48*b2*c2*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
48*a2*d2*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 48*b2*c2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
48*a2*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 70*a2*c2*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
35*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 70*a2*c2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
35*pow(b2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 96*a2*b2*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 
96*a2*b2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 63*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 
63*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 30*b2*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 96*b2*c2*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 
96*a2*d2*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 15*pow(c2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 30*b2*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 
96*b2*c2*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 96*a2*d2*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 15*pow(c2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 
30*b2*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 96*b2*c2*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 96*a2*d2*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 
15*pow(c2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) + 30*b2*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 96*b2*c2*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 
96*a2*d2*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 15*pow(c2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 
210*a2*c2*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 105*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
210*a2*c2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 105*pow(b2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
384*a2*b2*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 384*a2*b2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
315*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 315*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 
210*a2*c2*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 105*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
210*a2*c2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 105*pow(b2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
384*a2*b2*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 384*a2*b2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 
315*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 315*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 
96*b2*c2*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 96*a2*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 420*a2*c2*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 
210*k2*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 96*b2*c2*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 96*a2*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
420*a2*c2*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 210*k1*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 96*b2*c2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 
96*a2*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 420*a2*c2*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 
210*k1*pow(b2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 96*b2*c2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 
96*a2*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 420*a2*c2*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 
210*k2*pow(b2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 1152*a2*b2*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
1152*a2*b2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 1260*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 
1260*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 1152*a2*b2*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
1152*a2*b2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 1260*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 
1260*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 420*a2*c2*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 2304*a2*b2*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 
210*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 420*a2*c2*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 2304*a2*b2*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
210*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 420*a2*c2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 2304*a2*b2*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 
210*pow(b2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) + 420*a2*c2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 
2304*a2*b2*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 210*pow(b2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 
3780*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 3780*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 
3780*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 3780*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 
2304*a2*b2*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 7560*k2*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 2304*a2*b2*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 
7560*k1*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 2304*a2*b2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) + 
7560*k1*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 2304*a2*b2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) - 
7560*k2*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) - 7560*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) + 7560*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) - 
7560*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) + 7560*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7)))/2.;

return res;
}
double SDH_Integral_x_24_case_2_sub_1( double AP, double RP, double k1, double k2, double a2, double b2, double c2, double d2, double d )
{
double res =(3*AP*pow(d,-4)*pow(M_E,-((2*d + k2)*pow(RP,-1)))*pow(RP,2)*(pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1))) - 3*d*RP*(1 + pow(M_E,2*d*pow(RP,-1))) + 
3*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,2))*(10*c2*d2*k2*RP*pow(M_E,d*pow(RP,-1)) + k2*pow(d2,2)*pow(M_E,d*pow(RP,-1)) + 4*RP*pow(d2,2)*pow(M_E,d*pow(RP,-1)) - 
10*c2*d*d2*RP*pow(M_E,k2*pow(RP,-1)) - 2*c2*d2*pow(d,2)*pow(M_E,k2*pow(RP,-1)) - 12*b2*d2*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1)) - 
6*RP*pow(c2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1)) - 2*b2*d2*pow(d,3)*pow(M_E,k2*pow(RP,-1)) - 14*b2*c2*RP*pow(d,3)*pow(M_E,k2*pow(RP,-1)) - 
14*a2*d2*RP*pow(d,3)*pow(M_E,k2*pow(RP,-1)) - pow(c2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1)) - 2*b2*c2*pow(d,4)*pow(M_E,k2*pow(RP,-1)) - 
2*a2*d2*pow(d,4)*pow(M_E,k2*pow(RP,-1)) - 16*a2*c2*RP*pow(d,4)*pow(M_E,k2*pow(RP,-1)) - 8*RP*pow(b2,2)*pow(d,4)*pow(M_E,k2*pow(RP,-1)) - 
2*a2*c2*pow(d,5)*pow(M_E,k2*pow(RP,-1)) - 18*a2*b2*RP*pow(d,5)*pow(M_E,k2*pow(RP,-1)) - pow(b2,2)*pow(d,5)*pow(M_E,k2*pow(RP,-1)) - 
2*a2*b2*pow(d,6)*pow(M_E,k2*pow(RP,-1)) - 10*RP*pow(a2,2)*pow(d,6)*pow(M_E,k2*pow(RP,-1)) - pow(a2,2)*pow(d,7)*pow(M_E,k2*pow(RP,-1)) - 
d*pow(d2,2)*pow(M_E,k2*pow(RP,-1)) - 4*RP*pow(d2,2)*pow(M_E,k2*pow(RP,-1)) + 3*RP*gsl_sf_expint_Ei(-(d*pow(RP,-1)))*pow(d2,2)*pow(M_E,(d + k2)*pow(RP,-1)) - 
3*RP*gsl_sf_expint_Ei(-(k2*pow(RP,-1)))*pow(d2,2)*pow(M_E,(d + k2)*pow(RP,-1)) + 2*c2*d2*pow(M_E,d*pow(RP,-1))*pow(k2,2) + 12*b2*d2*RP*pow(M_E,d*pow(RP,-1))*pow(k2,2) + 
6*RP*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2) + 2*b2*d2*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 14*b2*c2*RP*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 
14*a2*d2*RP*pow(M_E,d*pow(RP,-1))*pow(k2,3) + pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 2*b2*c2*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 
2*a2*d2*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 16*a2*c2*RP*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 8*RP*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 
2*a2*c2*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 18*a2*b2*RP*pow(M_E,d*pow(RP,-1))*pow(k2,5) + pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 
2*a2*b2*pow(M_E,d*pow(RP,-1))*pow(k2,6) + 10*RP*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,6) + pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,7) + 
16*c2*d2*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 30*b2*d2*k2*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 15*k2*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 
16*c2*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 30*b2*d*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 15*d*pow(c2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 
48*b2*c2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 48*a2*d2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 70*a2*c2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 
35*pow(b2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 96*a2*b2*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 
63*pow(a2,2)*pow(d,5)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 48*b2*c2*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 48*a2*d2*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
70*a2*c2*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 35*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 96*a2*b2*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
63*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 30*b2*d2*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 96*b2*c2*k2*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 
96*a2*d2*k2*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 15*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 96*b2*c2*d*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 
30*b2*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 96*a2*d*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 15*pow(c2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 
210*a2*c2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 105*pow(b2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 384*a2*b2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 
315*pow(a2,2)*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 210*a2*c2*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
105*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 384*a2*b2*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
315*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 96*b2*c2*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 96*a2*d2*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 
420*a2*c2*k2*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 210*k2*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 96*b2*c2*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 
420*a2*c2*d*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 96*a2*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 210*d*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 
1152*a2*b2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 1260*pow(a2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
1152*a2*b2*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 1260*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 420*a2*c2*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 
2304*a2*b2*k2*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 210*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 420*a2*c2*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 
2304*a2*b2*d*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 210*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 3780*pow(a2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
3780*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 2304*a2*b2*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 7560*k2*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 
2304*a2*b2*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 7560*d*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 7560*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 
7560*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7)))/2.;

return res;
}
double SDH_Integral_x_24_case_2_sub_2( double AP, double RP, double k1, double k2, double a2, double b2, double c2, double d2, double d )
{
double res =(-3*AP*pow(d,-4)*pow(M_E,-((2*d + k1)*pow(RP,-1)))*pow(RP,2)*(3*d*RP + pow(d,2) + 3*pow(RP,2))*
(-10*c2*d2*k1*RP*pow(M_E,d*pow(RP,-1)) - k1*pow(d2,2)*pow(M_E,d*pow(RP,-1)) - 4*RP*pow(d2,2)*pow(M_E,d*pow(RP,-1)) + 10*c2*d*d2*RP*pow(M_E,k1*pow(RP,-1)) + 
2*c2*d2*pow(d,2)*pow(M_E,k1*pow(RP,-1)) + 12*b2*d2*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1)) + 6*RP*pow(c2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1)) + 
2*b2*d2*pow(d,3)*pow(M_E,k1*pow(RP,-1)) + 14*b2*c2*RP*pow(d,3)*pow(M_E,k1*pow(RP,-1)) + 14*a2*d2*RP*pow(d,3)*pow(M_E,k1*pow(RP,-1)) + 
pow(c2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1)) + 2*b2*c2*pow(d,4)*pow(M_E,k1*pow(RP,-1)) + 2*a2*d2*pow(d,4)*pow(M_E,k1*pow(RP,-1)) + 
16*a2*c2*RP*pow(d,4)*pow(M_E,k1*pow(RP,-1)) + 8*RP*pow(b2,2)*pow(d,4)*pow(M_E,k1*pow(RP,-1)) + 2*a2*c2*pow(d,5)*pow(M_E,k1*pow(RP,-1)) + 
18*a2*b2*RP*pow(d,5)*pow(M_E,k1*pow(RP,-1)) + pow(b2,2)*pow(d,5)*pow(M_E,k1*pow(RP,-1)) + 2*a2*b2*pow(d,6)*pow(M_E,k1*pow(RP,-1)) + 
10*RP*pow(a2,2)*pow(d,6)*pow(M_E,k1*pow(RP,-1)) + pow(a2,2)*pow(d,7)*pow(M_E,k1*pow(RP,-1)) + d*pow(d2,2)*pow(M_E,k1*pow(RP,-1)) + 
4*RP*pow(d2,2)*pow(M_E,k1*pow(RP,-1)) - 3*RP*gsl_sf_expint_Ei(-(d*pow(RP,-1)))*pow(d2,2)*pow(M_E,(d + k1)*pow(RP,-1)) + 
3*RP*gsl_sf_expint_Ei(d*pow(RP,-1))*pow(d2,2)*pow(M_E,(d + k1)*pow(RP,-1)) + 3*RP*gsl_sf_expint_Ei(-(k1*pow(RP,-1)))*pow(d2,2)*pow(M_E,(d + k1)*pow(RP,-1)) - 
3*RP*gsl_sf_expint_Ei(k1*pow(RP,-1))*pow(d2,2)*pow(M_E,(d + k1)*pow(RP,-1)) - 10*c2*d*d2*RP*pow(M_E,(2*d + k1)*pow(RP,-1)) + 
2*c2*d2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 12*b2*d2*RP*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 6*RP*pow(c2,2)*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 
2*b2*d2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 14*b2*c2*RP*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 14*a2*d2*RP*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 
pow(c2,2)*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 2*b2*c2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 2*a2*d2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 
16*a2*c2*RP*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 8*RP*pow(b2,2)*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 2*a2*c2*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 
18*a2*b2*RP*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1)) + pow(b2,2)*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 2*a2*b2*pow(d,6)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 
10*RP*pow(a2,2)*pow(d,6)*pow(M_E,(2*d + k1)*pow(RP,-1)) + pow(a2,2)*pow(d,7)*pow(M_E,(2*d + k1)*pow(RP,-1)) + d*pow(d2,2)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 
4*RP*pow(d2,2)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 10*c2*d2*k1*RP*pow(M_E,(d + 2*k1)*pow(RP,-1)) - k1*pow(d2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1)) + 
4*RP*pow(d2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1)) - 2*c2*d2*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 12*b2*d2*RP*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 
6*RP*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 2*c2*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) + 12*b2*d2*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) + 
6*RP*pow(c2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) - 2*b2*d2*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 14*b2*c2*RP*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 
14*a2*d2*RP*pow(M_E,d*pow(RP,-1))*pow(k1,3) - pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 2*b2*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) + 
14*b2*c2*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) + 14*a2*d2*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - pow(c2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - 
2*b2*c2*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 2*a2*d2*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 16*a2*c2*RP*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 
8*RP*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 2*b2*c2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) - 2*a2*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) + 
16*a2*c2*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) + 8*RP*pow(b2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) - 2*a2*c2*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 
18*a2*b2*RP*pow(M_E,d*pow(RP,-1))*pow(k1,5) - pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 2*a2*c2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) + 
18*a2*b2*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) - pow(b2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) - 2*a2*b2*pow(M_E,d*pow(RP,-1))*pow(k1,6) - 
10*RP*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,6) - 2*a2*b2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6) + 10*RP*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6) - 
pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,7) - pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,7) - 16*c2*d2*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 
30*b2*d2*k1*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 15*k1*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 16*c2*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
30*b2*d*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 15*d*pow(c2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 48*b2*c2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
48*a2*d2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 70*a2*c2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 35*pow(b2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
96*a2*b2*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 63*pow(a2,2)*pow(d,5)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 16*c2*d2*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 
30*b2*d*d2*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 15*d*pow(c2,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 
48*b2*c2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 48*a2*d2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 
70*a2*c2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 35*pow(b2,2)*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 
96*a2*b2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 63*pow(a2,2)*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) - 
16*c2*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) - 30*b2*d2*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) - 15*k1*pow(c2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) - 
48*b2*c2*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 48*a2*d2*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
48*b2*c2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 48*a2*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
70*a2*c2*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 35*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
70*a2*c2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 35*pow(b2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
96*a2*b2*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 96*a2*b2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
63*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 63*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 30*b2*d2*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 
96*b2*c2*k1*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 96*a2*d2*k1*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 15*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 
96*b2*c2*d*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 30*b2*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 96*a2*d*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 
15*pow(c2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 210*a2*c2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 105*pow(b2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 
384*a2*b2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 315*pow(a2,2)*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 96*b2*c2*d*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 
30*b2*d2*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 96*a2*d*d2*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 15*pow(c2,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 
210*a2*c2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 105*pow(b2,2)*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 
384*a2*b2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 315*pow(a2,2)*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) + 
30*b2*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 96*b2*c2*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 96*a2*d2*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 
15*pow(c2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 210*a2*c2*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
105*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 210*a2*c2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
105*pow(b2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 384*a2*b2*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
384*a2*b2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 315*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 
315*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 96*b2*c2*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 96*a2*d2*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
420*a2*c2*k1*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 210*k1*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 96*b2*c2*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 
420*a2*c2*d*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 96*a2*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 210*d*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 
1152*a2*b2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 1260*pow(a2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 96*b2*c2*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) + 
420*a2*c2*d*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) + 96*a2*d2*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) + 210*d*pow(b2,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) + 
1152*a2*b2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) + 1260*pow(a2,2)*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) - 
96*b2*c2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 96*a2*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 420*a2*c2*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 
210*k1*pow(b2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 1152*a2*b2*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
1152*a2*b2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 1260*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 
1260*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 420*a2*c2*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 2304*a2*b2*k1*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 
210*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 420*a2*c2*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 2304*a2*b2*d*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 
210*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 3780*pow(a2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 420*a2*c2*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) - 
2304*a2*b2*d*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) - 210*pow(b2,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) - 
3780*pow(a2,2)*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) + 420*a2*c2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 
2304*a2*b2*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 210*pow(b2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 
3780*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 3780*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 
2304*a2*b2*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 7560*k1*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 2304*a2*b2*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 
7560*d*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 2304*a2*b2*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,6) + 7560*d*pow(a2,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,6) - 
2304*a2*b2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 7560*k1*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 7560*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 
7560*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 7560*pow(a2,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,7) + 7560*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7)))/
2.;

return res;
}
double SDH_Integral_x_24_case_3( double AP, double RP, double k1, double k2, double a2, double b2, double c2, double d2, double d )
{
double res =(3*AP*pow(d,-4)*pow(M_E,-((d + k1 + k2)*pow(RP,-1)))*pow(RP,2)*(pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1))) - 3*d*RP*(1 + pow(M_E,2*d*pow(RP,-1))) + 
3*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,2))*(10*c2*d2*k2*RP*pow(M_E,k1*pow(RP,-1)) + k2*pow(d2,2)*pow(M_E,k1*pow(RP,-1)) + 4*RP*pow(d2,2)*pow(M_E,k1*pow(RP,-1)) - 
10*c2*d2*k1*RP*pow(M_E,k2*pow(RP,-1)) - k1*pow(d2,2)*pow(M_E,k2*pow(RP,-1)) - 4*RP*pow(d2,2)*pow(M_E,k2*pow(RP,-1)) + 
3*RP*gsl_sf_expint_Ei(-(k1*pow(RP,-1)))*pow(d2,2)*pow(M_E,(k1 + k2)*pow(RP,-1)) - 3*RP*gsl_sf_expint_Ei(-(k2*pow(RP,-1)))*pow(d2,2)*pow(M_E,(k1 + k2)*pow(RP,-1)) - 
2*c2*d2*pow(M_E,k2*pow(RP,-1))*pow(k1,2) - 12*b2*d2*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,2) - 6*RP*pow(c2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) - 
2*b2*d2*pow(M_E,k2*pow(RP,-1))*pow(k1,3) - 14*b2*c2*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,3) - 14*a2*d2*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,3) - 
pow(c2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) - 2*b2*c2*pow(M_E,k2*pow(RP,-1))*pow(k1,4) - 2*a2*d2*pow(M_E,k2*pow(RP,-1))*pow(k1,4) - 
16*a2*c2*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,4) - 8*RP*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) - 2*a2*c2*pow(M_E,k2*pow(RP,-1))*pow(k1,5) - 
18*a2*b2*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,5) - pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) - 2*a2*b2*pow(M_E,k2*pow(RP,-1))*pow(k1,6) - 
10*RP*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,6) - pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,7) + 2*c2*d2*pow(M_E,k1*pow(RP,-1))*pow(k2,2) + 
12*b2*d2*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,2) + 6*RP*pow(c2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) + 2*b2*d2*pow(M_E,k1*pow(RP,-1))*pow(k2,3) + 
14*b2*c2*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,3) + 14*a2*d2*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,3) + pow(c2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) + 
2*b2*c2*pow(M_E,k1*pow(RP,-1))*pow(k2,4) + 2*a2*d2*pow(M_E,k1*pow(RP,-1))*pow(k2,4) + 16*a2*c2*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,4) + 
8*RP*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) + 2*a2*c2*pow(M_E,k1*pow(RP,-1))*pow(k2,5) + 18*a2*b2*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,5) + 
pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) + 2*a2*b2*pow(M_E,k1*pow(RP,-1))*pow(k2,6) + 10*RP*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,6) + 
pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,7) + 16*c2*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 30*b2*d2*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
15*k2*pow(c2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 16*c2*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 30*b2*d2*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 
15*k1*pow(c2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 48*b2*c2*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 48*a2*d2*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
70*a2*c2*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 35*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 96*a2*b2*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
63*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 48*b2*c2*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 48*a2*d2*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
70*a2*c2*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 35*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 96*a2*b2*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
63*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 30*b2*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 96*b2*c2*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 
96*a2*d2*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 15*pow(c2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 30*b2*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 
96*b2*c2*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 96*a2*d2*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 15*pow(c2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 
210*a2*c2*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 105*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
384*a2*b2*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 315*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 
210*a2*c2*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 105*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
384*a2*b2*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 315*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 96*b2*c2*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 
96*a2*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 420*a2*c2*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 210*k2*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 
96*b2*c2*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 96*a2*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 420*a2*c2*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 
210*k1*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 1152*a2*b2*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
1260*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 1152*a2*b2*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 
1260*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 420*a2*c2*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 2304*a2*b2*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 
210*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 420*a2*c2*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 2304*a2*b2*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 
210*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 3780*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 
3780*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 2304*a2*b2*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 7560*k2*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 
2304*a2*b2*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 7560*k1*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 7560*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 
7560*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7)))/2.;

return res;
}









/* < S | dzH | S > */

double SDH_Integral_z_11_case_1( double AS, double RS, double k1, double k2, double a1, double b1, double c1, double d1, double d )
{
double res =(AS*RS*(d + RS)*pow(d,-2)*pow(M_E,-(d*pow(RS,-1)))*(-(pow(d1,2)*(RS + (k1 - RS)*pow(M_E,k1*pow(RS,-1)))) - 
pow(d1,2)*pow(M_E,-(k1*pow(RS,-1)))*(k1 + RS - RS*pow(M_E,k1*pow(RS,-1))) + pow(d1,2)*(RS + (k2 - RS)*pow(M_E,k2*pow(RS,-1))) + 
pow(d1,2)*pow(M_E,-(k2*pow(RS,-1)))*(k2 + RS - RS*pow(M_E,k2*pow(RS,-1))) - 2*c1*d1*pow(M_E,k1*pow(RS,-1))*(-2*k1*RS + pow(k1,2) + 2*pow(RS,2)) + 
2*c1*d1*pow(M_E,k2*pow(RS,-1))*(-2*k2*RS + pow(k2,2) + 2*pow(RS,2)) - 
2*c1*d1*pow(M_E,-(k1*pow(RS,-1)))*(2*k1*RS + pow(k1,2) - 2*(-1 + pow(M_E,k1*pow(RS,-1)))*pow(RS,2)) + 
2*c1*d1*pow(M_E,-(k2*pow(RS,-1)))*(2*k2*RS + pow(k2,2) - 2*(-1 + pow(M_E,k2*pow(RS,-1)))*pow(RS,2)) - 
2*b1*d1*pow(M_E,k1*pow(RS,-1))*(-3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) - 6*pow(RS,3)) - 
pow(c1,2)*pow(M_E,k1*pow(RS,-1))*(-3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) - 6*pow(RS,3)) + 
2*b1*d1*pow(M_E,k2*pow(RS,-1))*(-3*RS*pow(k2,2) + pow(k2,3) + 6*k2*pow(RS,2) - 6*pow(RS,3)) + 
pow(c1,2)*pow(M_E,k2*pow(RS,-1))*(-3*RS*pow(k2,2) + pow(k2,3) + 6*k2*pow(RS,2) - 6*pow(RS,3)) - 
2*b1*d1*pow(M_E,-(k1*pow(RS,-1)))*(3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) - 6*(-1 + pow(M_E,k1*pow(RS,-1)))*pow(RS,3)) - 
pow(c1,2)*pow(M_E,-(k1*pow(RS,-1)))*(3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) - 6*(-1 + pow(M_E,k1*pow(RS,-1)))*pow(RS,3)) + 
2*b1*d1*pow(M_E,-(k2*pow(RS,-1)))*(3*RS*pow(k2,2) + pow(k2,3) + 6*k2*pow(RS,2) - 6*(-1 + pow(M_E,k2*pow(RS,-1)))*pow(RS,3)) + 
pow(c1,2)*pow(M_E,-(k2*pow(RS,-1)))*(3*RS*pow(k2,2) + pow(k2,3) + 6*k2*pow(RS,2) - 6*(-1 + pow(M_E,k2*pow(RS,-1)))*pow(RS,3)) - 
2*b1*c1*pow(M_E,k1*pow(RS,-1))*(-4*RS*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RS,2) - 24*k1*pow(RS,3) + 24*pow(RS,4)) - 
2*a1*d1*pow(M_E,k1*pow(RS,-1))*(-4*RS*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RS,2) - 24*k1*pow(RS,3) + 24*pow(RS,4)) + 
2*b1*c1*pow(M_E,k2*pow(RS,-1))*(-4*RS*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RS,2) - 24*k2*pow(RS,3) + 24*pow(RS,4)) + 
2*a1*d1*pow(M_E,k2*pow(RS,-1))*(-4*RS*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RS,2) - 24*k2*pow(RS,3) + 24*pow(RS,4)) - 
2*b1*c1*pow(M_E,-(k1*pow(RS,-1)))*(4*RS*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RS,2) + 24*k1*pow(RS,3) - 24*(-1 + pow(M_E,k1*pow(RS,-1)))*pow(RS,4)) - 
2*a1*d1*pow(M_E,-(k1*pow(RS,-1)))*(4*RS*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RS,2) + 24*k1*pow(RS,3) - 24*(-1 + pow(M_E,k1*pow(RS,-1)))*pow(RS,4)) + 
2*b1*c1*pow(M_E,-(k2*pow(RS,-1)))*(4*RS*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RS,2) + 24*k2*pow(RS,3) - 24*(-1 + pow(M_E,k2*pow(RS,-1)))*pow(RS,4)) + 
2*a1*d1*pow(M_E,-(k2*pow(RS,-1)))*(4*RS*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RS,2) + 24*k2*pow(RS,3) - 24*(-1 + pow(M_E,k2*pow(RS,-1)))*pow(RS,4)) - 
2*a1*c1*pow(M_E,k1*pow(RS,-1))*(-5*RS*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RS,2) - 60*pow(k1,2)*pow(RS,3) + 120*k1*pow(RS,4) - 120*pow(RS,5)) - 
pow(b1,2)*pow(M_E,k1*pow(RS,-1))*(-5*RS*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RS,2) - 60*pow(k1,2)*pow(RS,3) + 120*k1*pow(RS,4) - 120*pow(RS,5)) + 
2*a1*c1*pow(M_E,k2*pow(RS,-1))*(-5*RS*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RS,2) - 60*pow(k2,2)*pow(RS,3) + 120*k2*pow(RS,4) - 120*pow(RS,5)) + 
pow(b1,2)*pow(M_E,k2*pow(RS,-1))*(-5*RS*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RS,2) - 60*pow(k2,2)*pow(RS,3) + 120*k2*pow(RS,4) - 120*pow(RS,5)) - 
2*a1*c1*pow(M_E,-(k1*pow(RS,-1)))*(5*RS*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RS,2) + 60*pow(k1,2)*pow(RS,3) + 120*k1*pow(RS,4) - 
120*(-1 + pow(M_E,k1*pow(RS,-1)))*pow(RS,5)) - pow(b1,2)*pow(M_E,-(k1*pow(RS,-1)))*
(5*RS*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RS,2) + 60*pow(k1,2)*pow(RS,3) + 120*k1*pow(RS,4) - 120*(-1 + pow(M_E,k1*pow(RS,-1)))*pow(RS,5)) + 
2*a1*c1*pow(M_E,-(k2*pow(RS,-1)))*(5*RS*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RS,2) + 60*pow(k2,2)*pow(RS,3) + 120*k2*pow(RS,4) - 
120*(-1 + pow(M_E,k2*pow(RS,-1)))*pow(RS,5)) + pow(b1,2)*pow(M_E,-(k2*pow(RS,-1)))*
(5*RS*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RS,2) + 60*pow(k2,2)*pow(RS,3) + 120*k2*pow(RS,4) - 120*(-1 + pow(M_E,k2*pow(RS,-1)))*pow(RS,5)) - 
2*a1*b1*pow(M_E,k1*pow(RS,-1))*(-6*RS*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RS,2) - 120*pow(k1,3)*pow(RS,3) + 360*pow(k1,2)*pow(RS,4) - 720*k1*pow(RS,5) + 
720*pow(RS,6)) + 2*a1*b1*pow(M_E,k2*pow(RS,-1))*(-6*RS*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RS,2) - 120*pow(k2,3)*pow(RS,3) + 360*pow(k2,2)*pow(RS,4) - 
720*k2*pow(RS,5) + 720*pow(RS,6)) - 2*a1*b1*pow(M_E,-(k1*pow(RS,-1)))*
(6*RS*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RS,2) + 120*pow(k1,3)*pow(RS,3) + 360*pow(k1,2)*pow(RS,4) + 720*k1*pow(RS,5) - 
720*(-1 + pow(M_E,k1*pow(RS,-1)))*pow(RS,6)) + 2*a1*b1*pow(M_E,-(k2*pow(RS,-1)))*
(6*RS*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RS,2) + 120*pow(k2,3)*pow(RS,3) + 360*pow(k2,2)*pow(RS,4) + 720*k2*pow(RS,5) - 
720*(-1 + pow(M_E,k2*pow(RS,-1)))*pow(RS,6)) - pow(a1,2)*pow(M_E,k1*pow(RS,-1))*
(-7*RS*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RS,2) - 210*pow(k1,4)*pow(RS,3) + 840*pow(k1,3)*pow(RS,4) - 2520*pow(k1,2)*pow(RS,5) + 5040*k1*pow(RS,6) - 
5040*pow(RS,7)) + pow(a1,2)*pow(M_E,k2*pow(RS,-1))*(-7*RS*pow(k2,6) + pow(k2,7) + 42*pow(k2,5)*pow(RS,2) - 210*pow(k2,4)*pow(RS,3) + 840*pow(k2,3)*pow(RS,4) - 
2520*pow(k2,2)*pow(RS,5) + 5040*k2*pow(RS,6) - 5040*pow(RS,7)) - 
pow(a1,2)*pow(M_E,-(k1*pow(RS,-1)))*(7*RS*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RS,2) + 210*pow(k1,4)*pow(RS,3) + 840*pow(k1,3)*pow(RS,4) + 
2520*pow(k1,2)*pow(RS,5) + 5040*k1*pow(RS,6) - 5040*(-1 + pow(M_E,k1*pow(RS,-1)))*pow(RS,7)) + 
pow(a1,2)*pow(M_E,-(k2*pow(RS,-1)))*(7*RS*pow(k2,6) + pow(k2,7) + 42*pow(k2,5)*pow(RS,2) + 210*pow(k2,4)*pow(RS,3) + 840*pow(k2,3)*pow(RS,4) + 
2520*pow(k2,2)*pow(RS,5) + 5040*k2*pow(RS,6) - 5040*(-1 + pow(M_E,k2*pow(RS,-1)))*pow(RS,7))))/2.;

return res;
}

double SDH_Integral_z_11_case_2_sub_1( double AS, double RS, double k1, double k2, double a1, double b1, double c1, double d1, double d )
{
double res =-(AS*RS*pow(d,-2)*(d + RS + d*pow(M_E,2*d*pow(RS,-1)) - RS*pow(M_E,2*d*pow(RS,-1)))*pow(M_E,-((2*d + k2)*pow(RS,-1)))*
(-(k2*pow(d1,2)*pow(M_E,d*pow(RS,-1))) - RS*pow(d1,2)*pow(M_E,d*pow(RS,-1)) + 8*a1*d1*RS*pow(d,3)*pow(M_E,k2*pow(RS,-1)) + 2*a1*d1*pow(d,4)*pow(M_E,k2*pow(RS,-1)) + 
7*RS*pow(a1,2)*pow(d,6)*pow(M_E,k2*pow(RS,-1)) + pow(a1,2)*pow(d,7)*pow(M_E,k2*pow(RS,-1)) + d*pow(d1,2)*pow(M_E,k2*pow(RS,-1)) + 
RS*pow(d1,2)*pow(M_E,k2*pow(RS,-1)) - 8*a1*d1*RS*pow(M_E,d*pow(RS,-1))*pow(k2,3) - 2*a1*d1*pow(M_E,d*pow(RS,-1))*pow(k2,4) - 
7*RS*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k2,6) - pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k2,7) + 24*a1*d1*pow(d,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,2) + 
42*pow(a1,2)*pow(d,5)*pow(M_E,k2*pow(RS,-1))*pow(RS,2) - 24*a1*d1*pow(M_E,d*pow(RS,-1))*pow(k2,2)*pow(RS,2) - 
42*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k2,5)*pow(RS,2) - 48*a1*d1*k2*pow(M_E,d*pow(RS,-1))*pow(RS,3) + 48*a1*d*d1*pow(M_E,k2*pow(RS,-1))*pow(RS,3) + 
210*pow(a1,2)*pow(d,4)*pow(M_E,k2*pow(RS,-1))*pow(RS,3) - 210*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k2,4)*pow(RS,3) + 
pow(c1,2)*(3*RS*pow(d,2)*pow(M_E,k2*pow(RS,-1)) + pow(d,3)*pow(M_E,k2*pow(RS,-1)) + 6*d*pow(M_E,k2*pow(RS,-1))*pow(RS,2) + 6*pow(M_E,k2*pow(RS,-1))*pow(RS,3) - 
pow(M_E,d*pow(RS,-1))*(3*RS*pow(k2,2) + pow(k2,3) + 6*k2*pow(RS,2) + 6*pow(RS,3))) - 48*a1*d1*pow(M_E,d*pow(RS,-1))*pow(RS,4) + 
48*a1*d1*pow(M_E,k2*pow(RS,-1))*pow(RS,4) + 840*pow(a1,2)*pow(d,3)*pow(M_E,k2*pow(RS,-1))*pow(RS,4) - 840*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k2,3)*pow(RS,4) + 
2520*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,5) - 2520*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k2,2)*pow(RS,5) + 
2*c1*(-2*d1*k2*RS*pow(M_E,d*pow(RS,-1)) + 2*d*d1*RS*pow(M_E,k2*pow(RS,-1)) + d1*pow(d,2)*pow(M_E,k2*pow(RS,-1)) + 5*a1*RS*pow(d,4)*pow(M_E,k2*pow(RS,-1)) + 
a1*pow(d,5)*pow(M_E,k2*pow(RS,-1)) - d1*pow(M_E,d*pow(RS,-1))*pow(k2,2) - 5*a1*RS*pow(M_E,d*pow(RS,-1))*pow(k2,4) - a1*pow(M_E,d*pow(RS,-1))*pow(k2,5) - 
2*d1*pow(M_E,d*pow(RS,-1))*pow(RS,2) + 2*d1*pow(M_E,k2*pow(RS,-1))*pow(RS,2) + 20*a1*pow(d,3)*pow(M_E,k2*pow(RS,-1))*pow(RS,2) - 
20*a1*pow(M_E,d*pow(RS,-1))*pow(k2,3)*pow(RS,2) + 60*a1*pow(d,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,3) - 60*a1*pow(M_E,d*pow(RS,-1))*pow(k2,2)*pow(RS,3) - 
120*a1*k2*pow(M_E,d*pow(RS,-1))*pow(RS,4) + 120*a1*d*pow(M_E,k2*pow(RS,-1))*pow(RS,4) + 
b1*(4*RS*pow(d,3)*pow(M_E,k2*pow(RS,-1)) + pow(d,4)*pow(M_E,k2*pow(RS,-1)) + 12*pow(d,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,2) + 24*d*pow(M_E,k2*pow(RS,-1))*pow(RS,3) + 
24*pow(M_E,k2*pow(RS,-1))*pow(RS,4) - pow(M_E,d*pow(RS,-1))*(4*RS*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RS,2) + 24*k2*pow(RS,3) + 24*pow(RS,4))) - 
120*a1*pow(M_E,d*pow(RS,-1))*pow(RS,5) + 120*a1*pow(M_E,k2*pow(RS,-1))*pow(RS,5)) + 
pow(b1,2)*(5*RS*pow(d,4)*pow(M_E,k2*pow(RS,-1)) + pow(d,5)*pow(M_E,k2*pow(RS,-1)) + 20*pow(d,3)*pow(M_E,k2*pow(RS,-1))*pow(RS,2) + 
60*pow(d,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,3) + 120*d*pow(M_E,k2*pow(RS,-1))*pow(RS,4) + 120*pow(M_E,k2*pow(RS,-1))*pow(RS,5) - 
pow(M_E,d*pow(RS,-1))*(5*RS*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RS,2) + 60*pow(k2,2)*pow(RS,3) + 120*k2*pow(RS,4) + 120*pow(RS,5))) - 
5040*k2*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(RS,6) + 5040*d*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,6) + 
2*b1*(d1*(3*RS*pow(d,2)*pow(M_E,k2*pow(RS,-1)) + pow(d,3)*pow(M_E,k2*pow(RS,-1)) + 6*d*pow(M_E,k2*pow(RS,-1))*pow(RS,2) + 6*pow(M_E,k2*pow(RS,-1))*pow(RS,3) - 
pow(M_E,d*pow(RS,-1))*(3*RS*pow(k2,2) + pow(k2,3) + 6*k2*pow(RS,2) + 6*pow(RS,3))) + 
a1*(6*RS*pow(d,5)*pow(M_E,k2*pow(RS,-1)) + pow(d,6)*pow(M_E,k2*pow(RS,-1)) + 30*pow(d,4)*pow(M_E,k2*pow(RS,-1))*pow(RS,2) + 
120*pow(d,3)*pow(M_E,k2*pow(RS,-1))*pow(RS,3) + 360*pow(d,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,4) + 720*d*pow(M_E,k2*pow(RS,-1))*pow(RS,5) + 
720*pow(M_E,k2*pow(RS,-1))*pow(RS,6) - pow(M_E,d*pow(RS,-1))*
(6*RS*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RS,2) + 120*pow(k2,3)*pow(RS,3) + 360*pow(k2,2)*pow(RS,4) + 720*k2*pow(RS,5) + 720*pow(RS,6)))) - 
5040*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(RS,7) + 5040*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,7)))/2.;

return res;
}

double SDH_Integral_z_11_case_2_sub_2( double AS, double RS, double k1, double k2, double a1, double b1, double c1, double d1, double d )
{
double res =(AS*(d + RS)*pow(d,-2)*pow(M_E,-(d*pow(RS,-1)))*(RS*pow(d1,2)*(d*pow(M_E,d*pow(RS,-1)) - RS*pow(M_E,d*pow(RS,-1)) + (-k1 + RS)*pow(M_E,k1*pow(RS,-1))) + 
RS*pow(d1,2)*(-((k1 + RS)*pow(M_E,d*pow(RS,-1))) + d*pow(M_E,k1*pow(RS,-1)) + RS*pow(M_E,k1*pow(RS,-1)))*pow(M_E,-((d + k1)*pow(RS,-1))) + 
2*c1*d1*RS*(-2*d*RS*pow(M_E,d*pow(RS,-1)) + pow(d,2)*pow(M_E,d*pow(RS,-1)) + 2*pow(M_E,d*pow(RS,-1))*pow(RS,2) - 
pow(M_E,k1*pow(RS,-1))*(-2*k1*RS + pow(k1,2) + 2*pow(RS,2))) + 
2*c1*d1*RS*pow(M_E,-((d + k1)*pow(RS,-1)))*(2*d*RS*pow(M_E,k1*pow(RS,-1)) + pow(d,2)*pow(M_E,k1*pow(RS,-1)) + 2*pow(M_E,k1*pow(RS,-1))*pow(RS,2) - 
pow(M_E,d*pow(RS,-1))*(2*k1*RS + pow(k1,2) + 2*pow(RS,2))) + 
2*b1*d1*RS*(-3*RS*pow(d,2)*pow(M_E,d*pow(RS,-1)) + pow(d,3)*pow(M_E,d*pow(RS,-1)) + 6*d*pow(M_E,d*pow(RS,-1))*pow(RS,2) - 
pow(M_E,k1*pow(RS,-1))*(-3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) - 6*pow(RS,3)) - 6*pow(M_E,d*pow(RS,-1))*pow(RS,3)) + 
RS*pow(c1,2)*(-3*RS*pow(d,2)*pow(M_E,d*pow(RS,-1)) + pow(d,3)*pow(M_E,d*pow(RS,-1)) + 6*d*pow(M_E,d*pow(RS,-1))*pow(RS,2) - 
pow(M_E,k1*pow(RS,-1))*(-3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) - 6*pow(RS,3)) - 6*pow(M_E,d*pow(RS,-1))*pow(RS,3)) + 
2*b1*d1*RS*pow(M_E,-((d + k1)*pow(RS,-1)))*(3*RS*pow(d,2)*pow(M_E,k1*pow(RS,-1)) + pow(d,3)*pow(M_E,k1*pow(RS,-1)) + 6*d*pow(M_E,k1*pow(RS,-1))*pow(RS,2) + 
6*pow(M_E,k1*pow(RS,-1))*pow(RS,3) - pow(M_E,d*pow(RS,-1))*(3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) + 6*pow(RS,3))) + 
RS*pow(c1,2)*pow(M_E,-((d + k1)*pow(RS,-1)))*(3*RS*pow(d,2)*pow(M_E,k1*pow(RS,-1)) + pow(d,3)*pow(M_E,k1*pow(RS,-1)) + 6*d*pow(M_E,k1*pow(RS,-1))*pow(RS,2) + 
6*pow(M_E,k1*pow(RS,-1))*pow(RS,3) - pow(M_E,d*pow(RS,-1))*(3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) + 6*pow(RS,3))) + 
2*b1*c1*RS*(-4*RS*pow(d,3)*pow(M_E,d*pow(RS,-1)) + pow(d,4)*pow(M_E,d*pow(RS,-1)) + 12*pow(d,2)*pow(M_E,d*pow(RS,-1))*pow(RS,2) - 
24*d*pow(M_E,d*pow(RS,-1))*pow(RS,3) + 24*pow(M_E,d*pow(RS,-1))*pow(RS,4) - 
pow(M_E,k1*pow(RS,-1))*(-4*RS*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RS,2) - 24*k1*pow(RS,3) + 24*pow(RS,4))) + 
2*a1*d1*RS*(-4*RS*pow(d,3)*pow(M_E,d*pow(RS,-1)) + pow(d,4)*pow(M_E,d*pow(RS,-1)) + 12*pow(d,2)*pow(M_E,d*pow(RS,-1))*pow(RS,2) - 
24*d*pow(M_E,d*pow(RS,-1))*pow(RS,3) + 24*pow(M_E,d*pow(RS,-1))*pow(RS,4) - 
pow(M_E,k1*pow(RS,-1))*(-4*RS*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RS,2) - 24*k1*pow(RS,3) + 24*pow(RS,4))) - 
2*b1*c1*(RS*pow(M_E,-(d*pow(RS,-1)))*(-4*RS*pow(d,3) - pow(d,4) - 12*pow(d,2)*pow(RS,2) - 24*d*pow(RS,3) + 24*(-1 + pow(M_E,d*pow(RS,-1)))*pow(RS,4)) + 
RS*pow(M_E,-(k1*pow(RS,-1)))*(4*RS*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RS,2) + 24*k1*pow(RS,3) - 24*(-1 + pow(M_E,k1*pow(RS,-1)))*pow(RS,4))) - 
2*a1*d1*(RS*pow(M_E,-(d*pow(RS,-1)))*(-4*RS*pow(d,3) - pow(d,4) - 12*pow(d,2)*pow(RS,2) - 24*d*pow(RS,3) + 24*(-1 + pow(M_E,d*pow(RS,-1)))*pow(RS,4)) + 
RS*pow(M_E,-(k1*pow(RS,-1)))*(4*RS*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RS,2) + 24*k1*pow(RS,3) - 24*(-1 + pow(M_E,k1*pow(RS,-1)))*pow(RS,4))) + 
2*a1*c1*RS*(-5*RS*pow(d,4)*pow(M_E,d*pow(RS,-1)) + pow(d,5)*pow(M_E,d*pow(RS,-1)) + 20*pow(d,3)*pow(M_E,d*pow(RS,-1))*pow(RS,2) - 
60*pow(d,2)*pow(M_E,d*pow(RS,-1))*pow(RS,3) + 120*d*pow(M_E,d*pow(RS,-1))*pow(RS,4) - 
pow(M_E,k1*pow(RS,-1))*(-5*RS*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RS,2) - 60*pow(k1,2)*pow(RS,3) + 120*k1*pow(RS,4) - 120*pow(RS,5)) - 
120*pow(M_E,d*pow(RS,-1))*pow(RS,5)) + RS*pow(b1,2)*(-5*RS*pow(d,4)*pow(M_E,d*pow(RS,-1)) + pow(d,5)*pow(M_E,d*pow(RS,-1)) + 
20*pow(d,3)*pow(M_E,d*pow(RS,-1))*pow(RS,2) - 60*pow(d,2)*pow(M_E,d*pow(RS,-1))*pow(RS,3) + 120*d*pow(M_E,d*pow(RS,-1))*pow(RS,4) - 
pow(M_E,k1*pow(RS,-1))*(-5*RS*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RS,2) - 60*pow(k1,2)*pow(RS,3) + 120*k1*pow(RS,4) - 120*pow(RS,5)) - 
120*pow(M_E,d*pow(RS,-1))*pow(RS,5)) - 2*a1*c1*(RS*pow(M_E,-(d*pow(RS,-1)))*
(-5*RS*pow(d,4) - pow(d,5) - 20*pow(d,3)*pow(RS,2) - 60*pow(d,2)*pow(RS,3) - 120*d*pow(RS,4) + 120*(-1 + pow(M_E,d*pow(RS,-1)))*pow(RS,5)) + 
RS*pow(M_E,-(k1*pow(RS,-1)))*(5*RS*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RS,2) + 60*pow(k1,2)*pow(RS,3) + 120*k1*pow(RS,4) - 
120*(-1 + pow(M_E,k1*pow(RS,-1)))*pow(RS,5))) - pow(b1,2)*
(RS*pow(M_E,-(d*pow(RS,-1)))*(-5*RS*pow(d,4) - pow(d,5) - 20*pow(d,3)*pow(RS,2) - 60*pow(d,2)*pow(RS,3) - 120*d*pow(RS,4) + 
120*(-1 + pow(M_E,d*pow(RS,-1)))*pow(RS,5)) + RS*pow(M_E,-(k1*pow(RS,-1)))*
(5*RS*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RS,2) + 60*pow(k1,2)*pow(RS,3) + 120*k1*pow(RS,4) - 120*(-1 + pow(M_E,k1*pow(RS,-1)))*pow(RS,5))) + 
2*a1*b1*RS*(-6*RS*pow(d,5)*pow(M_E,d*pow(RS,-1)) + pow(d,6)*pow(M_E,d*pow(RS,-1)) + 30*pow(d,4)*pow(M_E,d*pow(RS,-1))*pow(RS,2) - 
120*pow(d,3)*pow(M_E,d*pow(RS,-1))*pow(RS,3) + 360*pow(d,2)*pow(M_E,d*pow(RS,-1))*pow(RS,4) - 720*d*pow(M_E,d*pow(RS,-1))*pow(RS,5) + 
720*pow(M_E,d*pow(RS,-1))*pow(RS,6) - pow(M_E,k1*pow(RS,-1))*
(-6*RS*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RS,2) - 120*pow(k1,3)*pow(RS,3) + 360*pow(k1,2)*pow(RS,4) - 720*k1*pow(RS,5) + 720*pow(RS,6))) - 
2*a1*b1*(RS*pow(M_E,-(d*pow(RS,-1)))*(-6*RS*pow(d,5) - pow(d,6) - 30*pow(d,4)*pow(RS,2) - 120*pow(d,3)*pow(RS,3) - 360*pow(d,2)*pow(RS,4) - 720*d*pow(RS,5) + 
720*(-1 + pow(M_E,d*pow(RS,-1)))*pow(RS,6)) + RS*pow(M_E,-(k1*pow(RS,-1)))*
(6*RS*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RS,2) + 120*pow(k1,3)*pow(RS,3) + 360*pow(k1,2)*pow(RS,4) + 720*k1*pow(RS,5) - 
720*(-1 + pow(M_E,k1*pow(RS,-1)))*pow(RS,6))) + RS*pow(a1,2)*
(-7*RS*pow(d,6)*pow(M_E,d*pow(RS,-1)) + pow(d,7)*pow(M_E,d*pow(RS,-1)) + 42*pow(d,5)*pow(M_E,d*pow(RS,-1))*pow(RS,2) - 210*pow(d,4)*pow(M_E,d*pow(RS,-1))*pow(RS,3) + 
840*pow(d,3)*pow(M_E,d*pow(RS,-1))*pow(RS,4) - 2520*pow(d,2)*pow(M_E,d*pow(RS,-1))*pow(RS,5) + 5040*d*pow(M_E,d*pow(RS,-1))*pow(RS,6) - 
pow(M_E,k1*pow(RS,-1))*(-7*RS*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RS,2) - 210*pow(k1,4)*pow(RS,3) + 840*pow(k1,3)*pow(RS,4) - 2520*pow(k1,2)*pow(RS,5) + 
5040*k1*pow(RS,6) - 5040*pow(RS,7)) - 5040*pow(M_E,d*pow(RS,-1))*pow(RS,7)) - 
pow(a1,2)*(RS*pow(M_E,-(d*pow(RS,-1)))*(-7*RS*pow(d,6) - pow(d,7) - 42*pow(d,5)*pow(RS,2) - 210*pow(d,4)*pow(RS,3) - 840*pow(d,3)*pow(RS,4) - 
2520*pow(d,2)*pow(RS,5) - 5040*d*pow(RS,6) + 5040*(-1 + pow(M_E,d*pow(RS,-1)))*pow(RS,7)) + 
RS*pow(M_E,-(k1*pow(RS,-1)))*(7*RS*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RS,2) + 210*pow(k1,4)*pow(RS,3) + 840*pow(k1,3)*pow(RS,4) + 
2520*pow(k1,2)*pow(RS,5) + 5040*k1*pow(RS,6) - 5040*(-1 + pow(M_E,k1*pow(RS,-1)))*pow(RS,7)))))/2.;

return res;
}

double SDH_Integral_z_11_case_3( double AS, double RS, double k1, double k2, double a1, double b1, double c1, double d1, double d )
{
double res =-(AS*RS*pow(d,-2)*(d + RS + d*pow(M_E,2*d*pow(RS,-1)) - RS*pow(M_E,2*d*pow(RS,-1)))*pow(M_E,-((d + k1 + k2)*pow(RS,-1)))*
(pow(d1,2)*(-((k2 + RS)*pow(M_E,k1*pow(RS,-1))) + (k1 + RS)*pow(M_E,k2*pow(RS,-1))) + 5*RS*pow(b1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,4) + 
12*a1*b1*RS*pow(M_E,k2*pow(RS,-1))*pow(k1,5) + pow(b1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,5) + 2*a1*b1*pow(M_E,k2*pow(RS,-1))*pow(k1,6) + 
7*RS*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,6) + pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,7) - 5*RS*pow(b1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,4) - 
12*a1*b1*RS*pow(M_E,k1*pow(RS,-1))*pow(k2,5) - pow(b1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,5) - 2*a1*b1*pow(M_E,k1*pow(RS,-1))*pow(k2,6) - 
7*RS*pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,6) - pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,7) + 20*pow(b1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,3)*pow(RS,2) + 
60*a1*b1*pow(M_E,k2*pow(RS,-1))*pow(k1,4)*pow(RS,2) + 42*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,5)*pow(RS,2) - 
20*pow(b1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,3)*pow(RS,2) - 60*a1*b1*pow(M_E,k1*pow(RS,-1))*pow(k2,4)*pow(RS,2) - 
42*pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,5)*pow(RS,2) + 60*pow(b1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,2)*pow(RS,3) + 
240*a1*b1*pow(M_E,k2*pow(RS,-1))*pow(k1,3)*pow(RS,3) + 210*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,4)*pow(RS,3) - 
60*pow(b1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,2)*pow(RS,3) - 240*a1*b1*pow(M_E,k1*pow(RS,-1))*pow(k2,3)*pow(RS,3) - 
210*pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,4)*pow(RS,3) + pow(c1,2)*
(pow(M_E,k2*pow(RS,-1))*(3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) + 6*pow(RS,3)) - 
pow(M_E,k1*pow(RS,-1))*(3*RS*pow(k2,2) + pow(k2,3) + 6*k2*pow(RS,2) + 6*pow(RS,3))) - 120*k2*pow(b1,2)*pow(M_E,k1*pow(RS,-1))*pow(RS,4) + 
120*k1*pow(b1,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,4) + 720*a1*b1*pow(M_E,k2*pow(RS,-1))*pow(k1,2)*pow(RS,4) + 
840*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,3)*pow(RS,4) - 720*a1*b1*pow(M_E,k1*pow(RS,-1))*pow(k2,2)*pow(RS,4) - 
840*pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,3)*pow(RS,4) - 2*d1*
(-(c1*pow(M_E,k2*pow(RS,-1))*(2*k1*RS + pow(k1,2) + 2*pow(RS,2))) + c1*pow(M_E,k1*pow(RS,-1))*(2*k2*RS + pow(k2,2) + 2*pow(RS,2)) - 
b1*pow(M_E,k2*pow(RS,-1))*(3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) + 6*pow(RS,3)) + 
b1*pow(M_E,k1*pow(RS,-1))*(3*RS*pow(k2,2) + pow(k2,3) + 6*k2*pow(RS,2) + 6*pow(RS,3)) - 
a1*pow(M_E,k2*pow(RS,-1))*(4*RS*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RS,2) + 24*k1*pow(RS,3) + 24*pow(RS,4)) + 
a1*pow(M_E,k1*pow(RS,-1))*(4*RS*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RS,2) + 24*k2*pow(RS,3) + 24*pow(RS,4))) - 
1440*a1*b1*k2*pow(M_E,k1*pow(RS,-1))*pow(RS,5) - 120*pow(b1,2)*pow(M_E,k1*pow(RS,-1))*pow(RS,5) + 1440*a1*b1*k1*pow(M_E,k2*pow(RS,-1))*pow(RS,5) + 
120*pow(b1,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,5) + 2520*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,2)*pow(RS,5) - 
2520*pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,2)*pow(RS,5) - 2*c1*
(-(b1*pow(M_E,k2*pow(RS,-1))*(4*RS*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RS,2) + 24*k1*pow(RS,3) + 24*pow(RS,4))) + 
b1*pow(M_E,k1*pow(RS,-1))*(4*RS*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RS,2) + 24*k2*pow(RS,3) + 24*pow(RS,4)) - 
a1*pow(M_E,k2*pow(RS,-1))*(5*RS*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RS,2) + 60*pow(k1,2)*pow(RS,3) + 120*k1*pow(RS,4) + 120*pow(RS,5)) + 
a1*pow(M_E,k1*pow(RS,-1))*(5*RS*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RS,2) + 60*pow(k2,2)*pow(RS,3) + 120*k2*pow(RS,4) + 120*pow(RS,5))) - 
1440*a1*b1*pow(M_E,k1*pow(RS,-1))*pow(RS,6) - 5040*k2*pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(RS,6) + 1440*a1*b1*pow(M_E,k2*pow(RS,-1))*pow(RS,6) + 
5040*k1*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,6) - 5040*pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(RS,7) + 5040*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,7)))/2.;

return res;
}

/* < S | dzH | S > done */







/* < S | dzH | Z > */

double SDH_Integral_z_14_case_1( double ASP, double RSP, double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d )
{
double res =(ASP*RSP*pow(3,0.5)*pow(d,-3)*pow(M_E,-((d + k1 + k2)*pow(RSP,-1)))*(2*d*RSP + pow(d,2) + 2*pow(RSP,2))*
(4*b1*d2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,2) + 4*b1*d2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2) + b1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + 
5*b1*c2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + 5*a1*d2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) - b1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3) + 
5*b1*c2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3) + 5*a1*d2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3) + b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + 
a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + 6*b1*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + 6*a1*c2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) - 
b1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) - a1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) + 6*b1*b2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) + 
6*a1*c2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) + b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + 
7*a2*b1*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + 7*a1*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) - b1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) - 
a1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) + 7*a2*b1*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) + 7*a1*b2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) + 
a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) + a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) + 8*a1*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) - 
a2*b1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,6) - a1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,6) + 8*a1*a2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,6) + 
a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,7) - a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,7) - 4*b1*d2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,2) - 
4*b1*d2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2) - b1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - 5*b1*c2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - 
5*a1*d2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) + b1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3) - 5*b1*c2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3) - 
5*a1*d2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3) - b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 
6*b1*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 6*a1*c2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) + b1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) + 
a1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) - 6*b1*b2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) - 6*a1*c2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) - 
b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - 7*a2*b1*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - 
7*a1*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) + b1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) + a1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) - 
7*a2*b1*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) - 7*a1*b2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) - a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) - 
a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) - 8*a1*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) + a2*b1*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,6) + 
a1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,6) - 8*a1*a2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,6) - a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,7) + 
a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,7) - 8*b1*d2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) + 8*b1*d2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) - 
8*b1*d2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,2) + 8*b1*d2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,2) + 
15*b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 15*a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 
15*b1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 15*a1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 
24*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 24*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 
24*b1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 24*a1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 
35*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 35*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) - 
35*a2*b1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) - 35*a1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 
48*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) - 48*a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) - 
15*b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 15*a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 
15*b1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 15*a1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 
24*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 24*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 
24*b1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 24*a1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 
35*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) - 35*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) + 
35*a2*b1*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) + 35*a1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) - 
48*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) + 48*a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) - 8*b1*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
30*b1*c2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 30*a1*d2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) + 8*b1*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 
30*b1*c2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 30*a1*d2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 8*b1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,3) + 
30*b1*c2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,3) + 30*a1*d2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,3) - 
8*b1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,3) - 30*b1*c2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,3) - 
30*a1*d2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,3) + 72*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
72*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 72*b1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
72*a1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 140*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 
140*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 140*a2*b1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 
140*a1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 240*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) + 
240*a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) - 72*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 
72*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 72*b1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 
72*a1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 140*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) - 
140*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) - 140*a2*b1*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) - 
140*a1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) - 240*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) - 
240*a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) - 30*b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 30*a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 
144*b1*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 144*a1*c2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 30*b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 
30*a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 144*b1*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 144*a1*c2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 
30*b1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) - 30*a1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) - 
144*b1*b2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) - 144*a1*c2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) + 
30*b1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) + 30*a1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) + 
144*b1*b2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) + 144*a1*c2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) + 
420*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 420*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) - 
420*a2*b1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) - 420*a1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 
960*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) - 960*a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) - 
420*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) - 420*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 
420*a2*b1*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 420*a1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) - 
960*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) + 960*a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) + 
d1*(d2*(-((k2 + 2*RSP)*pow(M_E,k1*pow(RSP,-1))) + (k1 + 2*RSP)*pow(M_E,k2*pow(RSP,-1)) - (k1 - 2*RSP)*pow(M_E,(2*k1 + k2)*pow(RSP,-1)) + 
(k2 - 2*RSP)*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))) + 4*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,2) + 4*b2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2) + 
b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + 5*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) - b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3) + 
5*a2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3) + a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) - a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) - 
4*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,2) - 4*b2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2) - b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - 
5*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) + b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3) - 5*a2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3) - 
a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) + a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) - 8*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) + 
8*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) - 8*b2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,2) + 8*b2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,2) + 
15*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 15*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 
15*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 15*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 
c2*(-(pow(M_E,(2*k1 + k2)*pow(RSP,-1))*(-3*k1*RSP + pow(k1,2) + 3*pow(RSP,2))) + pow(M_E,k2*pow(RSP,-1))*(3*k1*RSP + pow(k1,2) + 3*pow(RSP,2)) + 
pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*(-3*k2*RSP + pow(k2,2) + 3*pow(RSP,2)) - pow(M_E,k1*pow(RSP,-1))*(3*k2*RSP + pow(k2,2) + 3*pow(RSP,2))) - 
8*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 30*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) + 8*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 
30*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 8*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,3) + 30*a2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,3) - 
8*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,3) - 30*a2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,3) - 30*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 
30*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 30*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) + 30*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4)) - 
144*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 144*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 840*a2*b1*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 
840*a1*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 144*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 144*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 
840*a2*b1*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 840*a1*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 144*b1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,5) + 
144*a1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,5) + 840*a2*b1*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,5) + 
840*a1*b2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,5) - 144*b1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,5) - 
144*a1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,5) - 840*a2*b1*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,5) - 
840*a1*b2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,5) + 2880*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) + 
2880*a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) - 2880*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) - 
2880*a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) + 
c1*(5*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + 5*b2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3) + b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + 
6*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) - b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) + 6*a2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) + 
a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) - a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) - 5*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - 
5*b2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3) - b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 6*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) + 
b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) - 6*a2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) - a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) + 
a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) + 15*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 15*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 
24*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 24*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 
15*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 15*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 
24*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 24*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 
d2*(-(pow(M_E,(2*k1 + k2)*pow(RSP,-1))*(-3*k1*RSP + pow(k1,2) + 3*pow(RSP,2))) + pow(M_E,k2*pow(RSP,-1))*(3*k1*RSP + pow(k1,2) + 3*pow(RSP,2)) + 
pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*(-3*k2*RSP + pow(k2,2) + 3*pow(RSP,2)) - pow(M_E,k1*pow(RSP,-1))*(3*k2*RSP + pow(k2,2) + 3*pow(RSP,2))) - 
30*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) + 30*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 30*b2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,3) - 
30*b2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,3) + 72*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
72*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) - 72*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 
72*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 
c2*(-(pow(M_E,(2*k1 + k2)*pow(RSP,-1))*(-4*RSP*pow(k1,2) + pow(k1,3) + 8*k1*pow(RSP,2) - 8*pow(RSP,3))) + 
pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*(-4*RSP*pow(k2,2) + pow(k2,3) + 8*k2*pow(RSP,2) - 8*pow(RSP,3)) + 
pow(M_E,k2*pow(RSP,-1))*(4*RSP*pow(k1,2) + pow(k1,3) + 8*k1*pow(RSP,2) + 8*pow(RSP,3)) - 
pow(M_E,k1*pow(RSP,-1))*(4*RSP*pow(k2,2) + pow(k2,3) + 8*k2*pow(RSP,2) + 8*pow(RSP,3))) - 30*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 
144*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 30*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 144*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 
30*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) - 144*a2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) + 30*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) + 
144*a2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) - 144*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 144*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 
144*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,5) - 144*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,5)) - 840*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) - 
840*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) - 5760*a1*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) + 840*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) + 
840*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) + 5760*a1*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) - 840*a2*b1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,6) - 
840*a1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,6) - 5760*a1*a2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,6) + 
840*a2*b1*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,6) + 840*a1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,6) + 
5760*a1*a2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,6) - 5760*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,7) + 5760*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,7) + 
5760*a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,7) - 5760*a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,7)))/2.;

return res;
}


double SDH_Integral_z_14_case_2_sub_1( double ASP, double RSP, double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d )
{
double res =-(ASP*RSP*pow(3,0.5)*pow(d,-3)*pow(M_E,-((2*d + k2)*pow(RSP,-1)))*(pow(d,2)*(-1 + pow(M_E,2*d*pow(RSP,-1))) - 2*d*RSP*(1 + pow(M_E,2*d*pow(RSP,-1))) + 
2*(-1 + pow(M_E,2*d*pow(RSP,-1)))*pow(RSP,2))*(-(d1*d2*k2*pow(M_E,d*pow(RSP,-1))) - 2*d1*d2*RSP*pow(M_E,d*pow(RSP,-1)) - 3*c2*d1*k2*RSP*pow(M_E,d*pow(RSP,-1)) + 
d*d1*d2*pow(M_E,k2*pow(RSP,-1)) + 3*c2*d*d1*RSP*pow(M_E,k2*pow(RSP,-1)) + 2*d1*d2*RSP*pow(M_E,k2*pow(RSP,-1)) + c2*d1*pow(d,2)*pow(M_E,k2*pow(RSP,-1)) + 
4*b2*d1*RSP*pow(d,2)*pow(M_E,k2*pow(RSP,-1)) + b2*d1*pow(d,3)*pow(M_E,k2*pow(RSP,-1)) + 5*a2*d1*RSP*pow(d,3)*pow(M_E,k2*pow(RSP,-1)) + 
5*a1*d2*RSP*pow(d,3)*pow(M_E,k2*pow(RSP,-1)) + a2*d1*pow(d,4)*pow(M_E,k2*pow(RSP,-1)) + a1*d2*pow(d,4)*pow(M_E,k2*pow(RSP,-1)) + 
6*a1*c2*RSP*pow(d,4)*pow(M_E,k2*pow(RSP,-1)) + a1*c2*pow(d,5)*pow(M_E,k2*pow(RSP,-1)) + 7*a1*b2*RSP*pow(d,5)*pow(M_E,k2*pow(RSP,-1)) + 
a1*b2*pow(d,6)*pow(M_E,k2*pow(RSP,-1)) + 8*a1*a2*RSP*pow(d,6)*pow(M_E,k2*pow(RSP,-1)) + a1*a2*pow(d,7)*pow(M_E,k2*pow(RSP,-1)) - 
c2*d1*pow(M_E,d*pow(RSP,-1))*pow(k2,2) - 4*b2*d1*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,2) - b2*d1*pow(M_E,d*pow(RSP,-1))*pow(k2,3) - 
5*a2*d1*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,3) - 5*a1*d2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,3) - a2*d1*pow(M_E,d*pow(RSP,-1))*pow(k2,4) - 
a1*d2*pow(M_E,d*pow(RSP,-1))*pow(k2,4) - 6*a1*c2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,4) - a1*c2*pow(M_E,d*pow(RSP,-1))*pow(k2,5) - 
7*a1*b2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,5) - a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k2,6) - 8*a1*a2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,6) - 
a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,7) - 3*c2*d1*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) - 8*b2*d1*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) + 
3*c2*d1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 8*b2*d*d1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 15*a2*d1*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 
15*a1*d2*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 24*a1*c2*pow(d,3)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 35*a1*b2*pow(d,4)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 
48*a1*a2*pow(d,5)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) - 15*a2*d1*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 15*a1*d2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 
24*a1*c2*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 35*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) - 48*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) - 
8*b2*d1*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) - 30*a2*d1*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) - 30*a1*d2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 
8*b2*d1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 30*a2*d*d1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 30*a1*d*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 
72*a1*c2*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 140*a1*b2*pow(d,3)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 
240*a1*a2*pow(d,4)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 72*a1*c2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 
140*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) - 240*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) - 30*a2*d1*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) - 
30*a1*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) - 144*a1*c2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 144*a1*c2*d*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 
30*a2*d1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 30*a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 420*a1*b2*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 
960*a1*a2*pow(d,3)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 420*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) - 
960*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) - 144*a1*c2*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) - 840*a1*b2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 
144*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 840*a1*b2*d*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 2880*a1*a2*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 
2880*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) + c1*(-3*d2*k2*RSP*pow(M_E,d*pow(RSP,-1)) + 3*d*d2*RSP*pow(M_E,k2*pow(RSP,-1)) + 
d2*pow(d,2)*pow(M_E,k2*pow(RSP,-1)) + 5*b2*RSP*pow(d,3)*pow(M_E,k2*pow(RSP,-1)) + b2*pow(d,4)*pow(M_E,k2*pow(RSP,-1)) + 6*a2*RSP*pow(d,4)*pow(M_E,k2*pow(RSP,-1)) + 
a2*pow(d,5)*pow(M_E,k2*pow(RSP,-1)) - d2*pow(M_E,d*pow(RSP,-1))*pow(k2,2) - 5*b2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,3) - b2*pow(M_E,d*pow(RSP,-1))*pow(k2,4) - 
6*a2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,4) - a2*pow(M_E,d*pow(RSP,-1))*pow(k2,5) - 3*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) + 
3*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 15*b2*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 24*a2*pow(d,3)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) - 
15*b2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 24*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 30*b2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 
30*b2*d*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 72*a2*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 72*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 
c2*(4*RSP*pow(d,2)*pow(M_E,k2*pow(RSP,-1)) + pow(d,3)*pow(M_E,k2*pow(RSP,-1)) + 8*d*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 8*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 
pow(M_E,d*pow(RSP,-1))*(4*RSP*pow(k2,2) + pow(k2,3) + 8*k2*pow(RSP,2) + 8*pow(RSP,3))) - 30*b2*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) - 
144*a2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 30*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 144*a2*d*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 
144*a2*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 144*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5)) - 840*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) - 
5760*a1*a2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) + 840*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) + 5760*a1*a2*d*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) + 
b1*(4*d2*RSP*pow(d,2)*pow(M_E,k2*pow(RSP,-1)) + d2*pow(d,3)*pow(M_E,k2*pow(RSP,-1)) + 6*b2*RSP*pow(d,4)*pow(M_E,k2*pow(RSP,-1)) + 
b2*pow(d,5)*pow(M_E,k2*pow(RSP,-1)) + 7*a2*RSP*pow(d,5)*pow(M_E,k2*pow(RSP,-1)) + a2*pow(d,6)*pow(M_E,k2*pow(RSP,-1)) - 4*d2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,2) - 
d2*pow(M_E,d*pow(RSP,-1))*pow(k2,3) - 6*b2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,4) - b2*pow(M_E,d*pow(RSP,-1))*pow(k2,5) - 7*a2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,5) - 
a2*pow(M_E,d*pow(RSP,-1))*pow(k2,6) - 8*d2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) + 8*d*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 
24*b2*pow(d,3)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 35*a2*pow(d,4)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) - 24*b2*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 
35*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) - 8*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 8*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 
72*b2*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 140*a2*pow(d,3)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 72*b2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 
140*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) - 144*b2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 144*b2*d*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 
420*a2*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 420*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 
c2*(5*RSP*pow(d,3)*pow(M_E,k2*pow(RSP,-1)) + pow(d,4)*pow(M_E,k2*pow(RSP,-1)) + 15*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 
30*d*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 30*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 
pow(M_E,d*pow(RSP,-1))*(5*RSP*pow(k2,3) + pow(k2,4) + 15*pow(k2,2)*pow(RSP,2) + 30*k2*pow(RSP,3) + 30*pow(RSP,4))) - 
144*b2*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) - 840*a2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 144*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 
840*a2*d*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 840*a2*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) + 840*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6)) - 
5760*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(RSP,7) + 5760*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,7)))/2.;

return res;
}


double SDH_Integral_z_14_case_2_sub_2( double ASP, double RSP, double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d )
{
double res =(ASP*RSP*pow(3,0.5)*pow(d,-3)*pow(M_E,-((2*d + k1)*pow(RSP,-1)))*(2*d*RSP + pow(d,2) + 2*pow(RSP,2))*
(d1*d2*k1*pow(M_E,d*pow(RSP,-1)) + 2*d1*d2*RSP*pow(M_E,d*pow(RSP,-1)) + 3*c2*d1*k1*RSP*pow(M_E,d*pow(RSP,-1)) - d*d1*d2*pow(M_E,k1*pow(RSP,-1)) - 
3*c2*d*d1*RSP*pow(M_E,k1*pow(RSP,-1)) - 2*d1*d2*RSP*pow(M_E,k1*pow(RSP,-1)) - c2*d1*pow(d,2)*pow(M_E,k1*pow(RSP,-1)) - 4*b2*d1*RSP*pow(d,2)*pow(M_E,k1*pow(RSP,-1)) - 
b2*d1*pow(d,3)*pow(M_E,k1*pow(RSP,-1)) - 5*a2*d1*RSP*pow(d,3)*pow(M_E,k1*pow(RSP,-1)) - 5*a1*d2*RSP*pow(d,3)*pow(M_E,k1*pow(RSP,-1)) - 
a2*d1*pow(d,4)*pow(M_E,k1*pow(RSP,-1)) - a1*d2*pow(d,4)*pow(M_E,k1*pow(RSP,-1)) - 6*a1*c2*RSP*pow(d,4)*pow(M_E,k1*pow(RSP,-1)) - 
a1*c2*pow(d,5)*pow(M_E,k1*pow(RSP,-1)) - 7*a1*b2*RSP*pow(d,5)*pow(M_E,k1*pow(RSP,-1)) - a1*b2*pow(d,6)*pow(M_E,k1*pow(RSP,-1)) - 
8*a1*a2*RSP*pow(d,6)*pow(M_E,k1*pow(RSP,-1)) - a1*a2*pow(d,7)*pow(M_E,k1*pow(RSP,-1)) + d*d1*d2*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 
3*c2*d*d1*RSP*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 2*d1*d2*RSP*pow(M_E,(2*d + k1)*pow(RSP,-1)) + c2*d1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 
4*b2*d1*RSP*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + b2*d1*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 5*a2*d1*RSP*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 
5*a1*d2*RSP*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + a2*d1*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + a1*d2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 
6*a1*c2*RSP*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + a1*c2*pow(d,5)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 7*a1*b2*RSP*pow(d,5)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + 
a1*b2*pow(d,6)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 8*a1*a2*RSP*pow(d,6)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + a1*a2*pow(d,7)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 
d1*d2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1)) + 2*d1*d2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1)) + 3*c2*d1*k1*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1)) + 
c2*d1*pow(M_E,d*pow(RSP,-1))*pow(k1,2) + 4*b2*d1*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,2) - c2*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2) + 
4*b2*d1*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2) + b2*d1*pow(M_E,d*pow(RSP,-1))*pow(k1,3) + 5*a2*d1*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,3) + 
5*a1*d2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,3) - b2*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3) + 5*a2*d1*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3) + 
5*a1*d2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3) + a2*d1*pow(M_E,d*pow(RSP,-1))*pow(k1,4) + a1*d2*pow(M_E,d*pow(RSP,-1))*pow(k1,4) + 
6*a1*c2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,4) - a2*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4) - a1*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4) + 
6*a1*c2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4) + a1*c2*pow(M_E,d*pow(RSP,-1))*pow(k1,5) + 7*a1*b2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,5) - 
a1*c2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5) + 7*a1*b2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5) + a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k1,6) + 
8*a1*a2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,6) - a1*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,6) + 8*a1*a2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,6) + 
a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,7) - a1*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,7) + 3*c2*d1*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) + 
8*b2*d1*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) - 3*c2*d1*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 8*b2*d*d1*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 
15*a2*d1*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 15*a1*d2*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 24*a1*c2*pow(d,3)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 
35*a1*b2*pow(d,4)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 48*a1*a2*pow(d,5)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) + 3*c2*d1*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) + 
8*b2*d*d1*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) + 15*a2*d1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) + 
15*a1*d2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) + 24*a1*c2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) + 
35*a1*b2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) + 48*a1*a2*pow(d,5)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) - 
3*c2*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,2) - 8*b2*d1*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,2) + 15*a2*d1*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 
15*a1*d2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 15*a2*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 
15*a1*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 24*a1*c2*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 
24*a1*c2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 35*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) - 
35*a1*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 48*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) - 
48*a1*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) + 8*b2*d1*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 30*a2*d1*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 
30*a1*d2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) - 8*b2*d1*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 30*a2*d*d1*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
30*a1*d*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 72*a1*c2*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 140*a1*b2*pow(d,3)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
240*a1*a2*pow(d,4)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 8*b2*d1*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) - 30*a2*d*d1*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) - 
30*a1*d*d2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) - 72*a1*c2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) - 
140*a1*b2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) - 240*a1*a2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) + 
8*b2*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,3) + 30*a2*d1*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,3) + 
30*a1*d2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,3) + 72*a1*c2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
72*a1*c2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 140*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 
140*a1*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 240*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) + 
240*a1*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) + 30*a2*d1*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 30*a1*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 
144*a1*c2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) - 144*a1*c2*d*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 30*a2*d1*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 
30*a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 420*a1*b2*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 960*a1*a2*pow(d,3)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 
144*a1*c2*d*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) + 30*a2*d1*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) + 30*a1*d2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) + 
420*a1*b2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) + 960*a1*a2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) - 
30*a2*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,4) - 30*a1*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,4) - 144*a1*c2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,4) + 
420*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) - 420*a1*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 
960*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) - 960*a1*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) + 
144*a1*c2*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 840*a1*b2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) - 144*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 
840*a1*b2*d*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 2880*a1*a2*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 144*a1*c2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,5) - 
840*a1*b2*d*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,5) - 2880*a1*a2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,5) + 
144*a1*c2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,5) + 840*a1*b2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,5) + 
2880*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) + 2880*a1*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) + 
c1*(3*d2*k1*RSP*pow(M_E,d*pow(RSP,-1)) - 3*d*d2*RSP*pow(M_E,k1*pow(RSP,-1)) - d2*pow(d,2)*pow(M_E,k1*pow(RSP,-1)) - 5*b2*RSP*pow(d,3)*pow(M_E,k1*pow(RSP,-1)) - 
b2*pow(d,4)*pow(M_E,k1*pow(RSP,-1)) - 6*a2*RSP*pow(d,4)*pow(M_E,k1*pow(RSP,-1)) - a2*pow(d,5)*pow(M_E,k1*pow(RSP,-1)) - 3*d*d2*RSP*pow(M_E,(2*d + k1)*pow(RSP,-1)) + 
d2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 5*b2*RSP*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + b2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 
6*a2*RSP*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + a2*pow(d,5)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + 3*d2*k1*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1)) + 
d2*pow(M_E,d*pow(RSP,-1))*pow(k1,2) - d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2) + 5*b2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,3) + 
5*b2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3) + b2*pow(M_E,d*pow(RSP,-1))*pow(k1,4) + 6*a2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,4) - 
b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4) + 6*a2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4) + a2*pow(M_E,d*pow(RSP,-1))*pow(k1,5) - 
a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5) + 3*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) - 3*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 
15*b2*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 24*a2*pow(d,3)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) + 3*d2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) + 
15*b2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) + 24*a2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) - 
3*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,2) + 15*b2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 15*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 
24*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 24*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 30*b2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) - 
30*b2*d*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 72*a2*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 30*b2*d*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) - 
72*a2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) + 30*b2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,3) + 
72*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 72*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
c2*(pow(d,3)*(-1 + pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1)) - 4*RSP*pow(d,2)*(1 + pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1)) + 
8*d*(-1 + pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 8*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 8*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) + 
pow(M_E,(d + 2*k1)*pow(RSP,-1))*(4*RSP*pow(k1,2) - pow(k1,3) - 8*k1*pow(RSP,2) + 8*pow(RSP,3)) + 
pow(M_E,d*pow(RSP,-1))*(4*RSP*pow(k1,2) + pow(k1,3) + 8*k1*pow(RSP,2) + 8*pow(RSP,3))) + 30*b2*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 
144*a2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) - 30*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 144*a2*d*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 
30*b2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) + 144*a2*d*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) - 30*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,4) - 
144*a2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,4) + 144*a2*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) - 144*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 
144*a2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,5) + 144*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,5)) + 840*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) + 
5760*a1*a2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) - 840*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) - 5760*a1*a2*d*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) + 
840*a1*b2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,6) + 5760*a1*a2*d*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,6) - 
840*a1*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,6) - 5760*a1*a2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,6) + 
b1*(-4*d2*RSP*pow(d,2)*pow(M_E,k1*pow(RSP,-1)) - d2*pow(d,3)*pow(M_E,k1*pow(RSP,-1)) - 6*b2*RSP*pow(d,4)*pow(M_E,k1*pow(RSP,-1)) - 
b2*pow(d,5)*pow(M_E,k1*pow(RSP,-1)) - 7*a2*RSP*pow(d,5)*pow(M_E,k1*pow(RSP,-1)) - a2*pow(d,6)*pow(M_E,k1*pow(RSP,-1)) - 
4*d2*RSP*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + d2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 6*b2*RSP*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + 
b2*pow(d,5)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 7*a2*RSP*pow(d,5)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + a2*pow(d,6)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + 
4*d2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,2) + 4*d2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2) + d2*pow(M_E,d*pow(RSP,-1))*pow(k1,3) - 
d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3) + 6*b2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,4) + 6*b2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4) + 
b2*pow(M_E,d*pow(RSP,-1))*pow(k1,5) + 7*a2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,5) - b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5) + 
7*a2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5) + a2*pow(M_E,d*pow(RSP,-1))*pow(k1,6) - a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,6) + 
8*d2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) - 8*d*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 24*b2*pow(d,3)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 
35*a2*pow(d,4)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) + 8*d*d2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) + 
24*b2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) + 35*a2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) - 
8*d2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,2) + 24*b2*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 
24*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 35*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) - 
35*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 8*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) - 8*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
72*b2*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 140*a2*pow(d,3)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 8*d2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) - 
72*b2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) - 140*a2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) + 
8*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,3) + 72*b2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 72*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
140*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 140*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 144*b2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) - 
144*b2*d*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 420*a2*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 144*b2*d*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) + 
420*a2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) - 144*b2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,4) + 
420*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) - 420*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 
c2*(pow(d,4)*(-1 + pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1)) - 5*RSP*pow(d,3)*(1 + pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1)) + 
15*pow(d,2)*(-1 + pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 30*d*(1 + pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
30*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 30*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) - 
pow(M_E,(d + 2*k1)*pow(RSP,-1))*(-5*RSP*pow(k1,3) + pow(k1,4) + 15*pow(k1,2)*pow(RSP,2) - 30*k1*pow(RSP,3) + 30*pow(RSP,4)) + 
pow(M_E,d*pow(RSP,-1))*(5*RSP*pow(k1,3) + pow(k1,4) + 15*pow(k1,2)*pow(RSP,2) + 30*k1*pow(RSP,3) + 30*pow(RSP,4))) + 
144*b2*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 840*a2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) - 144*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 
840*a2*d*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 144*b2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,5) - 840*a2*d*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,5) + 
144*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,5) + 840*a2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,5) + 840*a2*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) - 
840*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) + 840*a2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,6) - 840*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,6)) + 
5760*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(RSP,7) - 5760*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,7) - 5760*a1*a2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,7) + 
5760*a1*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,7)))/2.;

return res;
}


double SDH_Integral_z_14_case_3( double ASP, double RSP, double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d )
{
double res =-(ASP*RSP*pow(3,0.5)*pow(d,-3)*pow(M_E,-((d + k1 + k2)*pow(RSP,-1)))*
(pow(d,2)*(-1 + pow(M_E,2*d*pow(RSP,-1))) - 2*d*RSP*(1 + pow(M_E,2*d*pow(RSP,-1))) + 2*(-1 + pow(M_E,2*d*pow(RSP,-1)))*pow(RSP,2))*
(4*b1*d2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,2) + b1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + 5*b1*c2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + 
5*a1*d2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + 
6*b1*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + 6*a1*c2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + 
a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + 7*a2*b1*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + 7*a1*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + 
a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) + a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) + 8*a1*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) + 
a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,7) - 4*b1*d2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,2) - b1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - 
5*b1*c2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - 5*a1*d2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 
a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 6*b1*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 6*a1*c2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 
b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - 7*a2*b1*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - 
7*a1*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) - a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) - 
8*a1*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) - a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,7) - 8*b1*d2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) + 
8*b1*d2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 15*b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 15*a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 
24*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 24*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 
35*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 35*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 
48*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) - 15*b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 
15*a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 24*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 
24*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 35*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) - 
35*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) - 48*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) - 8*b1*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
30*b1*c2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 30*a1*d2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) + 8*b1*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 
30*b1*c2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 30*a1*d2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 72*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
72*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 140*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 
140*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 240*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) - 
72*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 72*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 
140*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) - 140*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) - 
240*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) - 30*b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 30*a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 
144*b1*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 144*a1*c2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 30*b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 
30*a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 144*b1*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 144*a1*c2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 
420*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 420*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 
960*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) - 420*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) - 
420*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) - 960*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) + 
d1*(-(d2*(k2 + 2*RSP)*pow(M_E,k1*pow(RSP,-1))) + d2*(k1 + 2*RSP)*pow(M_E,k2*pow(RSP,-1)) + 4*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,2) + 
b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + 5*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) - 
4*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,2) - b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - 5*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - 
a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 8*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) + 8*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 
15*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 15*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 
c2*pow(M_E,k2*pow(RSP,-1))*(3*k1*RSP + pow(k1,2) + 3*pow(RSP,2)) - c2*pow(M_E,k1*pow(RSP,-1))*(3*k2*RSP + pow(k2,2) + 3*pow(RSP,2)) - 
8*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 30*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) + 8*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 
30*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 30*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 30*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4)) - 
144*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 144*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 840*a2*b1*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 
840*a1*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 144*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 144*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 
840*a2*b1*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 840*a1*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 2880*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) - 
2880*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) + c1*
(5*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + 6*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + 
a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) - 5*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 
6*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) + 15*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 
24*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 15*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 24*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 
d2*pow(M_E,k2*pow(RSP,-1))*(3*k1*RSP + pow(k1,2) + 3*pow(RSP,2)) - d2*pow(M_E,k1*pow(RSP,-1))*(3*k2*RSP + pow(k2,2) + 3*pow(RSP,2)) - 
30*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) + 30*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 72*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) - 
72*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + c2*pow(M_E,k2*pow(RSP,-1))*(4*RSP*pow(k1,2) + pow(k1,3) + 8*k1*pow(RSP,2) + 8*pow(RSP,3)) - 
c2*pow(M_E,k1*pow(RSP,-1))*(4*RSP*pow(k2,2) + pow(k2,3) + 8*k2*pow(RSP,2) + 8*pow(RSP,3)) - 30*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 
144*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 30*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 144*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 
144*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 144*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5)) - 840*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) - 
840*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) - 5760*a1*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) + 840*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) + 
840*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) + 5760*a1*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) - 5760*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,7) + 
5760*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,7)))/2.;

return res;
}

/* < S | dzH | Z > done */




/* < X or Y | dzH | X or Y > */
double SDH_Integral_z_2233_case_1( double AP, double RP, double k1, double k2, double a2, double b2, double c2, double d2, double d )
{
double res =(3*AP*pow(d,-4)*pow(M_E,-((d + k1 + k2)*pow(RP,-1)))*pow(RP,2)*(-12*c2*d2*k2*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1)) - 2*c2*d2*k2*pow(d,3)*pow(M_E,k1*pow(RP,-1)) - 
4*c2*d2*RP*pow(d,3)*pow(M_E,k1*pow(RP,-1)) - 6*b2*d2*k2*RP*pow(d,3)*pow(M_E,k1*pow(RP,-1)) - 3*k2*RP*pow(c2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1)) - 
3*d*k2*RP*pow(d2,2)*pow(M_E,k1*pow(RP,-1)) - k2*pow(d,2)*pow(d2,2)*pow(M_E,k1*pow(RP,-1)) - 5*RP*pow(d,2)*pow(d2,2)*pow(M_E,k1*pow(RP,-1)) - 
pow(d,3)*pow(d2,2)*pow(M_E,k1*pow(RP,-1)) + 12*c2*d2*k1*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1)) + 2*c2*d2*k1*pow(d,3)*pow(M_E,k2*pow(RP,-1)) + 
4*c2*d2*RP*pow(d,3)*pow(M_E,k2*pow(RP,-1)) + 6*b2*d2*k1*RP*pow(d,3)*pow(M_E,k2*pow(RP,-1)) + 3*k1*RP*pow(c2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1)) + 
3*d*k1*RP*pow(d2,2)*pow(M_E,k2*pow(RP,-1)) + k1*pow(d,2)*pow(d2,2)*pow(M_E,k2*pow(RP,-1)) + 5*RP*pow(d,2)*pow(d2,2)*pow(M_E,k2*pow(RP,-1)) + 
pow(d,3)*pow(d2,2)*pow(M_E,k2*pow(RP,-1)) + 4*RP*gsl_sf_expint_Ei(-(k2*pow(RP,-1)))*pow(d,2)*pow(d2,2)*pow(M_E,(k1 + k2)*pow(RP,-1)) - 
4*RP*gsl_sf_expint_Ei(k2*pow(RP,-1))*pow(d,2)*pow(d2,2)*pow(M_E,(k1 + k2)*pow(RP,-1)) + 
gsl_sf_expint_Ei(-(k2*pow(RP,-1)))*pow(d,3)*pow(d2,2)*pow(M_E,(k1 + k2)*pow(RP,-1)) - gsl_sf_expint_Ei(k2*pow(RP,-1))*pow(d,3)*pow(d2,2)*pow(M_E,(k1 + k2)*pow(RP,-1)) - 
12*c2*d2*k1*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) - 2*c2*d2*k1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) + 4*c2*d2*RP*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) + 
6*b2*d2*k1*RP*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) + 3*k1*RP*pow(c2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) + 
3*d*k1*RP*pow(d2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) + k1*pow(d,2)*pow(d2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) - 
5*RP*pow(d,2)*pow(d2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) - pow(d,3)*pow(d2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) + 
12*c2*d2*k2*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) + 2*c2*d2*k2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) - 4*c2*d2*RP*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) - 
6*b2*d2*k2*RP*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) - 3*k2*RP*pow(c2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) - 
3*d*k2*RP*pow(d2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) - k2*pow(d,2)*pow(d2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) + 
5*RP*pow(d,2)*pow(d2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) + pow(d,3)*pow(d2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) + 6*c2*d*d2*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 
2*c2*d2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 14*b2*d2*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 
7*RP*pow(c2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 2*b2*d2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 8*b2*c2*RP*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 
8*a2*d2*RP*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + pow(c2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 6*c2*d*d2*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) + 
2*c2*d2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) - 14*b2*d2*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) - 
7*RP*pow(c2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) - 2*b2*d2*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) + 
8*b2*c2*RP*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) + 8*a2*d2*RP*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) - 
pow(c2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) + 6*b2*d*d2*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 3*d*RP*pow(c2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 
2*b2*d2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 16*b2*c2*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 16*a2*d2*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 
pow(c2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 2*b2*c2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 2*a2*d2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 
10*a2*c2*RP*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 5*RP*pow(b2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 
6*b2*d*d2*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + 3*d*RP*pow(c2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + 
2*b2*d2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) - 16*b2*c2*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) - 
16*a2*d2*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + pow(c2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) - 
2*b2*c2*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) - 2*a2*d2*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + 
10*a2*c2*RP*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + 5*RP*pow(b2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + 
6*b2*c2*d*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 6*a2*d*d2*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 2*b2*c2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 
2*a2*d2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 18*a2*c2*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 
9*RP*pow(b2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 2*a2*c2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 
12*a2*b2*RP*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + pow(b2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 
6*b2*c2*d*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) + 6*a2*d*d2*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) + 
2*b2*c2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) + 2*a2*d2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) - 
18*a2*c2*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) - 9*RP*pow(b2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) - 
2*a2*c2*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) + 12*a2*b2*RP*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) - 
pow(b2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) + 6*a2*c2*d*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 3*d*RP*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 
2*a2*c2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 20*a2*b2*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + pow(b2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 
2*a2*b2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 7*RP*pow(a2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 
6*a2*c2*d*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) + 3*d*RP*pow(b2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) + 
2*a2*c2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) - 20*a2*b2*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) + 
pow(b2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) - 2*a2*b2*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) + 
7*RP*pow(a2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) + 6*a2*b2*d*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,6) + 
2*a2*b2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,6) + 11*RP*pow(a2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,6) + pow(a2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,6) + 
6*a2*b2*d*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6) + 2*a2*b2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6) - 
11*RP*pow(a2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6) - pow(a2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6) + 
3*d*RP*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,7) + pow(a2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,7) + 
3*d*RP*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,7) + pow(a2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,7) - 
6*c2*d*d2*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 2*c2*d2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 14*b2*d2*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 
7*RP*pow(c2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 2*b2*d2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 8*b2*c2*RP*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 
8*a2*d2*RP*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - pow(c2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 6*c2*d*d2*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) - 
2*c2*d2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) + 14*b2*d2*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) + 
7*RP*pow(c2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) + 2*b2*d2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) - 
8*b2*c2*RP*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) - 8*a2*d2*RP*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) + 
pow(c2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) - 6*b2*d*d2*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 3*d*RP*pow(c2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 
2*b2*d2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 16*b2*c2*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 16*a2*d2*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 
pow(c2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 2*b2*c2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 2*a2*d2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 
10*a2*c2*RP*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 5*RP*pow(b2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 
6*b2*d*d2*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - 3*d*RP*pow(c2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - 
2*b2*d2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) + 16*b2*c2*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) + 
16*a2*d2*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - pow(c2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) + 
2*b2*c2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) + 2*a2*d2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - 
10*a2*c2*RP*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - 5*RP*pow(b2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - 
6*b2*c2*d*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 6*a2*d*d2*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 2*b2*c2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 
2*a2*d2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 18*a2*c2*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 
9*RP*pow(b2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 2*a2*c2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 
12*a2*b2*RP*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - pow(b2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 
6*b2*c2*d*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) - 6*a2*d*d2*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) - 
2*b2*c2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) - 2*a2*d2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) + 
18*a2*c2*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) + 9*RP*pow(b2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) + 
2*a2*c2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) - 12*a2*b2*RP*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) + 
pow(b2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) - 6*a2*c2*d*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 3*d*RP*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 
2*a2*c2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 20*a2*b2*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - pow(b2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 
2*a2*b2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 7*RP*pow(a2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 
6*a2*c2*d*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) - 3*d*RP*pow(b2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) - 
2*a2*c2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) + 20*a2*b2*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) - 
pow(b2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) + 2*a2*b2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) - 
7*RP*pow(a2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) - 6*a2*b2*d*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,6) - 
2*a2*b2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,6) - 11*RP*pow(a2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,6) - pow(a2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,6) - 
6*a2*b2*d*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6) - 2*a2*b2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6) + 
11*RP*pow(a2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6) + pow(a2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6) - 
3*d*RP*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,7) - pow(a2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,7) - 
3*d*RP*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,7) - pow(a2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,7) - 
30*c2*d*d2*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 20*c2*d2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 36*b2*d2*k2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 
18*k2*pow(c2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 6*b2*d2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 
16*b2*c2*k2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 16*a2*d2*k2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 
3*pow(c2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 12*d*pow(d2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 3*k2*pow(d2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
30*c2*d*d2*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 20*c2*d2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 36*b2*d2*k1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 
18*k1*pow(c2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 6*b2*d2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 
16*b2*c2*k1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 16*a2*d2*k1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 
3*pow(c2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 12*d*pow(d2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 3*k1*pow(d2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 
9*d*gsl_sf_expint_Ei(-(k2*pow(RP,-1)))*pow(d2,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*pow(RP,2) - 
9*d*gsl_sf_expint_Ei(k2*pow(RP,-1))*pow(d2,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*pow(RP,2) - 30*c2*d*d2*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) + 
20*c2*d2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) + 36*b2*d2*k1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) + 
18*k1*pow(c2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) - 6*b2*d2*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) - 
16*b2*c2*k1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) - 16*a2*d2*k1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) - 
3*pow(c2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) - 12*d*pow(d2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) + 
3*k1*pow(d2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) + 30*c2*d*d2*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) - 
20*c2*d2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) - 36*b2*d2*k2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) - 
18*k2*pow(c2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) + 6*b2*d2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) + 
16*b2*c2*k2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) + 16*a2*d2*k2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) + 
3*pow(c2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) + 12*d*pow(d2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) - 
3*k2*pow(d2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) + 6*c2*d2*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
36*b2*d*d2*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 18*d*pow(c2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
56*b2*c2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 56*a2*d2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
30*a2*c2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 15*pow(b2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
6*c2*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 36*b2*d*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
18*d*pow(c2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 56*b2*c2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
56*a2*d2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 30*a2*c2*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
15*pow(b2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 42*b2*c2*d*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
6*b2*d2*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 42*a2*d*d2*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 3*pow(c2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
80*a2*c2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 40*pow(b2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
48*a2*b2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 42*b2*c2*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
6*b2*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 42*a2*d*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
3*pow(c2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 80*a2*c2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
40*pow(b2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 48*a2*b2*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
6*b2*c2*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 48*a2*c2*d*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 6*a2*d2*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 
24*d*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 108*a2*b2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 
35*pow(a2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 6*b2*c2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
48*a2*c2*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 6*a2*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
24*d*pow(b2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 108*a2*b2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
35*pow(a2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 6*a2*c2*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 
54*a2*b2*d*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 3*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 
70*pow(a2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 6*a2*c2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 
54*a2*b2*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 3*pow(b2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 
70*pow(a2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 6*a2*b2*pow(M_E,k2*pow(RP,-1))*pow(k1,6)*pow(RP,2) + 
30*d*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,6)*pow(RP,2) + 6*a2*b2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6)*pow(RP,2) - 
30*d*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6)*pow(RP,2) + 3*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,7)*pow(RP,2) + 
3*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,7)*pow(RP,2) - 6*c2*d2*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
36*b2*d*d2*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 18*d*pow(c2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
56*b2*c2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 56*a2*d2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
30*a2*c2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 15*pow(b2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
6*c2*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 36*b2*d*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
18*d*pow(c2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 56*b2*c2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
56*a2*d2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 30*a2*c2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
15*pow(b2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 42*b2*c2*d*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
6*b2*d2*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 42*a2*d*d2*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 3*pow(c2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
80*a2*c2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 40*pow(b2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
48*a2*b2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 42*b2*c2*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
6*b2*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 42*a2*d*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
3*pow(c2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 80*a2*c2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
40*pow(b2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 48*a2*b2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
6*b2*c2*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 48*a2*c2*d*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 6*a2*d2*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 
24*d*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 108*a2*b2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 
35*pow(a2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 6*b2*c2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
48*a2*c2*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 6*a2*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
24*d*pow(b2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 108*a2*b2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
35*pow(a2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 6*a2*c2*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 
54*a2*b2*d*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 3*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 
70*pow(a2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 6*a2*c2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 
54*a2*b2*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 3*pow(b2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 
70*pow(a2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 6*a2*b2*pow(M_E,k1*pow(RP,-1))*pow(k2,6)*pow(RP,2) - 
30*d*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,6)*pow(RP,2) - 6*a2*b2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6)*pow(RP,2) + 
30*d*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6)*pow(RP,2) - 3*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,7)*pow(RP,2) - 
3*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,7)*pow(RP,2) - 48*c2*d*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 30*c2*d2*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 
90*b2*d*d2*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 45*d*k2*pow(c2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 36*b2*d2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 
112*b2*c2*k2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 112*a2*d2*k2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 
18*pow(c2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 16*b2*c2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 16*a2*d2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 
60*a2*c2*k2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 30*k2*pow(b2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 12*pow(d2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 
48*c2*d*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 30*c2*d2*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 90*b2*d*d2*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 
45*d*k1*pow(c2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 36*b2*d2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 112*b2*c2*k1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 
112*a2*d2*k1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 18*pow(c2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 
16*b2*c2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 16*a2*d2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 60*a2*c2*k1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 
30*k1*pow(b2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 12*pow(d2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 
9*gsl_sf_expint_Ei(-(k2*pow(RP,-1)))*pow(d2,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*pow(RP,3) - 
9*gsl_sf_expint_Ei(k2*pow(RP,-1))*pow(d2,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*pow(RP,3) + 48*c2*d*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 
30*c2*d2*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) + 90*b2*d*d2*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) + 
45*d*k1*pow(c2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 36*b2*d2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 
112*b2*c2*k1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 112*a2*d2*k1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 
18*pow(c2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) + 16*b2*c2*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) + 
16*a2*d2*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) + 60*a2*c2*k1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) + 
30*k1*pow(b2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 12*pow(d2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 
48*c2*d*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 30*c2*d2*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) - 
90*b2*d*d2*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) - 45*d*k2*pow(c2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 
36*b2*d2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 112*b2*c2*k2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 
112*a2*d2*k2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 18*pow(c2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) - 
16*b2*c2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) - 16*a2*d2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) - 
60*a2*c2*k2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) - 30*k2*pow(b2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 
12*pow(d2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 144*b2*c2*d*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
36*b2*d2*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 144*a2*d*d2*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
18*pow(c2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 240*a2*c2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
120*pow(b2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 144*a2*b2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
144*b2*c2*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 36*b2*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
144*a2*d*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 18*pow(c2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
240*a2*c2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 120*pow(b2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
144*a2*b2*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 42*b2*c2*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
210*a2*c2*d*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 42*a2*d2*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
105*d*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 432*a2*b2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
140*pow(a2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 42*b2*c2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
210*a2*c2*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 42*a2*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
105*d*pow(b2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 432*a2*b2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
140*pow(a2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 48*a2*c2*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 
288*a2*b2*d*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 24*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 
350*pow(a2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 48*a2*c2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 
288*a2*b2*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 24*pow(b2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 
350*pow(a2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 54*a2*b2*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,3) + 
189*d*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,3) - 54*a2*b2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,3) + 
189*d*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,3) + 30*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,6)*pow(RP,3) - 
30*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6)*pow(RP,3) - 144*b2*c2*d*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
36*b2*d2*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 144*a2*d*d2*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
18*pow(c2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 240*a2*c2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
120*pow(b2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 144*a2*b2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
144*b2*c2*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 36*b2*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
144*a2*d*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 18*pow(c2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
240*a2*c2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 120*pow(b2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
144*a2*b2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 42*b2*c2*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 
210*a2*c2*d*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 42*a2*d2*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 
105*d*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 432*a2*b2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 
140*pow(a2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 42*b2*c2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 
210*a2*c2*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 42*a2*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 
105*d*pow(b2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 432*a2*b2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 
140*pow(a2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 48*a2*c2*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 
288*a2*b2*d*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 24*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 
350*pow(a2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 48*a2*c2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 
288*a2*b2*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 24*pow(b2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 
350*pow(a2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 54*a2*b2*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 
189*d*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,3) + 54*a2*b2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 
189*d*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 30*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,6)*pow(RP,3) + 
30*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6)*pow(RP,3) - 
gsl_sf_expint_Ei(-(k1*pow(RP,-1)))*pow(d2,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*(4*RP*pow(d,2) + pow(d,3) + 9*d*pow(RP,2) + 9*pow(RP,3)) + 
gsl_sf_expint_Ei(k1*pow(RP,-1))*pow(d2,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*(4*RP*pow(d,2) + pow(d,3) + 9*d*pow(RP,2) + 9*pow(RP,3)) - 
48*c2*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 90*b2*d*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 288*b2*c2*d*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 
90*b2*d2*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 288*a2*d*d2*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 45*d*pow(c2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 
45*k2*pow(c2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 112*b2*c2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 112*a2*d2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 
480*a2*c2*k2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 240*k2*pow(b2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 
60*a2*c2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 288*a2*b2*k2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 
30*pow(b2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 48*c2*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 90*b2*d*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
288*b2*c2*d*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 90*b2*d2*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 288*a2*d*d2*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
45*d*pow(c2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 45*k1*pow(c2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 112*b2*c2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
112*a2*d2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 480*a2*c2*k1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
240*k1*pow(b2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 60*a2*c2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
288*a2*b2*k1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 30*pow(b2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
48*c2*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 90*b2*d*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 
288*b2*c2*d*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 90*b2*d2*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 
288*a2*d*d2*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 45*d*pow(c2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 
45*k1*pow(c2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 112*b2*c2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 
112*a2*d2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 480*a2*c2*k1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 
240*k1*pow(b2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 60*a2*c2*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 
288*a2*b2*k1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 30*pow(b2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 
48*c2*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 90*b2*d*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 
288*b2*c2*d*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 90*b2*d2*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 
288*a2*d*d2*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 45*d*pow(c2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 
45*k2*pow(c2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 112*b2*c2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 
112*a2*d2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 480*a2*c2*k2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 
240*k2*pow(b2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 60*a2*c2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 
288*a2*b2*k2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 30*pow(b2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 
144*b2*c2*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 630*a2*c2*d*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
144*a2*d2*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 315*d*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
1296*a2*b2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 420*pow(a2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
144*b2*c2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 630*a2*c2*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
144*a2*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 315*d*pow(b2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
1296*a2*b2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 420*pow(a2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
210*a2*c2*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 1152*a2*b2*d*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 
105*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 1400*pow(a2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 
210*a2*c2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 1152*a2*b2*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 
105*pow(b2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 1400*pow(a2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 
288*a2*b2*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,4) + 945*d*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,4) + 
288*a2*b2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,4) - 945*d*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,4) + 
189*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,4) + 189*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,4) - 
144*b2*c2*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 630*a2*c2*d*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
144*a2*d2*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 315*d*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
1296*a2*b2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 420*pow(a2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
144*b2*c2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 630*a2*c2*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
144*a2*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 315*d*pow(b2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
1296*a2*b2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 420*pow(a2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
210*a2*c2*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 1152*a2*b2*d*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 
105*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 1400*pow(a2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 
210*a2*c2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 1152*a2*b2*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 
105*pow(b2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 1400*pow(a2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 
288*a2*b2*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,4) - 945*d*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,4) - 
288*a2*b2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,4) + 945*d*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,4) - 
189*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,4) - 189*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,4) - 
288*b2*c2*d*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 90*b2*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 288*a2*d*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 
288*b2*c2*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 1260*a2*c2*d*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 288*a2*d2*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 
630*d*k2*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 45*pow(c2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 480*a2*c2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 
2592*a2*b2*k2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 240*pow(b2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 
288*a2*b2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 840*k2*pow(a2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 288*b2*c2*d*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
90*b2*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 288*a2*d*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 288*b2*c2*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
1260*a2*c2*d*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 288*a2*d2*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 630*d*k1*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
45*pow(c2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 480*a2*c2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 2592*a2*b2*k1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
240*pow(b2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 288*a2*b2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
840*k1*pow(a2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 288*b2*c2*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 
90*b2*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) + 288*a2*d*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 288*b2*c2*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) + 
1260*a2*c2*d*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 288*a2*d2*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) + 
630*d*k1*pow(b2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 45*pow(c2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 
480*a2*c2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 2592*a2*b2*k1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 
240*pow(b2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) + 288*a2*b2*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) + 
840*k1*pow(a2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 288*b2*c2*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 
90*b2*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) - 288*a2*d*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 288*b2*c2*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) - 
1260*a2*c2*d*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 288*a2*d2*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) - 
630*d*k2*pow(b2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 45*pow(c2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 
480*a2*c2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 2592*a2*b2*k2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 
240*pow(b2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) - 288*a2*b2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) - 
840*k2*pow(a2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 630*a2*c2*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 
3456*a2*b2*d*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 315*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 
4200*pow(a2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 630*a2*c2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 
3456*a2*b2*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 315*pow(b2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 
4200*pow(a2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 1152*a2*b2*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,5) + 
3780*d*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,5) - 1152*a2*b2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,5) + 
3780*d*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,5) + 945*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,5) - 
945*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,5) - 630*a2*c2*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 
3456*a2*b2*d*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 315*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 
4200*pow(a2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 630*a2*c2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 
3456*a2*b2*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 315*pow(b2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 
4200*pow(a2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 1152*a2*b2*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 
3780*d*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,5) + 1152*a2*b2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 
3780*d*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 945*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,5) + 
945*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,5) - 288*b2*c2*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 1260*a2*c2*d*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 
288*a2*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 1260*a2*c2*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 6912*a2*b2*d*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 
630*d*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 630*k2*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 2592*a2*b2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 
8400*k2*pow(a2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 840*pow(a2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 288*b2*c2*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 
1260*a2*c2*d*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 288*a2*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 1260*a2*c2*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 
6912*a2*b2*d*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 630*d*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 630*k1*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 
2592*a2*b2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 8400*k1*pow(a2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 
840*pow(a2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 288*b2*c2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 
1260*a2*c2*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) + 288*a2*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) + 
1260*a2*c2*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 6912*a2*b2*d*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 
630*d*pow(b2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) + 630*k1*pow(b2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) + 
2592*a2*b2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) + 8400*k1*pow(a2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 
840*pow(a2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 288*b2*c2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 
1260*a2*c2*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) - 288*a2*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) - 
1260*a2*c2*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 6912*a2*b2*d*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 
630*d*pow(b2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) - 630*k2*pow(b2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) - 
2592*a2*b2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) - 8400*k2*pow(a2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 
840*pow(a2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 3456*a2*b2*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,6) + 
11340*d*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,6) + 3456*a2*b2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,6) - 
11340*d*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,6) + 3780*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,6) + 
3780*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,6) - 3456*a2*b2*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,6) - 
11340*d*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,6) - 3456*a2*b2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,6) + 
11340*d*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,6) - 3780*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,6) - 
3780*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,6) - 1260*a2*c2*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 6912*a2*b2*d*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 
6912*a2*b2*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 22680*d*k2*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 630*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 
8400*pow(a2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) + 1260*a2*c2*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 6912*a2*b2*d*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 
6912*a2*b2*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 22680*d*k1*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 630*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 
8400*pow(a2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) - 1260*a2*c2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) + 
6912*a2*b2*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) - 6912*a2*b2*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) + 
22680*d*k1*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) - 630*pow(b2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) - 
8400*pow(a2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) + 1260*a2*c2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) - 
6912*a2*b2*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) + 6912*a2*b2*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) - 
22680*d*k2*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) + 630*pow(b2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) + 
8400*pow(a2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) + 11340*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,7) - 
11340*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,7) - 11340*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,7) + 
11340*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,7) - 6912*a2*b2*pow(M_E,k1*pow(RP,-1))*pow(RP,8) - 
22680*d*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,8) - 22680*k2*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,8) + 6912*a2*b2*pow(M_E,k2*pow(RP,-1))*pow(RP,8) + 
22680*d*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,8) + 22680*k1*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,8) + 6912*a2*b2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,8) - 
22680*d*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,8) + 22680*k1*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,8) - 
6912*a2*b2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,8) + 22680*d*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,8) - 
22680*k2*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,8) - 22680*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,9) + 22680*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,9) - 
22680*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,9) + 22680*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,9)))/2.;

return res;
}

double SDH_Integral_z_2233_case_2_sub_1( double AP, double RP, double k1, double k2, double a2, double b2, double c2, double d2, double d )
{
double res =(-3*AP*pow(d,-4)*pow(M_E,-((2*d + k2)*pow(RP,-1)))*pow(RP,2)*(12*c2*d2*k2*RP*pow(d,2)*pow(M_E,d*pow(RP,-1)) + 2*c2*d2*k2*pow(d,3)*pow(M_E,d*pow(RP,-1)) + 
4*c2*d2*RP*pow(d,3)*pow(M_E,d*pow(RP,-1)) + 6*b2*d2*k2*RP*pow(d,3)*pow(M_E,d*pow(RP,-1)) + 3*k2*RP*pow(c2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1)) + 
3*d*k2*RP*pow(d2,2)*pow(M_E,d*pow(RP,-1)) + k2*pow(d,2)*pow(d2,2)*pow(M_E,d*pow(RP,-1)) + 5*RP*pow(d,2)*pow(d2,2)*pow(M_E,d*pow(RP,-1)) + 
pow(d,3)*pow(d2,2)*pow(M_E,d*pow(RP,-1)) - 12*c2*d2*k2*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1)) + 2*c2*d2*k2*pow(d,3)*pow(M_E,3*d*pow(RP,-1)) + 
4*c2*d2*RP*pow(d,3)*pow(M_E,3*d*pow(RP,-1)) + 6*b2*d2*k2*RP*pow(d,3)*pow(M_E,3*d*pow(RP,-1)) + 3*k2*RP*pow(c2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1)) + 
3*d*k2*RP*pow(d2,2)*pow(M_E,3*d*pow(RP,-1)) - k2*pow(d,2)*pow(d2,2)*pow(M_E,3*d*pow(RP,-1)) - 5*RP*pow(d,2)*pow(d2,2)*pow(M_E,3*d*pow(RP,-1)) + 
pow(d,3)*pow(d2,2)*pow(M_E,3*d*pow(RP,-1)) - 22*c2*d2*RP*pow(d,3)*pow(M_E,k2*pow(RP,-1)) - 4*c2*d2*pow(d,4)*pow(M_E,k2*pow(RP,-1)) - 
26*b2*d2*RP*pow(d,4)*pow(M_E,k2*pow(RP,-1)) - 13*RP*pow(c2,2)*pow(d,4)*pow(M_E,k2*pow(RP,-1)) - 4*b2*d2*pow(d,5)*pow(M_E,k2*pow(RP,-1)) - 
30*b2*c2*RP*pow(d,5)*pow(M_E,k2*pow(RP,-1)) - 30*a2*d2*RP*pow(d,5)*pow(M_E,k2*pow(RP,-1)) - 2*pow(c2,2)*pow(d,5)*pow(M_E,k2*pow(RP,-1)) - 
4*b2*c2*pow(d,6)*pow(M_E,k2*pow(RP,-1)) - 4*a2*d2*pow(d,6)*pow(M_E,k2*pow(RP,-1)) - 34*a2*c2*RP*pow(d,6)*pow(M_E,k2*pow(RP,-1)) - 
17*RP*pow(b2,2)*pow(d,6)*pow(M_E,k2*pow(RP,-1)) - 4*a2*c2*pow(d,7)*pow(M_E,k2*pow(RP,-1)) - 38*a2*b2*RP*pow(d,7)*pow(M_E,k2*pow(RP,-1)) - 
2*pow(b2,2)*pow(d,7)*pow(M_E,k2*pow(RP,-1)) - 4*a2*b2*pow(d,8)*pow(M_E,k2*pow(RP,-1)) - 21*RP*pow(a2,2)*pow(d,8)*pow(M_E,k2*pow(RP,-1)) - 
2*pow(a2,2)*pow(d,9)*pow(M_E,k2*pow(RP,-1)) - 8*RP*pow(d,2)*pow(d2,2)*pow(M_E,k2*pow(RP,-1)) - 2*pow(d,3)*pow(d2,2)*pow(M_E,k2*pow(RP,-1)) + 
2*c2*d2*RP*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 2*b2*d2*RP*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1)) + RP*pow(c2,2)*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 
2*b2*c2*RP*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 2*a2*d2*RP*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 2*a2*c2*RP*pow(d,6)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 
RP*pow(b2,2)*pow(d,6)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 2*a2*b2*RP*pow(d,7)*pow(M_E,(2*d + k2)*pow(RP,-1)) + RP*pow(a2,2)*pow(d,8)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 
2*RP*pow(d,2)*pow(d2,2)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 6*c2*d*d2*RP*pow(M_E,d*pow(RP,-1))*pow(k2,2) + 2*c2*d2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2) + 
14*b2*d2*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2) + 7*RP*pow(c2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2) + 2*b2*d2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2) + 
8*b2*c2*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2) + 8*a2*d2*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2) + pow(c2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2) + 
6*c2*d*d2*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) - 2*c2*d2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) - 14*b2*d2*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) - 
7*RP*pow(c2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) + 2*b2*d2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) + 
8*b2*c2*RP*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) + 8*a2*d2*RP*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) + 
pow(c2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) + 6*b2*d*d2*RP*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 3*d*RP*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 
2*b2*d2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 16*b2*c2*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 16*a2*d2*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 
pow(c2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 2*b2*c2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 2*a2*d2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 
10*a2*c2*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 5*RP*pow(b2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 6*b2*d*d2*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) + 
3*d*RP*pow(c2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) - 2*b2*d2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) - 16*b2*c2*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) - 
16*a2*d2*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) - pow(c2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) + 2*b2*c2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) + 
2*a2*d2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) + 10*a2*c2*RP*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) + 
5*RP*pow(b2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) + 6*b2*c2*d*RP*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 6*a2*d*d2*RP*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 
2*b2*c2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 2*a2*d2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 18*a2*c2*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 
9*RP*pow(b2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 2*a2*c2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 12*a2*b2*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 
pow(b2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 6*b2*c2*d*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) + 6*a2*d*d2*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) - 
2*b2*c2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) - 2*a2*d2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) - 18*a2*c2*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) - 
9*RP*pow(b2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) + 2*a2*c2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) + 
12*a2*b2*RP*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) + pow(b2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) + 6*a2*c2*d*RP*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 
3*d*RP*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 2*a2*c2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 20*a2*b2*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 
pow(b2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 2*a2*b2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 7*RP*pow(a2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 
6*a2*c2*d*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) + 3*d*RP*pow(b2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) - 2*a2*c2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) - 
20*a2*b2*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) - pow(b2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) + 2*a2*b2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) + 
7*RP*pow(a2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) + 6*a2*b2*d*RP*pow(M_E,d*pow(RP,-1))*pow(k2,6) + 2*a2*b2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,6) + 
11*RP*pow(a2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,6) + pow(a2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,6) + 6*a2*b2*d*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,6) - 
2*a2*b2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6) - 11*RP*pow(a2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6) + 
pow(a2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6) + 3*d*RP*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,7) + pow(a2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,7) + 
3*d*RP*pow(a2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,7) - pow(a2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,7) + 30*c2*d*d2*k2*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 
20*c2*d2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 36*b2*d2*k2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 18*k2*pow(c2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 
6*b2*d2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 16*b2*c2*k2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 16*a2*d2*k2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 
3*pow(c2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 12*d*pow(d2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 3*k2*pow(d2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 
30*c2*d*d2*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) - 20*c2*d2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) - 36*b2*d2*k2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) - 
18*k2*pow(c2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) + 6*b2*d2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) + 
16*b2*c2*k2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) + 16*a2*d2*k2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) + 
3*pow(c2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) + 12*d*pow(d2,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) - 3*k2*pow(d2,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) - 
56*c2*d2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 84*b2*d2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 42*pow(c2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 
120*b2*c2*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 120*a2*d2*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 164*a2*c2*pow(d,5)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 
82*pow(b2,2)*pow(d,5)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 216*a2*b2*pow(d,6)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 
138*pow(a2,2)*pow(d,7)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 15*d*pow(d2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 
4*c2*d2*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 4*b2*c2*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 
4*a2*d2*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 8*a2*c2*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 
4*pow(b2,2)*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 12*a2*b2*pow(d,6)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 
8*pow(a2,2)*pow(d,7)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) - 9*d*pow(d2,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 
6*c2*d2*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 36*b2*d*d2*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 18*d*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
56*b2*c2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 56*a2*d2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
30*a2*c2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 15*pow(b2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
6*c2*d2*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 36*b2*d*d2*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
18*d*pow(c2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 56*b2*c2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
56*a2*d2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 30*a2*c2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
15*pow(b2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 42*b2*c2*d*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
6*b2*d2*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 42*a2*d*d2*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 3*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
80*a2*c2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 40*pow(b2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
48*a2*b2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 42*b2*c2*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
6*b2*d2*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 42*a2*d*d2*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
3*pow(c2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 80*a2*c2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
40*pow(b2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 48*a2*b2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
6*b2*c2*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 48*a2*c2*d*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 6*a2*d2*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
24*d*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 108*a2*b2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
35*pow(a2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 6*b2*c2*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
48*a2*c2*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 6*a2*d2*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
24*d*pow(b2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 108*a2*b2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
35*pow(a2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 6*a2*c2*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 
54*a2*b2*d*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 3*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 
70*pow(a2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 6*a2*c2*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 
54*a2*b2*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 3*pow(b2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 
70*pow(a2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 6*a2*b2*pow(M_E,d*pow(RP,-1))*pow(k2,6)*pow(RP,2) + 
30*d*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,6)*pow(RP,2) - 6*a2*b2*pow(M_E,3*d*pow(RP,-1))*pow(k2,6)*pow(RP,2) + 
30*d*pow(a2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6)*pow(RP,2) + 3*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,7)*pow(RP,2) - 
3*pow(a2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,7)*pow(RP,2) + 48*c2*d*d2*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 30*c2*d2*k2*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 
90*b2*d*d2*k2*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 45*d*k2*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 36*b2*d2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 
112*b2*c2*k2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 112*a2*d2*k2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 
18*pow(c2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 16*b2*c2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 16*a2*d2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 
60*a2*c2*k2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 30*k2*pow(b2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 12*pow(d2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 
48*c2*d*d2*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 30*c2*d2*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) + 90*b2*d*d2*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) + 
45*d*k2*pow(c2,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 36*b2*d2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 112*b2*c2*k2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 
112*a2*d2*k2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 18*pow(c2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) + 
16*b2*c2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) + 16*a2*d2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) + 60*a2*c2*k2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) + 
30*k2*pow(b2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 12*pow(d2,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 78*c2*d*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 
162*b2*d2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 81*pow(c2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 314*b2*c2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 
314*a2*d2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 558*a2*c2*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 279*pow(b2,2)*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 
918*a2*b2*pow(d,5)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 709*pow(a2,2)*pow(d,6)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 12*pow(d2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 
18*c2*d*d2*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) - 18*b2*d2*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) - 
9*pow(c2,2)*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) - 6*b2*c2*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) - 
6*a2*d2*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) + 18*a2*c2*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) + 
9*pow(b2,2)*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) + 54*a2*b2*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) + 
51*pow(a2,2)*pow(d,6)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) + 12*pow(d2,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) + 
144*b2*c2*d*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 36*b2*d2*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 144*a2*d*d2*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
18*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 240*a2*c2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
120*pow(b2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 144*a2*b2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
144*b2*c2*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 36*b2*d2*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
144*a2*d*d2*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 18*pow(c2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
240*a2*c2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 120*pow(b2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
144*a2*b2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 42*b2*c2*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
210*a2*c2*d*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 42*a2*d2*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
105*d*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 432*a2*b2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
140*pow(a2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 42*b2*c2*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
210*a2*c2*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 42*a2*d2*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
105*d*pow(b2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 432*a2*b2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
140*pow(a2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 48*a2*c2*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 
288*a2*b2*d*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 24*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 
350*pow(a2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 48*a2*c2*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 
288*a2*b2*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 24*pow(b2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 
350*pow(a2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 54*a2*b2*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,3) + 
189*d*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 54*a2*b2*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,3) + 
189*d*pow(a2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,3) + 30*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,6)*pow(RP,3) - 
30*pow(a2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6)*pow(RP,3) + gsl_sf_expint_Ei(-(d*pow(RP,-1)))*pow(d2,2)*pow(M_E,(d + k2)*pow(RP,-1))*
(-4*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1))) + pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1))) + 9*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,2) - 
9*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,3)) - gsl_sf_expint_Ei(-(k2*pow(RP,-1)))*pow(d2,2)*pow(M_E,(d + k2)*pow(RP,-1))*
(-4*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1))) + pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1))) + 9*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,2) - 
9*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,3)) + 48*c2*d2*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 90*b2*d*d2*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 
288*b2*c2*d*k2*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 90*b2*d2*k2*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 288*a2*d*d2*k2*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 
45*d*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 45*k2*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 112*b2*c2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 
112*a2*d2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 480*a2*c2*k2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 
240*k2*pow(b2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 60*a2*c2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 
288*a2*b2*k2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 30*pow(b2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 48*c2*d2*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 
90*b2*d*d2*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 288*b2*c2*d*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 90*b2*d2*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 
288*a2*d*d2*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 45*d*pow(c2,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 45*k2*pow(c2,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 
112*b2*c2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 112*a2*d2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 
480*a2*c2*k2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 240*k2*pow(b2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 
60*a2*c2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 288*a2*b2*k2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 
30*pow(b2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 48*c2*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 180*b2*d*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 
90*d*pow(c2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 544*b2*c2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 544*a2*d2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 
1380*a2*c2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 690*pow(b2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 
3024*a2*b2*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 2954*pow(a2,2)*pow(d,5)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 48*c2*d2*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) - 
32*b2*c2*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) - 32*a2*d2*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) + 
144*a2*b2*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) + 224*pow(a2,2)*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) + 
144*b2*c2*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 630*a2*c2*d*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 144*a2*d2*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 
315*d*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 1296*a2*b2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 
420*pow(a2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 144*b2*c2*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 
630*a2*c2*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 144*a2*d2*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 
315*d*pow(b2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 1296*a2*b2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 
420*pow(a2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 210*a2*c2*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 
1152*a2*b2*d*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 105*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 
1400*pow(a2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 210*a2*c2*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 
1152*a2*b2*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 105*pow(b2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 
1400*pow(a2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 288*a2*b2*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,4) + 
945*d*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,4) - 288*a2*b2*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,4) + 
945*d*pow(a2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,4) + 189*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,4) - 
189*pow(a2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,4) + 288*b2*c2*d*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 90*b2*d2*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 
288*a2*d*d2*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 288*b2*c2*k2*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 1260*a2*c2*d*k2*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 
288*a2*d2*k2*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 630*d*k2*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 45*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 
480*a2*c2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 2592*a2*b2*k2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 
240*pow(b2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 288*a2*b2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 
840*k2*pow(a2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 288*b2*c2*d*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 90*b2*d2*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 
288*a2*d*d2*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 288*b2*c2*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 1260*a2*c2*d*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 
288*a2*d2*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 630*d*k2*pow(b2,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 45*pow(c2,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 
480*a2*c2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 2592*a2*b2*k2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 
240*pow(b2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 288*a2*b2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 
840*k2*pow(a2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 576*b2*c2*d*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 90*b2*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 
576*a2*d*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 45*pow(c2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 2370*a2*c2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 
1185*pow(b2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 7488*a2*b2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 
9765*pow(a2,2)*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 90*b2*d2*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) + 
45*pow(c2,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) - 150*a2*c2*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) - 
75*pow(b2,2)*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) + 525*pow(a2,2)*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) + 
630*a2*c2*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 3456*a2*b2*d*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 
315*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 4200*pow(a2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 
630*a2*c2*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 3456*a2*b2*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 
315*pow(b2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 4200*pow(a2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 
1152*a2*b2*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,5) + 3780*d*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 
1152*a2*b2*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,5) + 3780*d*pow(a2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,5) + 
945*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,5) - 945*pow(a2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,5) + 288*b2*c2*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 
1260*a2*c2*d*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 288*a2*d2*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 1260*a2*c2*k2*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 
6912*a2*b2*d*k2*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 630*d*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 630*k2*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 
2592*a2*b2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 8400*k2*pow(a2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 
840*pow(a2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 288*b2*c2*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) + 1260*a2*c2*d*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 
288*a2*d2*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 1260*a2*c2*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) + 6912*a2*b2*d*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) + 
630*d*pow(b2,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 630*k2*pow(b2,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 2592*a2*b2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 
8400*k2*pow(a2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) + 840*pow(a2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 
288*b2*c2*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 2520*a2*c2*d*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 288*a2*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 
1260*d*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 12960*a2*b2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 
24360*pow(a2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 288*b2*c2*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,6) + 288*a2*d2*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,6) - 
864*a2*b2*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,6) + 3456*a2*b2*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,6) + 
11340*d*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,6) - 3456*a2*b2*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,6) + 
11340*d*pow(a2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,6) + 3780*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,6) - 
3780*pow(a2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,6) + 1260*a2*c2*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 6912*a2*b2*d*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 
6912*a2*b2*k2*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 22680*d*k2*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 630*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 
8400*pow(a2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 1260*a2*c2*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) + 6912*a2*b2*d*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) - 
6912*a2*b2*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) + 22680*d*k2*pow(a2,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) - 630*pow(b2,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) - 
8400*pow(a2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) - 1260*a2*c2*pow(M_E,k2*pow(RP,-1))*pow(RP,7) - 13824*a2*b2*d*pow(M_E,k2*pow(RP,-1))*pow(RP,7) - 
630*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) - 42420*pow(a2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 1260*a2*c2*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,7) + 
630*pow(b2,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,7) - 2940*pow(a2,2)*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,7) + 
11340*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,7) - 11340*pow(a2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,7) + 
6912*a2*b2*pow(M_E,d*pow(RP,-1))*pow(RP,8) + 22680*d*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,8) + 22680*k2*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,8) - 
6912*a2*b2*pow(M_E,3*d*pow(RP,-1))*pow(RP,8) + 22680*d*pow(a2,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,8) - 22680*k2*pow(a2,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,8) - 
6912*a2*b2*pow(M_E,k2*pow(RP,-1))*pow(RP,8) - 45360*d*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,8) + 6912*a2*b2*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,8) + 
22680*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,9) - 22680*pow(a2,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,9) - 22680*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,9) + 
22680*pow(a2,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,9)))/2.;

return res;
}

double SDH_Integral_z_2233_case_2_sub_2( double AP, double RP, double k1, double k2, double a2, double b2, double c2, double d2, double d )
{
double res =(-3*AP*pow(d,-4)*pow(M_E,-((2*d + k1)*pow(RP,-1)))*pow(RP,2)*(-12*c2*d2*k1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1)) - 2*c2*d2*k1*pow(d,3)*pow(M_E,d*pow(RP,-1)) - 
4*c2*d2*RP*pow(d,3)*pow(M_E,d*pow(RP,-1)) - 6*b2*d2*k1*RP*pow(d,3)*pow(M_E,d*pow(RP,-1)) - 3*k1*RP*pow(c2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1)) - 
3*d*k1*RP*pow(d2,2)*pow(M_E,d*pow(RP,-1)) - k1*pow(d,2)*pow(d2,2)*pow(M_E,d*pow(RP,-1)) - 5*RP*pow(d,2)*pow(d2,2)*pow(M_E,d*pow(RP,-1)) - 
pow(d,3)*pow(d2,2)*pow(M_E,d*pow(RP,-1)) + 22*c2*d2*RP*pow(d,3)*pow(M_E,k1*pow(RP,-1)) + 4*c2*d2*pow(d,4)*pow(M_E,k1*pow(RP,-1)) + 
26*b2*d2*RP*pow(d,4)*pow(M_E,k1*pow(RP,-1)) + 13*RP*pow(c2,2)*pow(d,4)*pow(M_E,k1*pow(RP,-1)) + 4*b2*d2*pow(d,5)*pow(M_E,k1*pow(RP,-1)) + 
30*b2*c2*RP*pow(d,5)*pow(M_E,k1*pow(RP,-1)) + 30*a2*d2*RP*pow(d,5)*pow(M_E,k1*pow(RP,-1)) + 2*pow(c2,2)*pow(d,5)*pow(M_E,k1*pow(RP,-1)) + 
4*b2*c2*pow(d,6)*pow(M_E,k1*pow(RP,-1)) + 4*a2*d2*pow(d,6)*pow(M_E,k1*pow(RP,-1)) + 34*a2*c2*RP*pow(d,6)*pow(M_E,k1*pow(RP,-1)) + 
17*RP*pow(b2,2)*pow(d,6)*pow(M_E,k1*pow(RP,-1)) + 4*a2*c2*pow(d,7)*pow(M_E,k1*pow(RP,-1)) + 38*a2*b2*RP*pow(d,7)*pow(M_E,k1*pow(RP,-1)) + 
2*pow(b2,2)*pow(d,7)*pow(M_E,k1*pow(RP,-1)) + 4*a2*b2*pow(d,8)*pow(M_E,k1*pow(RP,-1)) + 21*RP*pow(a2,2)*pow(d,8)*pow(M_E,k1*pow(RP,-1)) + 
2*pow(a2,2)*pow(d,9)*pow(M_E,k1*pow(RP,-1)) + 8*RP*pow(d,2)*pow(d2,2)*pow(M_E,k1*pow(RP,-1)) + 2*pow(d,3)*pow(d2,2)*pow(M_E,k1*pow(RP,-1)) + 
4*RP*gsl_sf_expint_Ei(-(k1*pow(RP,-1)))*pow(d,2)*pow(d2,2)*pow(M_E,(d + k1)*pow(RP,-1)) - 
4*RP*gsl_sf_expint_Ei(k1*pow(RP,-1))*pow(d,2)*pow(d2,2)*pow(M_E,(d + k1)*pow(RP,-1)) + 
gsl_sf_expint_Ei(-(k1*pow(RP,-1)))*pow(d,3)*pow(d2,2)*pow(M_E,(d + k1)*pow(RP,-1)) - gsl_sf_expint_Ei(k1*pow(RP,-1))*pow(d,3)*pow(d2,2)*pow(M_E,(d + k1)*pow(RP,-1)) - 
2*c2*d2*RP*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 2*b2*d2*RP*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1)) - RP*pow(c2,2)*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 
2*b2*c2*RP*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 2*a2*d2*RP*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 2*a2*c2*RP*pow(d,6)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 
RP*pow(b2,2)*pow(d,6)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 2*a2*b2*RP*pow(d,7)*pow(M_E,(2*d + k1)*pow(RP,-1)) - RP*pow(a2,2)*pow(d,8)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 
2*RP*pow(d,2)*pow(d2,2)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 12*c2*d2*k1*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1)) + 
2*c2*d2*k1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1)) - 4*c2*d2*RP*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1)) - 6*b2*d2*k1*RP*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1)) - 
3*k1*RP*pow(c2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1)) - 3*d*k1*RP*pow(d2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1)) - 
k1*pow(d,2)*pow(d2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1)) + 5*RP*pow(d,2)*pow(d2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1)) + pow(d,3)*pow(d2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1)) - 
6*c2*d*d2*RP*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 2*c2*d2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 14*b2*d2*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 
7*RP*pow(c2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 2*b2*d2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 8*b2*c2*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 
8*a2*d2*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - pow(c2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 6*c2*d*d2*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) - 
2*c2*d2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) + 14*b2*d2*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) + 
7*RP*pow(c2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) + 2*b2*d2*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) - 
8*b2*c2*RP*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) - 8*a2*d2*RP*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) + 
pow(c2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) - 6*b2*d*d2*RP*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 3*d*RP*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 
2*b2*d2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 16*b2*c2*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 16*a2*d2*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 
pow(c2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 2*b2*c2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 2*a2*d2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 
10*a2*c2*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 5*RP*pow(b2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 
6*b2*d*d2*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - 3*d*RP*pow(c2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - 
2*b2*d2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) + 16*b2*c2*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) + 
16*a2*d2*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - pow(c2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) + 
2*b2*c2*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) + 2*a2*d2*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - 
10*a2*c2*RP*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - 5*RP*pow(b2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - 
6*b2*c2*d*RP*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 6*a2*d*d2*RP*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 2*b2*c2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 
2*a2*d2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 18*a2*c2*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 9*RP*pow(b2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 
2*a2*c2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 12*a2*b2*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - pow(b2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 
6*b2*c2*d*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) - 6*a2*d*d2*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) - 
2*b2*c2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) - 2*a2*d2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) + 
18*a2*c2*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) + 9*RP*pow(b2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) + 
2*a2*c2*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) - 12*a2*b2*RP*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) + 
pow(b2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) - 6*a2*c2*d*RP*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 3*d*RP*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 
2*a2*c2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 20*a2*b2*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5) - pow(b2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 
2*a2*b2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 7*RP*pow(a2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 6*a2*c2*d*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) - 
3*d*RP*pow(b2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) - 2*a2*c2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) + 
20*a2*b2*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) - pow(b2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) + 
2*a2*b2*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) - 7*RP*pow(a2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) - 
6*a2*b2*d*RP*pow(M_E,d*pow(RP,-1))*pow(k1,6) - 2*a2*b2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,6) - 11*RP*pow(a2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,6) - 
pow(a2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,6) - 6*a2*b2*d*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6) - 
2*a2*b2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6) + 11*RP*pow(a2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6) + 
pow(a2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6) - 3*d*RP*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,7) - pow(a2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,7) - 
3*d*RP*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,7) - pow(a2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,7) - 
30*c2*d*d2*k1*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 20*c2*d2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 36*b2*d2*k1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 
18*k1*pow(c2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 6*b2*d2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 16*b2*c2*k1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 
16*a2*d2*k1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 3*pow(c2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 12*d*pow(d2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 
3*k1*pow(d2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 56*c2*d2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 84*b2*d2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
42*pow(c2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 120*b2*c2*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 120*a2*d2*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
164*a2*c2*pow(d,5)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 82*pow(b2,2)*pow(d,5)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 216*a2*b2*pow(d,6)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
138*pow(a2,2)*pow(d,7)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 15*d*pow(d2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
9*d*gsl_sf_expint_Ei(-(k1*pow(RP,-1)))*pow(d2,2)*pow(M_E,(d + k1)*pow(RP,-1))*pow(RP,2) - 
9*d*gsl_sf_expint_Ei(k1*pow(RP,-1))*pow(d2,2)*pow(M_E,(d + k1)*pow(RP,-1))*pow(RP,2) - 4*c2*d2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 
4*b2*c2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 4*a2*d2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 
8*a2*c2*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 4*pow(b2,2)*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 
12*a2*b2*pow(d,6)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 8*pow(a2,2)*pow(d,7)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) - 
9*d*pow(d2,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 30*c2*d*d2*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) - 
20*c2*d2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) - 36*b2*d2*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) - 
18*k1*pow(c2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) + 6*b2*d2*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) + 
16*b2*c2*k1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) + 16*a2*d2*k1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) + 
3*pow(c2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) + 12*d*pow(d2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) - 
3*k1*pow(d2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) - 6*c2*d2*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 36*b2*d*d2*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
18*d*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 56*b2*c2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
56*a2*d2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 30*a2*c2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
15*pow(b2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 6*c2*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
36*b2*d*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 18*d*pow(c2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
56*b2*c2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 56*a2*d2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
30*a2*c2*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 15*pow(b2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
42*b2*c2*d*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 6*b2*d2*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 42*a2*d*d2*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
3*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 80*a2*c2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
40*pow(b2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 48*a2*b2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
42*b2*c2*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 6*b2*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
42*a2*d*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 3*pow(c2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
80*a2*c2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 40*pow(b2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
48*a2*b2*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 6*b2*c2*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
48*a2*c2*d*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 6*a2*d2*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 24*d*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
108*a2*b2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 35*pow(a2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
6*b2*c2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 48*a2*c2*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
6*a2*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 24*d*pow(b2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
108*a2*b2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 35*pow(a2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
6*a2*c2*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 54*a2*b2*d*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 3*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 
70*pow(a2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 6*a2*c2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 
54*a2*b2*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 3*pow(b2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 
70*pow(a2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 6*a2*b2*pow(M_E,d*pow(RP,-1))*pow(k1,6)*pow(RP,2) - 
30*d*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,6)*pow(RP,2) - 6*a2*b2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6)*pow(RP,2) + 
30*d*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6)*pow(RP,2) - 3*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,7)*pow(RP,2) - 
3*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,7)*pow(RP,2) - 48*c2*d*d2*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 30*c2*d2*k1*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 
90*b2*d*d2*k1*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 45*d*k1*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 36*b2*d2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 
112*b2*c2*k1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 112*a2*d2*k1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 
18*pow(c2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 16*b2*c2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 16*a2*d2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 
60*a2*c2*k1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 30*k1*pow(b2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 12*pow(d2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 
78*c2*d*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 162*b2*d2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 81*pow(c2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 
314*b2*c2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 314*a2*d2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 558*a2*c2*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 
279*pow(b2,2)*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 918*a2*b2*pow(d,5)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 
709*pow(a2,2)*pow(d,6)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 12*pow(d2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 
9*gsl_sf_expint_Ei(-(k1*pow(RP,-1)))*pow(d2,2)*pow(M_E,(d + k1)*pow(RP,-1))*pow(RP,3) - 
9*gsl_sf_expint_Ei(k1*pow(RP,-1))*pow(d2,2)*pow(M_E,(d + k1)*pow(RP,-1))*pow(RP,3) + 18*c2*d*d2*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) + 
18*b2*d2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) + 9*pow(c2,2)*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) + 
6*b2*c2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) + 6*a2*d2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 
18*a2*c2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 9*pow(b2,2)*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 
54*a2*b2*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 51*pow(a2,2)*pow(d,6)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 
12*pow(d2,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 48*c2*d*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 30*c2*d2*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 
90*b2*d*d2*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 45*d*k1*pow(c2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 
36*b2*d2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 112*b2*c2*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 
112*a2*d2*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 18*pow(c2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 
16*b2*c2*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 16*a2*d2*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 
60*a2*c2*k1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 30*k1*pow(b2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 
12*pow(d2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 144*b2*c2*d*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 36*b2*d2*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
144*a2*d*d2*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 18*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
240*a2*c2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 120*pow(b2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
144*a2*b2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 144*b2*c2*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
36*b2*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 144*a2*d*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
18*pow(c2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 240*a2*c2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
120*pow(b2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 144*a2*b2*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
42*b2*c2*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 210*a2*c2*d*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 42*a2*d2*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
105*d*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 432*a2*b2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
140*pow(a2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 42*b2*c2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
210*a2*c2*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 42*a2*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
105*d*pow(b2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 432*a2*b2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
140*pow(a2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 48*a2*c2*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 
288*a2*b2*d*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 24*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 
350*pow(a2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 48*a2*c2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 
288*a2*b2*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 24*pow(b2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 
350*pow(a2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 54*a2*b2*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,3) - 
189*d*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,3) + 54*a2*b2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,3) - 
189*d*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,3) - 30*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,6)*pow(RP,3) + 
30*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6)*pow(RP,3) - 
gsl_sf_expint_Ei(-(d*pow(RP,-1)))*pow(d2,2)*pow(M_E,(d + k1)*pow(RP,-1))*(4*RP*pow(d,2) + pow(d,3) + 9*d*pow(RP,2) + 9*pow(RP,3)) + 
gsl_sf_expint_Ei(d*pow(RP,-1))*pow(d2,2)*pow(M_E,(d + k1)*pow(RP,-1))*(4*RP*pow(d,2) + pow(d,3) + 9*d*pow(RP,2) + 9*pow(RP,3)) - 
48*c2*d2*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 90*b2*d*d2*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 288*b2*c2*d*k1*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
90*b2*d2*k1*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 288*a2*d*d2*k1*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 45*d*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
45*k1*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 112*b2*c2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 112*a2*d2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
480*a2*c2*k1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 240*k1*pow(b2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
60*a2*c2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 288*a2*b2*k1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 30*pow(b2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 
48*c2*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 180*b2*d*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 90*d*pow(c2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 
544*b2*c2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 544*a2*d2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 1380*a2*c2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 
690*pow(b2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 3024*a2*b2*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 
2954*pow(a2,2)*pow(d,5)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 48*c2*d2*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) - 
32*b2*c2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) - 32*a2*d2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) + 
144*a2*b2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) + 224*pow(a2,2)*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) - 
48*c2*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 90*b2*d*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 288*b2*c2*d*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 
90*b2*d2*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 288*a2*d*d2*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 
45*d*pow(c2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 45*k1*pow(c2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 
112*b2*c2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 112*a2*d2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 
480*a2*c2*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 240*k1*pow(b2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 
60*a2*c2*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 288*a2*b2*k1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 
30*pow(b2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 144*b2*c2*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
630*a2*c2*d*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 144*a2*d2*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
315*d*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 1296*a2*b2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
420*pow(a2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 144*b2*c2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
630*a2*c2*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 144*a2*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
315*d*pow(b2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 1296*a2*b2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
420*pow(a2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 210*a2*c2*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 
1152*a2*b2*d*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 105*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 
1400*pow(a2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 210*a2*c2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 
1152*a2*b2*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 105*pow(b2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 
1400*pow(a2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 288*a2*b2*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,4) - 
945*d*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,4) - 288*a2*b2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,4) + 
945*d*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,4) - 189*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,4) - 
189*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,4) - 288*b2*c2*d*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 90*b2*d2*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 
288*a2*d*d2*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 288*b2*c2*k1*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 1260*a2*c2*d*k1*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 
288*a2*d2*k1*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 630*d*k1*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 45*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 
480*a2*c2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 2592*a2*b2*k1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 
240*pow(b2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 288*a2*b2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 
840*k1*pow(a2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 576*b2*c2*d*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 90*b2*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 
576*a2*d*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 45*pow(c2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 2370*a2*c2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 
1185*pow(b2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 7488*a2*b2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 
9765*pow(a2,2)*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 90*b2*d2*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) - 
45*pow(c2,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) + 150*a2*c2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) + 
75*pow(b2,2)*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) - 525*pow(a2,2)*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) - 
288*b2*c2*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 90*b2*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 288*a2*d*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 
288*b2*c2*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 1260*a2*c2*d*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 
288*a2*d2*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 630*d*k1*pow(b2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 
45*pow(c2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 480*a2*c2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 
2592*a2*b2*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 240*pow(b2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 
288*a2*b2*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 840*k1*pow(a2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 
630*a2*c2*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 3456*a2*b2*d*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 
315*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 4200*pow(a2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 
630*a2*c2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 3456*a2*b2*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 
315*pow(b2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 4200*pow(a2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 
1152*a2*b2*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,5) - 3780*d*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,5) + 
1152*a2*b2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,5) - 3780*d*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,5) - 
945*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,5) + 945*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,5) - 
288*b2*c2*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 1260*a2*c2*d*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 288*a2*d2*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 
1260*a2*c2*k1*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 6912*a2*b2*d*k1*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 630*d*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 
630*k1*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 2592*a2*b2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 8400*k1*pow(a2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 
840*pow(a2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 288*b2*c2*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 2520*a2*c2*d*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 
288*a2*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 1260*d*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 12960*a2*b2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 
24360*pow(a2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 288*b2*c2*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,6) + 288*a2*d2*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,6) - 
864*a2*b2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,6) - 288*b2*c2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) + 
1260*a2*c2*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 288*a2*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 1260*a2*c2*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) + 
6912*a2*b2*d*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) + 630*d*pow(b2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 
630*k1*pow(b2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 2592*a2*b2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 
8400*k1*pow(a2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) + 840*pow(a2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 
3456*a2*b2*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,6) - 11340*d*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,6) - 
3456*a2*b2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,6) + 11340*d*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,6) - 
3780*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,6) - 3780*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,6) - 
1260*a2*c2*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 6912*a2*b2*d*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 6912*a2*b2*k1*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 
22680*d*k1*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 630*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 8400*pow(a2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 
1260*a2*c2*pow(M_E,k1*pow(RP,-1))*pow(RP,7) + 13824*a2*b2*d*pow(M_E,k1*pow(RP,-1))*pow(RP,7) + 630*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) + 
42420*pow(a2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 1260*a2*c2*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,7) - 
630*pow(b2,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,7) + 2940*pow(a2,2)*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,7) + 
1260*a2*c2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) - 6912*a2*b2*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) + 
6912*a2*b2*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) - 22680*d*k1*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) + 
630*pow(b2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) + 8400*pow(a2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) - 
11340*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,7) + 11340*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,7) - 
6912*a2*b2*pow(M_E,d*pow(RP,-1))*pow(RP,8) - 22680*d*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,8) - 22680*k1*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,8) + 
6912*a2*b2*pow(M_E,k1*pow(RP,-1))*pow(RP,8) + 45360*d*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,8) + 6912*a2*b2*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,8) - 
6912*a2*b2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,8) + 22680*d*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,8) - 
22680*k1*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,8) - 22680*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,9) + 22680*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,9) - 
22680*pow(a2,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,9) + 22680*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,9)))/2.;

return res;
}

double SDH_Integral_z_2233_case_3( double AP, double RP, double k1, double k2, double a2, double b2, double c2, double d2, double d )
{
double res =(-3*AP*RP*pow(d,-4)*pow(M_E,-(d*pow(RP,-1)))*(RP*(gsl_sf_expint_Ei(-(k1*pow(RP,-1))) - gsl_sf_expint_Ei(-(k2*pow(RP,-1))))*pow(d,3)*pow(d2,2)*(1 + pow(M_E,2*d*pow(RP,-1))) + 
RP*pow(d,3)*pow(d2,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,k1*pow(RP,-1)) - pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1))) + 
2*c2*d2*RP*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1))) + 
RP*pow(d,2)*pow(d2,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1))) + 
4*(-gsl_sf_expint_Ei(-(k1*pow(RP,-1))) + gsl_sf_expint_Ei(-(k2*pow(RP,-1))))*pow(d,2)*pow(d2,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,2) + 
4*pow(d,2)*pow(d2,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,-(k1*pow(RP,-1))) - pow(M_E,-(k2*pow(RP,-1))))*pow(RP,2) - 
2*c2*d2*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,-(k1*pow(RP,-1))) - pow(M_E,-(k2*pow(RP,-1))))*pow(RP,2) + 
2*b2*d2*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2) + 
pow(c2,2)*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2) + 
3*d*pow(d2,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2) + 
8*c2*d2*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2) + 
2*c2*d2*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
8*b2*d2*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
4*pow(c2,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
2*b2*d2*RP*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
RP*pow(c2,2)*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
6*c2*d*d2*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
2*b2*c2*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
2*a2*d2*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
9*d*(gsl_sf_expint_Ei(-(k1*pow(RP,-1))) - gsl_sf_expint_Ei(-(k2*pow(RP,-1))))*pow(d2,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,3) + 
8*c2*d2*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,-(k1*pow(RP,-1))) - pow(M_E,-(k2*pow(RP,-1))))*pow(RP,3) + 
9*d*pow(d2,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,k1*pow(RP,-1)) - pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3) + 
18*c2*d*d2*(1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3) + 
8*b2*d2*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3) + 
4*pow(c2,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
pow(RP,3) + 3*pow(d2,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
pow(RP,3) + 6*c2*d2*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,3) + 
8*b2*c2*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,3) + 
8*a2*d2*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,3) + 
18*b2*d*d2*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,3) + 
9*d*pow(c2,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,3) + 
2*b2*d2*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
RP*pow(c2,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
8*b2*c2*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
8*a2*d2*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
6*b2*d2*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
3*pow(c2,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
8*a2*c2*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
4*pow(b2,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
2*b2*c2*RP*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
2*a2*d2*RP*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
6*b2*d*d2*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
3*d*pow(c2,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
2*a2*c2*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
pow(b2,2)*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
18*b2*c2*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
18*a2*d*d2*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
9*(-gsl_sf_expint_Ei(-(k1*pow(RP,-1))) + gsl_sf_expint_Ei(-(k2*pow(RP,-1))))*pow(d2,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,4) + 
9*pow(d2,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,-(k1*pow(RP,-1))) - pow(M_E,-(k2*pow(RP,-1))))*pow(RP,4) + 
18*c2*d*d2*(1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,k1*pow(RP,-1)) - pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4) + 
18*b2*d*d2*(1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4) + 
9*d*pow(c2,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4) + 
18*c2*d2*(-1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4) + 
18*b2*d2*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,4) + 
9*pow(c2,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,4) + 
18*b2*c2*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,4) + 
18*a2*d*d2*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,4) + 
18*b2*c2*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,4) + 
18*a2*d2*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,4) + 
18*a2*c2*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,4) + 
9*d*pow(b2,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,4) + 
2*b2*c2*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
2*a2*d2*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
8*a2*c2*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
4*pow(b2,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
6*b2*c2*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
6*a2*d2*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
8*a2*b2*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
18*a2*c2*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*
(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
9*pow(b2,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*
(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
2*a2*c2*RP*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
RP*pow(b2,2)*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
6*b2*c2*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
6*a2*d*d2*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
2*a2*b2*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
18*a2*c2*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
9*d*pow(b2,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
18*a2*b2*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
18*c2*d2*(-1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,-(k1*pow(RP,-1))) - pow(M_E,-(k2*pow(RP,-1))))*pow(RP,5) + 
18*b2*d2*(-1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,5) + 
9*pow(c2,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,5) + 
18*b2*c2*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,5) + 
18*a2*d2*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,5) + 
18*a2*c2*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,5) + 
9*pow(b2,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,5) + 
18*a2*b2*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4)))*pow(RP,5) + 
2*a2*c2*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5)) - 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
RP*pow(b2,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5)) - 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
8*a2*b2*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5)) - 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
6*a2*c2*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5)) - 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
3*pow(b2,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5)) - 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
4*pow(a2,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5)) - 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
18*a2*b2*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*
(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5)) - 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
9*pow(a2,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,5)*
(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5)) - 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
2*a2*b2*RP*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))) + 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
6*a2*c2*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))) + 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
3*d*pow(b2,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))) + 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
pow(a2,2)*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))) + 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
18*a2*b2*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))) + 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
9*d*pow(a2,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*
(-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))) + 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
2*a2*b2*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 
720*pow(RP,6)) - pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 
720*k2*pow(RP,5) + 720*pow(RP,6))) + 4*pow(a2,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 
720*pow(RP,6)) - pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 
720*k2*pow(RP,5) + 720*pow(RP,6))) + 6*a2*b2*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 
720*pow(RP,6)) - pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 
720*k2*pow(RP,5) + 720*pow(RP,6))) + 9*pow(a2,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*
(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 
720*pow(RP,6)) - pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 
720*k2*pow(RP,5) + 720*pow(RP,6))) + RP*pow(a2,2)*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 
720*pow(RP,6))) + pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 
720*k2*pow(RP,5) + 720*pow(RP,6))) + 6*a2*b2*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 
720*pow(RP,6))) + pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 
720*k2*pow(RP,5) + 720*pow(RP,6))) + 9*d*pow(a2,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(-(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 
720*pow(RP,6))) + pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 
720*k2*pow(RP,5) + 720*pow(RP,6))) + RP*pow(a2,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(7*RP*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RP,2) + 210*pow(k1,4)*pow(RP,3) + 840*pow(k1,3)*pow(RP,4) + 2520*pow(k1,2)*pow(RP,5) + 
5040*k1*pow(RP,6) + 5040*pow(RP,7)) - pow(M_E,k1*pow(RP,-1))*
(7*RP*pow(k2,6) + pow(k2,7) + 42*pow(k2,5)*pow(RP,2) + 210*pow(k2,4)*pow(RP,3) + 840*pow(k2,3)*pow(RP,4) + 2520*pow(k2,2)*pow(RP,5) + 5040*k2*pow(RP,6) + 
5040*pow(RP,7))) + 3*pow(a2,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(pow(M_E,k2*pow(RP,-1))*(7*RP*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RP,2) + 210*pow(k1,4)*pow(RP,3) + 840*pow(k1,3)*pow(RP,4) + 2520*pow(k1,2)*pow(RP,5) + 
5040*k1*pow(RP,6) + 5040*pow(RP,7)) - pow(M_E,k1*pow(RP,-1))*
(7*RP*pow(k2,6) + pow(k2,7) + 42*pow(k2,5)*pow(RP,2) + 210*pow(k2,4)*pow(RP,3) + 840*pow(k2,3)*pow(RP,4) + 2520*pow(k2,2)*pow(RP,5) + 5040*k2*pow(RP,6) + 
5040*pow(RP,7))) + 3*d*pow(a2,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(7*RP*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RP,2) + 210*pow(k1,4)*pow(RP,3) + 840*pow(k1,3)*pow(RP,4) + 2520*pow(k1,2)*pow(RP,5) + 
5040*k1*pow(RP,6) + 5040*pow(RP,7))) + pow(M_E,k1*pow(RP,-1))*
(7*RP*pow(k2,6) + pow(k2,7) + 42*pow(k2,5)*pow(RP,2) + 210*pow(k2,4)*pow(RP,3) + 840*pow(k2,3)*pow(RP,4) + 2520*pow(k2,2)*pow(RP,5) + 5040*k2*pow(RP,6) + 
5040*pow(RP,7)))))/2.;

return res;
}

/* < X or Y | dzH | X or Y >  done */










/* < Z | dzH | Z >  */
double SDH_Integral_z_44_case_1( double AP, double RP, double k1, double k2, double a2, double b2, double c2, double d2, double d )
{
double res =(-3*AP*RP*pow(d,-4)*pow(M_E,-((d + k1 + k2)*pow(RP,-1)))*(-8*c2*d2*k2*RP*pow(d,3)*pow(M_E,k1*pow(RP,-1)) - 3*k2*RP*pow(d,2)*pow(d2,2)*pow(M_E,k1*pow(RP,-1)) - 
k2*pow(d,3)*pow(d2,2)*pow(M_E,k1*pow(RP,-1)) - 3*RP*pow(d,3)*pow(d2,2)*pow(M_E,k1*pow(RP,-1)) + 8*c2*d2*k1*RP*pow(d,3)*pow(M_E,k2*pow(RP,-1)) + 
3*k1*RP*pow(d,2)*pow(d2,2)*pow(M_E,k2*pow(RP,-1)) + k1*pow(d,3)*pow(d2,2)*pow(M_E,k2*pow(RP,-1)) + 3*RP*pow(d,3)*pow(d2,2)*pow(M_E,k2*pow(RP,-1)) + 
2*RP*gsl_sf_expint_Ei(-(k2*pow(RP,-1)))*pow(d,3)*pow(d2,2)*pow(M_E,(k1 + k2)*pow(RP,-1)) - 
2*RP*gsl_sf_expint_Ei(k2*pow(RP,-1))*pow(d,3)*pow(d2,2)*pow(M_E,(k1 + k2)*pow(RP,-1)) - 8*c2*d2*k1*RP*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) + 
3*k1*RP*pow(d,2)*pow(d2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) + k1*pow(d,3)*pow(d2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) - 
3*RP*pow(d,3)*pow(d2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) + 8*c2*d2*k2*RP*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) - 
3*k2*RP*pow(d,2)*pow(d2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) - k2*pow(d,3)*pow(d2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) + 
3*RP*pow(d,3)*pow(d2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) + 6*c2*d2*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 2*c2*d2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 
10*b2*d2*RP*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 5*RP*pow(c2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 
6*c2*d2*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) + 2*c2*d2*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) - 
10*b2*d2*RP*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) - 5*RP*pow(c2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) + 
6*b2*d2*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 3*RP*pow(c2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 2*b2*d2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 
12*b2*c2*RP*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 12*a2*d2*RP*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + pow(c2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 
6*b2*d2*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + 3*RP*pow(c2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + 
2*b2*d2*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) - 12*b2*c2*RP*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) - 
12*a2*d2*RP*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + pow(c2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + 
6*b2*c2*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 6*a2*d2*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 2*b2*c2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 
2*a2*d2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 14*a2*c2*RP*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 
7*RP*pow(b2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 6*b2*c2*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) + 
6*a2*d2*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) + 2*b2*c2*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) + 
2*a2*d2*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) - 14*a2*c2*RP*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) - 
7*RP*pow(b2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) + 6*a2*c2*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 
3*RP*pow(b2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 2*a2*c2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 
16*a2*b2*RP*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + pow(b2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 
6*a2*c2*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) + 3*RP*pow(b2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) + 
2*a2*c2*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) - 16*a2*b2*RP*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) + 
pow(b2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) + 6*a2*b2*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,6) + 
2*a2*b2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,6) + 9*RP*pow(a2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,6) + 
6*a2*b2*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6) + 2*a2*b2*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6) - 
9*RP*pow(a2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6) + 3*RP*pow(a2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,7) + 
pow(a2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,7) + 3*RP*pow(a2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,7) + 
pow(a2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,7) - 6*c2*d2*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 
2*c2*d2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 10*b2*d2*RP*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 
5*RP*pow(c2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 6*c2*d2*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) - 
2*c2*d2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) + 10*b2*d2*RP*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) + 
5*RP*pow(c2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) - 6*b2*d2*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 
3*RP*pow(c2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 2*b2*d2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 
12*b2*c2*RP*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 12*a2*d2*RP*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - pow(c2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 
6*b2*d2*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - 3*RP*pow(c2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - 
2*b2*d2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) + 12*b2*c2*RP*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) + 
12*a2*d2*RP*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - pow(c2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - 
6*b2*c2*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 6*a2*d2*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 2*b2*c2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 
2*a2*d2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 14*a2*c2*RP*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 
7*RP*pow(b2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 6*b2*c2*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) - 
6*a2*d2*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) - 2*b2*c2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) - 
2*a2*d2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) + 14*a2*c2*RP*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) + 
7*RP*pow(b2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) - 6*a2*c2*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 
3*RP*pow(b2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 2*a2*c2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 
16*a2*b2*RP*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - pow(b2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 
6*a2*c2*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) - 3*RP*pow(b2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) - 
2*a2*c2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) + 16*a2*b2*RP*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) - 
pow(b2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) - 6*a2*b2*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,6) - 
2*a2*b2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,6) - 9*RP*pow(a2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,6) - 
6*a2*b2*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6) - 2*a2*b2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6) + 
9*RP*pow(a2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6) - 3*RP*pow(a2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,7) - 
pow(a2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,7) - 3*RP*pow(a2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,7) - 
pow(a2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,7) - 28*c2*d2*k2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 
12*c2*d2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 24*b2*d2*k2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 
12*k2*pow(c2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 6*d*k2*pow(d2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 
11*pow(d,2)*pow(d2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 28*c2*d2*k1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 12*c2*d2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 
24*b2*d2*k1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 12*k1*pow(c2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 
6*d*k1*pow(d2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 11*pow(d,2)*pow(d2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 
8*gsl_sf_expint_Ei(-(k2*pow(RP,-1)))*pow(d,2)*pow(d2,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*pow(RP,2) - 
8*gsl_sf_expint_Ei(k2*pow(RP,-1))*pow(d,2)*pow(d2,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*pow(RP,2) - 28*c2*d2*k1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) + 
12*c2*d2*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) + 24*b2*d2*k1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) + 
12*k1*pow(c2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) + 6*d*k1*pow(d2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) - 
11*pow(d,2)*pow(d2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) + 28*c2*d2*k2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) - 
12*c2*d2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) - 24*b2*d2*k2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) - 
12*k2*pow(c2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) - 6*d*k2*pow(d2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) + 
11*pow(d,2)*pow(d2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) + 12*c2*d*d2*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
34*b2*d2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 17*pow(c2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
40*b2*c2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 40*a2*d2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
12*c2*d*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 34*b2*d2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
17*pow(c2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 40*b2*c2*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
40*a2*d2*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 12*b2*d*d2*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
6*d*pow(c2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 40*b2*c2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
40*a2*d2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 60*a2*c2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
30*pow(b2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 12*b2*d*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
6*d*pow(c2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 40*b2*c2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
40*a2*d2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 60*a2*c2*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
30*pow(b2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 12*b2*c2*d*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 
12*a2*d*d2*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 46*a2*c2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 
23*pow(b2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 84*a2*b2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 
12*b2*c2*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 12*a2*d*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
46*a2*c2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 23*pow(b2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 
84*a2*b2*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 12*a2*c2*d*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 
6*d*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 52*a2*b2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 
56*pow(a2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 12*a2*c2*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 
6*d*pow(b2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 52*a2*b2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 
56*pow(a2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 12*a2*b2*d*pow(M_E,k2*pow(RP,-1))*pow(k1,6)*pow(RP,2) + 
29*pow(a2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,6)*pow(RP,2) + 12*a2*b2*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6)*pow(RP,2) - 
29*pow(a2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6)*pow(RP,2) + 6*d*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,7)*pow(RP,2) + 
6*d*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,7)*pow(RP,2) - 12*c2*d*d2*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
34*b2*d2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 17*pow(c2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
40*b2*c2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 40*a2*d2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
12*c2*d*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 34*b2*d2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
17*pow(c2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 40*b2*c2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
40*a2*d2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 12*b2*d*d2*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
6*d*pow(c2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 40*b2*c2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
40*a2*d2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 60*a2*c2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
30*pow(b2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 12*b2*d*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
6*d*pow(c2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 40*b2*c2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
40*a2*d2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 60*a2*c2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
30*pow(b2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 12*b2*c2*d*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 
12*a2*d*d2*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 46*a2*c2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 
23*pow(b2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 84*a2*b2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 
12*b2*c2*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 12*a2*d*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
46*a2*c2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 23*pow(b2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 
84*a2*b2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 12*a2*c2*d*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 
6*d*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 52*a2*b2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 
56*pow(a2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 12*a2*c2*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 
6*d*pow(b2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 52*a2*b2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 
56*pow(a2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 12*a2*b2*d*pow(M_E,k1*pow(RP,-1))*pow(k2,6)*pow(RP,2) - 
29*pow(a2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,6)*pow(RP,2) - 12*a2*b2*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6)*pow(RP,2) + 
29*pow(a2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6)*pow(RP,2) - 6*d*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,7)*pow(RP,2) - 
6*d*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,7)*pow(RP,2) - 60*c2*d*d2*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 
44*c2*d2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 84*b2*d2*k2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 
42*k2*pow(c2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 24*b2*d2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 
80*b2*c2*k2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 80*a2*d2*k2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 
12*pow(c2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 24*d*pow(d2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 6*k2*pow(d2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 
60*c2*d*d2*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 44*c2*d2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 84*b2*d2*k1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 
42*k1*pow(c2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 24*b2*d2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 
80*b2*c2*k1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 80*a2*d2*k1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 
12*pow(c2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 24*d*pow(d2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 6*k1*pow(d2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 
18*d*gsl_sf_expint_Ei(-(k2*pow(RP,-1)))*pow(d2,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*pow(RP,3) - 
18*d*gsl_sf_expint_Ei(k2*pow(RP,-1))*pow(d2,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*pow(RP,3) - 60*c2*d*d2*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) + 
44*c2*d2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) + 84*b2*d2*k1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) + 
42*k1*pow(c2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 24*b2*d2*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 
80*b2*c2*k1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 80*a2*d2*k1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 
12*pow(c2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 24*d*pow(d2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) + 
6*k1*pow(d2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) + 60*c2*d*d2*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) - 
44*c2*d2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) - 84*b2*d2*k2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) - 
42*k2*pow(c2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 24*b2*d2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 
80*b2*c2*k2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 80*a2*d2*k2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 
12*pow(c2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 24*d*pow(d2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) - 
6*k2*pow(d2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 12*c2*d2*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
72*b2*d*d2*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 36*d*pow(c2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
136*b2*c2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 136*a2*d2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
180*a2*c2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 90*pow(b2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
12*c2*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 72*b2*d*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
36*d*pow(c2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 136*b2*c2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
136*a2*d2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 180*a2*c2*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
90*pow(b2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 84*b2*c2*d*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
12*b2*d2*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 84*a2*d*d2*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
6*pow(c2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 200*a2*c2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
100*pow(b2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 336*a2*b2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
84*b2*c2*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 12*b2*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
84*a2*d*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 6*pow(c2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
200*a2*c2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 100*pow(b2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
336*a2*b2*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 12*b2*c2*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 
96*a2*c2*d*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 12*a2*d2*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 
48*d*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 276*a2*b2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 
280*pow(a2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 12*b2*c2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 
96*a2*c2*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 12*a2*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 
48*d*pow(b2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 276*a2*b2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 
280*pow(a2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 12*a2*c2*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,3) + 
108*a2*b2*d*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,3) + 6*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,3) + 
182*pow(a2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,3) + 12*a2*c2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,3) - 
108*a2*b2*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,3) + 6*pow(b2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,3) + 
182*pow(a2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,3) + 12*a2*b2*pow(M_E,k2*pow(RP,-1))*pow(k1,6)*pow(RP,3) + 
60*d*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,6)*pow(RP,3) + 12*a2*b2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6)*pow(RP,3) - 
60*d*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6)*pow(RP,3) + 6*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,7)*pow(RP,3) + 
6*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,7)*pow(RP,3) - 12*c2*d2*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
72*b2*d*d2*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 36*d*pow(c2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
136*b2*c2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 136*a2*d2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
180*a2*c2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 90*pow(b2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
12*c2*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 72*b2*d*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
36*d*pow(c2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 136*b2*c2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
136*a2*d2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 180*a2*c2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
90*pow(b2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 84*b2*c2*d*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 
12*b2*d2*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 84*a2*d*d2*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 
6*pow(c2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 200*a2*c2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 
100*pow(b2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 336*a2*b2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
84*b2*c2*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 12*b2*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
84*a2*d*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 6*pow(c2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 
200*a2*c2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 100*pow(b2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
336*a2*b2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 12*b2*c2*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 
96*a2*c2*d*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 12*a2*d2*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 
48*d*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 276*a2*b2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 
280*pow(a2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 12*b2*c2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 
96*a2*c2*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 12*a2*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 
48*d*pow(b2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 276*a2*b2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 
280*pow(a2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 12*a2*c2*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 
108*a2*b2*d*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 6*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 
182*pow(a2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 12*a2*c2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,3) + 
108*a2*b2*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 6*pow(b2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 
182*pow(a2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 12*a2*b2*pow(M_E,k1*pow(RP,-1))*pow(k2,6)*pow(RP,3) - 
60*d*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,6)*pow(RP,3) - 12*a2*b2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6)*pow(RP,3) + 
60*d*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6)*pow(RP,3) - 6*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,7)*pow(RP,3) - 
6*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,7)*pow(RP,3) - 
2*RP*gsl_sf_expint_Ei(-(k1*pow(RP,-1)))*pow(d2,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*(4*RP*pow(d,2) + pow(d,3) + 9*d*pow(RP,2) + 9*pow(RP,3)) + 
2*RP*gsl_sf_expint_Ei(k1*pow(RP,-1))*pow(d2,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*(4*RP*pow(d,2) + pow(d,3) + 9*d*pow(RP,2) + 9*pow(RP,3)) - 
96*c2*d*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 60*c2*d2*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 180*b2*d*d2*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 
90*d*k2*pow(c2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 84*b2*d2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 272*b2*c2*k2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 
272*a2*d2*k2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 42*pow(c2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 
80*b2*c2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 80*a2*d2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 360*a2*c2*k2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 
180*k2*pow(b2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 24*pow(d2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 96*c2*d*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
60*c2*d2*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 180*b2*d*d2*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 90*d*k1*pow(c2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
84*b2*d2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 272*b2*c2*k1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
272*a2*d2*k1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 42*pow(c2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
80*b2*c2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 80*a2*d2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 360*a2*c2*k1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
180*k1*pow(b2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 24*pow(d2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
18*gsl_sf_expint_Ei(-(k2*pow(RP,-1)))*pow(d2,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*pow(RP,4) - 
18*gsl_sf_expint_Ei(k2*pow(RP,-1))*pow(d2,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*pow(RP,4) + 96*c2*d*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 
60*c2*d2*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 180*b2*d*d2*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 
90*d*k1*pow(c2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 84*b2*d2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 
272*b2*c2*k1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 272*a2*d2*k1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 
42*pow(c2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 80*b2*c2*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 
80*a2*d2*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 360*a2*c2*k1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 
180*k1*pow(b2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 24*pow(d2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 
96*c2*d*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 60*c2*d2*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 
180*b2*d*d2*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 90*d*k2*pow(c2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 
84*b2*d2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 272*b2*c2*k2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 
272*a2*d2*k2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 42*pow(c2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 
80*b2*c2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 80*a2*d2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 
360*a2*c2*k2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 180*k2*pow(b2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 
24*pow(d2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 288*b2*c2*d*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
72*b2*d2*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 288*a2*d*d2*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
36*pow(c2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 600*a2*c2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
300*pow(b2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 1008*a2*b2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
288*b2*c2*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 72*b2*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
288*a2*d*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 36*pow(c2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
600*a2*c2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 300*pow(b2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
1008*a2*b2*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 84*b2*c2*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 
420*a2*c2*d*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 84*a2*d2*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 
210*d*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 1104*a2*b2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 
1120*pow(a2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 84*b2*c2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 
420*a2*c2*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 84*a2*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 
210*d*pow(b2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 1104*a2*b2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 
1120*pow(a2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 96*a2*c2*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,4) + 
576*a2*b2*d*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,4) + 48*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,4) + 
910*pow(a2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,4) - 96*a2*c2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,4) + 
576*a2*b2*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,4) - 48*pow(b2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,4) - 
910*pow(a2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,4) + 108*a2*b2*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,4) + 
378*d*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,4) - 108*a2*b2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,4) + 
378*d*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,4) + 60*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,6)*pow(RP,4) - 
60*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6)*pow(RP,4) - 288*b2*c2*d*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
72*b2*d2*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 288*a2*d*d2*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
36*pow(c2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 600*a2*c2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
300*pow(b2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 1008*a2*b2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
288*b2*c2*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 72*b2*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
288*a2*d*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 36*pow(c2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 
600*a2*c2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 300*pow(b2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
1008*a2*b2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 84*b2*c2*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 
420*a2*c2*d*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 84*a2*d2*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 
210*d*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 1104*a2*b2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 
1120*pow(a2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 84*b2*c2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 
420*a2*c2*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 84*a2*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 
210*d*pow(b2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 1104*a2*b2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 
1120*pow(a2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 96*a2*c2*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,4) - 
576*a2*b2*d*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,4) - 48*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,4) - 
910*pow(a2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,4) + 96*a2*c2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,4) - 
576*a2*b2*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,4) + 48*pow(b2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,4) + 
910*pow(a2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,4) - 108*a2*b2*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,4) - 
378*d*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,4) + 108*a2*b2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,4) - 
378*d*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,4) - 60*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,6)*pow(RP,4) + 
60*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6)*pow(RP,4) - 96*c2*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 180*b2*d*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 
576*b2*c2*d*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 180*b2*d2*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 576*a2*d*d2*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 
90*d*pow(c2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 90*k2*pow(c2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 272*b2*c2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 
272*a2*d2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 1200*a2*c2*k2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 
600*k2*pow(b2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 360*a2*c2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 
2016*a2*b2*k2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 180*pow(b2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 96*c2*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
180*b2*d*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 576*b2*c2*d*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 180*b2*d2*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
576*a2*d*d2*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 90*d*pow(c2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 90*k1*pow(c2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
272*b2*c2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 272*a2*d2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 1200*a2*c2*k1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
600*k1*pow(b2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 360*a2*c2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
2016*a2*b2*k1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 180*pow(b2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
96*c2*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 180*b2*d*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 
576*b2*c2*d*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) + 180*b2*d2*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 
576*a2*d*d2*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 90*d*pow(c2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) + 
90*k1*pow(c2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) + 272*b2*c2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) + 
272*a2*d2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) + 1200*a2*c2*k1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) + 
600*k1*pow(b2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 360*a2*c2*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 
2016*a2*b2*k1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 180*pow(b2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 
96*c2*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 180*b2*d*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 
576*b2*c2*d*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) - 180*b2*d2*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 
576*a2*d*d2*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 90*d*pow(c2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) - 
90*k2*pow(c2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) - 272*b2*c2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) - 
272*a2*d2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) - 1200*a2*c2*k2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) - 
600*k2*pow(b2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 360*a2*c2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 
2016*a2*b2*k2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 180*pow(b2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 
288*b2*c2*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 1260*a2*c2*d*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 
288*a2*d2*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 630*d*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 
3312*a2*b2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 3360*pow(a2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 
288*b2*c2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 1260*a2*c2*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 
288*a2*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 630*d*pow(b2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 
3312*a2*b2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 3360*pow(a2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 
420*a2*c2*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,5) + 2304*a2*b2*d*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,5) + 
210*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,5) + 3640*pow(a2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,5) + 
420*a2*c2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,5) - 2304*a2*b2*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,5) + 
210*pow(b2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,5) + 3640*pow(a2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,5) + 
576*a2*b2*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,5) + 1890*d*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,5) + 
576*a2*b2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,5) - 1890*d*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,5) + 
378*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,5) + 378*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,5) - 
288*b2*c2*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 1260*a2*c2*d*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 
288*a2*d2*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 630*d*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 
3312*a2*b2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 3360*pow(a2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 
288*b2*c2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 1260*a2*c2*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 
288*a2*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 630*d*pow(b2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 
3312*a2*b2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 3360*pow(a2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 
420*a2*c2*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 2304*a2*b2*d*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 
210*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 3640*pow(a2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 
420*a2*c2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,5) + 2304*a2*b2*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 
210*pow(b2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 3640*pow(a2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 
576*a2*b2*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,5) - 1890*d*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,5) - 
576*a2*b2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,5) + 1890*d*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,5) - 
378*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,5) - 378*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,5) - 
576*b2*c2*d*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 180*b2*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 576*a2*d*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 
576*b2*c2*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 2520*a2*c2*d*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 576*a2*d2*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 
1260*d*k2*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 90*pow(c2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 1200*a2*c2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 
6624*a2*b2*k2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 600*pow(b2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 
2016*a2*b2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 6720*k2*pow(a2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 576*b2*c2*d*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 
180*b2*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 576*a2*d*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 576*b2*c2*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 
2520*a2*c2*d*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 576*a2*d2*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 1260*d*k1*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 
90*pow(c2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 1200*a2*c2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 6624*a2*b2*k1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 
600*pow(b2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 2016*a2*b2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 
6720*k1*pow(a2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 576*b2*c2*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 
180*b2*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) + 576*a2*d*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 
576*b2*c2*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) + 2520*a2*c2*d*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 
576*a2*d2*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) + 1260*d*k1*pow(b2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 
90*pow(c2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 1200*a2*c2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 
6624*a2*b2*k1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 600*pow(b2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) + 
2016*a2*b2*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) + 6720*k1*pow(a2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 
576*b2*c2*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 180*b2*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) - 576*a2*d*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 
576*b2*c2*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) - 2520*a2*c2*d*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 
576*a2*d2*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) - 1260*d*k2*pow(b2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 
90*pow(c2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 1200*a2*c2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 
6624*a2*b2*k2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 600*pow(b2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) - 
2016*a2*b2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) - 6720*k2*pow(a2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 
1260*a2*c2*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,6) + 6912*a2*b2*d*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,6) + 
630*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,6) + 10920*pow(a2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,6) - 
1260*a2*c2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,6) + 6912*a2*b2*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,6) - 
630*pow(b2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,6) - 10920*pow(a2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,6) + 
2304*a2*b2*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,6) + 7560*d*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,6) - 
2304*a2*b2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,6) + 7560*d*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,6) + 
1890*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,6) - 1890*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,6) - 
1260*a2*c2*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,6) - 6912*a2*b2*d*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,6) - 
630*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,6) - 10920*pow(a2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,6) + 
1260*a2*c2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,6) - 6912*a2*b2*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,6) + 
630*pow(b2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,6) + 10920*pow(a2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,6) - 
2304*a2*b2*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,6) - 7560*d*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,6) + 
2304*a2*b2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,6) - 7560*d*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,6) - 
1890*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,6) + 1890*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,6) - 
576*b2*c2*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 2520*a2*c2*d*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 576*a2*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 
2520*a2*c2*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 13824*a2*b2*d*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 1260*d*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 
1260*k2*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 6624*a2*b2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 
21840*k2*pow(a2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 6720*pow(a2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) + 
576*b2*c2*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 2520*a2*c2*d*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 576*a2*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 
2520*a2*c2*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 13824*a2*b2*d*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 1260*d*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 
1260*k1*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 6624*a2*b2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 
21840*k1*pow(a2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 6720*pow(a2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 
576*b2*c2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) - 2520*a2*c2*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) + 576*a2*d2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) + 
2520*a2*c2*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) - 13824*a2*b2*d*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) - 
1260*d*pow(b2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) + 1260*k1*pow(b2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) + 
6624*a2*b2*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) + 21840*k1*pow(a2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) - 
6720*pow(a2,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) - 576*b2*c2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) + 
2520*a2*c2*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) - 576*a2*d2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) - 
2520*a2*c2*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) + 13824*a2*b2*d*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) + 
1260*d*pow(b2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) - 1260*k2*pow(b2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) - 
6624*a2*b2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) - 21840*k2*pow(a2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) + 
6720*pow(a2,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) + 6912*a2*b2*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,7) + 
22680*d*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,7) + 6912*a2*b2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,7) - 
22680*d*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,7) + 7560*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,7) + 
7560*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,7) - 6912*a2*b2*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,7) - 
22680*d*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,7) - 6912*a2*b2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,7) + 
22680*d*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,7) - 7560*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,7) - 
7560*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,7) - 2520*a2*c2*pow(M_E,k1*pow(RP,-1))*pow(RP,8) - 13824*a2*b2*d*pow(M_E,k1*pow(RP,-1))*pow(RP,8) - 
13824*a2*b2*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,8) - 45360*d*k2*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,8) - 1260*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,8) - 
21840*pow(a2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,8) + 2520*a2*c2*pow(M_E,k2*pow(RP,-1))*pow(RP,8) + 13824*a2*b2*d*pow(M_E,k2*pow(RP,-1))*pow(RP,8) + 
13824*a2*b2*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,8) + 45360*d*k1*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,8) + 1260*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,8) + 
21840*pow(a2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,8) - 2520*a2*c2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,8) + 
13824*a2*b2*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,8) - 13824*a2*b2*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,8) + 
45360*d*k1*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,8) - 1260*pow(b2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,8) - 
21840*pow(a2,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,8) + 2520*a2*c2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,8) - 
13824*a2*b2*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,8) + 13824*a2*b2*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,8) - 
45360*d*k2*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,8) + 1260*pow(b2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,8) + 
21840*pow(a2,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,8) + 22680*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,8) - 
22680*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,8) - 22680*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,8) + 
22680*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,8) - 13824*a2*b2*pow(M_E,k1*pow(RP,-1))*pow(RP,9) - 
45360*d*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,9) - 45360*k2*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,9) + 13824*a2*b2*pow(M_E,k2*pow(RP,-1))*pow(RP,9) + 
45360*d*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,9) + 45360*k1*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,9) + 13824*a2*b2*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,9) - 
45360*d*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,9) + 45360*k1*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,9) - 
13824*a2*b2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,9) + 45360*d*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,9) - 
45360*k2*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,9) - 45360*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,10) + 
45360*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,10) - 45360*pow(a2,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,10) + 
45360*pow(a2,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,10)))/2.;

return res;
}

double SDH_Integral_z_44_case_2_sub_1( double AP, double RP, double k1, double k2, double a2, double b2, double c2, double d2, double d )
{
double res =(3*AP*RP*pow(d,-4)*pow(M_E,-((2*d + k2)*pow(RP,-1)))*(8*c2*d2*k2*RP*pow(d,3)*pow(M_E,d*pow(RP,-1)) + 3*k2*RP*pow(d,2)*pow(d2,2)*pow(M_E,d*pow(RP,-1)) + 
k2*pow(d,3)*pow(d2,2)*pow(M_E,d*pow(RP,-1)) + 3*RP*pow(d,3)*pow(d2,2)*pow(M_E,d*pow(RP,-1)) + 8*c2*d2*k2*RP*pow(d,3)*pow(M_E,3*d*pow(RP,-1)) - 
3*k2*RP*pow(d,2)*pow(d2,2)*pow(M_E,3*d*pow(RP,-1)) + k2*pow(d,3)*pow(d2,2)*pow(M_E,3*d*pow(RP,-1)) + 3*RP*pow(d,3)*pow(d2,2)*pow(M_E,3*d*pow(RP,-1)) - 
14*c2*d2*RP*pow(d,4)*pow(M_E,k2*pow(RP,-1)) - 2*c2*d2*pow(d,5)*pow(M_E,k2*pow(RP,-1)) - 16*b2*d2*RP*pow(d,5)*pow(M_E,k2*pow(RP,-1)) - 
8*RP*pow(c2,2)*pow(d,5)*pow(M_E,k2*pow(RP,-1)) - 2*b2*d2*pow(d,6)*pow(M_E,k2*pow(RP,-1)) - 18*b2*c2*RP*pow(d,6)*pow(M_E,k2*pow(RP,-1)) - 
18*a2*d2*RP*pow(d,6)*pow(M_E,k2*pow(RP,-1)) - pow(c2,2)*pow(d,6)*pow(M_E,k2*pow(RP,-1)) - 2*b2*c2*pow(d,7)*pow(M_E,k2*pow(RP,-1)) - 
2*a2*d2*pow(d,7)*pow(M_E,k2*pow(RP,-1)) - 20*a2*c2*RP*pow(d,7)*pow(M_E,k2*pow(RP,-1)) - 10*RP*pow(b2,2)*pow(d,7)*pow(M_E,k2*pow(RP,-1)) - 
2*a2*c2*pow(d,8)*pow(M_E,k2*pow(RP,-1)) - 22*a2*b2*RP*pow(d,8)*pow(M_E,k2*pow(RP,-1)) - pow(b2,2)*pow(d,8)*pow(M_E,k2*pow(RP,-1)) - 
2*a2*b2*pow(d,9)*pow(M_E,k2*pow(RP,-1)) - 12*RP*pow(a2,2)*pow(d,9)*pow(M_E,k2*pow(RP,-1)) - pow(a2,2)*pow(d,10)*pow(M_E,k2*pow(RP,-1)) - 
6*RP*pow(d,3)*pow(d2,2)*pow(M_E,k2*pow(RP,-1)) - pow(d,4)*pow(d2,2)*pow(M_E,k2*pow(RP,-1)) - 2*c2*d2*RP*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1)) - 
2*c2*d2*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1)) - 4*b2*d2*RP*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1)) - 2*RP*pow(c2,2)*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1)) - 
2*b2*d2*pow(d,6)*pow(M_E,(2*d + k2)*pow(RP,-1)) - 6*b2*c2*RP*pow(d,6)*pow(M_E,(2*d + k2)*pow(RP,-1)) - 6*a2*d2*RP*pow(d,6)*pow(M_E,(2*d + k2)*pow(RP,-1)) - 
pow(c2,2)*pow(d,6)*pow(M_E,(2*d + k2)*pow(RP,-1)) - 2*b2*c2*pow(d,7)*pow(M_E,(2*d + k2)*pow(RP,-1)) - 2*a2*d2*pow(d,7)*pow(M_E,(2*d + k2)*pow(RP,-1)) - 
8*a2*c2*RP*pow(d,7)*pow(M_E,(2*d + k2)*pow(RP,-1)) - 4*RP*pow(b2,2)*pow(d,7)*pow(M_E,(2*d + k2)*pow(RP,-1)) - 2*a2*c2*pow(d,8)*pow(M_E,(2*d + k2)*pow(RP,-1)) - 
10*a2*b2*RP*pow(d,8)*pow(M_E,(2*d + k2)*pow(RP,-1)) - pow(b2,2)*pow(d,8)*pow(M_E,(2*d + k2)*pow(RP,-1)) - 2*a2*b2*pow(d,9)*pow(M_E,(2*d + k2)*pow(RP,-1)) - 
6*RP*pow(a2,2)*pow(d,9)*pow(M_E,(2*d + k2)*pow(RP,-1)) - pow(a2,2)*pow(d,10)*pow(M_E,(2*d + k2)*pow(RP,-1)) - pow(d,4)*pow(d2,2)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 
6*c2*d2*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2) + 2*c2*d2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2) + 10*b2*d2*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2) + 
5*RP*pow(c2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2) - 6*c2*d2*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) + 
2*c2*d2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) + 10*b2*d2*RP*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) + 
5*RP*pow(c2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) + 6*b2*d2*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 
3*RP*pow(c2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 2*b2*d2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 12*b2*c2*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 
12*a2*d2*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + pow(c2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,3) - 6*b2*d2*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) - 
3*RP*pow(c2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) + 2*b2*d2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) + 
12*b2*c2*RP*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) + 12*a2*d2*RP*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) + 
pow(c2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) + 6*b2*c2*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 6*a2*d2*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 
2*b2*c2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 2*a2*d2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 14*a2*c2*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 
7*RP*pow(b2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,4) - 6*b2*c2*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) - 
6*a2*d2*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) + 2*b2*c2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) + 2*a2*d2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) + 
14*a2*c2*RP*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) + 7*RP*pow(b2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) + 
6*a2*c2*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 3*RP*pow(b2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 2*a2*c2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 
16*a2*b2*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,5) + pow(b2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,5) - 6*a2*c2*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) - 
3*RP*pow(b2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) + 2*a2*c2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) + 
16*a2*b2*RP*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) + pow(b2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) + 6*a2*b2*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,6) + 
2*a2*b2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,6) + 9*RP*pow(a2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,6) - 6*a2*b2*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6) + 
2*a2*b2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6) + 9*RP*pow(a2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6) + 
3*RP*pow(a2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,7) + pow(a2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,7) - 
3*RP*pow(a2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,7) + pow(a2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,7) + 
28*c2*d2*k2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 12*c2*d2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 24*b2*d2*k2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 
12*k2*pow(c2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 6*d*k2*pow(d2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 11*pow(d,2)*pow(d2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 
28*c2*d2*k2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) + 12*c2*d2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) + 
24*b2*d2*k2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) + 12*k2*pow(c2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) + 
6*d*k2*pow(d2,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) - 11*pow(d,2)*pow(d2,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) - 52*c2*d2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 
70*b2*d2*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 35*pow(c2,2)*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 92*b2*c2*pow(d,5)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 
92*a2*d2*pow(d,5)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 118*a2*c2*pow(d,6)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 59*pow(b2,2)*pow(d,6)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 
148*a2*b2*pow(d,7)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 91*pow(a2,2)*pow(d,8)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 
17*pow(d,2)*pow(d2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 4*c2*d2*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) - 
2*b2*d2*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) - pow(c2,2)*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) - 
12*b2*c2*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) - 12*a2*d2*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) - 
26*a2*c2*pow(d,6)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) - 13*pow(b2,2)*pow(d,6)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) - 
44*a2*b2*pow(d,7)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) - 33*pow(a2,2)*pow(d,8)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 
5*pow(d,2)*pow(d2,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 12*c2*d*d2*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
34*b2*d2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 17*pow(c2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
40*b2*c2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 40*a2*d2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
12*c2*d*d2*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 34*b2*d2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
17*pow(c2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 40*b2*c2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
40*a2*d2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 12*b2*d*d2*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
6*d*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 40*b2*c2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
40*a2*d2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 60*a2*c2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
30*pow(b2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 12*b2*d*d2*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
6*d*pow(c2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 40*b2*c2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
40*a2*d2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 60*a2*c2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
30*pow(b2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 12*b2*c2*d*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
12*a2*d*d2*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 46*a2*c2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
23*pow(b2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 84*a2*b2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
12*b2*c2*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 12*a2*d*d2*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 
46*a2*c2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 23*pow(b2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
84*a2*b2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 12*a2*c2*d*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 
6*d*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 52*a2*b2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 
56*pow(a2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 12*a2*c2*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 
6*d*pow(b2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 52*a2*b2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 
56*pow(a2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 12*a2*b2*d*pow(M_E,d*pow(RP,-1))*pow(k2,6)*pow(RP,2) + 
29*pow(a2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,6)*pow(RP,2) + 12*a2*b2*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,6)*pow(RP,2) - 
29*pow(a2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6)*pow(RP,2) + 6*d*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,7)*pow(RP,2) + 
6*d*pow(a2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,7)*pow(RP,2) + 60*c2*d*d2*k2*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 44*c2*d2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 
84*b2*d2*k2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 42*k2*pow(c2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 24*b2*d2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 
80*b2*c2*k2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 80*a2*d2*k2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 12*pow(c2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 
24*d*pow(d2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 6*k2*pow(d2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 60*c2*d*d2*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 
44*c2*d2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 84*b2*d2*k2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 
42*k2*pow(c2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) + 24*b2*d2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) + 
80*b2*c2*k2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) + 80*a2*d2*k2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) + 
12*pow(c2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) + 24*d*pow(d2,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 6*k2*pow(d2,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 
116*c2*d2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 192*b2*d2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 96*pow(c2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 
312*b2*c2*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 312*a2*d2*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 488*a2*c2*pow(d,5)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 
244*pow(b2,2)*pow(d,5)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 732*a2*b2*pow(d,6)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 
528*pow(a2,2)*pow(d,7)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 30*d*pow(d2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 
4*c2*d2*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) - 16*b2*c2*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) - 
16*a2*d2*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) - 64*a2*c2*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) - 
32*pow(b2,2)*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) - 156*a2*b2*pow(d,6)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) - 
152*pow(a2,2)*pow(d,7)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) - 18*d*pow(d2,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) + 
12*c2*d2*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 72*b2*d*d2*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
36*d*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 136*b2*c2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
136*a2*d2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 180*a2*c2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
90*pow(b2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 12*c2*d2*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
72*b2*d*d2*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 36*d*pow(c2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
136*b2*c2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 136*a2*d2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
180*a2*c2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 90*pow(b2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
84*b2*c2*d*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 12*b2*d2*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 84*a2*d*d2*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
6*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 200*a2*c2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
100*pow(b2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 336*a2*b2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
84*b2*c2*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 12*b2*d2*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
84*a2*d*d2*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 6*pow(c2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 
200*a2*c2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 100*pow(b2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
336*a2*b2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 12*b2*c2*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 
96*a2*c2*d*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 12*a2*d2*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 
48*d*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 276*a2*b2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 
280*pow(a2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 12*b2*c2*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 
96*a2*c2*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 12*a2*d2*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 
48*d*pow(b2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 276*a2*b2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 
280*pow(a2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 12*a2*c2*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,3) + 
108*a2*b2*d*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,3) + 6*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,3) + 
182*pow(a2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 12*a2*c2*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,3) + 
108*a2*b2*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 6*pow(b2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 
182*pow(a2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,3) + 12*a2*b2*pow(M_E,d*pow(RP,-1))*pow(k2,6)*pow(RP,3) + 
60*d*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,6)*pow(RP,3) - 12*a2*b2*pow(M_E,3*d*pow(RP,-1))*pow(k2,6)*pow(RP,3) + 
60*d*pow(a2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6)*pow(RP,3) + 6*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,7)*pow(RP,3) - 
6*pow(a2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,7)*pow(RP,3) + 2*RP*gsl_sf_expint_Ei(-(d*pow(RP,-1)))*pow(d2,2)*pow(M_E,(d + k2)*pow(RP,-1))*
(-4*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1))) + pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1))) + 9*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,2) - 
9*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,3)) - 2*RP*gsl_sf_expint_Ei(-(k2*pow(RP,-1)))*pow(d2,2)*pow(M_E,(d + k2)*pow(RP,-1))*
(-4*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1))) + pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1))) + 9*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,2) - 
9*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,3)) + 96*c2*d*d2*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 60*c2*d2*k2*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 
180*b2*d*d2*k2*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 90*d*k2*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 84*b2*d2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 
272*b2*c2*k2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 272*a2*d2*k2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 
42*pow(c2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 80*b2*c2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 80*a2*d2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 
360*a2*c2*k2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 180*k2*pow(b2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 24*pow(d2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 
96*c2*d*d2*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 60*c2*d2*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 180*b2*d*d2*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 
90*d*k2*pow(c2,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 84*b2*d2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 272*b2*c2*k2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 
272*a2*d2*k2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 42*pow(c2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 
80*b2*c2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 80*a2*d2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 360*a2*c2*k2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 
180*k2*pow(b2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 24*pow(d2,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 156*c2*d*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 
336*b2*d2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 168*pow(c2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 724*b2*c2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 
724*a2*d2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 1476*a2*c2*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 
738*pow(b2,2)*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 2796*a2*b2*pow(d,5)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 
2468*pow(a2,2)*pow(d,6)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 24*pow(d2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 36*c2*d*d2*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) - 
24*b2*d2*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) - 12*pow(c2,2)*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) - 
12*b2*c2*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) - 12*a2*d2*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) - 
84*a2*c2*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) - 42*pow(b2,2)*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) - 
372*a2*b2*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) - 528*pow(a2,2)*pow(d,6)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) + 
24*pow(d2,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) + 288*b2*c2*d*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 72*b2*d2*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 
288*a2*d*d2*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 36*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 
600*a2*c2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 300*pow(b2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 
1008*a2*b2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 288*b2*c2*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
72*b2*d2*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 288*a2*d*d2*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
36*pow(c2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 600*a2*c2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
300*pow(b2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 1008*a2*b2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 
84*b2*c2*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 420*a2*c2*d*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 84*a2*d2*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 
210*d*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 1104*a2*b2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 
1120*pow(a2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 84*b2*c2*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 
420*a2*c2*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 84*a2*d2*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 
210*d*pow(b2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 1104*a2*b2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 
1120*pow(a2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 96*a2*c2*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,4) + 
576*a2*b2*d*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,4) + 48*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,4) + 
910*pow(a2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,4) - 96*a2*c2*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,4) + 
576*a2*b2*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,4) - 48*pow(b2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,4) - 
910*pow(a2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,4) + 108*a2*b2*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,4) + 
378*d*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,4) - 108*a2*b2*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,4) + 
378*d*pow(a2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,4) + 60*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,6)*pow(RP,4) - 
60*pow(a2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6)*pow(RP,4) + 96*c2*d2*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 180*b2*d*d2*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 
576*b2*c2*d*k2*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 180*b2*d2*k2*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 576*a2*d*d2*k2*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 
90*d*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 90*k2*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 272*b2*c2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 
272*a2*d2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 1200*a2*c2*k2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 
600*k2*pow(b2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 360*a2*c2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 
2016*a2*b2*k2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 180*pow(b2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 96*c2*d2*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 
180*b2*d*d2*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 576*b2*c2*d*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 180*b2*d2*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 
576*a2*d*d2*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 90*d*pow(c2,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 90*k2*pow(c2,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 
272*b2*c2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 272*a2*d2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 
1200*a2*c2*k2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 600*k2*pow(b2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 
360*a2*c2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 2016*a2*b2*k2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 
180*pow(b2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 96*c2*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 360*b2*d*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 
180*d*pow(c2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 1136*b2*c2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 1136*a2*d2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 
3240*a2*c2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 1620*pow(b2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 
8208*a2*b2*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 9268*pow(a2,2)*pow(d,5)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 96*c2*d2*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) - 
16*b2*c2*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) - 16*a2*d2*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) - 
432*a2*b2*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) - 1232*pow(a2,2)*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) + 
288*b2*c2*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 1260*a2*c2*d*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 288*a2*d2*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 
630*d*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 3312*a2*b2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 
3360*pow(a2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 288*b2*c2*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 
1260*a2*c2*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 288*a2*d2*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 
630*d*pow(b2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 3312*a2*b2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 
3360*pow(a2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 420*a2*c2*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,5) + 
2304*a2*b2*d*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,5) + 210*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,5) + 
3640*pow(a2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 420*a2*c2*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,5) + 
2304*a2*b2*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 210*pow(b2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 
3640*pow(a2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,5) + 576*a2*b2*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,5) + 
1890*d*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,5) - 576*a2*b2*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,5) + 
1890*d*pow(a2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,5) + 378*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,5) - 
378*pow(a2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,5) + 576*b2*c2*d*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 180*b2*d2*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 
576*a2*d*d2*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 576*b2*c2*k2*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 2520*a2*c2*d*k2*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 
576*a2*d2*k2*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 1260*d*k2*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 90*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 
1200*a2*c2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 6624*a2*b2*k2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 
600*pow(b2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 2016*a2*b2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 
6720*k2*pow(a2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 576*b2*c2*d*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 180*b2*d2*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) + 
576*a2*d*d2*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 576*b2*c2*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) + 2520*a2*c2*d*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 
576*a2*d2*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) + 1260*d*k2*pow(b2,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 90*pow(c2,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 
1200*a2*c2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 6624*a2*b2*k2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 
600*pow(b2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) + 2016*a2*b2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) + 
6720*k2*pow(a2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 1152*b2*c2*d*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 180*b2*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 
1152*a2*d*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 90*pow(c2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 4980*a2*c2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 
2490*pow(b2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 17856*a2*b2*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 
27090*pow(a2,2)*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 180*b2*d2*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,6) + 
90*pow(c2,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,6) - 60*a2*c2*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,6) - 
30*pow(b2,2)*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,6) - 1470*pow(a2,2)*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,6) + 
1260*a2*c2*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,6) + 6912*a2*b2*d*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,6) + 
630*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,6) + 10920*pow(a2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,6) - 
1260*a2*c2*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,6) + 6912*a2*b2*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,6) - 
630*pow(b2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,6) - 10920*pow(a2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,6) + 
2304*a2*b2*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,6) + 7560*d*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,6) - 
2304*a2*b2*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,6) + 7560*d*pow(a2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,6) + 
1890*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,6) - 1890*pow(a2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,6) + 576*b2*c2*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 
2520*a2*c2*d*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 576*a2*d2*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 2520*a2*c2*k2*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 
13824*a2*b2*d*k2*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 1260*d*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 1260*k2*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 
6624*a2*b2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 21840*k2*pow(a2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 
6720*pow(a2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 576*b2*c2*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) + 2520*a2*c2*d*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) - 
576*a2*d2*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) - 2520*a2*c2*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) + 13824*a2*b2*d*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) + 
1260*d*pow(b2,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) - 1260*k2*pow(b2,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) - 6624*a2*b2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) - 
21840*k2*pow(a2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) + 6720*pow(a2,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) - 
576*b2*c2*pow(M_E,k2*pow(RP,-1))*pow(RP,7) - 5040*a2*c2*d*pow(M_E,k2*pow(RP,-1))*pow(RP,7) - 576*a2*d2*pow(M_E,k2*pow(RP,-1))*pow(RP,7) - 
2520*d*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) - 27360*a2*b2*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) - 
58800*pow(a2,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 576*b2*c2*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,7) + 576*a2*d2*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,7) - 
288*a2*b2*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,7) + 6912*a2*b2*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,7) + 
22680*d*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,7) - 6912*a2*b2*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,7) + 
22680*d*pow(a2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,7) + 7560*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,7) - 
7560*pow(a2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,7) + 2520*a2*c2*pow(M_E,d*pow(RP,-1))*pow(RP,8) + 13824*a2*b2*d*pow(M_E,d*pow(RP,-1))*pow(RP,8) + 
13824*a2*b2*k2*pow(M_E,d*pow(RP,-1))*pow(RP,8) + 45360*d*k2*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,8) + 1260*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,8) + 
21840*pow(a2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,8) - 2520*a2*c2*pow(M_E,3*d*pow(RP,-1))*pow(RP,8) + 13824*a2*b2*d*pow(M_E,3*d*pow(RP,-1))*pow(RP,8) - 
13824*a2*b2*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,8) + 45360*d*k2*pow(a2,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,8) - 1260*pow(b2,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,8) - 
21840*pow(a2,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,8) - 2520*a2*c2*pow(M_E,k2*pow(RP,-1))*pow(RP,8) - 27648*a2*b2*d*pow(M_E,k2*pow(RP,-1))*pow(RP,8) - 
1260*pow(b2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,8) - 89880*pow(a2,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,8) + 2520*a2*c2*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,8) + 
1260*pow(b2,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,8) - 840*pow(a2,2)*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,8) + 
22680*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,8) - 22680*pow(a2,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,8) + 
13824*a2*b2*pow(M_E,d*pow(RP,-1))*pow(RP,9) + 45360*d*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,9) + 45360*k2*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,9) - 
13824*a2*b2*pow(M_E,3*d*pow(RP,-1))*pow(RP,9) + 45360*d*pow(a2,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,9) - 45360*k2*pow(a2,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,9) - 
13824*a2*b2*pow(M_E,k2*pow(RP,-1))*pow(RP,9) - 90720*d*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,9) + 13824*a2*b2*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,9) + 
45360*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,10) - 45360*pow(a2,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,10) - 45360*pow(a2,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,10) + 
45360*pow(a2,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,10)))/2.;

return res;
}

double SDH_Integral_z_44_case_2_sub_2( double AP, double RP, double k1, double k2, double a2, double b2, double c2, double d2, double d )
{
double res =(3*AP*RP*pow(d,-4)*pow(M_E,-((2*d + k1)*pow(RP,-1)))*(-8*c2*d2*k1*RP*pow(d,3)*pow(M_E,d*pow(RP,-1)) - 3*k1*RP*pow(d,2)*pow(d2,2)*pow(M_E,d*pow(RP,-1)) - 
k1*pow(d,3)*pow(d2,2)*pow(M_E,d*pow(RP,-1)) - 3*RP*pow(d,3)*pow(d2,2)*pow(M_E,d*pow(RP,-1)) + 14*c2*d2*RP*pow(d,4)*pow(M_E,k1*pow(RP,-1)) + 
2*c2*d2*pow(d,5)*pow(M_E,k1*pow(RP,-1)) + 16*b2*d2*RP*pow(d,5)*pow(M_E,k1*pow(RP,-1)) + 8*RP*pow(c2,2)*pow(d,5)*pow(M_E,k1*pow(RP,-1)) + 
2*b2*d2*pow(d,6)*pow(M_E,k1*pow(RP,-1)) + 18*b2*c2*RP*pow(d,6)*pow(M_E,k1*pow(RP,-1)) + 18*a2*d2*RP*pow(d,6)*pow(M_E,k1*pow(RP,-1)) + 
pow(c2,2)*pow(d,6)*pow(M_E,k1*pow(RP,-1)) + 2*b2*c2*pow(d,7)*pow(M_E,k1*pow(RP,-1)) + 2*a2*d2*pow(d,7)*pow(M_E,k1*pow(RP,-1)) + 
20*a2*c2*RP*pow(d,7)*pow(M_E,k1*pow(RP,-1)) + 10*RP*pow(b2,2)*pow(d,7)*pow(M_E,k1*pow(RP,-1)) + 2*a2*c2*pow(d,8)*pow(M_E,k1*pow(RP,-1)) + 
22*a2*b2*RP*pow(d,8)*pow(M_E,k1*pow(RP,-1)) + pow(b2,2)*pow(d,8)*pow(M_E,k1*pow(RP,-1)) + 2*a2*b2*pow(d,9)*pow(M_E,k1*pow(RP,-1)) + 
12*RP*pow(a2,2)*pow(d,9)*pow(M_E,k1*pow(RP,-1)) + pow(a2,2)*pow(d,10)*pow(M_E,k1*pow(RP,-1)) + 6*RP*pow(d,3)*pow(d2,2)*pow(M_E,k1*pow(RP,-1)) + 
pow(d,4)*pow(d2,2)*pow(M_E,k1*pow(RP,-1)) + 2*RP*gsl_sf_expint_Ei(-(k1*pow(RP,-1)))*pow(d,3)*pow(d2,2)*pow(M_E,(d + k1)*pow(RP,-1)) - 
2*RP*gsl_sf_expint_Ei(k1*pow(RP,-1))*pow(d,3)*pow(d2,2)*pow(M_E,(d + k1)*pow(RP,-1)) - 2*c2*d2*RP*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 
2*c2*d2*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 4*b2*d2*RP*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 2*RP*pow(c2,2)*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 
2*b2*d2*pow(d,6)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 6*b2*c2*RP*pow(d,6)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 6*a2*d2*RP*pow(d,6)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 
pow(c2,2)*pow(d,6)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 2*b2*c2*pow(d,7)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 2*a2*d2*pow(d,7)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 
8*a2*c2*RP*pow(d,7)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 4*RP*pow(b2,2)*pow(d,7)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 2*a2*c2*pow(d,8)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 
10*a2*b2*RP*pow(d,8)*pow(M_E,(2*d + k1)*pow(RP,-1)) + pow(b2,2)*pow(d,8)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 2*a2*b2*pow(d,9)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 
6*RP*pow(a2,2)*pow(d,9)*pow(M_E,(2*d + k1)*pow(RP,-1)) + pow(a2,2)*pow(d,10)*pow(M_E,(2*d + k1)*pow(RP,-1)) + pow(d,4)*pow(d2,2)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 
8*c2*d2*k1*RP*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1)) - 3*k1*RP*pow(d,2)*pow(d2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1)) - 
k1*pow(d,3)*pow(d2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1)) + 3*RP*pow(d,3)*pow(d2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1)) - 6*c2*d2*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 
2*c2*d2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 10*b2*d2*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 5*RP*pow(c2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 
6*c2*d2*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) - 2*c2*d2*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) + 
10*b2*d2*RP*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) + 5*RP*pow(c2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) - 
6*b2*d2*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 3*RP*pow(c2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 2*b2*d2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 
12*b2*c2*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 12*a2*d2*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - pow(c2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 
6*b2*d2*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - 3*RP*pow(c2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - 
2*b2*d2*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) + 12*b2*c2*RP*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) + 
12*a2*d2*RP*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - pow(c2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - 
6*b2*c2*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 6*a2*d2*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 2*b2*c2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 
2*a2*d2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 14*a2*c2*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 7*RP*pow(b2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 
6*b2*c2*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) - 6*a2*d2*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) - 
2*b2*c2*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) - 2*a2*d2*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) + 
14*a2*c2*RP*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) + 7*RP*pow(b2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) - 
6*a2*c2*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 3*RP*pow(b2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 2*a2*c2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 
16*a2*b2*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,5) - pow(b2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 
6*a2*c2*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) - 3*RP*pow(b2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) - 
2*a2*c2*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) + 16*a2*b2*RP*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) - 
pow(b2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) - 6*a2*b2*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,6) - 
2*a2*b2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,6) - 9*RP*pow(a2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,6) - 
6*a2*b2*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6) - 2*a2*b2*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6) + 
9*RP*pow(a2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6) - 3*RP*pow(a2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,7) - 
pow(a2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,7) - 3*RP*pow(a2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,7) - 
pow(a2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,7) - 28*c2*d2*k1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 
12*c2*d2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 24*b2*d2*k1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 12*k1*pow(c2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 
6*d*k1*pow(d2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 11*pow(d,2)*pow(d2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 52*c2*d2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
70*b2*d2*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 35*pow(c2,2)*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 92*b2*c2*pow(d,5)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
92*a2*d2*pow(d,5)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 118*a2*c2*pow(d,6)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 59*pow(b2,2)*pow(d,6)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
148*a2*b2*pow(d,7)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 91*pow(a2,2)*pow(d,8)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
17*pow(d,2)*pow(d2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 8*gsl_sf_expint_Ei(-(k1*pow(RP,-1)))*pow(d,2)*pow(d2,2)*pow(M_E,(d + k1)*pow(RP,-1))*pow(RP,2) - 
8*gsl_sf_expint_Ei(k1*pow(RP,-1))*pow(d,2)*pow(d2,2)*pow(M_E,(d + k1)*pow(RP,-1))*pow(RP,2) - 4*c2*d2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 
2*b2*d2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + pow(c2,2)*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 
12*b2*c2*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 12*a2*d2*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 
26*a2*c2*pow(d,6)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 13*pow(b2,2)*pow(d,6)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 
44*a2*b2*pow(d,7)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 33*pow(a2,2)*pow(d,8)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) - 
5*pow(d,2)*pow(d2,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 28*c2*d2*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) - 
12*c2*d2*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) - 24*b2*d2*k1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) - 
12*k1*pow(c2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) - 6*d*k1*pow(d2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) + 
11*pow(d,2)*pow(d2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) - 12*c2*d*d2*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
34*b2*d2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 17*pow(c2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
40*b2*c2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 40*a2*d2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
12*c2*d*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 34*b2*d2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
17*pow(c2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 40*b2*c2*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
40*a2*d2*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 12*b2*d*d2*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
6*d*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 40*b2*c2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
40*a2*d2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 60*a2*c2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
30*pow(b2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 12*b2*d*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
6*d*pow(c2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 40*b2*c2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
40*a2*d2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 60*a2*c2*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
30*pow(b2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 12*b2*c2*d*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
12*a2*d*d2*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 46*a2*c2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
23*pow(b2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 84*a2*b2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
12*b2*c2*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 12*a2*d*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 
46*a2*c2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 23*pow(b2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
84*a2*b2*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 12*a2*c2*d*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 
6*d*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 52*a2*b2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 
56*pow(a2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 12*a2*c2*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 
6*d*pow(b2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 52*a2*b2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 
56*pow(a2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 12*a2*b2*d*pow(M_E,d*pow(RP,-1))*pow(k1,6)*pow(RP,2) - 
29*pow(a2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,6)*pow(RP,2) - 12*a2*b2*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6)*pow(RP,2) + 
29*pow(a2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6)*pow(RP,2) - 6*d*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,7)*pow(RP,2) - 
6*d*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,7)*pow(RP,2) - 60*c2*d*d2*k1*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 
44*c2*d2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 84*b2*d2*k1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 42*k1*pow(c2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 
24*b2*d2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 80*b2*c2*k1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 80*a2*d2*k1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 
12*pow(c2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 24*d*pow(d2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 6*k1*pow(d2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 
116*c2*d2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 192*b2*d2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 96*pow(c2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 
312*b2*c2*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 312*a2*d2*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 488*a2*c2*pow(d,5)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 
244*pow(b2,2)*pow(d,5)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 732*a2*b2*pow(d,6)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 
528*pow(a2,2)*pow(d,7)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 30*d*pow(d2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 
18*d*gsl_sf_expint_Ei(-(k1*pow(RP,-1)))*pow(d2,2)*pow(M_E,(d + k1)*pow(RP,-1))*pow(RP,3) - 
18*d*gsl_sf_expint_Ei(k1*pow(RP,-1))*pow(d2,2)*pow(M_E,(d + k1)*pow(RP,-1))*pow(RP,3) - 4*c2*d2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 
16*b2*c2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 16*a2*d2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 
64*a2*c2*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 32*pow(b2,2)*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 
156*a2*b2*pow(d,6)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 152*pow(a2,2)*pow(d,7)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 
18*d*pow(d2,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) + 60*c2*d*d2*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 
44*c2*d2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 84*b2*d2*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 
42*k1*pow(c2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 24*b2*d2*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 
80*b2*c2*k1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 80*a2*d2*k1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 
12*pow(c2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 24*d*pow(d2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 
6*k1*pow(d2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 12*c2*d2*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 72*b2*d*d2*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
36*d*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 136*b2*c2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
136*a2*d2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 180*a2*c2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
90*pow(b2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 12*c2*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
72*b2*d*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 36*d*pow(c2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
136*b2*c2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 136*a2*d2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
180*a2*c2*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 90*pow(b2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
84*b2*c2*d*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 12*b2*d2*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 84*a2*d*d2*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
6*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 200*a2*c2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
100*pow(b2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 336*a2*b2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
84*b2*c2*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 12*b2*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
84*a2*d*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 6*pow(c2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
200*a2*c2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 100*pow(b2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
336*a2*b2*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 12*b2*c2*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 
96*a2*c2*d*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 12*a2*d2*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 
48*d*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 276*a2*b2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 
280*pow(a2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 12*b2*c2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 
96*a2*c2*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 12*a2*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 
48*d*pow(b2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 276*a2*b2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 
280*pow(a2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 12*a2*c2*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,3) - 
108*a2*b2*d*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,3) - 6*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,3) - 
182*pow(a2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,3) - 12*a2*c2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,3) + 
108*a2*b2*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,3) - 6*pow(b2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,3) - 
182*pow(a2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,3) - 12*a2*b2*pow(M_E,d*pow(RP,-1))*pow(k1,6)*pow(RP,3) - 
60*d*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,6)*pow(RP,3) - 12*a2*b2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6)*pow(RP,3) + 
60*d*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6)*pow(RP,3) - 6*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,7)*pow(RP,3) - 
6*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,7)*pow(RP,3) - 
2*RP*gsl_sf_expint_Ei(-(d*pow(RP,-1)))*pow(d2,2)*pow(M_E,(d + k1)*pow(RP,-1))*(4*RP*pow(d,2) + pow(d,3) + 9*d*pow(RP,2) + 9*pow(RP,3)) + 
2*RP*gsl_sf_expint_Ei(d*pow(RP,-1))*pow(d2,2)*pow(M_E,(d + k1)*pow(RP,-1))*(4*RP*pow(d,2) + pow(d,3) + 9*d*pow(RP,2) + 9*pow(RP,3)) - 
96*c2*d*d2*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 60*c2*d2*k1*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 180*b2*d*d2*k1*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
90*d*k1*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 84*b2*d2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 272*b2*c2*k1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
272*a2*d2*k1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 42*pow(c2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 80*b2*c2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
80*a2*d2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 360*a2*c2*k1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
180*k1*pow(b2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 24*pow(d2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 156*c2*d*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 
336*b2*d2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 168*pow(c2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 724*b2*c2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 
724*a2*d2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 1476*a2*c2*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 
738*pow(b2,2)*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 2796*a2*b2*pow(d,5)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 
2468*pow(a2,2)*pow(d,6)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 24*pow(d2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 
18*gsl_sf_expint_Ei(-(k1*pow(RP,-1)))*pow(d2,2)*pow(M_E,(d + k1)*pow(RP,-1))*pow(RP,4) - 
18*gsl_sf_expint_Ei(k1*pow(RP,-1))*pow(d2,2)*pow(M_E,(d + k1)*pow(RP,-1))*pow(RP,4) + 36*c2*d*d2*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) + 
24*b2*d2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) + 12*pow(c2,2)*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) + 
12*b2*c2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) + 12*a2*d2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) + 
84*a2*c2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) + 42*pow(b2,2)*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) + 
372*a2*b2*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) + 528*pow(a2,2)*pow(d,6)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) - 
24*pow(d2,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) - 96*c2*d*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 60*c2*d2*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 
180*b2*d*d2*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 90*d*k1*pow(c2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 
84*b2*d2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 272*b2*c2*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 
272*a2*d2*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 42*pow(c2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 
80*b2*c2*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 80*a2*d2*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 
360*a2*c2*k1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 180*k1*pow(b2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 
24*pow(d2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 288*b2*c2*d*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 72*b2*d2*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
288*a2*d*d2*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 36*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
600*a2*c2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 300*pow(b2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
1008*a2*b2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 288*b2*c2*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
72*b2*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 288*a2*d*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
36*pow(c2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 600*a2*c2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
300*pow(b2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 1008*a2*b2*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
84*b2*c2*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 420*a2*c2*d*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 84*a2*d2*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 
210*d*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 1104*a2*b2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 
1120*pow(a2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 84*b2*c2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 
420*a2*c2*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 84*a2*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 
210*d*pow(b2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 1104*a2*b2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 
1120*pow(a2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 96*a2*c2*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,4) - 
576*a2*b2*d*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,4) - 48*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,4) - 
910*pow(a2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,4) + 96*a2*c2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,4) - 
576*a2*b2*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,4) + 48*pow(b2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,4) + 
910*pow(a2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,4) - 108*a2*b2*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,4) - 
378*d*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,4) + 108*a2*b2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,4) - 
378*d*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,4) - 60*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,6)*pow(RP,4) + 
60*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6)*pow(RP,4) - 96*c2*d2*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 180*b2*d*d2*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 
576*b2*c2*d*k1*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 180*b2*d2*k1*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 576*a2*d*d2*k1*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 
90*d*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 90*k1*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 272*b2*c2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 
272*a2*d2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 1200*a2*c2*k1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 
600*k1*pow(b2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 360*a2*c2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 
2016*a2*b2*k1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 180*pow(b2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 96*c2*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 
360*b2*d*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 180*d*pow(c2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 1136*b2*c2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 
1136*a2*d2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 3240*a2*c2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 
1620*pow(b2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 8208*a2*b2*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 
9268*pow(a2,2)*pow(d,5)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 96*c2*d2*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) - 
16*b2*c2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) - 16*a2*d2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) - 
432*a2*b2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) - 1232*pow(a2,2)*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) - 
96*c2*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 180*b2*d*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 576*b2*c2*d*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 
180*b2*d2*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 576*a2*d*d2*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 
90*d*pow(c2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 90*k1*pow(c2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 
272*b2*c2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 272*a2*d2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 
1200*a2*c2*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 600*k1*pow(b2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 
360*a2*c2*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 2016*a2*b2*k1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 
180*pow(b2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 288*b2*c2*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 
1260*a2*c2*d*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 288*a2*d2*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 
630*d*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 3312*a2*b2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 
3360*pow(a2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 288*b2*c2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 
1260*a2*c2*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 288*a2*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 
630*d*pow(b2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 3312*a2*b2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 
3360*pow(a2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 420*a2*c2*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,5) - 
2304*a2*b2*d*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,5) - 210*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,5) - 
3640*pow(a2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,5) - 420*a2*c2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,5) + 
2304*a2*b2*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,5) - 210*pow(b2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,5) - 
3640*pow(a2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,5) - 576*a2*b2*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,5) - 
1890*d*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,5) - 576*a2*b2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,5) + 
1890*d*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,5) - 378*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,5) - 
378*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,5) - 576*b2*c2*d*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 180*b2*d2*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 
576*a2*d*d2*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 576*b2*c2*k1*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 2520*a2*c2*d*k1*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 
576*a2*d2*k1*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 1260*d*k1*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 90*pow(c2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 
1200*a2*c2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 6624*a2*b2*k1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 
600*pow(b2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 2016*a2*b2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 
6720*k1*pow(a2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 1152*b2*c2*d*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 180*b2*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 
1152*a2*d*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 90*pow(c2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 4980*a2*c2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 
2490*pow(b2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 17856*a2*b2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 
27090*pow(a2,2)*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 180*b2*d2*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,6) - 
90*pow(c2,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,6) + 60*a2*c2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,6) + 
30*pow(b2,2)*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,6) + 1470*pow(a2,2)*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,6) - 
576*b2*c2*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) + 180*b2*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 576*a2*d*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) + 
576*b2*c2*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 2520*a2*c2*d*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) + 
576*a2*d2*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 1260*d*k1*pow(b2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) + 
90*pow(c2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) + 1200*a2*c2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) + 
6624*a2*b2*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) + 600*pow(b2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 
2016*a2*b2*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 6720*k1*pow(a2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 
1260*a2*c2*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,6) - 6912*a2*b2*d*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,6) - 
630*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,6) - 10920*pow(a2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,6) + 
1260*a2*c2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,6) - 6912*a2*b2*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,6) + 
630*pow(b2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,6) + 10920*pow(a2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,6) - 
2304*a2*b2*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,6) - 7560*d*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,6) + 
2304*a2*b2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,6) - 7560*d*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,6) - 
1890*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,6) + 1890*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,6) - 
576*b2*c2*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 2520*a2*c2*d*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 576*a2*d2*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 
2520*a2*c2*k1*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 13824*a2*b2*d*k1*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 1260*d*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 
1260*k1*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 6624*a2*b2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 
21840*k1*pow(a2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 6720*pow(a2,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 576*b2*c2*pow(M_E,k1*pow(RP,-1))*pow(RP,7) + 
5040*a2*c2*d*pow(M_E,k1*pow(RP,-1))*pow(RP,7) + 576*a2*d2*pow(M_E,k1*pow(RP,-1))*pow(RP,7) + 2520*d*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) + 
27360*a2*b2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) + 58800*pow(a2,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) + 
576*b2*c2*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,7) + 576*a2*d2*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,7) - 
288*a2*b2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,7) - 576*b2*c2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) + 
2520*a2*c2*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) - 576*a2*d2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) - 2520*a2*c2*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) + 
13824*a2*b2*d*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) + 1260*d*pow(b2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) - 
1260*k1*pow(b2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) - 6624*a2*b2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) - 
21840*k1*pow(a2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) + 6720*pow(a2,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) - 
6912*a2*b2*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,7) - 22680*d*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,7) - 
6912*a2*b2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,7) + 22680*d*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,7) - 
7560*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,7) - 7560*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,7) - 
2520*a2*c2*pow(M_E,d*pow(RP,-1))*pow(RP,8) - 13824*a2*b2*d*pow(M_E,d*pow(RP,-1))*pow(RP,8) - 13824*a2*b2*k1*pow(M_E,d*pow(RP,-1))*pow(RP,8) - 
45360*d*k1*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,8) - 1260*pow(b2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,8) - 21840*pow(a2,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,8) + 
2520*a2*c2*pow(M_E,k1*pow(RP,-1))*pow(RP,8) + 27648*a2*b2*d*pow(M_E,k1*pow(RP,-1))*pow(RP,8) + 1260*pow(b2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,8) + 
89880*pow(a2,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,8) - 2520*a2*c2*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,8) - 
1260*pow(b2,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,8) + 840*pow(a2,2)*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,8) + 
2520*a2*c2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,8) - 13824*a2*b2*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,8) + 
13824*a2*b2*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,8) - 45360*d*k1*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,8) + 
1260*pow(b2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,8) + 21840*pow(a2,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,8) - 
22680*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,8) + 22680*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,8) - 
13824*a2*b2*pow(M_E,d*pow(RP,-1))*pow(RP,9) - 45360*d*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,9) - 45360*k1*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,9) + 
13824*a2*b2*pow(M_E,k1*pow(RP,-1))*pow(RP,9) + 90720*d*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,9) + 13824*a2*b2*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,9) - 
13824*a2*b2*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,9) + 45360*d*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,9) - 
45360*k1*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,9) - 45360*pow(a2,2)*pow(M_E,d*pow(RP,-1))*pow(RP,10) + 45360*pow(a2,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,10) - 
45360*pow(a2,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,10) + 45360*pow(a2,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,10)))/2.;

return res;
}

double SDH_Integral_z_44_case_3( double AP, double RP, double k1, double k2, double a2, double b2, double c2, double d2, double d )
{
double res =(-3*AP*pow(d,-4)*pow(M_E,-(d*pow(RP,-1)))*(RP*pow(d,3)*pow(d2,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*
pow(M_E,-((k1 + k2)*pow(RP,-1))) + 2*(-gsl_sf_expint_Ei(-(k1*pow(RP,-1))) + gsl_sf_expint_Ei(-(k2*pow(RP,-1))))*pow(d,3)*pow(d2,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*
pow(RP,2) + 2*pow(d,3)*pow(d2,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,-(k1*pow(RP,-1))) - pow(M_E,-(k2*pow(RP,-1))))*pow(RP,2) + 
3*pow(d,2)*pow(d2,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2) + 
4*c2*d2*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2) + 
2*c2*d2*RP*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
4*b2*d2*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
2*pow(c2,2)*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
6*c2*d2*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
8*(gsl_sf_expint_Ei(-(k1*pow(RP,-1))) - gsl_sf_expint_Ei(-(k2*pow(RP,-1))))*pow(d,2)*pow(d2,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,3) + 
4*c2*d2*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,-(k1*pow(RP,-1))) - pow(M_E,-(k2*pow(RP,-1))))*pow(RP,3) + 
8*pow(d,2)*pow(d2,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,k1*pow(RP,-1)) - pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3) + 
16*c2*d2*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3) + 
4*b2*d2*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3) + 
2*pow(c2,2)*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
pow(RP,3) + 6*d*pow(d2,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
pow(RP,3) + 12*c2*d*d2*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,3) + 
4*b2*c2*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,3) + 
4*a2*d2*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,3) + 
16*b2*d2*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,3) + 
8*pow(c2,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,3) + 
2*b2*d2*RP*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
RP*pow(c2,2)*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
4*b2*c2*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
4*a2*d2*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
12*b2*d*d2*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
6*d*pow(c2,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
4*a2*c2*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
2*pow(b2,2)*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
6*b2*d2*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
3*pow(c2,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
16*b2*c2*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
16*a2*d2*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
18*d*(-gsl_sf_expint_Ei(-(k1*pow(RP,-1))) + gsl_sf_expint_Ei(-(k2*pow(RP,-1))))*pow(d2,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,4) + 
18*d*pow(d2,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,-(k1*pow(RP,-1))) - pow(M_E,-(k2*pow(RP,-1))))*pow(RP,4) + 
16*c2*d2*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,k1*pow(RP,-1)) - pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4) + 
16*b2*d2*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4) + 
8*pow(c2,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4) + 
6*pow(d2,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4) + 
36*c2*d*d2*(1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4) + 
36*b2*d*d2*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,4) + 
18*d*pow(c2,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,4) + 
12*c2*d2*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,4) + 
16*b2*c2*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,4) + 
16*a2*d2*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,4) + 
36*b2*c2*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,4) + 
36*a2*d*d2*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,4) + 
12*b2*d2*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,4) + 
6*pow(c2,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,4) + 
16*a2*c2*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,4) + 
8*pow(b2,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,4) + 
2*b2*c2*RP*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
2*a2*d2*RP*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
4*a2*c2*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
2*pow(b2,2)*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
12*b2*c2*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
12*a2*d*d2*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
4*a2*b2*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
36*a2*c2*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*
(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
18*d*pow(b2,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*
(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
6*b2*c2*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
6*a2*d2*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
16*a2*c2*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
8*pow(b2,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
12*b2*c2*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
12*a2*d2*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
16*a2*b2*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
18*(gsl_sf_expint_Ei(-(k1*pow(RP,-1))) - gsl_sf_expint_Ei(-(k2*pow(RP,-1))))*pow(d2,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,5) + 
36*c2*d*d2*(1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,-(k1*pow(RP,-1))) - pow(M_E,-(k2*pow(RP,-1))))*pow(RP,5) + 
18*pow(d2,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,k1*pow(RP,-1)) - pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,5) + 
36*c2*d2*(-1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,5) + 
36*b2*d*d2*(1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,5) + 
18*d*pow(c2,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,5) + 
36*b2*c2*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,5) + 
36*a2*d*d2*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,5) + 
36*b2*d2*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,5) + 
18*pow(c2,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,5) + 
36*a2*c2*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,5) + 
18*d*pow(b2,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,5) + 
36*b2*c2*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,5) + 
36*a2*d2*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,5) + 
36*a2*b2*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4)))*pow(RP,5) + 
36*a2*c2*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4)))*pow(RP,5) + 
18*pow(b2,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4)))*pow(RP,5) + 
2*a2*c2*RP*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5)) - 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
RP*pow(b2,2)*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5)) - 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
4*a2*b2*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5)) - 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
12*a2*c2*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5)) - 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
6*d*pow(b2,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5)) - 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
2*pow(a2,2)*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5)) - 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
36*a2*b2*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*
(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5)) - 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
18*d*pow(a2,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,5)*
(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5)) - 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
6*a2*c2*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))) + 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
3*pow(b2,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))) + 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
16*a2*b2*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))) + 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
12*a2*c2*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*
(-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))) + 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
6*pow(b2,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*
(-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))) + 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
8*pow(a2,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*
(-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))) + 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
36*a2*b2*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,5)*
(-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))) + 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
36*c2*d2*(-1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,k1*pow(RP,-1)) - pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,6) + 
36*b2*d2*(-1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,6) + 
18*pow(c2,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,6) + 
36*b2*c2*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,6) + 
36*a2*d2*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,6) + 
36*a2*c2*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,6) + 
18*pow(b2,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,6) + 
36*a2*b2*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4)))*pow(RP,6) + 
18*pow(a2,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))) + 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5)))*pow(RP,6) + 
2*a2*b2*RP*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 
720*pow(RP,6)) - pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 
720*k2*pow(RP,5) + 720*pow(RP,6))) + 2*pow(a2,2)*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 
720*pow(RP,6)) - pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 
720*k2*pow(RP,5) + 720*pow(RP,6))) + 12*a2*b2*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 
720*pow(RP,6)) - pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 
720*k2*pow(RP,5) + 720*pow(RP,6))) + 18*d*pow(a2,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*
(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 
720*pow(RP,6)) - pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 
720*k2*pow(RP,5) + 720*pow(RP,6))) + 6*a2*b2*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 
720*pow(RP,6))) + pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 
720*k2*pow(RP,5) + 720*pow(RP,6))) + 8*pow(a2,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(-(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 
720*pow(RP,6))) + pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 
720*k2*pow(RP,5) + 720*pow(RP,6))) + 12*a2*b2*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*
(-(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 
720*pow(RP,6))) + pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 
720*k2*pow(RP,5) + 720*pow(RP,6))) + 18*pow(a2,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,5)*
(-(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 
720*pow(RP,6))) + pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 
720*k2*pow(RP,5) + 720*pow(RP,6))) + RP*pow(a2,2)*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(7*RP*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RP,2) + 210*pow(k1,4)*pow(RP,3) + 840*pow(k1,3)*pow(RP,4) + 2520*pow(k1,2)*pow(RP,5) + 
5040*k1*pow(RP,6) + 5040*pow(RP,7)) - pow(M_E,k1*pow(RP,-1))*
(7*RP*pow(k2,6) + pow(k2,7) + 42*pow(k2,5)*pow(RP,2) + 210*pow(k2,4)*pow(RP,3) + 840*pow(k2,3)*pow(RP,4) + 2520*pow(k2,2)*pow(RP,5) + 5040*k2*pow(RP,6) + 
5040*pow(RP,7))) + 6*d*pow(a2,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(pow(M_E,k2*pow(RP,-1))*(7*RP*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RP,2) + 210*pow(k1,4)*pow(RP,3) + 840*pow(k1,3)*pow(RP,4) + 2520*pow(k1,2)*pow(RP,5) + 
5040*k1*pow(RP,6) + 5040*pow(RP,7)) - pow(M_E,k1*pow(RP,-1))*
(7*RP*pow(k2,6) + pow(k2,7) + 42*pow(k2,5)*pow(RP,2) + 210*pow(k2,4)*pow(RP,3) + 840*pow(k2,3)*pow(RP,4) + 2520*pow(k2,2)*pow(RP,5) + 5040*k2*pow(RP,6) + 
5040*pow(RP,7))) + 3*pow(a2,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(7*RP*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RP,2) + 210*pow(k1,4)*pow(RP,3) + 840*pow(k1,3)*pow(RP,4) + 2520*pow(k1,2)*pow(RP,5) + 
5040*k1*pow(RP,6) + 5040*pow(RP,7))) + pow(M_E,k1*pow(RP,-1))*
(7*RP*pow(k2,6) + pow(k2,7) + 42*pow(k2,5)*pow(RP,2) + 210*pow(k2,4)*pow(RP,3) + 840*pow(k2,3)*pow(RP,4) + 2520*pow(k2,2)*pow(RP,5) + 5040*k2*pow(RP,6) + 
5040*pow(RP,7))) + 6*pow(a2,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*
(-(pow(M_E,k2*pow(RP,-1))*(7*RP*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RP,2) + 210*pow(k1,4)*pow(RP,3) + 840*pow(k1,3)*pow(RP,4) + 2520*pow(k1,2)*pow(RP,5) + 
5040*k1*pow(RP,6) + 5040*pow(RP,7))) + pow(M_E,k1*pow(RP,-1))*
(7*RP*pow(k2,6) + pow(k2,7) + 42*pow(k2,5)*pow(RP,2) + 210*pow(k2,4)*pow(RP,3) + 840*pow(k2,3)*pow(RP,4) + 2520*pow(k2,2)*pow(RP,5) + 5040*k2*pow(RP,6) + 
5040*pow(RP,7)))))/2.;
return res;
}
