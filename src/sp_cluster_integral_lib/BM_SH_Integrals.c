#include<stdio.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_sf_expint.h>

/* !!! CALL M_Ei function: 
*   *
*     * gsl_sf_expint_Ei( arg );
*       *
*         */

/* !!! Integral Arguments
*   *
*     * k1: a knot for the left  boundary for spline cubic polynomial
*       * k2: a knot for the right boundary for spline cubic polynomial
*         *
*           * d : the distance between sp-electrons and an external point charge
*             *
*               * a1, b1, c1, d1 : 1st cubic polynomial coefficients
*                 * a2, b2, c2, d2 : 2nd cubic polynomial coefficients
*                   *
*                     */



/* ******************************************************************************************     1. < S | Kernel_SH | S >        */

/* BM_SS SS */

// 1
double BMSH_Integral_11_case_1( double AS, double RS, double k1, double k2, double a1, double b1, double c1, double d1, double d )
{
double res = AS*pow(2*d,-1)*pow(M_E,-((d + k1 + k2)*pow(RS,-1)))*pow(RS,2)*
(-5*RS*pow(c1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,3) - 5*RS*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,3) - 
12*b1*c1*RS*pow(M_E,k2*pow(RS,-1))*pow(k1,4) - pow(c1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,4) - 
12*b1*c1*RS*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,4) + pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,4) - 
2*b1*c1*pow(M_E,k2*pow(RS,-1))*pow(k1,5) - 14*a1*c1*RS*pow(M_E,k2*pow(RS,-1))*pow(k1,5) - 
7*RS*pow(b1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,5) + 2*b1*c1*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,5) - 
14*a1*c1*RS*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,5) - 7*RS*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,5) - 
2*a1*c1*pow(M_E,k2*pow(RS,-1))*pow(k1,6) - 16*a1*b1*RS*pow(M_E,k2*pow(RS,-1))*pow(k1,6) - 
pow(b1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,6) + 2*a1*c1*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,6) - 
16*a1*b1*RS*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,6) + pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,6) - 
2*a1*b1*pow(M_E,k2*pow(RS,-1))*pow(k1,7) - 9*RS*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,7) + 
2*a1*b1*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,7) - 9*RS*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,7) - 
pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,8) + pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,8) + 
5*RS*pow(c1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,3) + 5*RS*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,3) + 
12*b1*c1*RS*pow(M_E,k1*pow(RS,-1))*pow(k2,4) + pow(c1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,4) + 
12*b1*c1*RS*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,4) - pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,4) + 
2*b1*c1*pow(M_E,k1*pow(RS,-1))*pow(k2,5) + 14*a1*c1*RS*pow(M_E,k1*pow(RS,-1))*pow(k2,5) + 
7*RS*pow(b1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,5) - 2*b1*c1*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,5) + 
14*a1*c1*RS*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,5) + 7*RS*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,5) + 
2*a1*c1*pow(M_E,k1*pow(RS,-1))*pow(k2,6) + 16*a1*b1*RS*pow(M_E,k1*pow(RS,-1))*pow(k2,6) + 
pow(b1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,6) - 2*a1*c1*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,6) + 
16*a1*b1*RS*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,6) - pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,6) + 
2*a1*b1*pow(M_E,k1*pow(RS,-1))*pow(k2,7) + 9*RS*pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,7) - 
2*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,7) + 9*RS*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,7) + 
pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,8) - pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,8) - 
15*pow(c1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,2)*pow(RS,2) + 
15*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,2)*pow(RS,2) - 48*b1*c1*pow(M_E,k2*pow(RS,-1))*pow(k1,3)*pow(RS,2) + 
48*b1*c1*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,3)*pow(RS,2) - 70*a1*c1*pow(M_E,k2*pow(RS,-1))*pow(k1,4)*pow(RS,2) - 
35*pow(b1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,4)*pow(RS,2) + 70*a1*c1*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,4)*pow(RS,2) + 
35*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,4)*pow(RS,2) - 96*a1*b1*pow(M_E,k2*pow(RS,-1))*pow(k1,5)*pow(RS,2) + 
96*a1*b1*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,5)*pow(RS,2) - 63*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,6)*pow(RS,2) + 
63*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,6)*pow(RS,2) + 
15*pow(c1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,2)*pow(RS,2) - 
15*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,2)*pow(RS,2) + 48*b1*c1*pow(M_E,k1*pow(RS,-1))*pow(k2,3)*pow(RS,2) - 
48*b1*c1*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,3)*pow(RS,2) + 70*a1*c1*pow(M_E,k1*pow(RS,-1))*pow(k2,4)*pow(RS,2) + 
35*pow(b1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,4)*pow(RS,2) - 70*a1*c1*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,4)*pow(RS,2) - 
35*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,4)*pow(RS,2) + 96*a1*b1*pow(M_E,k1*pow(RS,-1))*pow(k2,5)*pow(RS,2) - 
96*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,5)*pow(RS,2) + 63*pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,6)*pow(RS,2) - 
63*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,6)*pow(RS,2) + 
pow(d1,2)*(pow(M_E,(2*k1 + k2)*pow(RS,-1))*(-3*k1*RS + pow(k1,2) + 3*pow(RS,2)) - 
pow(M_E,k2*pow(RS,-1))*(3*k1*RS + pow(k1,2) + 3*pow(RS,2)) - 
pow(M_E,(k1 + 2*k2)*pow(RS,-1))*(-3*k2*RS + pow(k2,2) + 3*pow(RS,2)) + 
pow(M_E,k1*pow(RS,-1))*(3*k2*RS + pow(k2,2) + 3*pow(RS,2))) + 30*k2*pow(c1,2)*pow(M_E,k1*pow(RS,-1))*pow(RS,3) - 
30*k1*pow(c1,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,3) - 30*k1*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(RS,3) + 
30*k2*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(RS,3) - 144*b1*c1*pow(M_E,k2*pow(RS,-1))*pow(k1,2)*pow(RS,3) - 
144*b1*c1*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,2)*pow(RS,3) - 280*a1*c1*pow(M_E,k2*pow(RS,-1))*pow(k1,3)*pow(RS,3) - 
140*pow(b1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,3)*pow(RS,3) - 280*a1*c1*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,3)*pow(RS,3) - 
140*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,3)*pow(RS,3) - 480*a1*b1*pow(M_E,k2*pow(RS,-1))*pow(k1,4)*pow(RS,3) - 
480*a1*b1*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,4)*pow(RS,3) - 378*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,5)*pow(RS,3) - 
378*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,5)*pow(RS,3) + 144*b1*c1*pow(M_E,k1*pow(RS,-1))*pow(k2,2)*pow(RS,3) + 
144*b1*c1*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,2)*pow(RS,3) + 280*a1*c1*pow(M_E,k1*pow(RS,-1))*pow(k2,3)*pow(RS,3) + 
140*pow(b1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,3)*pow(RS,3) + 280*a1*c1*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,3)*pow(RS,3) + 
140*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,3)*pow(RS,3) + 480*a1*b1*pow(M_E,k1*pow(RS,-1))*pow(k2,4)*pow(RS,3) + 
480*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,4)*pow(RS,3) + 378*pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,5)*pow(RS,3) + 
378*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,5)*pow(RS,3) + 288*b1*c1*k2*pow(M_E,k1*pow(RS,-1))*pow(RS,4) + 
30*pow(c1,2)*pow(M_E,k1*pow(RS,-1))*pow(RS,4) - 288*b1*c1*k1*pow(M_E,k2*pow(RS,-1))*pow(RS,4) - 
30*pow(c1,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,4) + 288*b1*c1*k1*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(RS,4) + 
30*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(RS,4) - 288*b1*c1*k2*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(RS,4) - 
30*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(RS,4) - 840*a1*c1*pow(M_E,k2*pow(RS,-1))*pow(k1,2)*pow(RS,4) - 
420*pow(b1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,2)*pow(RS,4) + 840*a1*c1*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,2)*pow(RS,4) + 
420*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,2)*pow(RS,4) - 
1920*a1*b1*pow(M_E,k2*pow(RS,-1))*pow(k1,3)*pow(RS,4) + 1920*a1*b1*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,3)*pow(RS,4) - 
1890*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,4)*pow(RS,4) + 
1890*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,4)*pow(RS,4) + 
840*a1*c1*pow(M_E,k1*pow(RS,-1))*pow(k2,2)*pow(RS,4) + 420*pow(b1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,2)*pow(RS,4) - 
840*a1*c1*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,2)*pow(RS,4) - 
420*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,2)*pow(RS,4) + 
1920*a1*b1*pow(M_E,k1*pow(RS,-1))*pow(k2,3)*pow(RS,4) - 1920*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,3)*pow(RS,4) + 
1890*pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,4)*pow(RS,4) - 
1890*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,4)*pow(RS,4) + 288*b1*c1*pow(M_E,k1*pow(RS,-1))*pow(RS,5) + 
1680*a1*c1*k2*pow(M_E,k1*pow(RS,-1))*pow(RS,5) + 840*k2*pow(b1,2)*pow(M_E,k1*pow(RS,-1))*pow(RS,5) - 
288*b1*c1*pow(M_E,k2*pow(RS,-1))*pow(RS,5) - 1680*a1*c1*k1*pow(M_E,k2*pow(RS,-1))*pow(RS,5) - 
840*k1*pow(b1,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,5) - 288*b1*c1*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(RS,5) - 
1680*a1*c1*k1*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(RS,5) - 840*k1*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(RS,5) + 
288*b1*c1*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(RS,5) + 1680*a1*c1*k2*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(RS,5) + 
840*k2*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(RS,5) - 5760*a1*b1*pow(M_E,k2*pow(RS,-1))*pow(k1,2)*pow(RS,5) - 
5760*a1*b1*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,2)*pow(RS,5) - 
7560*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,3)*pow(RS,5) - 
7560*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,3)*pow(RS,5) + 
5760*a1*b1*pow(M_E,k1*pow(RS,-1))*pow(k2,2)*pow(RS,5) + 5760*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,2)*pow(RS,5) + 
7560*pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,3)*pow(RS,5) + 
7560*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,3)*pow(RS,5) + 
2*d1*(c1*(pow(M_E,(2*k1 + k2)*pow(RS,-1))*(-4*RS*pow(k1,2) + pow(k1,3) + 8*k1*pow(RS,2) - 8*pow(RS,3)) - 
pow(M_E,k2*pow(RS,-1))*(4*RS*pow(k1,2) + pow(k1,3) + 8*k1*pow(RS,2) + 8*pow(RS,3)) + 
pow(M_E,(k1 + 2*k2)*pow(RS,-1))*(4*RS*pow(k2,2) - pow(k2,3) - 8*k2*pow(RS,2) + 8*pow(RS,3)) + 
pow(M_E,k1*pow(RS,-1))*(4*RS*pow(k2,2) + pow(k2,3) + 8*k2*pow(RS,2) + 8*pow(RS,3))) + 
b1*(pow(M_E,(2*k1 + k2)*pow(RS,-1))*(-5*RS*pow(k1,3) + pow(k1,4) + 15*pow(k1,2)*pow(RS,2) - 30*k1*pow(RS,3) + 
30*pow(RS,4)) - pow(M_E,k2*pow(RS,-1))*
(5*RS*pow(k1,3) + pow(k1,4) + 15*pow(k1,2)*pow(RS,2) + 30*k1*pow(RS,3) + 30*pow(RS,4)) - 
pow(M_E,(k1 + 2*k2)*pow(RS,-1))*(-5*RS*pow(k2,3) + pow(k2,4) + 15*pow(k2,2)*pow(RS,2) - 30*k2*pow(RS,3) + 
30*pow(RS,4)) + pow(M_E,k1*pow(RS,-1))*
(5*RS*pow(k2,3) + pow(k2,4) + 15*pow(k2,2)*pow(RS,2) + 30*k2*pow(RS,3) + 30*pow(RS,4))) + 
a1*(pow(M_E,(2*k1 + k2)*pow(RS,-1))*(-6*RS*pow(k1,4) + pow(k1,5) + 24*pow(k1,3)*pow(RS,2) - 72*pow(k1,2)*pow(RS,3) + 
144*k1*pow(RS,4) - 144*pow(RS,5)) - 
pow(M_E,k2*pow(RS,-1))*(6*RS*pow(k1,4) + pow(k1,5) + 24*pow(k1,3)*pow(RS,2) + 72*pow(k1,2)*pow(RS,3) + 
144*k1*pow(RS,4) + 144*pow(RS,5)) + 
pow(M_E,(k1 + 2*k2)*pow(RS,-1))*(6*RS*pow(k2,4) - pow(k2,5) - 24*pow(k2,3)*pow(RS,2) + 72*pow(k2,2)*pow(RS,3) - 
144*k2*pow(RS,4) + 144*pow(RS,5)) + 
pow(M_E,k1*pow(RS,-1))*(6*RS*pow(k2,4) + pow(k2,5) + 24*pow(k2,3)*pow(RS,2) + 72*pow(k2,2)*pow(RS,3) + 
144*k2*pow(RS,4) + 144*pow(RS,5)))) + 1680*a1*c1*pow(M_E,k1*pow(RS,-1))*pow(RS,6) + 
11520*a1*b1*k2*pow(M_E,k1*pow(RS,-1))*pow(RS,6) + 840*pow(b1,2)*pow(M_E,k1*pow(RS,-1))*pow(RS,6) - 
1680*a1*c1*pow(M_E,k2*pow(RS,-1))*pow(RS,6) - 11520*a1*b1*k1*pow(M_E,k2*pow(RS,-1))*pow(RS,6) - 
840*pow(b1,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,6) + 1680*a1*c1*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(RS,6) + 
11520*a1*b1*k1*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(RS,6) + 840*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(RS,6) - 
1680*a1*c1*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(RS,6) - 11520*a1*b1*k2*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(RS,6) - 
840*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(RS,6) - 22680*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,2)*pow(RS,6) + 
22680*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,2)*pow(RS,6) + 
22680*pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,2)*pow(RS,6) - 
22680*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,2)*pow(RS,6) + 11520*a1*b1*pow(M_E,k1*pow(RS,-1))*pow(RS,7) + 
45360*k2*pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(RS,7) - 11520*a1*b1*pow(M_E,k2*pow(RS,-1))*pow(RS,7) - 
45360*k1*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,7) - 11520*a1*b1*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(RS,7) - 
45360*k1*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(RS,7) + 11520*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(RS,7) + 
45360*k2*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(RS,7) + 
d*(pow(d1,2)*((k2 + RS)*pow(M_E,k1*pow(RS,-1)) - (k1 + RS)*pow(M_E,k2*pow(RS,-1)) + 
(-k1 + RS)*pow(M_E,(2*k1 + k2)*pow(RS,-1)) + (k2 - RS)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))) - 
5*RS*pow(b1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,4) + 5*RS*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,4) - 
12*a1*b1*RS*pow(M_E,k2*pow(RS,-1))*pow(k1,5) - pow(b1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,5) + 
12*a1*b1*RS*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,5) - pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,5) - 
2*a1*b1*pow(M_E,k2*pow(RS,-1))*pow(k1,6) - 7*RS*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,6) - 
2*a1*b1*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,6) + 7*RS*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,6) - 
pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,7) - pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,7) + 
5*RS*pow(b1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,4) - 5*RS*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,4) + 
12*a1*b1*RS*pow(M_E,k1*pow(RS,-1))*pow(k2,5) + pow(b1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,5) - 
12*a1*b1*RS*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,5) + pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,5) + 
2*a1*b1*pow(M_E,k1*pow(RS,-1))*pow(k2,6) + 7*RS*pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,6) + 
2*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,6) - 7*RS*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,6) + 
pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,7) + pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,7) - 
20*pow(b1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,3)*pow(RS,2) - 
20*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,3)*pow(RS,2) - 
60*a1*b1*pow(M_E,k2*pow(RS,-1))*pow(k1,4)*pow(RS,2) - 60*a1*b1*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,4)*pow(RS,2) - 
42*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,5)*pow(RS,2) - 
42*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,5)*pow(RS,2) + 
20*pow(b1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,3)*pow(RS,2) + 
20*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,3)*pow(RS,2) + 
60*a1*b1*pow(M_E,k1*pow(RS,-1))*pow(k2,4)*pow(RS,2) + 60*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,4)*pow(RS,2) + 
42*pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,5)*pow(RS,2) + 
42*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,5)*pow(RS,2) - 
60*pow(b1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,2)*pow(RS,3) + 
60*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,2)*pow(RS,3) - 
240*a1*b1*pow(M_E,k2*pow(RS,-1))*pow(k1,3)*pow(RS,3) + 240*a1*b1*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,3)*pow(RS,3) - 
210*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,4)*pow(RS,3) + 
210*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,4)*pow(RS,3) + 
60*pow(b1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,2)*pow(RS,3) - 
60*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,2)*pow(RS,3) + 
240*a1*b1*pow(M_E,k1*pow(RS,-1))*pow(k2,3)*pow(RS,3) - 240*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,3)*pow(RS,3) + 
210*pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,4)*pow(RS,3) - 
210*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,4)*pow(RS,3) + 
pow(c1,2)*(-(pow(M_E,(2*k1 + k2)*pow(RS,-1))*(-3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) - 6*pow(RS,3))) + 
pow(M_E,(k1 + 2*k2)*pow(RS,-1))*(-3*RS*pow(k2,2) + pow(k2,3) + 6*k2*pow(RS,2) - 6*pow(RS,3)) - 
pow(M_E,k2*pow(RS,-1))*(3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) + 6*pow(RS,3)) + 
pow(M_E,k1*pow(RS,-1))*(3*RS*pow(k2,2) + pow(k2,3) + 6*k2*pow(RS,2) + 6*pow(RS,3))) + 
120*k2*pow(b1,2)*pow(M_E,k1*pow(RS,-1))*pow(RS,4) - 120*k1*pow(b1,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,4) - 
120*k1*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(RS,4) + 
120*k2*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(RS,4) - 720*a1*b1*pow(M_E,k2*pow(RS,-1))*pow(k1,2)*pow(RS,4) - 
720*a1*b1*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,2)*pow(RS,4) - 
840*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,3)*pow(RS,4) - 
840*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,3)*pow(RS,4) + 
720*a1*b1*pow(M_E,k1*pow(RS,-1))*pow(k2,2)*pow(RS,4) + 720*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,2)*pow(RS,4) + 
840*pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,3)*pow(RS,4) + 
840*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,3)*pow(RS,4) + 
2*d1*(c1*(-(pow(M_E,(2*k1 + k2)*pow(RS,-1))*(-2*k1*RS + pow(k1,2) + 2*pow(RS,2))) - 
pow(M_E,k2*pow(RS,-1))*(2*k1*RS + pow(k1,2) + 2*pow(RS,2)) + 
pow(M_E,(k1 + 2*k2)*pow(RS,-1))*(-2*k2*RS + pow(k2,2) + 2*pow(RS,2)) + 
pow(M_E,k1*pow(RS,-1))*(2*k2*RS + pow(k2,2) + 2*pow(RS,2))) + 
b1*(-(pow(M_E,(2*k1 + k2)*pow(RS,-1))*(-3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) - 6*pow(RS,3))) + 
pow(M_E,(k1 + 2*k2)*pow(RS,-1))*(-3*RS*pow(k2,2) + pow(k2,3) + 6*k2*pow(RS,2) - 6*pow(RS,3)) - 
pow(M_E,k2*pow(RS,-1))*(3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) + 6*pow(RS,3)) + 
pow(M_E,k1*pow(RS,-1))*(3*RS*pow(k2,2) + pow(k2,3) + 6*k2*pow(RS,2) + 6*pow(RS,3))) + 
a1*(-(pow(M_E,(2*k1 + k2)*pow(RS,-1))*(-4*RS*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RS,2) - 24*k1*pow(RS,3) + 
24*pow(RS,4))) - pow(M_E,k2*pow(RS,-1))*
(4*RS*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RS,2) + 24*k1*pow(RS,3) + 24*pow(RS,4)) + 
pow(M_E,(k1 + 2*k2)*pow(RS,-1))*(-4*RS*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RS,2) - 24*k2*pow(RS,3) + 
24*pow(RS,4)) + pow(M_E,k1*pow(RS,-1))*
(4*RS*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RS,2) + 24*k2*pow(RS,3) + 24*pow(RS,4)))) + 
1440*a1*b1*k2*pow(M_E,k1*pow(RS,-1))*pow(RS,5) + 120*pow(b1,2)*pow(M_E,k1*pow(RS,-1))*pow(RS,5) - 
1440*a1*b1*k1*pow(M_E,k2*pow(RS,-1))*pow(RS,5) - 120*pow(b1,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,5) + 
1440*a1*b1*k1*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(RS,5) + 120*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(RS,5) - 
1440*a1*b1*k2*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(RS,5) - 120*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(RS,5) - 
2520*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,2)*pow(RS,5) + 
2520*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(k1,2)*pow(RS,5) + 
2520*pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,2)*pow(RS,5) - 
2520*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(k2,2)*pow(RS,5) + 
2*c1*(b1*(-(pow(M_E,(2*k1 + k2)*pow(RS,-1))*
(-4*RS*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RS,2) - 24*k1*pow(RS,3) + 24*pow(RS,4))) - 
pow(M_E,k2*pow(RS,-1))*(4*RS*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RS,2) + 24*k1*pow(RS,3) + 24*pow(RS,4)) + 
pow(M_E,(k1 + 2*k2)*pow(RS,-1))*(-4*RS*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RS,2) - 24*k2*pow(RS,3) + 
24*pow(RS,4)) + pow(M_E,k1*pow(RS,-1))*
(4*RS*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RS,2) + 24*k2*pow(RS,3) + 24*pow(RS,4))) + 
a1*(-(pow(M_E,(2*k1 + k2)*pow(RS,-1))*(-5*RS*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RS,2) - 
60*pow(k1,2)*pow(RS,3) + 120*k1*pow(RS,4) - 120*pow(RS,5))) + 
pow(M_E,(k1 + 2*k2)*pow(RS,-1))*(-5*RS*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RS,2) - 
60*pow(k2,2)*pow(RS,3) + 120*k2*pow(RS,4) - 120*pow(RS,5)) - 
pow(M_E,k2*pow(RS,-1))*(5*RS*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RS,2) + 60*pow(k1,2)*pow(RS,3) + 
120*k1*pow(RS,4) + 120*pow(RS,5)) + 
pow(M_E,k1*pow(RS,-1))*(5*RS*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RS,2) + 60*pow(k2,2)*pow(RS,3) + 
120*k2*pow(RS,4) + 120*pow(RS,5)))) + 1440*a1*b1*pow(M_E,k1*pow(RS,-1))*pow(RS,6) + 
5040*k2*pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(RS,6) - 1440*a1*b1*pow(M_E,k2*pow(RS,-1))*pow(RS,6) - 
5040*k1*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,6) - 1440*a1*b1*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(RS,6) - 
5040*k1*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(RS,6) + 1440*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(RS,6) + 
5040*k2*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(RS,6) + 5040*pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(RS,7) - 
5040*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,7) + 5040*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(RS,7) - 
5040*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(RS,7)) + 45360*pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(RS,8) - 
45360*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,8) + 45360*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RS,-1))*pow(RS,8) - 
45360*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RS,-1))*pow(RS,8));

    return res;
}

//  2
double BMSH_Integral_11_case_2_sub_1( double AS, double RS, double k1, double k2, double a1, double b1, double c1, double d1, double d )
{
double res = AS*pow(2*d,-1)*pow(M_E,-((2*d + k2)*pow(RS,-1)))*pow(RS,2)*
(d*k2*pow(d1,2)*pow(M_E,d*pow(RS,-1)) + d*RS*pow(d1,2)*pow(M_E,d*pow(RS,-1)) + 3*k2*RS*pow(d1,2)*pow(M_E,d*pow(RS,-1)) + 
d*k2*pow(d1,2)*pow(M_E,3*d*pow(RS,-1)) + d*RS*pow(d1,2)*pow(M_E,3*d*pow(RS,-1)) - 
3*k2*RS*pow(d1,2)*pow(M_E,3*d*pow(RS,-1)) - 20*a1*d1*RS*pow(d,4)*pow(M_E,k2*pow(RS,-1)) - 
4*a1*d1*pow(d,5)*pow(M_E,k2*pow(RS,-1)) - 16*RS*pow(a1,2)*pow(d,7)*pow(M_E,k2*pow(RS,-1)) - 
2*pow(a1,2)*pow(d,8)*pow(M_E,k2*pow(RS,-1)) - 4*d*RS*pow(d1,2)*pow(M_E,k2*pow(RS,-1)) - 
2*pow(d,2)*pow(d1,2)*pow(M_E,k2*pow(RS,-1)) + 4*a1*d1*RS*pow(d,4)*pow(M_E,(2*d + k2)*pow(RS,-1)) + 
2*RS*pow(a1,2)*pow(d,7)*pow(M_E,(2*d + k2)*pow(RS,-1)) + 2*d*RS*pow(d1,2)*pow(M_E,(2*d + k2)*pow(RS,-1)) + 
pow(d1,2)*pow(M_E,d*pow(RS,-1))*pow(k2,2) - pow(d1,2)*pow(M_E,3*d*pow(RS,-1))*pow(k2,2) + 
8*a1*d*d1*RS*pow(M_E,d*pow(RS,-1))*pow(k2,3) + 8*a1*d*d1*RS*pow(M_E,3*d*pow(RS,-1))*pow(k2,3) + 
2*a1*d*d1*pow(M_E,d*pow(RS,-1))*pow(k2,4) + 12*a1*d1*RS*pow(M_E,d*pow(RS,-1))*pow(k2,4) + 
2*a1*d*d1*pow(M_E,3*d*pow(RS,-1))*pow(k2,4) - 12*a1*d1*RS*pow(M_E,3*d*pow(RS,-1))*pow(k2,4) + 
2*a1*d1*pow(M_E,d*pow(RS,-1))*pow(k2,5) - 2*a1*d1*pow(M_E,3*d*pow(RS,-1))*pow(k2,5) + 
7*d*RS*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k2,6) + 7*d*RS*pow(a1,2)*pow(M_E,3*d*pow(RS,-1))*pow(k2,6) + 
d*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k2,7) + 9*RS*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k2,7) + 
d*pow(a1,2)*pow(M_E,3*d*pow(RS,-1))*pow(k2,7) - 9*RS*pow(a1,2)*pow(M_E,3*d*pow(RS,-1))*pow(k2,7) + 
pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k2,8) - pow(a1,2)*pow(M_E,3*d*pow(RS,-1))*pow(k2,8) + 
3*pow(d1,2)*pow(M_E,d*pow(RS,-1))*pow(RS,2) - 3*pow(d1,2)*pow(M_E,3*d*pow(RS,-1))*pow(RS,2) - 
72*a1*d1*pow(d,3)*pow(M_E,k2*pow(RS,-1))*pow(RS,2) - 105*pow(a1,2)*pow(d,6)*pow(M_E,k2*pow(RS,-1))*pow(RS,2) - 
3*pow(d1,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,2) + 24*a1*d1*pow(d,3)*pow(M_E,(2*d + k2)*pow(RS,-1))*pow(RS,2) + 
21*pow(a1,2)*pow(d,6)*pow(M_E,(2*d + k2)*pow(RS,-1))*pow(RS,2) + 3*pow(d1,2)*pow(M_E,(2*d + k2)*pow(RS,-1))*pow(RS,2) + 
24*a1*d*d1*pow(M_E,d*pow(RS,-1))*pow(k2,2)*pow(RS,2) + 24*a1*d*d1*pow(M_E,3*d*pow(RS,-1))*pow(k2,2)*pow(RS,2) + 
48*a1*d1*pow(M_E,d*pow(RS,-1))*pow(k2,3)*pow(RS,2) - 48*a1*d1*pow(M_E,3*d*pow(RS,-1))*pow(k2,3)*pow(RS,2) + 
42*d*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k2,5)*pow(RS,2) + 42*d*pow(a1,2)*pow(M_E,3*d*pow(RS,-1))*pow(k2,5)*pow(RS,2) + 
63*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k2,6)*pow(RS,2) - 63*pow(a1,2)*pow(M_E,3*d*pow(RS,-1))*pow(k2,6)*pow(RS,2) + 
48*a1*d*d1*k2*pow(M_E,d*pow(RS,-1))*pow(RS,3) + 48*a1*d*d1*k2*pow(M_E,3*d*pow(RS,-1))*pow(RS,3) - 
192*a1*d1*pow(d,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,3) - 588*pow(a1,2)*pow(d,5)*pow(M_E,k2*pow(RS,-1))*pow(RS,3) + 
96*a1*d1*pow(d,2)*pow(M_E,(2*d + k2)*pow(RS,-1))*pow(RS,3) + 
168*pow(a1,2)*pow(d,5)*pow(M_E,(2*d + k2)*pow(RS,-1))*pow(RS,3) + 144*a1*d1*pow(M_E,d*pow(RS,-1))*pow(k2,2)*pow(RS,3) - 
144*a1*d1*pow(M_E,3*d*pow(RS,-1))*pow(k2,2)*pow(RS,3) + 210*d*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k2,4)*pow(RS,3) + 
210*d*pow(a1,2)*pow(M_E,3*d*pow(RS,-1))*pow(k2,4)*pow(RS,3) + 378*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k2,5)*pow(RS,3) - 
378*pow(a1,2)*pow(M_E,3*d*pow(RS,-1))*pow(k2,5)*pow(RS,3) + 48*a1*d*d1*pow(M_E,d*pow(RS,-1))*pow(RS,4) + 
288*a1*d1*k2*pow(M_E,d*pow(RS,-1))*pow(RS,4) + 48*a1*d*d1*pow(M_E,3*d*pow(RS,-1))*pow(RS,4) - 
288*a1*d1*k2*pow(M_E,3*d*pow(RS,-1))*pow(RS,4) - 336*a1*d*d1*pow(M_E,k2*pow(RS,-1))*pow(RS,4) - 
2730*pow(a1,2)*pow(d,4)*pow(M_E,k2*pow(RS,-1))*pow(RS,4) + 240*a1*d*d1*pow(M_E,(2*d + k2)*pow(RS,-1))*pow(RS,4) + 
1050*pow(a1,2)*pow(d,4)*pow(M_E,(2*d + k2)*pow(RS,-1))*pow(RS,4) + 
840*d*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k2,3)*pow(RS,4) + 840*d*pow(a1,2)*pow(M_E,3*d*pow(RS,-1))*pow(k2,3)*pow(RS,4) + 
1890*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k2,4)*pow(RS,4) - 1890*pow(a1,2)*pow(M_E,3*d*pow(RS,-1))*pow(k2,4)*pow(RS,4) + 
pow(c1,2)*(-2*pow(d,4)*pow(M_E,k2*pow(RS,-1)) + 2*RS*pow(d,3)*(-4 + pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k2*pow(RS,-1)) + 
3*pow(d,2)*(-7 + 3*pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k2*pow(RS,-1))*pow(RS,2) + 
d*(-36*pow(M_E,k2*pow(RS,-1))*pow(RS,3) + 24*pow(M_E,(2*d + k2)*pow(RS,-1))*pow(RS,3) + 
pow(M_E,d*pow(RS,-1))*(3*RS*pow(k2,2) + pow(k2,3) + 6*k2*pow(RS,2) + 6*pow(RS,3)) + 
pow(M_E,3*d*pow(RS,-1))*(3*RS*pow(k2,2) + pow(k2,3) + 6*k2*pow(RS,2) + 6*pow(RS,3))) - 
(-1 + pow(M_E,2*d*pow(RS,-1)))*(-30*pow(M_E,k2*pow(RS,-1))*pow(RS,4) + 
pow(M_E,d*pow(RS,-1))*(5*RS*pow(k2,3) + pow(k2,4) + 15*pow(k2,2)*pow(RS,2) + 30*k2*pow(RS,3) + 30*pow(RS,4)))) + 
288*a1*d1*pow(M_E,d*pow(RS,-1))*pow(RS,5) - 288*a1*d1*pow(M_E,3*d*pow(RS,-1))*pow(RS,5) - 
288*a1*d1*pow(M_E,k2*pow(RS,-1))*pow(RS,5) - 10080*pow(a1,2)*pow(d,3)*pow(M_E,k2*pow(RS,-1))*pow(RS,5) + 
288*a1*d1*pow(M_E,(2*d + k2)*pow(RS,-1))*pow(RS,5) + 5040*pow(a1,2)*pow(d,3)*pow(M_E,(2*d + k2)*pow(RS,-1))*pow(RS,5) + 
2520*d*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k2,2)*pow(RS,5) + 
2520*d*pow(a1,2)*pow(M_E,3*d*pow(RS,-1))*pow(k2,2)*pow(RS,5) + 7560*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k2,3)*pow(RS,5) - 
7560*pow(a1,2)*pow(M_E,3*d*pow(RS,-1))*pow(k2,3)*pow(RS,5) + 5040*d*k2*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(RS,6) + 
5040*d*k2*pow(a1,2)*pow(M_E,3*d*pow(RS,-1))*pow(RS,6) - 27720*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,6) + 
17640*pow(a1,2)*pow(d,2)*pow(M_E,(2*d + k2)*pow(RS,-1))*pow(RS,6) + 
22680*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k2,2)*pow(RS,6) - 22680*pow(a1,2)*pow(M_E,3*d*pow(RS,-1))*pow(k2,2)*pow(RS,6) + 
pow(b1,2)*(-2*pow(d,6)*pow(M_E,k2*pow(RS,-1)) + 2*RS*pow(d,5)*(-6 + pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k2*pow(RS,-1)) + 
5*pow(d,4)*(-11 + 3*pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k2*pow(RS,-1))*pow(RS,2) + 
40*pow(d,3)*(-5 + 2*pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k2*pow(RS,-1))*pow(RS,3) + 
60*pow(d,2)*(-9 + 5*pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k2*pow(RS,-1))*pow(RS,4) + 
d*(-960*pow(M_E,k2*pow(RS,-1))*pow(RS,5) + 720*pow(M_E,(2*d + k2)*pow(RS,-1))*pow(RS,5) + 
pow(M_E,d*pow(RS,-1))*(5*RS*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RS,2) + 60*pow(k2,2)*pow(RS,3) + 
120*k2*pow(RS,4) + 120*pow(RS,5)) + 
pow(M_E,3*d*pow(RS,-1))*(5*RS*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RS,2) + 60*pow(k2,2)*pow(RS,3) + 
120*k2*pow(RS,4) + 120*pow(RS,5))) - 
(-1 + pow(M_E,2*d*pow(RS,-1)))*(-840*pow(M_E,k2*pow(RS,-1))*pow(RS,6) + 
pow(M_E,d*pow(RS,-1))*(7*RS*pow(k2,5) + pow(k2,6) + 35*pow(k2,4)*pow(RS,2) + 140*pow(k2,3)*pow(RS,3) + 
420*pow(k2,2)*pow(RS,4) + 840*k2*pow(RS,5) + 840*pow(RS,6)))) - 
2*c1*(d1*(2*pow(d,3)*pow(M_E,k2*pow(RS,-1)) - 2*RS*pow(d,2)*(-3 + pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k2*pow(RS,-1)) - 
d*(-10*pow(M_E,k2*pow(RS,-1))*pow(RS,2) + 6*pow(M_E,(2*d + k2)*pow(RS,-1))*pow(RS,2) + 
pow(M_E,d*pow(RS,-1))*(2*k2*RS + pow(k2,2) + 2*pow(RS,2)) + 
pow(M_E,3*d*pow(RS,-1))*(2*k2*RS + pow(k2,2) + 2*pow(RS,2))) + 
(-1 + pow(M_E,2*d*pow(RS,-1)))*(-8*pow(M_E,k2*pow(RS,-1))*pow(RS,3) + 
pow(M_E,d*pow(RS,-1))*(4*RS*pow(k2,2) + pow(k2,3) + 8*k2*pow(RS,2) + 8*pow(RS,3)))) + 
b1*(2*pow(d,5)*pow(M_E,k2*pow(RS,-1)) - 2*RS*pow(d,4)*(-5 + pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k2*pow(RS,-1)) - 
12*pow(d,3)*(-3 + pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k2*pow(RS,-1))*pow(RS,2) - 
48*pow(d,2)*(-2 + pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k2*pow(RS,-1))*pow(RS,3) - 
d*(-168*pow(M_E,k2*pow(RS,-1))*pow(RS,4) + 120*pow(M_E,(2*d + k2)*pow(RS,-1))*pow(RS,4) + 
pow(M_E,d*pow(RS,-1))*(4*RS*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RS,2) + 24*k2*pow(RS,3) + 24*pow(RS,4)) + 
pow(M_E,3*d*pow(RS,-1))*(4*RS*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RS,2) + 24*k2*pow(RS,3) + 24*pow(RS,4)))\
+ (-1 + pow(M_E,2*d*pow(RS,-1)))*(-144*pow(M_E,k2*pow(RS,-1))*pow(RS,5) + 
pow(M_E,d*pow(RS,-1))*(6*RS*pow(k2,4) + pow(k2,5) + 24*pow(k2,3)*pow(RS,2) + 72*pow(k2,2)*pow(RS,3) + 
144*k2*pow(RS,4) + 144*pow(RS,5)))) + 
a1*(2*pow(d,6)*pow(M_E,k2*pow(RS,-1)) - 2*RS*pow(d,5)*(-6 + pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k2*pow(RS,-1)) - 
5*pow(d,4)*(-11 + 3*pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k2*pow(RS,-1))*pow(RS,2) - 
40*pow(d,3)*(-5 + 2*pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k2*pow(RS,-1))*pow(RS,3) - 
60*pow(d,2)*(-9 + 5*pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k2*pow(RS,-1))*pow(RS,4) - 
d*(-960*pow(M_E,k2*pow(RS,-1))*pow(RS,5) + 720*pow(M_E,(2*d + k2)*pow(RS,-1))*pow(RS,5) + 
pow(M_E,d*pow(RS,-1))*(5*RS*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RS,2) + 60*pow(k2,2)*pow(RS,3) + 
120*k2*pow(RS,4) + 120*pow(RS,5)) + 
pow(M_E,3*d*pow(RS,-1))*(5*RS*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RS,2) + 60*pow(k2,2)*pow(RS,3) + 
120*k2*pow(RS,4) + 120*pow(RS,5))) + 
(-1 + pow(M_E,2*d*pow(RS,-1)))*(-840*pow(M_E,k2*pow(RS,-1))*pow(RS,6) + 
pow(M_E,d*pow(RS,-1))*(7*RS*pow(k2,5) + pow(k2,6) + 35*pow(k2,4)*pow(RS,2) + 140*pow(k2,3)*pow(RS,3) + 
420*pow(k2,2)*pow(RS,4) + 840*k2*pow(RS,5) + 840*pow(RS,6))))) + 
5040*d*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(RS,7) + 45360*k2*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(RS,7) + 
5040*d*pow(a1,2)*pow(M_E,3*d*pow(RS,-1))*pow(RS,7) - 45360*k2*pow(a1,2)*pow(M_E,3*d*pow(RS,-1))*pow(RS,7) - 
50400*d*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,7) + 40320*d*pow(a1,2)*pow(M_E,(2*d + k2)*pow(RS,-1))*pow(RS,7) - 
2*b1*(d1*(2*pow(d,4)*pow(M_E,k2*pow(RS,-1)) - 2*RS*pow(d,3)*(-4 + pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k2*pow(RS,-1)) - 
3*pow(d,2)*(-7 + 3*pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k2*pow(RS,-1))*pow(RS,2) - 
d*(-36*pow(M_E,k2*pow(RS,-1))*pow(RS,3) + 24*pow(M_E,(2*d + k2)*pow(RS,-1))*pow(RS,3) + 
pow(M_E,d*pow(RS,-1))*(3*RS*pow(k2,2) + pow(k2,3) + 6*k2*pow(RS,2) + 6*pow(RS,3)) + 
pow(M_E,3*d*pow(RS,-1))*(3*RS*pow(k2,2) + pow(k2,3) + 6*k2*pow(RS,2) + 6*pow(RS,3))) + 
(-1 + pow(M_E,2*d*pow(RS,-1)))*(-30*pow(M_E,k2*pow(RS,-1))*pow(RS,4) + 
pow(M_E,d*pow(RS,-1))*(5*RS*pow(k2,3) + pow(k2,4) + 15*pow(k2,2)*pow(RS,2) + 30*k2*pow(RS,3) + 30*pow(RS,4))))\
+ a1*(2*pow(d,7)*pow(M_E,k2*pow(RS,-1)) - 2*RS*pow(d,6)*(-7 + pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k2*pow(RS,-1)) - 
6*pow(d,5)*(-13 + 3*pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k2*pow(RS,-1))*pow(RS,2) - 
120*pow(d,4)*(-3 + pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k2*pow(RS,-1))*pow(RS,3) - 
120*pow(d,3)*(-11 + 5*pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k2*pow(RS,-1))*pow(RS,4) - 
720*pow(d,2)*(-5 + 3*pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k2*pow(RS,-1))*pow(RS,5) - 
d*(-6480*pow(M_E,k2*pow(RS,-1))*pow(RS,6) + 5040*pow(M_E,(2*d + k2)*pow(RS,-1))*pow(RS,6) + 
pow(M_E,d*pow(RS,-1))*(6*RS*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RS,2) + 120*pow(k2,3)*pow(RS,3) + 
360*pow(k2,2)*pow(RS,4) + 720*k2*pow(RS,5) + 720*pow(RS,6)) + 
pow(M_E,3*d*pow(RS,-1))*(6*RS*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RS,2) + 120*pow(k2,3)*pow(RS,3) + 
360*pow(k2,2)*pow(RS,4) + 720*k2*pow(RS,5) + 720*pow(RS,6))) + 
(-1 + pow(M_E,2*d*pow(RS,-1)))*(-5760*pow(M_E,k2*pow(RS,-1))*pow(RS,7) + 
pow(M_E,d*pow(RS,-1))*(8*RS*pow(k2,6) + pow(k2,7) + 48*pow(k2,5)*pow(RS,2) + 240*pow(k2,4)*pow(RS,3) + 
960*pow(k2,3)*pow(RS,4) + 2880*pow(k2,2)*pow(RS,5) + 5760*k2*pow(RS,6) + 5760*pow(RS,7))))) + 
45360*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(RS,8) - 45360*pow(a1,2)*pow(M_E,3*d*pow(RS,-1))*pow(RS,8) - 
45360*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,8) + 45360*pow(a1,2)*pow(M_E,(2*d + k2)*pow(RS,-1))*pow(RS,8));

    return res;
}


// 3
double BMSH_Integral_11_case_2_sub_2( double AS, double RS, double k1, double k2, double a1, double b1, double c1, double d1, double d )
{
double res = AS*pow(2*d,-1)*pow(M_E,-((2*d + k1)*pow(RS,-1)))*pow(RS,2)*
(-(d*k1*pow(d1,2)*pow(M_E,d*pow(RS,-1))) - d*RS*pow(d1,2)*pow(M_E,d*pow(RS,-1)) - 3*k1*RS*pow(d1,2)*pow(M_E,d*pow(RS,-1)) + 
20*a1*d1*RS*pow(d,4)*pow(M_E,k1*pow(RS,-1)) + 4*a1*d1*pow(d,5)*pow(M_E,k1*pow(RS,-1)) + 
16*RS*pow(a1,2)*pow(d,7)*pow(M_E,k1*pow(RS,-1)) + 2*pow(a1,2)*pow(d,8)*pow(M_E,k1*pow(RS,-1)) + 
4*d*RS*pow(d1,2)*pow(M_E,k1*pow(RS,-1)) + 2*pow(d,2)*pow(d1,2)*pow(M_E,k1*pow(RS,-1)) + 
4*a1*d1*RS*pow(d,4)*pow(M_E,(2*d + k1)*pow(RS,-1)) + 2*RS*pow(a1,2)*pow(d,7)*pow(M_E,(2*d + k1)*pow(RS,-1)) + 
2*d*RS*pow(d1,2)*pow(M_E,(2*d + k1)*pow(RS,-1)) - d*k1*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RS,-1)) + 
d*RS*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RS,-1)) - 3*k1*RS*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RS,-1)) - 
pow(d1,2)*pow(M_E,d*pow(RS,-1))*pow(k1,2) + pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RS,-1))*pow(k1,2) - 
8*a1*d*d1*RS*pow(M_E,d*pow(RS,-1))*pow(k1,3) + 8*a1*d*d1*RS*pow(M_E,(d + 2*k1)*pow(RS,-1))*pow(k1,3) - 
2*a1*d*d1*pow(M_E,d*pow(RS,-1))*pow(k1,4) - 12*a1*d1*RS*pow(M_E,d*pow(RS,-1))*pow(k1,4) - 
2*a1*d*d1*pow(M_E,(d + 2*k1)*pow(RS,-1))*pow(k1,4) - 12*a1*d1*RS*pow(M_E,(d + 2*k1)*pow(RS,-1))*pow(k1,4) - 
2*a1*d1*pow(M_E,d*pow(RS,-1))*pow(k1,5) + 2*a1*d1*pow(M_E,(d + 2*k1)*pow(RS,-1))*pow(k1,5) - 
7*d*RS*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k1,6) + 7*d*RS*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RS,-1))*pow(k1,6) - 
d*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k1,7) - 9*RS*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k1,7) - 
d*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RS,-1))*pow(k1,7) - 9*RS*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RS,-1))*pow(k1,7) - 
pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k1,8) + pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RS,-1))*pow(k1,8) - 
3*pow(d1,2)*pow(M_E,d*pow(RS,-1))*pow(RS,2) + 72*a1*d1*pow(d,3)*pow(M_E,k1*pow(RS,-1))*pow(RS,2) + 
105*pow(a1,2)*pow(d,6)*pow(M_E,k1*pow(RS,-1))*pow(RS,2) + 3*pow(d1,2)*pow(M_E,k1*pow(RS,-1))*pow(RS,2) - 
24*a1*d1*pow(d,3)*pow(M_E,(2*d + k1)*pow(RS,-1))*pow(RS,2) - 
21*pow(a1,2)*pow(d,6)*pow(M_E,(2*d + k1)*pow(RS,-1))*pow(RS,2) - 3*pow(d1,2)*pow(M_E,(2*d + k1)*pow(RS,-1))*pow(RS,2) + 
3*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RS,-1))*pow(RS,2) - 24*a1*d*d1*pow(M_E,d*pow(RS,-1))*pow(k1,2)*pow(RS,2) - 
24*a1*d*d1*pow(M_E,(d + 2*k1)*pow(RS,-1))*pow(k1,2)*pow(RS,2) - 48*a1*d1*pow(M_E,d*pow(RS,-1))*pow(k1,3)*pow(RS,2) + 
48*a1*d1*pow(M_E,(d + 2*k1)*pow(RS,-1))*pow(k1,3)*pow(RS,2) - 42*d*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k1,5)*pow(RS,2) - 
42*d*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RS,-1))*pow(k1,5)*pow(RS,2) - 
63*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k1,6)*pow(RS,2) + 63*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RS,-1))*pow(k1,6)*pow(RS,2) - 
48*a1*d*d1*k1*pow(M_E,d*pow(RS,-1))*pow(RS,3) + 192*a1*d1*pow(d,2)*pow(M_E,k1*pow(RS,-1))*pow(RS,3) + 
588*pow(a1,2)*pow(d,5)*pow(M_E,k1*pow(RS,-1))*pow(RS,3) + 96*a1*d1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RS,-1))*pow(RS,3) + 
168*pow(a1,2)*pow(d,5)*pow(M_E,(2*d + k1)*pow(RS,-1))*pow(RS,3) + 48*a1*d*d1*k1*pow(M_E,(d + 2*k1)*pow(RS,-1))*pow(RS,3) - 
144*a1*d1*pow(M_E,d*pow(RS,-1))*pow(k1,2)*pow(RS,3) - 144*a1*d1*pow(M_E,(d + 2*k1)*pow(RS,-1))*pow(k1,2)*pow(RS,3) - 
210*d*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k1,4)*pow(RS,3) + 
210*d*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RS,-1))*pow(k1,4)*pow(RS,3) - 
378*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k1,5)*pow(RS,3) - 
378*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RS,-1))*pow(k1,5)*pow(RS,3) - 48*a1*d*d1*pow(M_E,d*pow(RS,-1))*pow(RS,4) - 
288*a1*d1*k1*pow(M_E,d*pow(RS,-1))*pow(RS,4) + 336*a1*d*d1*pow(M_E,k1*pow(RS,-1))*pow(RS,4) + 
2730*pow(a1,2)*pow(d,4)*pow(M_E,k1*pow(RS,-1))*pow(RS,4) - 240*a1*d*d1*pow(M_E,(2*d + k1)*pow(RS,-1))*pow(RS,4) - 
1050*pow(a1,2)*pow(d,4)*pow(M_E,(2*d + k1)*pow(RS,-1))*pow(RS,4) - 48*a1*d*d1*pow(M_E,(d + 2*k1)*pow(RS,-1))*pow(RS,4) + 
288*a1*d1*k1*pow(M_E,(d + 2*k1)*pow(RS,-1))*pow(RS,4) - 840*d*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k1,3)*pow(RS,4) - 
840*d*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RS,-1))*pow(k1,3)*pow(RS,4) - 
1890*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k1,4)*pow(RS,4) + 
1890*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RS,-1))*pow(k1,4)*pow(RS,4) - 
pow(c1,2)*(-2*pow(d,4)*pow(M_E,k1*pow(RS,-1)) - 2*RS*pow(d,3)*(4 + pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k1*pow(RS,-1)) + 
3*pow(d,2)*(-7 + 3*pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k1*pow(RS,-1))*pow(RS,2) + 
d*(pow(M_E,(d + 2*k1)*pow(RS,-1))*(-3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) - 6*pow(RS,3)) - 
36*pow(M_E,k1*pow(RS,-1))*pow(RS,3) - 24*pow(M_E,(2*d + k1)*pow(RS,-1))*pow(RS,3) + 
pow(M_E,d*pow(RS,-1))*(3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) + 6*pow(RS,3))) - 
30*pow(M_E,k1*pow(RS,-1))*pow(RS,4) + 30*pow(M_E,(2*d + k1)*pow(RS,-1))*pow(RS,4) - 
pow(M_E,(d + 2*k1)*pow(RS,-1))*(-5*RS*pow(k1,3) + pow(k1,4) + 15*pow(k1,2)*pow(RS,2) - 30*k1*pow(RS,3) + 
30*pow(RS,4)) + pow(M_E,d*pow(RS,-1))*(5*RS*pow(k1,3) + pow(k1,4) + 15*pow(k1,2)*pow(RS,2) + 30*k1*pow(RS,3) + 
30*pow(RS,4))) - 288*a1*d1*pow(M_E,d*pow(RS,-1))*pow(RS,5) + 288*a1*d1*pow(M_E,k1*pow(RS,-1))*pow(RS,5) + 
10080*pow(a1,2)*pow(d,3)*pow(M_E,k1*pow(RS,-1))*pow(RS,5) + 288*a1*d1*pow(M_E,(2*d + k1)*pow(RS,-1))*pow(RS,5) + 
5040*pow(a1,2)*pow(d,3)*pow(M_E,(2*d + k1)*pow(RS,-1))*pow(RS,5) - 288*a1*d1*pow(M_E,(d + 2*k1)*pow(RS,-1))*pow(RS,5) - 
2520*d*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k1,2)*pow(RS,5) + 
2520*d*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RS,-1))*pow(k1,2)*pow(RS,5) - 
7560*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k1,3)*pow(RS,5) - 
7560*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RS,-1))*pow(k1,3)*pow(RS,5) - 5040*d*k1*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(RS,6) + 
27720*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RS,-1))*pow(RS,6) - 
17640*pow(a1,2)*pow(d,2)*pow(M_E,(2*d + k1)*pow(RS,-1))*pow(RS,6) - 
5040*d*k1*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RS,-1))*pow(RS,6) - 22680*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k1,2)*pow(RS,6) + 
22680*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RS,-1))*pow(k1,2)*pow(RS,6) - 
pow(b1,2)*(-2*pow(d,6)*pow(M_E,k1*pow(RS,-1)) - 2*RS*pow(d,5)*(6 + pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k1*pow(RS,-1)) + 
5*pow(d,4)*(-11 + 3*pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k1*pow(RS,-1))*pow(RS,2) - 
40*pow(d,3)*(5 + 2*pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k1*pow(RS,-1))*pow(RS,3) + 
60*pow(d,2)*(-9 + 5*pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k1*pow(RS,-1))*pow(RS,4) + 
d*(pow(M_E,(d + 2*k1)*pow(RS,-1))*(-5*RS*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RS,2) - 60*pow(k1,2)*pow(RS,3) + 
120*k1*pow(RS,4) - 120*pow(RS,5)) - 960*pow(M_E,k1*pow(RS,-1))*pow(RS,5) - 
720*pow(M_E,(2*d + k1)*pow(RS,-1))*pow(RS,5) + 
pow(M_E,d*pow(RS,-1))*(5*RS*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RS,2) + 60*pow(k1,2)*pow(RS,3) + 
120*k1*pow(RS,4) + 120*pow(RS,5))) - 840*pow(M_E,k1*pow(RS,-1))*pow(RS,6) + 
840*pow(M_E,(2*d + k1)*pow(RS,-1))*pow(RS,6) - 
pow(M_E,(d + 2*k1)*pow(RS,-1))*(-7*RS*pow(k1,5) + pow(k1,6) + 35*pow(k1,4)*pow(RS,2) - 140*pow(k1,3)*pow(RS,3) + 
420*pow(k1,2)*pow(RS,4) - 840*k1*pow(RS,5) + 840*pow(RS,6)) + 
pow(M_E,d*pow(RS,-1))*(7*RS*pow(k1,5) + pow(k1,6) + 35*pow(k1,4)*pow(RS,2) + 140*pow(k1,3)*pow(RS,3) + 
420*pow(k1,2)*pow(RS,4) + 840*k1*pow(RS,5) + 840*pow(RS,6))) - 
2*c1*(d1*(-2*pow(d,3)*pow(M_E,k1*pow(RS,-1)) - 2*RS*pow(d,2)*(3 + pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k1*pow(RS,-1)) + 
d*(-10*pow(M_E,k1*pow(RS,-1))*pow(RS,2) + 6*pow(M_E,(2*d + k1)*pow(RS,-1))*pow(RS,2) + 
pow(M_E,(d + 2*k1)*pow(RS,-1))*(-2*k1*RS + pow(k1,2) + 2*pow(RS,2)) + 
pow(M_E,d*pow(RS,-1))*(2*k1*RS + pow(k1,2) + 2*pow(RS,2))) - 8*pow(M_E,k1*pow(RS,-1))*pow(RS,3) - 
8*pow(M_E,(2*d + k1)*pow(RS,-1))*pow(RS,3) + 
pow(M_E,(d + 2*k1)*pow(RS,-1))*(4*RS*pow(k1,2) - pow(k1,3) - 8*k1*pow(RS,2) + 8*pow(RS,3)) + 
pow(M_E,d*pow(RS,-1))*(4*RS*pow(k1,2) + pow(k1,3) + 8*k1*pow(RS,2) + 8*pow(RS,3))) + 
b1*(-2*pow(d,5)*pow(M_E,k1*pow(RS,-1)) - 2*RS*pow(d,4)*(5 + pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k1*pow(RS,-1)) + 
12*pow(d,3)*(-3 + pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k1*pow(RS,-1))*pow(RS,2) - 
48*pow(d,2)*(2 + pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k1*pow(RS,-1))*pow(RS,3) + 
d*(-168*pow(M_E,k1*pow(RS,-1))*pow(RS,4) + 120*pow(M_E,(2*d + k1)*pow(RS,-1))*pow(RS,4) + 
pow(M_E,(d + 2*k1)*pow(RS,-1))*(-4*RS*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RS,2) - 24*k1*pow(RS,3) + 
24*pow(RS,4)) + pow(M_E,d*pow(RS,-1))*
(4*RS*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RS,2) + 24*k1*pow(RS,3) + 24*pow(RS,4))) - 
144*pow(M_E,k1*pow(RS,-1))*pow(RS,5) - 144*pow(M_E,(2*d + k1)*pow(RS,-1))*pow(RS,5) + 
pow(M_E,(d + 2*k1)*pow(RS,-1))*(6*RS*pow(k1,4) - pow(k1,5) - 24*pow(k1,3)*pow(RS,2) + 72*pow(k1,2)*pow(RS,3) - 
144*k1*pow(RS,4) + 144*pow(RS,5)) + 
pow(M_E,d*pow(RS,-1))*(6*RS*pow(k1,4) + pow(k1,5) + 24*pow(k1,3)*pow(RS,2) + 72*pow(k1,2)*pow(RS,3) + 
144*k1*pow(RS,4) + 144*pow(RS,5))) + 
a1*(-2*pow(d,6)*pow(M_E,k1*pow(RS,-1)) - 2*RS*pow(d,5)*(6 + pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k1*pow(RS,-1)) + 
5*pow(d,4)*(-11 + 3*pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k1*pow(RS,-1))*pow(RS,2) - 
40*pow(d,3)*(5 + 2*pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k1*pow(RS,-1))*pow(RS,3) + 
60*pow(d,2)*(-9 + 5*pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k1*pow(RS,-1))*pow(RS,4) + 
d*(pow(M_E,(d + 2*k1)*pow(RS,-1))*(-5*RS*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RS,2) - 60*pow(k1,2)*pow(RS,3) + 
120*k1*pow(RS,4) - 120*pow(RS,5)) - 960*pow(M_E,k1*pow(RS,-1))*pow(RS,5) - 
720*pow(M_E,(2*d + k1)*pow(RS,-1))*pow(RS,5) + 
pow(M_E,d*pow(RS,-1))*(5*RS*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RS,2) + 60*pow(k1,2)*pow(RS,3) + 
120*k1*pow(RS,4) + 120*pow(RS,5))) - 840*pow(M_E,k1*pow(RS,-1))*pow(RS,6) + 
840*pow(M_E,(2*d + k1)*pow(RS,-1))*pow(RS,6) - 
pow(M_E,(d + 2*k1)*pow(RS,-1))*(-7*RS*pow(k1,5) + pow(k1,6) + 35*pow(k1,4)*pow(RS,2) - 140*pow(k1,3)*pow(RS,3) + 
420*pow(k1,2)*pow(RS,4) - 840*k1*pow(RS,5) + 840*pow(RS,6)) + 
pow(M_E,d*pow(RS,-1))*(7*RS*pow(k1,5) + pow(k1,6) + 35*pow(k1,4)*pow(RS,2) + 140*pow(k1,3)*pow(RS,3) + 
420*pow(k1,2)*pow(RS,4) + 840*k1*pow(RS,5) + 840*pow(RS,6)))) - 
5040*d*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(RS,7) - 45360*k1*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(RS,7) + 
50400*d*pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(RS,7) + 40320*d*pow(a1,2)*pow(M_E,(2*d + k1)*pow(RS,-1))*pow(RS,7) + 
5040*d*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RS,-1))*pow(RS,7) - 45360*k1*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RS,-1))*pow(RS,7) - 
2*b1*(d1*(-2*pow(d,4)*pow(M_E,k1*pow(RS,-1)) - 2*RS*pow(d,3)*(4 + pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k1*pow(RS,-1)) + 
3*pow(d,2)*(-7 + 3*pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k1*pow(RS,-1))*pow(RS,2) + 
d*(pow(M_E,(d + 2*k1)*pow(RS,-1))*(-3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) - 6*pow(RS,3)) - 
36*pow(M_E,k1*pow(RS,-1))*pow(RS,3) - 24*pow(M_E,(2*d + k1)*pow(RS,-1))*pow(RS,3) + 
pow(M_E,d*pow(RS,-1))*(3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) + 6*pow(RS,3))) - 
30*pow(M_E,k1*pow(RS,-1))*pow(RS,4) + 30*pow(M_E,(2*d + k1)*pow(RS,-1))*pow(RS,4) - 
pow(M_E,(d + 2*k1)*pow(RS,-1))*(-5*RS*pow(k1,3) + pow(k1,4) + 15*pow(k1,2)*pow(RS,2) - 30*k1*pow(RS,3) + 
30*pow(RS,4)) + pow(M_E,d*pow(RS,-1))*
(5*RS*pow(k1,3) + pow(k1,4) + 15*pow(k1,2)*pow(RS,2) + 30*k1*pow(RS,3) + 30*pow(RS,4))) + 
a1*(-2*pow(d,7)*pow(M_E,k1*pow(RS,-1)) - 2*RS*pow(d,6)*(7 + pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k1*pow(RS,-1)) + 
6*pow(d,5)*(-13 + 3*pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k1*pow(RS,-1))*pow(RS,2) - 
120*pow(d,4)*(3 + pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k1*pow(RS,-1))*pow(RS,3) + 
120*pow(d,3)*(-11 + 5*pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k1*pow(RS,-1))*pow(RS,4) - 
720*pow(d,2)*(5 + 3*pow(M_E,2*d*pow(RS,-1)))*pow(M_E,k1*pow(RS,-1))*pow(RS,5) + 
d*(-6480*pow(M_E,k1*pow(RS,-1))*pow(RS,6) + 5040*pow(M_E,(2*d + k1)*pow(RS,-1))*pow(RS,6) + 
pow(M_E,(d + 2*k1)*pow(RS,-1))*(-6*RS*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RS,2) - 
120*pow(k1,3)*pow(RS,3) + 360*pow(k1,2)*pow(RS,4) - 720*k1*pow(RS,5) + 720*pow(RS,6)) + 
pow(M_E,d*pow(RS,-1))*(6*RS*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RS,2) + 120*pow(k1,3)*pow(RS,3) + 
360*pow(k1,2)*pow(RS,4) + 720*k1*pow(RS,5) + 720*pow(RS,6))) - 5760*pow(M_E,k1*pow(RS,-1))*pow(RS,7) - 
5760*pow(M_E,(2*d + k1)*pow(RS,-1))*pow(RS,7) + 
pow(M_E,(d + 2*k1)*pow(RS,-1))*(8*RS*pow(k1,6) - pow(k1,7) - 48*pow(k1,5)*pow(RS,2) + 240*pow(k1,4)*pow(RS,3) - 
960*pow(k1,3)*pow(RS,4) + 2880*pow(k1,2)*pow(RS,5) - 5760*k1*pow(RS,6) + 5760*pow(RS,7)) + 
pow(M_E,d*pow(RS,-1))*(8*RS*pow(k1,6) + pow(k1,7) + 48*pow(k1,5)*pow(RS,2) + 240*pow(k1,4)*pow(RS,3) + 
960*pow(k1,3)*pow(RS,4) + 2880*pow(k1,2)*pow(RS,5) + 5760*k1*pow(RS,6) + 5760*pow(RS,7)))) - 
45360*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(RS,8) + 45360*pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(RS,8) - 
45360*pow(a1,2)*pow(M_E,(2*d + k1)*pow(RS,-1))*pow(RS,8) + 45360*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RS,-1))*pow(RS,8));

    return res;
}



// 4
double BMSH_Integral_11_case_3( double AS, double RS, double k1, double k2, double a1, double b1, double c1, double d1, double d )
{
double res = AS*pow(2*d,-1)*pow(M_E,-((d + k1 + k2)*pow(RS,-1)))*pow(RS,2)*
(d*pow(d1,2)*(1 + pow(M_E,2*d*pow(RS,-1)))*((k2 + RS)*pow(M_E,k1*pow(RS,-1)) - (k1 + RS)*pow(M_E,k2*pow(RS,-1))) + 
RS*pow(d1,2)*(-1 + pow(M_E,2*d*pow(RS,-1)))*(-((k2 + RS)*pow(M_E,k1*pow(RS,-1))) + (k1 + RS)*pow(M_E,k2*pow(RS,-1))) + 
2*c1*d1*RS*(-1 + pow(M_E,2*d*pow(RS,-1)))*(pow(M_E,k2*pow(RS,-1))*(2*k1*RS + pow(k1,2) + 2*pow(RS,2)) - 
pow(M_E,k1*pow(RS,-1))*(2*k2*RS + pow(k2,2) + 2*pow(RS,2))) + 
pow(d1,2)*(-1 + pow(M_E,2*d*pow(RS,-1)))*(pow(M_E,k2*pow(RS,-1))*(2*k1*RS + pow(k1,2) + 2*pow(RS,2)) - 
pow(M_E,k1*pow(RS,-1))*(2*k2*RS + pow(k2,2) + 2*pow(RS,2))) + 
2*c1*d*d1*(1 + pow(M_E,2*d*pow(RS,-1)))*(-(pow(M_E,k2*pow(RS,-1))*(2*k1*RS + pow(k1,2) + 2*pow(RS,2))) + 
pow(M_E,k1*pow(RS,-1))*(2*k2*RS + pow(k2,2) + 2*pow(RS,2))) + 
2*c1*d1*(-1 + pow(M_E,2*d*pow(RS,-1)))*(pow(M_E,k2*pow(RS,-1))*
(3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) + 6*pow(RS,3)) - 
pow(M_E,k1*pow(RS,-1))*(3*RS*pow(k2,2) + pow(k2,3) + 6*k2*pow(RS,2) + 6*pow(RS,3))) + 
2*b1*d1*RS*(-1 + pow(M_E,2*d*pow(RS,-1)))*(pow(M_E,k2*pow(RS,-1))*
(3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) + 6*pow(RS,3)) - 
pow(M_E,k1*pow(RS,-1))*(3*RS*pow(k2,2) + pow(k2,3) + 6*k2*pow(RS,2) + 6*pow(RS,3))) + 
RS*pow(c1,2)*(-1 + pow(M_E,2*d*pow(RS,-1)))*(pow(M_E,k2*pow(RS,-1))*
(3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) + 6*pow(RS,3)) - 
pow(M_E,k1*pow(RS,-1))*(3*RS*pow(k2,2) + pow(k2,3) + 6*k2*pow(RS,2) + 6*pow(RS,3))) + 
2*b1*d*d1*(1 + pow(M_E,2*d*pow(RS,-1)))*(-(pow(M_E,k2*pow(RS,-1))*
(3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) + 6*pow(RS,3))) + 
pow(M_E,k1*pow(RS,-1))*(3*RS*pow(k2,2) + pow(k2,3) + 6*k2*pow(RS,2) + 6*pow(RS,3))) + 
d*pow(c1,2)*(1 + pow(M_E,2*d*pow(RS,-1)))*(-(pow(M_E,k2*pow(RS,-1))*
(3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) + 6*pow(RS,3))) + 
pow(M_E,k1*pow(RS,-1))*(3*RS*pow(k2,2) + pow(k2,3) + 6*k2*pow(RS,2) + 6*pow(RS,3))) + 
2*b1*d1*(-1 + pow(M_E,2*d*pow(RS,-1)))*(pow(M_E,k2*pow(RS,-1))*
(4*RS*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RS,2) + 24*k1*pow(RS,3) + 24*pow(RS,4)) - 
pow(M_E,k1*pow(RS,-1))*(4*RS*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RS,2) + 24*k2*pow(RS,3) + 24*pow(RS,4))) + 
2*b1*c1*RS*(-1 + pow(M_E,2*d*pow(RS,-1)))*(pow(M_E,k2*pow(RS,-1))*
(4*RS*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RS,2) + 24*k1*pow(RS,3) + 24*pow(RS,4)) - 
pow(M_E,k1*pow(RS,-1))*(4*RS*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RS,2) + 24*k2*pow(RS,3) + 24*pow(RS,4))) + 
2*a1*d1*RS*(-1 + pow(M_E,2*d*pow(RS,-1)))*(pow(M_E,k2*pow(RS,-1))*
(4*RS*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RS,2) + 24*k1*pow(RS,3) + 24*pow(RS,4)) - 
pow(M_E,k1*pow(RS,-1))*(4*RS*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RS,2) + 24*k2*pow(RS,3) + 24*pow(RS,4))) + 
pow(c1,2)*(-1 + pow(M_E,2*d*pow(RS,-1)))*(pow(M_E,k2*pow(RS,-1))*
(4*RS*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RS,2) + 24*k1*pow(RS,3) + 24*pow(RS,4)) - 
pow(M_E,k1*pow(RS,-1))*(4*RS*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RS,2) + 24*k2*pow(RS,3) + 24*pow(RS,4))) + 
2*b1*c1*d*(1 + pow(M_E,2*d*pow(RS,-1)))*(-(pow(M_E,k2*pow(RS,-1))*
(4*RS*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RS,2) + 24*k1*pow(RS,3) + 24*pow(RS,4))) + 
pow(M_E,k1*pow(RS,-1))*(4*RS*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RS,2) + 24*k2*pow(RS,3) + 24*pow(RS,4))) + 
2*a1*d*d1*(1 + pow(M_E,2*d*pow(RS,-1)))*(-(pow(M_E,k2*pow(RS,-1))*
(4*RS*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RS,2) + 24*k1*pow(RS,3) + 24*pow(RS,4))) + 
pow(M_E,k1*pow(RS,-1))*(4*RS*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RS,2) + 24*k2*pow(RS,3) + 24*pow(RS,4))) + 
2*b1*c1*(-1 + pow(M_E,2*d*pow(RS,-1)))*(pow(M_E,k2*pow(RS,-1))*
(5*RS*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RS,2) + 60*pow(k1,2)*pow(RS,3) + 120*k1*pow(RS,4) + 120*pow(RS,5))\
- pow(M_E,k1*pow(RS,-1))*(5*RS*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RS,2) + 60*pow(k2,2)*pow(RS,3) + 
120*k2*pow(RS,4) + 120*pow(RS,5))) + 2*a1*d1*(-1 + pow(M_E,2*d*pow(RS,-1)))*
(pow(M_E,k2*pow(RS,-1))*(5*RS*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RS,2) + 60*pow(k1,2)*pow(RS,3) + 
120*k1*pow(RS,4) + 120*pow(RS,5)) - pow(M_E,k1*pow(RS,-1))*
(5*RS*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RS,2) + 60*pow(k2,2)*pow(RS,3) + 120*k2*pow(RS,4) + 120*pow(RS,5)))\
+ 2*a1*c1*RS*(-1 + pow(M_E,2*d*pow(RS,-1)))*(pow(M_E,k2*pow(RS,-1))*
(5*RS*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RS,2) + 60*pow(k1,2)*pow(RS,3) + 120*k1*pow(RS,4) + 120*pow(RS,5))\
- pow(M_E,k1*pow(RS,-1))*(5*RS*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RS,2) + 60*pow(k2,2)*pow(RS,3) + 
120*k2*pow(RS,4) + 120*pow(RS,5))) + RS*pow(b1,2)*(-1 + pow(M_E,2*d*pow(RS,-1)))*
(pow(M_E,k2*pow(RS,-1))*(5*RS*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RS,2) + 60*pow(k1,2)*pow(RS,3) + 
120*k1*pow(RS,4) + 120*pow(RS,5)) - pow(M_E,k1*pow(RS,-1))*
(5*RS*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RS,2) + 60*pow(k2,2)*pow(RS,3) + 120*k2*pow(RS,4) + 120*pow(RS,5)))\
+ 2*a1*c1*d*(1 + pow(M_E,2*d*pow(RS,-1)))*(-(pow(M_E,k2*pow(RS,-1))*
(5*RS*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RS,2) + 60*pow(k1,2)*pow(RS,3) + 120*k1*pow(RS,4) + 120*pow(RS,5))
) + pow(M_E,k1*pow(RS,-1))*(5*RS*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RS,2) + 60*pow(k2,2)*pow(RS,3) + 
120*k2*pow(RS,4) + 120*pow(RS,5))) + d*pow(b1,2)*(1 + pow(M_E,2*d*pow(RS,-1)))*
(-(pow(M_E,k2*pow(RS,-1))*(5*RS*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RS,2) + 60*pow(k1,2)*pow(RS,3) + 
120*k1*pow(RS,4) + 120*pow(RS,5))) + 
pow(M_E,k1*pow(RS,-1))*(5*RS*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RS,2) + 60*pow(k2,2)*pow(RS,3) + 
120*k2*pow(RS,4) + 120*pow(RS,5))) + 2*a1*c1*(-1 + pow(M_E,2*d*pow(RS,-1)))*
(pow(M_E,k2*pow(RS,-1))*(6*RS*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RS,2) + 120*pow(k1,3)*pow(RS,3) + 
360*pow(k1,2)*pow(RS,4) + 720*k1*pow(RS,5) + 720*pow(RS,6)) - 
pow(M_E,k1*pow(RS,-1))*(6*RS*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RS,2) + 120*pow(k2,3)*pow(RS,3) + 
360*pow(k2,2)*pow(RS,4) + 720*k2*pow(RS,5) + 720*pow(RS,6))) + 
2*a1*b1*RS*(-1 + pow(M_E,2*d*pow(RS,-1)))*(pow(M_E,k2*pow(RS,-1))*
(6*RS*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RS,2) + 120*pow(k1,3)*pow(RS,3) + 360*pow(k1,2)*pow(RS,4) + 
720*k1*pow(RS,5) + 720*pow(RS,6)) - pow(M_E,k1*pow(RS,-1))*
(6*RS*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RS,2) + 120*pow(k2,3)*pow(RS,3) + 360*pow(k2,2)*pow(RS,4) + 
720*k2*pow(RS,5) + 720*pow(RS,6))) + pow(b1,2)*(-1 + pow(M_E,2*d*pow(RS,-1)))*
(pow(M_E,k2*pow(RS,-1))*(6*RS*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RS,2) + 120*pow(k1,3)*pow(RS,3) + 
360*pow(k1,2)*pow(RS,4) + 720*k1*pow(RS,5) + 720*pow(RS,6)) - 
pow(M_E,k1*pow(RS,-1))*(6*RS*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RS,2) + 120*pow(k2,3)*pow(RS,3) + 
360*pow(k2,2)*pow(RS,4) + 720*k2*pow(RS,5) + 720*pow(RS,6))) + 
2*a1*b1*d*(1 + pow(M_E,2*d*pow(RS,-1)))*(-(pow(M_E,k2*pow(RS,-1))*
(6*RS*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RS,2) + 120*pow(k1,3)*pow(RS,3) + 360*pow(k1,2)*pow(RS,4) + 
720*k1*pow(RS,5) + 720*pow(RS,6))) + 
pow(M_E,k1*pow(RS,-1))*(6*RS*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RS,2) + 120*pow(k2,3)*pow(RS,3) + 
360*pow(k2,2)*pow(RS,4) + 720*k2*pow(RS,5) + 720*pow(RS,6))) + 
2*a1*b1*(-1 + pow(M_E,2*d*pow(RS,-1)))*(pow(M_E,k2*pow(RS,-1))*
(7*RS*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RS,2) + 210*pow(k1,4)*pow(RS,3) + 840*pow(k1,3)*pow(RS,4) + 
2520*pow(k1,2)*pow(RS,5) + 5040*k1*pow(RS,6) + 5040*pow(RS,7)) - 
pow(M_E,k1*pow(RS,-1))*(7*RS*pow(k2,6) + pow(k2,7) + 42*pow(k2,5)*pow(RS,2) + 210*pow(k2,4)*pow(RS,3) + 
840*pow(k2,3)*pow(RS,4) + 2520*pow(k2,2)*pow(RS,5) + 5040*k2*pow(RS,6) + 5040*pow(RS,7))) + 
RS*pow(a1,2)*(-1 + pow(M_E,2*d*pow(RS,-1)))*(pow(M_E,k2*pow(RS,-1))*
(7*RS*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RS,2) + 210*pow(k1,4)*pow(RS,3) + 840*pow(k1,3)*pow(RS,4) + 
2520*pow(k1,2)*pow(RS,5) + 5040*k1*pow(RS,6) + 5040*pow(RS,7)) - 
pow(M_E,k1*pow(RS,-1))*(7*RS*pow(k2,6) + pow(k2,7) + 42*pow(k2,5)*pow(RS,2) + 210*pow(k2,4)*pow(RS,3) + 
840*pow(k2,3)*pow(RS,4) + 2520*pow(k2,2)*pow(RS,5) + 5040*k2*pow(RS,6) + 5040*pow(RS,7))) + 
d*pow(a1,2)*(1 + pow(M_E,2*d*pow(RS,-1)))*(-(pow(M_E,k2*pow(RS,-1))*
(7*RS*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RS,2) + 210*pow(k1,4)*pow(RS,3) + 840*pow(k1,3)*pow(RS,4) + 
2520*pow(k1,2)*pow(RS,5) + 5040*k1*pow(RS,6) + 5040*pow(RS,7))) + 
pow(M_E,k1*pow(RS,-1))*(7*RS*pow(k2,6) + pow(k2,7) + 42*pow(k2,5)*pow(RS,2) + 210*pow(k2,4)*pow(RS,3) + 
840*pow(k2,3)*pow(RS,4) + 2520*pow(k2,2)*pow(RS,5) + 5040*k2*pow(RS,6) + 5040*pow(RS,7))) + 
pow(a1,2)*(-1 + pow(M_E,2*d*pow(RS,-1)))*(pow(M_E,k2*pow(RS,-1))*
(8*RS*pow(k1,7) + pow(k1,8) + 56*pow(k1,6)*pow(RS,2) + 336*pow(k1,5)*pow(RS,3) + 1680*pow(k1,4)*pow(RS,4) + 
6720*pow(k1,3)*pow(RS,5) + 20160*pow(k1,2)*pow(RS,6) + 40320*k1*pow(RS,7) + 40320*pow(RS,8)) - 
pow(M_E,k1*pow(RS,-1))*(8*RS*pow(k2,7) + pow(k2,8) + 56*pow(k2,6)*pow(RS,2) + 336*pow(k2,5)*pow(RS,3) + 
1680*pow(k2,4)*pow(RS,4) + 6720*pow(k2,3)*pow(RS,5) + 20160*pow(k2,2)*pow(RS,6) + 40320*k2*pow(RS,7) + 
40320*pow(RS,8))));

    return res;
}




/* BM_SS 14 SZ */

//1
double BMSH_Integral_14_case_1( double ASP, double RSP, double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d )
{
double res = -(ASP*pow(3,0.5)*pow(M_E,-((d + k1 + k2)*pow(RSP,-1)))*pow(RSP,2)*
(pow(d,2)*(-4*b1*d2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,2) - 4*b1*d2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2) - 
b1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) - 5*b1*c2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) - 
5*a1*d2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + b1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3) - 
5*b1*c2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3) - 5*a1*d2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3) - 
b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) - a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) - 
6*b1*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) - 6*a1*c2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + 
b1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) + a1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) - 
6*b1*b2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) - 6*a1*c2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) - 
b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) - a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) - 
7*a2*b1*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) - 7*a1*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + 
b1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) + a1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) - 
7*a2*b1*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) - 7*a1*b2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) - 
a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) - a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) - 
8*a1*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) + a2*b1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,6) + 
a1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,6) - 8*a1*a2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,6) - 
a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,7) + a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,7) + 
4*b1*d2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,2) + 4*b1*d2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2) + 
b1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) + 5*b1*c2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) + 
5*a1*d2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - b1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3) + 
5*b1*c2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3) + 5*a1*d2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3) + 
b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) + a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) + 
6*b1*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) + 6*a1*c2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 
b1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) - a1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) + 
6*b1*b2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) + 6*a1*c2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) + 
b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) + a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) + 
7*a2*b1*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) + 7*a1*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - 
b1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) - a1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) + 
7*a2*b1*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) + 7*a1*b2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) + 
a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) + a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) + 
8*a1*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) - a2*b1*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,6) - 
a1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,6) + 8*a1*a2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,6) + 
a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,7) - a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,7) + 
8*b1*d2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 8*b1*d2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 
8*b1*d2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,2) - 8*b1*d2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,2) - 
15*b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 15*a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 
15*b1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 
15*a1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 
24*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 24*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 
24*b1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 
24*a1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 
35*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) - 35*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 
35*a2*b1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 
35*a1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) - 
48*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) + 
48*a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) + 
15*b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 15*a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 
15*b1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 
15*a1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 
24*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 24*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 
24*b1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 
24*a1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 
35*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) + 35*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) - 
35*a2*b1*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) - 
35*a1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) + 
48*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) - 
48*a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) + 8*b1*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) + 
30*b1*c2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) + 30*a1*d2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
8*b1*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 30*b1*c2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 
30*a1*d2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 8*b1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,3) - 
30*b1*c2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,3) - 30*a1*d2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,3) + 
8*b1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,3) + 30*b1*c2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,3) + 
30*a1*d2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,3) - 72*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) - 
72*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) - 
72*b1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) - 
72*a1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) - 
140*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) - 140*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) - 
140*a2*b1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) - 
140*a1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) - 
240*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) - 
240*a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) + 
72*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 72*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 
72*b1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 
72*a1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 
140*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) + 140*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) + 
140*a2*b1*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) + 
140*a1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) + 
240*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) + 
240*a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) + 30*b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 
30*a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 144*b1*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 
144*a1*c2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 30*b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 
30*a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 144*b1*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 
144*a1*c2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 30*b1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) + 
30*a1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) + 144*b1*b2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) + 
144*a1*c2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) - 30*b1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) - 
30*a1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) - 144*b1*b2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) - 
144*a1*c2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) - 420*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) - 
420*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 
420*a2*b1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 
420*a1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) - 
960*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) + 
960*a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) + 
420*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 420*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) - 
420*a2*b1*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) - 
420*a1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 
960*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) - 
960*a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) + 
d1*(d2*((k2 + 2*RSP)*pow(M_E,k1*pow(RSP,-1)) - (k1 + 2*RSP)*pow(M_E,k2*pow(RSP,-1)) + 
(k1 - 2*RSP)*pow(M_E,(2*k1 + k2)*pow(RSP,-1)) + (-k2 + 2*RSP)*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))) - 
4*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,2) - 4*b2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2) - 
b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) - 5*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + 
b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3) - 5*a2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3) - 
a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) + 
4*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,2) + 4*b2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2) + 
b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) + 5*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - 
b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3) + 5*a2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3) + 
a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) + 
8*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 8*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 
8*b2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,2) - 8*b2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,2) - 
15*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 15*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 
15*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 15*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 
c2*(pow(M_E,(2*k1 + k2)*pow(RSP,-1))*(-3*k1*RSP + pow(k1,2) + 3*pow(RSP,2)) - 
pow(M_E,k2*pow(RSP,-1))*(3*k1*RSP + pow(k1,2) + 3*pow(RSP,2)) - 
pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*(-3*k2*RSP + pow(k2,2) + 3*pow(RSP,2)) + 
pow(M_E,k1*pow(RSP,-1))*(3*k2*RSP + pow(k2,2) + 3*pow(RSP,2))) + 8*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) + 
30*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 8*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 
30*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 8*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,3) - 
30*a2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,3) + 8*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,3) + 
30*a2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,3) + 30*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 
30*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 30*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) - 
30*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4)) + 144*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 
144*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 840*a2*b1*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 
840*a1*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 144*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 
144*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 840*a2*b1*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 
840*a1*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 144*b1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,5) - 
144*a1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,5) - 840*a2*b1*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,5) - 
840*a1*b2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,5) + 144*b1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,5) + 
144*a1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,5) + 840*a2*b1*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,5) + 
840*a1*b2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,5) - 2880*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) - 
2880*a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) + 
2880*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) + 
2880*a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) + 
c1*(-5*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) - 5*b2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3) - 
b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) - 6*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + 
b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) - 6*a2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) - 
a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) + 
5*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) + 5*b2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3) + 
b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) + 6*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 
b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) + 6*a2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) + 
a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) - 
15*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 15*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 
24*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 24*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 
15*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 15*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 
24*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 24*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 
d2*(pow(M_E,(2*k1 + k2)*pow(RSP,-1))*(-3*k1*RSP + pow(k1,2) + 3*pow(RSP,2)) - 
pow(M_E,k2*pow(RSP,-1))*(3*k1*RSP + pow(k1,2) + 3*pow(RSP,2)) - 
pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*(-3*k2*RSP + pow(k2,2) + 3*pow(RSP,2)) + 
pow(M_E,k1*pow(RSP,-1))*(3*k2*RSP + pow(k2,2) + 3*pow(RSP,2))) + 30*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
30*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 30*b2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,3) + 
30*b2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,3) - 72*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) - 
72*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 72*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 
72*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 
c2*(pow(M_E,(2*k1 + k2)*pow(RSP,-1))*(-4*RSP*pow(k1,2) + pow(k1,3) + 8*k1*pow(RSP,2) - 8*pow(RSP,3)) - 
pow(M_E,k2*pow(RSP,-1))*(4*RSP*pow(k1,2) + pow(k1,3) + 8*k1*pow(RSP,2) + 8*pow(RSP,3)) + 
pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*(4*RSP*pow(k2,2) - pow(k2,3) - 8*k2*pow(RSP,2) + 8*pow(RSP,3)) + 
pow(M_E,k1*pow(RSP,-1))*(4*RSP*pow(k2,2) + pow(k2,3) + 8*k2*pow(RSP,2) + 8*pow(RSP,3))) + 
30*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 144*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 
30*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 144*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 
30*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) + 144*a2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) - 
30*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) - 144*a2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) + 
144*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 144*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 
144*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,5) + 144*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,5)) + 
840*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) + 840*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) + 
5760*a1*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) - 840*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) - 
840*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) - 5760*a1*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) + 
840*a2*b1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,6) + 840*a1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,6) + 
5760*a1*a2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,6) - 840*a2*b1*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,6) - 
840*a1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,6) - 5760*a1*a2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,6) + 
5760*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,7) - 5760*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,7) - 
5760*a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,7) + 5760*a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,7)) + 
d*(-7*b1*d2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + 7*b1*d2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3) - 
b1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) - 8*b1*c2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) - 
8*a1*d2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) - b1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) + 
8*b1*c2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) + 8*a1*d2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) - 
b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) - a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) - 
9*b1*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) - 9*a1*c2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) - 
b1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) - a1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) + 
9*b1*b2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) + 9*a1*c2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) - 
b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) - a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) - 
10*a2*b1*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) - 10*a1*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) - 
b1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,6) - a1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,6) + 
10*a2*b1*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,6) + 10*a1*b2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,6) - 
a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,7) - a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,7) - 
11*a1*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,7) - a2*b1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,7) - 
a1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,7) + 11*a1*a2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,7) - 
a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,8) - a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,8) + 
7*b1*d2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - 7*b1*d2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3) + 
b1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) + 8*b1*c2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) + 
8*a1*d2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) + b1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) - 
8*b1*c2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) - 8*a1*d2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) + 
b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) + a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) + 
9*b1*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) + 9*a1*c2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) + 
b1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) + a1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) - 
9*b1*b2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) - 9*a1*c2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) + 
b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) + a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) + 
10*a2*b1*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) + 10*a1*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) + 
b1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,6) + a1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,6) - 
10*a2*b1*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,6) - 10*a1*b2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,6) + 
a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,7) + a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,7) + 
11*a1*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,7) + a2*b1*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,7) + 
a1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,7) - 11*a1*a2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,7) + 
a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,8) + a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,8) - 
24*b1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 
24*b1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 
35*b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 35*a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 
35*b1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 
35*a1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 
48*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) - 48*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) - 
48*b1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) - 
48*a1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) - 
63*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) - 63*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) - 
63*a2*b1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) - 
63*a1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) - 
80*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,6)*pow(RSP,2) - 
80*a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,6)*pow(RSP,2) + 
24*b1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 
24*b1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 
35*b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 35*a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 
35*b1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 
35*a1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 
48*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) + 48*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) + 
48*b1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) + 
48*a1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) + 
63*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) + 63*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) + 
63*a2*b1*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) + 
63*a1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) + 
80*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,6)*pow(RSP,2) + 
80*a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,6)*pow(RSP,2) + 48*b1*d2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
48*b1*d2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 48*b1*d2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,3) - 
48*b1*d2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,3) - 105*b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) - 
105*a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
105*b1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
105*a1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) - 
192*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) - 192*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 
192*b1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 
192*a1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) - 
315*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) - 315*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) + 
315*a2*b1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) + 
315*a1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) - 
480*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5)*pow(RSP,3) + 
480*a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5)*pow(RSP,3) + 
105*b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 105*a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 
105*b1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 
105*a1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 
192*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) + 192*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) - 
192*b1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) - 
192*a1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) + 
315*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) + 315*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) - 
315*a2*b1*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) - 
315*a1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) + 
480*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5)*pow(RSP,3) - 
480*a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5)*pow(RSP,3) + 48*b1*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 
210*b1*c2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 210*a1*d2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 
48*b1*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 210*b1*c2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 
210*a1*d2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 48*b1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) - 
210*b1*c2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) - 210*a1*d2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) + 
48*b1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) + 210*b1*c2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) + 
210*a1*d2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) - 576*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) - 
576*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) - 
576*b1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) - 
576*a1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) - 
1260*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) - 1260*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) - 
1260*a2*b1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) - 
1260*a1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) - 
2400*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,4) - 
2400*a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4)*pow(RSP,4) + 
576*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 576*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 
576*b1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 
576*a1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 
1260*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) + 1260*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) + 
1260*a2*b1*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) + 
1260*a1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) + 
2400*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,4) + 
2400*a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4)*pow(RSP,4) + 210*b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 
210*a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 1152*b1*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 
1152*a1*c2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 210*b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 
210*a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 1152*b1*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 
1152*a1*c2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 210*b1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,5) + 
210*a1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,5) + 1152*b1*b2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,5) + 
1152*a1*c2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,5) - 210*b1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,5) - 
210*a1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,5) - 1152*b1*b2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,5) - 
1152*a1*c2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,5) - 3780*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) - 
3780*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) + 
3780*a2*b1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) + 
3780*a1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) - 
9600*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,5) + 
9600*a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,5) + 
3780*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) + 3780*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) - 
3780*a2*b1*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) - 
3780*a1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) + 
9600*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,5) - 
9600*a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,5) + 
d1*(-7*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + 7*b2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3) - 
b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) - 8*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) - 
b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) + 8*a2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) - 
a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) - a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) + 
7*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - 7*b2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3) + 
b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) + 8*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) + 
b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) - 8*a2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) + 
a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) + a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) - 
24*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 24*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 
35*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 35*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 
24*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 24*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 
35*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 35*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 
d2*(-(pow(M_E,(2*k1 + k2)*pow(RSP,-1))*(-5*k1*RSP + pow(k1,2) + 8*pow(RSP,2))) - 
pow(M_E,k2*pow(RSP,-1))*(5*k1*RSP + pow(k1,2) + 8*pow(RSP,2)) + 
pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*(-5*k2*RSP + pow(k2,2) + 8*pow(RSP,2)) + 
pow(M_E,k1*pow(RSP,-1))*(5*k2*RSP + pow(k2,2) + 8*pow(RSP,2))) + 48*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
48*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 48*b2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,3) - 
48*b2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,3) - 105*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
105*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
105*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 
105*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 
c2*(-(pow(M_E,(2*k1 + k2)*pow(RSP,-1))*(-6*RSP*pow(k1,2) + pow(k1,3) + 15*k1*pow(RSP,2) - 15*pow(RSP,3))) + 
pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*(-6*RSP*pow(k2,2) + pow(k2,3) + 15*k2*pow(RSP,2) - 15*pow(RSP,3)) - 
pow(M_E,k2*pow(RSP,-1))*(6*RSP*pow(k1,2) + pow(k1,3) + 15*k1*pow(RSP,2) + 15*pow(RSP,3)) + 
pow(M_E,k1*pow(RSP,-1))*(6*RSP*pow(k2,2) + pow(k2,3) + 15*k2*pow(RSP,2) + 15*pow(RSP,3))) + 
48*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 210*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 
48*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 210*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 
48*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) - 210*a2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) + 
48*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) + 210*a2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) + 
210*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 210*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 
210*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,5) - 210*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,5)) + 
1152*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) + 1152*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) + 
7560*a2*b1*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) + 7560*a1*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) - 
1152*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) - 1152*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) - 
7560*a2*b1*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) - 7560*a1*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) - 
1152*b1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,6) - 1152*a1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,6) - 
7560*a2*b1*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,6) - 
7560*a1*b2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,6) + 1152*b1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,6) + 
1152*a1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,6) + 7560*a2*b1*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,6) + 
7560*a1*b2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,6) - 
28800*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,6) - 
28800*a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,6) + 
28800*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,6) + 
28800*a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,6) + 
c1*(-8*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + 8*b2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) - 
b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) - 9*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) - 
b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) + 9*a2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) - 
a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) - a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,6) + 
8*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 8*b2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) + 
b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) + 9*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) + 
b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) - 9*a2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) + 
a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) + a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,6) - 
35*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 35*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 
48*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) - 48*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 
35*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 35*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 
48*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) + 48*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) - 
105*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
105*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) - 
192*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 
192*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 
105*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 
105*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 
192*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) - 
192*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) + 
d2*(-(pow(M_E,(2*k1 + k2)*pow(RSP,-1))*(-6*RSP*pow(k1,2) + pow(k1,3) + 15*k1*pow(RSP,2) - 15*pow(RSP,3))) + 
pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*(-6*RSP*pow(k2,2) + pow(k2,3) + 15*k2*pow(RSP,2) - 15*pow(RSP,3)) - 
pow(M_E,k2*pow(RSP,-1))*(6*RSP*pow(k1,2) + pow(k1,3) + 15*k1*pow(RSP,2) + 15*pow(RSP,3)) + 
pow(M_E,k1*pow(RSP,-1))*(6*RSP*pow(k2,2) + pow(k2,3) + 15*k2*pow(RSP,2) + 15*pow(RSP,3))) + 
210*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 210*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 
210*b2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) + 210*b2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) - 
576*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) - 
576*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 
576*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 
576*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 
c2*(-(pow(M_E,(2*k1 + k2)*pow(RSP,-1))*
(-7*RSP*pow(k1,3) + pow(k1,4) + 24*pow(k1,2)*pow(RSP,2) - 48*k1*pow(RSP,3) + 48*pow(RSP,4))) - 
pow(M_E,k2*pow(RSP,-1))*(7*RSP*pow(k1,3) + pow(k1,4) + 24*pow(k1,2)*pow(RSP,2) + 48*k1*pow(RSP,3) + 
48*pow(RSP,4)) + pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*
(-7*RSP*pow(k2,3) + pow(k2,4) + 24*pow(k2,2)*pow(RSP,2) - 48*k2*pow(RSP,3) + 48*pow(RSP,4)) + 
pow(M_E,k1*pow(RSP,-1))*(7*RSP*pow(k2,3) + pow(k2,4) + 24*pow(k2,2)*pow(RSP,2) + 48*k2*pow(RSP,3) + 
48*pow(RSP,4))) + 210*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 
1152*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 210*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 
1152*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 210*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,5) + 
1152*a2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,5) - 210*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,5) - 
1152*a2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,5) + 1152*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) - 
1152*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) - 1152*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,6) + 
1152*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,6)) + 7560*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(RSP,7) + 
7560*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,7) + 57600*a1*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,7) - 
7560*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,7) - 7560*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,7) - 
57600*a1*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,7) + 7560*a2*b1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,7) + 
7560*a1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,7) + 57600*a1*a2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,7) - 
7560*a2*b1*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,7) - 7560*a1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,7) - 
57600*a1*a2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,7) + 57600*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,8) - 
57600*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,8) - 57600*a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,8) + 
57600*a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,8)) + 
RSP*(-7*b1*d2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + 7*b1*d2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3) - 
b1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) - 8*b1*c2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) - 
8*a1*d2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) - b1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) + 
8*b1*c2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) + 8*a1*d2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) - 
b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) - a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) - 
9*b1*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) - 9*a1*c2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) - 
b1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) - a1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) + 
9*b1*b2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) + 9*a1*c2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) - 
b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) - a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) - 
10*a2*b1*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) - 10*a1*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) - 
b1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,6) - a1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,6) + 
10*a2*b1*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,6) + 10*a1*b2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,6) - 
a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,7) - a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,7) - 
11*a1*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,7) - a2*b1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,7) - 
a1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,7) + 11*a1*a2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,7) - 
a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,8) - a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,8) + 
7*b1*d2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - 7*b1*d2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3) + 
b1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) + 8*b1*c2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) + 
8*a1*d2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) + b1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) - 
8*b1*c2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) - 8*a1*d2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) + 
b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) + a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) + 
9*b1*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) + 9*a1*c2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) + 
b1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) + a1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) - 
9*b1*b2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) - 9*a1*c2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) + 
b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) + a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) + 
10*a2*b1*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) + 10*a1*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) + 
b1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,6) + a1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,6) - 
10*a2*b1*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,6) - 10*a1*b2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,6) + 
a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,7) + a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,7) + 
11*a1*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,7) + a2*b1*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,7) + 
a1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,7) - 11*a1*a2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,7) + 
a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,8) + a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,8) - 
24*b1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 
24*b1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 
35*b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 35*a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 
35*b1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 
35*a1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 
48*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) - 48*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) - 
48*b1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) - 
48*a1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) - 
63*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) - 63*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) - 
63*a2*b1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) - 
63*a1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) - 
80*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,6)*pow(RSP,2) - 
80*a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,6)*pow(RSP,2) + 
24*b1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 
24*b1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 
35*b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 35*a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 
35*b1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 
35*a1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 
48*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) + 48*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) + 
48*b1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) + 
48*a1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) + 
63*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) + 63*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) + 
63*a2*b1*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) + 
63*a1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) + 
80*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,6)*pow(RSP,2) + 
80*a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,6)*pow(RSP,2) + 48*b1*d2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
48*b1*d2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 48*b1*d2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,3) - 
48*b1*d2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,3) - 105*b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) - 
105*a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
105*b1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
105*a1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) - 
192*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) - 192*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 
192*b1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 
192*a1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) - 
315*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) - 315*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) + 
315*a2*b1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) + 
315*a1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) - 
480*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5)*pow(RSP,3) + 
480*a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5)*pow(RSP,3) + 
105*b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 105*a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 
105*b1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 
105*a1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 
192*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) + 192*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) - 
192*b1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) - 
192*a1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) + 
315*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) + 315*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) - 
315*a2*b1*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) - 
315*a1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) + 
480*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5)*pow(RSP,3) - 
480*a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5)*pow(RSP,3) + 48*b1*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 
210*b1*c2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 210*a1*d2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 
48*b1*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 210*b1*c2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 
210*a1*d2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 48*b1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) - 
210*b1*c2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) - 210*a1*d2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) + 
48*b1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) + 210*b1*c2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) + 
210*a1*d2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) - 576*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) - 
576*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) - 
576*b1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) - 
576*a1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) - 
1260*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) - 1260*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) - 
1260*a2*b1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) - 
1260*a1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) - 
2400*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,4) - 
2400*a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4)*pow(RSP,4) + 
576*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 576*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 
576*b1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 
576*a1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 
1260*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) + 1260*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) + 
1260*a2*b1*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) + 
1260*a1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) + 
2400*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,4) + 
2400*a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4)*pow(RSP,4) + 210*b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 
210*a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 1152*b1*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 
1152*a1*c2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 210*b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 
210*a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 1152*b1*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 
1152*a1*c2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 210*b1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,5) + 
210*a1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,5) + 1152*b1*b2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,5) + 
1152*a1*c2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,5) - 210*b1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,5) - 
210*a1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,5) - 1152*b1*b2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,5) - 
1152*a1*c2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,5) - 3780*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) - 
3780*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) + 
3780*a2*b1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) + 
3780*a1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) - 
9600*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,5) + 
9600*a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,5) + 
3780*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) + 3780*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) - 
3780*a2*b1*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) - 
3780*a1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) + 
9600*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,5) - 
9600*a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,5) + 
d1*(-7*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + 7*b2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3) - 
b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) - 8*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) - 
b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) + 8*a2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) - 
a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) - a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) + 
7*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - 7*b2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3) + 
b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) + 8*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) + 
b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) - 8*a2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) + 
a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) + a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) - 
24*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 24*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 
35*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 35*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 
24*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 24*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 
35*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 35*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 
d2*(-(pow(M_E,(2*k1 + k2)*pow(RSP,-1))*(-5*k1*RSP + pow(k1,2) + 8*pow(RSP,2))) - 
pow(M_E,k2*pow(RSP,-1))*(5*k1*RSP + pow(k1,2) + 8*pow(RSP,2)) + 
pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*(-5*k2*RSP + pow(k2,2) + 8*pow(RSP,2)) + 
pow(M_E,k1*pow(RSP,-1))*(5*k2*RSP + pow(k2,2) + 8*pow(RSP,2))) + 48*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
48*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 48*b2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,3) - 
48*b2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,3) - 105*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
105*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
105*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 
105*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 
c2*(-(pow(M_E,(2*k1 + k2)*pow(RSP,-1))*(-6*RSP*pow(k1,2) + pow(k1,3) + 15*k1*pow(RSP,2) - 15*pow(RSP,3))) + 
pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*(-6*RSP*pow(k2,2) + pow(k2,3) + 15*k2*pow(RSP,2) - 15*pow(RSP,3)) - 
pow(M_E,k2*pow(RSP,-1))*(6*RSP*pow(k1,2) + pow(k1,3) + 15*k1*pow(RSP,2) + 15*pow(RSP,3)) + 
pow(M_E,k1*pow(RSP,-1))*(6*RSP*pow(k2,2) + pow(k2,3) + 15*k2*pow(RSP,2) + 15*pow(RSP,3))) + 
48*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 210*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 
48*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 210*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 
48*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) - 210*a2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) + 
48*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) + 210*a2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) + 
210*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 210*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 
210*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,5) - 210*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,5)) + 
1152*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) + 1152*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) + 
7560*a2*b1*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) + 7560*a1*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) - 
1152*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) - 1152*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) - 
7560*a2*b1*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) - 7560*a1*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) - 
1152*b1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,6) - 1152*a1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,6) - 
7560*a2*b1*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,6) - 
7560*a1*b2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,6) + 1152*b1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,6) + 
1152*a1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,6) + 7560*a2*b1*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,6) + 
7560*a1*b2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,6) - 
28800*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,6) - 
28800*a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,6) + 
28800*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,6) + 
28800*a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,6) + 
c1*(-8*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + 8*b2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) - 
b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) - 9*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) - 
b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) + 9*a2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) - 
a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) - a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,6) + 
8*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 8*b2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) + 
b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) + 9*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) + 
b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) - 9*a2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) + 
a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) + a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,6) - 
35*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 35*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 
48*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) - 48*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 
35*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 35*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 
48*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) + 48*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) - 
105*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
105*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) - 
192*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 
192*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 
105*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 
105*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 
192*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) - 
192*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) + 
d2*(-(pow(M_E,(2*k1 + k2)*pow(RSP,-1))*(-6*RSP*pow(k1,2) + pow(k1,3) + 15*k1*pow(RSP,2) - 15*pow(RSP,3))) + 
pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*(-6*RSP*pow(k2,2) + pow(k2,3) + 15*k2*pow(RSP,2) - 15*pow(RSP,3)) - 
pow(M_E,k2*pow(RSP,-1))*(6*RSP*pow(k1,2) + pow(k1,3) + 15*k1*pow(RSP,2) + 15*pow(RSP,3)) + 
pow(M_E,k1*pow(RSP,-1))*(6*RSP*pow(k2,2) + pow(k2,3) + 15*k2*pow(RSP,2) + 15*pow(RSP,3))) + 
210*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 210*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 
210*b2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) + 210*b2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) - 
576*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) - 
576*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 
576*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 
576*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 
c2*(-(pow(M_E,(2*k1 + k2)*pow(RSP,-1))*
(-7*RSP*pow(k1,3) + pow(k1,4) + 24*pow(k1,2)*pow(RSP,2) - 48*k1*pow(RSP,3) + 48*pow(RSP,4))) - 
pow(M_E,k2*pow(RSP,-1))*(7*RSP*pow(k1,3) + pow(k1,4) + 24*pow(k1,2)*pow(RSP,2) + 48*k1*pow(RSP,3) + 
48*pow(RSP,4)) + pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*
(-7*RSP*pow(k2,3) + pow(k2,4) + 24*pow(k2,2)*pow(RSP,2) - 48*k2*pow(RSP,3) + 48*pow(RSP,4)) + 
pow(M_E,k1*pow(RSP,-1))*(7*RSP*pow(k2,3) + pow(k2,4) + 24*pow(k2,2)*pow(RSP,2) + 48*k2*pow(RSP,3) + 
48*pow(RSP,4))) + 210*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 
1152*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 210*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 
1152*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 210*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,5) + 
1152*a2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,5) - 210*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,5) - 
1152*a2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,5) + 1152*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) - 
1152*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) - 1152*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,6) + 
1152*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,6)) + 7560*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(RSP,7) + 
7560*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,7) + 57600*a1*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,7) - 
7560*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,7) - 7560*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,7) - 
57600*a1*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,7) + 7560*a2*b1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,7) + 
7560*a1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,7) + 57600*a1*a2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,7) - 
7560*a2*b1*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,7) - 7560*a1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,7) - 
57600*a1*a2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,7) + 57600*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,8) - 
57600*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,8) - 57600*a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,8) + 
57600*a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,8)))*pow(2*pow(d,2),-1));

    return res;
}


//  2
double BMSH_Integral_14_case_2_sub_1( double ASP, double RSP, double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d )
{
double res = -(ASP*pow(3,0.5)*pow(M_E,-((2*d + k2)*pow(RSP,-1)))*pow(RSP,2)*
(5*d*d1*d2*k2*RSP*pow(M_E,d*pow(RSP,-1)) + d1*d2*k2*pow(d,2)*pow(M_E,d*pow(RSP,-1)) + 
2*d1*d2*RSP*pow(d,2)*pow(M_E,d*pow(RSP,-1)) + 3*c2*d1*k2*RSP*pow(d,2)*pow(M_E,d*pow(RSP,-1)) + 
5*d*d1*d2*k2*RSP*pow(M_E,3*d*pow(RSP,-1)) - d1*d2*k2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1)) - 
2*d1*d2*RSP*pow(d,2)*pow(M_E,3*d*pow(RSP,-1)) - 3*c2*d1*k2*RSP*pow(d,2)*pow(M_E,3*d*pow(RSP,-1)) - 
8*d1*d2*RSP*pow(d,2)*pow(M_E,k2*pow(RSP,-1)) - 2*d1*d2*pow(d,3)*pow(M_E,k2*pow(RSP,-1)) - 
10*c2*d1*RSP*pow(d,3)*pow(M_E,k2*pow(RSP,-1)) - 2*c2*d1*pow(d,4)*pow(M_E,k2*pow(RSP,-1)) - 
12*b2*d1*RSP*pow(d,4)*pow(M_E,k2*pow(RSP,-1)) - 2*b2*d1*pow(d,5)*pow(M_E,k2*pow(RSP,-1)) - 
14*a2*d1*RSP*pow(d,5)*pow(M_E,k2*pow(RSP,-1)) - 14*a1*d2*RSP*pow(d,5)*pow(M_E,k2*pow(RSP,-1)) - 
2*a2*d1*pow(d,6)*pow(M_E,k2*pow(RSP,-1)) - 2*a1*d2*pow(d,6)*pow(M_E,k2*pow(RSP,-1)) - 
16*a1*c2*RSP*pow(d,6)*pow(M_E,k2*pow(RSP,-1)) - 2*a1*c2*pow(d,7)*pow(M_E,k2*pow(RSP,-1)) - 
18*a1*b2*RSP*pow(d,7)*pow(M_E,k2*pow(RSP,-1)) - 2*a1*b2*pow(d,8)*pow(M_E,k2*pow(RSP,-1)) - 
20*a1*a2*RSP*pow(d,8)*pow(M_E,k2*pow(RSP,-1)) - 2*a1*a2*pow(d,9)*pow(M_E,k2*pow(RSP,-1)) - 
2*d1*d2*RSP*pow(d,2)*pow(M_E,(2*d + k2)*pow(RSP,-1)) - 2*c2*d1*RSP*pow(d,3)*pow(M_E,(2*d + k2)*pow(RSP,-1)) - 
2*b2*d1*RSP*pow(d,4)*pow(M_E,(2*d + k2)*pow(RSP,-1)) - 2*a2*d1*RSP*pow(d,5)*pow(M_E,(2*d + k2)*pow(RSP,-1)) - 
2*a1*d2*RSP*pow(d,5)*pow(M_E,(2*d + k2)*pow(RSP,-1)) - 2*a1*c2*RSP*pow(d,6)*pow(M_E,(2*d + k2)*pow(RSP,-1)) - 
2*a1*b2*RSP*pow(d,7)*pow(M_E,(2*d + k2)*pow(RSP,-1)) - 2*a1*a2*RSP*pow(d,8)*pow(M_E,(2*d + k2)*pow(RSP,-1)) + 
d*d1*d2*pow(M_E,d*pow(RSP,-1))*pow(k2,2) + 6*c2*d*d1*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,2) + 
d1*d2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,2) + c2*d1*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,2) + 
4*b2*d1*RSP*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,2) + d*d1*d2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2) + 
6*c2*d*d1*RSP*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2) - d1*d2*RSP*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2) - 
c2*d1*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2) - 4*b2*d1*RSP*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2) + 
c2*d*d1*pow(M_E,d*pow(RSP,-1))*pow(k2,3) + c2*d1*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,3) + 
7*b2*d*d1*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,3) + b2*d1*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,3) + 
5*a2*d1*RSP*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,3) + 5*a1*d2*RSP*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,3) + 
c2*d*d1*pow(M_E,3*d*pow(RSP,-1))*pow(k2,3) - c2*d1*RSP*pow(M_E,3*d*pow(RSP,-1))*pow(k2,3) + 
7*b2*d*d1*RSP*pow(M_E,3*d*pow(RSP,-1))*pow(k2,3) - b2*d1*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,3) - 
5*a2*d1*RSP*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,3) - 5*a1*d2*RSP*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,3) + 
b2*d*d1*pow(M_E,d*pow(RSP,-1))*pow(k2,4) + b2*d1*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,4) + 
8*a2*d*d1*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,4) + 8*a1*d*d2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,4) + 
a2*d1*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,4) + a1*d2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,4) + 
6*a1*c2*RSP*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,4) + b2*d*d1*pow(M_E,3*d*pow(RSP,-1))*pow(k2,4) - 
b2*d1*RSP*pow(M_E,3*d*pow(RSP,-1))*pow(k2,4) + 8*a2*d*d1*RSP*pow(M_E,3*d*pow(RSP,-1))*pow(k2,4) + 
8*a1*d*d2*RSP*pow(M_E,3*d*pow(RSP,-1))*pow(k2,4) - a2*d1*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,4) - 
a1*d2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,4) - 6*a1*c2*RSP*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,4) + 
a2*d*d1*pow(M_E,d*pow(RSP,-1))*pow(k2,5) + a1*d*d2*pow(M_E,d*pow(RSP,-1))*pow(k2,5) + 
9*a1*c2*d*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,5) + a2*d1*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,5) + 
a1*d2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,5) + a1*c2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,5) + 
7*a1*b2*RSP*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,5) + a2*d*d1*pow(M_E,3*d*pow(RSP,-1))*pow(k2,5) + 
a1*d*d2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,5) + 9*a1*c2*d*RSP*pow(M_E,3*d*pow(RSP,-1))*pow(k2,5) - 
a2*d1*RSP*pow(M_E,3*d*pow(RSP,-1))*pow(k2,5) - a1*d2*RSP*pow(M_E,3*d*pow(RSP,-1))*pow(k2,5) - 
a1*c2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,5) - 7*a1*b2*RSP*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,5) + 
a1*c2*d*pow(M_E,d*pow(RSP,-1))*pow(k2,6) + a1*c2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,6) + 
10*a1*b2*d*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,6) + a1*b2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,6) + 
8*a1*a2*RSP*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,6) + a1*c2*d*pow(M_E,3*d*pow(RSP,-1))*pow(k2,6) - 
a1*c2*RSP*pow(M_E,3*d*pow(RSP,-1))*pow(k2,6) + 10*a1*b2*d*RSP*pow(M_E,3*d*pow(RSP,-1))*pow(k2,6) - 
a1*b2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,6) - 8*a1*a2*RSP*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,6) + 
a1*b2*d*pow(M_E,d*pow(RSP,-1))*pow(k2,7) + a1*b2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,7) + 
11*a1*a2*d*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,7) + a1*a2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,7) + 
a1*b2*d*pow(M_E,3*d*pow(RSP,-1))*pow(k2,7) - a1*b2*RSP*pow(M_E,3*d*pow(RSP,-1))*pow(k2,7) + 
11*a1*a2*d*RSP*pow(M_E,3*d*pow(RSP,-1))*pow(k2,7) - a1*a2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,7) + 
a1*a2*d*pow(M_E,d*pow(RSP,-1))*pow(k2,8) + a1*a2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,8) + 
a1*a2*d*pow(M_E,3*d*pow(RSP,-1))*pow(k2,8) - a1*a2*RSP*pow(M_E,3*d*pow(RSP,-1))*pow(k2,8) + 
8*d*d1*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) + 15*c2*d*d1*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) + 
5*d1*d2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) + 3*c2*d1*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) + 
8*b2*d1*k2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) + 8*d*d1*d2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,2) + 
15*c2*d*d1*k2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,2) - 5*d1*d2*k2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,2) - 
3*c2*d1*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,2) - 8*b2*d1*k2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,2) - 
13*d*d1*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) - 24*c2*d1*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) - 
39*b2*d1*pow(d,3)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) - 58*a2*d1*pow(d,4)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) - 
58*a1*d2*pow(d,4)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) - 81*a1*c2*pow(d,5)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) - 
108*a1*b2*pow(d,6)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) - 139*a1*a2*pow(d,7)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) - 
3*d*d1*d2*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,2) - 6*c2*d1*pow(d,2)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,2) - 
9*b2*d1*pow(d,3)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,2) - 
12*a2*d1*pow(d,4)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,2) - 
12*a1*d2*pow(d,4)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,2) - 
15*a1*c2*pow(d,5)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,2) - 
18*a1*b2*pow(d,6)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,2) - 
21*a1*a2*pow(d,7)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,2) + 6*c2*d1*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 
24*b2*d*d1*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 15*a2*d1*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 
15*a1*d2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 6*c2*d1*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 
24*b2*d*d1*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 
15*a2*d1*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 
15*a1*d2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 7*b2*d1*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 
35*a2*d*d1*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 35*a1*d*d2*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 
24*a1*c2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 7*b2*d1*pow(M_E,3*d*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 
35*a2*d*d1*pow(M_E,3*d*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 35*a1*d*d2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 
24*a1*c2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 
48*a1*c2*d*pow(M_E,d*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) + 8*a2*d1*pow(M_E,d*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) + 
8*a1*d2*pow(M_E,d*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) + 35*a1*b2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) + 
48*a1*c2*d*pow(M_E,3*d*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) - 8*a2*d1*pow(M_E,3*d*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) - 
8*a1*d2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) - 
35*a1*b2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) + 9*a1*c2*pow(M_E,d*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) + 
63*a1*b2*d*pow(M_E,d*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) + 48*a1*a2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) - 
9*a1*c2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) + 63*a1*b2*d*pow(M_E,3*d*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) - 
48*a1*a2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) + 10*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k2,6)*pow(RSP,2) + 
80*a1*a2*d*pow(M_E,d*pow(RSP,-1))*pow(k2,6)*pow(RSP,2) - 10*a1*b2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,6)*pow(RSP,2) + 
80*a1*a2*d*pow(M_E,3*d*pow(RSP,-1))*pow(k2,6)*pow(RSP,2) + 11*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,7)*pow(RSP,2) - 
11*a1*a2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,7)*pow(RSP,2) + 15*c2*d*d1*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 
8*d1*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 15*c2*d1*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 
48*b2*d*d1*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 8*b2*d1*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 
30*a2*d1*k2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 30*a1*d2*k2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 
15*c2*d*d1*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,3) - 8*d1*d2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,3) - 
15*c2*d1*k2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,3) + 48*b2*d*d1*k2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,3) - 
8*b2*d1*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,3) - 30*a2*d1*k2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,3) - 
30*a1*d2*k2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,3) - 30*c2*d*d1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 
8*d1*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 80*b2*d1*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 
170*a2*d1*pow(d,3)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 170*a1*d2*pow(d,3)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 
312*a1*c2*pow(d,4)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 518*a1*b2*pow(d,5)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 
800*a1*a2*pow(d,6)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 8*d1*d2*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,3) - 
16*b2*d1*pow(d,2)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,3) - 
40*a2*d1*pow(d,3)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,3) - 
40*a1*d2*pow(d,3)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,3) - 
72*a1*c2*pow(d,4)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,3) - 
112*a1*b2*pow(d,5)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,3) - 
160*a1*a2*pow(d,6)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,3) + 24*b2*d1*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 
105*a2*d*d1*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 105*a1*d*d2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 
72*a1*c2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 24*b2*d1*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 
105*a2*d*d1*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 105*a1*d*d2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 
72*a1*c2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 
192*a1*c2*d*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) + 35*a2*d1*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) + 
35*a1*d2*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) + 140*a1*b2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) + 
192*a1*c2*d*pow(M_E,3*d*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) - 35*a2*d1*pow(M_E,3*d*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) - 
35*a1*d2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) - 
140*a1*b2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) + 
48*a1*c2*pow(M_E,d*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) + 315*a1*b2*d*pow(M_E,d*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) + 
240*a1*a2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) - 
48*a1*c2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) + 315*a1*b2*d*pow(M_E,3*d*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) - 
240*a1*a2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) + 
63*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k2,5)*pow(RSP,3) + 480*a1*a2*d*pow(M_E,d*pow(RSP,-1))*pow(k2,5)*pow(RSP,3) - 
63*a1*b2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,5)*pow(RSP,3) + 480*a1*a2*d*pow(M_E,3*d*pow(RSP,-1))*pow(k2,5)*pow(RSP,3) + 
80*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,6)*pow(RSP,3) - 80*a1*a2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,6)*pow(RSP,3) + 
15*c2*d1*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 48*b2*d*d1*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 
48*b2*d1*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 210*a2*d*d1*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 
210*a1*d*d2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 30*a2*d1*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 
30*a1*d2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 144*a1*c2*k2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) - 
15*c2*d1*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,4) + 48*b2*d*d1*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,4) - 
48*b2*d1*k2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,4) + 210*a2*d*d1*k2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,4) + 
210*a1*d*d2*k2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,4) - 30*a2*d1*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,4) - 
30*a1*d2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,4) - 144*a1*c2*k2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,4) - 
15*c2*d1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 96*b2*d*d1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 
345*a2*d1*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 345*a1*d2*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 
912*a1*c2*pow(d,3)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 1995*a1*b2*pow(d,4)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 
3840*a1*a2*pow(d,5)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 15*c2*d1*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,4) - 
75*a2*d1*pow(d,2)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,4) - 
75*a1*d2*pow(d,2)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,4) - 
240*a1*c2*pow(d,3)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,4) - 
525*a1*b2*pow(d,4)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,4) - 
960*a1*a2*pow(d,5)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,4) + 
576*a1*c2*d*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 105*a2*d1*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 
105*a1*d2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 420*a1*b2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 
576*a1*c2*d*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) - 105*a2*d1*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) - 
105*a1*d2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) - 
420*a1*b2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 
192*a1*c2*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) + 1260*a1*b2*d*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) + 
960*a1*a2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) - 
192*a1*c2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) + 1260*a1*b2*d*pow(M_E,3*d*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) - 
960*a1*a2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) + 
315*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k2,4)*pow(RSP,4) + 2400*a1*a2*d*pow(M_E,d*pow(RSP,-1))*pow(k2,4)*pow(RSP,4) - 
315*a1*b2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,4)*pow(RSP,4) + 2400*a1*a2*d*pow(M_E,3*d*pow(RSP,-1))*pow(k2,4)*pow(RSP,4) + 
480*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,5)*pow(RSP,4) - 480*a1*a2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,5)*pow(RSP,4) + 
48*b2*d1*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 210*a2*d*d1*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 
210*a1*d*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 1152*a1*c2*d*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 
210*a2*d1*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 210*a1*d2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 
144*a1*c2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 840*a1*b2*k2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) - 
48*b2*d1*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,5) + 210*a2*d*d1*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,5) + 
210*a1*d*d2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,5) + 1152*a1*c2*d*k2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,5) - 
210*a2*d1*k2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,5) - 210*a1*d2*k2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,5) - 
144*a1*c2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,5) - 840*a1*b2*k2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,5) - 
48*b2*d1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 420*a2*d*d1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 
420*a1*d*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 1872*a1*c2*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 
5880*a1*b2*pow(d,3)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 14880*a1*a2*pow(d,4)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 
48*b2*d1*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,5) - 432*a1*c2*pow(d,2)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,5) - 
1680*a1*b2*pow(d,3)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,5) - 
4320*a1*a2*pow(d,4)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,5) + 576*a1*c2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) + 
3780*a1*b2*d*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) + 
2880*a1*a2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) - 
576*a1*c2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) + 3780*a1*b2*d*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) - 
2880*a1*a2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) + 
1260*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,5) + 9600*a1*a2*d*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,5) - 
1260*a1*b2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,3)*pow(RSP,5) + 9600*a1*a2*d*pow(M_E,3*d*pow(RSP,-1))*pow(k2,3)*pow(RSP,5) + 
2400*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,4)*pow(RSP,5) - 2400*a1*a2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,4)*pow(RSP,5) + 
1152*a1*c2*d*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) + 210*a2*d1*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) + 
210*a1*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) + 1152*a1*c2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) + 
7560*a1*b2*d*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) + 840*a1*b2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) + 
5760*a1*a2*k2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) + 1152*a1*c2*d*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,6) - 
210*a2*d1*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,6) - 210*a1*d2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,6) - 
1152*a1*c2*k2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,6) + 7560*a1*b2*d*k2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,6) - 
840*a1*b2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,6) - 5760*a1*a2*k2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,6) - 
2304*a1*c2*d*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) - 210*a2*d1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) - 
210*a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) - 12180*a1*b2*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) - 
44160*a1*a2*pow(d,3)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) + 210*a2*d1*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,6) + 
210*a1*d2*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,6) - 2940*a1*b2*pow(d,2)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,6) - 
13440*a1*a2*pow(d,3)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,6) + 
3780*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,6) + 28800*a1*a2*d*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,6) - 
3780*a1*b2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2)*pow(RSP,6) + 28800*a1*a2*d*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2)*pow(RSP,6) + 
9600*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,6) - 9600*a1*a2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,3)*pow(RSP,6) + 
1152*a1*c2*pow(M_E,d*pow(RSP,-1))*pow(RSP,7) + 7560*a1*b2*d*pow(M_E,d*pow(RSP,-1))*pow(RSP,7) + 
7560*a1*b2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,7) + 57600*a1*a2*d*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,7) + 
5760*a1*a2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,7) - 1152*a1*c2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,7) + 
7560*a1*b2*d*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,7) - 7560*a1*b2*k2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,7) + 
57600*a1*a2*d*k2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,7) - 5760*a1*a2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,7) - 
1152*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,7) - 15120*a1*b2*d*pow(M_E,k2*pow(RSP,-1))*pow(RSP,7) - 
92160*a1*a2*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,7) + 1152*a1*c2*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,7) - 
23040*a1*a2*pow(d,2)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,7) + 
28800*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,7) - 28800*a1*a2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2)*pow(RSP,7) + 
c1*(3*d2*k2*RSP*pow(d,2)*pow(M_E,d*pow(RSP,-1)) - 3*d2*k2*RSP*pow(d,2)*pow(M_E,3*d*pow(RSP,-1)) - 
10*d2*RSP*pow(d,3)*pow(M_E,k2*pow(RSP,-1)) - 2*d2*pow(d,4)*pow(M_E,k2*pow(RSP,-1)) - 
16*a2*RSP*pow(d,6)*pow(M_E,k2*pow(RSP,-1)) - 2*a2*pow(d,7)*pow(M_E,k2*pow(RSP,-1)) - 
2*d2*RSP*pow(d,3)*pow(M_E,(2*d + k2)*pow(RSP,-1)) - 2*a2*RSP*pow(d,6)*pow(M_E,(2*d + k2)*pow(RSP,-1)) + 
6*d*d2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,2) + d2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,2) + 
6*d*d2*RSP*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2) - d2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2) + 
d*d2*pow(M_E,d*pow(RSP,-1))*pow(k2,3) + d2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,3) + 
d*d2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,3) - d2*RSP*pow(M_E,3*d*pow(RSP,-1))*pow(k2,3) + 
6*a2*RSP*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,4) - 6*a2*RSP*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,4) + 
9*a2*d*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,5) + a2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,5) + 
9*a2*d*RSP*pow(M_E,3*d*pow(RSP,-1))*pow(k2,5) - a2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,5) + 
a2*d*pow(M_E,d*pow(RSP,-1))*pow(k2,6) + a2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,6) + 
a2*d*pow(M_E,3*d*pow(RSP,-1))*pow(k2,6) - a2*RSP*pow(M_E,3*d*pow(RSP,-1))*pow(k2,6) + 
15*d*d2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) + 3*d2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) + 
15*d*d2*k2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,2) - 3*d2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,2) - 
24*d2*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) - 81*a2*pow(d,5)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) - 
6*d2*pow(d,2)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,2) - 
15*a2*pow(d,5)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,2) + 6*d2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 
6*d2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 24*a2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 
24*a2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 48*a2*d*pow(M_E,d*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) + 
48*a2*d*pow(M_E,3*d*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) + 9*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) - 
9*a2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) + 15*d*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 
15*d2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 15*d*d2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,3) - 
15*d2*k2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,3) - 30*d*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 
312*a2*pow(d,4)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 72*a2*pow(d,4)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,3) + 
72*a2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 
72*a2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 192*a2*d*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) + 
192*a2*d*pow(M_E,3*d*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) + 48*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) - 
48*a2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) + 15*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 
144*a2*k2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) - 15*d2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,4) - 
144*a2*k2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,4) - 15*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 
912*a2*pow(d,3)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 15*d2*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,4) - 
240*a2*pow(d,3)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,4) + 576*a2*d*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 
576*a2*d*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 192*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) - 
192*a2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) + 
c2*(-2*pow(d,5)*pow(M_E,k2*pow(RSP,-1)) - 2*RSP*pow(d,4)*(6 + pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k2*pow(RSP,-1)) - 
3*pow(d,3)*(13 + 3*pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) - 
pow(d,2)*(80*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 16*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,3) - 
pow(M_E,d*pow(RSP,-1))*(4*RSP*pow(k2,2) + pow(k2,3) + 8*k2*pow(RSP,2) + 8*pow(RSP,3)) + 
pow(M_E,3*d*pow(RSP,-1))*(4*RSP*pow(k2,2) + pow(k2,3) + 8*k2*pow(RSP,2) + 8*pow(RSP,3))) - 
RSP*(-1 + pow(M_E,2*d*pow(RSP,-1)))*(-48*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 
pow(M_E,d*pow(RSP,-1))*(7*RSP*pow(k2,3) + pow(k2,4) + 24*pow(k2,2)*pow(RSP,2) + 48*k2*pow(RSP,3) + 
48*pow(RSP,4))) + d*(-96*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 
pow(M_E,d*pow(RSP,-1))*(7*RSP*pow(k2,3) + pow(k2,4) + 24*pow(k2,2)*pow(RSP,2) + 48*k2*pow(RSP,3) + 
48*pow(RSP,4)) + pow(M_E,3*d*pow(RSP,-1))*
(7*RSP*pow(k2,3) + pow(k2,4) + 24*pow(k2,2)*pow(RSP,2) + 48*k2*pow(RSP,3) + 48*pow(RSP,4)))) + 
1152*a2*d*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 144*a2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 
1152*a2*d*k2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,5) - 144*a2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,5) - 
1872*a2*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 432*a2*pow(d,2)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,5) + 
576*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) - 576*a2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) + 
b2*(-2*pow(d,6)*pow(M_E,k2*pow(RSP,-1)) - 2*RSP*pow(d,5)*(7 + pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k2*pow(RSP,-1)) - 
2*pow(d,4)*(29 + 6*pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) - 
10*pow(d,3)*(17 + 4*pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 
pow(d,2)*(345*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 75*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,4) - 
pow(M_E,d*pow(RSP,-1))*(5*RSP*pow(k2,3) + pow(k2,4) + 15*pow(k2,2)*pow(RSP,2) + 30*k2*pow(RSP,3) + 
30*pow(RSP,4)) + pow(M_E,3*d*pow(RSP,-1))*
(5*RSP*pow(k2,3) + pow(k2,4) + 15*pow(k2,2)*pow(RSP,2) + 30*k2*pow(RSP,3) + 30*pow(RSP,4))) - 
RSP*(-1 + pow(M_E,2*d*pow(RSP,-1)))*(-210*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 
pow(M_E,d*pow(RSP,-1))*(8*RSP*pow(k2,4) + pow(k2,5) + 35*pow(k2,3)*pow(RSP,2) + 105*pow(k2,2)*pow(RSP,3) + 
210*k2*pow(RSP,4) + 210*pow(RSP,5))) + 
d*(-420*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 
pow(M_E,d*pow(RSP,-1))*(8*RSP*pow(k2,4) + pow(k2,5) + 35*pow(k2,3)*pow(RSP,2) + 105*pow(k2,2)*pow(RSP,3) + 
210*k2*pow(RSP,4) + 210*pow(RSP,5)) + 
pow(M_E,3*d*pow(RSP,-1))*(8*RSP*pow(k2,4) + pow(k2,5) + 35*pow(k2,3)*pow(RSP,2) + 105*pow(k2,2)*pow(RSP,3) + 
210*k2*pow(RSP,4) + 210*pow(RSP,5)))) + 1152*a2*d*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) + 
1152*a2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) + 1152*a2*d*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,6) - 
1152*a2*k2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,6) - 2304*a2*d*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) + 
1152*a2*pow(M_E,d*pow(RSP,-1))*pow(RSP,7) - 1152*a2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,7) - 
1152*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,7) + 1152*a2*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,7)) + 
7560*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(RSP,8) + 57600*a1*a2*d*pow(M_E,d*pow(RSP,-1))*pow(RSP,8) + 
57600*a1*a2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,8) - 7560*a1*b2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,8) + 
57600*a1*a2*d*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,8) - 57600*a1*a2*k2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,8) - 
7560*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,8) - 115200*a1*a2*d*pow(M_E,k2*pow(RSP,-1))*pow(RSP,8) + 
7560*a1*b2*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,8) + 
b1*(-12*d2*RSP*pow(d,4)*pow(M_E,k2*pow(RSP,-1)) - 2*d2*pow(d,5)*pow(M_E,k2*pow(RSP,-1)) - 
18*a2*RSP*pow(d,7)*pow(M_E,k2*pow(RSP,-1)) - 2*a2*pow(d,8)*pow(M_E,k2*pow(RSP,-1)) - 
2*d2*RSP*pow(d,4)*pow(M_E,(2*d + k2)*pow(RSP,-1)) - 2*a2*RSP*pow(d,7)*pow(M_E,(2*d + k2)*pow(RSP,-1)) + 
4*d2*RSP*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,2) - 4*d2*RSP*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2) + 
7*d*d2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,3) + d2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,3) + 
7*d*d2*RSP*pow(M_E,3*d*pow(RSP,-1))*pow(k2,3) - d2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,3) + 
d*d2*pow(M_E,d*pow(RSP,-1))*pow(k2,4) + d2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,4) + 
d*d2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,4) - d2*RSP*pow(M_E,3*d*pow(RSP,-1))*pow(k2,4) + 
7*a2*RSP*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,5) - 7*a2*RSP*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,5) + 
10*a2*d*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,6) + a2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,6) + 
10*a2*d*RSP*pow(M_E,3*d*pow(RSP,-1))*pow(k2,6) - a2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,6) + 
a2*d*pow(M_E,d*pow(RSP,-1))*pow(k2,7) + a2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,7) + 
a2*d*pow(M_E,3*d*pow(RSP,-1))*pow(k2,7) - a2*RSP*pow(M_E,3*d*pow(RSP,-1))*pow(k2,7) + 
8*d2*k2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) - 8*d2*k2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,2) - 
39*d2*pow(d,3)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) - 108*a2*pow(d,6)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) - 
9*d2*pow(d,3)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,2) - 
18*a2*pow(d,6)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,2) + 24*d*d2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 
24*d*d2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 7*d2*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 
7*d2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 35*a2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) - 
35*a2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) + 63*a2*d*pow(M_E,d*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) + 
63*a2*d*pow(M_E,3*d*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) + 10*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,6)*pow(RSP,2) - 
10*a2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,6)*pow(RSP,2) + 48*d*d2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 
8*d2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 48*d*d2*k2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,3) - 
8*d2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,3) - 80*d2*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 
518*a2*pow(d,5)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 16*d2*pow(d,2)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,3) - 
112*a2*pow(d,5)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,3) + 24*d2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 
24*d2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 140*a2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) - 
140*a2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) + 
315*a2*d*pow(M_E,d*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) + 315*a2*d*pow(M_E,3*d*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) + 
63*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,5)*pow(RSP,3) - 63*a2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,5)*pow(RSP,3) + 
48*d*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 48*d2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 
48*d*d2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,4) - 48*d2*k2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,4) - 
96*d*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 1995*a2*pow(d,4)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 
525*a2*pow(d,4)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,4) + 
420*a2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) - 
420*a2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 
1260*a2*d*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) + 1260*a2*d*pow(M_E,3*d*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) + 
315*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,4)*pow(RSP,4) - 315*a2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,4)*pow(RSP,4) + 
48*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 840*a2*k2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) - 
48*d2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,5) - 840*a2*k2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,5) - 
48*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 5880*a2*pow(d,3)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 
48*d2*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,5) - 1680*a2*pow(d,3)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,5) + 
3780*a2*d*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) + 3780*a2*d*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) + 
1260*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,5) - 1260*a2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,3)*pow(RSP,5) + 
c2*(-2*pow(d,6)*pow(M_E,k2*pow(RSP,-1)) - 2*RSP*pow(d,5)*(7 + pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k2*pow(RSP,-1)) - 
2*pow(d,4)*(29 + 6*pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) - 
10*pow(d,3)*(17 + 4*pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 
pow(d,2)*(345*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 75*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,4) - 
pow(M_E,d*pow(RSP,-1))*(5*RSP*pow(k2,3) + pow(k2,4) + 15*pow(k2,2)*pow(RSP,2) + 30*k2*pow(RSP,3) + 
30*pow(RSP,4)) + pow(M_E,3*d*pow(RSP,-1))*
(5*RSP*pow(k2,3) + pow(k2,4) + 15*pow(k2,2)*pow(RSP,2) + 30*k2*pow(RSP,3) + 30*pow(RSP,4))) - 
RSP*(-1 + pow(M_E,2*d*pow(RSP,-1)))*(-210*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 
pow(M_E,d*pow(RSP,-1))*(8*RSP*pow(k2,4) + pow(k2,5) + 35*pow(k2,3)*pow(RSP,2) + 105*pow(k2,2)*pow(RSP,3) + 
210*k2*pow(RSP,4) + 210*pow(RSP,5))) + 
d*(-420*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 
pow(M_E,d*pow(RSP,-1))*(8*RSP*pow(k2,4) + pow(k2,5) + 35*pow(k2,3)*pow(RSP,2) + 105*pow(k2,2)*pow(RSP,3) + 
210*k2*pow(RSP,4) + 210*pow(RSP,5)) + 
pow(M_E,3*d*pow(RSP,-1))*(8*RSP*pow(k2,4) + pow(k2,5) + 35*pow(k2,3)*pow(RSP,2) + 105*pow(k2,2)*pow(RSP,3) + 
210*k2*pow(RSP,4) + 210*pow(RSP,5)))) + 7560*a2*d*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) + 
840*a2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) + 7560*a2*d*k2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,6) - 
840*a2*pow(d,2)*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,6) - 12180*a2*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) - 
2940*a2*pow(d,2)*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,6) + 3780*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,6) - 
3780*a2*pow(M_E,3*d*pow(RSP,-1))*pow(k2,2)*pow(RSP,6) + 
b2*(-2*pow(d,7)*pow(M_E,k2*pow(RSP,-1)) - 2*RSP*pow(d,6)*(8 + pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k2*pow(RSP,-1)) - 
3*pow(d,5)*(27 + 5*pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) - 
24*pow(d,4)*(13 + 3*pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 
48*pow(d,3)*(19 + 5*pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 
pow(d,2)*(1872*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 432*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,5) - 
pow(M_E,d*pow(RSP,-1))*(6*RSP*pow(k2,4) + pow(k2,5) + 24*pow(k2,3)*pow(RSP,2) + 72*pow(k2,2)*pow(RSP,3) + 
144*k2*pow(RSP,4) + 144*pow(RSP,5)) + 
pow(M_E,3*d*pow(RSP,-1))*(6*RSP*pow(k2,4) + pow(k2,5) + 24*pow(k2,3)*pow(RSP,2) + 72*pow(k2,2)*pow(RSP,3) + 
144*k2*pow(RSP,4) + 144*pow(RSP,5))) - 
RSP*(-1 + pow(M_E,2*d*pow(RSP,-1)))*(-1152*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) + 
pow(M_E,d*pow(RSP,-1))*(9*RSP*pow(k2,5) + pow(k2,6) + 48*pow(k2,4)*pow(RSP,2) + 192*pow(k2,3)*pow(RSP,3) + 
576*pow(k2,2)*pow(RSP,4) + 1152*k2*pow(RSP,5) + 1152*pow(RSP,6))) + 
d*(-2304*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) + 
pow(M_E,d*pow(RSP,-1))*(9*RSP*pow(k2,5) + pow(k2,6) + 48*pow(k2,4)*pow(RSP,2) + 192*pow(k2,3)*pow(RSP,3) + 
576*pow(k2,2)*pow(RSP,4) + 1152*k2*pow(RSP,5) + 1152*pow(RSP,6)) + 
pow(M_E,3*d*pow(RSP,-1))*(9*RSP*pow(k2,5) + pow(k2,6) + 48*pow(k2,4)*pow(RSP,2) + 192*pow(k2,3)*pow(RSP,3) + 
576*pow(k2,2)*pow(RSP,4) + 1152*k2*pow(RSP,5) + 1152*pow(RSP,6)))) + 
7560*a2*d*pow(M_E,d*pow(RSP,-1))*pow(RSP,7) + 7560*a2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,7) + 
7560*a2*d*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,7) - 7560*a2*k2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,7) - 
15120*a2*d*pow(M_E,k2*pow(RSP,-1))*pow(RSP,7) + 7560*a2*pow(M_E,d*pow(RSP,-1))*pow(RSP,8) - 
7560*a2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,8) - 7560*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,8) + 
7560*a2*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,8)) + 57600*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(RSP,9) - 
57600*a1*a2*pow(M_E,3*d*pow(RSP,-1))*pow(RSP,9) - 57600*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,9) + 
57600*a1*a2*pow(M_E,(2*d + k2)*pow(RSP,-1))*pow(RSP,9))*pow(2*pow(d,2),-1));

    return res;
}

//  3
double BMSH_Integral_14_case_2_sub_2( double ASP, double RSP, double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d )
{
double res = ASP*pow(3,0.5)*pow(M_E,-((2*d + k1)*pow(RSP,-1)))*pow(RSP,2)*
(5*d*d1*d2*k1*RSP*pow(M_E,d*pow(RSP,-1)) + d1*d2*k1*pow(d,2)*pow(M_E,d*pow(RSP,-1)) + 
2*d1*d2*RSP*pow(d,2)*pow(M_E,d*pow(RSP,-1)) + 3*c2*d1*k1*RSP*pow(d,2)*pow(M_E,d*pow(RSP,-1)) - 
8*d1*d2*RSP*pow(d,2)*pow(M_E,k1*pow(RSP,-1)) - 2*d1*d2*pow(d,3)*pow(M_E,k1*pow(RSP,-1)) - 
10*c2*d1*RSP*pow(d,3)*pow(M_E,k1*pow(RSP,-1)) - 2*c2*d1*pow(d,4)*pow(M_E,k1*pow(RSP,-1)) - 
12*b2*d1*RSP*pow(d,4)*pow(M_E,k1*pow(RSP,-1)) - 2*b2*d1*pow(d,5)*pow(M_E,k1*pow(RSP,-1)) - 
14*a2*d1*RSP*pow(d,5)*pow(M_E,k1*pow(RSP,-1)) - 14*a1*d2*RSP*pow(d,5)*pow(M_E,k1*pow(RSP,-1)) - 
2*a2*d1*pow(d,6)*pow(M_E,k1*pow(RSP,-1)) - 2*a1*d2*pow(d,6)*pow(M_E,k1*pow(RSP,-1)) - 
16*a1*c2*RSP*pow(d,6)*pow(M_E,k1*pow(RSP,-1)) - 2*a1*c2*pow(d,7)*pow(M_E,k1*pow(RSP,-1)) - 
18*a1*b2*RSP*pow(d,7)*pow(M_E,k1*pow(RSP,-1)) - 2*a1*b2*pow(d,8)*pow(M_E,k1*pow(RSP,-1)) - 
20*a1*a2*RSP*pow(d,8)*pow(M_E,k1*pow(RSP,-1)) - 2*a1*a2*pow(d,9)*pow(M_E,k1*pow(RSP,-1)) + 
2*d1*d2*RSP*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + 2*c2*d1*RSP*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + 
2*b2*d1*RSP*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + 2*a2*d1*RSP*pow(d,5)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + 
2*a1*d2*RSP*pow(d,5)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + 2*a1*c2*RSP*pow(d,6)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + 
2*a1*b2*RSP*pow(d,7)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + 2*a1*a2*RSP*pow(d,8)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 
5*d*d1*d2*k1*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1)) - d1*d2*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1)) + 
2*d1*d2*RSP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1)) + 3*c2*d1*k1*RSP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1)) + 
d*d1*d2*pow(M_E,d*pow(RSP,-1))*pow(k1,2) + 6*c2*d*d1*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,2) + 
d1*d2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,2) + c2*d1*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,2) + 
4*b2*d1*RSP*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,2) + d*d1*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2) - 
6*c2*d*d1*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2) + d1*d2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2) - 
c2*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2) + 
4*b2*d1*RSP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2) + c2*d*d1*pow(M_E,d*pow(RSP,-1))*pow(k1,3) + 
c2*d1*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,3) + 7*b2*d*d1*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,3) + 
b2*d1*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,3) + 5*a2*d1*RSP*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,3) + 
5*a1*d2*RSP*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,3) + c2*d*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3) + 
c2*d1*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3) - 7*b2*d*d1*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3) - 
b2*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3) + 
5*a2*d1*RSP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3) + 
5*a1*d2*RSP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3) + b2*d*d1*pow(M_E,d*pow(RSP,-1))*pow(k1,4) + 
b2*d1*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,4) + 8*a2*d*d1*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,4) + 
8*a1*d*d2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,4) + a2*d1*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,4) + 
a1*d2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,4) + 6*a1*c2*RSP*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,4) + 
b2*d*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4) + b2*d1*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4) - 
8*a2*d*d1*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4) - 8*a1*d*d2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4) - 
a2*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4) - a1*d2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4) + 
6*a1*c2*RSP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4) + a2*d*d1*pow(M_E,d*pow(RSP,-1))*pow(k1,5) + 
a1*d*d2*pow(M_E,d*pow(RSP,-1))*pow(k1,5) + 9*a1*c2*d*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,5) + 
a2*d1*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,5) + a1*d2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,5) + 
a1*c2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,5) + 7*a1*b2*RSP*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,5) + 
a2*d*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5) + a1*d*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5) - 
9*a1*c2*d*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5) + a2*d1*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5) + 
a1*d2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5) - a1*c2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5) + 
7*a1*b2*RSP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5) + a1*c2*d*pow(M_E,d*pow(RSP,-1))*pow(k1,6) + 
a1*c2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,6) + 10*a1*b2*d*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,6) + 
a1*b2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,6) + 8*a1*a2*RSP*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,6) + 
a1*c2*d*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,6) + a1*c2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,6) - 
10*a1*b2*d*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,6) - a1*b2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,6) + 
8*a1*a2*RSP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,6) + a1*b2*d*pow(M_E,d*pow(RSP,-1))*pow(k1,7) + 
a1*b2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,7) + 11*a1*a2*d*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,7) + 
a1*a2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,7) + a1*b2*d*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,7) + 
a1*b2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,7) - 11*a1*a2*d*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,7) - 
a1*a2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,7) + a1*a2*d*pow(M_E,d*pow(RSP,-1))*pow(k1,8) + 
a1*a2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,8) + a1*a2*d*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,8) + 
a1*a2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,8) + 8*d*d1*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) + 
15*c2*d*d1*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) + 5*d1*d2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) + 
3*c2*d1*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) + 8*b2*d1*k1*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) - 
13*d*d1*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 24*c2*d1*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 
39*b2*d1*pow(d,3)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 58*a2*d1*pow(d,4)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 
58*a1*d2*pow(d,4)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 81*a1*c2*pow(d,5)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 
108*a1*b2*pow(d,6)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 139*a1*a2*pow(d,7)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 
3*d*d1*d2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) - 6*c2*d1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) - 
9*b2*d1*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) - 
12*a2*d1*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) - 
12*a1*d2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) - 
15*a1*c2*pow(d,5)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) - 
18*a1*b2*pow(d,6)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) - 
21*a1*a2*pow(d,7)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) + 8*d*d1*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,2) + 
15*c2*d*d1*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,2) - 5*d1*d2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,2) - 
3*c2*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,2) - 
8*b2*d1*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,2) + 6*c2*d1*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 
24*b2*d*d1*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 15*a2*d1*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 
15*a1*d2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 
6*c2*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 
24*b2*d*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 
15*a2*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 
15*a1*d2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 
7*b2*d1*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 35*a2*d*d1*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 
35*a1*d*d2*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 24*a1*c2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 
7*b2*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 
35*a2*d*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 
35*a1*d*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 
24*a1*c2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 
48*a1*c2*d*pow(M_E,d*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 8*a2*d1*pow(M_E,d*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 
8*a1*d2*pow(M_E,d*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 35*a1*b2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 
48*a1*c2*d*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) - 
8*a2*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) - 
8*a1*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) - 
35*a1*b2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 
9*a1*c2*pow(M_E,d*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) + 63*a1*b2*d*pow(M_E,d*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) + 
48*a1*a2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) - 
9*a1*c2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) + 
63*a1*b2*d*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) - 
48*a1*a2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) + 
10*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k1,6)*pow(RSP,2) + 80*a1*a2*d*pow(M_E,d*pow(RSP,-1))*pow(k1,6)*pow(RSP,2) - 
10*a1*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,6)*pow(RSP,2) + 
80*a1*a2*d*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,6)*pow(RSP,2) + 11*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,7)*pow(RSP,2) - 
11*a1*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,7)*pow(RSP,2) + 15*c2*d*d1*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 
8*d1*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 15*c2*d1*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 
48*b2*d*d1*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 8*b2*d1*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 
30*a2*d1*k1*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 30*a1*d2*k1*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) - 
30*c2*d*d1*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 8*d1*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
80*b2*d1*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 170*a2*d1*pow(d,3)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
170*a1*d2*pow(d,3)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 312*a1*c2*pow(d,4)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
518*a1*b2*pow(d,5)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 800*a1*a2*pow(d,6)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
8*d1*d2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) + 16*b2*d1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) + 
40*a2*d1*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) + 
40*a1*d2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) + 
72*a1*c2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) + 
112*a1*b2*pow(d,5)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) + 
160*a1*a2*pow(d,6)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) - 15*c2*d*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,3) + 
8*d1*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,3) + 15*c2*d1*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,3) - 
48*b2*d*d1*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,3) + 8*b2*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,3) + 
30*a2*d1*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,3) + 
30*a1*d2*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,3) + 24*b2*d1*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
105*a2*d*d1*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 105*a1*d*d2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
72*a1*c2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
24*b2*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) - 
105*a2*d*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) - 
105*a1*d*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
72*a1*c2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
192*a1*c2*d*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 35*a2*d1*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 
35*a1*d2*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 140*a1*b2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) - 
192*a1*c2*d*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 
35*a2*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 
35*a1*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 
140*a1*b2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 
48*a1*c2*pow(M_E,d*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) + 315*a1*b2*d*pow(M_E,d*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) + 
240*a1*a2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) + 
48*a1*c2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) - 
315*a1*b2*d*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) + 
240*a1*a2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) + 
63*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k1,5)*pow(RSP,3) + 480*a1*a2*d*pow(M_E,d*pow(RSP,-1))*pow(k1,5)*pow(RSP,3) + 
63*a1*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5)*pow(RSP,3) - 
480*a1*a2*d*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5)*pow(RSP,3) + 80*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,6)*pow(RSP,3) + 
80*a1*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,6)*pow(RSP,3) + 15*c2*d1*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 
48*b2*d*d1*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 48*b2*d1*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 
210*a2*d*d1*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 210*a1*d*d2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 
30*a2*d1*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 30*a1*d2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 
144*a1*c2*k1*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) - 15*c2*d1*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 
96*b2*d*d1*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 345*a2*d1*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 
345*a1*d2*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 912*a1*c2*pow(d,3)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 
1995*a1*b2*pow(d,4)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 3840*a1*a2*pow(d,5)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 
15*c2*d1*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) - 75*a2*d1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) - 
75*a1*d2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) - 
240*a1*c2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) - 
525*a1*b2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) - 
960*a1*a2*pow(d,5)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) - 15*c2*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,4) + 
48*b2*d*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,4) - 48*b2*d1*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,4) + 
210*a2*d*d1*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,4) + 210*a1*d*d2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,4) - 
30*a2*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,4) - 
30*a1*d2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,4) - 
144*a1*c2*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,4) + 
576*a1*c2*d*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 105*a2*d1*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 
105*a1*d2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 420*a1*b2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 
576*a1*c2*d*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) - 
105*a2*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) - 
105*a1*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) - 
420*a1*b2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 
192*a1*c2*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) + 1260*a1*b2*d*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) + 
960*a1*a2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) - 
192*a1*c2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) + 
1260*a1*b2*d*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) - 
960*a1*a2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) + 
315*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k1,4)*pow(RSP,4) + 2400*a1*a2*d*pow(M_E,d*pow(RSP,-1))*pow(k1,4)*pow(RSP,4) - 
315*a1*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4)*pow(RSP,4) + 
2400*a1*a2*d*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4)*pow(RSP,4) + 
480*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,5)*pow(RSP,4) - 480*a1*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5)*pow(RSP,4) + 
48*b2*d1*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 210*a2*d*d1*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 
210*a1*d*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 1152*a1*c2*d*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 
210*a2*d1*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 210*a1*d2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 
144*a1*c2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 840*a1*b2*k1*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) - 
48*b2*d1*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 420*a2*d*d1*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 
420*a1*d*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 1872*a1*c2*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 
5880*a1*b2*pow(d,3)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 14880*a1*a2*pow(d,4)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 
48*b2*d1*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,5) + 432*a1*c2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,5) + 
1680*a1*b2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,5) + 
4320*a1*a2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,5) + 48*b2*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,5) - 
210*a2*d*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,5) - 210*a1*d*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,5) - 
1152*a1*c2*d*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,5) + 210*a2*d1*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,5) + 
210*a1*d2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,5) + 144*a1*c2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,5) + 
840*a1*b2*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,5) + 576*a1*c2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) + 
3780*a1*b2*d*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) + 
2880*a1*a2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) + 
576*a1*c2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) - 
3780*a1*b2*d*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) + 
2880*a1*a2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) + 
1260*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,5) + 9600*a1*a2*d*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,5) + 
1260*a1*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,5) - 
9600*a1*a2*d*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,5) + 
2400*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,4)*pow(RSP,5) + 2400*a1*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4)*pow(RSP,5) + 
1152*a1*c2*d*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) + 210*a2*d1*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) + 
210*a1*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) + 1152*a1*c2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) + 
7560*a1*b2*d*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) + 840*a1*b2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) + 
5760*a1*a2*k1*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) - 2304*a1*c2*d*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) - 
210*a2*d1*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) - 210*a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) - 
12180*a1*b2*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) - 44160*a1*a2*pow(d,3)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) + 
210*a2*d1*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,6) + 210*a1*d2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,6) - 
2940*a1*b2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,6) - 
13440*a1*a2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,6) + 
1152*a1*c2*d*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,6) - 210*a2*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,6) - 
210*a1*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,6) - 1152*a1*c2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,6) + 
7560*a1*b2*d*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,6) - 
840*a1*b2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,6) - 
5760*a1*a2*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,6) + 
3780*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,6) + 28800*a1*a2*d*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,6) - 
3780*a1*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,6) + 
28800*a1*a2*d*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,6) + 
9600*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,6) - 9600*a1*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,6) + 
1152*a1*c2*pow(M_E,d*pow(RSP,-1))*pow(RSP,7) + 7560*a1*b2*d*pow(M_E,d*pow(RSP,-1))*pow(RSP,7) + 
7560*a1*b2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,7) + 57600*a1*a2*d*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,7) + 
5760*a1*a2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,7) - 1152*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,7) - 
15120*a1*b2*d*pow(M_E,k1*pow(RSP,-1))*pow(RSP,7) - 92160*a1*a2*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,7) - 
1152*a1*c2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,7) + 23040*a1*a2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,7) + 
1152*a1*c2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,7) - 7560*a1*b2*d*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,7) + 
7560*a1*b2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,7) - 57600*a1*a2*d*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,7) + 
5760*a1*a2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,7) + 28800*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,7) + 
28800*a1*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,7) + 
c1*(3*d2*k1*RSP*pow(d,2)*pow(M_E,d*pow(RSP,-1)) - 10*d2*RSP*pow(d,3)*pow(M_E,k1*pow(RSP,-1)) - 
2*d2*pow(d,4)*pow(M_E,k1*pow(RSP,-1)) - 16*a2*RSP*pow(d,6)*pow(M_E,k1*pow(RSP,-1)) - 
2*a2*pow(d,7)*pow(M_E,k1*pow(RSP,-1)) + 2*d2*RSP*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + 
2*a2*RSP*pow(d,6)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + 3*d2*k1*RSP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1)) + 
6*d*d2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,2) + d2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,2) - 
6*d*d2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2) - d2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2) + 
d*d2*pow(M_E,d*pow(RSP,-1))*pow(k1,3) + d2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,3) + 
d*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3) + d2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3) + 
6*a2*RSP*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,4) + 6*a2*RSP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4) + 
9*a2*d*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,5) + a2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,5) - 
9*a2*d*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5) - a2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5) + 
a2*d*pow(M_E,d*pow(RSP,-1))*pow(k1,6) + a2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,6) + 
a2*d*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,6) + a2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,6) + 
15*d*d2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) + 3*d2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) - 
24*d2*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 81*a2*pow(d,5)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 
6*d2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) - 15*a2*pow(d,5)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) + 
15*d*d2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,2) - 3*d2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,2) + 
6*d2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 6*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 
24*a2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 
24*a2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 
48*a2*d*pow(M_E,d*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 48*a2*d*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 
9*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) - 9*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) + 
15*d*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 15*d2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) - 
30*d*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 312*a2*pow(d,4)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) + 
72*a2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) - 15*d*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,3) + 
15*d2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,3) + 72*a2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
72*a2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
192*a2*d*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) - 192*a2*d*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 
48*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) + 48*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) + 
15*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 144*a2*k1*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) - 
15*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 912*a2*pow(d,3)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 
15*d2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) - 240*a2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) - 
15*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,4) - 144*a2*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,4) + 
576*a2*d*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 576*a2*d*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 
192*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) - 192*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) + 
c2*(-2*pow(d,5)*pow(M_E,k1*pow(RSP,-1)) + 2*RSP*pow(d,4)*(-6 + pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1)) - 
3*pow(d,3)*(13 + 3*pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) + 
pow(d,2)*(-80*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) + 16*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) + 
pow(M_E,(d + 2*k1)*pow(RSP,-1))*(4*RSP*pow(k1,2) - pow(k1,3) - 8*k1*pow(RSP,2) + 8*pow(RSP,3)) + 
pow(M_E,d*pow(RSP,-1))*(4*RSP*pow(k1,2) + pow(k1,3) + 8*k1*pow(RSP,2) + 8*pow(RSP,3))) + 
d*(-96*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 
pow(M_E,(d + 2*k1)*pow(RSP,-1))*(-7*RSP*pow(k1,3) + pow(k1,4) + 24*pow(k1,2)*pow(RSP,2) - 48*k1*pow(RSP,3) + 
48*pow(RSP,4)) + pow(M_E,d*pow(RSP,-1))*
(7*RSP*pow(k1,3) + pow(k1,4) + 24*pow(k1,2)*pow(RSP,2) + 48*k1*pow(RSP,3) + 48*pow(RSP,4))) + 
RSP*(-48*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 48*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) + 
pow(M_E,(d + 2*k1)*pow(RSP,-1))*(-7*RSP*pow(k1,3) + pow(k1,4) + 24*pow(k1,2)*pow(RSP,2) - 48*k1*pow(RSP,3) + 
48*pow(RSP,4)) + pow(M_E,d*pow(RSP,-1))*
(7*RSP*pow(k1,3) + pow(k1,4) + 24*pow(k1,2)*pow(RSP,2) + 48*k1*pow(RSP,3) + 48*pow(RSP,4)))) + 
1152*a2*d*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 144*a2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) - 
1872*a2*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 432*a2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,5) - 
1152*a2*d*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,5) + 144*a2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,5) + 
576*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) + 576*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) + 
b2*(-2*pow(d,6)*pow(M_E,k1*pow(RSP,-1)) + 2*RSP*pow(d,5)*(-7 + pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1)) - 
2*pow(d,4)*(29 + 6*pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) + 
10*pow(d,3)*(-17 + 4*pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) + 
pow(d,2)*(-345*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 75*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) - 
pow(M_E,(d + 2*k1)*pow(RSP,-1))*(-5*RSP*pow(k1,3) + pow(k1,4) + 15*pow(k1,2)*pow(RSP,2) - 30*k1*pow(RSP,3) + 
30*pow(RSP,4)) + pow(M_E,d*pow(RSP,-1))*
(5*RSP*pow(k1,3) + pow(k1,4) + 15*pow(k1,2)*pow(RSP,2) + 30*k1*pow(RSP,3) + 30*pow(RSP,4))) + 
d*(pow(M_E,(d + 2*k1)*pow(RSP,-1))*(-8*RSP*pow(k1,4) + pow(k1,5) + 35*pow(k1,3)*pow(RSP,2) - 
105*pow(k1,2)*pow(RSP,3) + 210*k1*pow(RSP,4) - 210*pow(RSP,5)) - 420*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 
pow(M_E,d*pow(RSP,-1))*(8*RSP*pow(k1,4) + pow(k1,5) + 35*pow(k1,3)*pow(RSP,2) + 105*pow(k1,2)*pow(RSP,3) + 
210*k1*pow(RSP,4) + 210*pow(RSP,5))) + 
RSP*(pow(M_E,(d + 2*k1)*pow(RSP,-1))*(-8*RSP*pow(k1,4) + pow(k1,5) + 35*pow(k1,3)*pow(RSP,2) - 
105*pow(k1,2)*pow(RSP,3) + 210*k1*pow(RSP,4) - 210*pow(RSP,5)) - 210*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 
210*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,5) + 
pow(M_E,d*pow(RSP,-1))*(8*RSP*pow(k1,4) + pow(k1,5) + 35*pow(k1,3)*pow(RSP,2) + 105*pow(k1,2)*pow(RSP,3) + 
210*k1*pow(RSP,4) + 210*pow(RSP,5)))) + 1152*a2*d*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) + 
1152*a2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) - 2304*a2*d*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) + 
1152*a2*d*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,6) - 1152*a2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,6) + 
1152*a2*pow(M_E,d*pow(RSP,-1))*pow(RSP,7) - 1152*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,7) - 
1152*a2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,7) + 1152*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,7)) + 
7560*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(RSP,8) + 57600*a1*a2*d*pow(M_E,d*pow(RSP,-1))*pow(RSP,8) + 
57600*a1*a2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,8) - 7560*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,8) - 
115200*a1*a2*d*pow(M_E,k1*pow(RSP,-1))*pow(RSP,8) + 7560*a1*b2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,8) - 
7560*a1*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,8) + 57600*a1*a2*d*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,8) - 
57600*a1*a2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,8) + 
b1*(-12*d2*RSP*pow(d,4)*pow(M_E,k1*pow(RSP,-1)) - 2*d2*pow(d,5)*pow(M_E,k1*pow(RSP,-1)) - 
18*a2*RSP*pow(d,7)*pow(M_E,k1*pow(RSP,-1)) - 2*a2*pow(d,8)*pow(M_E,k1*pow(RSP,-1)) + 
2*d2*RSP*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + 2*a2*RSP*pow(d,7)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + 
4*d2*RSP*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,2) + 4*d2*RSP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2) + 
7*d*d2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,3) + d2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,3) - 
7*d*d2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3) - d2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3) + 
d*d2*pow(M_E,d*pow(RSP,-1))*pow(k1,4) + d2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,4) + 
d*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4) + d2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4) + 
7*a2*RSP*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,5) + 7*a2*RSP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5) + 
10*a2*d*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,6) + a2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,6) - 
10*a2*d*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,6) - a2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,6) + 
a2*d*pow(M_E,d*pow(RSP,-1))*pow(k1,7) + a2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,7) + 
a2*d*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,7) + a2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,7) + 
8*d2*k1*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) - 39*d2*pow(d,3)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 
108*a2*pow(d,6)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 9*d2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) - 
18*a2*pow(d,6)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) - 
8*d2*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,2) + 24*d*d2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 
24*d*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 7*d2*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 
7*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 
35*a2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) - 
35*a2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 
63*a2*d*pow(M_E,d*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) + 63*a2*d*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) + 
10*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,6)*pow(RSP,2) - 10*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,6)*pow(RSP,2) + 
48*d*d2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 8*d2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) - 
80*d2*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 518*a2*pow(d,5)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) + 
16*d2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) + 
112*a2*pow(d,5)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) - 48*d*d2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,3) + 
8*d2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,3) + 24*d2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
24*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
140*a2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 
140*a2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 
315*a2*d*pow(M_E,d*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) - 315*a2*d*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) + 
63*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,5)*pow(RSP,3) + 63*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5)*pow(RSP,3) + 
48*d*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 48*d2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) - 
96*d*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 1995*a2*pow(d,4)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 
525*a2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) + 48*d*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,4) - 
48*d2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,4) + 420*a2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) - 
420*a2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 
1260*a2*d*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) + 
1260*a2*d*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) + 315*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,4)*pow(RSP,4) - 
315*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4)*pow(RSP,4) + 48*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 
840*a2*k1*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) - 48*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 
5880*a2*pow(d,3)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 48*d2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,5) + 
1680*a2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,5) + 48*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,5) + 
840*a2*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,5) + 3780*a2*d*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) - 
3780*a2*d*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) + 1260*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,5) + 
1260*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,5) + 
c2*(-2*pow(d,6)*pow(M_E,k1*pow(RSP,-1)) + 2*RSP*pow(d,5)*(-7 + pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1)) - 
2*pow(d,4)*(29 + 6*pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) + 
10*pow(d,3)*(-17 + 4*pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) + 
pow(d,2)*(-345*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 75*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) - 
pow(M_E,(d + 2*k1)*pow(RSP,-1))*(-5*RSP*pow(k1,3) + pow(k1,4) + 15*pow(k1,2)*pow(RSP,2) - 30*k1*pow(RSP,3) + 
30*pow(RSP,4)) + pow(M_E,d*pow(RSP,-1))*
(5*RSP*pow(k1,3) + pow(k1,4) + 15*pow(k1,2)*pow(RSP,2) + 30*k1*pow(RSP,3) + 30*pow(RSP,4))) + 
d*(pow(M_E,(d + 2*k1)*pow(RSP,-1))*(-8*RSP*pow(k1,4) + pow(k1,5) + 35*pow(k1,3)*pow(RSP,2) - 
105*pow(k1,2)*pow(RSP,3) + 210*k1*pow(RSP,4) - 210*pow(RSP,5)) - 420*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 
pow(M_E,d*pow(RSP,-1))*(8*RSP*pow(k1,4) + pow(k1,5) + 35*pow(k1,3)*pow(RSP,2) + 105*pow(k1,2)*pow(RSP,3) + 
210*k1*pow(RSP,4) + 210*pow(RSP,5))) + 
RSP*(pow(M_E,(d + 2*k1)*pow(RSP,-1))*(-8*RSP*pow(k1,4) + pow(k1,5) + 35*pow(k1,3)*pow(RSP,2) - 
105*pow(k1,2)*pow(RSP,3) + 210*k1*pow(RSP,4) - 210*pow(RSP,5)) - 210*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 
210*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,5) + 
pow(M_E,d*pow(RSP,-1))*(8*RSP*pow(k1,4) + pow(k1,5) + 35*pow(k1,3)*pow(RSP,2) + 105*pow(k1,2)*pow(RSP,3) + 
210*k1*pow(RSP,4) + 210*pow(RSP,5)))) + 7560*a2*d*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) + 
840*a2*pow(d,2)*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) - 12180*a2*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) - 
2940*a2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,6) + 7560*a2*d*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,6) - 
840*a2*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,6) + 3780*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,6) - 
3780*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,6) + 
b2*(-2*pow(d,7)*pow(M_E,k1*pow(RSP,-1)) + 2*RSP*pow(d,6)*(-8 + pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1)) - 
3*pow(d,5)*(27 + 5*pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) + 
24*pow(d,4)*(-13 + 3*pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
48*pow(d,3)*(19 + 5*pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 
pow(d,2)*(-1872*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 432*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,5) + 
pow(M_E,(d + 2*k1)*pow(RSP,-1))*(6*RSP*pow(k1,4) - pow(k1,5) - 24*pow(k1,3)*pow(RSP,2) + 
72*pow(k1,2)*pow(RSP,3) - 144*k1*pow(RSP,4) + 144*pow(RSP,5)) + 
pow(M_E,d*pow(RSP,-1))*(6*RSP*pow(k1,4) + pow(k1,5) + 24*pow(k1,3)*pow(RSP,2) + 72*pow(k1,2)*pow(RSP,3) + 
144*k1*pow(RSP,4) + 144*pow(RSP,5))) + 
d*(-2304*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) + 
pow(M_E,(d + 2*k1)*pow(RSP,-1))*(-9*RSP*pow(k1,5) + pow(k1,6) + 48*pow(k1,4)*pow(RSP,2) - 
192*pow(k1,3)*pow(RSP,3) + 576*pow(k1,2)*pow(RSP,4) - 1152*k1*pow(RSP,5) + 1152*pow(RSP,6)) + 
pow(M_E,d*pow(RSP,-1))*(9*RSP*pow(k1,5) + pow(k1,6) + 48*pow(k1,4)*pow(RSP,2) + 192*pow(k1,3)*pow(RSP,3) + 
576*pow(k1,2)*pow(RSP,4) + 1152*k1*pow(RSP,5) + 1152*pow(RSP,6))) + 
RSP*(-1152*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) - 1152*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,6) + 
pow(M_E,(d + 2*k1)*pow(RSP,-1))*(-9*RSP*pow(k1,5) + pow(k1,6) + 48*pow(k1,4)*pow(RSP,2) - 
192*pow(k1,3)*pow(RSP,3) + 576*pow(k1,2)*pow(RSP,4) - 1152*k1*pow(RSP,5) + 1152*pow(RSP,6)) + 
pow(M_E,d*pow(RSP,-1))*(9*RSP*pow(k1,5) + pow(k1,6) + 48*pow(k1,4)*pow(RSP,2) + 192*pow(k1,3)*pow(RSP,3) + 
576*pow(k1,2)*pow(RSP,4) + 1152*k1*pow(RSP,5) + 1152*pow(RSP,6)))) + 
7560*a2*d*pow(M_E,d*pow(RSP,-1))*pow(RSP,7) + 7560*a2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,7) - 
15120*a2*d*pow(M_E,k1*pow(RSP,-1))*pow(RSP,7) - 7560*a2*d*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,7) + 
7560*a2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,7) + 7560*a2*pow(M_E,d*pow(RSP,-1))*pow(RSP,8) - 
7560*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,8) + 7560*a2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,8) - 
7560*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,8)) + 57600*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(RSP,9) - 
57600*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,9) - 57600*a1*a2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,9) + 
57600*a1*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,9))*pow(2*pow(d,2),-1);

    return res;
}

//  4
double BMSH_Integral_14_case_3( double ASP, double RSP, double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d )
{
double res = ASP*pow(3,0.5)*pow(M_E,-((d + k1 + k2)*pow(RSP,-1)))*pow(RSP,2)*
(-(pow(d,2)*(-1 + pow(M_E,2*d*pow(RSP,-1)))*(4*b1*d2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,2) + 
b1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + 5*b1*c2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + 
5*a1*d2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + 
a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + 6*b1*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + 
6*a1*c2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + 
a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + 7*a2*b1*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + 
7*a1*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) + 
a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) + 8*a1*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) + 
a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,7) - 4*b1*d2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,2) - 
b1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - 5*b1*c2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - 
5*a1*d2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 
a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 6*b1*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 
6*a1*c2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - 
a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - 7*a2*b1*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - 
7*a1*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) - 
a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) - 8*a1*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) - 
a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,7) - 8*b1*d2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) + 
8*b1*d2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 15*b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 
15*a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 24*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 
24*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 35*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 
35*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 48*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) - 
15*b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 15*a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 
24*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 24*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 
35*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) - 35*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) - 
48*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) - 8*b1*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
30*b1*c2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 30*a1*d2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) + 
8*b1*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 30*b1*c2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 
30*a1*d2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 72*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
72*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 140*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 
140*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 240*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) - 
72*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 72*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 
140*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) - 140*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) - 
240*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) - 30*b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 
30*a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 144*b1*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 
144*a1*c2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 30*b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 
30*a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 144*b1*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 
144*a1*c2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 420*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 
420*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 960*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) - 
420*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) - 420*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) - 
960*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) + 
d1*(-(d2*(k2 + 2*RSP)*pow(M_E,k1*pow(RSP,-1))) + d2*(k1 + 2*RSP)*pow(M_E,k2*pow(RSP,-1)) + 
4*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,2) + b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + 
5*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) - 
4*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,2) - b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - 
5*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 
8*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) + 8*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 
15*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 15*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 
c2*pow(M_E,k2*pow(RSP,-1))*(3*k1*RSP + pow(k1,2) + 3*pow(RSP,2)) - 
c2*pow(M_E,k1*pow(RSP,-1))*(3*k2*RSP + pow(k2,2) + 3*pow(RSP,2)) - 8*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
30*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) + 8*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 
30*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 30*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 
30*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4)) - 144*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 
144*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 840*a2*b1*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 
840*a1*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 144*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 
144*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 840*a2*b1*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 
840*a1*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 2880*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) - 
2880*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) + 
c1*(5*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + 
6*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) - 
5*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 
6*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) + 
15*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 24*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 
15*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 24*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 
d2*pow(M_E,k2*pow(RSP,-1))*(3*k1*RSP + pow(k1,2) + 3*pow(RSP,2)) - 
d2*pow(M_E,k1*pow(RSP,-1))*(3*k2*RSP + pow(k2,2) + 3*pow(RSP,2)) - 30*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) + 
30*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 72*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) - 
72*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 
c2*pow(M_E,k2*pow(RSP,-1))*(4*RSP*pow(k1,2) + pow(k1,3) + 8*k1*pow(RSP,2) + 8*pow(RSP,3)) - 
c2*pow(M_E,k1*pow(RSP,-1))*(4*RSP*pow(k2,2) + pow(k2,3) + 8*k2*pow(RSP,2) + 8*pow(RSP,3)) - 
30*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 144*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 
30*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 144*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 
144*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 144*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5)) - 
840*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) - 840*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) - 
5760*a1*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) + 840*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) + 
840*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) + 5760*a1*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) - 
5760*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,7) + 5760*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,7))) - 
RSP*(-1 + pow(M_E,2*d*pow(RSP,-1)))*(7*b1*d2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + 
b1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + 8*b1*c2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + 
8*a1*d2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + 
a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + 9*b1*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + 
9*a1*c2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) + 
a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) + 10*a2*b1*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) + 
10*a1*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) + a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,7) + 
a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,7) + 11*a1*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,7) + 
a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,8) - 7*b1*d2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - 
b1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 8*b1*c2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 
8*a1*d2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - 
a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - 9*b1*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - 
9*a1*c2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) - 
a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) - 10*a2*b1*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) - 
10*a1*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) - a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,7) - 
a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,7) - 11*a1*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,7) - 
a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,8) + 24*b1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 
35*b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 35*a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 
48*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 48*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 
63*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) + 63*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) + 
80*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,6)*pow(RSP,2) - 24*b1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 
35*b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 35*a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 
48*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) - 48*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) - 
63*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) - 63*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) - 
80*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,6)*pow(RSP,2) - 48*b1*d2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) + 
48*b1*d2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 105*b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
105*a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 192*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 
192*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 315*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) + 
315*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) + 480*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5)*pow(RSP,3) - 
105*b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 105*a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 
192*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) - 192*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) - 
315*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) - 315*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) - 
480*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5)*pow(RSP,3) - 48*b1*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 
210*b1*c2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 210*a1*d2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 
48*b1*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 210*b1*c2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 
210*a1*d2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 576*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 
576*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 1260*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) + 
1260*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) + 2400*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,4) - 
576*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) - 576*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) - 
1260*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) - 1260*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) - 
2400*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,4) - 210*b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 
210*a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 1152*b1*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 
1152*a1*c2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 210*b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 
210*a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 1152*b1*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 
1152*a1*c2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 3780*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) + 
3780*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) + 9600*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,5) - 
3780*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) - 3780*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) - 
9600*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,5) + 
d1*(7*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + 
8*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) - 
7*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 
8*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) + 
24*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 35*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 
24*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 35*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 
d2*pow(M_E,k2*pow(RSP,-1))*(5*k1*RSP + pow(k1,2) + 8*pow(RSP,2)) - 
d2*pow(M_E,k1*pow(RSP,-1))*(5*k2*RSP + pow(k2,2) + 8*pow(RSP,2)) - 48*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) + 
48*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 105*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) - 
105*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 
c2*pow(M_E,k2*pow(RSP,-1))*(6*RSP*pow(k1,2) + pow(k1,3) + 15*k1*pow(RSP,2) + 15*pow(RSP,3)) - 
c2*pow(M_E,k1*pow(RSP,-1))*(6*RSP*pow(k2,2) + pow(k2,3) + 15*k2*pow(RSP,2) + 15*pow(RSP,3)) - 
48*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 210*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 
48*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 210*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 
210*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 210*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5)) - 
1152*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) - 1152*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) - 
7560*a2*b1*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) - 7560*a1*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) + 
1152*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) + 1152*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) + 
7560*a2*b1*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) + 7560*a1*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) + 
28800*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,6) - 28800*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,6) + 
c1*(8*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + 
9*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) - 
8*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - 
9*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) + 
35*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 48*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) - 
35*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 48*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) + 
105*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 192*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) - 
105*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 192*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) + 
d2*pow(M_E,k2*pow(RSP,-1))*(6*RSP*pow(k1,2) + pow(k1,3) + 15*k1*pow(RSP,2) + 15*pow(RSP,3)) - 
d2*pow(M_E,k1*pow(RSP,-1))*(6*RSP*pow(k2,2) + pow(k2,3) + 15*k2*pow(RSP,2) + 15*pow(RSP,3)) - 
210*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 210*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 
576*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) - 576*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 
c2*pow(M_E,k2*pow(RSP,-1))*(7*RSP*pow(k1,3) + pow(k1,4) + 24*pow(k1,2)*pow(RSP,2) + 48*k1*pow(RSP,3) + 
48*pow(RSP,4)) - c2*pow(M_E,k1*pow(RSP,-1))*
(7*RSP*pow(k2,3) + pow(k2,4) + 24*pow(k2,2)*pow(RSP,2) + 48*k2*pow(RSP,3) + 48*pow(RSP,4)) - 
210*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 1152*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 
210*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 1152*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 
1152*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) + 1152*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6)) - 
7560*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(RSP,7) - 7560*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,7) - 
57600*a1*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,7) + 7560*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,7) + 
7560*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,7) + 57600*a1*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,7) - 
57600*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,8) + 57600*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,8)) + 
d*(1 + pow(M_E,2*d*pow(RSP,-1)))*(7*b1*d2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + b1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + 
8*b1*c2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + 8*a1*d2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + 
b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + 
9*b1*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + 9*a1*c2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + 
b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) + a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) + 
10*a2*b1*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) + 10*a1*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) + 
a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,7) + a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,7) + 
11*a1*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,7) + a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,8) - 
7*b1*d2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - b1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 
8*b1*c2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 8*a1*d2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 
b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - 
9*b1*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - 9*a1*c2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - 
b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) - a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) - 
10*a2*b1*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) - 10*a1*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) - 
a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,7) - a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,7) - 
11*a1*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,7) - a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,8) + 
24*b1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 35*b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 
35*a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 48*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 
48*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 63*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) + 
63*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) + 80*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,6)*pow(RSP,2) - 
24*b1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 35*b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 
35*a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 48*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) - 
48*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) - 63*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) - 
63*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) - 80*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,6)*pow(RSP,2) - 
48*b1*d2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) + 48*b1*d2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 
105*b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 105*a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
192*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 192*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 
315*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) + 315*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) + 
480*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5)*pow(RSP,3) - 105*b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 
105*a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 192*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) - 
192*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) - 315*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) - 
315*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) - 480*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5)*pow(RSP,3) - 
48*b1*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 210*b1*c2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 
210*a1*d2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 48*b1*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 
210*b1*c2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 210*a1*d2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 
576*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 576*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 
1260*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) + 1260*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) + 
2400*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,4) - 576*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) - 
576*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) - 1260*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) - 
1260*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) - 2400*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,4) - 
210*b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 210*a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 
1152*b1*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 1152*a1*c2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 
210*b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 210*a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 
1152*b1*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 1152*a1*c2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 
3780*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) + 3780*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) + 
9600*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,5) - 3780*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) - 
3780*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) - 9600*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,5) + 
d1*(7*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + 
8*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) - 
7*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 
8*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) + 
24*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 35*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 
24*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 35*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 
d2*pow(M_E,k2*pow(RSP,-1))*(5*k1*RSP + pow(k1,2) + 8*pow(RSP,2)) - 
d2*pow(M_E,k1*pow(RSP,-1))*(5*k2*RSP + pow(k2,2) + 8*pow(RSP,2)) - 48*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) + 
48*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 105*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) - 
105*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 
c2*pow(M_E,k2*pow(RSP,-1))*(6*RSP*pow(k1,2) + pow(k1,3) + 15*k1*pow(RSP,2) + 15*pow(RSP,3)) - 
c2*pow(M_E,k1*pow(RSP,-1))*(6*RSP*pow(k2,2) + pow(k2,3) + 15*k2*pow(RSP,2) + 15*pow(RSP,3)) - 
48*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 210*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 
48*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 210*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 
210*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 210*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5)) - 
1152*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) - 1152*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) - 
7560*a2*b1*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) - 7560*a1*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) + 
1152*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) + 1152*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) + 
7560*a2*b1*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) + 7560*a1*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) + 
28800*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,6) - 28800*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,6) + 
c1*(8*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + 
9*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) - 
8*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - 
9*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) + 
35*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 48*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) - 
35*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 48*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) + 
105*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 192*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) - 
105*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 192*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) + 
d2*pow(M_E,k2*pow(RSP,-1))*(6*RSP*pow(k1,2) + pow(k1,3) + 15*k1*pow(RSP,2) + 15*pow(RSP,3)) - 
d2*pow(M_E,k1*pow(RSP,-1))*(6*RSP*pow(k2,2) + pow(k2,3) + 15*k2*pow(RSP,2) + 15*pow(RSP,3)) - 
210*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 210*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 
576*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) - 576*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 
c2*pow(M_E,k2*pow(RSP,-1))*(7*RSP*pow(k1,3) + pow(k1,4) + 24*pow(k1,2)*pow(RSP,2) + 48*k1*pow(RSP,3) + 
48*pow(RSP,4)) - c2*pow(M_E,k1*pow(RSP,-1))*
(7*RSP*pow(k2,3) + pow(k2,4) + 24*pow(k2,2)*pow(RSP,2) + 48*k2*pow(RSP,3) + 48*pow(RSP,4)) - 
210*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 1152*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 
210*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 1152*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 
1152*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) + 1152*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6)) - 
7560*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(RSP,7) - 7560*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,7) - 
57600*a1*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,7) + 7560*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,7) + 
7560*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,7) + 57600*a1*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,7) - 
57600*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,8) + 57600*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,8)))*pow(2*pow(d,2),-1);

    return res;
}




/* BM SH 2233 XX YY */

//  1
 double BMSH_Integral_2233_case_1( double AP, double RP, double k1, double k2, double a1, double b1, double c1, double d1, double d )
{
double res = 3*AP*pow(M_E,-((d + k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(-20*c1*d1*k2*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1)) - 2*c1*d1*k2*pow(d,3)*pow(M_E,k1*pow(RP,-1)) - 
4*c1*d1*RP*pow(d,3)*pow(M_E,k1*pow(RP,-1)) - 6*b1*d1*k2*RP*pow(d,3)*pow(M_E,k1*pow(RP,-1)) - 
3*k2*RP*pow(c1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1)) - 8*d*k2*RP*pow(d1,2)*pow(M_E,k1*pow(RP,-1)) - 
2*k2*pow(d,2)*pow(d1,2)*pow(M_E,k1*pow(RP,-1)) - 8*RP*pow(d,2)*pow(d1,2)*pow(M_E,k1*pow(RP,-1)) - 
pow(d,3)*pow(d1,2)*pow(M_E,k1*pow(RP,-1)) + 20*c1*d1*k1*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1)) + 
2*c1*d1*k1*pow(d,3)*pow(M_E,k2*pow(RP,-1)) + 4*c1*d1*RP*pow(d,3)*pow(M_E,k2*pow(RP,-1)) + 
6*b1*d1*k1*RP*pow(d,3)*pow(M_E,k2*pow(RP,-1)) + 3*k1*RP*pow(c1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1)) + 
8*d*k1*RP*pow(d1,2)*pow(M_E,k2*pow(RP,-1)) + 2*k1*pow(d,2)*pow(d1,2)*pow(M_E,k2*pow(RP,-1)) + 
8*RP*pow(d,2)*pow(d1,2)*pow(M_E,k2*pow(RP,-1)) + pow(d,3)*pow(d1,2)*pow(M_E,k2*pow(RP,-1)) + 
6*RP*gsl_sf_expint_Ei(-(k2*pow(RP,-1)))*pow(d,2)*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1)) - 
6*RP*gsl_sf_expint_Ei(k2*pow(RP,-1))*pow(d,2)*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1)) + 
gsl_sf_expint_Ei(-(k2*pow(RP,-1)))*pow(d,3)*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1)) - 
gsl_sf_expint_Ei(k2*pow(RP,-1))*pow(d,3)*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1)) - 
20*c1*d1*k1*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) - 2*c1*d1*k1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) + 
4*c1*d1*RP*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) + 6*b1*d1*k1*RP*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) + 
3*k1*RP*pow(c1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) + 8*d*k1*RP*pow(d1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) + 
2*k1*pow(d,2)*pow(d1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) - 8*RP*pow(d,2)*pow(d1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) - 
pow(d,3)*pow(d1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) + 20*c1*d1*k2*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) + 
2*c1*d1*k2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) - 4*c1*d1*RP*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) - 
6*b1*d1*k2*RP*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) - 3*k2*RP*pow(c1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) - 
8*d*k2*RP*pow(d1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) - 2*k2*pow(d,2)*pow(d1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) + 
8*RP*pow(d,2)*pow(d1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) + pow(d,3)*pow(d1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) + 
18*c1*d*d1*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 4*c1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 
24*b1*d1*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 12*RP*pow(c1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 
2*b1*d1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 8*b1*c1*RP*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 
8*a1*d1*RP*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + pow(c1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 
d*pow(d1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + RP*pow(d1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 
18*c1*d*d1*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) + 4*c1*d1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) - 
24*b1*d1*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) - 
12*RP*pow(c1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) - 
2*b1*d1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) + 
8*b1*c1*RP*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) + 
8*a1*d1*RP*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) - 
pow(c1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) - d*pow(d1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) - 
RP*pow(d1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) + 2*c1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 
2*c1*d1*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 20*b1*d*d1*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 
10*d*RP*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 4*b1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 
28*b1*c1*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 28*a1*d1*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 
2*pow(c1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 2*b1*c1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 
2*a1*d1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 10*a1*c1*RP*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 
5*RP*pow(b1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) - 2*c1*d*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) - 
2*c1*d1*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + 20*b1*d*d1*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + 
10*d*RP*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + 4*b1*d1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) - 
28*b1*c1*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) - 
28*a1*d1*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + 
2*pow(c1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) - 
2*b1*c1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) - 2*a1*d1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + 
10*a1*c1*RP*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + 
5*RP*pow(b1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + 2*b1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 
22*b1*c1*d*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 2*b1*d1*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 
22*a1*d*d1*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + d*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 
RP*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 4*b1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 
4*a1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 32*a1*c1*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 
16*RP*pow(b1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 2*a1*c1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 
12*a1*b1*RP*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + pow(b1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) - 
2*b1*d*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) + 22*b1*c1*d*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) - 
2*b1*d1*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) + 22*a1*d*d1*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) - 
d*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) - RP*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) + 
4*b1*c1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) + 4*a1*d1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) - 
32*a1*c1*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) - 
16*RP*pow(b1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) - 
2*a1*c1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) + 
12*a1*b1*RP*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) - 
pow(b1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) + 2*b1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 
2*a1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 2*b1*c1*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 
24*a1*c1*d*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 2*a1*d1*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 
12*d*RP*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 4*a1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 
36*a1*b1*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 2*pow(b1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 
2*a1*b1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 7*RP*pow(a1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) - 
2*b1*c1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) - 2*a1*d*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) - 
2*b1*c1*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) + 24*a1*c1*d*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) - 
2*a1*d1*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) + 12*d*RP*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) + 
4*a1*c1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) - 
36*a1*b1*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) + 
2*pow(b1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) - 
2*a1*b1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) + 
7*RP*pow(a1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) + 2*a1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,6) + 
2*a1*c1*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,6) + 26*a1*b1*d*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,6) + 
d*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,6) + RP*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,6) + 
4*a1*b1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,6) + 20*RP*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,6) + 
pow(a1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,6) - 2*a1*c1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6) - 
2*a1*c1*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6) + 26*a1*b1*d*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6) - 
d*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6) - RP*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6) + 
4*a1*b1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6) - 
20*RP*pow(a1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6) - 
pow(a1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6) + 2*a1*b1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,7) + 
2*a1*b1*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,7) + 14*d*RP*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,7) + 
2*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,7) - 2*a1*b1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,7) - 
2*a1*b1*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,7) + 14*d*RP*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,7) + 
2*pow(a1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,7) + d*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,8) + 
RP*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,8) - d*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,8) - 
RP*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,8) - 18*c1*d*d1*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 
4*c1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 24*b1*d1*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 
12*RP*pow(c1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 2*b1*d1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 
8*b1*c1*RP*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 8*a1*d1*RP*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 
pow(c1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - d*pow(d1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 
RP*pow(d1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 18*c1*d*d1*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) - 
4*c1*d1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) + 
24*b1*d1*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) + 
12*RP*pow(c1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) + 
2*b1*d1*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) - 
8*b1*c1*RP*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) - 
8*a1*d1*RP*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) + 
pow(c1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) + d*pow(d1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) + 
RP*pow(d1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) - 2*c1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 
2*c1*d1*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 20*b1*d*d1*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 
10*d*RP*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 4*b1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 
28*b1*c1*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 28*a1*d1*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 
2*pow(c1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 2*b1*c1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 
2*a1*d1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 10*a1*c1*RP*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 
5*RP*pow(b1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) + 2*c1*d*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) + 
2*c1*d1*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - 20*b1*d*d1*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - 
10*d*RP*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - 4*b1*d1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) + 
28*b1*c1*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) + 
28*a1*d1*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - 
2*pow(c1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) + 
2*b1*c1*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) + 2*a1*d1*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - 
10*a1*c1*RP*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - 
5*RP*pow(b1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - 2*b1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 
22*b1*c1*d*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 2*b1*d1*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 
22*a1*d*d1*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - d*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 
RP*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 4*b1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 
4*a1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 32*a1*c1*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 
16*RP*pow(b1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 2*a1*c1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 
12*a1*b1*RP*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - pow(b1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) + 
2*b1*d*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) - 22*b1*c1*d*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) + 
2*b1*d1*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) - 22*a1*d*d1*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) + 
d*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) + RP*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) - 
4*b1*c1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) - 4*a1*d1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) + 
32*a1*c1*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) + 
16*RP*pow(b1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) + 
2*a1*c1*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) - 
12*a1*b1*RP*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) + 
pow(b1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) - 2*b1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 
2*a1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 2*b1*c1*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 
24*a1*c1*d*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 2*a1*d1*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 
12*d*RP*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 4*a1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 
36*a1*b1*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 2*pow(b1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 
2*a1*b1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 7*RP*pow(a1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) + 
2*b1*c1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) + 2*a1*d*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) + 
2*b1*c1*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) - 24*a1*c1*d*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) + 
2*a1*d1*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) - 12*d*RP*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) - 
4*a1*c1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) + 
36*a1*b1*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) - 
2*pow(b1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) + 
2*a1*b1*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) - 
7*RP*pow(a1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) - 2*a1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,6) - 
2*a1*c1*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,6) - 26*a1*b1*d*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,6) - 
d*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,6) - RP*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,6) - 
4*a1*b1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,6) - 20*RP*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,6) - 
pow(a1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,6) + 2*a1*c1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6) + 
2*a1*c1*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6) - 26*a1*b1*d*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6) + 
d*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6) + RP*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6) - 
4*a1*b1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6) + 
20*RP*pow(a1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6) + 
pow(a1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6) - 2*a1*b1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,7) - 
2*a1*b1*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,7) - 14*d*RP*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,7) - 
2*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,7) + 2*a1*b1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,7) + 
2*a1*b1*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,7) - 14*d*RP*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,7) - 
2*pow(a1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,7) - d*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,8) - 
RP*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,8) + d*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,8) + 
RP*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,8) - 66*c1*d*d1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 
32*c1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 60*b1*d1*k2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 
30*k2*pow(c1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 6*b1*d1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 
16*b1*c1*k2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 16*a1*d1*k2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 
3*pow(c1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 23*d*pow(d1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 
8*k2*pow(d1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 66*c1*d*d1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 
32*c1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 60*b1*d1*k1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 
30*k1*pow(c1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 6*b1*d1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 
16*b1*c1*k1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 16*a1*d1*k1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 
3*pow(c1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 23*d*pow(d1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 
8*k1*pow(d1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 
15*d*gsl_sf_expint_Ei(-(k2*pow(RP,-1)))*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*pow(RP,2) - 
15*d*gsl_sf_expint_Ei(k2*pow(RP,-1))*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*pow(RP,2) - 
66*c1*d*d1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) + 32*c1*d1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) + 
60*b1*d1*k1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) + 
30*k1*pow(c1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) - 
6*b1*d1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) - 
16*b1*c1*k1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) - 
16*a1*d1*k1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) - 
3*pow(c1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) - 
23*d*pow(d1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) + 8*k1*pow(d1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) + 
66*c1*d*d1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) - 32*c1*d1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) - 
60*b1*d1*k2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) - 
30*k2*pow(c1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) + 
6*b1*d1*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) + 
16*b1*c1*k2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) + 
16*a1*d1*k2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) + 
3*pow(c1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) + 
23*d*pow(d1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) - 8*k2*pow(d1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) + 
18*c1*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 90*b1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
45*d*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 96*b1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
96*a1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
30*a1*c1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
15*pow(b1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
18*c1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
90*b1*d*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
45*d*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
96*b1*c1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
96*a1*d1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
30*a1*c1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
15*pow(b1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
118*b1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 20*b1*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
118*a1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 10*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
140*a1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
70*pow(b1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
48*a1*b1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
118*b1*c1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
20*b1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
118*a1*d*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
10*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
140*a1*c1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
70*pow(b1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
48*a1*b1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
22*b1*c1*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 150*a1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 
22*a1*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 75*d*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 
192*a1*b1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 
35*pow(a1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 
22*b1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
150*a1*c1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 
22*a1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
75*d*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 
192*a1*b1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
35*pow(a1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 
24*a1*c1*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 186*a1*b1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 
12*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 
126*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 
24*a1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 
186*a1*b1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 
12*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 
126*pow(a1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 
26*a1*b1*pow(M_E,k2*pow(RP,-1))*pow(k1,6)*pow(RP,2) + 113*d*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,6)*pow(RP,2) + 
26*a1*b1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6)*pow(RP,2) - 
113*d*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6)*pow(RP,2) + 
14*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,7)*pow(RP,2) + 
14*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,7)*pow(RP,2) - 18*c1*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
90*b1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 45*d*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
96*b1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
96*a1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
30*a1*c1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
15*pow(b1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
18*c1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
90*b1*d*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
45*d*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
96*b1*c1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
96*a1*d1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
30*a1*c1*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
15*pow(b1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
118*b1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 20*b1*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
118*a1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 10*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
140*a1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
70*pow(b1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
48*a1*b1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
118*b1*c1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
20*b1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
118*a1*d*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
10*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
140*a1*c1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
70*pow(b1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
48*a1*b1*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
22*b1*c1*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 150*a1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 
22*a1*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 75*d*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 
192*a1*b1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 
35*pow(a1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 
22*b1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
150*a1*c1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 
22*a1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
75*d*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 
192*a1*b1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
35*pow(a1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 
24*a1*c1*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 186*a1*b1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 
12*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 
126*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 
24*a1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 
186*a1*b1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 
12*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 
126*pow(a1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 
26*a1*b1*pow(M_E,k1*pow(RP,-1))*pow(k2,6)*pow(RP,2) - 113*d*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,6)*pow(RP,2) - 
26*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6)*pow(RP,2) + 
113*d*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6)*pow(RP,2) - 
14*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,7)*pow(RP,2) - 
14*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,7)*pow(RP,2) - 96*c1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 
66*c1*d1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 210*b1*d*d1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 
105*d*k2*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 60*b1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 
192*b1*c1*k2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 192*a1*d1*k2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 
30*pow(c1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 16*b1*c1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 
16*a1*d1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 60*a1*c1*k2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 
30*k2*pow(b1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 23*pow(d1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 
96*c1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 66*c1*d1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 
210*b1*d*d1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 105*d*k1*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 
60*b1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 192*b1*c1*k1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 
192*a1*d1*k1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 30*pow(c1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 
16*b1*c1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 16*a1*d1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 
60*a1*c1*k1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 30*k1*pow(b1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 
23*pow(d1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 
15*gsl_sf_expint_Ei(-(k2*pow(RP,-1)))*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*pow(RP,3) - 
15*gsl_sf_expint_Ei(k2*pow(RP,-1))*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*pow(RP,3) + 
96*c1*d*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 66*c1*d1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) + 
210*b1*d*d1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) + 105*d*k1*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 
60*b1*d1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 
192*b1*c1*k1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 
192*a1*d1*k1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 
30*pow(c1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) + 
16*b1*c1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) + 
16*a1*d1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) + 
60*a1*c1*k1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) + 
30*k1*pow(b1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 
23*pow(d1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 96*c1*d*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 
66*c1*d1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) - 210*b1*d*d1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) - 
105*d*k2*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 
60*b1*d1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 
192*b1*c1*k2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 
192*a1*d1*k2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 
30*pow(c1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) - 
16*b1*c1*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) - 
16*a1*d1*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) - 
60*a1*c1*k2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) - 
30*k2*pow(b1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 
23*pow(d1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 384*b1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
90*b1*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 384*a1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
45*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 420*a1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
210*pow(b1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
144*a1*b1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
384*b1*c1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
90*b1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
384*a1*d*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
45*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
420*a1*c1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
210*pow(b1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
144*a1*b1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
118*b1*c1*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 630*a1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
118*a1*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 315*d*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
768*a1*b1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
140*pow(a1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
118*b1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
630*a1*c1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
118*a1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
315*d*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
768*a1*b1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
140*pow(a1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
150*a1*c1*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 960*a1*b1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 
75*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 
630*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 
150*a1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 
960*a1*b1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 
75*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 
630*pow(a1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 
186*a1*b1*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,3) + 693*d*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,3) - 
186*a1*b1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,3) + 
693*d*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,3) + 
113*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,6)*pow(RP,3) - 
113*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6)*pow(RP,3) - 
384*b1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 90*b1*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
384*a1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 45*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
420*a1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
210*pow(b1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
144*a1*b1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
384*b1*c1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
90*b1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
384*a1*d*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
45*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
420*a1*c1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
210*pow(b1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
144*a1*b1*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
118*b1*c1*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 630*a1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 
118*a1*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 315*d*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 
768*a1*b1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 
140*pow(a1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
118*b1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 
630*a1*c1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
118*a1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 
315*d*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
768*a1*b1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 
140*pow(a1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 
150*a1*c1*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 960*a1*b1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 
75*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 
630*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 
150*a1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 
960*a1*b1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 
75*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 
630*pow(a1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 
186*a1*b1*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 693*d*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,3) + 
186*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 
693*d*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 
113*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,6)*pow(RP,3) + 
113*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6)*pow(RP,3) - 
gsl_sf_expint_Ei(-(k1*pow(RP,-1)))*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*
(6*RP*pow(d,2) + pow(d,3) + 15*d*pow(RP,2) + 15*pow(RP,3)) + 
gsl_sf_expint_Ei(k1*pow(RP,-1))*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*
(6*RP*pow(d,2) + pow(d,3) + 15*d*pow(RP,2) + 15*pow(RP,3)) - 96*c1*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 
210*b1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 768*b1*c1*d*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 
210*b1*d1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 768*a1*d*d1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 
105*d*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 105*k2*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 
192*b1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 192*a1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 
840*a1*c1*k2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 420*k2*pow(b1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 
60*a1*c1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 288*a1*b1*k2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 
30*pow(b1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 96*c1*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
210*b1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 768*b1*c1*d*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
210*b1*d1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 768*a1*d*d1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
105*d*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 105*k1*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
192*b1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 192*a1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
840*a1*c1*k1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 420*k1*pow(b1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
60*a1*c1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 288*a1*b1*k1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
30*pow(b1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 96*c1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 
210*b1*d*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 768*b1*c1*d*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 
210*b1*d1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 768*a1*d*d1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 
105*d*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 105*k1*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 
192*b1*c1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 
192*a1*d1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 
840*a1*c1*k1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 
420*k1*pow(b1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 
60*a1*c1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 
288*a1*b1*k1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 
30*pow(b1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 96*c1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 
210*b1*d*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 768*b1*c1*d*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 
210*b1*d1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 768*a1*d*d1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 
105*d*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 105*k2*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 
192*b1*c1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 
192*a1*d1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 
840*a1*c1*k2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 
420*k2*pow(b1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 
60*a1*c1*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 
288*a1*b1*k2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 
30*pow(b1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 384*b1*c1*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
1890*a1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 384*a1*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
945*d*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
2304*a1*b1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
420*pow(a1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
384*b1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
1890*a1*c1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
384*a1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
945*d*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
2304*a1*b1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
420*pow(a1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
630*a1*c1*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 3840*a1*b1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 
315*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 
2520*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 
630*a1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 
3840*a1*b1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 
315*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 
2520*pow(a1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 
960*a1*b1*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,4) + 3465*d*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,4) + 
960*a1*b1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,4) - 
3465*d*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,4) + 
693*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,4) + 
693*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,4) - 384*b1*c1*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
1890*a1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 384*a1*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
945*d*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
2304*a1*b1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
420*pow(a1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
384*b1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 
1890*a1*c1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
384*a1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 
945*d*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
2304*a1*b1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 
420*pow(a1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
630*a1*c1*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 3840*a1*b1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 
315*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 
2520*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 
630*a1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 
3840*a1*b1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 
315*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 
2520*pow(a1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 
960*a1*b1*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,4) - 3465*d*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,4) - 
960*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,4) + 
3465*d*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,4) - 
693*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,4) - 
693*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,4) - 768*b1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 
210*b1*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 768*a1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 
768*b1*c1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 3780*a1*c1*d*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 
768*a1*d1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 1890*d*k2*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 
105*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 840*a1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 
4608*a1*b1*k2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 420*pow(b1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 
288*a1*b1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 840*k2*pow(a1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 
768*b1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 210*b1*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
768*a1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 768*b1*c1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
3780*a1*c1*d*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 768*a1*d1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
1890*d*k1*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 105*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
840*a1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 4608*a1*b1*k1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
420*pow(b1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 288*a1*b1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
840*k1*pow(a1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 768*b1*c1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 
210*b1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) + 768*a1*d*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 
768*b1*c1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) + 3780*a1*c1*d*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 
768*a1*d1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) + 1890*d*k1*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 
105*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 840*a1*c1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 
4608*a1*b1*k1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 
420*pow(b1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) + 
288*a1*b1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) + 
840*k1*pow(a1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 
768*b1*c1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 210*b1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) - 
768*a1*d*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 768*b1*c1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) - 
3780*a1*c1*d*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 768*a1*d1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) - 
1890*d*k2*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 105*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 
840*a1*c1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 
4608*a1*b1*k2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 
420*pow(b1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) - 
288*a1*b1*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) - 
840*k2*pow(a1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 
1890*a1*c1*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 11520*a1*b1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 
945*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 
7560*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 
1890*a1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 
11520*a1*b1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 
945*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 
7560*pow(a1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 
3840*a1*b1*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,5) + 13860*d*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,5) - 
3840*a1*b1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,5) + 
13860*d*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,5) + 
3465*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,5) - 
3465*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,5) - 
1890*a1*c1*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 11520*a1*b1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 
945*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 
7560*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 
1890*a1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 
11520*a1*b1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 
945*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 
7560*pow(a1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 
3840*a1*b1*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 13860*d*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,5) + 
3840*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 
13860*d*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 
3465*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,5) + 
3465*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,5) - 768*b1*c1*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 
3780*a1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 768*a1*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 
3780*a1*c1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 23040*a1*b1*d*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 
1890*d*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 1890*k2*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 
4608*a1*b1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 15120*k2*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 
840*pow(a1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 768*b1*c1*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 
3780*a1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 768*a1*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 
3780*a1*c1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 23040*a1*b1*d*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 
1890*d*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 1890*k1*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 
4608*a1*b1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 15120*k1*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 
840*pow(a1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 768*b1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 
3780*a1*c1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) + 768*a1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) + 
3780*a1*c1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 23040*a1*b1*d*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 
1890*d*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) + 1890*k1*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) + 
4608*a1*b1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) + 
15120*k1*pow(a1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 
840*pow(a1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 768*b1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 
3780*a1*c1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) - 768*a1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) - 
3780*a1*c1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 23040*a1*b1*d*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 
1890*d*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) - 1890*k2*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) - 
4608*a1*b1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) - 
15120*k2*pow(a1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 
840*pow(a1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 
11520*a1*b1*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,6) + 41580*d*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,6) + 
11520*a1*b1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,6) - 
41580*d*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,6) + 
13860*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,6) + 
13860*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,6) - 
11520*a1*b1*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,6) - 41580*d*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,6) - 
11520*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,6) + 
41580*d*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,6) - 
13860*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,6) - 
13860*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,6) - 3780*a1*c1*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 
23040*a1*b1*d*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 23040*a1*b1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 
83160*d*k2*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 1890*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 
15120*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) + 3780*a1*c1*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 
23040*a1*b1*d*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 23040*a1*b1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 
83160*d*k1*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 1890*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 
15120*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) - 3780*a1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) + 
23040*a1*b1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) - 23040*a1*b1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) + 
83160*d*k1*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) - 
1890*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) - 
15120*pow(a1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) + 
3780*a1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) - 23040*a1*b1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) + 
23040*a1*b1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) - 
83160*d*k2*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) + 
1890*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) + 
15120*pow(a1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) + 
41580*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,7) - 
41580*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,7) - 
41580*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,7) + 
41580*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,7) - 23040*a1*b1*pow(M_E,k1*pow(RP,-1))*pow(RP,8) - 
83160*d*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,8) - 83160*k2*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,8) + 
23040*a1*b1*pow(M_E,k2*pow(RP,-1))*pow(RP,8) + 83160*d*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,8) + 
83160*k1*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,8) + 23040*a1*b1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,8) - 
83160*d*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,8) + 
83160*k1*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,8) - 23040*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,8) + 
83160*d*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,8) - 
83160*k2*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,8) - 83160*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,9) + 
83160*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,9) - 83160*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,9) + 
83160*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,9))*pow(2*pow(d,3),-1);

    return res;
}


// 2
double BMSH_Integral_2233_case_2_sub_1( double AP, double RP, double k1, double k2, double a1, double b1, double c1, double d1, double d )
{
double res = -3*AP*pow(M_E,-((2*d + k2)*pow(RP,-1)))*pow(RP,3)*(20*c1*d1*k2*RP*pow(d,2)*pow(M_E,d*pow(RP,-1)) + 
2*c1*d1*k2*pow(d,3)*pow(M_E,d*pow(RP,-1)) + 4*c1*d1*RP*pow(d,3)*pow(M_E,d*pow(RP,-1)) + 
6*b1*d1*k2*RP*pow(d,3)*pow(M_E,d*pow(RP,-1)) + 3*k2*RP*pow(c1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1)) + 
8*d*k2*RP*pow(d1,2)*pow(M_E,d*pow(RP,-1)) + 2*k2*pow(d,2)*pow(d1,2)*pow(M_E,d*pow(RP,-1)) + 
8*RP*pow(d,2)*pow(d1,2)*pow(M_E,d*pow(RP,-1)) + pow(d,3)*pow(d1,2)*pow(M_E,d*pow(RP,-1)) - 
20*c1*d1*k2*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1)) + 2*c1*d1*k2*pow(d,3)*pow(M_E,3*d*pow(RP,-1)) + 
4*c1*d1*RP*pow(d,3)*pow(M_E,3*d*pow(RP,-1)) + 6*b1*d1*k2*RP*pow(d,3)*pow(M_E,3*d*pow(RP,-1)) + 
3*k2*RP*pow(c1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1)) + 8*d*k2*RP*pow(d1,2)*pow(M_E,3*d*pow(RP,-1)) - 
2*k2*pow(d,2)*pow(d1,2)*pow(M_E,3*d*pow(RP,-1)) - 8*RP*pow(d,2)*pow(d1,2)*pow(M_E,3*d*pow(RP,-1)) + 
pow(d,3)*pow(d1,2)*pow(M_E,3*d*pow(RP,-1)) - 44*c1*d1*RP*pow(d,3)*pow(M_E,k2*pow(RP,-1)) - 
8*c1*d1*pow(d,4)*pow(M_E,k2*pow(RP,-1)) - 52*b1*d1*RP*pow(d,4)*pow(M_E,k2*pow(RP,-1)) - 
26*RP*pow(c1,2)*pow(d,4)*pow(M_E,k2*pow(RP,-1)) - 8*b1*d1*pow(d,5)*pow(M_E,k2*pow(RP,-1)) - 
60*b1*c1*RP*pow(d,5)*pow(M_E,k2*pow(RP,-1)) - 60*a1*d1*RP*pow(d,5)*pow(M_E,k2*pow(RP,-1)) - 
4*pow(c1,2)*pow(d,5)*pow(M_E,k2*pow(RP,-1)) - 8*b1*c1*pow(d,6)*pow(M_E,k2*pow(RP,-1)) - 
8*a1*d1*pow(d,6)*pow(M_E,k2*pow(RP,-1)) - 68*a1*c1*RP*pow(d,6)*pow(M_E,k2*pow(RP,-1)) - 
34*RP*pow(b1,2)*pow(d,6)*pow(M_E,k2*pow(RP,-1)) - 8*a1*c1*pow(d,7)*pow(M_E,k2*pow(RP,-1)) - 
76*a1*b1*RP*pow(d,7)*pow(M_E,k2*pow(RP,-1)) - 4*pow(b1,2)*pow(d,7)*pow(M_E,k2*pow(RP,-1)) - 
8*a1*b1*pow(d,8)*pow(M_E,k2*pow(RP,-1)) - 42*RP*pow(a1,2)*pow(d,8)*pow(M_E,k2*pow(RP,-1)) - 
4*pow(a1,2)*pow(d,9)*pow(M_E,k2*pow(RP,-1)) - 17*RP*pow(d,2)*pow(d1,2)*pow(M_E,k2*pow(RP,-1)) - 
4*pow(d,3)*pow(d1,2)*pow(M_E,k2*pow(RP,-1)) + RP*pow(d,2)*pow(d1,2)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 
18*c1*d*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k2,2) + 4*c1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2) + 
24*b1*d1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2) + 12*RP*pow(c1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2) + 
2*b1*d1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2) + 8*b1*c1*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2) + 
8*a1*d1*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2) + pow(c1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2) + 
d*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2) + RP*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2) + 
18*c1*d*d1*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) - 4*c1*d1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) - 
24*b1*d1*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) - 12*RP*pow(c1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) + 
2*b1*d1*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) + 8*b1*c1*RP*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) + 
8*a1*d1*RP*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) + pow(c1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) + 
d*pow(d1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) - RP*pow(d1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) + 
2*c1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 2*c1*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 
20*b1*d*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 10*d*RP*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 
4*b1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 28*b1*c1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 
28*a1*d1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 2*pow(c1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 
2*b1*c1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 2*a1*d1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 
10*a1*c1*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 5*RP*pow(b1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 
2*c1*d*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) - 2*c1*d1*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) + 
20*b1*d*d1*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) + 10*d*RP*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) - 
4*b1*d1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) - 28*b1*c1*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) - 
28*a1*d1*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) - 2*pow(c1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) + 
2*b1*c1*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) + 2*a1*d1*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) + 
10*a1*c1*RP*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) + 5*RP*pow(b1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) + 
2*b1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 22*b1*c1*d*RP*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 
2*b1*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 22*a1*d*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 
d*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + RP*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 
4*b1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 4*a1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 
32*a1*c1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 16*RP*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 
2*a1*c1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 12*a1*b1*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 
pow(b1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 2*b1*d*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) + 
22*b1*c1*d*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) - 2*b1*d1*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) + 
22*a1*d*d1*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) + d*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) - 
RP*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) - 4*b1*c1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) - 
4*a1*d1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) - 32*a1*c1*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) - 
16*RP*pow(b1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) + 2*a1*c1*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) + 
12*a1*b1*RP*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) + pow(b1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) + 
2*b1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 2*a1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 
2*b1*c1*RP*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 24*a1*c1*d*RP*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 
2*a1*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 12*d*RP*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 
4*a1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 36*a1*b1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 
2*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 2*a1*b1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 
7*RP*pow(a1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 2*b1*c1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) + 
2*a1*d*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) - 2*b1*c1*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) + 
24*a1*c1*d*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) - 2*a1*d1*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) + 
12*d*RP*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) - 4*a1*c1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) - 
36*a1*b1*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) - 2*pow(b1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) + 
2*a1*b1*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) + 7*RP*pow(a1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) + 
2*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k2,6) + 2*a1*c1*RP*pow(M_E,d*pow(RP,-1))*pow(k2,6) + 
26*a1*b1*d*RP*pow(M_E,d*pow(RP,-1))*pow(k2,6) + d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,6) + 
RP*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,6) + 4*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,6) + 
20*RP*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,6) + pow(a1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,6) + 
2*a1*c1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,6) - 2*a1*c1*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,6) + 
26*a1*b1*d*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,6) + d*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6) - 
RP*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6) - 4*a1*b1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6) - 
20*RP*pow(a1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6) + pow(a1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6) + 
2*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k2,7) + 2*a1*b1*RP*pow(M_E,d*pow(RP,-1))*pow(k2,7) + 
14*d*RP*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,7) + 2*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,7) + 
2*a1*b1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,7) - 2*a1*b1*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,7) + 
14*d*RP*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,7) - 2*pow(a1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,7) + 
d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,8) + RP*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,8) + 
d*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,8) - RP*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,8) + 
66*c1*d*d1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 32*c1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 
60*b1*d1*k2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 30*k2*pow(c1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 
6*b1*d1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 16*b1*c1*k2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 
16*a1*d1*k2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 3*pow(c1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 
23*d*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 8*k2*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 
66*c1*d*d1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) - 32*c1*d1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) - 
60*b1*d1*k2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) - 30*k2*pow(c1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) + 
6*b1*d1*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) + 16*b1*c1*k2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) + 
16*a1*d1*k2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) + 3*pow(c1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) + 
23*d*pow(d1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) - 8*k2*pow(d1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) - 
116*c1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 176*b1*d1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 
88*pow(c1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 252*b1*c1*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 
252*a1*d1*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 344*a1*c1*pow(d,5)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 
172*pow(b1,2)*pow(d,5)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 452*a1*b1*pow(d,6)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 
288*pow(a1,2)*pow(d,7)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 31*d*pow(d1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 
16*c1*d1*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) - 16*b1*d1*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) - 
8*pow(c1,2)*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) - 
16*b1*c1*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) - 16*a1*d1*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) - 
16*a1*c1*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) - 
8*pow(b1,2)*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) - 
16*a1*b1*pow(d,6)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) - 
8*pow(a1,2)*pow(d,7)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) - 15*d*pow(d1,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 
18*c1*d1*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 90*b1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
45*d*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 96*b1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
96*a1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
30*a1*c1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
15*pow(b1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 18*c1*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
90*b1*d*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 45*d*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
96*b1*c1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
96*a1*d1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
30*a1*c1*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
15*pow(b1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
118*b1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 20*b1*d1*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
118*a1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 10*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
140*a1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
70*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
48*a1*b1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 118*b1*c1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
20*b1*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 118*a1*d*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
10*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
140*a1*c1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
70*pow(b1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
48*a1*b1*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 22*b1*c1*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
150*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 22*a1*d1*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
75*d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 192*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
35*pow(a1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 22*b1*c1*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
150*a1*c1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 22*a1*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
75*d*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 
192*a1*b1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
35*pow(a1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 24*a1*c1*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 
186*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 12*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 
126*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 24*a1*c1*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 
186*a1*b1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 12*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 
126*pow(a1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 26*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k2,6)*pow(RP,2) + 
113*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,6)*pow(RP,2) - 26*a1*b1*pow(M_E,3*d*pow(RP,-1))*pow(k2,6)*pow(RP,2) + 
113*d*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6)*pow(RP,2) + 14*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,7)*pow(RP,2) - 
14*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,7)*pow(RP,2) + 96*c1*d*d1*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 
66*c1*d1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 210*b1*d*d1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 
105*d*k2*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 60*b1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 
192*b1*c1*k2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 192*a1*d1*k2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 
30*pow(c1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 16*b1*c1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 
16*a1*d1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 60*a1*c1*k2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 
30*k2*pow(b1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 23*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 
96*c1*d*d1*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 66*c1*d1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) + 
210*b1*d*d1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) + 105*d*k2*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 
60*b1*d1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 192*b1*c1*k2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 
192*a1*d1*k2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 30*pow(c1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) + 
16*b1*c1*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) + 16*a1*d1*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) + 
60*a1*c1*k2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) + 30*k2*pow(b1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 
23*pow(d1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 162*c1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 
360*b1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 180*pow(c1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 
710*b1*c1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 710*a1*d1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 
1260*a1*c1*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 630*pow(b1,2)*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 
2058*a1*b1*pow(d,5)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 1576*pow(a1,2)*pow(d,6)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 
23*pow(d1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 30*c1*d*d1*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) - 
60*b1*d1*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) - 
30*pow(c1,2)*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) - 
90*b1*c1*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) - 90*a1*d1*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) - 
120*a1*c1*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) - 
60*pow(b1,2)*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) - 
150*a1*b1*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) - 
90*pow(a1,2)*pow(d,6)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) + 23*pow(d1,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) + 
384*b1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 90*b1*d1*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
384*a1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 45*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
420*a1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
210*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
144*a1*b1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 384*b1*c1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
90*b1*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 384*a1*d*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
45*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
420*a1*c1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
210*pow(b1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
144*a1*b1*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 118*b1*c1*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
630*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 118*a1*d1*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
315*d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 768*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
140*pow(a1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 118*b1*c1*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
630*a1*c1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 118*a1*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
315*d*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 
768*a1*b1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
140*pow(a1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 150*a1*c1*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 
960*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 75*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 
630*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 150*a1*c1*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 
960*a1*b1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 75*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 
630*pow(a1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 186*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,3) + 
693*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 186*a1*b1*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,3) + 
693*d*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,3) + 113*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,6)*pow(RP,3) - 
113*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6)*pow(RP,3) + 
gsl_sf_expint_Ei(-(d*pow(RP,-1)))*pow(d1,2)*pow(M_E,(d + k2)*pow(RP,-1))*
(-6*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1))) + pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1))) + 
15*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,2) - 15*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,3)) - 
gsl_sf_expint_Ei(-(k2*pow(RP,-1)))*pow(d1,2)*pow(M_E,(d + k2)*pow(RP,-1))*
(-6*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1))) + pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1))) + 
15*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,2) - 15*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,3)) + 
96*c1*d1*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 210*b1*d*d1*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 
768*b1*c1*d*k2*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 210*b1*d1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 
768*a1*d*d1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 105*d*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 
105*k2*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 192*b1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 
192*a1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 840*a1*c1*k2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 
420*k2*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 60*a1*c1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 
288*a1*b1*k2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 30*pow(b1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
96*c1*d1*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 210*b1*d*d1*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 
768*b1*c1*d*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 210*b1*d1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 
768*a1*d*d1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 105*d*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 
105*k2*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 192*b1*c1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 
192*a1*d1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 840*a1*c1*k2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 
420*k2*pow(b1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 60*a1*c1*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 
288*a1*b1*k2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 30*pow(b1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 
96*c1*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 420*b1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 
210*d*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 1344*b1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 
1344*a1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 3420*a1*c1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 
1710*pow(b1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 7392*a1*b1*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 
7098*pow(a1,2)*pow(d,5)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 96*c1*d1*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) - 
192*b1*c1*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) - 
192*a1*d1*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) - 
480*a1*c1*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) - 
240*pow(b1,2)*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) - 
864*a1*b1*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) - 
672*pow(a1,2)*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) + 384*b1*c1*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 
1890*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 384*a1*d1*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 
945*d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 
2304*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 
420*pow(a1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 384*b1*c1*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 
1890*a1*c1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 384*a1*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 
945*d*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
2304*a1*b1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 
420*pow(a1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 630*a1*c1*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 
3840*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 315*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 
2520*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 
630*a1*c1*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 3840*a1*b1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 
315*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 
2520*pow(a1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 
960*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,4) + 3465*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,4) - 
960*a1*b1*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,4) + 3465*d*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,4) + 
693*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,4) - 693*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,4) + 
768*b1*c1*d*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 210*b1*d1*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 
768*a1*d*d1*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 768*b1*c1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 
3780*a1*c1*d*k2*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 768*a1*d1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 
1890*d*k2*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 105*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 
840*a1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 4608*a1*b1*k2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 
420*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 288*a1*b1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 
840*k2*pow(a1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 768*b1*c1*d*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 
210*b1*d1*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 768*a1*d*d1*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 
768*b1*c1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 3780*a1*c1*d*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 
768*a1*d1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 1890*d*k2*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 
105*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 840*a1*c1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 
4608*a1*b1*k2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 420*pow(b1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 
288*a1*b1*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 840*k2*pow(a1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 
1536*b1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 210*b1*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 
1536*a1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 105*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 
6510*a1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 3255*pow(b1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 
20256*a1*b1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 25725*pow(a1,2)*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
210*b1*d1*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) + 105*pow(c1,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) - 
1050*a1*c1*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) - 
525*pow(b1,2)*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) - 
3360*a1*b1*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) - 
3675*pow(a1,2)*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) + 1890*a1*c1*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 
11520*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 945*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 
7560*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 
1890*a1*c1*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 11520*a1*b1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 
945*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 
7560*pow(a1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 
3840*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,5) + 13860*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 
3840*a1*b1*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,5) + 13860*d*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,5) + 
3465*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,5) - 3465*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,5) + 
768*b1*c1*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 3780*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 
768*a1*d1*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 3780*a1*c1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 
23040*a1*b1*d*k2*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 1890*d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 
1890*k2*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 4608*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 
15120*k2*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 840*pow(a1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 
768*b1*c1*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) + 3780*a1*c1*d*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 
768*a1*d1*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 3780*a1*c1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) + 
23040*a1*b1*d*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) + 1890*d*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 
1890*k2*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 4608*a1*b1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 
15120*k2*pow(a1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) + 840*pow(a1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 
768*b1*c1*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 7560*a1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 
768*a1*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 3780*d*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 
39168*a1*b1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 71400*pow(a1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 
768*b1*c1*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,6) + 768*a1*d1*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,6) - 
6912*a1*b1*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,6) - 
13440*pow(a1,2)*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,6) + 
11520*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,6) + 41580*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,6) - 
11520*a1*b1*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,6) + 41580*d*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,6) + 
13860*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,6) - 13860*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,6) + 
3780*a1*c1*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 23040*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 
23040*a1*b1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 83160*d*k2*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 
1890*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 15120*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 
3780*a1*c1*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) + 23040*a1*b1*d*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) - 
23040*a1*b1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) + 83160*d*k2*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) - 
1890*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) - 15120*pow(a1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) - 
3780*a1*c1*pow(M_E,k2*pow(RP,-1))*pow(RP,7) - 46080*a1*b1*d*pow(M_E,k2*pow(RP,-1))*pow(RP,7) - 
1890*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) - 139860*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 
3780*a1*c1*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,7) + 1890*pow(b1,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,7) - 
26460*pow(a1,2)*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,7) + 
41580*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,7) - 41580*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,7) + 
23040*a1*b1*pow(M_E,d*pow(RP,-1))*pow(RP,8) + 83160*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,8) + 
83160*k2*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,8) - 23040*a1*b1*pow(M_E,3*d*pow(RP,-1))*pow(RP,8) + 
83160*d*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,8) - 83160*k2*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,8) - 
23040*a1*b1*pow(M_E,k2*pow(RP,-1))*pow(RP,8) - 166320*d*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,8) + 
23040*a1*b1*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,8) + 83160*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,9) - 
83160*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,9) - 83160*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,9) + 
83160*pow(a1,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,9))*pow(2*pow(d,3),-1);

    return res;
}


//  3
double BMSH_Integral_2233_case_2_sub_2( double AP, double RP, double k1, double k2, double a1, double b1, double c1, double d1, double d )
{
double res = -3*AP*pow(M_E,-((2*d + k1)*pow(RP,-1)))*pow(RP,3)*(-20*c1*d1*k1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1)) - 
2*c1*d1*k1*pow(d,3)*pow(M_E,d*pow(RP,-1)) - 4*c1*d1*RP*pow(d,3)*pow(M_E,d*pow(RP,-1)) - 
6*b1*d1*k1*RP*pow(d,3)*pow(M_E,d*pow(RP,-1)) - 3*k1*RP*pow(c1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1)) - 
8*d*k1*RP*pow(d1,2)*pow(M_E,d*pow(RP,-1)) - 2*k1*pow(d,2)*pow(d1,2)*pow(M_E,d*pow(RP,-1)) - 
8*RP*pow(d,2)*pow(d1,2)*pow(M_E,d*pow(RP,-1)) - pow(d,3)*pow(d1,2)*pow(M_E,d*pow(RP,-1)) + 
44*c1*d1*RP*pow(d,3)*pow(M_E,k1*pow(RP,-1)) + 8*c1*d1*pow(d,4)*pow(M_E,k1*pow(RP,-1)) + 
52*b1*d1*RP*pow(d,4)*pow(M_E,k1*pow(RP,-1)) + 26*RP*pow(c1,2)*pow(d,4)*pow(M_E,k1*pow(RP,-1)) + 
8*b1*d1*pow(d,5)*pow(M_E,k1*pow(RP,-1)) + 60*b1*c1*RP*pow(d,5)*pow(M_E,k1*pow(RP,-1)) + 
60*a1*d1*RP*pow(d,5)*pow(M_E,k1*pow(RP,-1)) + 4*pow(c1,2)*pow(d,5)*pow(M_E,k1*pow(RP,-1)) + 
8*b1*c1*pow(d,6)*pow(M_E,k1*pow(RP,-1)) + 8*a1*d1*pow(d,6)*pow(M_E,k1*pow(RP,-1)) + 
68*a1*c1*RP*pow(d,6)*pow(M_E,k1*pow(RP,-1)) + 34*RP*pow(b1,2)*pow(d,6)*pow(M_E,k1*pow(RP,-1)) + 
8*a1*c1*pow(d,7)*pow(M_E,k1*pow(RP,-1)) + 76*a1*b1*RP*pow(d,7)*pow(M_E,k1*pow(RP,-1)) + 
4*pow(b1,2)*pow(d,7)*pow(M_E,k1*pow(RP,-1)) + 8*a1*b1*pow(d,8)*pow(M_E,k1*pow(RP,-1)) + 
42*RP*pow(a1,2)*pow(d,8)*pow(M_E,k1*pow(RP,-1)) + 4*pow(a1,2)*pow(d,9)*pow(M_E,k1*pow(RP,-1)) + 
17*RP*pow(d,2)*pow(d1,2)*pow(M_E,k1*pow(RP,-1)) + 4*pow(d,3)*pow(d1,2)*pow(M_E,k1*pow(RP,-1)) + 
6*RP*gsl_sf_expint_Ei(-(k1*pow(RP,-1)))*pow(d,2)*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1)) - 
6*RP*gsl_sf_expint_Ei(k1*pow(RP,-1))*pow(d,2)*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1)) + 
gsl_sf_expint_Ei(-(k1*pow(RP,-1)))*pow(d,3)*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1)) - 
gsl_sf_expint_Ei(k1*pow(RP,-1))*pow(d,3)*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1)) - 
RP*pow(d,2)*pow(d1,2)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 20*c1*d1*k1*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1)) + 
2*c1*d1*k1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1)) - 4*c1*d1*RP*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1)) - 
6*b1*d1*k1*RP*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1)) - 3*k1*RP*pow(c1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1)) - 
8*d*k1*RP*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1)) - 2*k1*pow(d,2)*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1)) + 
8*RP*pow(d,2)*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1)) + pow(d,3)*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1)) - 
18*c1*d*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 4*c1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 
24*b1*d1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 12*RP*pow(c1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 
2*b1*d1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 8*b1*c1*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 
8*a1*d1*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - pow(c1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 
d*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - RP*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 
18*c1*d*d1*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) - 4*c1*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) + 
24*b1*d1*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) + 
12*RP*pow(c1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) + 
2*b1*d1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) - 8*b1*c1*RP*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) - 
8*a1*d1*RP*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) + 
pow(c1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) + d*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) + 
RP*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) - 2*c1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 
2*c1*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 20*b1*d*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 
10*d*RP*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 4*b1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 
28*b1*c1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 28*a1*d1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 
2*pow(c1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 2*b1*c1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 
2*a1*d1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 10*a1*c1*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 
5*RP*pow(b1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,3) + 2*c1*d*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) + 
2*c1*d1*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - 20*b1*d*d1*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - 
10*d*RP*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - 4*b1*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) + 
28*b1*c1*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) + 
28*a1*d1*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - 
2*pow(c1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) + 
2*b1*c1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) + 2*a1*d1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - 
10*a1*c1*RP*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - 
5*RP*pow(b1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - 2*b1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 
22*b1*c1*d*RP*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 2*b1*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 
22*a1*d*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k1,4) - d*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 
RP*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 4*b1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 
4*a1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 32*a1*c1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 
16*RP*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 2*a1*c1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 
12*a1*b1*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - pow(b1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,4) + 
2*b1*d*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) - 22*b1*c1*d*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) + 
2*b1*d1*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) - 22*a1*d*d1*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) + 
d*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) + RP*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) - 
4*b1*c1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) - 4*a1*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) + 
32*a1*c1*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) + 
16*RP*pow(b1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) + 
2*a1*c1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) - 
12*a1*b1*RP*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) + 
pow(b1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) - 2*b1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 
2*a1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 2*b1*c1*RP*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 
24*a1*c1*d*RP*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 2*a1*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 
12*d*RP*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 4*a1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 
36*a1*b1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 2*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 
2*a1*b1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 7*RP*pow(a1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,5) + 
2*b1*c1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) + 2*a1*d*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) + 
2*b1*c1*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) - 24*a1*c1*d*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) + 
2*a1*d1*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) - 12*d*RP*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) - 
4*a1*c1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) + 
36*a1*b1*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) - 
2*pow(b1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) + 
2*a1*b1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) - 
7*RP*pow(a1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) - 2*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k1,6) - 
2*a1*c1*RP*pow(M_E,d*pow(RP,-1))*pow(k1,6) - 26*a1*b1*d*RP*pow(M_E,d*pow(RP,-1))*pow(k1,6) - 
d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,6) - RP*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,6) - 
4*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,6) - 20*RP*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,6) - 
pow(a1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,6) + 2*a1*c1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6) + 
2*a1*c1*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6) - 26*a1*b1*d*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6) + 
d*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6) + RP*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6) - 
4*a1*b1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6) + 
20*RP*pow(a1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6) + 
pow(a1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6) - 2*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k1,7) - 
2*a1*b1*RP*pow(M_E,d*pow(RP,-1))*pow(k1,7) - 14*d*RP*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,7) - 
2*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,7) + 2*a1*b1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,7) + 
2*a1*b1*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,7) - 14*d*RP*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,7) - 
2*pow(a1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,7) - d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,8) - 
RP*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,8) + d*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,8) + 
RP*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,8) - 66*c1*d*d1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 
32*c1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 60*b1*d1*k1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 
30*k1*pow(c1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 6*b1*d1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 
16*b1*c1*k1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 16*a1*d1*k1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 
3*pow(c1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 23*d*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 
8*k1*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 116*c1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
176*b1*d1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 88*pow(c1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
252*b1*c1*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 252*a1*d1*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
344*a1*c1*pow(d,5)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 172*pow(b1,2)*pow(d,5)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
452*a1*b1*pow(d,6)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 288*pow(a1,2)*pow(d,7)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
31*d*pow(d1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
15*d*gsl_sf_expint_Ei(-(k1*pow(RP,-1)))*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1))*pow(RP,2) - 
15*d*gsl_sf_expint_Ei(k1*pow(RP,-1))*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1))*pow(RP,2) - 
16*c1*d1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) - 16*b1*d1*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) - 
8*pow(c1,2)*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) - 
16*b1*c1*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) - 16*a1*d1*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) - 
16*a1*c1*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) - 
8*pow(b1,2)*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) - 
16*a1*b1*pow(d,6)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) - 
8*pow(a1,2)*pow(d,7)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) - 15*d*pow(d1,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 
66*c1*d*d1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) - 32*c1*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) - 
60*b1*d1*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) - 
30*k1*pow(c1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) + 
6*b1*d1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) + 
16*b1*c1*k1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) + 
16*a1*d1*k1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) + 
3*pow(c1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) + 23*d*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) - 
8*k1*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) - 18*c1*d1*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
90*b1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 45*d*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
96*b1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
96*a1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
30*a1*c1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
15*pow(b1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
18*c1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
90*b1*d*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
45*d*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
96*b1*c1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
96*a1*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
30*a1*c1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
15*pow(b1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
118*b1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 20*b1*d1*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
118*a1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 10*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
140*a1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
70*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
48*a1*b1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
118*b1*c1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
20*b1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
118*a1*d*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
10*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
140*a1*c1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
70*pow(b1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
48*a1*b1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
22*b1*c1*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 150*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
22*a1*d1*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 75*d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
192*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
35*pow(a1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
22*b1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 
150*a1*c1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
22*a1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 
75*d*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
192*a1*b1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 
35*pow(a1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
24*a1*c1*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 186*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 
12*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 
126*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 
24*a1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 
186*a1*b1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 
12*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 
126*pow(a1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 
26*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k1,6)*pow(RP,2) - 113*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,6)*pow(RP,2) - 
26*a1*b1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6)*pow(RP,2) + 
113*d*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6)*pow(RP,2) - 
14*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,7)*pow(RP,2) - 14*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,7)*pow(RP,2) - 
96*c1*d*d1*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 66*c1*d1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 
210*b1*d*d1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 105*d*k1*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 
60*b1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 192*b1*c1*k1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 
192*a1*d1*k1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 30*pow(c1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 
16*b1*c1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 16*a1*d1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 
60*a1*c1*k1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 30*k1*pow(b1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 
23*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 162*c1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 
360*b1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 180*pow(c1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 
710*b1*c1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 710*a1*d1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 
1260*a1*c1*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 630*pow(b1,2)*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 
2058*a1*b1*pow(d,5)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 1576*pow(a1,2)*pow(d,6)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 
23*pow(d1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 
15*gsl_sf_expint_Ei(-(k1*pow(RP,-1)))*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1))*pow(RP,3) - 
15*gsl_sf_expint_Ei(k1*pow(RP,-1))*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1))*pow(RP,3) + 
30*c1*d*d1*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) + 60*b1*d1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) + 
30*pow(c1,2)*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) + 
90*b1*c1*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) + 90*a1*d1*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) + 
120*a1*c1*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) + 
60*pow(b1,2)*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) + 
150*a1*b1*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) + 
90*pow(a1,2)*pow(d,6)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 23*pow(d1,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 
96*c1*d*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 66*c1*d1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 
210*b1*d*d1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 105*d*k1*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 
60*b1*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 
192*b1*c1*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 
192*a1*d1*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 
30*pow(c1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 
16*b1*c1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 16*a1*d1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 
60*a1*c1*k1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 
30*k1*pow(b1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 
23*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 384*b1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
90*b1*d1*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 384*a1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
45*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 420*a1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
210*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
144*a1*b1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
384*b1*c1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
90*b1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
384*a1*d*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
45*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
420*a1*c1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
210*pow(b1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
144*a1*b1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
118*b1*c1*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 630*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
118*a1*d1*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 315*d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
768*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
140*pow(a1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
118*b1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
630*a1*c1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
118*a1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
315*d*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
768*a1*b1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
140*pow(a1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
150*a1*c1*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 960*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 
75*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 
630*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 
150*a1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 
960*a1*b1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 
75*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 
630*pow(a1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 
186*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,3) - 693*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,3) + 
186*a1*b1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,3) - 
693*d*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,3) - 
113*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,6)*pow(RP,3) + 
113*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6)*pow(RP,3) - 
gsl_sf_expint_Ei(-(d*pow(RP,-1)))*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1))*
(6*RP*pow(d,2) + pow(d,3) + 15*d*pow(RP,2) + 15*pow(RP,3)) + 
gsl_sf_expint_Ei(d*pow(RP,-1))*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1))*
(6*RP*pow(d,2) + pow(d,3) + 15*d*pow(RP,2) + 15*pow(RP,3)) - 96*c1*d1*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
210*b1*d*d1*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 768*b1*c1*d*k1*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
210*b1*d1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 768*a1*d*d1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
105*d*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 105*k1*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
192*b1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 192*a1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
840*a1*c1*k1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 420*k1*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
60*a1*c1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 288*a1*b1*k1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
30*pow(b1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 96*c1*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 
420*b1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 210*d*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 
1344*b1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 1344*a1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 
3420*a1*c1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 1710*pow(b1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 
7392*a1*b1*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 7098*pow(a1,2)*pow(d,5)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 
96*c1*d1*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) - 192*b1*c1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) - 
192*a1*d1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) - 
480*a1*c1*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) - 
240*pow(b1,2)*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) - 
864*a1*b1*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) - 
672*pow(a1,2)*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) - 96*c1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 
210*b1*d*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 768*b1*c1*d*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 
210*b1*d1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 768*a1*d*d1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 
105*d*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 105*k1*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 
192*b1*c1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 
192*a1*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 
840*a1*c1*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 
420*k1*pow(b1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 
60*a1*c1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 
288*a1*b1*k1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 
30*pow(b1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 384*b1*c1*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
1890*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 384*a1*d1*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
945*d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
2304*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
420*pow(a1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
384*b1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
1890*a1*c1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
384*a1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
945*d*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
2304*a1*b1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
420*pow(a1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
630*a1*c1*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 3840*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 
315*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 
2520*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 
630*a1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 
3840*a1*b1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 
315*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 
2520*pow(a1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 
960*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,4) - 3465*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,4) - 
960*a1*b1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,4) + 
3465*d*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,4) - 
693*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,4) - 
693*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,4) - 768*b1*c1*d*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 
210*b1*d1*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 768*a1*d*d1*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 
768*b1*c1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 3780*a1*c1*d*k1*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 
768*a1*d1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 1890*d*k1*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 
105*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 840*a1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 
4608*a1*b1*k1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 420*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 
288*a1*b1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 840*k1*pow(a1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 
1536*b1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 210*b1*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 
1536*a1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 105*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 
6510*a1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 3255*pow(b1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 
20256*a1*b1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 25725*pow(a1,2)*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 
210*b1*d1*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) - 105*pow(c1,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) + 
1050*a1*c1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) + 
525*pow(b1,2)*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) + 
3360*a1*b1*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) + 
3675*pow(a1,2)*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) - 768*b1*c1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 
210*b1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 768*a1*d*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 
768*b1*c1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 3780*a1*c1*d*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 
768*a1*d1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 1890*d*k1*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 
105*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 840*a1*c1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 
4608*a1*b1*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 
420*pow(b1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 
288*a1*b1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 
840*k1*pow(a1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 
1890*a1*c1*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 11520*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 
945*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 
7560*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 
1890*a1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 
11520*a1*b1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 
945*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 
7560*pow(a1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 
3840*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,5) - 13860*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,5) + 
3840*a1*b1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,5) - 
13860*d*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,5) - 
3465*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,5) + 
3465*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,5) - 768*b1*c1*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 
3780*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 768*a1*d1*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 
3780*a1*c1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 23040*a1*b1*d*k1*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 
1890*d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 1890*k1*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 
4608*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 15120*k1*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 
840*pow(a1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 768*b1*c1*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 
7560*a1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 768*a1*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 
3780*d*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 39168*a1*b1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 
71400*pow(a1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 768*b1*c1*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,6) + 
768*a1*d1*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,6) - 6912*a1*b1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,6) - 
13440*pow(a1,2)*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,6) - 768*b1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) + 
3780*a1*c1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 768*a1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 
3780*a1*c1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) + 23040*a1*b1*d*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) + 
1890*d*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 1890*k1*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 
4608*a1*b1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 
15120*k1*pow(a1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) + 
840*pow(a1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 11520*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,6) - 
41580*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,6) - 
11520*a1*b1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,6) + 
41580*d*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,6) - 
13860*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,6) - 
13860*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,6) - 3780*a1*c1*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 
23040*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 23040*a1*b1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 
83160*d*k1*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 1890*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 
15120*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 3780*a1*c1*pow(M_E,k1*pow(RP,-1))*pow(RP,7) + 
46080*a1*b1*d*pow(M_E,k1*pow(RP,-1))*pow(RP,7) + 1890*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) + 
139860*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 3780*a1*c1*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,7) - 
1890*pow(b1,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,7) + 
26460*pow(a1,2)*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,7) + 3780*a1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) - 
23040*a1*b1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) + 23040*a1*b1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) - 
83160*d*k1*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) + 1890*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) + 
15120*pow(a1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) - 
41580*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,7) + 
41580*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,7) - 23040*a1*b1*pow(M_E,d*pow(RP,-1))*pow(RP,8) - 
83160*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,8) - 83160*k1*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,8) + 
23040*a1*b1*pow(M_E,k1*pow(RP,-1))*pow(RP,8) + 166320*d*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,8) + 
23040*a1*b1*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,8) - 23040*a1*b1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,8) + 
83160*d*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,8) - 83160*k1*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,8) - 
83160*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,9) + 83160*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,9) - 
83160*pow(a1,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,9) + 83160*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,9))*
pow(2*pow(d,3),-1);

    return res;
}



//  4
double BMSH_Integral_2233_case_3( double AP, double RP, double k1, double k2, double a1, double b1, double c1, double d1, double d )
{
double res = 3*AP*pow(M_E,-(d*pow(RP,-1)))*pow(RP,2)*(RP*(-gsl_sf_expint_Ei(-(k1*pow(RP,-1))) + gsl_sf_expint_Ei(-(k2*pow(RP,-1))))*pow(d,3)*
pow(d1,2)*(1 + pow(M_E,2*d*pow(RP,-1))) + RP*pow(d,3)*pow(d1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*
(pow(M_E,-(k1*pow(RP,-1))) - pow(M_E,-(k2*pow(RP,-1)))) + 
2*RP*pow(d,2)*pow(d1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*
pow(M_E,-((k1 + k2)*pow(RP,-1))) + 2*c1*d1*RP*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*
(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1))) + 
6*(gsl_sf_expint_Ei(-(k1*pow(RP,-1))) - gsl_sf_expint_Ei(-(k2*pow(RP,-1))))*pow(d,2)*pow(d1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*
pow(RP,2) - 6*pow(d,2)*pow(d1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,-(k1*pow(RP,-1))) - pow(M_E,-(k2*pow(RP,-1))))*
pow(RP,2) + 2*c1*d1*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,-(k1*pow(RP,-1))) - pow(M_E,-(k2*pow(RP,-1))))*
pow(RP,2) + 12*c1*d1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*
((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2) + 
2*b1*d1*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*
pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2) + pow(c1,2)*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*
(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2) + 
6*d*pow(d1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*
pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2) + 2*b1*d1*RP*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*
pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - 
pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
RP*pow(c1,2)*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))\
+ d*RP*pow(d1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))\
+ 12*c1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))\
+ 2*b1*c1*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))\
+ 2*a1*d1*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))\
+ 4*c1*d1*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + 
pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
12*b1*d1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + 
pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
6*pow(c1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + 
pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
pow(d1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + 
pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
15*d*(-gsl_sf_expint_Ei(-(k1*pow(RP,-1))) + gsl_sf_expint_Ei(-(k2*pow(RP,-1))))*pow(d1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*
pow(RP,3) + 15*d*pow(d1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,-(k1*pow(RP,-1))) - pow(M_E,-(k2*pow(RP,-1))))*
pow(RP,3) + 12*c1*d1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,k1*pow(RP,-1)) - pow(M_E,k2*pow(RP,-1)))*
pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3) + 12*b1*d1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*
((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3) + 
6*pow(c1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*
pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3) + 6*pow(d1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*
((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3) + 
30*c1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*
pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3) + 30*b1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*
pow(RP,3) + 15*d*pow(c1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*
pow(RP,3) + 12*c1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + 
pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,3) + 
12*b1*c1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + 
pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,3) + 
12*a1*d1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + 
pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,3) + 
2*c1*d*d1*RP*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
2*b1*c1*RP*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
2*a1*d1*RP*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
12*b1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
6*d*pow(c1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
2*a1*c1*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
pow(b1,2)*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
30*b1*c1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
30*a1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
4*b1*d1*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
2*RP*pow(c1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
2*c1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
12*b1*c1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
12*a1*d1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
12*b1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
6*pow(c1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
12*a1*c1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
6*pow(b1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
15*(gsl_sf_expint_Ei(-(k1*pow(RP,-1))) - gsl_sf_expint_Ei(-(k2*pow(RP,-1))))*pow(d1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*
pow(RP,4) + 30*c1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,-(k1*pow(RP,-1))) - pow(M_E,-(k2*pow(RP,-1))))*pow(RP,4) + 
15*pow(d1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,k1*pow(RP,-1)) - pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
pow(RP,4) + 30*c1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*
pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4) + 30*b1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*
(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4) + 
15*d*pow(c1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*
pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4) + 30*b1*c1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*
pow(RP,4) + 30*a1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*
pow(RP,4) + 30*b1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + 
pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,4) + 
15*pow(c1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + 
pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,4) + 
30*a1*c1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,4) + 
15*d*pow(b1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,4) + 
30*b1*c1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,4) + 
30*a1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,4) + 
2*b1*d*d1*RP*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
d*RP*pow(c1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
2*a1*c1*RP*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
RP*pow(b1,2)*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
12*b1*c1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
12*a1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
2*a1*b1*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
30*a1*c1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
15*d*pow(b1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
30*a1*b1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*
(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
4*b1*c1*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
4*a1*d1*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
2*b1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
pow(c1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
12*a1*c1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
6*pow(b1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
12*b1*c1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
12*a1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
12*a1*b1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
30*a1*c1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
15*pow(b1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
30*c1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,k1*pow(RP,-1)) - pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
pow(RP,5) + 30*b1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*
pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,5) + 15*pow(c1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*
((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,5) + 
30*b1*c1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + 
pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,5) + 
30*a1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + 
pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,5) + 
30*a1*c1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,5) + 
15*pow(b1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,5) + 
30*a1*b1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4)))*
pow(RP,5) + 2*b1*c1*d*RP*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 
120*k1*pow(RP,4) + 120*pow(RP,5)) - pow(M_E,k1*pow(RP,-1))*
(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5)))\
+ 2*a1*d*d1*RP*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 
120*k1*pow(RP,4) + 120*pow(RP,5)) - pow(M_E,k1*pow(RP,-1))*
(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5)))\
+ 2*a1*b1*RP*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 
120*k1*pow(RP,4) + 120*pow(RP,5)) - pow(M_E,k1*pow(RP,-1))*
(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5)))\
+ 12*a1*c1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 
120*k1*pow(RP,4) + 120*pow(RP,5)) - pow(M_E,k1*pow(RP,-1))*
(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5)))\
+ 6*d*pow(b1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 
120*k1*pow(RP,4) + 120*pow(RP,5)) - pow(M_E,k1*pow(RP,-1))*
(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5)))\
+ pow(a1,2)*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 
120*k1*pow(RP,4) + 120*pow(RP,5)) - pow(M_E,k1*pow(RP,-1))*
(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5)))\
+ 30*a1*b1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 
120*k1*pow(RP,4) + 120*pow(RP,5)) - pow(M_E,k1*pow(RP,-1))*
(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5)))\
+ 15*d*pow(a1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*
(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 
120*k1*pow(RP,4) + 120*pow(RP,5)) - pow(M_E,k1*pow(RP,-1))*
(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5)))\
+ 4*a1*c1*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 
120*k1*pow(RP,4) + 120*pow(RP,5))) + 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 
120*k2*pow(RP,4) + 120*pow(RP,5))) + 2*RP*pow(b1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*
pow(M_E,-((k1 + k2)*pow(RP,-1)))*(-(pow(M_E,k2*pow(RP,-1))*
(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))
) + pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 
120*k2*pow(RP,4) + 120*pow(RP,5))) + 2*b1*c1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
pow(RP,2)*(-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 
120*k1*pow(RP,4) + 120*pow(RP,5))) + 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 
120*k2*pow(RP,4) + 120*pow(RP,5))) + 2*a1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
pow(RP,2)*(-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 
120*k1*pow(RP,4) + 120*pow(RP,5))) + 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 
120*k2*pow(RP,4) + 120*pow(RP,5))) + 12*a1*b1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*
pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*(-(pow(M_E,k2*pow(RP,-1))*
(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))
) + pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 
120*k2*pow(RP,4) + 120*pow(RP,5))) + 12*a1*c1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
pow(RP,3)*(-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 
120*k1*pow(RP,4) + 120*pow(RP,5))) + 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 
120*k2*pow(RP,4) + 120*pow(RP,5))) + 6*pow(b1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
pow(RP,3)*(-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 
120*k1*pow(RP,4) + 120*pow(RP,5))) + 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 
120*k2*pow(RP,4) + 120*pow(RP,5))) + 6*pow(a1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*
pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*(-(pow(M_E,k2*pow(RP,-1))*
(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))
) + pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 
120*k2*pow(RP,4) + 120*pow(RP,5))) + 30*a1*b1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
pow(RP,4)*(-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 
120*k1*pow(RP,4) + 120*pow(RP,5))) + 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 
120*k2*pow(RP,4) + 120*pow(RP,5))) + 15*pow(a1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
pow(RP,5)*(-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 
120*k1*pow(RP,4) + 120*pow(RP,5))) + 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 
120*k2*pow(RP,4) + 120*pow(RP,5))) + 2*a1*c1*d*RP*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 
360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 720*pow(RP,6)) - 
pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 
360*pow(k2,2)*pow(RP,4) + 720*k2*pow(RP,5) + 720*pow(RP,6))) + 
d*RP*pow(b1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 
360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 720*pow(RP,6)) - 
pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 
360*pow(k2,2)*pow(RP,4) + 720*k2*pow(RP,5) + 720*pow(RP,6))) + 
RP*pow(a1,2)*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 
360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 720*pow(RP,6)) - 
pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 
360*pow(k2,2)*pow(RP,4) + 720*k2*pow(RP,5) + 720*pow(RP,6))) + 
12*a1*b1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 
360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 720*pow(RP,6)) - 
pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 
360*pow(k2,2)*pow(RP,4) + 720*k2*pow(RP,5) + 720*pow(RP,6))) + 
15*d*pow(a1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 
360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 720*pow(RP,6)) - 
pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 
360*pow(k2,2)*pow(RP,4) + 720*k2*pow(RP,5) + 720*pow(RP,6))) + 
4*a1*b1*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 
360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 720*pow(RP,6))) + 
pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 
360*pow(k2,2)*pow(RP,4) + 720*k2*pow(RP,5) + 720*pow(RP,6))) + 
2*a1*c1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 
360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 720*pow(RP,6))) + 
pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 
360*pow(k2,2)*pow(RP,4) + 720*k2*pow(RP,5) + 720*pow(RP,6))) + 
pow(b1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 
360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 720*pow(RP,6))) + 
pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 
360*pow(k2,2)*pow(RP,4) + 720*k2*pow(RP,5) + 720*pow(RP,6))) + 
6*pow(a1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 
360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 720*pow(RP,6))) + 
pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 
360*pow(k2,2)*pow(RP,4) + 720*k2*pow(RP,5) + 720*pow(RP,6))) + 
12*a1*b1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(-(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 
360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 720*pow(RP,6))) + 
pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 
360*pow(k2,2)*pow(RP,4) + 720*k2*pow(RP,5) + 720*pow(RP,6))) + 
15*pow(a1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*
(-(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 
360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 720*pow(RP,6))) + 
pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 
360*pow(k2,2)*pow(RP,4) + 720*k2*pow(RP,5) + 720*pow(RP,6))) + 
2*a1*b1*d*RP*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(7*RP*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RP,2) + 210*pow(k1,4)*pow(RP,3) + 
840*pow(k1,3)*pow(RP,4) + 2520*pow(k1,2)*pow(RP,5) + 5040*k1*pow(RP,6) + 5040*pow(RP,7)) - 
pow(M_E,k1*pow(RP,-1))*(7*RP*pow(k2,6) + pow(k2,7) + 42*pow(k2,5)*pow(RP,2) + 210*pow(k2,4)*pow(RP,3) + 
840*pow(k2,3)*pow(RP,4) + 2520*pow(k2,2)*pow(RP,5) + 5040*k2*pow(RP,6) + 5040*pow(RP,7))) + 
6*d*pow(a1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(pow(M_E,k2*pow(RP,-1))*(7*RP*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RP,2) + 210*pow(k1,4)*pow(RP,3) + 
840*pow(k1,3)*pow(RP,4) + 2520*pow(k1,2)*pow(RP,5) + 5040*k1*pow(RP,6) + 5040*pow(RP,7)) - 
pow(M_E,k1*pow(RP,-1))*(7*RP*pow(k2,6) + pow(k2,7) + 42*pow(k2,5)*pow(RP,2) + 210*pow(k2,4)*pow(RP,3) + 
840*pow(k2,3)*pow(RP,4) + 2520*pow(k2,2)*pow(RP,5) + 5040*k2*pow(RP,6) + 5040*pow(RP,7))) + 
2*RP*pow(a1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(7*RP*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RP,2) + 210*pow(k1,4)*pow(RP,3) + 
840*pow(k1,3)*pow(RP,4) + 2520*pow(k1,2)*pow(RP,5) + 5040*k1*pow(RP,6) + 5040*pow(RP,7))) + 
pow(M_E,k1*pow(RP,-1))*(7*RP*pow(k2,6) + pow(k2,7) + 42*pow(k2,5)*pow(RP,2) + 210*pow(k2,4)*pow(RP,3) + 
840*pow(k2,3)*pow(RP,4) + 2520*pow(k2,2)*pow(RP,5) + 5040*k2*pow(RP,6) + 5040*pow(RP,7))) + 
2*a1*b1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(7*RP*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RP,2) + 210*pow(k1,4)*pow(RP,3) + 
840*pow(k1,3)*pow(RP,4) + 2520*pow(k1,2)*pow(RP,5) + 5040*k1*pow(RP,6) + 5040*pow(RP,7))) + 
pow(M_E,k1*pow(RP,-1))*(7*RP*pow(k2,6) + pow(k2,7) + 42*pow(k2,5)*pow(RP,2) + 210*pow(k2,4)*pow(RP,3) + 
840*pow(k2,3)*pow(RP,4) + 2520*pow(k2,2)*pow(RP,5) + 5040*k2*pow(RP,6) + 5040*pow(RP,7))) + 
6*pow(a1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(-(pow(M_E,k2*pow(RP,-1))*(7*RP*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RP,2) + 210*pow(k1,4)*pow(RP,3) + 
840*pow(k1,3)*pow(RP,4) + 2520*pow(k1,2)*pow(RP,5) + 5040*k1*pow(RP,6) + 5040*pow(RP,7))) + 
pow(M_E,k1*pow(RP,-1))*(7*RP*pow(k2,6) + pow(k2,7) + 42*pow(k2,5)*pow(RP,2) + 210*pow(k2,4)*pow(RP,3) + 
840*pow(k2,3)*pow(RP,4) + 2520*pow(k2,2)*pow(RP,5) + 5040*k2*pow(RP,6) + 5040*pow(RP,7))) + 
d*RP*pow(a1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(pow(M_E,k2*pow(RP,-1))*(8*RP*pow(k1,7) + pow(k1,8) + 56*pow(k1,6)*pow(RP,2) + 336*pow(k1,5)*pow(RP,3) + 
1680*pow(k1,4)*pow(RP,4) + 6720*pow(k1,3)*pow(RP,5) + 20160*pow(k1,2)*pow(RP,6) + 40320*k1*pow(RP,7) + 
40320*pow(RP,8)) - pow(M_E,k1*pow(RP,-1))*
(8*RP*pow(k2,7) + pow(k2,8) + 56*pow(k2,6)*pow(RP,2) + 336*pow(k2,5)*pow(RP,3) + 1680*pow(k2,4)*pow(RP,4) + 
6720*pow(k2,3)*pow(RP,5) + 20160*pow(k2,2)*pow(RP,6) + 40320*k2*pow(RP,7) + 40320*pow(RP,8))) + 
pow(a1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(8*RP*pow(k1,7) + pow(k1,8) + 56*pow(k1,6)*pow(RP,2) + 336*pow(k1,5)*pow(RP,3) + 
1680*pow(k1,4)*pow(RP,4) + 6720*pow(k1,3)*pow(RP,5) + 20160*pow(k1,2)*pow(RP,6) + 40320*k1*pow(RP,7) + 
40320*pow(RP,8))) + pow(M_E,k1*pow(RP,-1))*
(8*RP*pow(k2,7) + pow(k2,8) + 56*pow(k2,6)*pow(RP,2) + 336*pow(k2,5)*pow(RP,3) + 1680*pow(k2,4)*pow(RP,4) + 
6720*pow(k2,3)*pow(RP,5) + 20160*pow(k2,2)*pow(RP,6) + 40320*k2*pow(RP,7) + 40320*pow(RP,8))))*pow(2*pow(d,3),-1);

    return res;
}



/* BM SR 44 ZZ */

//  1
double BMSH_Integral_44_case_1( double AP, double RP, double k1, double k2, double a1, double b1, double c1, double d1, double d )
{
double res = -3*AP*pow(M_E,-((d + k1 + k2)*pow(RP,-1)))*pow(RP,2)*(-8*c1*d1*k2*RP*pow(d,3)*pow(M_E,k1*pow(RP,-1)) - 7*k2*RP*pow(d,2)*pow(d1,2)*pow(M_E,k1*pow(RP,-1)) - k2*pow(d,3)*pow(d1,2)*pow(M_E,k1*pow(RP,-1)) - 3*RP*pow(d,3)*pow(d1,2)*pow(M_E,k1*pow(RP,-1)) + 8*c1*d1*k1*RP*pow(d,3)*pow(M_E,k2*pow(RP,-1)) + 
7*k1*RP*pow(d,2)*pow(d1,2)*pow(M_E,k2*pow(RP,-1)) + k1*pow(d,3)*pow(d1,2)*pow(M_E,k2*pow(RP,-1)) + 3*RP*pow(d,3)*pow(d1,2)*pow(M_E,k2*pow(RP,-1)) + 2*RP*gsl_sf_expint_Ei(-(k2*pow(RP,-1)))*pow(d,3)*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1)) - 
2*RP*gsl_sf_expint_Ei(k2*pow(RP,-1))*pow(d,3)*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1)) - 8*c1*d1*k1*RP*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) + 7*k1*RP*pow(d,2)*pow(d1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) + k1*pow(d,3)*pow(d1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) - 3*RP*pow(d,3)*pow(d1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) + 
8*c1*d1*k2*RP*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) - 7*k2*RP*pow(d,2)*pow(d1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) - k2*pow(d,3)*pow(d1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) + 3*RP*pow(d,3)*pow(d1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) + 16*c1*d1*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 
2*c1*d1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 10*b1*d1*RP*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 5*RP*pow(c1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 2*d*RP*pow(d1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + pow(d,2)*pow(d1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 
16*c1*d1*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) + 2*c1*d1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) - 10*b1*d1*RP*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) - 5*RP*pow(c1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) - 
2*d*RP*pow(d1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) - pow(d,2)*pow(d1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) + 4*c1*d*d1*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 2*c1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 18*b1*d1*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 
9*RP*pow(c1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 2*b1*d1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 12*b1*c1*RP*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 12*a1*d1*RP*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + pow(c1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) - 
4*c1*d*d1*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) - 2*c1*d1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + 18*b1*d1*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + 9*RP*pow(c1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + 2*b1*d1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) - 
12*b1*c1*RP*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) - 12*a1*d1*RP*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + pow(c1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + 4*b1*d*d1*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 2*d*RP*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 
2*b1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 20*b1*c1*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 20*a1*d1*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + pow(c1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 2*b1*c1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 
2*a1*d1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 14*a1*c1*RP*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 7*RP*pow(b1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) - 4*b1*d*d1*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) - 2*d*RP*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) - 
2*b1*d1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) + 20*b1*c1*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) + 20*a1*d1*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) - pow(c1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) + 2*b1*c1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) + 
2*a1*d1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) - 14*a1*c1*RP*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) - 7*RP*pow(b1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) + 4*b1*c1*d*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 4*a1*d*d1*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 
2*b1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 2*a1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 22*a1*c1*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 11*RP*pow(b1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 2*a1*c1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 
16*a1*b1*RP*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + pow(b1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) - 4*b1*c1*d*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) - 4*a1*d*d1*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) - 2*b1*c1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) - 
2*a1*d1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) + 22*a1*c1*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) + 11*RP*pow(b1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) + 2*a1*c1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) - 
16*a1*b1*RP*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) + pow(b1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) + 4*a1*c1*d*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,6) + 2*d*RP*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,6) + 2*a1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,6) + 
24*a1*b1*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,6) + pow(b1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,6) + 2*a1*b1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,6) + 9*RP*pow(a1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,6) - 4*a1*c1*d*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6) - 
2*d*RP*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6) - 2*a1*c1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6) + 24*a1*b1*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6) - pow(b1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6) + 2*a1*b1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6) - 
9*RP*pow(a1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6) + 4*a1*b1*d*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,7) + 2*a1*b1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,7) + 13*RP*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,7) + pow(a1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,7) - 
4*a1*b1*d*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,7) - 2*a1*b1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,7) + 13*RP*pow(a1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,7) + pow(a1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,7) + 2*d*RP*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,8) + 
pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,8) - 2*d*RP*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,8) - pow(a1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,8) - 16*c1*d1*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 2*c1*d1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 
10*b1*d1*RP*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 5*RP*pow(c1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 2*d*RP*pow(d1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - pow(d,2)*pow(d1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 16*c1*d1*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) - 
2*c1*d1*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) + 10*b1*d1*RP*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) + 5*RP*pow(c1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) + 2*d*RP*pow(d1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) + pow(d,2)*pow(d1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) - 
4*c1*d*d1*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 2*c1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 18*b1*d1*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 9*RP*pow(c1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 2*b1*d1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 
12*b1*c1*RP*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 12*a1*d1*RP*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - pow(c1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) + 4*c1*d*d1*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) + 2*c1*d1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - 
18*b1*d1*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - 9*RP*pow(c1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - 2*b1*d1*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) + 12*b1*c1*RP*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) + 
12*a1*d1*RP*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - pow(c1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - 4*b1*d*d1*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 2*d*RP*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 2*b1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 
20*b1*c1*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 20*a1*d1*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - pow(c1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 2*b1*c1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 2*a1*d1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 
14*a1*c1*RP*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 7*RP*pow(b1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) + 4*b1*d*d1*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) + 2*d*RP*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) + 2*b1*d1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) - 
20*b1*c1*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) - 20*a1*d1*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) + pow(c1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) - 2*b1*c1*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) - 2*a1*d1*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) + 
14*a1*c1*RP*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) + 7*RP*pow(b1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) - 4*b1*c1*d*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 4*a1*d*d1*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 2*b1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 
2*a1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 22*a1*c1*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 11*RP*pow(b1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 2*a1*c1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 16*a1*b1*RP*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 
pow(b1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) + 4*b1*c1*d*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) + 4*a1*d*d1*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) + 2*b1*c1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) + 2*a1*d1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) - 
22*a1*c1*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) - 11*RP*pow(b1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) - 2*a1*c1*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) + 16*a1*b1*RP*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) - 
pow(b1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) - 4*a1*c1*d*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,6) - 2*d*RP*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,6) - 2*a1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,6) - 24*a1*b1*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,6) - 
pow(b1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,6) - 2*a1*b1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,6) - 9*RP*pow(a1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,6) + 4*a1*c1*d*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6) + 2*d*RP*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6) + 
2*a1*c1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6) - 24*a1*b1*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6) + pow(b1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6) - 2*a1*b1*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6) + 9*RP*pow(a1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6) - 
4*a1*b1*d*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,7) - 2*a1*b1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,7) - 13*RP*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,7) - pow(a1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,7) + 4*a1*b1*d*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,7) + 
2*a1*b1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,7) - 13*RP*pow(a1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,7) - pow(a1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,7) - 2*d*RP*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,8) - pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,8) + 
2*d*RP*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,8) + pow(a1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,8) - 56*c1*d1*k2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 12*c1*d1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 24*b1*d1*k2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 
12*k2*pow(c1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 16*d*k2*pow(d1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 19*pow(d,2)*pow(d1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 56*c1*d1*k1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 12*c1*d1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 
24*b1*d1*k1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 12*k1*pow(c1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 16*d*k1*pow(d1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 19*pow(d,2)*pow(d1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 
12*gsl_sf_expint_Ei(-(k2*pow(RP,-1)))*pow(d,2)*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*pow(RP,2) - 12*gsl_sf_expint_Ei(k2*pow(RP,-1))*pow(d,2)*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*pow(RP,2) - 56*c1*d1*k1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) + 12*c1*d1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) + 
24*b1*d1*k1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) + 12*k1*pow(c1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) + 16*d*k1*pow(d1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) - 19*pow(d,2)*pow(d1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) + 
56*c1*d1*k2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) - 12*c1*d1*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) - 24*b1*d1*k2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) - 12*k2*pow(c1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) - 
16*d*k2*pow(d1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) + 19*pow(d,2)*pow(d1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) + 36*c1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 78*b1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 39*pow(c1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
40*b1*c1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 40*a1*d1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 2*pow(d1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 36*c1*d*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
78*b1*d1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 39*pow(c1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 40*b1*c1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 40*a1*d1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
2*pow(d1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 4*c1*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 40*b1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 20*d*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 104*b1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
104*a1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 60*a1*c1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 30*pow(b1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 4*c1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
40*b1*d*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 20*d*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 104*b1*c1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 104*a1*d1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
60*a1*c1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 30*pow(b1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 44*b1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 4*b1*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 44*a1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 
2*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 134*a1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 67*pow(b1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 84*a1*b1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 44*b1*c1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
4*b1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 44*a1*d*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 2*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 134*a1*c1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
67*pow(b1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 84*a1*b1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 4*b1*c1*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 48*a1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 4*a1*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 
24*d*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 168*a1*b1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 56*pow(a1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 4*b1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 48*a1*c1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 
4*a1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 24*d*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 168*a1*b1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 56*pow(a1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 
4*a1*c1*pow(M_E,k2*pow(RP,-1))*pow(k1,6)*pow(RP,2) + 52*a1*b1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,6)*pow(RP,2) + 2*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,6)*pow(RP,2) + 103*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,6)*pow(RP,2) - 4*a1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6)*pow(RP,2) + 
52*a1*b1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6)*pow(RP,2) - 2*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6)*pow(RP,2) - 103*pow(a1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6)*pow(RP,2) + 4*a1*b1*pow(M_E,k2*pow(RP,-1))*pow(k1,7)*pow(RP,2) + 
28*d*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,7)*pow(RP,2) - 4*a1*b1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,7)*pow(RP,2) + 28*d*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,7)*pow(RP,2) + 2*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,8)*pow(RP,2) - 2*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,8)*pow(RP,2) - 
36*c1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 78*b1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 39*pow(c1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 40*b1*c1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 40*a1*d1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
2*pow(d1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 36*c1*d*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 78*b1*d1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 39*pow(c1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
40*b1*c1*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 40*a1*d1*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 2*pow(d1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 4*c1*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 40*b1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
20*d*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 104*b1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 104*a1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 60*a1*c1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 30*pow(b1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
4*c1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 40*b1*d*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 20*d*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 104*b1*c1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
104*a1*d1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 60*a1*c1*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 30*pow(b1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 44*b1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 
4*b1*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 44*a1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 2*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 134*a1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 67*pow(b1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 
84*a1*b1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 44*b1*c1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 4*b1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 44*a1*d*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
2*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 134*a1*c1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 67*pow(b1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 84*a1*b1*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 
4*b1*c1*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 48*a1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 4*a1*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 24*d*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 168*a1*b1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 
56*pow(a1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 4*b1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 48*a1*c1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 4*a1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 
24*d*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 168*a1*b1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 56*pow(a1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 4*a1*c1*pow(M_E,k1*pow(RP,-1))*pow(k2,6)*pow(RP,2) - 
52*a1*b1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,6)*pow(RP,2) - 2*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,6)*pow(RP,2) - 103*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,6)*pow(RP,2) + 4*a1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6)*pow(RP,2) - 52*a1*b1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6)*pow(RP,2) + 
2*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6)*pow(RP,2) + 103*pow(a1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6)*pow(RP,2) - 4*a1*b1*pow(M_E,k1*pow(RP,-1))*pow(k2,7)*pow(RP,2) - 28*d*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,7)*pow(RP,2) + 4*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,7)*pow(RP,2) - 
28*d*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,7)*pow(RP,2) - 2*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,8)*pow(RP,2) + 2*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,8)*pow(RP,2) - 132*c1*d*d1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 80*c1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 
180*b1*d1*k2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 90*k2*pow(c1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 24*b1*d1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 80*b1*c1*k2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 80*a1*d1*k2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 
12*pow(c1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 46*d*pow(d1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 16*k2*pow(d1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 132*c1*d*d1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 80*c1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 180*b1*d1*k1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 
90*k1*pow(c1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 24*b1*d1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 80*b1*c1*k1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 80*a1*d1*k1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 12*pow(c1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 
46*d*pow(d1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 16*k1*pow(d1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 30*d*gsl_sf_expint_Ei(-(k2*pow(RP,-1)))*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*pow(RP,3) - 30*d*gsl_sf_expint_Ei(k2*pow(RP,-1))*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*pow(RP,3) - 
132*c1*d*d1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) + 80*c1*d1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) + 180*b1*d1*k1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) + 90*k1*pow(c1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 24*b1*d1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 
80*b1*c1*k1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 80*a1*d1*k1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 12*pow(c1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 46*d*pow(d1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) + 16*k1*pow(d1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) + 
132*c1*d*d1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) - 80*c1*d1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) - 180*b1*d1*k2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) - 90*k2*pow(c1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 24*b1*d1*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 
80*b1*c1*k2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 80*a1*d1*k2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 12*pow(c1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 46*d*pow(d1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) - 16*k2*pow(d1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 
36*c1*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 180*b1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 90*d*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 336*b1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 336*a1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
180*a1*c1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 90*pow(b1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 36*c1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 180*b1*d*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
90*d*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 336*b1*c1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 336*a1*d1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 180*a1*c1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
90*pow(b1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 236*b1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 40*b1*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 236*a1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 20*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
560*a1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 280*pow(b1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 336*a1*b1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 236*b1*c1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
40*b1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 236*a1*d*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 20*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 560*a1*c1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
280*pow(b1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 336*a1*b1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 44*b1*c1*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 300*a1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 44*a1*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 
150*d*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 864*a1*b1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 280*pow(a1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 44*b1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 
300*a1*c1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 44*a1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 150*d*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 864*a1*b1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 
280*pow(a1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 48*a1*c1*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,3) + 372*a1*b1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,3) + 24*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,3) + 630*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,3) + 
48*a1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,3) - 372*a1*b1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,3) + 24*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,3) + 630*pow(a1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,3) + 
52*a1*b1*pow(M_E,k2*pow(RP,-1))*pow(k1,6)*pow(RP,3) + 226*d*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,6)*pow(RP,3) + 52*a1*b1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6)*pow(RP,3) - 226*d*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6)*pow(RP,3) + 28*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,7)*pow(RP,3) + 
28*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,7)*pow(RP,3) - 36*c1*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 180*b1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 90*d*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 336*b1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
336*a1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 180*a1*c1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 90*pow(b1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 36*c1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
180*b1*d*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 90*d*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 336*b1*c1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 336*a1*d1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
180*a1*c1*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 90*pow(b1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 236*b1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 40*b1*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 236*a1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 
20*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 560*a1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 280*pow(b1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 336*a1*b1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
236*b1*c1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 40*b1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 236*a1*d*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 20*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 
560*a1*c1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 280*pow(b1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 336*a1*b1*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 44*b1*c1*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 
300*a1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 44*a1*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 150*d*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 864*a1*b1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 280*pow(a1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 
44*b1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 300*a1*c1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 44*a1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 150*d*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 
864*a1*b1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 280*pow(a1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 48*a1*c1*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 372*a1*b1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 24*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 
630*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 48*a1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,3) + 372*a1*b1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 24*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 
630*pow(a1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 52*a1*b1*pow(M_E,k1*pow(RP,-1))*pow(k2,6)*pow(RP,3) - 226*d*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,6)*pow(RP,3) - 52*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6)*pow(RP,3) + 
226*d*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6)*pow(RP,3) - 28*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,7)*pow(RP,3) - 28*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,7)*pow(RP,3) - 
2*RP*gsl_sf_expint_Ei(-(k1*pow(RP,-1)))*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*(6*RP*pow(d,2) + pow(d,3) + 15*d*pow(RP,2) + 15*pow(RP,3)) + 2*RP*gsl_sf_expint_Ei(k1*pow(RP,-1))*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*(6*RP*pow(d,2) + pow(d,3) + 15*d*pow(RP,2) + 15*pow(RP,3)) - 
192*c1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 132*c1*d1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 420*b1*d*d1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 210*d*k2*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 180*b1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 672*b1*c1*k2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 
672*a1*d1*k2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 90*pow(c1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 80*b1*c1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 80*a1*d1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 360*a1*c1*k2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 
180*k2*pow(b1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 46*pow(d1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 192*c1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 132*c1*d1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 420*b1*d*d1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 210*d*k1*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
180*b1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 672*b1*c1*k1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 672*a1*d1*k1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 90*pow(c1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 80*b1*c1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
80*a1*d1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 360*a1*c1*k1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 180*k1*pow(b1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 46*pow(d1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 30*gsl_sf_expint_Ei(-(k2*pow(RP,-1)))*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*pow(RP,4) - 
30*gsl_sf_expint_Ei(k2*pow(RP,-1))*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*pow(RP,4) + 192*c1*d*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 132*c1*d1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 420*b1*d*d1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 210*d*k1*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 
180*b1*d1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 672*b1*c1*k1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 672*a1*d1*k1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 90*pow(c1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 
80*b1*c1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 80*a1*d1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 360*a1*c1*k1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 180*k1*pow(b1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 46*pow(d1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 
192*c1*d*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 132*c1*d1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 420*b1*d*d1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 210*d*k2*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 180*b1*d1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 
672*b1*c1*k2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 672*a1*d1*k2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 90*pow(c1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 80*b1*c1*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 
80*a1*d1*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 360*a1*c1*k2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 180*k2*pow(b1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 46*pow(d1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 768*b1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
180*b1*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 768*a1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 90*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 1680*a1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 840*pow(b1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
1008*a1*b1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 768*b1*c1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 180*b1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 768*a1*d*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
90*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 1680*a1*c1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 840*pow(b1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 1008*a1*b1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
236*b1*c1*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 1260*a1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 236*a1*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 630*d*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 3456*a1*b1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 
1120*pow(a1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 236*b1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 1260*a1*c1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 236*a1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 
630*d*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 3456*a1*b1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 1120*pow(a1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 300*a1*c1*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,4) + 
1920*a1*b1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,4) + 150*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,4) + 3150*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,4) - 300*a1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,4) + 1920*a1*b1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,4) - 
150*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,4) - 3150*pow(a1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,4) + 372*a1*b1*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,4) + 1386*d*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,4) - 
372*a1*b1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,4) + 1386*d*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,4) + 226*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,6)*pow(RP,4) - 226*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6)*pow(RP,4) - 768*b1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
180*b1*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 768*a1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 90*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 1680*a1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 840*pow(b1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
1008*a1*b1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 768*b1*c1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 180*b1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 768*a1*d*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 
90*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 1680*a1*c1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 840*pow(b1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 1008*a1*b1*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
236*b1*c1*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 1260*a1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 236*a1*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 630*d*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 3456*a1*b1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 
1120*pow(a1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 236*b1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 1260*a1*c1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 236*a1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 
630*d*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 3456*a1*b1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 1120*pow(a1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 300*a1*c1*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,4) - 
1920*a1*b1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,4) - 150*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,4) - 3150*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,4) + 300*a1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,4) - 1920*a1*b1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,4) + 
150*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,4) + 3150*pow(a1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,4) - 372*a1*b1*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,4) - 1386*d*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,4) + 
372*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,4) - 1386*d*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,4) - 226*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,6)*pow(RP,4) + 226*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6)*pow(RP,4) - 192*c1*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 
420*b1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 1536*b1*c1*d*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 420*b1*d1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 1536*a1*d*d1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 210*d*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 210*k2*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 
672*b1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 672*a1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 3360*a1*c1*k2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 1680*k2*pow(b1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 360*a1*c1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 
2016*a1*b1*k2*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 180*pow(b1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 192*c1*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 420*b1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 1536*b1*c1*d*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 420*b1*d1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
1536*a1*d*d1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 210*d*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 210*k1*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 672*b1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 672*a1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 3360*a1*c1*k1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
1680*k1*pow(b1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 360*a1*c1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 2016*a1*b1*k1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 180*pow(b1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 192*c1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 
420*b1*d*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 1536*b1*c1*d*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) + 420*b1*d1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 1536*a1*d*d1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 210*d*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) + 
210*k1*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) + 672*b1*c1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) + 672*a1*d1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) + 3360*a1*c1*k1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) + 
1680*k1*pow(b1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 360*a1*c1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 2016*a1*b1*k1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 180*pow(b1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 
192*c1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 420*b1*d*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 1536*b1*c1*d*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) - 420*b1*d1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 1536*a1*d*d1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 
210*d*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) - 210*k2*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) - 672*b1*c1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) - 672*a1*d1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) - 3360*a1*c1*k2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) - 
1680*k2*pow(b1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 360*a1*c1*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 2016*a1*b1*k2*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 180*pow(b1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 
768*b1*c1*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 3780*a1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 768*a1*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 1890*d*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 10368*a1*b1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 
3360*pow(a1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 768*b1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 3780*a1*c1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 768*a1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 
1890*d*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 10368*a1*b1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 3360*pow(a1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 1260*a1*c1*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,5) + 
7680*a1*b1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,5) + 630*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,5) + 12600*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,5) + 1260*a1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,5) - 7680*a1*b1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,5) + 
630*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,5) + 12600*pow(a1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,5) + 1920*a1*b1*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,5) + 6930*d*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,5) + 
1920*a1*b1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,5) - 6930*d*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,5) + 1386*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,5) + 1386*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,5) - 
768*b1*c1*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 3780*a1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 768*a1*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 1890*d*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 10368*a1*b1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 
3360*pow(a1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 768*b1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 3780*a1*c1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 768*a1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 
1890*d*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 10368*a1*b1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 3360*pow(a1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 1260*a1*c1*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 
7680*a1*b1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 630*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 12600*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 1260*a1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,5) + 7680*a1*b1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 
630*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 12600*pow(a1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 1920*a1*b1*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,5) - 6930*d*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,5) - 
1920*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,5) + 6930*d*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,5) - 1386*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,5) - 1386*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,5) - 1536*b1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 
420*b1*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 1536*a1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 1536*b1*c1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 7560*a1*c1*d*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 1536*a1*d1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 3780*d*k2*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 
210*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 3360*a1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 20736*a1*b1*k2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 1680*pow(b1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 2016*a1*b1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 
6720*k2*pow(a1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 1536*b1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 420*b1*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 1536*a1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 1536*b1*c1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 7560*a1*c1*d*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 
1536*a1*d1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 3780*d*k1*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 210*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 3360*a1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 20736*a1*b1*k1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 
1680*pow(b1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 2016*a1*b1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 6720*k1*pow(a1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 1536*b1*c1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 420*b1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) + 
1536*a1*d*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 1536*b1*c1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) + 7560*a1*c1*d*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 1536*a1*d1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) + 3780*d*k1*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 
210*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 3360*a1*c1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 20736*a1*b1*k1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 1680*pow(b1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) + 
2016*a1*b1*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) + 6720*k1*pow(a1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 1536*b1*c1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 420*b1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) - 1536*a1*d*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 
1536*b1*c1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) - 7560*a1*c1*d*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 1536*a1*d1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) - 3780*d*k2*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 210*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 
3360*a1*c1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 20736*a1*b1*k2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 1680*pow(b1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) - 2016*a1*b1*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) - 
6720*k2*pow(a1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 3780*a1*c1*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,6) + 23040*a1*b1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,6) + 1890*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,6) + 37800*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,6) - 
3780*a1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,6) + 23040*a1*b1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,6) - 1890*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,6) - 37800*pow(a1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,6) + 
7680*a1*b1*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,6) + 27720*d*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,6) - 7680*a1*b1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,6) + 27720*d*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,6) + 6930*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,6) - 
6930*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,6) - 3780*a1*c1*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,6) - 23040*a1*b1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,6) - 1890*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,6) - 37800*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,6) + 
3780*a1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,6) - 23040*a1*b1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,6) + 1890*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,6) + 37800*pow(a1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,6) - 
7680*a1*b1*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,6) - 27720*d*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,6) + 7680*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,6) - 27720*d*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,6) - 6930*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,6) + 
6930*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,6) - 1536*b1*c1*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 7560*a1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 1536*a1*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 7560*a1*c1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 46080*a1*b1*d*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 
3780*d*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 3780*k2*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 20736*a1*b1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 75600*k2*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 6720*pow(a1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) + 
1536*b1*c1*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 7560*a1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 1536*a1*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 7560*a1*c1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 46080*a1*b1*d*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 3780*d*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 
3780*k1*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 20736*a1*b1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 75600*k1*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 6720*pow(a1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 1536*b1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) - 
7560*a1*c1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) + 1536*a1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) + 7560*a1*c1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) - 46080*a1*b1*d*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) - 3780*d*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) + 
3780*k1*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) + 20736*a1*b1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) + 75600*k1*pow(a1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) - 6720*pow(a1,2)*pow(d,3)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) - 
1536*b1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) + 7560*a1*c1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) - 1536*a1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) - 7560*a1*c1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) + 46080*a1*b1*d*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) + 
3780*d*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) - 3780*k2*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) - 20736*a1*b1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) - 75600*k2*pow(a1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) + 
6720*pow(a1,2)*pow(d,3)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) + 23040*a1*b1*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,7) + 83160*d*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,7) + 23040*a1*b1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,7) - 
83160*d*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,7) + 27720*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,7) + 27720*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,7) - 23040*a1*b1*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,7) - 
83160*d*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,7) - 23040*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,7) + 83160*d*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,7) - 27720*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,7) - 
27720*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,7) - 7560*a1*c1*pow(M_E,k1*pow(RP,-1))*pow(RP,8) - 46080*a1*b1*d*pow(M_E,k1*pow(RP,-1))*pow(RP,8) - 46080*a1*b1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,8) - 166320*d*k2*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,8) - 
3780*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,8) - 75600*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,8) + 7560*a1*c1*pow(M_E,k2*pow(RP,-1))*pow(RP,8) + 46080*a1*b1*d*pow(M_E,k2*pow(RP,-1))*pow(RP,8) + 46080*a1*b1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,8) + 166320*d*k1*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,8) + 
3780*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,8) + 75600*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,8) - 7560*a1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,8) + 46080*a1*b1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,8) - 46080*a1*b1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,8) + 
166320*d*k1*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,8) - 3780*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,8) - 75600*pow(a1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,8) + 7560*a1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,8) - 46080*a1*b1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,8) + 
46080*a1*b1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,8) - 166320*d*k2*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,8) + 3780*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,8) + 75600*pow(a1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,8) + 83160*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,8) - 
83160*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,8) - 83160*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,8) + 83160*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,8) - 46080*a1*b1*pow(M_E,k1*pow(RP,-1))*pow(RP,9) - 166320*d*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,9) - 
166320*k2*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,9) + 46080*a1*b1*pow(M_E,k2*pow(RP,-1))*pow(RP,9) + 166320*d*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,9) + 166320*k1*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,9) + 46080*a1*b1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,9) - 
166320*d*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,9) + 166320*k1*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,9) - 46080*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,9) + 166320*d*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,9) - 166320*k2*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,9) - 
166320*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,10) + 166320*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,10) - 166320*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,10) + 166320*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,10))*pow(2*pow(d,3),-1);

    return res;
}


//  2
double BMSH_Integral_44_case_2_sub_1( double AP, double RP, double k1, double k2, double a1, double b1, double c1, double d1, double d )
{
double res = 3*AP*pow(M_E,-((2*d + k2)*pow(RP,-1)))*pow(RP,2)*(8*c1*d1*k2*RP*pow(d,3)*pow(M_E,d*pow(RP,-1)) + 7*k2*RP*pow(d,2)*pow(d1,2)*pow(M_E,d*pow(RP,-1)) + k2*pow(d,3)*pow(d1,2)*pow(M_E,d*pow(RP,-1)) + 3*RP*pow(d,3)*pow(d1,2)*pow(M_E,d*pow(RP,-1)) + 8*c1*d1*k2*RP*pow(d,3)*pow(M_E,3*d*pow(RP,-1)) - 
7*k2*RP*pow(d,2)*pow(d1,2)*pow(M_E,3*d*pow(RP,-1)) + k2*pow(d,3)*pow(d1,2)*pow(M_E,3*d*pow(RP,-1)) + 3*RP*pow(d,3)*pow(d1,2)*pow(M_E,3*d*pow(RP,-1)) - 28*c1*d1*RP*pow(d,4)*pow(M_E,k2*pow(RP,-1)) - 4*c1*d1*pow(d,5)*pow(M_E,k2*pow(RP,-1)) - 32*b1*d1*RP*pow(d,5)*pow(M_E,k2*pow(RP,-1)) - 
16*RP*pow(c1,2)*pow(d,5)*pow(M_E,k2*pow(RP,-1)) - 4*b1*d1*pow(d,6)*pow(M_E,k2*pow(RP,-1)) - 36*b1*c1*RP*pow(d,6)*pow(M_E,k2*pow(RP,-1)) - 36*a1*d1*RP*pow(d,6)*pow(M_E,k2*pow(RP,-1)) - 2*pow(c1,2)*pow(d,6)*pow(M_E,k2*pow(RP,-1)) - 4*b1*c1*pow(d,7)*pow(M_E,k2*pow(RP,-1)) - 4*a1*d1*pow(d,7)*pow(M_E,k2*pow(RP,-1)) - 
40*a1*c1*RP*pow(d,7)*pow(M_E,k2*pow(RP,-1)) - 20*RP*pow(b1,2)*pow(d,7)*pow(M_E,k2*pow(RP,-1)) - 4*a1*c1*pow(d,8)*pow(M_E,k2*pow(RP,-1)) - 44*a1*b1*RP*pow(d,8)*pow(M_E,k2*pow(RP,-1)) - 2*pow(b1,2)*pow(d,8)*pow(M_E,k2*pow(RP,-1)) - 4*a1*b1*pow(d,9)*pow(M_E,k2*pow(RP,-1)) - 24*RP*pow(a1,2)*pow(d,9)*pow(M_E,k2*pow(RP,-1)) - 
2*pow(a1,2)*pow(d,10)*pow(M_E,k2*pow(RP,-1)) - 12*RP*pow(d,3)*pow(d1,2)*pow(M_E,k2*pow(RP,-1)) - 2*pow(d,4)*pow(d1,2)*pow(M_E,k2*pow(RP,-1)) + 4*c1*d1*RP*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 4*b1*d1*RP*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 2*RP*pow(c1,2)*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 
4*b1*c1*RP*pow(d,6)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 4*a1*d1*RP*pow(d,6)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 4*a1*c1*RP*pow(d,7)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 2*RP*pow(b1,2)*pow(d,7)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 4*a1*b1*RP*pow(d,8)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 
2*RP*pow(a1,2)*pow(d,9)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 2*RP*pow(d,3)*pow(d1,2)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 16*c1*d1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2) + 2*c1*d1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2) + 10*b1*d1*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2) + 
5*RP*pow(c1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2) + 2*d*RP*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2) + pow(d,2)*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2) - 16*c1*d1*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) + 2*c1*d1*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) + 
10*b1*d1*RP*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) + 5*RP*pow(c1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) + 2*d*RP*pow(d1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) - pow(d,2)*pow(d1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) + 4*c1*d*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 
2*c1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 18*b1*d1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 9*RP*pow(c1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 2*b1*d1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 12*b1*c1*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 
12*a1*d1*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + pow(c1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 4*c1*d*d1*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) - 2*c1*d1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) - 18*b1*d1*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) - 
9*RP*pow(c1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) + 2*b1*d1*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) + 12*b1*c1*RP*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) + 12*a1*d1*RP*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) + pow(c1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) + 
4*b1*d*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 2*d*RP*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 2*b1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 20*b1*c1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 20*a1*d1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + pow(c1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 
2*b1*c1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 2*a1*d1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 14*a1*c1*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 7*RP*pow(b1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 4*b1*d*d1*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) + 2*d*RP*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) - 
2*b1*d1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) - 20*b1*c1*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) - 20*a1*d1*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) - pow(c1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) + 2*b1*c1*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) + 
2*a1*d1*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) + 14*a1*c1*RP*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) + 7*RP*pow(b1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) + 4*b1*c1*d*RP*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 4*a1*d*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 2*b1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 
2*a1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 22*a1*c1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 11*RP*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 2*a1*c1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 16*a1*b1*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 
pow(b1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 4*b1*c1*d*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) + 4*a1*d*d1*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) - 2*b1*c1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) - 2*a1*d1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) - 22*a1*c1*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) - 
11*RP*pow(b1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) + 2*a1*c1*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) + 16*a1*b1*RP*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) + pow(b1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) + 4*a1*c1*d*RP*pow(M_E,d*pow(RP,-1))*pow(k2,6) + 
2*d*RP*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,6) + 2*a1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,6) + 24*a1*b1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,6) + pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,6) + 2*a1*b1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,6) + 9*RP*pow(a1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,6) + 
4*a1*c1*d*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,6) + 2*d*RP*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6) - 2*a1*c1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6) - 24*a1*b1*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6) - pow(b1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6) + 2*a1*b1*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6) + 
9*RP*pow(a1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6) + 4*a1*b1*d*RP*pow(M_E,d*pow(RP,-1))*pow(k2,7) + 2*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,7) + 13*RP*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,7) + pow(a1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,7) + 4*a1*b1*d*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,7) - 
2*a1*b1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,7) - 13*RP*pow(a1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,7) + pow(a1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,7) + 2*d*RP*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,8) + pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,8) + 
2*d*RP*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,8) - pow(a1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,8) + 56*c1*d1*k2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 12*c1*d1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 24*b1*d1*k2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 
12*k2*pow(c1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 16*d*k2*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 19*pow(d,2)*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 56*c1*d1*k2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) + 12*c1*d1*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) + 
24*b1*d1*k2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) + 12*k2*pow(c1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) + 16*d*k2*pow(d1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) - 19*pow(d,2)*pow(d1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) - 108*c1*d1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 
146*b1*d1*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 73*pow(c1,2)*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 192*b1*c1*pow(d,5)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 192*a1*d1*pow(d,5)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 246*a1*c1*pow(d,6)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 
123*pow(b1,2)*pow(d,6)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 308*a1*b1*pow(d,7)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 189*pow(a1,2)*pow(d,8)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 37*pow(d,2)*pow(d1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 12*c1*d1*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 
18*b1*d1*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 9*pow(c1,2)*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 24*b1*c1*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 24*a1*d1*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 30*a1*c1*pow(d,6)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 
15*pow(b1,2)*pow(d,6)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 36*a1*b1*pow(d,7)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 21*pow(a1,2)*pow(d,8)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 5*pow(d,2)*pow(d1,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 36*c1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
78*b1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 39*pow(c1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 40*b1*c1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 40*a1*d1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 2*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
36*c1*d*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 78*b1*d1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 39*pow(c1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 40*b1*c1*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 40*a1*d1*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
2*pow(d1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 4*c1*d1*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 40*b1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 20*d*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 104*b1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
104*a1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 60*a1*c1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 30*pow(b1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 4*c1*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 40*b1*d*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
20*d*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 104*b1*c1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 104*a1*d1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 60*a1*c1*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
30*pow(b1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 44*b1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 4*b1*d1*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 44*a1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 2*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
134*a1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 67*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 84*a1*b1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 44*b1*c1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 4*b1*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
44*a1*d*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 2*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 134*a1*c1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 67*pow(b1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 84*a1*b1*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
4*b1*c1*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 48*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 4*a1*d1*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 24*d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 168*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 
56*pow(a1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 4*b1*c1*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 48*a1*c1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 4*a1*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 24*d*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 
168*a1*b1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 56*pow(a1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 4*a1*c1*pow(M_E,d*pow(RP,-1))*pow(k2,6)*pow(RP,2) + 52*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k2,6)*pow(RP,2) + 2*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,6)*pow(RP,2) + 
103*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,6)*pow(RP,2) - 4*a1*c1*pow(M_E,3*d*pow(RP,-1))*pow(k2,6)*pow(RP,2) + 52*a1*b1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,6)*pow(RP,2) - 2*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6)*pow(RP,2) - 103*pow(a1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6)*pow(RP,2) + 
4*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k2,7)*pow(RP,2) + 28*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,7)*pow(RP,2) - 4*a1*b1*pow(M_E,3*d*pow(RP,-1))*pow(k2,7)*pow(RP,2) + 28*d*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,7)*pow(RP,2) + 2*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,8)*pow(RP,2) - 
2*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,8)*pow(RP,2) + 132*c1*d*d1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 80*c1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 180*b1*d1*k2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 90*k2*pow(c1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 
24*b1*d1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 80*b1*c1*k2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 80*a1*d1*k2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 12*pow(c1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 46*d*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 16*k2*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 
132*c1*d*d1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 80*c1*d1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 180*b1*d1*k2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 90*k2*pow(c1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) + 24*b1*d1*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) + 
80*b1*c1*k2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) + 80*a1*d1*k2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) + 12*pow(c1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) + 46*d*pow(d1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 16*k2*pow(d1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 
248*c1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 424*b1*d1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 212*pow(c1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 696*b1*c1*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 696*a1*d1*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 
1088*a1*c1*pow(d,5)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 544*pow(b1,2)*pow(d,5)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 1624*a1*b1*pow(d,6)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 1164*pow(a1,2)*pow(d,7)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 62*d*pow(d1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 
16*c1*d1*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) + 16*b1*d1*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) + 8*pow(c1,2)*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) + 64*b1*c1*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) + 64*a1*d1*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) + 
128*a1*c1*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) + 64*pow(b1,2)*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) + 208*a1*b1*pow(d,6)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) + 152*pow(a1,2)*pow(d,7)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) - 30*d*pow(d1,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) + 
36*c1*d1*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 180*b1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 90*d*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 336*b1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 336*a1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
180*a1*c1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 90*pow(b1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 36*c1*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 180*b1*d*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 90*d*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
336*b1*c1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 336*a1*d1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 180*a1*c1*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 90*pow(b1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 236*b1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
40*b1*d1*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 236*a1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 20*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 560*a1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 280*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
336*a1*b1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 236*b1*c1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 40*b1*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 236*a1*d*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 20*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 
560*a1*c1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 280*pow(b1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 336*a1*b1*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 44*b1*c1*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 300*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 
44*a1*d1*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 150*d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 864*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 280*pow(a1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 44*b1*c1*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 
300*a1*c1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 44*a1*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 150*d*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 864*a1*b1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 280*pow(a1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 
48*a1*c1*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,3) + 372*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,3) + 24*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,3) + 630*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 48*a1*c1*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,3) + 
372*a1*b1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 24*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 630*pow(a1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,3) + 52*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k2,6)*pow(RP,3) + 226*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,6)*pow(RP,3) - 
52*a1*b1*pow(M_E,3*d*pow(RP,-1))*pow(k2,6)*pow(RP,3) + 226*d*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6)*pow(RP,3) + 28*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,7)*pow(RP,3) - 28*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,7)*pow(RP,3) + 
2*RP*gsl_sf_expint_Ei(-(d*pow(RP,-1)))*pow(d1,2)*pow(M_E,(d + k2)*pow(RP,-1))*(-6*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1))) + pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1))) + 15*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,2) - 15*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,3)) - 
2*RP*gsl_sf_expint_Ei(-(k2*pow(RP,-1)))*pow(d1,2)*pow(M_E,(d + k2)*pow(RP,-1))*(-6*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1))) + pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1))) + 15*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,2) - 15*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,3)) + 192*c1*d*d1*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 
132*c1*d1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 420*b1*d*d1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 210*d*k2*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 180*b1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 672*b1*c1*k2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 672*a1*d1*k2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 
90*pow(c1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 80*b1*c1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 80*a1*d1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 360*a1*c1*k2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 180*k2*pow(b1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 
46*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 192*c1*d*d1*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 132*c1*d1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 420*b1*d*d1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 210*d*k2*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 180*b1*d1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 
672*b1*c1*k2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 672*a1*d1*k2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 90*pow(c1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 80*b1*c1*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 80*a1*d1*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 
360*a1*c1*k2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 180*k2*pow(b1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 46*pow(d1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 324*c1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 780*b1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 
390*pow(c1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 1756*b1*c1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 1756*a1*d1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 3600*a1*c1*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 1800*pow(b1,2)*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 
6756*a1*b1*pow(d,5)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 5882*pow(a1,2)*pow(d,6)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 46*pow(d1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 60*c1*d*d1*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) - 60*b1*d1*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) - 
30*pow(c1,2)*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) + 60*b1*c1*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) + 60*a1*d1*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) + 360*a1*c1*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) + 180*pow(b1,2)*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) + 
900*a1*b1*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) + 870*pow(a1,2)*pow(d,6)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) + 46*pow(d1,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) + 768*b1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 180*b1*d1*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 
768*a1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 90*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 1680*a1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 840*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 1008*a1*b1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 
768*b1*c1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 180*b1*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 768*a1*d*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 90*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 1680*a1*c1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
840*pow(b1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 1008*a1*b1*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 236*b1*c1*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 1260*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 236*a1*d1*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 
630*d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 3456*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 1120*pow(a1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 236*b1*c1*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 1260*a1*c1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 
236*a1*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 630*d*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 3456*a1*b1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 1120*pow(a1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 300*a1*c1*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,4) + 
1920*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,4) + 150*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,4) + 3150*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,4) - 300*a1*c1*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,4) + 1920*a1*b1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,4) - 
150*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,4) - 3150*pow(a1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,4) + 372*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,4) + 1386*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,4) - 372*a1*b1*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,4) + 
1386*d*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,4) + 226*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,6)*pow(RP,4) - 226*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6)*pow(RP,4) + 192*c1*d1*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 420*b1*d*d1*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 1536*b1*c1*d*k2*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 
420*b1*d1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 1536*a1*d*d1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 210*d*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 210*k2*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 672*b1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 672*a1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 
3360*a1*c1*k2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 1680*k2*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 360*a1*c1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 2016*a1*b1*k2*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 180*pow(b1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 
192*c1*d1*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 420*b1*d*d1*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 1536*b1*c1*d*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 420*b1*d1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 1536*a1*d*d1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 210*d*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 
210*k2*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 672*b1*c1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 672*a1*d1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 3360*a1*c1*k2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 1680*k2*pow(b1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 
360*a1*c1*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 2016*a1*b1*k2*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 180*pow(b1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 192*c1*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 840*b1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 420*d*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 
2976*b1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 2976*a1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 8760*a1*c1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 4380*pow(b1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 21984*a1*b1*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 
24276*pow(a1,2)*pow(d,5)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 192*c1*d1*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) - 96*b1*c1*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) - 96*a1*d1*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) + 480*a1*c1*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) + 
240*pow(b1,2)*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) + 2592*a1*b1*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) + 3696*pow(a1,2)*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) + 768*b1*c1*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 3780*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 
768*a1*d1*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 1890*d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 10368*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 3360*pow(a1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 768*b1*c1*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 
3780*a1*c1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 768*a1*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 1890*d*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 10368*a1*b1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 3360*pow(a1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 
1260*a1*c1*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,5) + 7680*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,5) + 630*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,5) + 12600*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 1260*a1*c1*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,5) + 
7680*a1*b1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 630*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 12600*pow(a1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,5) + 1920*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,5) + 6930*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,5) - 
1920*a1*b1*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,5) + 6930*d*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,5) + 1386*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,5) - 1386*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,5) + 1536*b1*c1*d*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 
420*b1*d1*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 1536*a1*d*d1*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 1536*b1*c1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 7560*a1*c1*d*k2*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 1536*a1*d1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 3780*d*k2*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 
210*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 3360*a1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 20736*a1*b1*k2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 1680*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 2016*a1*b1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 
6720*k2*pow(a1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 1536*b1*c1*d*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 420*b1*d1*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) + 1536*a1*d*d1*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 1536*b1*c1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) + 7560*a1*c1*d*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 
1536*a1*d1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) + 3780*d*k2*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 210*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 3360*a1*c1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 20736*a1*b1*k2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 
1680*pow(b1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) + 2016*a1*b1*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) + 6720*k2*pow(a1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 3072*b1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 420*b1*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 3072*a1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 
210*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 14700*a1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 7350*pow(b1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 53472*a1*b1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 79170*pow(a1,2)*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 
420*b1*d1*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,6) + 210*pow(c1,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,6) - 420*a1*c1*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,6) - 210*pow(b1,2)*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,6) + 3360*a1*b1*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,6) + 
10290*pow(a1,2)*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,6) + 3780*a1*c1*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,6) + 23040*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,6) + 1890*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,6) + 37800*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,6) - 
3780*a1*c1*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,6) + 23040*a1*b1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,6) - 1890*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,6) - 37800*pow(a1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,6) + 7680*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,6) + 
27720*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,6) - 7680*a1*b1*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,6) + 27720*d*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,6) + 6930*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,6) - 6930*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,6) + 
1536*b1*c1*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 7560*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 1536*a1*d1*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 7560*a1*c1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 46080*a1*b1*d*k2*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 3780*d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 
3780*k2*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 20736*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 75600*k2*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 6720*pow(a1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 1536*b1*c1*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) + 
7560*a1*c1*d*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) - 1536*a1*d1*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) - 7560*a1*c1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) + 46080*a1*b1*d*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) + 3780*d*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) - 3780*k2*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) - 
20736*a1*b1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) - 75600*k2*pow(a1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) + 6720*pow(a1,2)*pow(d,3)*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) - 1536*b1*c1*pow(M_E,k2*pow(RP,-1))*pow(RP,7) - 15120*a1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(RP,7) - 
1536*a1*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,7) - 7560*d*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) - 89856*a1*b1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) - 193200*pow(a1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 1536*b1*c1*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,7) + 
1536*a1*d1*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,7) - 2304*a1*b1*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,7) + 13440*pow(a1,2)*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,7) + 23040*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,7) + 83160*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,7) - 
23040*a1*b1*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,7) + 83160*d*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,7) + 27720*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,7) - 27720*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,7) + 7560*a1*c1*pow(M_E,d*pow(RP,-1))*pow(RP,8) + 
46080*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(RP,8) + 46080*a1*b1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,8) + 166320*d*k2*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,8) + 3780*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,8) + 75600*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,8) - 7560*a1*c1*pow(M_E,3*d*pow(RP,-1))*pow(RP,8) + 
46080*a1*b1*d*pow(M_E,3*d*pow(RP,-1))*pow(RP,8) - 46080*a1*b1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,8) + 166320*d*k2*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,8) - 3780*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,8) - 75600*pow(a1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,8) - 7560*a1*c1*pow(M_E,k2*pow(RP,-1))*pow(RP,8) - 
92160*a1*b1*d*pow(M_E,k2*pow(RP,-1))*pow(RP,8) - 3780*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,8) - 325080*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,8) + 7560*a1*c1*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,8) + 3780*pow(b1,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,8) - 
7560*pow(a1,2)*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,8) + 83160*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,8) - 83160*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,8) + 46080*a1*b1*pow(M_E,d*pow(RP,-1))*pow(RP,9) + 166320*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,9) + 
166320*k2*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,9) - 46080*a1*b1*pow(M_E,3*d*pow(RP,-1))*pow(RP,9) + 166320*d*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,9) - 166320*k2*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,9) - 46080*a1*b1*pow(M_E,k2*pow(RP,-1))*pow(RP,9) - 332640*d*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,9) + 
46080*a1*b1*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,9) + 166320*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,10) - 166320*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,10) - 166320*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,10) + 166320*pow(a1,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,10))*pow(2*pow(d,3),-1);

    return res;
}


//  3
double BMSH_Integral_44_case_2_sub_2( double AP, double RP, double k1, double k2, double a1, double b1, double c1, double d1, double d )
{
double res = 3*AP*pow(M_E,-((2*d + k1)*pow(RP,-1)))*pow(RP,2)*(-8*c1*d1*k1*RP*pow(d,3)*pow(M_E,d*pow(RP,-1)) - 7*k1*RP*pow(d,2)*pow(d1,2)*pow(M_E,d*pow(RP,-1)) - k1*pow(d,3)*pow(d1,2)*pow(M_E,d*pow(RP,-1)) - 3*RP*pow(d,3)*pow(d1,2)*pow(M_E,d*pow(RP,-1)) + 28*c1*d1*RP*pow(d,4)*pow(M_E,k1*pow(RP,-1)) + 
4*c1*d1*pow(d,5)*pow(M_E,k1*pow(RP,-1)) + 32*b1*d1*RP*pow(d,5)*pow(M_E,k1*pow(RP,-1)) + 16*RP*pow(c1,2)*pow(d,5)*pow(M_E,k1*pow(RP,-1)) + 4*b1*d1*pow(d,6)*pow(M_E,k1*pow(RP,-1)) + 36*b1*c1*RP*pow(d,6)*pow(M_E,k1*pow(RP,-1)) + 36*a1*d1*RP*pow(d,6)*pow(M_E,k1*pow(RP,-1)) + 2*pow(c1,2)*pow(d,6)*pow(M_E,k1*pow(RP,-1)) + 
4*b1*c1*pow(d,7)*pow(M_E,k1*pow(RP,-1)) + 4*a1*d1*pow(d,7)*pow(M_E,k1*pow(RP,-1)) + 40*a1*c1*RP*pow(d,7)*pow(M_E,k1*pow(RP,-1)) + 20*RP*pow(b1,2)*pow(d,7)*pow(M_E,k1*pow(RP,-1)) + 4*a1*c1*pow(d,8)*pow(M_E,k1*pow(RP,-1)) + 44*a1*b1*RP*pow(d,8)*pow(M_E,k1*pow(RP,-1)) + 2*pow(b1,2)*pow(d,8)*pow(M_E,k1*pow(RP,-1)) + 
4*a1*b1*pow(d,9)*pow(M_E,k1*pow(RP,-1)) + 24*RP*pow(a1,2)*pow(d,9)*pow(M_E,k1*pow(RP,-1)) + 2*pow(a1,2)*pow(d,10)*pow(M_E,k1*pow(RP,-1)) + 12*RP*pow(d,3)*pow(d1,2)*pow(M_E,k1*pow(RP,-1)) + 2*pow(d,4)*pow(d1,2)*pow(M_E,k1*pow(RP,-1)) + 
2*RP*gsl_sf_expint_Ei(-(k1*pow(RP,-1)))*pow(d,3)*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1)) - 2*RP*gsl_sf_expint_Ei(k1*pow(RP,-1))*pow(d,3)*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1)) + 4*c1*d1*RP*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 4*b1*d1*RP*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 
2*RP*pow(c1,2)*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 4*b1*c1*RP*pow(d,6)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 4*a1*d1*RP*pow(d,6)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 4*a1*c1*RP*pow(d,7)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 2*RP*pow(b1,2)*pow(d,7)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 
4*a1*b1*RP*pow(d,8)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 2*RP*pow(a1,2)*pow(d,9)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 2*RP*pow(d,3)*pow(d1,2)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 8*c1*d1*k1*RP*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1)) - 7*k1*RP*pow(d,2)*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1)) - 
k1*pow(d,3)*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1)) + 3*RP*pow(d,3)*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1)) - 16*c1*d1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 2*c1*d1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 10*b1*d1*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 
5*RP*pow(c1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 2*d*RP*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - pow(d,2)*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 16*c1*d1*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) - 2*c1*d1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) + 
10*b1*d1*RP*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) + 5*RP*pow(c1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) + 2*d*RP*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) + pow(d,2)*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) - 4*c1*d*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 
2*c1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 18*b1*d1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 9*RP*pow(c1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 2*b1*d1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 12*b1*c1*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 
12*a1*d1*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - pow(c1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,3) + 4*c1*d*d1*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) + 2*c1*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - 18*b1*d1*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - 
9*RP*pow(c1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - 2*b1*d1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) + 12*b1*c1*RP*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) + 12*a1*d1*RP*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - pow(c1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - 
4*b1*d*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 2*d*RP*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 2*b1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 20*b1*c1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 20*a1*d1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - pow(c1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 
2*b1*c1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 2*a1*d1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 14*a1*c1*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 7*RP*pow(b1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,4) + 4*b1*d*d1*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) + 
2*d*RP*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) + 2*b1*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) - 20*b1*c1*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) - 20*a1*d1*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) + pow(c1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) - 
2*b1*c1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) - 2*a1*d1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) + 14*a1*c1*RP*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) + 7*RP*pow(b1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) - 4*b1*c1*d*RP*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 
4*a1*d*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 2*b1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 2*a1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 22*a1*c1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 11*RP*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 2*a1*c1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 
16*a1*b1*RP*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,5) - pow(b1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,5) + 4*b1*c1*d*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) + 4*a1*d*d1*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) + 2*b1*c1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) + 
2*a1*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) - 22*a1*c1*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) - 11*RP*pow(b1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) - 2*a1*c1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) + 16*a1*b1*RP*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) - 
pow(b1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) - 4*a1*c1*d*RP*pow(M_E,d*pow(RP,-1))*pow(k1,6) - 2*d*RP*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,6) - 2*a1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,6) - 24*a1*b1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,6) - pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,6) - 
2*a1*b1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,6) - 9*RP*pow(a1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,6) + 4*a1*c1*d*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6) + 2*d*RP*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6) + 2*a1*c1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6) - 
24*a1*b1*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6) + pow(b1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6) - 2*a1*b1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6) + 9*RP*pow(a1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6) - 4*a1*b1*d*RP*pow(M_E,d*pow(RP,-1))*pow(k1,7) - 
2*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,7) - 13*RP*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,7) - pow(a1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,7) + 4*a1*b1*d*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,7) + 2*a1*b1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,7) - 
13*RP*pow(a1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,7) - pow(a1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,7) - 2*d*RP*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,8) - pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,8) + 2*d*RP*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,8) + 
pow(a1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,8) - 56*c1*d1*k1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 12*c1*d1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 24*b1*d1*k1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 12*k1*pow(c1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 
16*d*k1*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 19*pow(d,2)*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 108*c1*d1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 146*b1*d1*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 73*pow(c1,2)*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
192*b1*c1*pow(d,5)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 192*a1*d1*pow(d,5)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 246*a1*c1*pow(d,6)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 123*pow(b1,2)*pow(d,6)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 308*a1*b1*pow(d,7)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
189*pow(a1,2)*pow(d,8)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 37*pow(d,2)*pow(d1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 12*gsl_sf_expint_Ei(-(k1*pow(RP,-1)))*pow(d,2)*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1))*pow(RP,2) - 12*gsl_sf_expint_Ei(k1*pow(RP,-1))*pow(d,2)*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1))*pow(RP,2) - 
12*c1*d1*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) - 18*b1*d1*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) - 9*pow(c1,2)*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) - 24*b1*c1*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) - 24*a1*d1*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) - 
30*a1*c1*pow(d,6)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) - 15*pow(b1,2)*pow(d,6)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) - 36*a1*b1*pow(d,7)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) - 21*pow(a1,2)*pow(d,8)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) - 5*pow(d,2)*pow(d1,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 
56*c1*d1*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) - 12*c1*d1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) - 24*b1*d1*k1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) - 12*k1*pow(c1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) - 16*d*k1*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) + 
19*pow(d,2)*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) - 36*c1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 78*b1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 39*pow(c1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 40*b1*c1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
40*a1*d1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 2*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 36*c1*d*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 78*b1*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
39*pow(c1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 40*b1*c1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 40*a1*d1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 2*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
4*c1*d1*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 40*b1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 20*d*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 104*b1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 104*a1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
60*a1*c1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 30*pow(b1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 4*c1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 40*b1*d*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
20*d*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 104*b1*c1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 104*a1*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 60*a1*c1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
30*pow(b1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 44*b1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 4*b1*d1*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 44*a1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 2*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
134*a1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 67*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 84*a1*b1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 44*b1*c1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 4*b1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
44*a1*d*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 2*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 134*a1*c1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 67*pow(b1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
84*a1*b1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 4*b1*c1*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 48*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 4*a1*d1*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 24*d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 
168*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 56*pow(a1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 4*b1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 48*a1*c1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 4*a1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 
24*d*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 168*a1*b1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 56*pow(a1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 4*a1*c1*pow(M_E,d*pow(RP,-1))*pow(k1,6)*pow(RP,2) - 
52*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k1,6)*pow(RP,2) - 2*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,6)*pow(RP,2) - 103*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,6)*pow(RP,2) + 4*a1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6)*pow(RP,2) - 52*a1*b1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6)*pow(RP,2) + 
2*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6)*pow(RP,2) + 103*pow(a1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6)*pow(RP,2) - 4*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k1,7)*pow(RP,2) - 28*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,7)*pow(RP,2) + 4*a1*b1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,7)*pow(RP,2) - 
28*d*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,7)*pow(RP,2) - 2*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,8)*pow(RP,2) + 2*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,8)*pow(RP,2) - 132*c1*d*d1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 80*c1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 
180*b1*d1*k1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 90*k1*pow(c1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 24*b1*d1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 80*b1*c1*k1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 80*a1*d1*k1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 
12*pow(c1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 46*d*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 16*k1*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 248*c1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 424*b1*d1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 212*pow(c1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 
696*b1*c1*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 696*a1*d1*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 1088*a1*c1*pow(d,5)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 544*pow(b1,2)*pow(d,5)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 1624*a1*b1*pow(d,6)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 
1164*pow(a1,2)*pow(d,7)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 62*d*pow(d1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 30*d*gsl_sf_expint_Ei(-(k1*pow(RP,-1)))*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1))*pow(RP,3) - 30*d*gsl_sf_expint_Ei(k1*pow(RP,-1))*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1))*pow(RP,3) - 
16*c1*d1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) + 16*b1*d1*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) + 8*pow(c1,2)*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) + 64*b1*c1*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) + 64*a1*d1*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) + 
128*a1*c1*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) + 64*pow(b1,2)*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) + 208*a1*b1*pow(d,6)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) + 152*pow(a1,2)*pow(d,7)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 30*d*pow(d1,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) + 
132*c1*d*d1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 80*c1*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 180*b1*d1*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 90*k1*pow(c1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 24*b1*d1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 
80*b1*c1*k1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 80*a1*d1*k1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 12*pow(c1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 46*d*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 16*k1*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 
36*c1*d1*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 180*b1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 90*d*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 336*b1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 336*a1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
180*a1*c1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 90*pow(b1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 36*c1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 180*b1*d*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
90*d*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 336*b1*c1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 336*a1*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 180*a1*c1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
90*pow(b1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 236*b1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 40*b1*d1*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 236*a1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 20*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
560*a1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 280*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 336*a1*b1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 236*b1*c1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
40*b1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 236*a1*d*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 20*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 560*a1*c1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
280*pow(b1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 336*a1*b1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 44*b1*c1*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 300*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 44*a1*d1*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 
150*d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 864*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 280*pow(a1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 44*b1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 300*a1*c1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 
44*a1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 150*d*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 864*a1*b1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 280*pow(a1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 
48*a1*c1*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,3) - 372*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,3) - 24*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,3) - 630*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,3) - 48*a1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,3) + 
372*a1*b1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,3) - 24*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,3) - 630*pow(a1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,3) - 52*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k1,6)*pow(RP,3) - 
226*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,6)*pow(RP,3) - 52*a1*b1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6)*pow(RP,3) + 226*d*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6)*pow(RP,3) - 28*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,7)*pow(RP,3) - 28*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,7)*pow(RP,3) - 
2*RP*gsl_sf_expint_Ei(-(d*pow(RP,-1)))*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1))*(6*RP*pow(d,2) + pow(d,3) + 15*d*pow(RP,2) + 15*pow(RP,3)) + 2*RP*gsl_sf_expint_Ei(d*pow(RP,-1))*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1))*(6*RP*pow(d,2) + pow(d,3) + 15*d*pow(RP,2) + 15*pow(RP,3)) - 192*c1*d*d1*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
132*c1*d1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 420*b1*d*d1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 210*d*k1*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 180*b1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 672*b1*c1*k1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 672*a1*d1*k1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
90*pow(c1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 80*b1*c1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 80*a1*d1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 360*a1*c1*k1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 180*k1*pow(b1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
46*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 324*c1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 780*b1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 390*pow(c1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 1756*b1*c1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 1756*a1*d1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 
3600*a1*c1*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 1800*pow(b1,2)*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 6756*a1*b1*pow(d,5)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 5882*pow(a1,2)*pow(d,6)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 46*pow(d1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 
30*gsl_sf_expint_Ei(-(k1*pow(RP,-1)))*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1))*pow(RP,4) - 30*gsl_sf_expint_Ei(k1*pow(RP,-1))*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1))*pow(RP,4) + 60*c1*d*d1*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) + 60*b1*d1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) + 
30*pow(c1,2)*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) - 60*b1*c1*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) - 60*a1*d1*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) - 360*a1*c1*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) - 180*pow(b1,2)*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) - 
900*a1*b1*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) - 870*pow(a1,2)*pow(d,6)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) - 46*pow(d1,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) - 192*c1*d*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 132*c1*d1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 
420*b1*d*d1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 210*d*k1*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 180*b1*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 672*b1*c1*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 672*a1*d1*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 
90*pow(c1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 80*b1*c1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 80*a1*d1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 360*a1*c1*k1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 
180*k1*pow(b1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 46*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 768*b1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 180*b1*d1*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 768*a1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
90*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 1680*a1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 840*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 1008*a1*b1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 768*b1*c1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
180*b1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 768*a1*d*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 90*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 1680*a1*c1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
840*pow(b1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 1008*a1*b1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 236*b1*c1*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 1260*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 236*a1*d1*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 
630*d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 3456*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 1120*pow(a1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 236*b1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 
1260*a1*c1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 236*a1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 630*d*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 3456*a1*b1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 
1120*pow(a1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 300*a1*c1*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,4) - 1920*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,4) - 150*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,4) - 3150*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,4) + 
300*a1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,4) - 1920*a1*b1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,4) + 150*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,4) + 3150*pow(a1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,4) - 
372*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,4) - 1386*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,4) + 372*a1*b1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,4) - 1386*d*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,4) - 226*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,6)*pow(RP,4) + 
226*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6)*pow(RP,4) - 192*c1*d1*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 420*b1*d*d1*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 1536*b1*c1*d*k1*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 420*b1*d1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 1536*a1*d*d1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 
210*d*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 210*k1*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 672*b1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 672*a1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 3360*a1*c1*k1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 
1680*k1*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 360*a1*c1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 2016*a1*b1*k1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 180*pow(b1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 192*c1*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 840*b1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 
420*d*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 2976*b1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 2976*a1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 8760*a1*c1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 4380*pow(b1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 
21984*a1*b1*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 24276*pow(a1,2)*pow(d,5)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 192*c1*d1*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) - 96*b1*c1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) - 96*a1*d1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) + 
480*a1*c1*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) + 240*pow(b1,2)*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) + 2592*a1*b1*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) + 3696*pow(a1,2)*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) - 192*c1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 
420*b1*d*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 1536*b1*c1*d*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 420*b1*d1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 1536*a1*d*d1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 210*d*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 
210*k1*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 672*b1*c1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 672*a1*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 3360*a1*c1*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 
1680*k1*pow(b1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 360*a1*c1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 2016*a1*b1*k1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 180*pow(b1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 768*b1*c1*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 
3780*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 768*a1*d1*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 1890*d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 10368*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 3360*pow(a1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 
768*b1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 3780*a1*c1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 768*a1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 1890*d*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 
10368*a1*b1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 3360*pow(a1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 1260*a1*c1*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,5) - 7680*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,5) - 
630*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,5) - 12600*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,5) - 1260*a1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,5) + 7680*a1*b1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,5) - 
630*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,5) - 12600*pow(a1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,5) - 1920*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,5) - 6930*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,5) - 
1920*a1*b1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,5) + 6930*d*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,5) - 1386*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,5) - 1386*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,5) - 1536*b1*c1*d*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 
420*b1*d1*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 1536*a1*d*d1*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 1536*b1*c1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 7560*a1*c1*d*k1*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 1536*a1*d1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 3780*d*k1*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 
210*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 3360*a1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 20736*a1*b1*k1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 1680*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 2016*a1*b1*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 
6720*k1*pow(a1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 3072*b1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 420*b1*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 3072*a1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 210*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 14700*a1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 
7350*pow(b1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 53472*a1*b1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 79170*pow(a1,2)*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 420*b1*d1*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,6) - 210*pow(c1,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,6) + 
420*a1*c1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,6) + 210*pow(b1,2)*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,6) - 3360*a1*b1*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,6) - 10290*pow(a1,2)*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,6) - 1536*b1*c1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) + 
420*b1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 1536*a1*d*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) + 1536*b1*c1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 7560*a1*c1*d*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) + 1536*a1*d1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 
3780*d*k1*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) + 210*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) + 3360*a1*c1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) + 20736*a1*b1*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) + 1680*pow(b1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 
2016*a1*b1*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 6720*k1*pow(a1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 3780*a1*c1*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,6) - 23040*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,6) - 1890*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,6) - 
37800*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,6) + 3780*a1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,6) - 23040*a1*b1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,6) + 1890*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,6) + 
37800*pow(a1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,6) - 7680*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,6) - 27720*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,6) + 7680*a1*b1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,6) - 
27720*d*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,6) - 6930*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,6) + 6930*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,6) - 1536*b1*c1*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 7560*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 
1536*a1*d1*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 7560*a1*c1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 46080*a1*b1*d*k1*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 3780*d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 3780*k1*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 20736*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 
75600*k1*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 6720*pow(a1,2)*pow(d,3)*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 1536*b1*c1*pow(M_E,k1*pow(RP,-1))*pow(RP,7) + 15120*a1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(RP,7) + 1536*a1*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,7) + 7560*d*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) + 
89856*a1*b1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) + 193200*pow(a1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) + 1536*b1*c1*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,7) + 1536*a1*d1*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,7) - 2304*a1*b1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,7) + 
13440*pow(a1,2)*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,7) - 1536*b1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) + 7560*a1*c1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) - 1536*a1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) - 7560*a1*c1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) + 
46080*a1*b1*d*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) + 3780*d*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) - 3780*k1*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) - 20736*a1*b1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) - 75600*k1*pow(a1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) + 
6720*pow(a1,2)*pow(d,3)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) - 23040*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,7) - 83160*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,7) - 23040*a1*b1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,7) + 
83160*d*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,7) - 27720*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,7) - 27720*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,7) - 7560*a1*c1*pow(M_E,d*pow(RP,-1))*pow(RP,8) - 46080*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(RP,8) - 
46080*a1*b1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,8) - 166320*d*k1*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,8) - 3780*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,8) - 75600*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,8) + 7560*a1*c1*pow(M_E,k1*pow(RP,-1))*pow(RP,8) + 92160*a1*b1*d*pow(M_E,k1*pow(RP,-1))*pow(RP,8) + 
3780*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,8) + 325080*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,8) - 7560*a1*c1*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,8) - 3780*pow(b1,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,8) + 7560*pow(a1,2)*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,8) + 
7560*a1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,8) - 46080*a1*b1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,8) + 46080*a1*b1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,8) - 166320*d*k1*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,8) + 3780*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,8) + 
75600*pow(a1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,8) - 83160*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,8) + 83160*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,8) - 46080*a1*b1*pow(M_E,d*pow(RP,-1))*pow(RP,9) - 166320*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,9) - 
166320*k1*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,9) + 46080*a1*b1*pow(M_E,k1*pow(RP,-1))*pow(RP,9) + 332640*d*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,9) + 46080*a1*b1*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,9) - 46080*a1*b1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,9) + 
166320*d*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,9) - 166320*k1*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,9) - 166320*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,10) + 166320*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,10) - 166320*pow(a1,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,10) + 
166320*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,10))*pow(2*pow(d,3),-1);

    return res;
}


// 4
double BMSH_Integral_44_case_3( double AP, double RP, double k1, double k2, double a1, double b1, double c1, double d1, double d )
{
double res = 3*AP*RP*pow(M_E,-(d*pow(RP,-1)))*(RP*pow(d,3)*pow(d1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1))) + 2*(gsl_sf_expint_Ei(-(k1*pow(RP,-1))) - gsl_sf_expint_Ei(-(k2*pow(RP,-1))))*pow(d,3)*pow(d1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,2) - 
2*pow(d,3)*pow(d1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,-(k1*pow(RP,-1))) - pow(M_E,-(k2*pow(RP,-1))))*pow(RP,2) + 4*c1*d1*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2) + 
5*pow(d,2)*pow(d1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2) + 
RP*pow(d,2)*pow(d1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
10*c1*d1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
2*c1*d1*RP*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
4*b1*d1*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
2*pow(c1,2)*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
2*d*pow(d1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
12*(-gsl_sf_expint_Ei(-(k1*pow(RP,-1))) + gsl_sf_expint_Ei(-(k2*pow(RP,-1))))*pow(d,2)*pow(d1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,3) + 12*pow(d,2)*pow(d1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,-(k1*pow(RP,-1))) - pow(M_E,-(k2*pow(RP,-1))))*pow(RP,3) + 
4*c1*d1*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,k1*pow(RP,-1)) - pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3) + 4*b1*d1*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3) + 
2*pow(c1,2)*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3) + 
12*d*pow(d1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3) + 
24*c1*d1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3) + 
24*b1*d1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,3) + 
12*pow(c1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,3) + 
2*pow(d1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,3) + 
24*c1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,3) + 
4*b1*c1*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,3) + 
4*a1*d1*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,3) + 
2*c1*d1*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
10*b1*d1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
5*pow(c1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
4*c1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
24*b1*c1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
24*a1*d1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
2*b1*d1*RP*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
RP*pow(c1,2)*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
4*c1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
4*b1*c1*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
4*a1*d1*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
24*b1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
12*d*pow(c1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
4*a1*c1*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
2*pow(b1,2)*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
30*d*(gsl_sf_expint_Ei(-(k1*pow(RP,-1))) - gsl_sf_expint_Ei(-(k2*pow(RP,-1))))*pow(d1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,4) + 24*c1*d1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,-(k1*pow(RP,-1))) - pow(M_E,-(k2*pow(RP,-1))))*pow(RP,4) + 
30*d*pow(d1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,k1*pow(RP,-1)) - pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4) + 60*c1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4) + 
24*b1*d1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4) + 
12*pow(c1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4) + 
12*pow(d1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4) + 
24*c1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,4) + 
24*b1*c1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,4) + 
24*a1*d1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,4) + 
60*b1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,4) + 
30*d*pow(c1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,4) + 
24*b1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,4) + 
12*pow(c1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,4) + 
24*a1*c1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,4) + 
12*pow(b1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,4) + 
60*b1*c1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,4) + 
60*a1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,4) + 
2*b1*d1*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
RP*pow(c1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
10*b1*c1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 10*a1*d1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
4*b1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
2*pow(c1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
24*a1*c1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 12*pow(b1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
24*b1*c1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
24*a1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
24*a1*b1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 2*b1*c1*RP*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
2*a1*d1*RP*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
4*b1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
2*d*pow(c1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
4*a1*c1*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 2*pow(b1,2)*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
24*b1*c1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
24*a1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
4*a1*b1*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 60*a1*c1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*
(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
30*d*pow(b1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
30*(-gsl_sf_expint_Ei(-(k1*pow(RP,-1))) + gsl_sf_expint_Ei(-(k2*pow(RP,-1))))*pow(d1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,5) + 30*pow(d1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,-(k1*pow(RP,-1))) - pow(M_E,-(k2*pow(RP,-1))))*pow(RP,5) + 
60*c1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,k1*pow(RP,-1)) - pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,5) + 60*b1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,5) + 
30*d*pow(c1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,5) + 60*c1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,5) + 
60*b1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,5) + 
30*pow(c1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,5) + 
60*b1*c1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,5) + 
60*a1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,5) + 
60*b1*c1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,5) + 
60*a1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,5) + 
60*a1*c1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,5) + 
30*d*pow(b1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,5) + 
60*a1*c1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4)))*pow(RP,5) + 
30*pow(b1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4)))*pow(RP,5) + 
60*a1*b1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4)))*pow(RP,5) + 
2*b1*c1*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5)) - 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
2*a1*d1*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5)) - 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
10*a1*c1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5)) - 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
5*pow(b1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5)) - 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
4*b1*c1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5)) - 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
4*a1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5)) - 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
24*a1*b1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5)) - 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
24*a1*c1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5)) - 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
12*pow(b1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5)) - 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
12*pow(a1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5)) - 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
60*a1*b1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,5)*(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5)) - 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
2*a1*c1*RP*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))) + 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
RP*pow(b1,2)*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))) + 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
4*b1*c1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*(-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))) + 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
4*a1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*(-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))) + 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
4*a1*b1*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*(-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))) + 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
24*a1*c1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*(-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))) + 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
12*d*pow(b1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*(-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))) + 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
2*pow(a1,2)*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*(-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))) + 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
60*a1*b1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*(-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))) + 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
30*d*pow(a1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,5)*(-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))) + 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 60*c1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,-(k1*pow(RP,-1))) - pow(M_E,-(k2*pow(RP,-1))))*pow(RP,6) + 
60*b1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,6) + 30*pow(c1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,6) + 
60*b1*c1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,6) + 
60*a1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,6) + 
60*a1*c1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,6) + 
30*pow(b1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,6) + 
60*a1*b1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4)))*pow(RP,6) + 
30*pow(a1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5)) - 
pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5)))*pow(RP,6) + 
2*a1*c1*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 720*pow(RP,6)) - 
pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 720*k2*pow(RP,5) + 720*pow(RP,6))) + 
RP*pow(b1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 720*pow(RP,6)) - 
pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 720*k2*pow(RP,5) + 720*pow(RP,6))) + 
10*a1*b1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 720*pow(RP,6)) - 
pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 720*k2*pow(RP,5) + 720*pow(RP,6))) + 
4*a1*c1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 720*pow(RP,6)) - 
pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 720*k2*pow(RP,5) + 720*pow(RP,6))) + 
2*pow(b1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 720*pow(RP,6)) - 
pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 720*k2*pow(RP,5) + 720*pow(RP,6))) + 
12*pow(a1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 720*pow(RP,6)) - 
pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 720*k2*pow(RP,5) + 720*pow(RP,6))) + 
24*a1*b1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 720*pow(RP,6)) - 
pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 720*k2*pow(RP,5) + 720*pow(RP,6))) + 
30*pow(a1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,5)*(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 720*pow(RP,6)) - 
pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 720*k2*pow(RP,5) + 720*pow(RP,6))) + 
2*a1*b1*RP*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(-(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 720*pow(RP,6))) + 
pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 720*k2*pow(RP,5) + 720*pow(RP,6))) + 
4*a1*c1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*(-(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 720*pow(RP,6))) + 
pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 720*k2*pow(RP,5) + 720*pow(RP,6))) + 
2*d*pow(b1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*(-(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 720*pow(RP,6))) + 
pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 720*k2*pow(RP,5) + 720*pow(RP,6))) + 
2*pow(a1,2)*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*(-(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 720*pow(RP,6))) + 
pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 720*k2*pow(RP,5) + 720*pow(RP,6))) + 
24*a1*b1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*(-(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 720*pow(RP,6))) + 
pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 720*k2*pow(RP,5) + 720*pow(RP,6))) + 
30*d*pow(a1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*(-(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 720*pow(RP,6))) + 
pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 720*k2*pow(RP,5) + 720*pow(RP,6))) + 
2*a1*b1*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(7*RP*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RP,2) + 210*pow(k1,4)*pow(RP,3) + 840*pow(k1,3)*pow(RP,4) + 2520*pow(k1,2)*pow(RP,5) + 5040*k1*pow(RP,6) + 5040*pow(RP,7)) - 
pow(M_E,k1*pow(RP,-1))*(7*RP*pow(k2,6) + pow(k2,7) + 42*pow(k2,5)*pow(RP,2) + 210*pow(k2,4)*pow(RP,3) + 840*pow(k2,3)*pow(RP,4) + 2520*pow(k2,2)*pow(RP,5) + 5040*k2*pow(RP,6) + 5040*pow(RP,7))) + 
5*pow(a1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*(pow(M_E,k2*pow(RP,-1))*(7*RP*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RP,2) + 210*pow(k1,4)*pow(RP,3) + 840*pow(k1,3)*pow(RP,4) + 2520*pow(k1,2)*pow(RP,5) + 5040*k1*pow(RP,6) + 5040*pow(RP,7)) - 
pow(M_E,k1*pow(RP,-1))*(7*RP*pow(k2,6) + pow(k2,7) + 42*pow(k2,5)*pow(RP,2) + 210*pow(k2,4)*pow(RP,3) + 840*pow(k2,3)*pow(RP,4) + 2520*pow(k2,2)*pow(RP,5) + 5040*k2*pow(RP,6) + 5040*pow(RP,7))) + 
4*a1*b1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*(pow(M_E,k2*pow(RP,-1))*(7*RP*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RP,2) + 210*pow(k1,4)*pow(RP,3) + 840*pow(k1,3)*pow(RP,4) + 2520*pow(k1,2)*pow(RP,5) + 5040*k1*pow(RP,6) + 5040*pow(RP,7)) - 
pow(M_E,k1*pow(RP,-1))*(7*RP*pow(k2,6) + pow(k2,7) + 42*pow(k2,5)*pow(RP,2) + 210*pow(k2,4)*pow(RP,3) + 840*pow(k2,3)*pow(RP,4) + 2520*pow(k2,2)*pow(RP,5) + 5040*k2*pow(RP,6) + 5040*pow(RP,7))) + 
12*pow(a1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*(pow(M_E,k2*pow(RP,-1))*(7*RP*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RP,2) + 210*pow(k1,4)*pow(RP,3) + 840*pow(k1,3)*pow(RP,4) + 2520*pow(k1,2)*pow(RP,5) + 5040*k1*pow(RP,6) + 5040*pow(RP,7)) - 
pow(M_E,k1*pow(RP,-1))*(7*RP*pow(k2,6) + pow(k2,7) + 42*pow(k2,5)*pow(RP,2) + 210*pow(k2,4)*pow(RP,3) + 840*pow(k2,3)*pow(RP,4) + 2520*pow(k2,2)*pow(RP,5) + 5040*k2*pow(RP,6) + 5040*pow(RP,7))) + 
RP*pow(a1,2)*pow(d,3)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(-(pow(M_E,k2*pow(RP,-1))*(7*RP*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RP,2) + 210*pow(k1,4)*pow(RP,3) + 840*pow(k1,3)*pow(RP,4) + 2520*pow(k1,2)*pow(RP,5) + 5040*k1*pow(RP,6) + 5040*pow(RP,7))) + 
pow(M_E,k1*pow(RP,-1))*(7*RP*pow(k2,6) + pow(k2,7) + 42*pow(k2,5)*pow(RP,2) + 210*pow(k2,4)*pow(RP,3) + 840*pow(k2,3)*pow(RP,4) + 2520*pow(k2,2)*pow(RP,5) + 5040*k2*pow(RP,6) + 5040*pow(RP,7))) + 
4*a1*b1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*(-(pow(M_E,k2*pow(RP,-1))*(7*RP*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RP,2) + 210*pow(k1,4)*pow(RP,3) + 840*pow(k1,3)*pow(RP,4) + 2520*pow(k1,2)*pow(RP,5) + 5040*k1*pow(RP,6) + 5040*pow(RP,7))) + 
pow(M_E,k1*pow(RP,-1))*(7*RP*pow(k2,6) + pow(k2,7) + 42*pow(k2,5)*pow(RP,2) + 210*pow(k2,4)*pow(RP,3) + 840*pow(k2,3)*pow(RP,4) + 2520*pow(k2,2)*pow(RP,5) + 5040*k2*pow(RP,6) + 5040*pow(RP,7))) + 
12*d*pow(a1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*(-(pow(M_E,k2*pow(RP,-1))*(7*RP*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RP,2) + 210*pow(k1,4)*pow(RP,3) + 840*pow(k1,3)*pow(RP,4) + 2520*pow(k1,2)*pow(RP,5) + 5040*k1*pow(RP,6) + 5040*pow(RP,7))) + 
pow(M_E,k1*pow(RP,-1))*(7*RP*pow(k2,6) + pow(k2,7) + 42*pow(k2,5)*pow(RP,2) + 210*pow(k2,4)*pow(RP,3) + 840*pow(k2,3)*pow(RP,4) + 2520*pow(k2,2)*pow(RP,5) + 5040*k2*pow(RP,6) + 5040*pow(RP,7))) + 
RP*pow(a1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*(pow(M_E,k2*pow(RP,-1))*(8*RP*pow(k1,7) + pow(k1,8) + 56*pow(k1,6)*pow(RP,2) + 336*pow(k1,5)*pow(RP,3) + 1680*pow(k1,4)*pow(RP,4) + 6720*pow(k1,3)*pow(RP,5) + 20160*pow(k1,2)*pow(RP,6) + 40320*k1*pow(RP,7) + 40320*pow(RP,8)) - 
pow(M_E,k1*pow(RP,-1))*(8*RP*pow(k2,7) + pow(k2,8) + 56*pow(k2,6)*pow(RP,2) + 336*pow(k2,5)*pow(RP,3) + 1680*pow(k2,4)*pow(RP,4) + 6720*pow(k2,3)*pow(RP,5) + 20160*pow(k2,2)*pow(RP,6) + 40320*k2*pow(RP,7) + 40320*pow(RP,8))) + 
2*pow(a1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*(pow(M_E,k2*pow(RP,-1))*(8*RP*pow(k1,7) + pow(k1,8) + 56*pow(k1,6)*pow(RP,2) + 336*pow(k1,5)*pow(RP,3) + 1680*pow(k1,4)*pow(RP,4) + 6720*pow(k1,3)*pow(RP,5) + 20160*pow(k1,2)*pow(RP,6) + 40320*k1*pow(RP,7) + 40320*pow(RP,8)) - 
pow(M_E,k1*pow(RP,-1))*(8*RP*pow(k2,7) + pow(k2,8) + 56*pow(k2,6)*pow(RP,2) + 336*pow(k2,5)*pow(RP,3) + 1680*pow(k2,4)*pow(RP,4) + 6720*pow(k2,3)*pow(RP,5) + 20160*pow(k2,2)*pow(RP,6) + 40320*k2*pow(RP,7) + 40320*pow(RP,8))) + 
2*d*pow(a1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*(-(pow(M_E,k2*pow(RP,-1))*(8*RP*pow(k1,7) + pow(k1,8) + 56*pow(k1,6)*pow(RP,2) + 336*pow(k1,5)*pow(RP,3) + 1680*pow(k1,4)*pow(RP,4) + 6720*pow(k1,3)*pow(RP,5) + 20160*pow(k1,2)*pow(RP,6) + 40320*k1*pow(RP,7) + 
40320*pow(RP,8))) + pow(M_E,k1*pow(RP,-1))*(8*RP*pow(k2,7) + pow(k2,8) + 56*pow(k2,6)*pow(RP,2) + 336*pow(k2,5)*pow(RP,3) + 1680*pow(k2,4)*pow(RP,4) + 6720*pow(k2,3)*pow(RP,5) + 20160*pow(k2,2)*pow(RP,6) + 40320*k2*pow(RP,7) + 40320*pow(RP,8))))*pow(2*pow(d,3),-1);

    return res;
}
