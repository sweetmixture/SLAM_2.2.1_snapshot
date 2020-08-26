#include<stdio.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_sf_expint.h>

/* !!! CALL Ei function: 
 *
 * gsl_sf_expint_Ei( arg );
 *
 */

/* !!! Integral Arguments
 *
 * k1: a knot for the left  boundary for spline cubic polynomial
 * k2: a knot for the right boundary for spline cubic polynomial
 *
 * d : the distance between sp-electrons and an external point charge
 *
 * a1, b1, c1, d1 : 1st cubic polynomial coefficients
 * a2, b2, c2, d2 : 2nd cubic polynomial coefficients
 *
 */



/* ******************************************************************************************     1. < S | Kernel_SH | S >        */

/*      Case_1:    k1 < k2 < d          */
double SH_Integral_11_case_1( double AS, double RS, double k1, double k2, double a1, double b1, double c1, double d1, double d )
{   double res=AS*RS*pow(2*d,-1)*pow(M_E,-(d*pow(RS,-1)))*(RS*pow(d1,2)*((-k1 + RS)*pow(M_E,k1*pow(RS,-1)) + (k2 - RS)*pow(M_E,k2*pow(RS,-1))) + 
        RS*pow(d1,2)*((k2 + RS)*pow(M_E,k1*pow(RS,-1)) - (k1 + RS)*pow(M_E,k2*pow(RS,-1)))*pow(M_E,-((k1 + k2)*pow(RS,-1))) + 
        2*c1*d1*RS*(-(pow(M_E,k1*pow(RS,-1))*(-2*k1*RS + pow(k1,2) + 2*pow(RS,2))) + 
        pow(M_E,k2*pow(RS,-1))*(-2*k2*RS + pow(k2,2) + 2*pow(RS,2))) + 
        2*c1*d1*RS*pow(M_E,-((k1 + k2)*pow(RS,-1)))*(-(pow(M_E,k2*pow(RS,-1))*(2*k1*RS + pow(k1,2) + 2*pow(RS,2))) + 
        pow(M_E,k1*pow(RS,-1))*(2*k2*RS + pow(k2,2) + 2*pow(RS,2))) + 
        2*b1*d1*RS*(-(pow(M_E,k1*pow(RS,-1))*(-3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) - 6*pow(RS,3))) + 
        pow(M_E,k2*pow(RS,-1))*(-3*RS*pow(k2,2) + pow(k2,3) + 6*k2*pow(RS,2) - 6*pow(RS,3))) + 
        RS*pow(c1,2)*(-(pow(M_E,k1*pow(RS,-1))*(-3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) - 6*pow(RS,3))) + 
        pow(M_E,k2*pow(RS,-1))*(-3*RS*pow(k2,2) + pow(k2,3) + 6*k2*pow(RS,2) - 6*pow(RS,3))) + 
        2*b1*d1*RS*pow(M_E,-((k1 + k2)*pow(RS,-1)))*(-(pow(M_E,k2*pow(RS,-1))*(3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) + 6*pow(RS,3))) + 
        pow(M_E,k1*pow(RS,-1))*(3*RS*pow(k2,2) + pow(k2,3) + 6*k2*pow(RS,2) + 6*pow(RS,3))) + 
        RS*pow(c1,2)*pow(M_E,-((k1 + k2)*pow(RS,-1)))*(-(pow(M_E,k2*pow(RS,-1))*
        (3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) + 6*pow(RS,3))) + 
        pow(M_E,k1*pow(RS,-1))*(3*RS*pow(k2,2) + pow(k2,3) + 6*k2*pow(RS,2) + 6*pow(RS,3))) + 
        2*b1*c1*RS*(-(pow(M_E,k1*pow(RS,-1))*(-4*RS*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RS,2) - 24*k1*pow(RS,3) + 24*pow(RS,4))) + 
        pow(M_E,k2*pow(RS,-1))*(-4*RS*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RS,2) - 24*k2*pow(RS,3) + 24*pow(RS,4))) + 
        2*a1*d1*RS*(-(pow(M_E,k1*pow(RS,-1))*(-4*RS*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RS,2) - 24*k1*pow(RS,3) + 24*pow(RS,4))) + 
        pow(M_E,k2*pow(RS,-1))*(-4*RS*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RS,2) - 24*k2*pow(RS,3) + 24*pow(RS,4))) + 
        2*b1*c1*RS*pow(M_E,-((k1 + k2)*pow(RS,-1)))*(-(pow(M_E,k2*pow(RS,-1))*
        (4*RS*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RS,2) + 24*k1*pow(RS,3) + 24*pow(RS,4))) + 
        pow(M_E,k1*pow(RS,-1))*(4*RS*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RS,2) + 24*k2*pow(RS,3) + 24*pow(RS,4))) + 
        2*a1*d1*RS*pow(M_E,-((k1 + k2)*pow(RS,-1)))*(-(pow(M_E,k2*pow(RS,-1))*
        (4*RS*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RS,2) + 24*k1*pow(RS,3) + 24*pow(RS,4))) + 
        pow(M_E,k1*pow(RS,-1))*(4*RS*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RS,2) + 24*k2*pow(RS,3) + 24*pow(RS,4))) + 
        2*a1*c1*RS*(-(pow(M_E,k1*pow(RS,-1))*(-5*RS*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RS,2) - 60*pow(k1,2)*pow(RS,3) + 
        120*k1*pow(RS,4) - 120*pow(RS,5))) + pow(M_E,k2*pow(RS,-1))*
        (-5*RS*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RS,2) - 60*pow(k2,2)*pow(RS,3) + 120*k2*pow(RS,4) - 120*pow(RS,5))) + 
        RS*pow(b1,2)*(-(pow(M_E,k1*pow(RS,-1))*(-5*RS*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RS,2) - 60*pow(k1,2)*pow(RS,3) + 
        120*k1*pow(RS,4) - 120*pow(RS,5))) + pow(M_E,k2*pow(RS,-1))*
        (-5*RS*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RS,2) - 60*pow(k2,2)*pow(RS,3) + 120*k2*pow(RS,4) - 120*pow(RS,5))) + 
        2*a1*c1*RS*pow(M_E,-((k1 + k2)*pow(RS,-1)))*(-(pow(M_E,k2*pow(RS,-1))*
        (5*RS*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RS,2) + 60*pow(k1,2)*pow(RS,3) + 120*k1*pow(RS,4) + 120*pow(RS,5))) + 
        pow(M_E,k1*pow(RS,-1))*(5*RS*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RS,2) + 60*pow(k2,2)*pow(RS,3) + 120*k2*pow(RS,4) + 
        120*pow(RS,5))) + RS*pow(b1,2)*pow(M_E,-((k1 + k2)*pow(RS,-1)))*
        (-(pow(M_E,k2*pow(RS,-1))*(5*RS*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RS,2) + 60*pow(k1,2)*pow(RS,3) + 120*k1*pow(RS,4) + 
        120*pow(RS,5))) + pow(M_E,k1*pow(RS,-1))*(5*RS*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RS,2) + 60*pow(k2,2)*pow(RS,3) + 
        120*k2*pow(RS,4) + 120*pow(RS,5))) + 2*a1*b1*RS*
        (-(pow(M_E,k1*pow(RS,-1))*(-6*RS*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RS,2) - 120*pow(k1,3)*pow(RS,3) + 
        360*pow(k1,2)*pow(RS,4) - 720*k1*pow(RS,5) + 720*pow(RS,6))) + 
        pow(M_E,k2*pow(RS,-1))*(-6*RS*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RS,2) - 120*pow(k2,3)*pow(RS,3) + 
        360*pow(k2,2)*pow(RS,4) - 720*k2*pow(RS,5) + 720*pow(RS,6))) + 
        2*a1*b1*RS*pow(M_E,-((k1 + k2)*pow(RS,-1)))*(-(pow(M_E,k2*pow(RS,-1))*
        (6*RS*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RS,2) + 120*pow(k1,3)*pow(RS,3) + 360*pow(k1,2)*pow(RS,4) + 
        720*k1*pow(RS,5) + 720*pow(RS,6))) + pow(M_E,k1*pow(RS,-1))*
        (6*RS*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RS,2) + 120*pow(k2,3)*pow(RS,3) + 360*pow(k2,2)*pow(RS,4) + 
        720*k2*pow(RS,5) + 720*pow(RS,6))) + RS*pow(a1,2)*
        (-(pow(M_E,k1*pow(RS,-1))*(-7*RS*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RS,2) - 210*pow(k1,4)*pow(RS,3) + 
        840*pow(k1,3)*pow(RS,4) - 2520*pow(k1,2)*pow(RS,5) + 5040*k1*pow(RS,6) - 5040*pow(RS,7))) + 
        pow(M_E,k2*pow(RS,-1))*(-7*RS*pow(k2,6) + pow(k2,7) + 42*pow(k2,5)*pow(RS,2) - 210*pow(k2,4)*pow(RS,3) + 
        840*pow(k2,3)*pow(RS,4) - 2520*pow(k2,2)*pow(RS,5) + 5040*k2*pow(RS,6) - 5040*pow(RS,7))) + 
        RS*pow(a1,2)*pow(M_E,-((k1 + k2)*pow(RS,-1)))*(-(pow(M_E,k2*pow(RS,-1))*
        (7*RS*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RS,2) + 210*pow(k1,4)*pow(RS,3) + 840*pow(k1,3)*pow(RS,4) + 
        2520*pow(k1,2)*pow(RS,5) + 5040*k1*pow(RS,6) + 5040*pow(RS,7))) + 
        pow(M_E,k1*pow(RS,-1))*(7*RS*pow(k2,6) + pow(k2,7) + 42*pow(k2,5)*pow(RS,2) + 210*pow(k2,4)*pow(RS,3) + 
        840*pow(k2,3)*pow(RS,4) + 2520*pow(k2,2)*pow(RS,5) + 5040*k2*pow(RS,6) + 5040*pow(RS,7))));
    return res;
}
/* *********** CASE 1 DONE *********** */


/*      Case_2:     k1 < d < k2         */
/*      Case_2_Sub_1:   d < r < k2      */
double SH_Integral_11_case_2_sub_1( double AS, double RS, double k1, double k2, double a1, double b1, double c1, double d1, double d )
{
    double res=AS*pow(2*d,-1)*(-1 + pow(M_E,2*d*pow(RS,-1)))*pow(M_E,-((2*d + k2)*pow(RS,-1)))*pow(RS,2)*
        (-(k2*pow(d1,2)*pow(M_E,d*pow(RS,-1))) - RS*pow(d1,2)*pow(M_E,d*pow(RS,-1)) + 8*a1*d1*RS*pow(d,3)*pow(M_E,k2*pow(RS,-1)) + 
        2*a1*d1*pow(d,4)*pow(M_E,k2*pow(RS,-1)) + 7*RS*pow(a1,2)*pow(d,6)*pow(M_E,k2*pow(RS,-1)) + pow(a1,2)*pow(d,7)*pow(M_E,k2*pow(RS,-1)) + 
        d*pow(d1,2)*pow(M_E,k2*pow(RS,-1)) + RS*pow(d1,2)*pow(M_E,k2*pow(RS,-1)) - 8*a1*d1*RS*pow(M_E,d*pow(RS,-1))*pow(k2,3) - 
        2*a1*d1*pow(M_E,d*pow(RS,-1))*pow(k2,4) - 7*RS*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k2,6) - pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k2,7) + 
        24*a1*d1*pow(d,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,2) + 42*pow(a1,2)*pow(d,5)*pow(M_E,k2*pow(RS,-1))*pow(RS,2) - 
        24*a1*d1*pow(M_E,d*pow(RS,-1))*pow(k2,2)*pow(RS,2) - 42*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k2,5)*pow(RS,2) - 
        48*a1*d1*k2*pow(M_E,d*pow(RS,-1))*pow(RS,3) + 48*a1*d*d1*pow(M_E,k2*pow(RS,-1))*pow(RS,3) + 
        210*pow(a1,2)*pow(d,4)*pow(M_E,k2*pow(RS,-1))*pow(RS,3) - 210*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k2,4)*pow(RS,3) + 
        pow(c1,2)*(3*RS*pow(d,2)*pow(M_E,k2*pow(RS,-1)) + pow(d,3)*pow(M_E,k2*pow(RS,-1)) + 6*d*pow(M_E,k2*pow(RS,-1))*pow(RS,2) + 
        6*pow(M_E,k2*pow(RS,-1))*pow(RS,3) - pow(M_E,d*pow(RS,-1))*(3*RS*pow(k2,2) + pow(k2,3) + 6*k2*pow(RS,2) + 6*pow(RS,3))) - 
        48*a1*d1*pow(M_E,d*pow(RS,-1))*pow(RS,4) + 48*a1*d1*pow(M_E,k2*pow(RS,-1))*pow(RS,4) + 
        840*pow(a1,2)*pow(d,3)*pow(M_E,k2*pow(RS,-1))*pow(RS,4) - 840*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k2,3)*pow(RS,4) + 
        2520*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,5) - 2520*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(k2,2)*pow(RS,5) + 
        2*c1*(-2*d1*k2*RS*pow(M_E,d*pow(RS,-1)) + 2*d*d1*RS*pow(M_E,k2*pow(RS,-1)) + d1*pow(d,2)*pow(M_E,k2*pow(RS,-1)) + 
        5*a1*RS*pow(d,4)*pow(M_E,k2*pow(RS,-1)) + a1*pow(d,5)*pow(M_E,k2*pow(RS,-1)) - d1*pow(M_E,d*pow(RS,-1))*pow(k2,2) - 
        5*a1*RS*pow(M_E,d*pow(RS,-1))*pow(k2,4) - a1*pow(M_E,d*pow(RS,-1))*pow(k2,5) - 2*d1*pow(M_E,d*pow(RS,-1))*pow(RS,2) + 
        2*d1*pow(M_E,k2*pow(RS,-1))*pow(RS,2) + 20*a1*pow(d,3)*pow(M_E,k2*pow(RS,-1))*pow(RS,2) - 
        20*a1*pow(M_E,d*pow(RS,-1))*pow(k2,3)*pow(RS,2) + 60*a1*pow(d,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,3) - 
        60*a1*pow(M_E,d*pow(RS,-1))*pow(k2,2)*pow(RS,3) - 120*a1*k2*pow(M_E,d*pow(RS,-1))*pow(RS,4) + 
        120*a1*d*pow(M_E,k2*pow(RS,-1))*pow(RS,4) + b1*(4*RS*pow(d,3)*pow(M_E,k2*pow(RS,-1)) + pow(d,4)*pow(M_E,k2*pow(RS,-1)) + 
        12*pow(d,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,2) + 24*d*pow(M_E,k2*pow(RS,-1))*pow(RS,3) + 24*pow(M_E,k2*pow(RS,-1))*pow(RS,4) - 
        pow(M_E,d*pow(RS,-1))*(4*RS*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RS,2) + 24*k2*pow(RS,3) + 24*pow(RS,4))) - 
        120*a1*pow(M_E,d*pow(RS,-1))*pow(RS,5) + 120*a1*pow(M_E,k2*pow(RS,-1))*pow(RS,5)) + 
        pow(b1,2)*(5*RS*pow(d,4)*pow(M_E,k2*pow(RS,-1)) + pow(d,5)*pow(M_E,k2*pow(RS,-1)) + 20*pow(d,3)*pow(M_E,k2*pow(RS,-1))*pow(RS,2) + 
        60*pow(d,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,3) + 120*d*pow(M_E,k2*pow(RS,-1))*pow(RS,4) + 120*pow(M_E,k2*pow(RS,-1))*pow(RS,5) - 
        pow(M_E,d*pow(RS,-1))*(5*RS*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RS,2) + 60*pow(k2,2)*pow(RS,3) + 120*k2*pow(RS,4) + 
        120*pow(RS,5))) - 5040*k2*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(RS,6) + 5040*d*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,6) + 
        2*b1*(d1*(3*RS*pow(d,2)*pow(M_E,k2*pow(RS,-1)) + pow(d,3)*pow(M_E,k2*pow(RS,-1)) + 6*d*pow(M_E,k2*pow(RS,-1))*pow(RS,2) + 
        6*pow(M_E,k2*pow(RS,-1))*pow(RS,3) - pow(M_E,d*pow(RS,-1))*(3*RS*pow(k2,2) + pow(k2,3) + 6*k2*pow(RS,2) + 6*pow(RS,3))) + 
        a1*(6*RS*pow(d,5)*pow(M_E,k2*pow(RS,-1)) + pow(d,6)*pow(M_E,k2*pow(RS,-1)) + 30*pow(d,4)*pow(M_E,k2*pow(RS,-1))*pow(RS,2) + 
        120*pow(d,3)*pow(M_E,k2*pow(RS,-1))*pow(RS,3) + 360*pow(d,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,4) + 
        720*d*pow(M_E,k2*pow(RS,-1))*pow(RS,5) + 720*pow(M_E,k2*pow(RS,-1))*pow(RS,6) - 
        pow(M_E,d*pow(RS,-1))*(6*RS*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RS,2) + 120*pow(k2,3)*pow(RS,3) + 
        360*pow(k2,2)*pow(RS,4) + 720*k2*pow(RS,5) + 720*pow(RS,6)))) - 5040*pow(a1,2)*pow(M_E,d*pow(RS,-1))*pow(RS,7) + 
        5040*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,7));
    return res;
}
/*      Case_2_Sub_2:   k1 < r < d      */
double SH_Integral_11_case_2_sub_2( double AS, double RS, double k1, double k2, double a1, double b1, double c1, double d1, double d )
{
    double res=AS*RS*pow(2*d,-1)*pow(M_E,-(d*pow(RS,-1)))*(RS*pow(d1,2)*
        (d*pow(M_E,d*pow(RS,-1)) - RS*pow(M_E,d*pow(RS,-1)) + (-k1 + RS)*pow(M_E,k1*pow(RS,-1))) + 
        RS*pow(d1,2)*(-((k1 + RS)*pow(M_E,d*pow(RS,-1))) + d*pow(M_E,k1*pow(RS,-1)) + RS*pow(M_E,k1*pow(RS,-1)))*
        pow(M_E,-((d + k1)*pow(RS,-1))) + 2*c1*d1*RS*(-2*d*RS*pow(M_E,d*pow(RS,-1)) + pow(d,2)*pow(M_E,d*pow(RS,-1)) + 
        2*pow(M_E,d*pow(RS,-1))*pow(RS,2) - pow(M_E,k1*pow(RS,-1))*(-2*k1*RS + pow(k1,2) + 2*pow(RS,2))) + 
        2*c1*d1*RS*pow(M_E,-((d + k1)*pow(RS,-1)))*(2*d*RS*pow(M_E,k1*pow(RS,-1)) + pow(d,2)*pow(M_E,k1*pow(RS,-1)) + 
        2*pow(M_E,k1*pow(RS,-1))*pow(RS,2) - pow(M_E,d*pow(RS,-1))*(2*k1*RS + pow(k1,2) + 2*pow(RS,2))) + 
        2*b1*d1*RS*(-3*RS*pow(d,2)*pow(M_E,d*pow(RS,-1)) + pow(d,3)*pow(M_E,d*pow(RS,-1)) + 6*d*pow(M_E,d*pow(RS,-1))*pow(RS,2) - 
        pow(M_E,k1*pow(RS,-1))*(-3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) - 6*pow(RS,3)) - 6*pow(M_E,d*pow(RS,-1))*pow(RS,3)) + 
        RS*pow(c1,2)*(-3*RS*pow(d,2)*pow(M_E,d*pow(RS,-1)) + pow(d,3)*pow(M_E,d*pow(RS,-1)) + 6*d*pow(M_E,d*pow(RS,-1))*pow(RS,2) - 
        pow(M_E,k1*pow(RS,-1))*(-3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) - 6*pow(RS,3)) - 6*pow(M_E,d*pow(RS,-1))*pow(RS,3)) + 
        2*b1*d1*RS*pow(M_E,-((d + k1)*pow(RS,-1)))*(3*RS*pow(d,2)*pow(M_E,k1*pow(RS,-1)) + pow(d,3)*pow(M_E,k1*pow(RS,-1)) + 
        6*d*pow(M_E,k1*pow(RS,-1))*pow(RS,2) + 6*pow(M_E,k1*pow(RS,-1))*pow(RS,3) - 
        pow(M_E,d*pow(RS,-1))*(3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) + 6*pow(RS,3))) + 
        RS*pow(c1,2)*pow(M_E,-((d + k1)*pow(RS,-1)))*(3*RS*pow(d,2)*pow(M_E,k1*pow(RS,-1)) + pow(d,3)*pow(M_E,k1*pow(RS,-1)) + 
        6*d*pow(M_E,k1*pow(RS,-1))*pow(RS,2) + 6*pow(M_E,k1*pow(RS,-1))*pow(RS,3) - 
        pow(M_E,d*pow(RS,-1))*(3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) + 6*pow(RS,3))) + 
        2*b1*c1*RS*(-4*RS*pow(d,3)*pow(M_E,d*pow(RS,-1)) + pow(d,4)*pow(M_E,d*pow(RS,-1)) + 12*pow(d,2)*pow(M_E,d*pow(RS,-1))*pow(RS,2) - 
        24*d*pow(M_E,d*pow(RS,-1))*pow(RS,3) + 24*pow(M_E,d*pow(RS,-1))*pow(RS,4) - 
        pow(M_E,k1*pow(RS,-1))*(-4*RS*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RS,2) - 24*k1*pow(RS,3) + 24*pow(RS,4))) + 
        2*a1*d1*RS*(-4*RS*pow(d,3)*pow(M_E,d*pow(RS,-1)) + pow(d,4)*pow(M_E,d*pow(RS,-1)) + 12*pow(d,2)*pow(M_E,d*pow(RS,-1))*pow(RS,2) - 
        24*d*pow(M_E,d*pow(RS,-1))*pow(RS,3) + 24*pow(M_E,d*pow(RS,-1))*pow(RS,4) - 
        pow(M_E,k1*pow(RS,-1))*(-4*RS*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RS,2) - 24*k1*pow(RS,3) + 24*pow(RS,4))) - 
        2*b1*c1*(RS*pow(M_E,-(d*pow(RS,-1)))*(-4*RS*pow(d,3) - pow(d,4) - 12*pow(d,2)*pow(RS,2) - 24*d*pow(RS,3) + 
        24*(-1 + pow(M_E,d*pow(RS,-1)))*pow(RS,4)) + 
        RS*pow(M_E,-(k1*pow(RS,-1)))*(4*RS*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RS,2) + 24*k1*pow(RS,3) - 
        24*(-1 + pow(M_E,k1*pow(RS,-1)))*pow(RS,4))) - 
        2*a1*d1*(RS*pow(M_E,-(d*pow(RS,-1)))*(-4*RS*pow(d,3) - pow(d,4) - 12*pow(d,2)*pow(RS,2) - 24*d*pow(RS,3) + 
        24*(-1 + pow(M_E,d*pow(RS,-1)))*pow(RS,4)) + 
        RS*pow(M_E,-(k1*pow(RS,-1)))*(4*RS*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RS,2) + 24*k1*pow(RS,3) - 
        24*(-1 + pow(M_E,k1*pow(RS,-1)))*pow(RS,4))) + 
        2*a1*c1*RS*(-5*RS*pow(d,4)*pow(M_E,d*pow(RS,-1)) + pow(d,5)*pow(M_E,d*pow(RS,-1)) + 20*pow(d,3)*pow(M_E,d*pow(RS,-1))*pow(RS,2) - 
        60*pow(d,2)*pow(M_E,d*pow(RS,-1))*pow(RS,3) + 120*d*pow(M_E,d*pow(RS,-1))*pow(RS,4) - 
        pow(M_E,k1*pow(RS,-1))*(-5*RS*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RS,2) - 60*pow(k1,2)*pow(RS,3) + 120*k1*pow(RS,4) - 
        120*pow(RS,5)) - 120*pow(M_E,d*pow(RS,-1))*pow(RS,5)) + 
        RS*pow(b1,2)*(-5*RS*pow(d,4)*pow(M_E,d*pow(RS,-1)) + pow(d,5)*pow(M_E,d*pow(RS,-1)) + 20*pow(d,3)*pow(M_E,d*pow(RS,-1))*pow(RS,2) - 
        60*pow(d,2)*pow(M_E,d*pow(RS,-1))*pow(RS,3) + 120*d*pow(M_E,d*pow(RS,-1))*pow(RS,4) - 
        pow(M_E,k1*pow(RS,-1))*(-5*RS*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RS,2) - 60*pow(k1,2)*pow(RS,3) + 120*k1*pow(RS,4) - 
        120*pow(RS,5)) - 120*pow(M_E,d*pow(RS,-1))*pow(RS,5)) - 
        2*a1*c1*(RS*pow(M_E,-(d*pow(RS,-1)))*(-5*RS*pow(d,4) - pow(d,5) - 20*pow(d,3)*pow(RS,2) - 60*pow(d,2)*pow(RS,3) - 
        120*d*pow(RS,4) + 120*(-1 + pow(M_E,d*pow(RS,-1)))*pow(RS,5)) + 
        RS*pow(M_E,-(k1*pow(RS,-1)))*(5*RS*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RS,2) + 60*pow(k1,2)*pow(RS,3) + 120*k1*pow(RS,4) - 
        120*(-1 + pow(M_E,k1*pow(RS,-1)))*pow(RS,5))) - 
        pow(b1,2)*(RS*pow(M_E,-(d*pow(RS,-1)))*(-5*RS*pow(d,4) - pow(d,5) - 20*pow(d,3)*pow(RS,2) - 60*pow(d,2)*pow(RS,3) - 
        120*d*pow(RS,4) + 120*(-1 + pow(M_E,d*pow(RS,-1)))*pow(RS,5)) + 
        RS*pow(M_E,-(k1*pow(RS,-1)))*(5*RS*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RS,2) + 60*pow(k1,2)*pow(RS,3) + 120*k1*pow(RS,4) - 
        120*(-1 + pow(M_E,k1*pow(RS,-1)))*pow(RS,5))) + 
        2*a1*b1*RS*(-6*RS*pow(d,5)*pow(M_E,d*pow(RS,-1)) + pow(d,6)*pow(M_E,d*pow(RS,-1)) + 30*pow(d,4)*pow(M_E,d*pow(RS,-1))*pow(RS,2) - 
        120*pow(d,3)*pow(M_E,d*pow(RS,-1))*pow(RS,3) + 360*pow(d,2)*pow(M_E,d*pow(RS,-1))*pow(RS,4) - 
        720*d*pow(M_E,d*pow(RS,-1))*pow(RS,5) + 720*pow(M_E,d*pow(RS,-1))*pow(RS,6) - 
        pow(M_E,k1*pow(RS,-1))*(-6*RS*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RS,2) - 120*pow(k1,3)*pow(RS,3) + 
        360*pow(k1,2)*pow(RS,4) - 720*k1*pow(RS,5) + 720*pow(RS,6))) - 
        2*a1*b1*(RS*pow(M_E,-(d*pow(RS,-1)))*(-6*RS*pow(d,5) - pow(d,6) - 30*pow(d,4)*pow(RS,2) - 120*pow(d,3)*pow(RS,3) - 
        360*pow(d,2)*pow(RS,4) - 720*d*pow(RS,5) + 720*(-1 + pow(M_E,d*pow(RS,-1)))*pow(RS,6)) + 
        RS*pow(M_E,-(k1*pow(RS,-1)))*(6*RS*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RS,2) + 120*pow(k1,3)*pow(RS,3) + 
        360*pow(k1,2)*pow(RS,4) + 720*k1*pow(RS,5) - 720*(-1 + pow(M_E,k1*pow(RS,-1)))*pow(RS,6))) + 
        RS*pow(a1,2)*(-7*RS*pow(d,6)*pow(M_E,d*pow(RS,-1)) + pow(d,7)*pow(M_E,d*pow(RS,-1)) + 42*pow(d,5)*pow(M_E,d*pow(RS,-1))*pow(RS,2) - 
        210*pow(d,4)*pow(M_E,d*pow(RS,-1))*pow(RS,3) + 840*pow(d,3)*pow(M_E,d*pow(RS,-1))*pow(RS,4) - 
        2520*pow(d,2)*pow(M_E,d*pow(RS,-1))*pow(RS,5) + 5040*d*pow(M_E,d*pow(RS,-1))*pow(RS,6) - 
        pow(M_E,k1*pow(RS,-1))*(-7*RS*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RS,2) - 210*pow(k1,4)*pow(RS,3) + 
        840*pow(k1,3)*pow(RS,4) - 2520*pow(k1,2)*pow(RS,5) + 5040*k1*pow(RS,6) - 5040*pow(RS,7)) - 
        5040*pow(M_E,d*pow(RS,-1))*pow(RS,7)) - pow(a1,2)*
        (RS*pow(M_E,-(d*pow(RS,-1)))*(-7*RS*pow(d,6) - pow(d,7) - 42*pow(d,5)*pow(RS,2) - 210*pow(d,4)*pow(RS,3) - 
        840*pow(d,3)*pow(RS,4) - 2520*pow(d,2)*pow(RS,5) - 5040*d*pow(RS,6) + 5040*(-1 + pow(M_E,d*pow(RS,-1)))*pow(RS,7)) + 
        RS*pow(M_E,-(k1*pow(RS,-1)))*(7*RS*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RS,2) + 210*pow(k1,4)*pow(RS,3) + 
        840*pow(k1,3)*pow(RS,4) + 2520*pow(k1,2)*pow(RS,5) + 5040*k1*pow(RS,6) - 5040*(-1 + pow(M_E,k1*pow(RS,-1)))*pow(RS,7))));
    return res;   
}
/* *********** CASE 2 DONE *********** */


/*      Case_3:     d < k1 < k2         */
double SH_Integral_11_case_3( double AS, double RS, double k1, double k2, double a1, double b1, double c1, double d1, double d )
{   
    double res=AS*pow(2*d,-1)*(-1 + pow(M_E,2*d*pow(RS,-1)))*pow(M_E,-((d + k1 + k2)*pow(RS,-1)))*pow(RS,2)*
        (pow(d1,2)*(-((k2 + RS)*pow(M_E,k1*pow(RS,-1))) + (k1 + RS)*pow(M_E,k2*pow(RS,-1))) + 5*RS*pow(b1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,4) + 
        12*a1*b1*RS*pow(M_E,k2*pow(RS,-1))*pow(k1,5) + pow(b1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,5) + 2*a1*b1*pow(M_E,k2*pow(RS,-1))*pow(k1,6) + 
        7*RS*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,6) + pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,7) - 
        5*RS*pow(b1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,4) - 12*a1*b1*RS*pow(M_E,k1*pow(RS,-1))*pow(k2,5) - 
        pow(b1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,5) - 2*a1*b1*pow(M_E,k1*pow(RS,-1))*pow(k2,6) - 
        7*RS*pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,6) - pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,7) + 
        20*pow(b1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,3)*pow(RS,2) + 60*a1*b1*pow(M_E,k2*pow(RS,-1))*pow(k1,4)*pow(RS,2) + 
        42*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,5)*pow(RS,2) - 20*pow(b1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,3)*pow(RS,2) - 
        60*a1*b1*pow(M_E,k1*pow(RS,-1))*pow(k2,4)*pow(RS,2) - 42*pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,5)*pow(RS,2) + 
        60*pow(b1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,2)*pow(RS,3) + 240*a1*b1*pow(M_E,k2*pow(RS,-1))*pow(k1,3)*pow(RS,3) + 
        210*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,4)*pow(RS,3) - 60*pow(b1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,2)*pow(RS,3) - 
        240*a1*b1*pow(M_E,k1*pow(RS,-1))*pow(k2,3)*pow(RS,3) - 210*pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,4)*pow(RS,3) + 
        pow(c1,2)*(pow(M_E,k2*pow(RS,-1))*(3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) + 6*pow(RS,3)) - 
        pow(M_E,k1*pow(RS,-1))*(3*RS*pow(k2,2) + pow(k2,3) + 6*k2*pow(RS,2) + 6*pow(RS,3))) - 
        120*k2*pow(b1,2)*pow(M_E,k1*pow(RS,-1))*pow(RS,4) + 120*k1*pow(b1,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,4) + 
        720*a1*b1*pow(M_E,k2*pow(RS,-1))*pow(k1,2)*pow(RS,4) + 840*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,3)*pow(RS,4) - 
        720*a1*b1*pow(M_E,k1*pow(RS,-1))*pow(k2,2)*pow(RS,4) - 840*pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,3)*pow(RS,4) - 
        2*d1*(c1*(-(pow(M_E,k2*pow(RS,-1))*(2*k1*RS + pow(k1,2) + 2*pow(RS,2))) + 
        pow(M_E,k1*pow(RS,-1))*(2*k2*RS + pow(k2,2) + 2*pow(RS,2))) + 
        b1*(-(pow(M_E,k2*pow(RS,-1))*(3*RS*pow(k1,2) + pow(k1,3) + 6*k1*pow(RS,2) + 6*pow(RS,3))) + 
        pow(M_E,k1*pow(RS,-1))*(3*RS*pow(k2,2) + pow(k2,3) + 6*k2*pow(RS,2) + 6*pow(RS,3))) + 
        a1*(-(pow(M_E,k2*pow(RS,-1))*(4*RS*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RS,2) + 24*k1*pow(RS,3) + 24*pow(RS,4))) + 
        pow(M_E,k1*pow(RS,-1))*(4*RS*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RS,2) + 24*k2*pow(RS,3) + 24*pow(RS,4)))) - 
        1440*a1*b1*k2*pow(M_E,k1*pow(RS,-1))*pow(RS,5) - 120*pow(b1,2)*pow(M_E,k1*pow(RS,-1))*pow(RS,5) + 
        1440*a1*b1*k1*pow(M_E,k2*pow(RS,-1))*pow(RS,5) + 120*pow(b1,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,5) + 
        2520*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(k1,2)*pow(RS,5) - 2520*pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(k2,2)*pow(RS,5) - 
        2*c1*(b1*(-(pow(M_E,k2*pow(RS,-1))*(4*RS*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RS,2) + 24*k1*pow(RS,3) + 24*pow(RS,4))) + 
        pow(M_E,k1*pow(RS,-1))*(4*RS*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RS,2) + 24*k2*pow(RS,3) + 24*pow(RS,4))) + 
        a1*(-(pow(M_E,k2*pow(RS,-1))*(5*RS*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RS,2) + 60*pow(k1,2)*pow(RS,3) + 120*k1*pow(RS,4) + 
        120*pow(RS,5))) + pow(M_E,k1*pow(RS,-1))*
        (5*RS*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RS,2) + 60*pow(k2,2)*pow(RS,3) + 120*k2*pow(RS,4) + 120*pow(RS,5)))) - 
        1440*a1*b1*pow(M_E,k1*pow(RS,-1))*pow(RS,6) - 5040*k2*pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(RS,6) + 
        1440*a1*b1*pow(M_E,k2*pow(RS,-1))*pow(RS,6) + 5040*k1*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,6) - 
        5040*pow(a1,2)*pow(M_E,k1*pow(RS,-1))*pow(RS,7) + 5040*pow(a1,2)*pow(M_E,k2*pow(RS,-1))*pow(RS,7));
    return res;
}
/* *********** CASE 3 DONE *********** */

/* ******************************************************************************************     1. < S | Kernel_SH | S >        */


        /*      *       *       *       *       *       *       *       *       *       *       *       *       */


/* ******************************************************************************************     2. < S | Kernel_SH | Z >        */

/*      Case_1:    k1 < k2 < d          */
double SH_Integral_14_case_1( double ASP, double RSP, double k1, double k2, double a1, double b1, double c1, double d1,
        double a2, double b2, double c2, double d2, double d )
{   
    double res=-(ASP*(d + RSP)*pow(3,0.5)*pow(M_E,-((d + k1 + k2)*pow(RSP,-1)))*pow(RSP,2)*
        (-4*b1*d2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,2) - 4*b1*d2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2) - 
        b1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) - 5*b1*c2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) - 
        5*a1*d2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + b1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3) - 
        5*b1*c2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3) - 5*a1*d2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3) - 
        b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) - a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) - 6*b1*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) - 
        6*a1*c2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + b1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) + 
        a1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) - 6*b1*b2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) - 
        6*a1*c2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4) - b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) - 
        a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) - 7*a2*b1*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) - 
        7*a1*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + b1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) + 
        a1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) - 7*a2*b1*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) - 
        7*a1*b2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5) - a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) - 
        a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) - 8*a1*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) + 
        a2*b1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,6) + a1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,6) - 
        8*a1*a2*RSP*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,6) - a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,7) + 
        a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,7) + 4*b1*d2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,2) + 
        4*b1*d2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2) + b1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) + 
        5*b1*c2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) + 5*a1*d2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - 
        b1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3) + 5*b1*c2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3) + 
        5*a1*d2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3) + b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) + 
        a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) + 6*b1*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) + 
        6*a1*c2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - b1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) - 
        a1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) + 6*b1*b2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) + 
        6*a1*c2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4) + b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) + 
        a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) + 7*a2*b1*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) + 
        7*a1*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - b1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) - 
        a1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) + 7*a2*b1*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) + 
        7*a1*b2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5) + a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) + 
        a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) + 8*a1*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) - 
        a2*b1*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,6) - a1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,6) + 
        8*a1*a2*RSP*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,6) + a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,7) - 
        a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,7) + 8*b1*d2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 
        8*b1*d2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 8*b1*d2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,2) - 
        8*b1*d2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,2) - 15*b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 
        15*a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 15*b1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 
        15*a1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 24*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 
        24*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 24*b1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 
        24*a1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 35*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) - 
        35*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 35*a2*b1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 
        35*a1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) - 48*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) + 
        48*a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) + 15*b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 
        15*a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 15*b1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 
        15*a1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 24*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 
        24*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 24*b1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 
        24*a1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 35*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) + 
        35*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) - 35*a2*b1*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) - 
        35*a1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) + 48*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) - 
        48*a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) + 8*b1*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) + 
        30*b1*c2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) + 30*a1*d2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
        8*b1*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 30*b1*c2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 
        30*a1*d2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 8*b1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,3) - 
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
        240*a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) + 30*b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 
        30*a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 144*b1*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 
        144*a1*c2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 30*b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 
        30*a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 144*b1*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 
        144*a1*c2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 30*b1*c2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) + 
        30*a1*d2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) + 144*b1*b2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) + 
        144*a1*c2*k1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,4) - 30*b1*c2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) - 
        30*a1*d2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) - 144*b1*b2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) - 
        144*a1*c2*k2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,4) - 420*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) - 
        420*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 420*a2*b1*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 
        420*a1*b2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) - 960*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) + 
        960*a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) + 420*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 
        420*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) - 420*a2*b1*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) - 
        420*a1*b2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 960*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) - 
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
        2880*a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) + 2880*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) + 
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
        5760*a1*a2*pow(M_E,(2*k1 + k2)*pow(RSP,-1))*pow(RSP,7) + 5760*a1*a2*pow(M_E,(k1 + 2*k2)*pow(RSP,-1))*pow(RSP,7))*pow(2*pow(d,2),-1));
    return res;
}
/* *********** CASE 1 DONE *********** */

/*      Case_2:     k1 < d < k2         */
/*      Case_2_Sub_1:   d < r < k2      */
double SH_Integral_14_case_2_sub_1( double ASP, double RSP, double k1, double k2, double a1, double b1, double c1, double d1,
        double a2, double b2, double c2, double d2, double d )
{
    double res=ASP*pow(3,0.5)*(d + RSP + d*pow(M_E,2*d*pow(RSP,-1)) - RSP*pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,-((2*d + k2)*pow(RSP,-1)))*pow(RSP,2)*
        (-(d1*d2*k2*pow(M_E,d*pow(RSP,-1))) - 2*d1*d2*RSP*pow(M_E,d*pow(RSP,-1)) - 3*c2*d1*k2*RSP*pow(M_E,d*pow(RSP,-1)) + 
        d*d1*d2*pow(M_E,k2*pow(RSP,-1)) + 3*c2*d*d1*RSP*pow(M_E,k2*pow(RSP,-1)) + 2*d1*d2*RSP*pow(M_E,k2*pow(RSP,-1)) + 
        c2*d1*pow(d,2)*pow(M_E,k2*pow(RSP,-1)) + 4*b2*d1*RSP*pow(d,2)*pow(M_E,k2*pow(RSP,-1)) + b2*d1*pow(d,3)*pow(M_E,k2*pow(RSP,-1)) + 
        5*a2*d1*RSP*pow(d,3)*pow(M_E,k2*pow(RSP,-1)) + 5*a1*d2*RSP*pow(d,3)*pow(M_E,k2*pow(RSP,-1)) + a2*d1*pow(d,4)*pow(M_E,k2*pow(RSP,-1)) + 
        a1*d2*pow(d,4)*pow(M_E,k2*pow(RSP,-1)) + 6*a1*c2*RSP*pow(d,4)*pow(M_E,k2*pow(RSP,-1)) + a1*c2*pow(d,5)*pow(M_E,k2*pow(RSP,-1)) + 
        7*a1*b2*RSP*pow(d,5)*pow(M_E,k2*pow(RSP,-1)) + a1*b2*pow(d,6)*pow(M_E,k2*pow(RSP,-1)) + 8*a1*a2*RSP*pow(d,6)*pow(M_E,k2*pow(RSP,-1)) + 
        a1*a2*pow(d,7)*pow(M_E,k2*pow(RSP,-1)) - c2*d1*pow(M_E,d*pow(RSP,-1))*pow(k2,2) - 4*b2*d1*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,2) - 
        b2*d1*pow(M_E,d*pow(RSP,-1))*pow(k2,3) - 5*a2*d1*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,3) - 5*a1*d2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,3) - 
        a2*d1*pow(M_E,d*pow(RSP,-1))*pow(k2,4) - a1*d2*pow(M_E,d*pow(RSP,-1))*pow(k2,4) - 6*a1*c2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,4) - 
        a1*c2*pow(M_E,d*pow(RSP,-1))*pow(k2,5) - 7*a1*b2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,5) - a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k2,6) - 
        8*a1*a2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,6) - a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,7) - 3*c2*d1*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) - 
        8*b2*d1*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) + 3*c2*d1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 
        8*b2*d*d1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 15*a2*d1*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 
        15*a1*d2*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 24*a1*c2*pow(d,3)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 
        35*a1*b2*pow(d,4)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 48*a1*a2*pow(d,5)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) - 
        15*a2*d1*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 15*a1*d2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 
        24*a1*c2*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 35*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) - 
        48*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,5)*pow(RSP,2) - 8*b2*d1*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) - 
        30*a2*d1*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) - 30*a1*d2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 
        8*b2*d1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 30*a2*d*d1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 
        30*a1*d*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 72*a1*c2*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 
        140*a1*b2*pow(d,3)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 240*a1*a2*pow(d,4)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 
        72*a1*c2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 140*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) - 
        240*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,4)*pow(RSP,3) - 30*a2*d1*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) - 
        30*a1*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) - 144*a1*c2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 
        144*a1*c2*d*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 30*a2*d1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 
        30*a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 420*a1*b2*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 
        960*a1*a2*pow(d,3)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 420*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) - 
        960*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,4) - 144*a1*c2*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) - 
        840*a1*b2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 144*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 
        840*a1*b2*d*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 2880*a1*a2*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 
        2880*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) + 
        c1*(-3*d2*k2*RSP*pow(M_E,d*pow(RSP,-1)) + 3*d*d2*RSP*pow(M_E,k2*pow(RSP,-1)) + d2*pow(d,2)*pow(M_E,k2*pow(RSP,-1)) + 
        5*b2*RSP*pow(d,3)*pow(M_E,k2*pow(RSP,-1)) + b2*pow(d,4)*pow(M_E,k2*pow(RSP,-1)) + 6*a2*RSP*pow(d,4)*pow(M_E,k2*pow(RSP,-1)) + 
        a2*pow(d,5)*pow(M_E,k2*pow(RSP,-1)) - d2*pow(M_E,d*pow(RSP,-1))*pow(k2,2) - 5*b2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,3) - 
        b2*pow(M_E,d*pow(RSP,-1))*pow(k2,4) - 6*a2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,4) - a2*pow(M_E,d*pow(RSP,-1))*pow(k2,5) - 
        3*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) + 3*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 
        15*b2*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 24*a2*pow(d,3)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) - 
        15*b2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 24*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 
        30*b2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 30*b2*d*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 
        72*a2*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 72*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 
        c2*(4*RSP*pow(d,2)*pow(M_E,k2*pow(RSP,-1)) + pow(d,3)*pow(M_E,k2*pow(RSP,-1)) + 8*d*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 
        8*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - pow(M_E,d*pow(RSP,-1))*(4*RSP*pow(k2,2) + pow(k2,3) + 8*k2*pow(RSP,2) + 8*pow(RSP,3)))\
        - 30*b2*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) - 144*a2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 
        30*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 144*a2*d*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 144*a2*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 
        144*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5)) - 840*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) - 
        5760*a1*a2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) + 840*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) + 
        5760*a1*a2*d*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) + b1*
        (4*d2*RSP*pow(d,2)*pow(M_E,k2*pow(RSP,-1)) + d2*pow(d,3)*pow(M_E,k2*pow(RSP,-1)) + 6*b2*RSP*pow(d,4)*pow(M_E,k2*pow(RSP,-1)) + 
        b2*pow(d,5)*pow(M_E,k2*pow(RSP,-1)) + 7*a2*RSP*pow(d,5)*pow(M_E,k2*pow(RSP,-1)) + a2*pow(d,6)*pow(M_E,k2*pow(RSP,-1)) - 
        4*d2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,2) - d2*pow(M_E,d*pow(RSP,-1))*pow(k2,3) - 6*b2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,4) - 
        b2*pow(M_E,d*pow(RSP,-1))*pow(k2,5) - 7*a2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k2,5) - a2*pow(M_E,d*pow(RSP,-1))*pow(k2,6) - 
        8*d2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) + 8*d*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 
        24*b2*pow(d,3)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 35*a2*pow(d,4)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) - 
        24*b2*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) - 35*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,4)*pow(RSP,2) - 
        8*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 8*d2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 
        72*b2*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 140*a2*pow(d,3)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 
        72*b2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) - 140*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,3)*pow(RSP,3) - 
        144*b2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 144*b2*d*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 
        420*a2*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 420*a2*pow(M_E,d*pow(RSP,-1))*pow(k2,2)*pow(RSP,4) + 
        c2*(5*RSP*pow(d,3)*pow(M_E,k2*pow(RSP,-1)) + pow(d,4)*pow(M_E,k2*pow(RSP,-1)) + 15*pow(d,2)*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 
        30*d*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 30*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 
        pow(M_E,d*pow(RSP,-1))*(5*RSP*pow(k2,3) + pow(k2,4) + 15*pow(k2,2)*pow(RSP,2) + 30*k2*pow(RSP,3) + 30*pow(RSP,4))) - 
        144*b2*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) - 840*a2*k2*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 
        144*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 840*a2*d*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) - 
        840*a2*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) + 840*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6)) - 
        5760*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(RSP,7) + 5760*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,7))*pow(2*pow(d,2),-1);
    return res;
}
/*      Case_2_Sub_2:   k1 < r < d      */
double SH_Integral_14_case_2_sub_2( double ASP, double RSP, double k1, double k2, double a1, double b1, double c1, double d1,
        double a2, double b2, double c2, double d2, double d )
{
    double res=ASP*(d + RSP)*pow(3,0.5)*pow(M_E,-((2*d + k1)*pow(RSP,-1)))*pow(RSP,2)*
        (d1*d2*k1*pow(M_E,d*pow(RSP,-1)) + 2*d1*d2*RSP*pow(M_E,d*pow(RSP,-1)) + 3*c2*d1*k1*RSP*pow(M_E,d*pow(RSP,-1)) - 
        d*d1*d2*pow(M_E,k1*pow(RSP,-1)) - 3*c2*d*d1*RSP*pow(M_E,k1*pow(RSP,-1)) - 2*d1*d2*RSP*pow(M_E,k1*pow(RSP,-1)) - 
        c2*d1*pow(d,2)*pow(M_E,k1*pow(RSP,-1)) - 4*b2*d1*RSP*pow(d,2)*pow(M_E,k1*pow(RSP,-1)) - b2*d1*pow(d,3)*pow(M_E,k1*pow(RSP,-1)) - 
        5*a2*d1*RSP*pow(d,3)*pow(M_E,k1*pow(RSP,-1)) - 5*a1*d2*RSP*pow(d,3)*pow(M_E,k1*pow(RSP,-1)) - a2*d1*pow(d,4)*pow(M_E,k1*pow(RSP,-1)) - 
        a1*d2*pow(d,4)*pow(M_E,k1*pow(RSP,-1)) - 6*a1*c2*RSP*pow(d,4)*pow(M_E,k1*pow(RSP,-1)) - a1*c2*pow(d,5)*pow(M_E,k1*pow(RSP,-1)) - 
        7*a1*b2*RSP*pow(d,5)*pow(M_E,k1*pow(RSP,-1)) - a1*b2*pow(d,6)*pow(M_E,k1*pow(RSP,-1)) - 8*a1*a2*RSP*pow(d,6)*pow(M_E,k1*pow(RSP,-1)) - 
        a1*a2*pow(d,7)*pow(M_E,k1*pow(RSP,-1)) + d*d1*d2*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 3*c2*d*d1*RSP*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 
        2*d1*d2*RSP*pow(M_E,(2*d + k1)*pow(RSP,-1)) + c2*d1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 
        4*b2*d1*RSP*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + b2*d1*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 
        5*a2*d1*RSP*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 5*a1*d2*RSP*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + 
        a2*d1*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + a1*d2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 
        6*a1*c2*RSP*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + a1*c2*pow(d,5)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 
        7*a1*b2*RSP*pow(d,5)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + a1*b2*pow(d,6)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 
        8*a1*a2*RSP*pow(d,6)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + a1*a2*pow(d,7)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 
        d1*d2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1)) + 2*d1*d2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1)) + 
        3*c2*d1*k1*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1)) + c2*d1*pow(M_E,d*pow(RSP,-1))*pow(k1,2) + 
        4*b2*d1*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,2) - c2*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2) + 
        4*b2*d1*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2) + b2*d1*pow(M_E,d*pow(RSP,-1))*pow(k1,3) + 
        5*a2*d1*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,3) + 5*a1*d2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,3) - 
        b2*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3) + 5*a2*d1*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3) + 
        5*a1*d2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3) + a2*d1*pow(M_E,d*pow(RSP,-1))*pow(k1,4) + 
        a1*d2*pow(M_E,d*pow(RSP,-1))*pow(k1,4) + 6*a1*c2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,4) - 
        a2*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4) - a1*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4) + 
        6*a1*c2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4) + a1*c2*pow(M_E,d*pow(RSP,-1))*pow(k1,5) + 
        7*a1*b2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,5) - a1*c2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5) + 
        7*a1*b2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5) + a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k1,6) + 
        8*a1*a2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,6) - a1*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,6) + 
        8*a1*a2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,6) + a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,7) - 
        a1*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,7) + 3*c2*d1*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) + 
        8*b2*d1*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) - 3*c2*d1*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 
        8*b2*d*d1*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 15*a2*d1*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 
        15*a1*d2*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 24*a1*c2*pow(d,3)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 
        35*a1*b2*pow(d,4)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 48*a1*a2*pow(d,5)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) + 
        3*c2*d1*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) + 8*b2*d*d1*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) + 
        15*a2*d1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) + 15*a1*d2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) + 
        24*a1*c2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) + 35*a1*b2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) + 
        48*a1*a2*pow(d,5)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) - 3*c2*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,2) - 
        8*b2*d1*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,2) + 15*a2*d1*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 
        15*a1*d2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 15*a2*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 
        15*a1*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 24*a1*c2*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 
        24*a1*c2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 35*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) - 
        35*a1*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 48*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) - 
        48*a1*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5)*pow(RSP,2) + 8*b2*d1*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 
        30*a2*d1*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) + 30*a1*d2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) - 
        8*b2*d1*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 30*a2*d*d1*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
        30*a1*d*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 72*a1*c2*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
        140*a1*b2*pow(d,3)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 240*a1*a2*pow(d,4)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
        8*b2*d1*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) - 30*a2*d*d1*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) - 
        30*a1*d*d2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) - 72*a1*c2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) - 
        140*a1*b2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) - 240*a1*a2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) + 
        8*b2*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,3) + 30*a2*d1*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,3) + 
        30*a1*d2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,3) + 72*a1*c2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
        72*a1*c2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 140*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 
        140*a1*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 240*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) + 
        240*a1*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4)*pow(RSP,3) + 30*a2*d1*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 
        30*a1*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 144*a1*c2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) - 
        144*a1*c2*d*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 30*a2*d1*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 
        30*a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 420*a1*b2*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 
        960*a1*a2*pow(d,3)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 144*a1*c2*d*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) + 
        30*a2*d1*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) + 30*a1*d2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) + 
        420*a1*b2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) + 960*a1*a2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) - 
        30*a2*d1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,4) - 30*a1*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,4) - 
        144*a1*c2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,4) + 420*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) - 
        420*a1*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 960*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) - 
        960*a1*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,4) + 144*a1*c2*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 
        840*a1*b2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) - 144*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 
        840*a1*b2*d*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 2880*a1*a2*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 
        144*a1*c2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,5) - 840*a1*b2*d*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,5) - 
        2880*a1*a2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,5) + 144*a1*c2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,5) + 
        840*a1*b2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,5) + 2880*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) + 
        2880*a1*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) + 
        c1*(3*d2*k1*RSP*pow(M_E,d*pow(RSP,-1)) - 3*d*d2*RSP*pow(M_E,k1*pow(RSP,-1)) - d2*pow(d,2)*pow(M_E,k1*pow(RSP,-1)) - 
        5*b2*RSP*pow(d,3)*pow(M_E,k1*pow(RSP,-1)) - b2*pow(d,4)*pow(M_E,k1*pow(RSP,-1)) - 6*a2*RSP*pow(d,4)*pow(M_E,k1*pow(RSP,-1)) - 
        a2*pow(d,5)*pow(M_E,k1*pow(RSP,-1)) - 3*d*d2*RSP*pow(M_E,(2*d + k1)*pow(RSP,-1)) + d2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 
        5*b2*RSP*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + b2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 
        6*a2*RSP*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + a2*pow(d,5)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + 
        3*d2*k1*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1)) + d2*pow(M_E,d*pow(RSP,-1))*pow(k1,2) - d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2) + 
        5*b2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,3) + 5*b2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3) + 
        b2*pow(M_E,d*pow(RSP,-1))*pow(k1,4) + 6*a2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,4) - b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4) + 
        6*a2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4) + a2*pow(M_E,d*pow(RSP,-1))*pow(k1,5) - 
        a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5) + 3*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) - 3*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 
        15*b2*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 24*a2*pow(d,3)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) + 
        3*d2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) + 15*b2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) + 
        24*a2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) - 3*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,2) + 
        15*b2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 15*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 
        24*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 24*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 
        30*b2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) - 30*b2*d*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
        72*a2*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 30*b2*d*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) - 
        72*a2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) + 30*b2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,3) + 
        72*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 72*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
        c2*(pow(d,3)*(-1 + pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1)) - 
        4*RSP*pow(d,2)*(1 + pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1)) + 
        8*d*(-1 + pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 8*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
        8*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) + 
        pow(M_E,(d + 2*k1)*pow(RSP,-1))*(4*RSP*pow(k1,2) - pow(k1,3) - 8*k1*pow(RSP,2) + 8*pow(RSP,3)) + 
        pow(M_E,d*pow(RSP,-1))*(4*RSP*pow(k1,2) + pow(k1,3) + 8*k1*pow(RSP,2) + 8*pow(RSP,3))) + 
        30*b2*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) + 144*a2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) - 30*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 
        144*a2*d*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 30*b2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) + 
        144*a2*d*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) - 30*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,4) - 
        144*a2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,4) + 144*a2*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) - 
        144*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 144*a2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,5) + 
        144*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,5)) + 840*a1*b2*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) + 
        5760*a1*a2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) - 840*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) - 
        5760*a1*a2*d*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) + 840*a1*b2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,6) + 
        5760*a1*a2*d*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,6) - 840*a1*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,6) - 
        5760*a1*a2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,6) + 
        b1*(-4*d2*RSP*pow(d,2)*pow(M_E,k1*pow(RSP,-1)) - d2*pow(d,3)*pow(M_E,k1*pow(RSP,-1)) - 6*b2*RSP*pow(d,4)*pow(M_E,k1*pow(RSP,-1)) - 
        b2*pow(d,5)*pow(M_E,k1*pow(RSP,-1)) - 7*a2*RSP*pow(d,5)*pow(M_E,k1*pow(RSP,-1)) - a2*pow(d,6)*pow(M_E,k1*pow(RSP,-1)) - 
        4*d2*RSP*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + d2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 
        6*b2*RSP*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + b2*pow(d,5)*pow(M_E,(2*d + k1)*pow(RSP,-1)) - 
        7*a2*RSP*pow(d,5)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + a2*pow(d,6)*pow(M_E,(2*d + k1)*pow(RSP,-1)) + 
        4*d2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,2) + 4*d2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2) + 
        d2*pow(M_E,d*pow(RSP,-1))*pow(k1,3) - d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3) + 6*b2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,4) + 
        6*b2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4) + b2*pow(M_E,d*pow(RSP,-1))*pow(k1,5) + 
        7*a2*RSP*pow(M_E,d*pow(RSP,-1))*pow(k1,5) - b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5) + 
        7*a2*RSP*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,5) + a2*pow(M_E,d*pow(RSP,-1))*pow(k1,6) - 
        a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,6) + 8*d2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,2) - 
        8*d*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 24*b2*pow(d,3)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 
        35*a2*pow(d,4)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) + 8*d*d2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) + 
        24*b2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) + 35*a2*pow(d,4)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,2) - 
        8*d2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,2) + 24*b2*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 
        24*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) + 35*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) - 
        35*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,4)*pow(RSP,2) + 8*d2*pow(M_E,d*pow(RSP,-1))*pow(RSP,3) - 
        8*d2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 72*b2*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
        140*a2*pow(d,3)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 8*d2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) - 
        72*b2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) - 140*a2*pow(d,3)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,3) + 
        8*d2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,3) + 72*b2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 
        72*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) + 140*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 
        140*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,3)*pow(RSP,3) + 144*b2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,4) - 
        144*b2*d*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 420*a2*pow(d,2)*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 
        144*b2*d*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) + 420*a2*pow(d,2)*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) - 
        144*b2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,4) + 420*a2*pow(M_E,d*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) - 
        420*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(k1,2)*pow(RSP,4) + 
        c2*(pow(d,4)*(-1 + pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1)) - 
        5*RSP*pow(d,3)*(1 + pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1)) + 
        15*pow(d,2)*(-1 + pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) - 
        30*d*(1 + pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 30*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 
        30*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,4) - 
        pow(M_E,(d + 2*k1)*pow(RSP,-1))*(-5*RSP*pow(k1,3) + pow(k1,4) + 15*pow(k1,2)*pow(RSP,2) - 30*k1*pow(RSP,3) + 
        30*pow(RSP,4)) + pow(M_E,d*pow(RSP,-1))*(5*RSP*pow(k1,3) + pow(k1,4) + 15*pow(k1,2)*pow(RSP,2) + 30*k1*pow(RSP,3) + 
        30*pow(RSP,4))) + 144*b2*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) + 840*a2*k1*pow(M_E,d*pow(RSP,-1))*pow(RSP,5) - 
        144*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 840*a2*d*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 
        144*b2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,5) - 840*a2*d*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,5) + 
        144*b2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,5) + 840*a2*k1*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,5) + 
        840*a2*pow(M_E,d*pow(RSP,-1))*pow(RSP,6) - 840*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) + 
        840*a2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,6) - 840*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,6)) + 
        5760*a1*a2*pow(M_E,d*pow(RSP,-1))*pow(RSP,7) - 5760*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,7) - 
        5760*a1*a2*pow(M_E,(2*d + k1)*pow(RSP,-1))*pow(RSP,7) + 5760*a1*a2*pow(M_E,(d + 2*k1)*pow(RSP,-1))*pow(RSP,7))*pow(2*pow(d,2),-1);
return res;
}
/* *********** CASE 2 DONE *********** */

/*      Case_3:     d < k1 < k2         */
double SH_Integral_14_case_3( double ASP, double RSP, double k1, double k2, double a1, double b1, double c1, double d1,
        double a2, double b2, double c2, double d2, double d )
{
    double res=ASP*pow(3,0.5)*(d + RSP + d*pow(M_E,2*d*pow(RSP,-1)) - RSP*pow(M_E,2*d*pow(RSP,-1)))*pow(M_E,-((d + k1 + k2)*pow(RSP,-1)))*pow(RSP,2)*
        (4*b1*d2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,2) + b1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + 
        5*b1*c2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + 5*a1*d2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + 
        b1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + a1*d2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + 6*b1*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + 
        6*a1*c2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + 
        7*a2*b1*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + 7*a1*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) + 
        a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) + a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) + 8*a1*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,6) + 
        a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,7) - 4*b1*d2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,2) - b1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - 
        5*b1*c2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - 5*a1*d2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - 
        b1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - a1*d2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 6*b1*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 
        6*a1*c2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - 
        7*a2*b1*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - 7*a1*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) - 
        a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) - a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) - 8*a1*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,6) - 
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
        d1*(d2*(-((k2 + 2*RSP)*pow(M_E,k1*pow(RSP,-1))) + (k1 + 2*RSP)*pow(M_E,k2*pow(RSP,-1))) + 4*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,2) + 
        b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + 5*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) - 
        4*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,2) - b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - 5*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - 
        a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 8*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,2) + 8*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,2) + 
        15*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) - 15*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) + 
        c2*(pow(M_E,k2*pow(RSP,-1))*(3*k1*RSP + pow(k1,2) + 3*pow(RSP,2)) - 
        pow(M_E,k1*pow(RSP,-1))*(3*k2*RSP + pow(k2,2) + 3*pow(RSP,2))) - 8*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) - 
        30*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) + 8*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 
        30*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) - 30*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 30*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4))\
        - 144*b1*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 144*a1*c2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 
        840*a2*b1*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) - 840*a1*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 
        144*b1*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 144*a1*c2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 
        840*a2*b1*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 840*a1*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5) + 
        2880*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,5) - 2880*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,5) + 
        c1*(5*b2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,3) + b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + 6*a2*RSP*pow(M_E,k2*pow(RSP,-1))*pow(k1,4) + 
        a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,5) - 5*b2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,3) - b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - 
        6*a2*RSP*pow(M_E,k1*pow(RSP,-1))*pow(k2,4) - a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,5) + 
        15*b2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,2) + 24*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,3)*pow(RSP,2) - 
        15*b2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,2) - 24*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,3)*pow(RSP,2) + 
        d2*(pow(M_E,k2*pow(RSP,-1))*(3*k1*RSP + pow(k1,2) + 3*pow(RSP,2)) - 
        pow(M_E,k1*pow(RSP,-1))*(3*k2*RSP + pow(k2,2) + 3*pow(RSP,2))) - 30*b2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,3) + 
        30*b2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,3) + 72*a2*pow(M_E,k2*pow(RSP,-1))*pow(k1,2)*pow(RSP,3) - 
        72*a2*pow(M_E,k1*pow(RSP,-1))*pow(k2,2)*pow(RSP,3) + 
        c2*(pow(M_E,k2*pow(RSP,-1))*(4*RSP*pow(k1,2) + pow(k1,3) + 8*k1*pow(RSP,2) + 8*pow(RSP,3)) - 
        pow(M_E,k1*pow(RSP,-1))*(4*RSP*pow(k2,2) + pow(k2,3) + 8*k2*pow(RSP,2) + 8*pow(RSP,3))) - 
        30*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) - 144*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,4) + 
        30*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) + 144*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,4) - 
        144*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,5) + 144*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,5)) - 
        840*a2*b1*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) - 840*a1*b2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) - 
        5760*a1*a2*k2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,6) + 840*a2*b1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) + 
        840*a1*b2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) + 5760*a1*a2*k1*pow(M_E,k2*pow(RSP,-1))*pow(RSP,6) - 
        5760*a1*a2*pow(M_E,k1*pow(RSP,-1))*pow(RSP,7) + 5760*a1*a2*pow(M_E,k2*pow(RSP,-1))*pow(RSP,7))*pow(2*pow(d,2),-1);
    return res;
}
/* *********** CASE 3 DONE *********** */

/* ******************************************************************************************     2. < S | Kernel_SH | Z >        */


/*      *       *       *       *       *       *       *       *       *       *       *       *       */


/* ******************************************************************************************     3. < XorY | Kernel_SH | XorY >        */

/*      Case_1:    k1 < k2 < d          */
double SH_Integral_2233_case_1( double AP, double RP, double k1, double k2, double a1, double b1, double c1, double d1, double d )
{
    double res=3*AP*pow(M_E,-((d + k1 + k2)*pow(RP,-1)))*pow(RP,3)*(-10*c1*d*d1*k2*RP*pow(M_E,k1*pow(RP,-1)) - 2*c1*d1*k2*pow(d,2)*pow(M_E,k1*pow(RP,-1)) - 
        4*c1*d1*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1)) - 6*b1*d1*k2*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1)) - 3*k2*RP*pow(c1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1)) - 
        d*k2*pow(d1,2)*pow(M_E,k1*pow(RP,-1)) - 4*d*RP*pow(d1,2)*pow(M_E,k1*pow(RP,-1)) - k2*RP*pow(d1,2)*pow(M_E,k1*pow(RP,-1)) - 
        pow(d,2)*pow(d1,2)*pow(M_E,k1*pow(RP,-1)) + 10*c1*d*d1*k1*RP*pow(M_E,k2*pow(RP,-1)) + 2*c1*d1*k1*pow(d,2)*pow(M_E,k2*pow(RP,-1)) + 
        4*c1*d1*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1)) + 6*b1*d1*k1*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1)) + 3*k1*RP*pow(c1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1)) + 
        d*k1*pow(d1,2)*pow(M_E,k2*pow(RP,-1)) + 4*d*RP*pow(d1,2)*pow(M_E,k2*pow(RP,-1)) + k1*RP*pow(d1,2)*pow(M_E,k2*pow(RP,-1)) + 
        pow(d,2)*pow(d1,2)*pow(M_E,k2*pow(RP,-1)) + 3*d*RP*gsl_sf_expint_Ei(-(k2*pow(RP,-1)))*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1)) - 
        3*d*RP*gsl_sf_expint_Ei(k2*pow(RP,-1))*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1)) + gsl_sf_expint_Ei(-(k2*pow(RP,-1)))*pow(d,2)*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1)) - 
        gsl_sf_expint_Ei(k2*pow(RP,-1))*pow(d,2)*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1)) - 10*c1*d*d1*k1*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1)) - 
        2*c1*d1*k1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) + 4*c1*d1*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) + 
        6*b1*d1*k1*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) + 3*k1*RP*pow(c1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) + 
        d*k1*pow(d1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) - 4*d*RP*pow(d1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) + k1*RP*pow(d1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) - 
        pow(d,2)*pow(d1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) + 10*c1*d*d1*k2*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) + 2*c1*d1*k2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) - 
        4*c1*d1*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) - 6*b1*d1*k2*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) - 
        3*k2*RP*pow(c1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) - d*k2*pow(d1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) + 4*d*RP*pow(d1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) - 
        k2*RP*pow(d1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) + pow(d,2)*pow(d1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) + 2*c1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 
        2*c1*d1*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 12*b1*d*d1*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 6*d*RP*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 
        2*b1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 8*b1*c1*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 8*a1*d1*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 
        pow(c1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 2*c1*d*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) + 2*c1*d1*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) - 
        12*b1*d*d1*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) - 6*d*RP*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) - 
        2*b1*d1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) + 8*b1*c1*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) + 
        8*a1*d1*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) - pow(c1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) + 
        2*b1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 14*b1*c1*d*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 2*b1*d1*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 
        14*a1*d*d1*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + d*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + RP*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 
        2*b1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 2*a1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 10*a1*c1*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 
        5*RP*pow(b1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 2*b1*d*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) - 
        14*b1*c1*d*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + 2*b1*d1*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) - 
        14*a1*d*d1*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + d*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + 
        RP*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) - 2*b1*c1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) - 
        2*a1*d1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + 10*a1*c1*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + 
        5*RP*pow(b1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + 2*b1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 2*a1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 
        2*b1*c1*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 16*a1*c1*d*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 2*a1*d1*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 
        8*d*RP*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 2*a1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 12*a1*b1*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 
        pow(b1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 2*b1*c1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) + 2*a1*d*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) + 
        2*b1*c1*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) - 16*a1*c1*d*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) + 
        2*a1*d1*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) - 8*d*RP*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) - 
        2*a1*c1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) + 12*a1*b1*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) - 
        pow(b1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) + 2*a1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 2*a1*c1*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 
        18*a1*b1*d*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + d*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + RP*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 
        2*a1*b1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 7*RP*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 
        2*a1*c1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) + 2*a1*c1*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) - 
        18*a1*b1*d*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) + d*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) + 
        RP*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) - 2*a1*b1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) + 
        7*RP*pow(a1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) + 2*a1*b1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,6) + 2*a1*b1*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,6) + 
        10*d*RP*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,6) + pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,6) + 2*a1*b1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6) + 
        2*a1*b1*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6) - 10*d*RP*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6) - 
        pow(a1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6) + d*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,7) + RP*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,7) + 
        d*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,7) + RP*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,7) - 2*c1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 
        2*c1*d1*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 12*b1*d*d1*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 6*d*RP*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 
        2*b1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 8*b1*c1*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 8*a1*d1*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 
        pow(c1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 2*c1*d*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) - 2*c1*d1*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) + 
        12*b1*d*d1*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) + 6*d*RP*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) + 
        2*b1*d1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) - 8*b1*c1*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) - 
        8*a1*d1*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) + pow(c1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) - 
        2*b1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 14*b1*c1*d*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 2*b1*d1*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 
        14*a1*d*d1*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - d*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - RP*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 
        2*b1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 2*a1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 10*a1*c1*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 
        5*RP*pow(b1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 2*b1*d*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) + 
        14*b1*c1*d*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - 2*b1*d1*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) + 
        14*a1*d*d1*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - d*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - 
        RP*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) + 2*b1*c1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) + 
        2*a1*d1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - 10*a1*c1*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - 
        5*RP*pow(b1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - 2*b1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 2*a1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 
        2*b1*c1*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 16*a1*c1*d*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 2*a1*d1*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 
        8*d*RP*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 2*a1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 12*a1*b1*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 
        pow(b1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 2*b1*c1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) - 2*a1*d*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) - 
        2*b1*c1*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) + 16*a1*c1*d*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) - 
        2*a1*d1*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) + 8*d*RP*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) + 
        2*a1*c1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) - 12*a1*b1*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) + 
        pow(b1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) - 2*a1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 2*a1*c1*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 
        18*a1*b1*d*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - d*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - RP*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 
        2*a1*b1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 7*RP*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 
        2*a1*c1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) - 2*a1*c1*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) + 
        18*a1*b1*d*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) - d*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) - 
        RP*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) + 2*a1*b1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) - 
        7*RP*pow(a1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) - 2*a1*b1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,6) - 2*a1*b1*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,6) - 
        10*d*RP*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,6) - pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,6) - 2*a1*b1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6) - 
        2*a1*b1*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6) + 10*d*RP*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6) + 
        pow(a1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6) - d*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,7) - RP*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,7) - 
        d*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,7) - RP*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,7) - 16*c1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 
        10*c1*d1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 30*b1*d*d1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 15*d*k2*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 
        6*b1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 16*b1*c1*k2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 
        16*a1*d1*k2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 3*pow(c1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 4*pow(d1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
        16*c1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 10*c1*d1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 30*b1*d*d1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 
        15*d*k1*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 6*b1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 16*b1*c1*k1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 
        16*a1*d1*k1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 3*pow(c1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 4*pow(d1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 
        3*gsl_sf_expint_Ei(-(k2*pow(RP,-1)))*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*pow(RP,2) - 
        3*gsl_sf_expint_Ei(k2*pow(RP,-1))*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*pow(RP,2) + 16*c1*d*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) - 
        10*c1*d1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) + 30*b1*d*d1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) + 
        15*d*k1*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) - 6*b1*d1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) - 
        16*b1*c1*k1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) - 16*a1*d1*k1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) - 
        3*pow(c1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) - 4*pow(d1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) - 
        16*c1*d*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) + 10*c1*d1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) - 
        30*b1*d*d1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) - 15*d*k2*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) + 
        6*b1*d1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) + 16*b1*c1*k2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) + 
        16*a1*d1*k2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) + 3*pow(c1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) + 
        4*pow(d1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) + 48*b1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
        12*b1*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 48*a1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
        6*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 30*a1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
        15*pow(b1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 48*b1*c1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
        12*b1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 48*a1*d*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
        6*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 30*a1*c1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
        15*pow(b1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 14*b1*c1*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
        70*a1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 14*a1*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
        35*d*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 48*a1*b1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
        14*b1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 70*a1*c1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
        14*a1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 35*d*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
        48*a1*b1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 16*a1*c1*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 
        96*a1*b1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 8*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 
        35*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 16*a1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 
        96*a1*b1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 8*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
        35*pow(a1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 18*a1*b1*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 
        63*d*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 18*a1*b1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 
        63*d*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 10*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,6)*pow(RP,2) - 
        10*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6)*pow(RP,2) - 48*b1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
        12*b1*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 48*a1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
        6*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 30*a1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
        15*pow(b1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 48*b1*c1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
        12*b1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 48*a1*d*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
        6*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 30*a1*c1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
        15*pow(b1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 14*b1*c1*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
        70*a1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 14*a1*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
        35*d*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 48*a1*b1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
        14*b1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 70*a1*c1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
        14*a1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 35*d*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
        48*a1*b1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 16*a1*c1*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 
        96*a1*b1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 8*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 
        35*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 16*a1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 
        96*a1*b1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 8*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
        35*pow(a1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 18*a1*b1*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 
        63*d*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 18*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 
        63*d*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 10*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,6)*pow(RP,2) + 
        10*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6)*pow(RP,2) - 
        gsl_sf_expint_Ei(-(k1*pow(RP,-1)))*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*(3*d*RP + pow(d,2) + 3*pow(RP,2)) + 
        gsl_sf_expint_Ei(k1*pow(RP,-1))*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*(3*d*RP + pow(d,2) + 3*pow(RP,2)) - 16*c1*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 
        30*b1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 96*b1*c1*d*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 30*b1*d1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 
        96*a1*d*d1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 15*d*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 15*k2*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 
        16*b1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 16*a1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 60*a1*c1*k2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 
        30*k2*pow(b1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 16*c1*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 30*b1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 
        96*b1*c1*d*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 30*b1*d1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 96*a1*d*d1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 
        15*d*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 15*k1*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 16*b1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 
        16*a1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 60*a1*c1*k1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 
        30*k1*pow(b1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 16*c1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 
        30*b1*d*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 96*b1*c1*d*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) + 
        30*b1*d1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 96*a1*d*d1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 
        15*d*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) + 15*k1*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) + 
        16*b1*c1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) + 16*a1*d1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) + 
        60*a1*c1*k1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) + 30*k1*pow(b1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 
        16*c1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 30*b1*d*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 
        96*b1*c1*d*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) - 30*b1*d1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 
        96*a1*d*d1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 15*d*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) - 
        15*k2*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) - 16*b1*c1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) - 
        16*a1*d1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) - 60*a1*c1*k2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) - 
        30*k2*pow(b1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 48*b1*c1*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
        210*a1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 48*a1*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
        105*d*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 144*a1*b1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
        48*b1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 210*a1*c1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
        48*a1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 105*d*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
        144*a1*b1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 70*a1*c1*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
        384*a1*b1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 35*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
        140*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 70*a1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
        384*a1*b1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 35*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
        140*pow(a1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 96*a1*b1*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 
        315*d*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 96*a1*b1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 
        315*d*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 63*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,3) + 
        63*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,3) - 48*b1*c1*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
        210*a1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 48*a1*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
        105*d*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 144*a1*b1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
        48*b1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 210*a1*c1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
        48*a1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 105*d*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
        144*a1*b1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 70*a1*c1*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 
        384*a1*b1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 35*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 
        140*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 70*a1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
        384*a1*b1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 35*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 
        140*pow(a1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 96*a1*b1*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 
        315*d*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 96*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 
        315*d*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 63*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 
        63*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 96*b1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 30*b1*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 
        96*a1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 96*b1*c1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 420*a1*c1*d*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 
        96*a1*d1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 210*d*k2*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 15*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 
        60*a1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 288*a1*b1*k2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 
        30*pow(b1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 96*b1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 30*b1*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
        96*a1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 96*b1*c1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 420*a1*c1*d*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
        96*a1*d1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 210*d*k1*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 15*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
        60*a1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 288*a1*b1*k1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
        30*pow(b1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 96*b1*c1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 
        30*b1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 96*a1*d*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 
        96*b1*c1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 420*a1*c1*d*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 
        96*a1*d1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 210*d*k1*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 
        15*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 60*a1*c1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 
        288*a1*b1*k1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 30*pow(b1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 
        96*b1*c1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 30*b1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 96*a1*d*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 
        96*b1*c1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 420*a1*c1*d*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 
        96*a1*d1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 210*d*k2*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 
        15*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 60*a1*c1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 
        288*a1*b1*k2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 30*pow(b1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 
        210*a1*c1*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 1152*a1*b1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
        105*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 420*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
        210*a1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 1152*a1*b1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
        105*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 420*pow(a1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
        384*a1*b1*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 1260*d*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 
        384*a1*b1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 1260*d*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 
        315*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,4) - 315*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,4) - 
        210*a1*c1*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 1152*a1*b1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
        105*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 420*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 
        210*a1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 1152*a1*b1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 
        105*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 420*pow(a1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
        384*a1*b1*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 1260*d*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 
        384*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 1260*d*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 
        315*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,4) + 315*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,4) - 
        96*b1*c1*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 420*a1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 96*a1*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 
        420*a1*c1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 2304*a1*b1*d*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 210*d*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 
        210*k2*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 288*a1*b1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 
        840*k2*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 96*b1*c1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 420*a1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
        96*a1*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 420*a1*c1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 2304*a1*b1*d*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
        210*d*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 210*k1*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 288*a1*b1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
        840*k1*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 96*b1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 
        420*a1*c1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) + 96*a1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) + 
        420*a1*c1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 2304*a1*b1*d*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 
        210*d*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) + 210*k1*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) + 
        288*a1*b1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) + 840*k1*pow(a1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 
        96*b1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 420*a1*c1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) - 96*a1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) - 
        420*a1*c1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 2304*a1*b1*d*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 
        210*d*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) - 210*k2*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) - 
        288*a1*b1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) - 840*k2*pow(a1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 
        1152*a1*b1*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 3780*d*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 
        1152*a1*b1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 3780*d*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 
        1260*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,5) + 1260*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,5) - 
        1152*a1*b1*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 3780*d*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 
        1152*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 3780*d*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 
        1260*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 1260*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 
        420*a1*c1*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 2304*a1*b1*d*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 2304*a1*b1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 
        7560*d*k2*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 210*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 840*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 
        420*a1*c1*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 2304*a1*b1*d*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 2304*a1*b1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 
        7560*d*k1*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 210*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 840*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 
        420*a1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) + 2304*a1*b1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 
        2304*a1*b1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) + 7560*d*k1*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 
        210*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 840*pow(a1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) + 
        420*a1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) - 2304*a1*b1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 
        2304*a1*b1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) - 7560*d*k2*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 
        210*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 840*pow(a1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 
        3780*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,6) - 3780*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,6) - 
        3780*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,6) + 3780*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,6) - 
        2304*a1*b1*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 7560*d*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 7560*k2*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) + 
        2304*a1*b1*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 7560*d*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 7560*k1*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 
        2304*a1*b1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) - 7560*d*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) + 
        7560*k1*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) - 2304*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) + 
        7560*d*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) - 7560*k2*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) - 
        7560*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,8) + 7560*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,8) - 7560*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,8) + 
        7560*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,8))*pow(2*pow(d,3),-1);
    return res;
}
/* *********** CASE 1 DONE *********** */

/*      Case_2:     k1 < d < k2         */
/*      Case_2_Sub_1:   d < r < k2      */
double SH_Integral_2233_case_2_sub_1( double AP, double RP, double k1, double k2, double a1, double b1, double c1, double d1, double d )
{
    double res=3*AP*pow(M_E,-((2*d + k2)*pow(RP,-1)))*pow(RP,3)*(-10*c1*d*d1*k2*RP*pow(M_E,d*pow(RP,-1)) - 2*c1*d1*k2*pow(d,2)*pow(M_E,d*pow(RP,-1)) - 
        4*c1*d1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1)) - 6*b1*d1*k2*RP*pow(d,2)*pow(M_E,d*pow(RP,-1)) - 3*k2*RP*pow(c1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1)) - 
        d*k2*pow(d1,2)*pow(M_E,d*pow(RP,-1)) - 4*d*RP*pow(d1,2)*pow(M_E,d*pow(RP,-1)) - k2*RP*pow(d1,2)*pow(M_E,d*pow(RP,-1)) - pow(d,2)*pow(d1,2)*pow(M_E,d*pow(RP,-1)) - 
        10*c1*d*d1*k2*RP*pow(M_E,3*d*pow(RP,-1)) + 2*c1*d1*k2*pow(d,2)*pow(M_E,3*d*pow(RP,-1)) + 4*c1*d1*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1)) + 
        6*b1*d1*k2*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1)) + 3*k2*RP*pow(c1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1)) - d*k2*pow(d1,2)*pow(M_E,3*d*pow(RP,-1)) - 
        4*d*RP*pow(d1,2)*pow(M_E,3*d*pow(RP,-1)) + k2*RP*pow(d1,2)*pow(M_E,3*d*pow(RP,-1)) + pow(d,2)*pow(d1,2)*pow(M_E,3*d*pow(RP,-1)) + 
        16*c1*d1*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1)) + 4*c1*d1*pow(d,3)*pow(M_E,k2*pow(RP,-1)) + 20*b1*d1*RP*pow(d,3)*pow(M_E,k2*pow(RP,-1)) + 
        10*RP*pow(c1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1)) + 4*b1*d1*pow(d,4)*pow(M_E,k2*pow(RP,-1)) + 24*b1*c1*RP*pow(d,4)*pow(M_E,k2*pow(RP,-1)) + 
        24*a1*d1*RP*pow(d,4)*pow(M_E,k2*pow(RP,-1)) + 2*pow(c1,2)*pow(d,4)*pow(M_E,k2*pow(RP,-1)) + 4*b1*c1*pow(d,5)*pow(M_E,k2*pow(RP,-1)) + 
        4*a1*d1*pow(d,5)*pow(M_E,k2*pow(RP,-1)) + 28*a1*c1*RP*pow(d,5)*pow(M_E,k2*pow(RP,-1)) + 14*RP*pow(b1,2)*pow(d,5)*pow(M_E,k2*pow(RP,-1)) + 
        4*a1*c1*pow(d,6)*pow(M_E,k2*pow(RP,-1)) + 32*a1*b1*RP*pow(d,6)*pow(M_E,k2*pow(RP,-1)) + 2*pow(b1,2)*pow(d,6)*pow(M_E,k2*pow(RP,-1)) + 
        4*a1*b1*pow(d,7)*pow(M_E,k2*pow(RP,-1)) + 18*RP*pow(a1,2)*pow(d,7)*pow(M_E,k2*pow(RP,-1)) + 2*pow(a1,2)*pow(d,8)*pow(M_E,k2*pow(RP,-1)) + 
        5*d*RP*pow(d1,2)*pow(M_E,k2*pow(RP,-1)) + 2*pow(d,2)*pow(d1,2)*pow(M_E,k2*pow(RP,-1)) + 4*c1*d1*RP*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 
        4*b1*d1*RP*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 2*RP*pow(c1,2)*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 4*b1*c1*RP*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 
        4*a1*d1*RP*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 4*a1*c1*RP*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 2*RP*pow(b1,2)*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 
        4*a1*b1*RP*pow(d,6)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 2*RP*pow(a1,2)*pow(d,7)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 3*d*RP*pow(d1,2)*pow(M_E,(2*d + k2)*pow(RP,-1)) - 
        2*c1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k2,2) - 2*c1*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k2,2) - 12*b1*d*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k2,2) - 
        6*d*RP*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2) - 2*b1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2) - 8*b1*c1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2) - 
        8*a1*d1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2) - pow(c1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2) - 2*c1*d*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) + 
        2*c1*d1*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) - 12*b1*d*d1*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) - 6*d*RP*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) + 
        2*b1*d1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) + 8*b1*c1*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) + 
        8*a1*d1*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) + pow(c1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) - 2*b1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k2,3) - 
        14*b1*c1*d*RP*pow(M_E,d*pow(RP,-1))*pow(k2,3) - 2*b1*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k2,3) - 14*a1*d*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k2,3) - 
        d*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3) - RP*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3) - 2*b1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3) - 
        2*a1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3) - 10*a1*c1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3) - 
        5*RP*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3) - 2*b1*d*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) - 14*b1*c1*d*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) + 
        2*b1*d1*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) - 14*a1*d*d1*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) - d*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) + 
        RP*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) + 2*b1*c1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) + 2*a1*d1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) + 
        10*a1*c1*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) + 5*RP*pow(b1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) - 2*b1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k2,4) - 
        2*a1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k2,4) - 2*b1*c1*RP*pow(M_E,d*pow(RP,-1))*pow(k2,4) - 16*a1*c1*d*RP*pow(M_E,d*pow(RP,-1))*pow(k2,4) - 
        2*a1*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k2,4) - 8*d*RP*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4) - 2*a1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4) - 
        12*a1*b1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4) - pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4) - 2*b1*c1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) - 
        2*a1*d*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) + 2*b1*c1*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) - 16*a1*c1*d*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) + 
        2*a1*d1*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) - 8*d*RP*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) + 2*a1*c1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) + 
        12*a1*b1*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) + pow(b1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) - 2*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k2,5) - 
        2*a1*c1*RP*pow(M_E,d*pow(RP,-1))*pow(k2,5) - 18*a1*b1*d*RP*pow(M_E,d*pow(RP,-1))*pow(k2,5) - d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5) - 
        RP*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5) - 2*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5) - 7*RP*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5) - 
        2*a1*c1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) + 2*a1*c1*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) - 18*a1*b1*d*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) - 
        d*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) + RP*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) + 2*a1*b1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) + 
        7*RP*pow(a1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) - 2*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k2,6) - 2*a1*b1*RP*pow(M_E,d*pow(RP,-1))*pow(k2,6) - 
        10*d*RP*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,6) - pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,6) - 2*a1*b1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,6) + 
        2*a1*b1*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,6) - 10*d*RP*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6) + pow(a1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6) - 
        d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,7) - RP*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,7) - d*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,7) + 
        RP*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,7) - 16*c1*d*d1*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 10*c1*d1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 
        30*b1*d*d1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 15*d*k2*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 6*b1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 
        16*b1*c1*k2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 16*a1*d1*k2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 
        3*pow(c1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 4*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 16*c1*d*d1*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) + 
        10*c1*d1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) - 30*b1*d*d1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) - 15*d*k2*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) + 
        6*b1*d1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) + 16*b1*c1*k2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) + 
        16*a1*d1*k2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) + 3*pow(c1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) + 4*pow(d1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) + 
        26*c1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 48*b1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 24*pow(c1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 
        78*b1*c1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 78*a1*d1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 116*a1*c1*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 
        58*pow(b1,2)*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 162*a1*b1*pow(d,5)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 
        108*pow(a1,2)*pow(d,6)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 4*pow(d1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 6*c1*d*d1*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 
        12*b1*d1*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 6*pow(c1,2)*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 
        18*b1*c1*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 18*a1*d1*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 
        24*a1*c1*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 12*pow(b1,2)*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 
        30*a1*b1*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 18*pow(a1,2)*pow(d,6)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) - 
        4*pow(d1,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) - 48*b1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 12*b1*d1*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
        48*a1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 6*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
        30*a1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 15*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
        48*b1*c1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 12*b1*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
        48*a1*d*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 6*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
        30*a1*c1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 15*pow(b1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
        14*b1*c1*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 70*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 14*a1*d1*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
        35*d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 48*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
        14*b1*c1*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 70*a1*c1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
        14*a1*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 35*d*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
        48*a1*b1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 16*a1*c1*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 
        96*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 8*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 
        35*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 16*a1*c1*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 
        96*a1*b1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 8*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
        35*pow(a1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 18*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 
        63*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 18*a1*b1*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 
        63*d*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 10*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,6)*pow(RP,2) + 
        10*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6)*pow(RP,2) + gsl_sf_expint_Ei(-(d*pow(RP,-1)))*pow(d1,2)*pow(M_E,(d + k2)*pow(RP,-1))*
        (pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1))) - 3*d*RP*(1 + pow(M_E,2*d*pow(RP,-1))) + 3*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,2)) - 
        gsl_sf_expint_Ei(-(k2*pow(RP,-1)))*pow(d1,2)*pow(M_E,(d + k2)*pow(RP,-1))*
        (pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1))) - 3*d*RP*(1 + pow(M_E,2*d*pow(RP,-1))) + 3*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,2)) - 
        16*c1*d1*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 30*b1*d*d1*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 96*b1*c1*d*k2*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 
        30*b1*d1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 96*a1*d*d1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 15*d*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 
        15*k2*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 16*b1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 16*a1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 
        60*a1*c1*k2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 30*k2*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 16*c1*d1*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 
        30*b1*d*d1*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 96*b1*c1*d*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) + 30*b1*d1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 
        96*a1*d*d1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 15*d*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) + 15*k2*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) + 
        16*b1*c1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) + 16*a1*d1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) + 
        60*a1*c1*k2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) + 30*k2*pow(b1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) + 16*c1*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 
        60*b1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 30*d*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 160*b1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 
        160*a1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 340*a1*c1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 
        170*pow(b1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 624*a1*b1*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 
        518*pow(a1,2)*pow(d,5)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 16*c1*d1*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) + 
        32*b1*c1*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) + 32*a1*d1*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) + 
        80*a1*c1*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) + 40*pow(b1,2)*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) + 
        144*a1*b1*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) + 112*pow(a1,2)*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) - 
        48*b1*c1*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 210*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 48*a1*d1*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
        105*d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 144*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
        48*b1*c1*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 210*a1*c1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
        48*a1*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 105*d*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
        144*a1*b1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 70*a1*c1*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 
        384*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 35*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 
        140*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 70*a1*c1*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 
        384*a1*b1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 35*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
        140*pow(a1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 96*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 
        315*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 96*a1*b1*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 
        315*d*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 63*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,3) + 
        63*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 96*b1*c1*d*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 30*b1*d1*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
        96*a1*d*d1*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 96*b1*c1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 420*a1*c1*d*k2*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
        96*a1*d1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 210*d*k2*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 15*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
        60*a1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 288*a1*b1*k2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
        30*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 96*b1*c1*d*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 30*b1*d1*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 
        96*a1*d*d1*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 96*b1*c1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 420*a1*c1*d*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 
        96*a1*d1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 210*d*k2*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 15*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 
        60*a1*c1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 288*a1*b1*k2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 
        30*pow(b1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 192*b1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 30*b1*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
        192*a1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 15*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 690*a1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
        345*pow(b1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 1824*a1*b1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
        1995*pow(a1,2)*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 30*b1*d1*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) - 
        15*pow(c1,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) + 150*a1*c1*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) + 
        75*pow(b1,2)*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) + 480*a1*b1*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) + 
        525*pow(a1,2)*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) - 210*a1*c1*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
        1152*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 105*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
        420*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 210*a1*c1*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
        1152*a1*b1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 105*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 
        420*pow(a1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 384*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 
        1260*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 384*a1*b1*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 
        1260*d*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 315*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,4) + 
        315*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,4) - 96*b1*c1*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 420*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 
        96*a1*d1*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 420*a1*c1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 2304*a1*b1*d*k2*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 
        210*d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 210*k2*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 288*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 
        840*k2*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 96*b1*c1*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 420*a1*c1*d*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 
        96*a1*d1*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 420*a1*c1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 2304*a1*b1*d*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 
        210*d*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 210*k2*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 288*a1*b1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 
        840*k2*pow(a1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 96*b1*c1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 840*a1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
        96*a1*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 420*d*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 3744*a1*b1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
        5880*pow(a1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 96*b1*c1*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) - 96*a1*d1*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) + 
        864*a1*b1*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) + 1680*pow(a1,2)*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) - 
        1152*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 3780*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 
        1152*a1*b1*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 3780*d*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 
        1260*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,5) + 1260*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 
        420*a1*c1*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 2304*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 2304*a1*b1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 
        7560*d*k2*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 210*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 840*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 
        420*a1*c1*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 2304*a1*b1*d*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) + 2304*a1*b1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 
        7560*d*k2*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) + 210*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) + 
        840*pow(a1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) + 420*a1*c1*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 4608*a1*b1*d*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 
        210*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 12180*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 420*a1*c1*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,6) - 
        210*pow(b1,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,6) + 2940*pow(a1,2)*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,6) - 
        3780*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,6) + 3780*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,6) - 
        2304*a1*b1*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 7560*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 7560*k2*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 
        2304*a1*b1*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) - 7560*d*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) + 7560*k2*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) + 
        2304*a1*b1*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 15120*d*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) - 2304*a1*b1*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,7) - 
        7560*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,8) + 7560*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,8) + 7560*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,8) - 
        7560*pow(a1,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,8))*pow(2*pow(d,3),-1);
    return res;
}
/*      Case_2_Sub_2:   k1 < r < d      */
double SH_Integral_2233_case_2_sub_2( double AP, double RP, double k1, double k2, double a1, double b1, double c1, double d1, double d )
{
    double res=-3*AP*pow(M_E,-((2*d + k1)*pow(RP,-1)))*pow(RP,3)*(-10*c1*d*d1*k1*RP*pow(M_E,d*pow(RP,-1)) - 2*c1*d1*k1*pow(d,2)*pow(M_E,d*pow(RP,-1)) - 
        4*c1*d1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1)) - 6*b1*d1*k1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1)) - 3*k1*RP*pow(c1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1)) - 
        d*k1*pow(d1,2)*pow(M_E,d*pow(RP,-1)) - 4*d*RP*pow(d1,2)*pow(M_E,d*pow(RP,-1)) - k1*RP*pow(d1,2)*pow(M_E,d*pow(RP,-1)) - pow(d,2)*pow(d1,2)*pow(M_E,d*pow(RP,-1)) + 
        16*c1*d1*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1)) + 4*c1*d1*pow(d,3)*pow(M_E,k1*pow(RP,-1)) + 20*b1*d1*RP*pow(d,3)*pow(M_E,k1*pow(RP,-1)) + 
        10*RP*pow(c1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1)) + 4*b1*d1*pow(d,4)*pow(M_E,k1*pow(RP,-1)) + 24*b1*c1*RP*pow(d,4)*pow(M_E,k1*pow(RP,-1)) + 
        24*a1*d1*RP*pow(d,4)*pow(M_E,k1*pow(RP,-1)) + 2*pow(c1,2)*pow(d,4)*pow(M_E,k1*pow(RP,-1)) + 4*b1*c1*pow(d,5)*pow(M_E,k1*pow(RP,-1)) + 
        4*a1*d1*pow(d,5)*pow(M_E,k1*pow(RP,-1)) + 28*a1*c1*RP*pow(d,5)*pow(M_E,k1*pow(RP,-1)) + 14*RP*pow(b1,2)*pow(d,5)*pow(M_E,k1*pow(RP,-1)) + 
        4*a1*c1*pow(d,6)*pow(M_E,k1*pow(RP,-1)) + 32*a1*b1*RP*pow(d,6)*pow(M_E,k1*pow(RP,-1)) + 2*pow(b1,2)*pow(d,6)*pow(M_E,k1*pow(RP,-1)) + 
        4*a1*b1*pow(d,7)*pow(M_E,k1*pow(RP,-1)) + 18*RP*pow(a1,2)*pow(d,7)*pow(M_E,k1*pow(RP,-1)) + 2*pow(a1,2)*pow(d,8)*pow(M_E,k1*pow(RP,-1)) + 
        5*d*RP*pow(d1,2)*pow(M_E,k1*pow(RP,-1)) + 2*pow(d,2)*pow(d1,2)*pow(M_E,k1*pow(RP,-1)) + 
        3*d*RP*gsl_sf_expint_Ei(-(k1*pow(RP,-1)))*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1)) - 3*d*RP*gsl_sf_expint_Ei(k1*pow(RP,-1))*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1)) + 
        gsl_sf_expint_Ei(-(k1*pow(RP,-1)))*pow(d,2)*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1)) - gsl_sf_expint_Ei(k1*pow(RP,-1))*pow(d,2)*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1)) - 
        4*c1*d1*RP*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 4*b1*d1*RP*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 2*RP*pow(c1,2)*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 
        4*b1*c1*RP*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 4*a1*d1*RP*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 4*a1*c1*RP*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 
        2*RP*pow(b1,2)*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 4*a1*b1*RP*pow(d,6)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 
        2*RP*pow(a1,2)*pow(d,7)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 3*d*RP*pow(d1,2)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 10*c1*d*d1*k1*RP*pow(M_E,(d + 2*k1)*pow(RP,-1)) + 
        2*c1*d1*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1)) - 4*c1*d1*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1)) - 6*b1*d1*k1*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1)) - 
        3*k1*RP*pow(c1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1)) - d*k1*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1)) + 4*d*RP*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1)) - 
        k1*RP*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1)) + pow(d,2)*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1)) - 2*c1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 
        2*c1*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 12*b1*d*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 6*d*RP*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 
        2*b1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 8*b1*c1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 8*a1*d1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 
        pow(c1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 2*c1*d*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) - 2*c1*d1*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) + 
        12*b1*d*d1*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) + 6*d*RP*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) + 
        2*b1*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) - 8*b1*c1*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) - 
        8*a1*d1*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) + pow(c1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) - 
        2*b1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 14*b1*c1*d*RP*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 2*b1*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 
        14*a1*d*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k1,3) - d*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - RP*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 
        2*b1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 2*a1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 10*a1*c1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 
        5*RP*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 2*b1*d*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) + 
        14*b1*c1*d*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - 2*b1*d1*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) + 
        14*a1*d*d1*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - d*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - 
        RP*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) + 2*b1*c1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) + 
        2*a1*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - 10*a1*c1*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - 
        5*RP*pow(b1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - 2*b1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 2*a1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 
        2*b1*c1*RP*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 16*a1*c1*d*RP*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 2*a1*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 
        8*d*RP*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 2*a1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 12*a1*b1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 
        pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 2*b1*c1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) - 2*a1*d*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) - 
        2*b1*c1*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) + 16*a1*c1*d*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) - 
        2*a1*d1*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) + 8*d*RP*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) + 
        2*a1*c1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) - 12*a1*b1*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) + 
        pow(b1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) - 2*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 2*a1*c1*RP*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 
        18*a1*b1*d*RP*pow(M_E,d*pow(RP,-1))*pow(k1,5) - d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5) - RP*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 
        2*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 7*RP*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 2*a1*c1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) - 
        2*a1*c1*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) + 18*a1*b1*d*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) - 
        d*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) - RP*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) + 
        2*a1*b1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) - 7*RP*pow(a1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) - 
        2*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k1,6) - 2*a1*b1*RP*pow(M_E,d*pow(RP,-1))*pow(k1,6) - 10*d*RP*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,6) - 
        pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,6) - 2*a1*b1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6) - 2*a1*b1*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6) + 
        10*d*RP*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6) + pow(a1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6) - 
        d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,7) - RP*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,7) - d*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,7) - 
        RP*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,7) - 16*c1*d*d1*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 10*c1*d1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 
        30*b1*d*d1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 15*d*k1*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 6*b1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 
        16*b1*c1*k1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 16*a1*d1*k1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 
        3*pow(c1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 4*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 26*c1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
        48*b1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 24*pow(c1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 78*b1*c1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
        78*a1*d1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 116*a1*c1*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
        58*pow(b1,2)*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 162*a1*b1*pow(d,5)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
        108*pow(a1,2)*pow(d,6)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 4*pow(d1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
        3*gsl_sf_expint_Ei(-(k1*pow(RP,-1)))*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1))*pow(RP,2) - 
        3*gsl_sf_expint_Ei(k1*pow(RP,-1))*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1))*pow(RP,2) + 6*c1*d*d1*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 
        12*b1*d1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 6*pow(c1,2)*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 
        18*b1*c1*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 18*a1*d1*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 
        24*a1*c1*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 12*pow(b1,2)*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 
        30*a1*b1*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 18*pow(a1,2)*pow(d,6)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) - 
        4*pow(d1,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) - 16*c1*d*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) + 
        10*c1*d1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) - 30*b1*d*d1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) - 
        15*d*k1*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) + 6*b1*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) + 
        16*b1*c1*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) + 16*a1*d1*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) + 
        3*pow(c1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) + 4*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) - 
        48*b1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 12*b1*d1*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 48*a1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
        6*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 30*a1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
        15*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 48*b1*c1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
        12*b1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 48*a1*d*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
        6*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 30*a1*c1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
        15*pow(b1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 14*b1*c1*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
        70*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 14*a1*d1*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
        35*d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 48*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
        14*b1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 70*a1*c1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
        14*a1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 35*d*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
        48*a1*b1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 16*a1*c1*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
        96*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 8*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
        35*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 16*a1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
        96*a1*b1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 8*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 
        35*pow(a1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 18*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 
        63*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 18*a1*b1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 
        63*d*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 10*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,6)*pow(RP,2) + 
        10*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6)*pow(RP,2) - 
        gsl_sf_expint_Ei(-(d*pow(RP,-1)))*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1))*(3*d*RP + pow(d,2) + 3*pow(RP,2)) + 
        gsl_sf_expint_Ei(d*pow(RP,-1))*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1))*(3*d*RP + pow(d,2) + 3*pow(RP,2)) - 16*c1*d1*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 
        30*b1*d*d1*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 96*b1*c1*d*k1*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 30*b1*d1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 
        96*a1*d*d1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 15*d*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 15*k1*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 
        16*b1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 16*a1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 60*a1*c1*k1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 
        30*k1*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 16*c1*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 60*b1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 
        30*d*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 160*b1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 160*a1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 
        340*a1*c1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 170*pow(b1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 
        624*a1*b1*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 518*pow(a1,2)*pow(d,5)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 
        16*c1*d1*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 32*b1*c1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 
        32*a1*d1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 80*a1*c1*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 
        40*pow(b1,2)*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 144*a1*b1*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 
        112*pow(a1,2)*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 16*c1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 
        30*b1*d*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 96*b1*c1*d*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 
        30*b1*d1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 96*a1*d*d1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 
        15*d*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 15*k1*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 
        16*b1*c1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 16*a1*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 
        60*a1*c1*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 30*k1*pow(b1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 
        48*b1*c1*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 210*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 48*a1*d1*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
        105*d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 144*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
        48*b1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 210*a1*c1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
        48*a1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 105*d*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
        144*a1*b1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 70*a1*c1*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
        384*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 35*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
        140*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 70*a1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
        384*a1*b1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 35*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
        140*pow(a1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 96*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 
        315*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 96*a1*b1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 
        315*d*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 63*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,3) - 
        63*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,3) - 96*b1*c1*d*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 30*b1*d1*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
        96*a1*d*d1*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 96*b1*c1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 420*a1*c1*d*k1*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
        96*a1*d1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 210*d*k1*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 15*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
        60*a1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 288*a1*b1*k1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
        30*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 192*b1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 30*b1*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 
        192*a1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 15*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 690*a1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 
        345*pow(b1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 1824*a1*b1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 
        1995*pow(a1,2)*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 30*b1*d1*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) - 
        15*pow(c1,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) + 150*a1*c1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) + 
        75*pow(b1,2)*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) + 480*a1*b1*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) + 
        525*pow(a1,2)*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) - 96*b1*c1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 
        30*b1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 96*a1*d*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 96*b1*c1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 
        420*a1*c1*d*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 96*a1*d1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 
        210*d*k1*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 15*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 
        60*a1*c1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 288*a1*b1*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 
        30*pow(b1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 210*a1*c1*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
        1152*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 105*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
        420*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 210*a1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
        1152*a1*b1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 105*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
        420*pow(a1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 384*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 
        1260*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 384*a1*b1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 
        1260*d*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 315*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,4) + 
        315*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,4) - 96*b1*c1*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 420*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 
        96*a1*d1*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 420*a1*c1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 2304*a1*b1*d*k1*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 
        210*d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 210*k1*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 288*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 
        840*k1*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 96*b1*c1*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 840*a1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 
        96*a1*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 420*d*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 3744*a1*b1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 
        5880*pow(a1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 96*b1*c1*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) + 96*a1*d1*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) - 
        864*a1*b1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) - 1680*pow(a1,2)*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) - 
        96*b1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 420*a1*c1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 96*a1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 
        420*a1*c1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 2304*a1*b1*d*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 
        210*d*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 210*k1*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 
        288*a1*b1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 840*k1*pow(a1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 
        1152*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 3780*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 
        1152*a1*b1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 3780*d*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 
        1260*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,5) - 1260*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,5) - 
        420*a1*c1*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 2304*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 2304*a1*b1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 
        7560*d*k1*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 210*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 840*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 
        420*a1*c1*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 4608*a1*b1*d*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 210*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 
        12180*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 420*a1*c1*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,6) - 
        210*pow(b1,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,6) + 2940*pow(a1,2)*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,6) + 
        420*a1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 2304*a1*b1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) + 
        2304*a1*b1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 7560*d*k1*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) + 
        210*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) + 840*pow(a1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 
        3780*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,6) + 3780*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,6) - 
        2304*a1*b1*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 7560*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 7560*k1*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 
        2304*a1*b1*pow(M_E,k1*pow(RP,-1))*pow(RP,7) + 15120*d*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) + 2304*a1*b1*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,7) - 
        2304*a1*b1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) + 7560*d*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) - 
        7560*k1*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) - 7560*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,8) + 7560*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,8) - 
        7560*pow(a1,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,8) + 7560*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,8))*pow(2*pow(d,3),-1);
    return res;
}
/* *********** CASE 2 DONE *********** */

/*      Case_3:     d < k1 < k2         */
double SH_Integral_2233_case_3( double AP, double RP, double k1, double k2, double a1, double b1, double c1, double d1, double d )
{
    double res=3*AP*pow(M_E,-(d*pow(RP,-1)))*pow(RP,2)*(RP*(gsl_sf_expint_Ei(-(k1*pow(RP,-1))) - gsl_sf_expint_Ei(-(k2*pow(RP,-1))))*pow(d,2)*pow(d1,2)*(-1 + pow(M_E,2*d*pow(RP,-1))) + 
        RP*pow(d,2)*pow(d1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,k1*pow(RP,-1)) - pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1))) + 
        2*c1*d1*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1))) + 
        d*RP*pow(d1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1))) + 
        3*d*(-gsl_sf_expint_Ei(-(k1*pow(RP,-1))) + gsl_sf_expint_Ei(-(k2*pow(RP,-1))))*pow(d1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,2) - 
        2*c1*d1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,-(k1*pow(RP,-1))) - pow(M_E,-(k2*pow(RP,-1))))*pow(RP,2) + 
        3*d*pow(d1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,-(k1*pow(RP,-1))) - pow(M_E,-(k2*pow(RP,-1))))*pow(RP,2) + 
        2*b1*d1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2) + 
        pow(c1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2) + 
        pow(d1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2) + 
        6*c1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2) + 
        2*c1*d*d1*RP*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
        6*b1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
        3*d*pow(c1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
        2*b1*d1*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
        RP*pow(c1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
        2*c1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
        2*b1*c1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
        2*a1*d1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
        3*(gsl_sf_expint_Ei(-(k1*pow(RP,-1))) - gsl_sf_expint_Ei(-(k2*pow(RP,-1))))*pow(d1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,3) + 
        6*c1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,-(k1*pow(RP,-1))) - pow(M_E,-(k2*pow(RP,-1))))*pow(RP,3) + 
        3*pow(d1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,k1*pow(RP,-1)) - pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3) + 
        6*c1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3) + 
        6*b1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3) + 
        3*d*pow(c1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3) + 
        6*b1*c1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,3) + 
        6*a1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,3) + 
        6*b1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,3) + 
        3*pow(c1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,3) + 
        2*b1*d*d1*RP*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
        d*RP*pow(c1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
        6*b1*c1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
        6*a1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
        6*a1*c1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
        (pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
        3*d*pow(b1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
        (pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
        2*b1*c1*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
        2*a1*d1*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
        2*b1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
        pow(c1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
        2*a1*c1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
        pow(b1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
        6*b1*c1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
        (-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
        6*a1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
        (-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
        6*c1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,k1*pow(RP,-1)) - pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4) + 
        6*b1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4) + 
        3*pow(c1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4) + 
        6*b1*c1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,4) + 
        6*a1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,4) + 
        6*a1*c1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,4) + 
        3*pow(b1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,4) + 
        2*b1*c1*d*RP*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
        pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
        2*a1*d*d1*RP*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
        pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
        6*a1*c1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
        pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
        3*d*pow(b1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
        pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
        6*a1*b1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
        (pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
        pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
        2*a1*c1*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
        pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
        RP*pow(b1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
        pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
        2*b1*c1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
        pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
        2*a1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
        pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
        2*a1*b1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
        pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
        6*a1*c1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
        (-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
        pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
        3*pow(b1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
        (-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
        pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
        6*a1*b1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*
        (-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
        pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
        2*a1*c1*d*RP*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5)) - 
        pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
        d*RP*pow(b1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5)) - 
        pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
        6*a1*b1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5)) - 
        pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
        3*d*pow(a1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
        (pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5)) - 
        pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
        2*a1*b1*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))) + 
        pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
        2*a1*c1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))) + 
        pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
        pow(b1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))) + 
        pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
        pow(a1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))) + 
        pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
        6*a1*b1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
        (-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))) + 
        pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
        3*pow(a1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*
        (-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 120*pow(RP,5))) + 
        pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
        2*a1*b1*d*RP*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 
        720*pow(RP,6)) - pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 
        720*k2*pow(RP,5) + 720*pow(RP,6))) + 3*d*pow(a1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 
        720*pow(RP,6)) - pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 
        720*k2*pow(RP,5) + 720*pow(RP,6))) + RP*pow(a1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (-(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 
        720*pow(RP,6))) + pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 
        720*k2*pow(RP,5) + 720*pow(RP,6))) + 2*a1*b1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (-(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 
        720*pow(RP,6))) + pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 
        720*k2*pow(RP,5) + 720*pow(RP,6))) + 3*pow(a1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
        (-(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 
        720*pow(RP,6))) + pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 
        720*k2*pow(RP,5) + 720*pow(RP,6))) + d*RP*pow(a1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (pow(M_E,k2*pow(RP,-1))*(7*RP*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RP,2) + 210*pow(k1,4)*pow(RP,3) + 840*pow(k1,3)*pow(RP,4) + 2520*pow(k1,2)*pow(RP,5) + 
        5040*k1*pow(RP,6) + 5040*pow(RP,7)) - pow(M_E,k1*pow(RP,-1))*
        (7*RP*pow(k2,6) + pow(k2,7) + 42*pow(k2,5)*pow(RP,2) + 210*pow(k2,4)*pow(RP,3) + 840*pow(k2,3)*pow(RP,4) + 2520*pow(k2,2)*pow(RP,5) + 
        5040*k2*pow(RP,6) + 5040*pow(RP,7))) + pow(a1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (-(pow(M_E,k2*pow(RP,-1))*(7*RP*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RP,2) + 210*pow(k1,4)*pow(RP,3) + 840*pow(k1,3)*pow(RP,4) + 
        2520*pow(k1,2)*pow(RP,5) + 5040*k1*pow(RP,6) + 5040*pow(RP,7))) + 
        pow(M_E,k1*pow(RP,-1))*(7*RP*pow(k2,6) + pow(k2,7) + 42*pow(k2,5)*pow(RP,2) + 210*pow(k2,4)*pow(RP,3) + 840*pow(k2,3)*pow(RP,4) + 
        2520*pow(k2,2)*pow(RP,5) + 5040*k2*pow(RP,6) + 5040*pow(RP,7))))*pow(2*pow(d,3),-1);
    return res;
}
/* *********** CASE 3 DONE *********** */
/* ******************************************************************************************     3. < XorY | Kernel_SH | XorY >        */

/*      *       *       *       *       *       *       *       *       *       *       *       *       */

/* ******************************************************************************************     4. < Z | Kernel_SH | Z >        */

/*      Case_1:    k1 < k2 < d          */
double SH_Integral_44_case_1( double AP, double RP, double k1, double k2, double a1, double b1, double c1, double d1, double d )
{
    double res=-3*AP*pow(M_E,-((d + k1 + k2)*pow(RP,-1)))*pow(RP,2)*(-8*c1*d1*k2*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1)) - 
        2*d*k2*RP*pow(d1,2)*pow(M_E,k1*pow(RP,-1)) - k2*pow(d,2)*pow(d1,2)*pow(M_E,k1*pow(RP,-1)) - 
        3*RP*pow(d,2)*pow(d1,2)*pow(M_E,k1*pow(RP,-1)) + 8*c1*d1*k1*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1)) + 
        2*d*k1*RP*pow(d1,2)*pow(M_E,k2*pow(RP,-1)) + k1*pow(d,2)*pow(d1,2)*pow(M_E,k2*pow(RP,-1)) + 
        3*RP*pow(d,2)*pow(d1,2)*pow(M_E,k2*pow(RP,-1)) + 2*RP*gsl_sf_expint_Ei(-(k2*pow(RP,-1)))*pow(d,2)*pow(d1,2)*
        pow(M_E,(k1 + k2)*pow(RP,-1)) - 2*RP*gsl_sf_expint_Ei(k2*pow(RP,-1))*pow(d,2)*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1)) - 
        8*c1*d1*k1*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) + 2*d*k1*RP*pow(d1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) + 
        k1*pow(d,2)*pow(d1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) - 3*RP*pow(d,2)*pow(d1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1)) + 
        8*c1*d1*k2*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) - 2*d*k2*RP*pow(d1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) - 
        k2*pow(d,2)*pow(d1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) + 3*RP*pow(d,2)*pow(d1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1)) + 
        4*c1*d*d1*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 2*c1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 
        10*b1*d1*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 5*RP*pow(c1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2) + 
        4*c1*d*d1*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) + 2*c1*d1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) - 
        10*b1*d1*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) - 5*RP*pow(c1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2) + 
        4*b1*d*d1*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 2*d*RP*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 
        2*b1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 12*b1*c1*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 
        12*a1*d1*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + pow(c1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3) + 
        4*b1*d*d1*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + 2*d*RP*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + 
        2*b1*d1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) - 12*b1*c1*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) - 
        12*a1*d1*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + pow(c1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3) + 
        4*b1*c1*d*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 4*a1*d*d1*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 
        2*b1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 2*a1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 
        14*a1*c1*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 7*RP*pow(b1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4) + 
        4*b1*c1*d*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) + 4*a1*d*d1*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) + 
        2*b1*c1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) + 2*a1*d1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) - 
        14*a1*c1*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) - 7*RP*pow(b1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4) + 
        4*a1*c1*d*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 2*d*RP*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 
        2*a1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 16*a1*b1*RP*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 
        pow(b1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5) + 4*a1*c1*d*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) + 
        2*d*RP*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) + 2*a1*c1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) - 
        16*a1*b1*RP*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) + pow(b1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5) + 
        4*a1*b1*d*RP*pow(M_E,k2*pow(RP,-1))*pow(k1,6) + 2*a1*b1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,6) + 
        9*RP*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,6) + 4*a1*b1*d*RP*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6) + 
        2*a1*b1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6) - 9*RP*pow(a1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6) + 
        2*d*RP*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,7) + pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,7) + 
        2*d*RP*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,7) + pow(a1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,7) - 
        4*c1*d*d1*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 2*c1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 
        10*b1*d1*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 5*RP*pow(c1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2) - 
        4*c1*d*d1*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) - 2*c1*d1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) + 
        10*b1*d1*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) + 5*RP*pow(c1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2) - 
        4*b1*d*d1*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 2*d*RP*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 
        2*b1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 12*b1*c1*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 
        12*a1*d1*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - pow(c1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3) - 
        4*b1*d*d1*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - 2*d*RP*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - 
        2*b1*d1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) + 12*b1*c1*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) + 
        12*a1*d1*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - pow(c1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3) - 
        4*b1*c1*d*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 4*a1*d*d1*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 
        2*b1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 2*a1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 
        14*a1*c1*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 7*RP*pow(b1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4) - 
        4*b1*c1*d*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) - 4*a1*d*d1*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) - 
        2*b1*c1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) - 2*a1*d1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) + 
        14*a1*c1*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) + 7*RP*pow(b1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4) - 
        4*a1*c1*d*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 2*d*RP*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 
        2*a1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 16*a1*b1*RP*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 
        pow(b1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5) - 4*a1*c1*d*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) - 
        2*d*RP*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) - 2*a1*c1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) + 
        16*a1*b1*RP*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) - pow(b1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5) - 
        4*a1*b1*d*RP*pow(M_E,k1*pow(RP,-1))*pow(k2,6) - 2*a1*b1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,6) - 
        9*RP*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,6) - 4*a1*b1*d*RP*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6) - 
        2*a1*b1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6) + 9*RP*pow(a1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6) - 
        2*d*RP*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,7) - pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,7) - 
        2*d*RP*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,7) - pow(a1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,7) - 
        20*c1*d*d1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 12*c1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 
        24*b1*d1*k2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 12*k2*pow(c1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 
        8*d*pow(d1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) - 2*k2*pow(d1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
        20*c1*d*d1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 12*c1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 
        24*b1*d1*k1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 12*k1*pow(c1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 
        8*d*pow(d1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 2*k1*pow(d1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) + 
        6*d*gsl_sf_expint_Ei(-(k2*pow(RP,-1)))*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*pow(RP,2) - 
        6*d*gsl_sf_expint_Ei(k2*pow(RP,-1))*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*pow(RP,2) - 
        20*c1*d*d1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) + 12*c1*d1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) + 
        24*b1*d1*k1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) + 12*k1*pow(c1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) - 
        8*d*pow(d1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) + 2*k1*pow(d1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,2) + 
        20*c1*d*d1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) - 12*c1*d1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) - 
        24*b1*d1*k2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) - 12*k2*pow(c1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) + 
        8*d*pow(d1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) - 2*k2*pow(d1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,2) + 
        4*c1*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 24*b1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
        12*d*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 40*b1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
        40*a1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 4*c1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
        24*b1*d*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 12*d*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
        40*b1*c1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
        40*a1*d1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 28*b1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
        4*b1*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 28*a1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
        2*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 60*a1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
        30*pow(b1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 28*b1*c1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
        4*b1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 28*a1*d*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
        2*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
        60*a1*c1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 
        30*pow(b1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 4*b1*c1*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 
        32*a1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 4*a1*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 
        16*d*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 84*a1*b1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 
        4*b1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 32*a1*c1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 
        4*a1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 16*d*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 
        84*a1*b1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 4*a1*c1*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 
        36*a1*b1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 2*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 
        56*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 4*a1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 
        36*a1*b1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 2*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 
        56*pow(a1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 4*a1*b1*pow(M_E,k2*pow(RP,-1))*pow(k1,6)*pow(RP,2) + 
        20*d*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,6)*pow(RP,2) + 4*a1*b1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6)*pow(RP,2) - 
        20*d*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6)*pow(RP,2) + 2*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,7)*pow(RP,2) + 
        2*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,7)*pow(RP,2) - 4*c1*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
        24*b1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 12*d*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
        40*b1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 40*a1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
        4*c1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 24*b1*d*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
        12*d*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
        40*b1*c1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
        40*a1*d1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 28*b1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
        4*b1*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 28*a1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
        2*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 60*a1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
        30*pow(b1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 28*b1*c1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
        4*b1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 28*a1*d*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
        2*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
        60*a1*c1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
        30*pow(b1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 4*b1*c1*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 
        32*a1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 4*a1*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 
        16*d*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 84*a1*b1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 
        4*b1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 32*a1*c1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 
        4*a1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 16*d*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 
        84*a1*b1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 4*a1*c1*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 
        36*a1*b1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 2*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 
        56*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 4*a1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 
        36*a1*b1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 2*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 
        56*pow(a1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 4*a1*b1*pow(M_E,k1*pow(RP,-1))*pow(k2,6)*pow(RP,2) - 
        20*d*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,6)*pow(RP,2) - 4*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6)*pow(RP,2) + 
        20*d*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6)*pow(RP,2) - 2*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,7)*pow(RP,2) - 
        2*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,7)*pow(RP,2) - 
        2*RP*gsl_sf_expint_Ei(-(k1*pow(RP,-1)))*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*(3*d*RP + pow(d,2) + 3*pow(RP,2)) + 
        2*RP*gsl_sf_expint_Ei(k1*pow(RP,-1))*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*(3*d*RP + pow(d,2) + 3*pow(RP,2)) - 
        32*c1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 20*c1*d1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 
        60*b1*d*d1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 30*d*k2*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 
        24*b1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 80*b1*c1*k2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 
        80*a1*d1*k2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 12*pow(c1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) - 
        8*pow(d1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 32*c1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 
        20*c1*d1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 60*b1*d*d1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 
        30*d*k1*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 24*b1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 
        80*b1*c1*k1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 80*a1*d1*k1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 
        12*pow(c1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 8*pow(d1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) + 
        6*gsl_sf_expint_Ei(-(k2*pow(RP,-1)))*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*pow(RP,3) - 
        6*gsl_sf_expint_Ei(k2*pow(RP,-1))*pow(d1,2)*pow(M_E,(k1 + k2)*pow(RP,-1))*pow(RP,3) + 
        32*c1*d*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 20*c1*d1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) + 
        60*b1*d*d1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) + 30*d*k1*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 
        24*b1*d1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 80*b1*c1*k1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 
        80*a1*d1*k1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 12*pow(c1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 
        8*pow(d1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,3) - 32*c1*d*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 
        20*c1*d1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) - 60*b1*d*d1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) - 
        30*d*k2*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 24*b1*d1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 
        80*b1*c1*k2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 80*a1*d1*k2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 
        12*pow(c1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 8*pow(d1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,3) + 
        96*b1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 24*b1*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
        96*a1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 12*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
        180*a1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 90*pow(b1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
        96*b1*c1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 24*b1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
        96*a1*d*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 12*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
        180*a1*c1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
        90*pow(b1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 28*b1*c1*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
        140*a1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 28*a1*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
        70*d*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 336*a1*b1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
        28*b1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 140*a1*c1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
        28*a1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 70*d*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
        336*a1*b1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 32*a1*c1*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 
        192*a1*b1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 16*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 
        280*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 32*a1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 
        192*a1*b1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 16*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 
        280*pow(a1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 36*a1*b1*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,3) + 
        126*d*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,3) - 36*a1*b1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,3) + 
        126*d*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,3) + 20*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,6)*pow(RP,3) - 
        20*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,6)*pow(RP,3) - 96*b1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
        24*b1*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 96*a1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
        12*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 180*a1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
        90*pow(b1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 96*b1*c1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
        24*b1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 96*a1*d*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
        12*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
        180*a1*c1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
        90*pow(b1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 28*b1*c1*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 
        140*a1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 28*a1*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 
        70*d*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 336*a1*b1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
        28*b1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 140*a1*c1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
        28*a1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 70*d*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
        336*a1*b1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 32*a1*c1*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 
        192*a1*b1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 16*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 
        280*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 32*a1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 
        192*a1*b1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 16*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 
        280*pow(a1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 36*a1*b1*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 
        126*d*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,3) + 36*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 
        126*d*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 20*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,6)*pow(RP,3) + 
        20*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,6)*pow(RP,3) - 32*c1*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 
        60*b1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 192*b1*c1*d*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 
        60*b1*d1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 192*a1*d*d1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 
        30*d*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 30*k2*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 
        80*b1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 80*a1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 
        360*a1*c1*k2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) - 180*k2*pow(b1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 
        32*c1*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 60*b1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
        192*b1*c1*d*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 60*b1*d1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
        192*a1*d*d1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 30*d*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
        30*k1*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 80*b1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
        80*a1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 360*a1*c1*k1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
        180*k1*pow(b1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 32*c1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 
        60*b1*d*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 192*b1*c1*d*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 
        60*b1*d1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 192*a1*d*d1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 
        30*d*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 30*k1*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 
        80*b1*c1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 80*a1*d1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 
        360*a1*c1*k1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) + 
        180*k1*pow(b1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,4) - 32*c1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 
        60*b1*d*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 192*b1*c1*d*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 
        60*b1*d1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 192*a1*d*d1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 
        30*d*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 30*k2*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 
        80*b1*c1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 80*a1*d1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 
        360*a1*c1*k2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) - 
        180*k2*pow(b1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,4) + 96*b1*c1*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
        420*a1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 96*a1*d1*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
        210*d*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 1008*a1*b1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
        96*b1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 420*a1*c1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
        96*a1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 210*d*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 
        1008*a1*b1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 140*a1*c1*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 
        768*a1*b1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 70*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 
        1120*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 140*a1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 
        768*a1*b1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 70*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 
        1120*pow(a1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 192*a1*b1*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,4) + 
        630*d*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,4) + 192*a1*b1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,4) - 
        630*d*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,4) + 126*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,5)*pow(RP,4) + 
        126*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,5)*pow(RP,4) - 96*b1*c1*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
        420*a1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 96*a1*d1*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
        210*d*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 1008*a1*b1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
        96*b1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 420*a1*c1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
        96*a1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 210*d*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
        1008*a1*b1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 140*a1*c1*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 
        768*a1*b1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 70*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 
        1120*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 140*a1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 
        768*a1*b1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 70*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 
        1120*pow(a1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 192*a1*b1*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,4) - 
        630*d*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,4) - 192*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,4) + 
        630*d*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,4) - 126*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,5)*pow(RP,4) - 
        126*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,5)*pow(RP,4) - 192*b1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 
        60*b1*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 192*a1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 
        192*b1*c1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 840*a1*c1*d*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 
        192*a1*d1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 420*d*k2*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 
        30*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 360*a1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 
        2016*a1*b1*k2*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 180*pow(b1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 
        192*b1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 60*b1*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
        192*a1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 192*b1*c1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
        840*a1*c1*d*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 192*a1*d1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
        420*d*k1*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 30*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
        360*a1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 2016*a1*b1*k1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
        180*pow(b1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 192*b1*c1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 
        60*b1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) + 192*a1*d*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 
        192*b1*c1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) + 840*a1*c1*d*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 
        192*a1*d1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) + 420*d*k1*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 
        30*pow(c1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 360*a1*c1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 
        2016*a1*b1*k1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 180*pow(b1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,5) - 
        192*b1*c1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 60*b1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) - 
        192*a1*d*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 192*b1*c1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) - 
        840*a1*c1*d*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 192*a1*d1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) - 
        420*d*k2*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 30*pow(c1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 
        360*a1*c1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 2016*a1*b1*k2*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 
        180*pow(b1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,5) + 420*a1*c1*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 
        2304*a1*b1*d*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 210*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 
        3360*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 420*a1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 
        2304*a1*b1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 210*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 
        3360*pow(a1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 768*a1*b1*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,5) + 
        2520*d*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,5) - 768*a1*b1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,5) + 
        2520*d*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,5) + 630*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,4)*pow(RP,5) - 
        630*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,4)*pow(RP,5) - 420*a1*c1*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 
        2304*a1*b1*d*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 210*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 
        3360*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 420*a1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 
        2304*a1*b1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 210*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 
        3360*pow(a1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 768*a1*b1*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 
        2520*d*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,5) + 768*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 
        2520*d*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 630*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,4)*pow(RP,5) + 
        630*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,4)*pow(RP,5) - 192*b1*c1*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 
        840*a1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 192*a1*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 
        840*a1*c1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 4608*a1*b1*d*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 
        420*d*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 420*k2*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 
        2016*a1*b1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) - 6720*k2*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 
        192*b1*c1*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 840*a1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 192*a1*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 
        840*a1*c1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 4608*a1*b1*d*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 
        420*d*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 420*k1*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 
        2016*a1*b1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 6720*k1*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 
        192*b1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 840*a1*c1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) + 
        192*a1*d1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) + 840*a1*c1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 
        4608*a1*b1*d*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 420*d*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) + 
        420*k1*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) + 2016*a1*b1*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) + 
        6720*k1*pow(a1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,6) - 192*b1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 
        840*a1*c1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) - 192*a1*d1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) - 
        840*a1*c1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 4608*a1*b1*d*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 
        420*d*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) - 420*k2*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) - 
        2016*a1*b1*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) - 6720*k2*pow(a1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,6) + 
        2304*a1*b1*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,6) + 7560*d*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,6) + 
        2304*a1*b1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,6) - 
        7560*d*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,6) + 2520*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,3)*pow(RP,6) + 
        2520*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,3)*pow(RP,6) - 2304*a1*b1*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,6) - 
        7560*d*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,6) - 2304*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,6) + 
        7560*d*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,6) - 2520*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,3)*pow(RP,6) - 
        2520*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,3)*pow(RP,6) - 840*a1*c1*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 
        4608*a1*b1*d*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 4608*a1*b1*k2*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 
        15120*d*k2*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 420*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 
        6720*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) + 840*a1*c1*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 
        4608*a1*b1*d*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 4608*a1*b1*k1*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 
        15120*d*k1*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 420*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 
        6720*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) - 840*a1*c1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) + 
        4608*a1*b1*d*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) - 4608*a1*b1*k1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) + 
        15120*d*k1*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) - 420*pow(b1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) - 
        6720*pow(a1,2)*pow(d,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,7) + 840*a1*c1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) - 
        4608*a1*b1*d*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) + 4608*a1*b1*k2*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) - 
        15120*d*k2*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) + 420*pow(b1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) + 
        6720*pow(a1,2)*pow(d,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,7) + 7560*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(k1,2)*pow(RP,7) - 
        7560*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(k1,2)*pow(RP,7) - 7560*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(k2,2)*pow(RP,7) + 
        7560*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(k2,2)*pow(RP,7) - 4608*a1*b1*pow(M_E,k1*pow(RP,-1))*pow(RP,8) - 
        15120*d*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,8) - 15120*k2*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,8) + 
        4608*a1*b1*pow(M_E,k2*pow(RP,-1))*pow(RP,8) + 15120*d*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,8) + 
        15120*k1*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,8) + 4608*a1*b1*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,8) - 
        15120*d*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,8) + 15120*k1*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,8) - 
        4608*a1*b1*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,8) + 15120*d*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,8) - 
        15120*k2*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,8) - 15120*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,9) + 
        15120*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,9) - 15120*pow(a1,2)*pow(M_E,(2*k1 + k2)*pow(RP,-1))*pow(RP,9) + 
        15120*pow(a1,2)*pow(M_E,(k1 + 2*k2)*pow(RP,-1))*pow(RP,9))*pow(2*pow(d,3),-1);
    return res;
}
/* *********** CASE 1 DONE *********** */

/*      Case_2:     k1 < d < k2         */
/*      Case_2_Sub_1:   d < r < k2      */
double SH_Integral_44_case_2_sub_1( double AP, double RP, double k1, double k2, double a1, double b1, double c1, double d1, double d )
{
    double res=3*AP*pow(M_E,-((2*d + k2)*pow(RP,-1)))*pow(RP,2)*(8*c1*d1*k2*RP*pow(d,2)*pow(M_E,d*pow(RP,-1)) + 2*d*k2*RP*pow(d1,2)*pow(M_E,d*pow(RP,-1)) + 
        k2*pow(d,2)*pow(d1,2)*pow(M_E,d*pow(RP,-1)) + 3*RP*pow(d,2)*pow(d1,2)*pow(M_E,d*pow(RP,-1)) - 
        8*c1*d1*k2*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1)) + 2*d*k2*RP*pow(d1,2)*pow(M_E,3*d*pow(RP,-1)) - 
        k2*pow(d,2)*pow(d1,2)*pow(M_E,3*d*pow(RP,-1)) - 3*RP*pow(d,2)*pow(d1,2)*pow(M_E,3*d*pow(RP,-1)) - 
        12*c1*d1*RP*pow(d,3)*pow(M_E,k2*pow(RP,-1)) - 2*c1*d1*pow(d,4)*pow(M_E,k2*pow(RP,-1)) - 14*b1*d1*RP*pow(d,4)*pow(M_E,k2*pow(RP,-1)) - 
        7*RP*pow(c1,2)*pow(d,4)*pow(M_E,k2*pow(RP,-1)) - 2*b1*d1*pow(d,5)*pow(M_E,k2*pow(RP,-1)) - 16*b1*c1*RP*pow(d,5)*pow(M_E,k2*pow(RP,-1)) - 
        16*a1*d1*RP*pow(d,5)*pow(M_E,k2*pow(RP,-1)) - pow(c1,2)*pow(d,5)*pow(M_E,k2*pow(RP,-1)) - 2*b1*c1*pow(d,6)*pow(M_E,k2*pow(RP,-1)) - 
        2*a1*d1*pow(d,6)*pow(M_E,k2*pow(RP,-1)) - 18*a1*c1*RP*pow(d,6)*pow(M_E,k2*pow(RP,-1)) - 9*RP*pow(b1,2)*pow(d,6)*pow(M_E,k2*pow(RP,-1)) - 
        2*a1*c1*pow(d,7)*pow(M_E,k2*pow(RP,-1)) - 20*a1*b1*RP*pow(d,7)*pow(M_E,k2*pow(RP,-1)) - pow(b1,2)*pow(d,7)*pow(M_E,k2*pow(RP,-1)) - 
        2*a1*b1*pow(d,8)*pow(M_E,k2*pow(RP,-1)) - 11*RP*pow(a1,2)*pow(d,8)*pow(M_E,k2*pow(RP,-1)) - pow(a1,2)*pow(d,9)*pow(M_E,k2*pow(RP,-1)) - 
        5*RP*pow(d,2)*pow(d1,2)*pow(M_E,k2*pow(RP,-1)) - pow(d,3)*pow(d1,2)*pow(M_E,k2*pow(RP,-1)) + 
        4*c1*d1*RP*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 2*c1*d1*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 
        6*b1*d1*RP*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 3*RP*pow(c1,2)*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 
        2*b1*d1*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 8*b1*c1*RP*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 
        8*a1*d1*RP*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1)) + pow(c1,2)*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 
        2*b1*c1*pow(d,6)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 2*a1*d1*pow(d,6)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 
        10*a1*c1*RP*pow(d,6)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 5*RP*pow(b1,2)*pow(d,6)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 
        2*a1*c1*pow(d,7)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 12*a1*b1*RP*pow(d,7)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 
        pow(b1,2)*pow(d,7)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 2*a1*b1*pow(d,8)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 
        7*RP*pow(a1,2)*pow(d,8)*pow(M_E,(2*d + k2)*pow(RP,-1)) + pow(a1,2)*pow(d,9)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 
        RP*pow(d,2)*pow(d1,2)*pow(M_E,(2*d + k2)*pow(RP,-1)) + pow(d,3)*pow(d1,2)*pow(M_E,(2*d + k2)*pow(RP,-1)) + 
        4*c1*d*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k2,2) + 2*c1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2) + 
        10*b1*d1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2) + 5*RP*pow(c1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2) + 
        4*c1*d*d1*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) - 2*c1*d1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) - 
        10*b1*d1*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) - 5*RP*pow(c1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2) + 
        4*b1*d*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 2*d*RP*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 
        2*b1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 12*b1*c1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 
        12*a1*d1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + pow(c1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3) + 
        4*b1*d*d1*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) + 2*d*RP*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) - 
        2*b1*d1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) - 12*b1*c1*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) - 
        12*a1*d1*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) - pow(c1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3) + 
        4*b1*c1*d*RP*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 4*a1*d*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 
        2*b1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 2*a1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 
        14*a1*c1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 7*RP*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4) + 
        4*b1*c1*d*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) + 4*a1*d*d1*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) - 
        2*b1*c1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) - 2*a1*d1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) - 
        14*a1*c1*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) - 7*RP*pow(b1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4) + 
        4*a1*c1*d*RP*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 2*d*RP*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 
        2*a1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 16*a1*b1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 
        pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5) + 4*a1*c1*d*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) + 
        2*d*RP*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) - 2*a1*c1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) - 
        16*a1*b1*RP*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) - pow(b1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5) + 
        4*a1*b1*d*RP*pow(M_E,d*pow(RP,-1))*pow(k2,6) + 2*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,6) + 
        9*RP*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,6) + 4*a1*b1*d*RP*pow(M_E,3*d*pow(RP,-1))*pow(k2,6) - 
        2*a1*b1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6) - 9*RP*pow(a1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6) + 
        2*d*RP*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,7) + pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,7) + 
        2*d*RP*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,7) - pow(a1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,7) + 
        20*c1*d*d1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 12*c1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 
        24*b1*d1*k2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 12*k2*pow(c1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 
        8*d*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 2*k2*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 
        20*c1*d*d1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) - 12*c1*d1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) - 
        24*b1*d1*k2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) - 12*k2*pow(c1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) + 
        8*d*pow(d1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) - 2*k2*pow(d1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,2) - 
        36*c1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 52*b1*d1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 
        26*pow(c1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 72*b1*c1*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 
        72*a1*d1*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 96*a1*c1*pow(d,5)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 
        48*pow(b1,2)*pow(d,5)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 124*a1*b1*pow(d,6)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 
        78*pow(a1,2)*pow(d,7)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 10*d*pow(d1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,2) - 
        4*c1*d1*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 4*b1*d1*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 
        2*pow(c1,2)*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 16*b1*c1*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 
        16*a1*d1*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 32*a1*c1*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 
        16*pow(b1,2)*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 52*a1*b1*pow(d,6)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 
        38*pow(a1,2)*pow(d,7)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) - 6*d*pow(d1,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,2) + 
        4*c1*d1*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 24*b1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
        12*d*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 40*b1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
        40*a1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 4*c1*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
        24*b1*d*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 12*d*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 
        40*b1*c1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) - 40*a1*d1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,2) + 
        28*b1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 4*b1*d1*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
        28*a1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 2*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
        60*a1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 30*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
        28*b1*c1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 4*b1*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
        28*a1*d*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 2*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 
        60*a1*c1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) - 30*pow(b1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,2) + 
        4*b1*c1*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 32*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
        4*a1*d1*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 16*d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
        84*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 4*b1*c1*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
        32*a1*c1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 4*a1*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
        16*d*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) - 84*a1*b1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,2) + 
        4*a1*c1*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 36*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 
        2*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 56*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 
        4*a1*c1*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 36*a1*b1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 
        2*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,2) - 56*pow(a1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,2) + 
        4*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k2,6)*pow(RP,2) + 20*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,6)*pow(RP,2) - 
        4*a1*b1*pow(M_E,3*d*pow(RP,-1))*pow(k2,6)*pow(RP,2) + 20*d*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6)*pow(RP,2) + 
        2*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,7)*pow(RP,2) - 2*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,7)*pow(RP,2) - 
        2*RP*gsl_sf_expint_Ei(-(d*pow(RP,-1)))*pow(d1,2)*pow(M_E,(d + k2)*pow(RP,-1))*
        (pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1))) - 3*d*RP*(1 + pow(M_E,2*d*pow(RP,-1))) + 3*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,2)) + 
        2*RP*gsl_sf_expint_Ei(-(k2*pow(RP,-1)))*pow(d1,2)*pow(M_E,(d + k2)*pow(RP,-1))*
        (pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1))) - 3*d*RP*(1 + pow(M_E,2*d*pow(RP,-1))) + 3*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,2)) + 
        32*c1*d*d1*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 20*c1*d1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 
        60*b1*d*d1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 30*d*k2*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 
        24*b1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 80*b1*c1*k2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 
        80*a1*d1*k2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 12*pow(c1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 
        8*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 32*c1*d*d1*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 
        20*c1*d1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) + 60*b1*d*d1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) + 
        30*d*k2*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 24*b1*d1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 
        80*b1*c1*k2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 80*a1*d1*k2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 
        12*pow(c1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 8*pow(d1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,3) - 
        52*c1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 108*b1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 
        54*pow(c1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 204*b1*c1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 
        204*a1*d1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 352*a1*c1*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 
        176*pow(b1,2)*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 564*a1*b1*pow(d,5)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 
        426*pow(a1,2)*pow(d,6)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 8*pow(d1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,3) - 
        12*c1*d*d1*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) - 12*b1*d1*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) - 
        6*pow(c1,2)*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) + 12*b1*c1*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) + 
        12*a1*d1*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) + 72*a1*c1*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) + 
        36*pow(b1,2)*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) + 180*a1*b1*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) + 
        174*pow(a1,2)*pow(d,6)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) + 8*pow(d1,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,3) + 
        96*b1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 24*b1*d1*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
        96*a1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 12*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
        180*a1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 90*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
        96*b1*c1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 24*b1*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
        96*a1*d*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 12*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 
        180*a1*c1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) - 90*pow(b1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,3) + 
        28*b1*c1*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 140*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
        28*a1*d1*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 70*d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
        336*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 28*b1*c1*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
        140*a1*c1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 28*a1*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
        70*d*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) - 336*a1*b1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,3) + 
        32*a1*c1*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 192*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 
        16*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 280*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 
        32*a1*c1*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 192*a1*b1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 
        16*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,3) - 280*pow(a1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,3) + 
        36*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,3) + 126*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,3) - 
        36*a1*b1*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,3) + 126*d*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,3) + 
        20*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,6)*pow(RP,3) - 20*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,6)*pow(RP,3) + 
        32*c1*d1*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 60*b1*d*d1*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 192*b1*c1*d*k2*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 
        60*b1*d1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 192*a1*d*d1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 
        30*d*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 30*k2*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 
        80*b1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 80*a1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 
        360*a1*c1*k2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 180*k2*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
        32*c1*d1*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 60*b1*d*d1*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 
        192*b1*c1*d*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 60*b1*d1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 
        192*a1*d*d1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) + 30*d*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 
        30*k2*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 80*b1*c1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 
        80*a1*d1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 360*a1*c1*k2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 
        180*k2*pow(b1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,4) - 32*c1*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 
        120*b1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 60*d*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 
        368*b1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 368*a1*d1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 
        920*a1*c1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 460*pow(b1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 
        1968*a1*b1*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) - 1876*pow(a1,2)*pow(d,5)*pow(M_E,k2*pow(RP,-1))*pow(RP,4) + 
        32*c1*d1*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) - 16*b1*c1*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) - 
        16*a1*d1*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) + 80*a1*c1*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) + 
        40*pow(b1,2)*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) + 432*a1*b1*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) + 
        616*pow(a1,2)*pow(d,5)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,4) + 96*b1*c1*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 
        420*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 96*a1*d1*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 
        210*d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 1008*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
        96*b1*c1*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 420*a1*c1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
        96*a1*d1*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 210*d*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) - 
        1008*a1*b1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,4) + 140*a1*c1*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 
        768*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 70*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 
        1120*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 140*a1*c1*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 
        768*a1*b1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 70*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,4) - 
        1120*pow(a1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,4) + 192*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,4) + 
        630*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,4) - 192*a1*b1*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,4) + 
        630*d*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,4) + 126*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,5)*pow(RP,4) - 
        126*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,5)*pow(RP,4) + 192*b1*c1*d*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 
        60*b1*d1*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 192*a1*d*d1*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 192*b1*c1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 
        840*a1*c1*d*k2*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 192*a1*d1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 
        420*d*k2*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 30*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 
        360*a1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 2016*a1*b1*k2*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 
        180*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 192*b1*c1*d*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 
        60*b1*d1*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 192*a1*d*d1*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 
        192*b1*c1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 840*a1*c1*d*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 
        192*a1*d1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) + 420*d*k2*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 
        30*pow(c1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 360*a1*c1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 
        2016*a1*b1*k2*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 180*pow(b1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,5) - 
        384*b1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 60*b1*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 
        384*a1*d*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 30*pow(c1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 
        1620*a1*c1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 810*pow(b1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 
        5088*a1*b1*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) - 6510*pow(a1,2)*pow(d,4)*pow(M_E,k2*pow(RP,-1))*pow(RP,5) + 
        60*b1*d1*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) + 30*pow(c1,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) - 
        60*a1*c1*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) - 30*pow(b1,2)*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) + 
        480*a1*b1*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) + 1470*pow(a1,2)*pow(d,4)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,5) + 
        420*a1*c1*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 2304*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 
        210*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 3360*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 
        420*a1*c1*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 2304*a1*b1*d*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 
        210*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,5) - 3360*pow(a1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,5) + 
        768*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,5) + 2520*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,5) - 
        768*a1*b1*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,5) + 2520*d*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,5) + 
        630*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,4)*pow(RP,5) - 630*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,4)*pow(RP,5) + 
        192*b1*c1*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 840*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 192*a1*d1*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 
        840*a1*c1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 4608*a1*b1*d*k2*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 
        420*d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 420*k2*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 
        2016*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 6720*k2*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 
        192*b1*c1*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) + 840*a1*c1*d*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 
        192*a1*d1*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 840*a1*c1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) + 
        4608*a1*b1*d*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) + 420*d*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 
        420*k2*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 2016*a1*b1*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 
        6720*k2*pow(a1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,6) - 192*b1*c1*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 
        1680*a1*c1*d*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 192*a1*d1*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 
        840*d*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 8928*a1*b1*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) - 
        16800*pow(a1,2)*pow(d,3)*pow(M_E,k2*pow(RP,-1))*pow(RP,6) + 192*b1*c1*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,6) + 
        192*a1*d1*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,6) - 288*a1*b1*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,6) + 
        1680*pow(a1,2)*pow(d,3)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,6) + 2304*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,6) + 
        7560*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,6) - 2304*a1*b1*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,6) + 
        7560*d*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,6) + 2520*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,3)*pow(RP,6) - 
        2520*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,3)*pow(RP,6) + 840*a1*c1*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 
        4608*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 4608*a1*b1*k2*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 
        15120*d*k2*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 420*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 
        6720*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 840*a1*c1*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) + 
        4608*a1*b1*d*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) - 4608*a1*b1*k2*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) + 
        15120*d*k2*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) - 420*pow(b1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) - 
        6720*pow(a1,2)*pow(d,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,7) - 840*a1*c1*pow(M_E,k2*pow(RP,-1))*pow(RP,7) - 
        9216*a1*b1*d*pow(M_E,k2*pow(RP,-1))*pow(RP,7) - 420*pow(b1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) - 
        29400*pow(a1,2)*pow(d,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,7) + 840*a1*c1*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,7) + 
        420*pow(b1,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,7) - 840*pow(a1,2)*pow(d,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,7) + 
        7560*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k2,2)*pow(RP,7) - 7560*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(k2,2)*pow(RP,7) + 
        4608*a1*b1*pow(M_E,d*pow(RP,-1))*pow(RP,8) + 15120*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,8) + 
        15120*k2*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,8) - 4608*a1*b1*pow(M_E,3*d*pow(RP,-1))*pow(RP,8) + 
        15120*d*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,8) - 15120*k2*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,8) - 
        4608*a1*b1*pow(M_E,k2*pow(RP,-1))*pow(RP,8) - 30240*d*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,8) + 
        4608*a1*b1*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,8) + 15120*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,9) - 
        15120*pow(a1,2)*pow(M_E,3*d*pow(RP,-1))*pow(RP,9) - 15120*pow(a1,2)*pow(M_E,k2*pow(RP,-1))*pow(RP,9) + 
        15120*pow(a1,2)*pow(M_E,(2*d + k2)*pow(RP,-1))*pow(RP,9))*pow(2*pow(d,3),-1);
    return res;
}
/*      Case_2_Sub_2:   k1 < r < d      */
double SH_Integral_44_case_2_sub_2( double AP, double RP, double k1, double k2, double a1, double b1, double c1, double d1, double d )
{
    double res=3*AP*pow(M_E,-((2*d + k1)*pow(RP,-1)))*pow(RP,2)*(-8*c1*d1*k1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1)) - 2*d*k1*RP*pow(d1,2)*pow(M_E,d*pow(RP,-1)) - 
        k1*pow(d,2)*pow(d1,2)*pow(M_E,d*pow(RP,-1)) - 3*RP*pow(d,2)*pow(d1,2)*pow(M_E,d*pow(RP,-1)) + 
        12*c1*d1*RP*pow(d,3)*pow(M_E,k1*pow(RP,-1)) + 2*c1*d1*pow(d,4)*pow(M_E,k1*pow(RP,-1)) + 14*b1*d1*RP*pow(d,4)*pow(M_E,k1*pow(RP,-1)) + 
        7*RP*pow(c1,2)*pow(d,4)*pow(M_E,k1*pow(RP,-1)) + 2*b1*d1*pow(d,5)*pow(M_E,k1*pow(RP,-1)) + 16*b1*c1*RP*pow(d,5)*pow(M_E,k1*pow(RP,-1)) + 
        16*a1*d1*RP*pow(d,5)*pow(M_E,k1*pow(RP,-1)) + pow(c1,2)*pow(d,5)*pow(M_E,k1*pow(RP,-1)) + 2*b1*c1*pow(d,6)*pow(M_E,k1*pow(RP,-1)) + 
        2*a1*d1*pow(d,6)*pow(M_E,k1*pow(RP,-1)) + 18*a1*c1*RP*pow(d,6)*pow(M_E,k1*pow(RP,-1)) + 9*RP*pow(b1,2)*pow(d,6)*pow(M_E,k1*pow(RP,-1)) + 
        2*a1*c1*pow(d,7)*pow(M_E,k1*pow(RP,-1)) + 20*a1*b1*RP*pow(d,7)*pow(M_E,k1*pow(RP,-1)) + pow(b1,2)*pow(d,7)*pow(M_E,k1*pow(RP,-1)) + 
        2*a1*b1*pow(d,8)*pow(M_E,k1*pow(RP,-1)) + 11*RP*pow(a1,2)*pow(d,8)*pow(M_E,k1*pow(RP,-1)) + pow(a1,2)*pow(d,9)*pow(M_E,k1*pow(RP,-1)) + 
        5*RP*pow(d,2)*pow(d1,2)*pow(M_E,k1*pow(RP,-1)) + pow(d,3)*pow(d1,2)*pow(M_E,k1*pow(RP,-1)) + 
        2*RP*gsl_sf_expint_Ei(-(k1*pow(RP,-1)))*pow(d,2)*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1)) - 
        2*RP*gsl_sf_expint_Ei(k1*pow(RP,-1))*pow(d,2)*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1)) - 
        4*c1*d1*RP*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 2*c1*d1*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 
        6*b1*d1*RP*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 3*RP*pow(c1,2)*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 
        2*b1*d1*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 8*b1*c1*RP*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 
        8*a1*d1*RP*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1)) + pow(c1,2)*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 
        2*b1*c1*pow(d,6)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 2*a1*d1*pow(d,6)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 
        10*a1*c1*RP*pow(d,6)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 5*RP*pow(b1,2)*pow(d,6)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 
        2*a1*c1*pow(d,7)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 12*a1*b1*RP*pow(d,7)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 
        pow(b1,2)*pow(d,7)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 2*a1*b1*pow(d,8)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 
        7*RP*pow(a1,2)*pow(d,8)*pow(M_E,(2*d + k1)*pow(RP,-1)) + pow(a1,2)*pow(d,9)*pow(M_E,(2*d + k1)*pow(RP,-1)) - 
        RP*pow(d,2)*pow(d1,2)*pow(M_E,(2*d + k1)*pow(RP,-1)) + pow(d,3)*pow(d1,2)*pow(M_E,(2*d + k1)*pow(RP,-1)) + 
        8*c1*d1*k1*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1)) - 2*d*k1*RP*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1)) - 
        k1*pow(d,2)*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1)) + 3*RP*pow(d,2)*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1)) - 
        4*c1*d*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 2*c1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 
        10*b1*d1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 5*RP*pow(c1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2) - 
        4*c1*d*d1*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) - 2*c1*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) + 
        10*b1*d1*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) + 5*RP*pow(c1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2) - 
        4*b1*d*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 2*d*RP*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 
        2*b1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 12*b1*c1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 
        12*a1*d1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - pow(c1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3) - 
        4*b1*d*d1*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - 2*d*RP*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - 
        2*b1*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) + 12*b1*c1*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) + 
        12*a1*d1*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - pow(c1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3) - 
        4*b1*c1*d*RP*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 4*a1*d*d1*RP*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 
        2*b1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 2*a1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 
        14*a1*c1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 7*RP*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4) - 
        4*b1*c1*d*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) - 4*a1*d*d1*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) - 
        2*b1*c1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) - 2*a1*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) + 
        14*a1*c1*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) + 7*RP*pow(b1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4) - 
        4*a1*c1*d*RP*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 2*d*RP*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 
        2*a1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 16*a1*b1*RP*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 
        pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5) - 4*a1*c1*d*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) - 
        2*d*RP*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) - 2*a1*c1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) + 
        16*a1*b1*RP*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) - pow(b1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5) - 
        4*a1*b1*d*RP*pow(M_E,d*pow(RP,-1))*pow(k1,6) - 2*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,6) - 
        9*RP*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,6) - 4*a1*b1*d*RP*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6) - 
        2*a1*b1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6) + 9*RP*pow(a1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6) - 
        2*d*RP*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,7) - pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,7) - 
        2*d*RP*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,7) - pow(a1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,7) - 
        20*c1*d*d1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 12*c1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 
        24*b1*d1*k1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 12*k1*pow(c1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 
        8*d*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) - 2*k1*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,2) + 
        36*c1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 52*b1*d1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
        26*pow(c1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 72*b1*c1*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
        72*a1*d1*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 96*a1*c1*pow(d,5)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
        48*pow(b1,2)*pow(d,5)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 124*a1*b1*pow(d,6)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
        78*pow(a1,2)*pow(d,7)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 10*d*pow(d1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,2) + 
        6*d*gsl_sf_expint_Ei(-(k1*pow(RP,-1)))*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1))*pow(RP,2) - 
        6*d*gsl_sf_expint_Ei(k1*pow(RP,-1))*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1))*pow(RP,2) - 
        4*c1*d1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 4*b1*d1*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 
        2*pow(c1,2)*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 16*b1*c1*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 
        16*a1*d1*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 32*a1*c1*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 
        16*pow(b1,2)*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 52*a1*b1*pow(d,6)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 
        38*pow(a1,2)*pow(d,7)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) - 6*d*pow(d1,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,2) + 
        20*c1*d*d1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) - 12*c1*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) - 
        24*b1*d1*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) - 12*k1*pow(c1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) + 
        8*d*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) - 2*k1*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,2) - 
        4*c1*d1*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 24*b1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
        12*d*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 40*b1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
        40*a1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 4*c1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 
        24*b1*d*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) + 12*d*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
        40*b1*c1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 
        40*a1*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,2) - 28*b1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
        4*b1*d1*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 28*a1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
        2*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 60*a1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
        30*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 28*b1*c1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
        4*b1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) + 28*a1*d*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
        2*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
        60*a1*c1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 
        30*pow(b1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,2) - 4*b1*c1*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
        32*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 4*a1*d1*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
        16*d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 84*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
        4*b1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 32*a1*c1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
        4*a1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) + 16*d*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 
        84*a1*b1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,2) - 4*a1*c1*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 
        36*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 2*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 
        56*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 4*a1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,2) + 
        36*a1*b1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 2*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 
        56*pow(a1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,2) - 4*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k1,6)*pow(RP,2) - 
        20*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,6)*pow(RP,2) - 4*a1*b1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6)*pow(RP,2) + 
        20*d*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6)*pow(RP,2) - 2*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,7)*pow(RP,2) - 
        2*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,7)*pow(RP,2) - 
        2*RP*gsl_sf_expint_Ei(-(d*pow(RP,-1)))*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1))*(3*d*RP + pow(d,2) + 3*pow(RP,2)) + 
        2*RP*gsl_sf_expint_Ei(d*pow(RP,-1))*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1))*(3*d*RP + pow(d,2) + 3*pow(RP,2)) - 
        32*c1*d*d1*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 20*c1*d1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 
        60*b1*d*d1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 30*d*k1*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 
        24*b1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 80*b1*c1*k1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 
        80*a1*d1*k1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 12*pow(c1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) - 
        8*pow(d1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,3) + 52*c1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 
        108*b1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 54*pow(c1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 
        204*b1*c1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 204*a1*d1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 
        352*a1*c1*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 176*pow(b1,2)*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 
        564*a1*b1*pow(d,5)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 426*pow(a1,2)*pow(d,6)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 
        8*pow(d1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,3) + 6*gsl_sf_expint_Ei(-(k1*pow(RP,-1)))*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1))*pow(RP,3) - 
        6*gsl_sf_expint_Ei(k1*pow(RP,-1))*pow(d1,2)*pow(M_E,(d + k1)*pow(RP,-1))*pow(RP,3) + 
        12*c1*d*d1*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) + 12*b1*d1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) + 
        6*pow(c1,2)*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 12*b1*c1*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 
        12*a1*d1*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 72*a1*c1*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 
        36*pow(b1,2)*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 180*a1*b1*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 
        174*pow(a1,2)*pow(d,6)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 8*pow(d1,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,3) - 
        32*c1*d*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 20*c1*d1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 
        60*b1*d*d1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 30*d*k1*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 
        24*b1*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 80*b1*c1*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 
        80*a1*d1*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 12*pow(c1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) + 
        8*pow(d1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,3) - 96*b1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
        24*b1*d1*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 96*a1*d*d1*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
        12*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 180*a1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 
        90*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 96*b1*c1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
        24*b1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 96*a1*d*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
        12*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
        180*a1*c1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) + 
        90*pow(b1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,3) - 28*b1*c1*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
        140*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 28*a1*d1*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 
        70*d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 336*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
        28*b1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 140*a1*c1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
        28*a1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 70*d*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) + 
        336*a1*b1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,3) - 32*a1*c1*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 
        192*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 16*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 
        280*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 32*a1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 
        192*a1*b1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 16*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,3) + 
        280*pow(a1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,3) - 36*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,3) - 
        126*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,3) + 36*a1*b1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,3) - 
        126*d*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,3) - 20*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,6)*pow(RP,3) + 
        20*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,6)*pow(RP,3) - 32*c1*d1*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
        60*b1*d*d1*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 192*b1*c1*d*k1*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
        60*b1*d1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 192*a1*d*d1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
        30*d*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 30*k1*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
        80*b1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 80*a1*d1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 
        360*a1*c1*k1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) - 180*k1*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,4) + 
        32*c1*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 120*b1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 
        60*d*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 368*b1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 
        368*a1*d1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 920*a1*c1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 
        460*pow(b1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 1968*a1*b1*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 
        1876*pow(a1,2)*pow(d,5)*pow(M_E,k1*pow(RP,-1))*pow(RP,4) + 32*c1*d1*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) - 
        16*b1*c1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) - 16*a1*d1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) + 
        80*a1*c1*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) + 40*pow(b1,2)*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) + 
        432*a1*b1*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) + 616*pow(a1,2)*pow(d,5)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,4) - 
        32*c1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 60*b1*d*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 
        192*b1*c1*d*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 60*b1*d1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 
        192*a1*d*d1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) + 30*d*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 
        30*k1*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 80*b1*c1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 
        80*a1*d1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 360*a1*c1*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 
        180*k1*pow(b1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,4) - 96*b1*c1*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
        420*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 96*a1*d1*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
        210*d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 1008*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
        96*b1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 420*a1*c1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
        96*a1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) + 210*d*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 
        1008*a1*b1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,4) - 140*a1*c1*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 
        768*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 70*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 
        1120*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 140*a1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,4) + 
        768*a1*b1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 70*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 
        1120*pow(a1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,4) - 192*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,4) - 
        630*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,4) - 192*a1*b1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,4) + 
        630*d*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,4) - 126*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,5)*pow(RP,4) - 
        126*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,5)*pow(RP,4) - 192*b1*c1*d*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 
        60*b1*d1*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 192*a1*d*d1*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 192*b1*c1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 
        840*a1*c1*d*k1*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 192*a1*d1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 
        420*d*k1*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 30*pow(c1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 
        360*a1*c1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 2016*a1*b1*k1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) - 
        180*pow(b1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,5) + 384*b1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 
        60*b1*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 384*a1*d*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 
        30*pow(c1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 1620*a1*c1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 
        810*pow(b1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 5088*a1*b1*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) + 
        6510*pow(a1,2)*pow(d,4)*pow(M_E,k1*pow(RP,-1))*pow(RP,5) - 60*b1*d1*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) - 
        30*pow(c1,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) + 60*a1*c1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) + 
        30*pow(b1,2)*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) - 480*a1*b1*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) - 
        1470*pow(a1,2)*pow(d,4)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,5) - 192*b1*c1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 
        60*b1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 192*a1*d*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 
        192*b1*c1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 840*a1*c1*d*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 
        192*a1*d1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 420*d*k1*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 
        30*pow(c1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 360*a1*c1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 
        2016*a1*b1*k1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) + 180*pow(b1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,5) - 
        420*a1*c1*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 2304*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 
        210*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 3360*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 
        420*a1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 2304*a1*b1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 
        210*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,5) + 
        3360*pow(a1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,5) - 768*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,5) - 
        2520*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,5) + 768*a1*b1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,5) - 
        2520*d*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,5) - 630*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,4)*pow(RP,5) + 
        630*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,4)*pow(RP,5) - 192*b1*c1*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 
        840*a1*c1*d*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 192*a1*d1*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 840*a1*c1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 
        4608*a1*b1*d*k1*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 420*d*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 
        420*k1*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 2016*a1*b1*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) - 
        6720*k1*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,6) + 192*b1*c1*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 
        1680*a1*c1*d*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 192*a1*d1*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 
        840*d*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 8928*a1*b1*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 
        16800*pow(a1,2)*pow(d,3)*pow(M_E,k1*pow(RP,-1))*pow(RP,6) + 192*b1*c1*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,6) + 
        192*a1*d1*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,6) - 288*a1*b1*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,6) + 
        1680*pow(a1,2)*pow(d,3)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,6) - 192*b1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) + 
        840*a1*c1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 192*a1*d1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 
        840*a1*c1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) + 4608*a1*b1*d*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) + 
        420*d*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 420*k1*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 
        2016*a1*b1*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 6720*k1*pow(a1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,6) - 
        2304*a1*b1*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,6) - 7560*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,6) - 
        2304*a1*b1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,6) + 7560*d*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,6) - 
        2520*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,3)*pow(RP,6) - 2520*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,3)*pow(RP,6) - 
        840*a1*c1*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 4608*a1*b1*d*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 
        4608*a1*b1*k1*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 15120*d*k1*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 
        420*pow(b1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) - 6720*pow(a1,2)*pow(d,2)*pow(M_E,d*pow(RP,-1))*pow(RP,7) + 
        840*a1*c1*pow(M_E,k1*pow(RP,-1))*pow(RP,7) + 9216*a1*b1*d*pow(M_E,k1*pow(RP,-1))*pow(RP,7) + 
        420*pow(b1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) + 29400*pow(a1,2)*pow(d,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,7) - 
        840*a1*c1*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,7) - 420*pow(b1,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,7) + 
        840*pow(a1,2)*pow(d,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,7) + 840*a1*c1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) - 
        4608*a1*b1*d*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) + 4608*a1*b1*k1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) - 
        15120*d*k1*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) + 420*pow(b1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) + 
        6720*pow(a1,2)*pow(d,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,7) - 7560*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(k1,2)*pow(RP,7) + 
        7560*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(k1,2)*pow(RP,7) - 4608*a1*b1*pow(M_E,d*pow(RP,-1))*pow(RP,8) - 
        15120*d*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,8) - 15120*k1*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,8) + 
        4608*a1*b1*pow(M_E,k1*pow(RP,-1))*pow(RP,8) + 30240*d*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,8) + 
        4608*a1*b1*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,8) - 4608*a1*b1*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,8) + 
        15120*d*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,8) - 15120*k1*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,8) - 
        15120*pow(a1,2)*pow(M_E,d*pow(RP,-1))*pow(RP,9) + 15120*pow(a1,2)*pow(M_E,k1*pow(RP,-1))*pow(RP,9) - 
        15120*pow(a1,2)*pow(M_E,(2*d + k1)*pow(RP,-1))*pow(RP,9) + 15120*pow(a1,2)*pow(M_E,(d + 2*k1)*pow(RP,-1))*pow(RP,9))*pow(2*pow(d,3),-1);
    return res;
}
/* *********** CASE 2 DONE *********** */

/*      Case_3:     d < k1 < k2         */
double SH_Integral_44_case_3( double AP, double RP, double k1, double k2, double a1, double b1, double c1, double d1, double d )
{
    double res=3*AP*RP*pow(M_E,-(d*pow(RP,-1)))*(RP*pow(d,2)*pow(d1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*
        (-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1))) + 
        2*(-gsl_sf_expint_Ei(-(k1*pow(RP,-1))) + gsl_sf_expint_Ei(-(k2*pow(RP,-1))))*pow(d,2)*pow(d1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,2) + 
        2*pow(d,2)*pow(d1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,-(k1*pow(RP,-1))) - pow(M_E,-(k2*pow(RP,-1))))*pow(RP,2) + 
        2*d*pow(d1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*
        pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2) + 4*c1*d1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*
        (-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2) + 
        2*c1*d1*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
        4*b1*d1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
        2*pow(c1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
        4*c1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2))) + 
        6*d*(gsl_sf_expint_Ei(-(k1*pow(RP,-1))) - gsl_sf_expint_Ei(-(k2*pow(RP,-1))))*pow(d1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,3) + 
        4*c1*d1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,-(k1*pow(RP,-1))) - pow(M_E,-(k2*pow(RP,-1))))*pow(RP,3) + 
        6*d*pow(d1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,k1*pow(RP,-1)) - pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3) + 
        12*c1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*
        pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3) + 4*b1*d1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*
        (-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3) + 
        2*pow(c1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*
        pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3) + 2*pow(d1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*
        (-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3) + 
        4*c1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,3) + 
        4*b1*c1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,3) + 
        4*a1*d1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,3) + 
        12*b1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,3)
        + 6*d*pow(c1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,3)
        + 2*b1*d1*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
        RP*pow(c1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
        4*b1*c1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
        4*a1*d1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
        4*b1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
        (pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
        2*pow(c1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
        (pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
        4*a1*c1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
        (pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
        2*pow(b1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
        (pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
        4*b1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
        2*d*pow(c1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
        12*b1*c1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
        (-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
        12*a1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
        (-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3))) + 
        6*(-gsl_sf_expint_Ei(-(k1*pow(RP,-1))) + gsl_sf_expint_Ei(-(k2*pow(RP,-1))))*pow(d1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(RP,4) + 
        6*pow(d1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,-(k1*pow(RP,-1))) - pow(M_E,-(k2*pow(RP,-1))))*pow(RP,4) + 
        12*c1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,k1*pow(RP,-1)) - pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4) + 
        12*b1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*
        pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4) + 6*d*pow(c1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*
        ((k2 + RP)*pow(M_E,k1*pow(RP,-1)) - (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4) + 
        12*c1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*
        pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4) + 12*b1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,4) + 
        6*pow(c1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,4) + 
        12*b1*c1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,4)
        + 12*a1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (-(pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2))) + pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,4)
        + 12*b1*c1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,4) + 
        12*a1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,4) + 
        12*a1*c1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,4) + 
        6*d*pow(b1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (-(pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3))) + 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,4) + 
        2*b1*c1*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
        pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
        2*a1*d1*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
        pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
        4*a1*c1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
        pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
        2*pow(b1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
        pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
        4*b1*c1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
        (pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
        pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
        4*a1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
        (pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
        pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
        4*a1*b1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
        (pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
        pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
        12*a1*c1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*
        (pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
        pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
        6*pow(b1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*
        (pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
        pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
        4*b1*c1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
        pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
        4*a1*d*d1*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
        pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
        12*a1*c1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
        (-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
        pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
        6*d*pow(b1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
        (-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
        pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
        12*a1*b1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*
        (-(pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4))) + 
        pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4))) + 
        12*c1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*(pow(M_E,-(k1*pow(RP,-1))) - pow(M_E,-(k2*pow(RP,-1))))*pow(RP,5) + 
        12*b1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*(-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*
        pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,5) + 6*pow(c1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*
        (-((k2 + RP)*pow(M_E,k1*pow(RP,-1))) + (k1 + RP)*pow(M_E,k2*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,5) + 
        12*b1*c1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,5) + 
        12*a1*d1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (pow(M_E,k2*pow(RP,-1))*(2*k1*RP + pow(k1,2) + 2*pow(RP,2)) - pow(M_E,k1*pow(RP,-1))*(2*k2*RP + pow(k2,2) + 2*pow(RP,2)))*pow(RP,5) + 
        12*a1*c1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,5) + 
        6*pow(b1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (pow(M_E,k2*pow(RP,-1))*(3*RP*pow(k1,2) + pow(k1,3) + 6*k1*pow(RP,2) + 6*pow(RP,3)) - 
        pow(M_E,k1*pow(RP,-1))*(3*RP*pow(k2,2) + pow(k2,3) + 6*k2*pow(RP,2) + 6*pow(RP,3)))*pow(RP,5) + 
        12*a1*b1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (pow(M_E,k2*pow(RP,-1))*(4*RP*pow(k1,3) + pow(k1,4) + 12*pow(k1,2)*pow(RP,2) + 24*k1*pow(RP,3) + 24*pow(RP,4)) - 
        pow(M_E,k1*pow(RP,-1))*(4*RP*pow(k2,3) + pow(k2,4) + 12*pow(k2,2)*pow(RP,2) + 24*k2*pow(RP,3) + 24*pow(RP,4)))*pow(RP,5) + 
        2*a1*c1*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 
        120*pow(RP,5)) - pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 
        120*k2*pow(RP,4) + 120*pow(RP,5))) + RP*pow(b1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 
        120*pow(RP,5)) - pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 
        120*k2*pow(RP,4) + 120*pow(RP,5))) + 4*a1*b1*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 
        120*pow(RP,5)) - pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 
        120*k2*pow(RP,4) + 120*pow(RP,5))) + 4*a1*c1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
        (pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 
        120*pow(RP,5)) - pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 
        120*k2*pow(RP,4) + 120*pow(RP,5))) + 2*pow(b1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
        (pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 
        120*pow(RP,5)) - pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 
        120*k2*pow(RP,4) + 120*pow(RP,5))) + 2*pow(a1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        pow(RP,3)*(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 
        120*k1*pow(RP,4) + 120*pow(RP,5)) - pow(M_E,k1*pow(RP,-1))*
        (5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 120*k2*pow(RP,4) + 120*pow(RP,5))) + 
        12*a1*b1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*
        (pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 
        120*pow(RP,5)) - pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 
        120*k2*pow(RP,4) + 120*pow(RP,5))) + 6*pow(a1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,5)*
        (pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 
        120*pow(RP,5)) - pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 
        120*k2*pow(RP,4) + 120*pow(RP,5))) + 4*a1*c1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 
        120*pow(RP,5))) + pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 
        120*k2*pow(RP,4) + 120*pow(RP,5))) + 2*d*pow(b1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 
        120*pow(RP,5))) + pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 
        120*k2*pow(RP,4) + 120*pow(RP,5))) + 12*a1*b1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
        (-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 
        120*pow(RP,5))) + pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 
        120*k2*pow(RP,4) + 120*pow(RP,5))) + 6*d*pow(a1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*
        (-(pow(M_E,k2*pow(RP,-1))*(5*RP*pow(k1,4) + pow(k1,5) + 20*pow(k1,3)*pow(RP,2) + 60*pow(k1,2)*pow(RP,3) + 120*k1*pow(RP,4) + 
        120*pow(RP,5))) + pow(M_E,k1*pow(RP,-1))*(5*RP*pow(k2,4) + pow(k2,5) + 20*pow(k2,3)*pow(RP,2) + 60*pow(k2,2)*pow(RP,3) + 
        120*k2*pow(RP,4) + 120*pow(RP,5))) + 2*a1*b1*RP*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 
        720*k1*pow(RP,5) + 720*pow(RP,6)) - pow(M_E,k1*pow(RP,-1))*
        (6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 720*k2*pow(RP,5) + 
        720*pow(RP,6))) + 2*pow(a1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 
        720*k1*pow(RP,5) + 720*pow(RP,6)) - pow(M_E,k1*pow(RP,-1))*
        (6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 720*k2*pow(RP,5) + 
        720*pow(RP,6))) + 4*a1*b1*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
        (pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 
        720*k1*pow(RP,5) + 720*pow(RP,6)) - pow(M_E,k1*pow(RP,-1))*
        (6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 720*k2*pow(RP,5) + 
        720*pow(RP,6))) + 6*pow(a1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,4)*
        (pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 360*pow(k1,2)*pow(RP,4) + 
        720*k1*pow(RP,5) + 720*pow(RP,6)) - pow(M_E,k1*pow(RP,-1))*
        (6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 720*k2*pow(RP,5) + 
        720*pow(RP,6))) + 4*a1*b1*d*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (-(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 
        360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 720*pow(RP,6))) + 
        pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 
        720*k2*pow(RP,5) + 720*pow(RP,6))) + 6*d*pow(a1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
        (-(pow(M_E,k2*pow(RP,-1))*(6*RP*pow(k1,5) + pow(k1,6) + 30*pow(k1,4)*pow(RP,2) + 120*pow(k1,3)*pow(RP,3) + 
        360*pow(k1,2)*pow(RP,4) + 720*k1*pow(RP,5) + 720*pow(RP,6))) + 
        pow(M_E,k1*pow(RP,-1))*(6*RP*pow(k2,5) + pow(k2,6) + 30*pow(k2,4)*pow(RP,2) + 120*pow(k2,3)*pow(RP,3) + 360*pow(k2,2)*pow(RP,4) + 
        720*k2*pow(RP,5) + 720*pow(RP,6))) + RP*pow(a1,2)*pow(d,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*
        (pow(M_E,k2*pow(RP,-1))*(7*RP*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RP,2) + 210*pow(k1,4)*pow(RP,3) + 840*pow(k1,3)*pow(RP,4) + 
        2520*pow(k1,2)*pow(RP,5) + 5040*k1*pow(RP,6) + 5040*pow(RP,7)) - 
        pow(M_E,k1*pow(RP,-1))*(7*RP*pow(k2,6) + pow(k2,7) + 42*pow(k2,5)*pow(RP,2) + 210*pow(k2,4)*pow(RP,3) + 840*pow(k2,3)*pow(RP,4) + 
        2520*pow(k2,2)*pow(RP,5) + 5040*k2*pow(RP,6) + 5040*pow(RP,7))) + 
        2*pow(a1,2)*(-1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,3)*
        (pow(M_E,k2*pow(RP,-1))*(7*RP*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RP,2) + 210*pow(k1,4)*pow(RP,3) + 840*pow(k1,3)*pow(RP,4) + 
        2520*pow(k1,2)*pow(RP,5) + 5040*k1*pow(RP,6) + 5040*pow(RP,7)) - 
        pow(M_E,k1*pow(RP,-1))*(7*RP*pow(k2,6) + pow(k2,7) + 42*pow(k2,5)*pow(RP,2) + 210*pow(k2,4)*pow(RP,3) + 840*pow(k2,3)*pow(RP,4) + 
        2520*pow(k2,2)*pow(RP,5) + 5040*k2*pow(RP,6) + 5040*pow(RP,7))) + 
        2*d*pow(a1,2)*(1 + pow(M_E,2*d*pow(RP,-1)))*pow(M_E,-((k1 + k2)*pow(RP,-1)))*pow(RP,2)*
        (-(pow(M_E,k2*pow(RP,-1))*(7*RP*pow(k1,6) + pow(k1,7) + 42*pow(k1,5)*pow(RP,2) + 210*pow(k1,4)*pow(RP,3) + 
        840*pow(k1,3)*pow(RP,4) + 2520*pow(k1,2)*pow(RP,5) + 5040*k1*pow(RP,6) + 5040*pow(RP,7))) + 
        pow(M_E,k1*pow(RP,-1))*(7*RP*pow(k2,6) + pow(k2,7) + 42*pow(k2,5)*pow(RP,2) + 210*pow(k2,4)*pow(RP,3) + 840*pow(k2,3)*pow(RP,4) + 
        2520*pow(k2,2)*pow(RP,5) + 5040*k2*pow(RP,6) + 5040*pow(RP,7))))*pow(2*pow(d,3),-1);
    return res;
}
/* *********** CASE 3 DONE *********** */

/* ******************************************************************************************     3. < Z | Kernel_SH | Z >        */




