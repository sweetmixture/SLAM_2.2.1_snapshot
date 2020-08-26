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



/* ******************************************************************************************     1. < S | Kernel_CH | S >        */

/*      Case_1:    k1 < k2 < d          */
double CH_Integral_11_case_1( double k1, double k2, double a1, double b1, double c1, double d1, double d )
{   double res=pow(1260*d,-1)*(420*pow(d1,2)*(-pow(k1,3) + pow(k2,3)) + 252*pow(c1,2)*(-pow(k1,5) + pow(k2,5)) - 
            42*d1*(15*c1*(pow(k1,4) - pow(k2,4)) + 2*(6*b1*(pow(k1,5) - pow(k2,5)) + 5*a1*(pow(k1,6) - pow(k2,6)))) - 
            60*c1*(7*b1*(pow(k1,6) - pow(k2,6)) + 6*a1*(pow(k1,7) - pow(k2,7))) + 5*(36*pow(b1,2)*(-pow(k1,7) + pow(k2,7)) + 
            63*a1*b1*(-pow(k1,8) + pow(k2,8)) + 28*pow(a1,2)*(-pow(k1,9) + pow(k2,9))));
    return res;
}
/* *********** CASE 1 DONE *********** */


/*      Case_2:     k1 < d < k2         */
/*      Case_2_Sub_1:   d < r < k2      */
double CH_Integral_11_case_2_sub_1( double k1, double k2, double a1, double b1, double c1, double d1, double d )
{   double res=-(pow(4,-1)*pow(c1,2)*(pow(d,4) - pow(k2,4))) - c1*pow(15,-1)*
        (5*(pow(d,3) - pow(k2,3))*(2*d1 + a1*(pow(d,3) + pow(k2,3))) + 6*b1*(pow(d,5) - pow(k2,5))) - 
        pow(6,-1)*pow(b1,2)*(pow(d,6) - pow(k2,6)) - b1*pow(14,-1)*(7*d1*(pow(d,4) - pow(k2,4)) + 4*a1*(pow(d,7) - pow(k2,7))) + 
        pow(40,-1)*(20*pow(d1,2)*(-pow(d,2) + pow(k2,2)) - 16*a1*d1*(pow(d,5) - pow(k2,5)) - 5*pow(a1,2)*(pow(d,8) - pow(k2,8)));
    return res;
}
/*      Case_2_Sub_2:   k1 < r < d      */
double CH_Integral_11_case_2_sub_2( double k1, double k2, double a1, double b1, double c1, double d1, double d )
{   double res=pow(1260*d,-1)*(252*pow(c1,2)*(pow(d,5) - pow(k1,5)) + 180*pow(b1,2)*(pow(d,7) - pow(k1,7)) + 
        30*c1*(14*b1*(pow(d,6) - pow(k1,6)) + 3*(7*d1*pow(d,4) + 4*a1*pow(d,7) - 7*d1*pow(k1,4) - 4*a1*pow(k1,7))) + 
        63*b1*(8*d1*(pow(d,5) - pow(k1,5)) + 5*a1*(pow(d,8) - pow(k1,8))) + 
        140*(3*pow(d1,2)*(pow(d,3) - pow(k1,3)) + 3*a1*d1*(pow(d,6) - pow(k1,6)) + pow(a1,2)*(pow(d,9) - pow(k1,9))));
    return res;   
}
/* *********** CASE 2 DONE *********** */


/*      Case_3:     d < k1 < k2         */
double CH_Integral_11_case_3( double k1, double k2, double a1, double b1, double c1, double d1, double d )
{   double res=pow(840,-1)*(-420*pow(d1,2)*(pow(k1,2) - pow(k2,2)) - 210*pow(c1,2)*(pow(k1,4) - pow(k2,4)) - 
        28*d1*(20*c1*(pow(k1,3) - pow(k2,3)) + 3*(5*b1*pow(k1,4) + 4*a1*pow(k1,5) - 5*b1*pow(k2,4) - 4*a1*pow(k2,5))) - 
        56*c1*(6*b1*(pow(k1,5) - pow(k2,5)) + 5*a1*(pow(k1,6) - pow(k2,6))) - 
        5*(28*pow(b1,2)*(pow(k1,6) - pow(k2,6)) + 48*a1*b1*(pow(k1,7) - pow(k2,7)) + 21*pow(a1,2)*(pow(k1,8) - pow(k2,8))));
    return res;
}
/* *********** CASE 3 DONE *********** */

/* ******************************************************************************************     1. < S | Kernel_CH | S >        */


        /*      *       *       *       *       *       *       *       *       *       *       *       *       */


/* ******************************************************************************************     2. < S | Kernel_CH | Z >        */

/*      Case_1:    k1 < k2 < d          */
double CH_Integral_14_case_1( double k1, double k2, double a1, double b1, double c1, double d1,
        double a2, double b2, double c2, double d2, double d )
{   double res=(-(d1*d2*pow(2,-1)*pow(k1,4)) - 2*(c2*d1 + c1*d2)*pow(5,-1)*pow(k1,5) - (c1*c2 + b2*d1 + b1*d2)*pow(3,-1)*pow(k1,6) -
        2*(b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(7,-1)*pow(k1,7) - (b1*b2 + a2*c1 + a1*c2)*pow(4,-1)*pow(k1,8) -
        2*(a2*b1 + a1*b2)*pow(9,-1)*pow(k1,9) - a1*a2*pow(5,-1)*pow(k1,10) + d1*d2*pow(2,-1)*pow(k2,4) +
        2*(c2*d1 + c1*d2)*pow(5,-1)*pow(k2,5) + (c1*c2 + b2*d1 + b1*d2)*pow(3,-1)*pow(k2,6) +
        2*(b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(7,-1)*pow(k2,7) + (b1*b2 + a2*c1 + a1*c2)*pow(4,-1)*pow(k2,8) +
        2*(a2*b1 + a1*b2)*pow(9,-1)*pow(k2,9) + a1*a2*pow(5,-1)*pow(k2,10))*pow(2*pow(3,0.5)*pow(d,2),-1);
    return res;
}
/* *********** CASE 1 DONE *********** */

/*      Case_2:     k1 < d < k2         */
/*      Case_2_Sub_1:   d < r < k2      */
double CH_Integral_14_case_2_sub_1( double k1, double k2, double a1, double b1, double c1, double d1,
        double a2, double b2, double c2, double d2, double d )
{   double res=d*(420*d1*d2*(-d + k2) + 210*c2*d1*(-pow(d,2) + pow(k2,2)) + 210*c1*d2*(-pow(d,2) + pow(k2,2)) + 
        140*c1*c2*(-pow(d,3) + pow(k2,3)) + 140*b2*d1*(-pow(d,3) + pow(k2,3)) + 140*b1*d2*(-pow(d,3) + pow(k2,3)) + 
        105*b2*c1*(-pow(d,4) + pow(k2,4)) + 105*b1*c2*(-pow(d,4) + pow(k2,4)) + 105*a2*d1*(-pow(d,4) + pow(k2,4)) + 
        105*a1*d2*(-pow(d,4) + pow(k2,4)) + 84*b1*b2*(-pow(d,5) + pow(k2,5)) + 84*a2*c1*(-pow(d,5) + pow(k2,5)) + 
        84*a1*c2*(-pow(d,5) + pow(k2,5)) + 70*a2*b1*(-pow(d,6) + pow(k2,6)) + 70*a1*b2*(-pow(d,6) + pow(k2,6)) + 
        60*a1*a2*(-pow(d,7) + pow(k2,7)))*pow(420*pow(3,0.5),-1);
    return res;
}
/*      Case_2_Sub_2:   k1 < r < d      */
double CH_Integral_14_case_2_sub_2( double k1, double k2, double a1, double b1, double c1, double d1,
        double a2, double b2, double c2, double d2, double d )
{   double res=(d1*d2*pow(2,-1)*pow(d,4) + 2*c2*d1*pow(5,-1)*pow(d,5) + b2*d1*pow(3,-1)*pow(d,6) + 
        b1*pow(252,-1)*(72*c2*d + 84*d2 + 7*(9*b2 + 8*a2*d)*pow(d,2))*pow(d,6) + 2*a2*d1*pow(7,-1)*pow(d,7) + 
        2*a1*d2*pow(7,-1)*pow(d,7) + a1*c2*pow(4,-1)*pow(d,8) + c1*(2*d2*pow(5,-1)*pow(d,5) + c2*pow(3,-1)*pow(d,6) + 
        2*b2*pow(7,-1)*pow(d,7) + a2*pow(4,-1)*pow(d,8)) + 2*a1*b2*pow(9,-1)*pow(d,9) + a1*a2*pow(5,-1)*pow(d,10) - 
        d1*d2*pow(2,-1)*pow(k1,4) - 2*(c2*d1 + c1*d2)*pow(5,-1)*pow(k1,5) - (c1*c2 + b2*d1 + b1*d2)*pow(3,-1)*pow(k1,6) - 
        2*(b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(7,-1)*pow(k1,7) - (b1*b2 + a2*c1 + a1*c2)*pow(4,-1)*pow(k1,8) - 
        2*(a2*b1 + a1*b2)*pow(9,-1)*pow(k1,9) - a1*a2*pow(5,-1)*pow(k1,10))*pow(2*pow(3,0.5)*pow(d,2),-1);
    return res;
}
/* *********** CASE 2 DONE *********** */

/*      Case_3:     d < k1 < k2         */
double CH_Integral_14_case_3( double k1, double k2, double a1, double b1, double c1, double d1,
        double a2, double b2, double c2, double d2, double d )
{   double res=d*(420*d1*d2*k2 + k1*(-420*d1*d2 - 210*(c2*d1 + c1*d2)*k1 - 140*(c1*c2 + b2*d1 + b1*d2)*pow(k1,2) - 
        105*(b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(k1,3) - 84*(b1*b2 + a2*c1 + a1*c2)*pow(k1,4) - 70*(a2*b1 + a1*b2)*pow(k1,5) - 
        60*a1*a2*pow(k1,6)) + 210*(c2*d1 + c1*d2)*pow(k2,2) + 140*(c1*c2 + b2*d1 + b1*d2)*pow(k2,3) + 
        105*(b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(k2,4) + 84*(b1*b2 + a2*c1 + a1*c2)*pow(k2,5) + 70*(a2*b1 + a1*b2)*pow(k2,6) + 
        60*a1*a2*pow(k2,7))*pow(420*pow(3,0.5),-1);
    return res;
}
/* *********** CASE 3 DONE *********** */

/* ******************************************************************************************     2. < S | Kernel_CH | Z >        */


        /*      *       *       *       *       *       *       *       *       *       *       *       *       */


/* ******************************************************************************************     3. < XorY | Kernel_CH | XorY >        */

/*      Case_1:    k1 < k2 < d          */
double CH_Integral_2233_case_1( double k1, double k2, double a1, double b1, double c1, double d1, double d )
{   double res=(1980*pow(c1,2)*pow(k1,7) + 3465*b1*c1*pow(k1,8) + 3080*a1*c1*pow(k1,9) + 1540*pow(b1,2)*pow(k1,9) + 
        2772*a1*b1*pow(k1,10) + 1260*pow(a1,2)*pow(k1,11) + 2772*pow(d1,2)*(pow(k1,5) - pow(k2,5)) - 
        1980*pow(c1,2)*pow(k2,7) + 165*d1*(28*c1*(pow(k1,6) - pow(k2,6)) + 3*(8*b1*(pow(k1,7) - pow(k2,7)) + 
        7*a1*(pow(k1,8) - pow(k2,8)))) - 3465*b1*c1*pow(k2,8) - 55*pow(d,2)*(420*pow(d1,2)*(pow(k1,3) - 
        pow(k2,3)) + 252*pow(c1,2)*(pow(k1,5) - pow(k2,5)) + 42*d1*(15*c1*(pow(k1,4) - pow(k2,4)) + 
        2*(6*b1*pow(k1,5) + 5*a1*pow(k1,6) - 6*b1*pow(k2,5) - 5*a1*pow(k2,6))) + 60*c1*(7*b1*(pow(k1,6) - 
        pow(k2,6)) + 6*a1*(pow(k1,7) - pow(k2,7))) + 5*(36*pow(b1,2)*(pow(k1,7) - pow(k2,7)) + 
        63*a1*b1*(pow(k1,8) - pow(k2,8)) + 28*pow(a1,2)*(pow(k1,9) - pow(k2,9)))) - 3080*a1*c1*pow(k2,9) - 
        1540*pow(b1,2)*pow(k2,9) - 2772*a1*b1*pow(k2,10) - 1260*pow(a1,2)*pow(k2,11))*pow(69300*pow(d,3),-1);
    return res;
}
/* *********** CASE 1 DONE *********** */

/*      Case_2:     k1 < d < k2         */
/*      Case_2_Sub_1:   d < r < k2      */
double CH_Integral_2233_case_2_sub_1( double k1, double k2, double a1, double b1, double c1, double d1, double d )
{   double res=pow(5,-1)*(-2*c1*d1*k2*pow(d,2) - 4*c1*d1*pow(3,-1)*pow(d,3) - 3*b1*d1*pow(2,-1)*pow(d,4) - 
        3*pow(4,-1)*pow(c1,2)*pow(d,4) - 4*b1*c1*pow(3,-1)*pow(d,5) - 4*a1*d1*pow(3,-1)*pow(d,5) - 
        7*a1*c1*pow(6,-1)*pow(d,6) - 7*pow(12,-1)*pow(b1,2)*pow(d,6) - 36*a1*b1*pow(35,-1)*pow(d,7) - 
        11*pow(24,-1)*pow(a1,2)*pow(d,8) + log(d)*pow(d,2)*pow(d1,2) - log(k2)*pow(d,2)*pow(d1,2) - 
        5*pow(2,-1)*pow(d,2)*pow(d1,2) - b1*d1*pow(d,2)*pow(k2,2) - pow(2,-1)*pow(c1,2)*pow(d,2)*pow(k2,2) + 
        5*pow(2,-1)*pow(d1,2)*pow(k2,2) + 10*c1*d1*pow(3,-1)*pow(k2,3) - 2*b1*c1*pow(3,-1)*pow(d,2)*pow(k2,3) - 
        2*a1*d1*pow(3,-1)*pow(d,2)*pow(k2,3) + 5*b1*d1*pow(2,-1)*pow(k2,4) + 5*pow(4,-1)*pow(c1,2)*pow(k2,4) - 
        a1*c1*pow(2,-1)*pow(d,2)*pow(k2,4) - pow(4,-1)*pow(b1,2)*pow(d,2)*pow(k2,4) + 2*b1*c1*pow(k2,5) + 
        2*a1*d1*pow(k2,5) - 2*a1*b1*pow(5,-1)*pow(d,2)*pow(k2,5) + 5*a1*c1*pow(3,-1)*pow(k2,6) + 
        5*pow(6,-1)*pow(b1,2)*pow(k2,6) - pow(6,-1)*pow(a1,2)*pow(d,2)*pow(k2,6) + 10*a1*b1*pow(7,-1)*pow(k2,7) + 
        5*pow(8,-1)*pow(a1,2)*pow(k2,8));
    return res;
}
/*      Case_2_Sub_2:   k1 < r < d      */
double CH_Integral_2233_case_2_sub_2( double k1, double k2, double a1, double b1, double c1, double d1, double d )
{   double res=(1980*pow(c1,2)*(6*pow(d,7) - 7*pow(d,2)*pow(k1,5) + pow(k1,7)) + 220*pow(b1,2)*(38*pow(d,9) - 
        45*pow(d,2)*pow(k1,7) + 7*pow(k1,9)) + 55*c1*(42*d1*(13*pow(d,6) - 15*pow(d,2)*pow(k1,4) + 
        2*pow(k1,6)) + 21*b1*(17*pow(d,8) - 20*pow(d,2)*pow(k1,6) + 3*pow(k1,8)) + 8*a1*(38*pow(d,9) - 
        45*pow(d,2)*pow(k1,7) + 7*pow(k1,9))) + 99*b1*(40*d1*(6*pow(d,7) - 7*pow(d,2)*pow(k1,5) + 
        pow(k1,7)) + 7*a1*(21*pow(d,10) - 25*pow(d,2)*pow(k1,8) + 4*pow(k1,10))) + 
        7*(132*pow(d1,2)*(22*pow(d,5) - 25*pow(d,2)*pow(k1,3) + 3*pow(k1,5)) + 165*a1*d1*(17*pow(d,8) - 
        20*pow(d,2)*pow(k1,6) + 3*pow(k1,8)) + 20*pow(a1,2)*(46*pow(d,11) - 55*pow(d,2)*pow(k1,9) + 
        9*pow(k1,11))))*pow(69300*pow(d,3),-1);
    return res;
}
/* *********** CASE 2 DONE *********** */

/*      Case_3:     d < k1 < k2         */
double CH_Integral_2233_case_3( double k1, double k2, double a1, double b1, double c1, double d1, double d )
{   double res=pow(5,-1)*(2*c1*d1*k1*pow(d,2) - 2*c1*d1*k2*pow(d,2) + log(k1)*pow(d,2)*pow(d1,2) - 
        log(k2)*pow(d,2)*pow(d1,2) + b1*d1*pow(d,2)*pow(k1,2) + pow(2,-1)*pow(c1,2)*pow(d,2)*pow(k1,2) - 
        5*pow(2,-1)*pow(d1,2)*pow(k1,2) - 10*c1*d1*pow(3,-1)*pow(k1,3) + 2*b1*c1*pow(3,-1)*pow(d,2)*pow(k1,3) + 
        2*a1*d1*pow(3,-1)*pow(d,2)*pow(k1,3) - 5*b1*d1*pow(2,-1)*pow(k1,4) - 5*pow(4,-1)*pow(c1,2)*pow(k1,4) + 
        a1*c1*pow(2,-1)*pow(d,2)*pow(k1,4) + pow(4,-1)*pow(b1,2)*pow(d,2)*pow(k1,4) - 2*b1*c1*pow(k1,5) - 
        2*a1*d1*pow(k1,5) + 2*a1*b1*pow(5,-1)*pow(d,2)*pow(k1,5) - 5*a1*c1*pow(3,-1)*pow(k1,6) - 
        5*pow(6,-1)*pow(b1,2)*pow(k1,6) + pow(6,-1)*pow(a1,2)*pow(d,2)*pow(k1,6) - 10*a1*b1*pow(7,-1)*pow(k1,7) - 
        5*pow(8,-1)*pow(a1,2)*pow(k1,8) - b1*d1*pow(d,2)*pow(k2,2) - pow(2,-1)*pow(c1,2)*pow(d,2)*pow(k2,2) + 
        5*pow(2,-1)*pow(d1,2)*pow(k2,2) + 10*c1*d1*pow(3,-1)*pow(k2,3) - 2*b1*c1*pow(3,-1)*pow(d,2)*pow(k2,3) - 
        2*a1*d1*pow(3,-1)*pow(d,2)*pow(k2,3) + 5*b1*d1*pow(2,-1)*pow(k2,4) + 5*pow(4,-1)*pow(c1,2)*pow(k2,4) - 
        a1*c1*pow(2,-1)*pow(d,2)*pow(k2,4) - pow(4,-1)*pow(b1,2)*pow(d,2)*pow(k2,4) + 2*b1*c1*pow(k2,5) + 
        2*a1*d1*pow(k2,5) - 2*a1*b1*pow(5,-1)*pow(d,2)*pow(k2,5) + 5*a1*c1*pow(3,-1)*pow(k2,6) + 
        5*pow(6,-1)*pow(b1,2)*pow(k2,6) - pow(6,-1)*pow(a1,2)*pow(d,2)*pow(k2,6) + 
        10*a1*b1*pow(7,-1)*pow(k2,7) + 5*pow(8,-1)*pow(a1,2)*pow(k2,8));
    return res;
}
/* *********** CASE 3 DONE *********** */


/* ******************************************************************************************     3. < XorY | Kernel_CH | XorY >        */

        /*      *       *       *       *       *       *       *       *       *       *       *       *       */

/* ******************************************************************************************     4. < Z | Kernel_CH | Z >        */

/*      Case_1:    k1 < k2 < d          */
double CH_Integral_44_case_1( double k1, double k2, double a1, double b1, double c1, double d1, double d )
{   double res=-((3960*pow(c1,2)*pow(k1,7) + 6930*b1*c1*pow(k1,8) + 6160*a1*c1*pow(k1,9) + 3080*pow(b1,2)*pow(k1,9) + 
        5544*a1*b1*pow(k1,10) + 2520*pow(a1,2)*pow(k1,11) + 5544*pow(d1,2)*(pow(k1,5) - pow(k2,5)) - 3960*pow(c1,2)*pow(k2,7) + 
        330*d1*(28*c1*(pow(k1,6) - pow(k2,6)) + 3*(8*b1*(pow(k1,7) - pow(k2,7)) + 7*a1*(pow(k1,8) - pow(k2,8)))) - 6930*b1*c1*pow(k2,8) + 
        55*pow(d,2)*(420*pow(d1,2)*(pow(k1,3) - pow(k2,3)) + 252*pow(c1,2)*(pow(k1,5) - pow(k2,5)) + 42*d1*(15*c1*(pow(k1,4) - 
        pow(k2,4)) + 2*(6*b1*pow(k1,5) + 5*a1*pow(k1,6) - 6*b1*pow(k2,5) - 5*a1*pow(k2,6))) + 60*c1*(7*b1*(pow(k1,6) - 
        pow(k2,6)) + 6*a1*(pow(k1,7) - pow(k2,7))) + 5*(36*pow(b1,2)*(pow(k1,7) - pow(k2,7)) + 63*a1*b1*(pow(k1,8) - 
        pow(k2,8)) + 28*pow(a1,2)*(pow(k1,9) - pow(k2,9)))) - 6160*a1*c1*pow(k2,9) - 3080*pow(b1,2)*pow(k2,9) - 
        5544*a1*b1*pow(k2,10) - 2520*pow(a1,2)*pow(k2,11))*pow(69300*pow(d,3),-1));
    return res;
}
/* *********** CASE 1 DONE *********** */

/*      Case_2:     k1 < d < k2         */
/*      Case_2_Sub_1:   d < r < k2      */
double CH_Integral_44_case_2_sub_1( double k1, double k2, double a1, double b1, double c1, double d1, double d )
{   double res=pow(5,-1)*(4*c1*d1*k2*pow(d,2) - 22*c1*d1*pow(3,-1)*pow(d,3) - 9*b1*d1*pow(2,-1)*pow(d,4) - 
        9*pow(4,-1)*pow(c1,2)*pow(d,4) - 10*b1*c1*pow(3,-1)*pow(d,5) - 10*a1*d1*pow(3,-1)*pow(d,5) - 
        8*a1*c1*pow(3,-1)*pow(d,6) - 4*pow(3,-1)*pow(b1,2)*pow(d,6) - 78*a1*b1*pow(35,-1)*pow(d,7) - 
        23*pow(24,-1)*pow(a1,2)*pow(d,8) - 2*log(d)*pow(d,2)*pow(d1,2) + 2*log(k2)*pow(d,2)*pow(d1,2) - 
        5*pow(2,-1)*pow(d,2)*pow(d1,2) + 2*b1*d1*pow(d,2)*pow(k2,2) + pow(c1,2)*pow(d,2)*pow(k2,2) + 
        5*pow(2,-1)*pow(d1,2)*pow(k2,2) + 10*c1*d1*pow(3,-1)*pow(k2,3) + 4*b1*c1*pow(3,-1)*pow(d,2)*pow(k2,3) + 
        4*a1*d1*pow(3,-1)*pow(d,2)*pow(k2,3) + 5*b1*d1*pow(2,-1)*pow(k2,4) + 5*pow(4,-1)*pow(c1,2)*pow(k2,4) + 
        a1*c1*pow(d,2)*pow(k2,4) + pow(2,-1)*pow(b1,2)*pow(d,2)*pow(k2,4) + 2*b1*c1*pow(k2,5) + 
        2*a1*d1*pow(k2,5) + 4*a1*b1*pow(5,-1)*pow(d,2)*pow(k2,5) + 5*a1*c1*pow(3,-1)*pow(k2,6) + 
        5*pow(6,-1)*pow(b1,2)*pow(k2,6) + pow(3,-1)*pow(a1,2)*pow(d,2)*pow(k2,6) + 
        10*a1*b1*pow(7,-1)*pow(k2,7) + 5*pow(8,-1)*pow(a1,2)*pow(k2,8));
    return res;
}
/*      Case_2_Sub_2:   k1 < r < d      */
double CH_Integral_44_case_2_sub_2( double k1, double k2, double a1, double b1, double c1, double d1, double d )
{   double res=(1980*pow(c1,2)*(9*pow(d,7) - 7*pow(d,2)*pow(k1,5) - 2*pow(k1,7)) + 110*c1*(21*d1*(19*pow(d,6) - 
        15*pow(d,2)*pow(k1,4) - 4*pow(k1,6)) + 21*b1*(13*pow(d,8) - 10*pow(d,2)*pow(k1,6) - 3*pow(k1,8)) + 
        4*a1*(59*pow(d,9) - 45*pow(d,2)*pow(k1,7) - 14*pow(k1,9))) + 220*pow(b1,2)*(59*pow(d,9) - 
        45*pow(d,2)*pow(k1,7) - 14*pow(k1,9)) + 99*b1*(40*d1*(9*pow(d,7) - 7*pow(d,2)*pow(k1,5) - 
        2*pow(k1,7)) + 7*a1*(33*pow(d,10) - 25*pow(d,2)*pow(k1,8) - 8*pow(k1,10))) + 
        14*(66*pow(d1,2)*(31*pow(d,5) - 25*pow(d,2)*pow(k1,3) - 6*pow(k1,5)) + 
        165*a1*d1*(13*pow(d,8) - 10*pow(d,2)*pow(k1,6) - 3*pow(k1,8)) + 
        10*pow(a1,2)*(73*pow(d,11) - 55*pow(d,2)*pow(k1,9) - 
        18*pow(k1,11))))*pow(69300*pow(d,3),-1);
    return res;
}
/* *********** CASE 2 DONE *********** */

/*      Case_3:     d < k1 < k2         */
double CH_Integral_44_case_3( double k1, double k2, double a1, double b1, double c1, double d1, double d )
{   double res=pow(5,-1)*(-4*c1*d1*k1*pow(d,2) + 4*c1*d1*k2*pow(d,2) - 2*log(k1)*pow(d,2)*pow(d1,2) + 
        2*log(k2)*pow(d,2)*pow(d1,2) - 2*b1*d1*pow(d,2)*pow(k1,2) - pow(c1,2)*pow(d,2)*pow(k1,2) - 
        5*pow(2,-1)*pow(d1,2)*pow(k1,2) - 10*c1*d1*pow(3,-1)*pow(k1,3) - 4*b1*c1*pow(3,-1)*pow(d,2)*pow(k1,3) - 
        4*a1*d1*pow(3,-1)*pow(d,2)*pow(k1,3) - 5*b1*d1*pow(2,-1)*pow(k1,4) - 5*pow(4,-1)*pow(c1,2)*pow(k1,4) - 
        a1*c1*pow(d,2)*pow(k1,4) - pow(2,-1)*pow(b1,2)*pow(d,2)*pow(k1,4) - 2*b1*c1*pow(k1,5) - 2*a1*d1*pow(k1,5) - 
        4*a1*b1*pow(5,-1)*pow(d,2)*pow(k1,5) - 5*a1*c1*pow(3,-1)*pow(k1,6) - 5*pow(6,-1)*pow(b1,2)*pow(k1,6) - 
        pow(3,-1)*pow(a1,2)*pow(d,2)*pow(k1,6) - 10*a1*b1*pow(7,-1)*pow(k1,7) - 5*pow(8,-1)*pow(a1,2)*pow(k1,8) + 
        2*b1*d1*pow(d,2)*pow(k2,2) + pow(c1,2)*pow(d,2)*pow(k2,2) + 5*pow(2,-1)*pow(d1,2)*pow(k2,2) + 
        10*c1*d1*pow(3,-1)*pow(k2,3) + 4*b1*c1*pow(3,-1)*pow(d,2)*pow(k2,3) + 4*a1*d1*pow(3,-1)*pow(d,2)*pow(k2,3) + 
        5*b1*d1*pow(2,-1)*pow(k2,4) + 5*pow(4,-1)*pow(c1,2)*pow(k2,4) + a1*c1*pow(d,2)*pow(k2,4) + 
        pow(2,-1)*pow(b1,2)*pow(d,2)*pow(k2,4) + 2*b1*c1*pow(k2,5) + 2*a1*d1*pow(k2,5) + 
        4*a1*b1*pow(5,-1)*pow(d,2)*pow(k2,5) + 5*a1*c1*pow(3,-1)*pow(k2,6) + 5*pow(6,-1)*pow(b1,2)*pow(k2,6) + 
        pow(3,-1)*pow(a1,2)*pow(d,2)*pow(k2,6) + 10*a1*b1*pow(7,-1)*pow(k2,7) + 5*pow(8,-1)*pow(a1,2)*pow(k2,8));
    return res;
}
/* *********** CASE 3 DONE *********** */

/* ******************************************************************************************     3. < Z | Kernel_CH | Z >        */

















/* x,y and z integrals 
 *
 * x = r sin[the] cos[pi]
 * y = r sin[the] sin[pi]
 * z = r cos[the]
 *
 *  Note that the integral results are in Bohr Unit!
 * 
 */


/* x */

// SX
double Integral_x_12( double k1, double k2, double a1, double b1, double c1, double d1,
double a2, double b2, double c2, double d2 )
{
double res =(-(d1*d2*pow(4,-1)*pow(k1,4)) - (c2*d1 + c1*d2)*pow(5,-1)*pow(k1,5) - (c1*c2 + b2*d1 + b1*d2)*pow(6,-1)*pow(k1,6) - 
(b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(7,-1)*pow(k1,7) - (b1*b2 + a2*c1 + a1*c2)*pow(8,-1)*pow(k1,8) - 
(a2*b1 + a1*b2)*pow(9,-1)*pow(k1,9) - a1*a2*pow(10,-1)*pow(k1,10) + d1*d2*pow(4,-1)*pow(k2,4) + 
(c2*d1 + c1*d2)*pow(5,-1)*pow(k2,5) + (c1*c2 + b2*d1 + b1*d2)*pow(6,-1)*pow(k2,6) + 
(b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(7,-1)*pow(k2,7) + (b1*b2 + a2*c1 + a1*c2)*pow(8,-1)*pow(k2,8) + 
(a2*b1 + a1*b2)*pow(9,-1)*pow(k2,9) + a1*a2*pow(10,-1)*pow(k2,10))*pow(pow(3,0.5),-1);

return res;
}

/* y */

// SY
double Integral_y_13( double k1, double k2, double a1, double b1, double c1, double d1,
double a2, double b2, double c2, double d2 )
{
double res =(-(d1*d2*pow(4,-1)*pow(k1,4)) - (c2*d1 + c1*d2)*pow(5,-1)*pow(k1,5) - (c1*c2 + b2*d1 + b1*d2)*pow(6,-1)*pow(k1,6) - 
(b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(7,-1)*pow(k1,7) - (b1*b2 + a2*c1 + a1*c2)*pow(8,-1)*pow(k1,8) - 
(a2*b1 + a1*b2)*pow(9,-1)*pow(k1,9) - a1*a2*pow(10,-1)*pow(k1,10) + d1*d2*pow(4,-1)*pow(k2,4) + 
(c2*d1 + c1*d2)*pow(5,-1)*pow(k2,5) + (c1*c2 + b2*d1 + b1*d2)*pow(6,-1)*pow(k2,6) + 
(b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(7,-1)*pow(k2,7) + (b1*b2 + a2*c1 + a1*c2)*pow(8,-1)*pow(k2,8) + 
(a2*b1 + a1*b2)*pow(9,-1)*pow(k2,9) + a1*a2*pow(10,-1)*pow(k2,10))*pow(pow(3,0.5),-1);

return res;
}

/* z */
double Integral_z_14( double k1, double k2, double a1, double b1, double c1, double d1,
double a2, double b2, double c2, double d2 )
{
double res =(-(d1*d2*pow(4,-1)*pow(k1,4)) - (c2*d1 + c1*d2)*pow(5,-1)*pow(k1,5) - (c1*c2 + b2*d1 + b1*d2)*pow(6,-1)*pow(k1,6) - 
(b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(7,-1)*pow(k1,7) - (b1*b2 + a2*c1 + a1*c2)*pow(8,-1)*pow(k1,8) - 
(a2*b1 + a1*b2)*pow(9,-1)*pow(k1,9) - a1*a2*pow(10,-1)*pow(k1,10) + d1*d2*pow(4,-1)*pow(k2,4) + 
(c2*d1 + c1*d2)*pow(5,-1)*pow(k2,5) + (c1*c2 + b2*d1 + b1*d2)*pow(6,-1)*pow(k2,6) + 
(b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(7,-1)*pow(k2,7) + (b1*b2 + a2*c1 + a1*c2)*pow(8,-1)*pow(k2,8) + 
(a2*b1 + a1*b2)*pow(9,-1)*pow(k2,9) + a1*a2*pow(10,-1)*pow(k2,10))*pow(pow(3,0.5),-1);

return res;
}

