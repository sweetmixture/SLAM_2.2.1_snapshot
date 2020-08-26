#ifndef __CH_INTEGRALS__
    #define __CH_INTEGRALS__

/* Functions Calculating CH Integrals 
 * by the kernel with spline cubinc polynomial set
 * of radial wave-functions
 */

/* !!! Integral Arguments
 *
 * k1: a knot for the left  boundary for spline cubic polynomial
 * k2: a knot for the right boundary for spline cubic polynomial
 *
 * d : the distance between sp-electrons and an external point charge
 *
 * a1, b1, c1, d1 : 1st cubic polynomial coefficients
 * a2, b2, c2, d2 : 2nd cubic polynomial coefficients if there is radial function of p state
 *
 */


//SS
double CH_Integral_11_case_1( double k1, double k2, double a1, double b1, double c1, double d1, double d );

double CH_Integral_11_case_2_sub_1( double k1, double k2, double a1, double b1, double c1, double d1, double d );

double CH_Integral_11_case_2_sub_2( double k1, double k2, double a1, double b1, double c1, double d1, double d );

double CH_Integral_11_case_3( double k1, double k2, double a1, double b1, double c1, double d1, double d );

//SZ
double CH_Integral_14_case_1( double k1, double k2, double a1, double b1, double c1, double d1,
        double a2, double b2, double c2, double d2, double d );
#define CH_Integral_41_case_1 CH_Integral_14_case_1

double CH_Integral_14_case_2_sub_1( double k1, double k2, double a1, double b1, double c1, double d1,
        double a2, double b2, double c2, double d2, double d );
#define CH_Integral_41_case_2_sub_1 CH_Integral_14_case_2_sub_1

double CH_Integral_14_case_2_sub_2( double k1, double k2, double a1, double b1, double c1, double d1,
        double a2, double b2, double c2, double d2, double d );
#define CH_Integral_41_case_2_sub_2 CH_Integral_14_case_2_sub_2

double CH_Integral_14_case_3( double k1, double k2, double a1, double b1, double c1, double d1,
        double a2, double b2, double c2, double d2, double d );
#define CH_Integral_41_case_3 CH_Integral_14_case_3

//XX YY
double CH_Integral_2233_case_1( double k1, double k2, double a1, double b1, double c1, double d1, double d );

double CH_Integral_2233_case_2_sub_1( double k1, double k2, double a1, double b1, double c1, double d1, double d );

double CH_Integral_2233_case_2_sub_2( double k1, double k2, double a1, double b1, double c1, double d1, double d );

double CH_Integral_2233_case_3( double k1, double k2, double a1, double b1, double c1, double d1, double d );

//ZZ
double CH_Integral_44_case_1( double k1, double k2, double a1, double b1, double c1, double d1, double d );

double CH_Integral_44_case_2_sub_1( double k1, double k2, double a1, double b1, double c1, double d1, double d );

double CH_Integral_44_case_2_sub_2( double k1, double k2, double a1, double b1, double c1, double d1, double d );

double CH_Integral_44_case_3( double k1, double k2, double a1, double b1, double c1, double d1, double d );








// x y and z integral

double Integral_x_12( double k1, double k2, double a1, double b1, double c1, double d1,
        double a2, double b2, double c2, double d2 );

double Integral_y_13( double k1, double k2, double a1, double b1, double c1, double d1,
        double a2, double b2, double c2, double d2 );

double Integral_z_14( double k1, double k2, double a1, double b1, double c1, double d1,
        double a2, double b2, double c2, double d2 );

#endif
