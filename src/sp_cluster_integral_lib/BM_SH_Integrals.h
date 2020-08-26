#ifndef __BM_SH_INTEGRALS__
    #define __BM_SH_INTEGRALS__

/* Functions Calculating SH Integrals 
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
double BMSH_Integral_11_case_1( double AS, double RS, double k1, double k2, double a1, double b1, double c1, double d1, double d );

double BMSH_Integral_11_case_2_sub_1( double AS, double RS, double k1, double k2, double a1, double b1, double c1, double d1, double d );

double BMSH_Integral_11_case_2_sub_2( double AS, double RS, double k1, double k2, double a1, double b1, double c1, double d1, double d );

double BMSH_Integral_11_case_3( double AS, double RS, double k1, double k2, double a1, double b1, double c1, double d1, double d );

//SZ
double BMSH_Integral_14_case_1( double ASP, double RSP, double k1, double k2, double a1, double b1, double c1, double d1,
        double a2, double b2, double c2, double d2, double d );
#define BMSH_Integral_41_case_1 BMSH_Integral_14_case_1

double BMSH_Integral_14_case_2_sub_1( double ASP, double RSP, double k1, double k2, double a1, double b1, double c1, double d1,
        double a2, double b2, double c2, double d2, double d );
#define BMSH_Integral_41_case_2_sub_1 BMSH_Integral_14_case_2_sub_1

double BMSH_Integral_14_case_2_sub_2( double ASP, double RSP, double k1, double k2, double a1, double b1, double c1, double d1,
        double a2, double b2, double c2, double d2, double d );
#define BMSH_Integral_41_case_2_sub_2 BMSH_Integral_14_case_2_sub_2

double BMSH_Integral_14_case_3( double ASP, double RSP, double k1, double k2, double a1, double b1, double c1, double d1,
        double a2, double b2, double c2, double d2, double d );
#define BMSH_Integral_41_case_3 BMSH_Integral_14_case_3

//XX YY
double BMSH_Integral_2233_case_1( double AP, double RP, double k1, double k2, double a1, double b1, double c1, double d1, double d );

double BMSH_Integral_2233_case_2_sub_1( double AP, double RP, double k1, double k2, double a1, double b1, double c1, double d1, double d );

double BMSH_Integral_2233_case_2_sub_2( double AP, double RP, double k1, double k2, double a1, double b1, double c1, double d1, double d );

double BMSH_Integral_2233_case_3( double AP, double RP, double k1, double k2, double a1, double b1, double c1, double d1, double d );

//ZZ
double BMSH_Integral_44_case_1( double AP, double RP, double k1, double k2, double a1, double b1, double c1, double d1, double d );

double BMSH_Integral_44_case_2_sub_1( double AP, double RP, double k1, double k2, double a1, double b1, double c1, double d1, double d );

double BMSH_Integral_44_case_2_sub_2( double AP, double RP, double k1, double k2, double a1, double b1, double c1, double d1, double d );

double BMSH_Integral_44_case_3( double AP, double RP, double k1, double k2, double a1, double b1, double c1, double d1, double d );

#endif
