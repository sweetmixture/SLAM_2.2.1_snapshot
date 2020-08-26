#ifndef __CDDH_INTEGRAL__
    #define __CDDH_INTEGRAL__



/* < S | dxxH | S > = < S | dyyH | S > */

double CDDH_Integral_xx_11_case_1( double k1, double k2, double a1, double b1, double c1, double d1, double d );
#define CDDH_Integral_yy_11_case_1  CDDH_Integral_xx_11_case_1
double CDDH_Integral_xx_11_case_2_sub_1( double k1, double k2, double a1, double b1, double c1, double d1, double d );
#define CDDH_Integral_yy_11_case_2_sub_1  CDDH_Integral_xx_11_case_2_sub_1
double CDDH_Integral_xx_11_case_2_sub_2( double k1, double k2, double a1, double b1, double c1, double d1, double d );
#define CDDH_Integral_yy_11_case_2_sub_2  CDDH_Integral_xx_11_case_2_sub_2
double CDDH_Integral_xx_11_case_3( double k1, double k2, double a1, double b1, double c1, double d1, double d );
#define CDDH_Integral_yy_11_case_3  CDDH_Integral_xx_11_case_3

/* < S | dxxH | S > done */




/* < S | dxxH | Z > = < S | dyy | Z > */

double CDDH_Integral_xx_14_case_1( double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d );
#define CDDH_Integral_yy_14_case_1  CDDH_Integral_xx_14_case_1
double CDDH_Integral_xx_14_case_2_sub_1( double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d );
#define CDDH_Integral_yy_14_case_2_sub_1  CDDH_Integral_xx_14_case_2_sub_1
double CDDH_Integral_xx_14_case_2_sub_2( double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d );
#define CDDH_Integral_yy_14_case_2_sub_2  CDDH_Integral_xx_14_case_2_sub_2
double CDDH_Integral_xx_14_case_3( double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d );
#define CDDH_Integral_yy_14_case_3  CDDH_Integral_xx_14_case_3

/* < S | dxxH | Z > done */


/* < X | dxxH | X > = < Y | dyyH | Y > */
double CDDH_Integral_xx_22_case_1( double k1, double k2, double a2, double b2, double c2, double d2, double d );
#define CDDH_Integral_yy_33_case_1  CDDH_Integral_xx_22_case_1
double CDDH_Integral_xx_22_case_2_sub_1( double k1, double k2, double a2, double b2, double c2, double d2, double d );
#define CDDH_Integral_yy_33_case_2_sub_1  CDDH_Integral_xx_22_case_2_sub_1
double CDDH_Integral_xx_22_case_2_sub_2( double k1, double k2, double a2, double b2, double c2, double d2, double d );
#define CDDH_Integral_yy_33_case_2_sub_2  CDDH_Integral_xx_22_case_2_sub_2
double CDDH_Integral_xx_22_case_3( double k1, double k2, double a2, double b2, double c2, double d2, double d );
#define CDDH_Integral_yy_33_case_3  CDDH_Integral_xx_22_case_3

/* < X | dxxH | X >  done */


/* < Y | dxxH | Y > = < X | dyyH | X > */
double CDDH_Integral_xx_33_case_1( double k1, double k2, double a2, double b2, double c2, double d2, double d );
#define CDDH_Integral_yy_22_case_1  CDDH_Integral_xx_33_case_1
double CDDH_Integral_xx_33_case_2_sub_1( double k1, double k2, double a2, double b2, double c2, double d2, double d );
#define CDDH_Integral_yy_22_case_2_sub_1  CDDH_Integral_xx_33_case_2_sub_1
double CDDH_Integral_xx_33_case_2_sub_2( double k1, double k2, double a2, double b2, double c2, double d2, double d );
#define CDDH_Integral_yy_22_case_2_sub_2  CDDH_Integral_xx_33_case_2_sub_12
double CDDH_Integral_xx_33_case_3( double k1, double k2, double a2, double b2, double c2, double d2, double d );
#define CDDH_Integral_yy_22_case_3  CDDH_Integral_xx_33_case_3

/* < Y | dxxH | Y >  done */



/* < Z | dxxH | Z >  == < Z | dyyH | Z > */
double CDDH_Integral_xx_44_case_1( double k1, double k2, double a2, double b2, double c2, double d2, double d );
#define CDDH_Integral_yy_44_case_1  CDDH_Integral_xx_44_case_1
double CDDH_Integral_xx_44_case_2_sub_1( double k1, double k2, double a2, double b2, double c2, double d2, double d );
#define CDDH_Integral_yy_44_case_2_sub_1  CDDH_Integral_xx_44_case_2_sub_1
double CDDH_Integral_xx_44_case_2_sub_2( double k1, double k2, double a2, double b2, double c2, double d2, double d );
#define CDDH_Integral_yy_44_case_2_sub_2  CDDH_Integral_xx_44_case_2_sub_2
double CDDH_Integral_xx_44_case_3( double k1, double k2, double a2, double b2, double c2, double d2, double d );
#define CDDH_Integral_yy_44_case_3  CDDH_Integral_xx_44_case_3

/////////////////////////////////////////////////////// dxxH & dyyH




/* *************************************************** dxyH */

/* < X | dxyH | Y > */
double CDDH_Integral_xy_23_case_1( double k1, double k2, double a2, double b2, double c2, double d2, double d );

double CDDH_Integral_xy_23_case_2_sub_1( double k1, double k2, double a2, double b2, double c2, double d2, double d );

double CDDH_Integral_xy_23_case_2_sub_2( double k1, double k2, double a2, double b2, double c2, double d2, double d );

double CDDH_Integral_xy_23_case_3( double k1, double k2, double a2, double b2, double c2, double d2, double d );
/* < X | dxyH | Y >  done */




/* ************************************************** dxzH */

/* < S | dxzH | X >  = < S | dyz | Y > */
double CDDH_Integral_xz_12_case_1( double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d );
#define CDDH_Integral_yz_13_case_1  CDDH_Integral_xz_12_case_1
double CDDH_Integral_xz_12_case_2_sub_1( double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d );
#define CDDH_Integral_yz_13_case_2_sub_1  CDDH_Integral_xz_12_case_2_sub_1
double CDDH_Integral_xz_12_case_2_sub_2( double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d );
#define CDDH_Integral_yz_13_case_2_sub_2  CDDH_Integral_xz_12_case_2_sub_2
double CDDH_Integral_xz_12_case_3( double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d );
#define CDDH_Integral_yz_13_case_3  CDDH_Integral_xz_12_case_3
/* < S | dxzH | X >  done */


/* < X | dxzH | Z >  = < Y | dyz | Z > */
double CDDH_Integral_xz_24_case_1( double k1, double k2, double a2, double b2, double c2, double d2, double d );
#define CDDH_Integral_yz_34_case_1  CDDH_Integral_xz_24_case_1
double CDDH_Integral_xz_24_case_2_sub_1( double k1, double k2, double a2, double b2, double c2, double d2, double d );
#define CDDH_Integral_yz_34_case_2_sub_1  CDDH_Integral_xz_24_case_2_sub_1
double CDDH_Integral_xz_24_case_2_sub_2( double k1, double k2, double a2, double b2, double c2, double d2, double d );
#define CDDH_Integral_yz_34_case_2_sub_2  CDDH_Integral_xz_24_case_2_sub_2
double CDDH_Integral_xz_24_case_3( double k1, double k2, double a2, double b2, double c2, double d2, double d );
#define CDDH_Integral_yz_34_case_3  CDDH_Integral_xz_24_case_3
/* < X | dxzH | Z >  done */





/* < S | dzzH | S > */

double CDDH_Integral_zz_11_case_1( double k1, double k2, double a1, double b1, double c1, double d1, double d );

double CDDH_Integral_zz_11_case_2_sub_1( double k1, double k2, double a1, double b1, double c1, double d1, double d );

double CDDH_Integral_zz_11_case_2_sub_2( double k1, double k2, double a1, double b1, double c1, double d1, double d );

double CDDH_Integral_zz_11_case_3( double k1, double k2, double a1, double b1, double c1, double d1, double d );

/* < S | dzzH | S > done */




/* < S | dzzH | Z > */

double CDDH_Integral_zz_14_case_1( double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d );

double CDDH_Integral_zz_14_case_2_sub_1( double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d );

double CDDH_Integral_zz_14_case_2_sub_2( double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d );

double CDDH_Integral_zz_14_case_3( double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d );

/* < S | dzzH | Z > done */


/* < X | dzzH | X > = < Y | dzzH | Y > */
double CDDH_Integral_zz_2233_case_1( double k1, double k2, double a2, double b2, double c2, double d2, double d );

double CDDH_Integral_zz_2233_case_2_sub_1( double k1, double k2, double a2, double b2, double c2, double d2, double d );

double CDDH_Integral_zz_2233_case_2_sub_2( double k1, double k2, double a2, double b2, double c2, double d2, double d );

double CDDH_Integral_zz_2233_case_3( double k1, double k2, double a2, double b2, double c2, double d2, double d );

/* < X | dzzH | X >  done */


/* < Z | dzzH | Z >  */
double CDDH_Integral_zz_44_case_1( double k1, double k2, double a2, double b2, double c2, double d2, double d );

double CDDH_Integral_zz_44_case_2_sub_1( double k1, double k2, double a2, double b2, double c2, double d2, double d );

double CDDH_Integral_zz_44_case_2_sub_2( double k1, double k2, double a2, double b2, double c2, double d2, double d );

double CDDH_Integral_zz_44_case_3( double k1, double k2, double a2, double b2, double c2, double d2, double d );

#endif
