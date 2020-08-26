#ifndef __BM_SDH_INTEGRAL__
    #define __BM_SDH_INTEGRAL__

// dxH SX == dyH SY
double BMSDH_Integral_x_12_case_1( double ASP, double RSP, double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d );
#define BMSDH_Integral_x_21_case_1    BMSDH_Integral_x_12_case_1
#define BMSDH_Integral_y_13_case_1    BMSDH_Integral_x_12_case_1
#define BMSDH_Integral_y_31_case_1    BMSDH_Integral_x_12_case_1

double BMSDH_Integral_x_12_case_2_sub_1( double ASP, double RSP, double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d );
#define BMSDH_Integral_x_21_case_2_sub_1  BMSDH_Integral_x_12_case_2_sub_1
#define BMSDH_Integral_y_13_case_2_sub_1  BMSDH_Integral_x_12_case_2_sub_1
#define BMSDH_Integral_y_31_case_2_sub_1  BMSDH_Integral_x_12_case_2_sub_1

double BMSDH_Integral_x_12_case_2_sub_2( double ASP, double RSP, double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d );
#define BMSDH_Integral_x_21_case_2_sub_2  BMSDH_Integral_x_12_case_2_sub_2
#define BMSDH_Integral_y_13_case_2_sub_2  BMSDH_Integral_x_12_case_2_sub_2
#define BMSDH_Integral_y_31_case_2_sub_2  BMSDH_Integral_x_12_case_2_sub_2

double BMSDH_Integral_x_12_case_3( double ASP, double RSP, double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d );
#define BMSDH_Integral_x_21_case_3    BMSDH_Integral_x_12_case_3
#define BMSDH_Integral_y_13_case_3    BMSDH_Integral_x_12_case_3
#define BMSDH_Integral_y_31_case_3    BMSDH_Integral_x_12_case_3


// dxH XZ == dyH YZ
double BMSDH_Integral_x_24_case_1( double AP, double RP, double k1, double k2, double a2, double b2, double c2, double d2, double d );
#define BMSDH_Integral_x_42_case_1    BMSDH_Integral_x_24_case_1
#define BMSDH_Integral_y_34_case_1    BMSDH_Integral_x_24_case_1
#define BMSDH_Integral_y_43_case_1    BMSDH_Integral_x_24_case_1

double BMSDH_Integral_x_24_case_2_sub_1( double AP, double RP, double k1, double k2, double a2, double b2, double c2, double d2, double d );
#define BMSDH_Integral_x_42_case_2_sub_1  BMSDH_Integral_x_24_case_2_sub_1
#define BMSDH_Integral_y_34_case_2_sub_1  BMSDH_Integral_x_24_case_2_sub_1
#define BMSDH_Integral_y_43_case_2_sub_1  BMSDH_Integral_x_24_case_2_sub_1

double BMSDH_Integral_x_24_case_2_sub_2( double AP, double RP, double k1, double k2, double a2, double b2, double c2, double d2, double d );
#define BMSDH_Integral_x_42_case_2_sub_2  BMSDH_Integral_x_24_case_2_sub_2
#define BMSDH_Integral_y_34_case_2_sub_2  BMSDH_Integral_x_24_case_2_sub_2
#define BMSDH_Integral_y_43_case_2_sub_2  BMSDH_Integral_x_24_case_2_sub_2

double BMSDH_Integral_x_24_case_3( double AP, double RP, double k1, double k2, double a2, double b2, double c2, double d2, double d );
#define BMSDH_Integral_x_42_case_3    BMSDH_Integral_x_24_case_3
#define BMSDH_Integral_y_34_case_3    BMSDH_Integral_x_24_case_3
#define BMSDH_Integral_y_43_case_3    BMSDH_Integral_x_24_case_3


// dzH SS
double BMSDH_Integral_z_11_case_1( double AS, double RS, double k1, double k2, double a1, double b1, double c1, double d1, double d );

double BMSDH_Integral_z_11_case_2_sub_1( double AS, double RS, double k1, double k2, double a1, double b1, double c1, double d1, double d );

double BMSDH_Integral_z_11_case_2_sub_2( double AS, double RS, double k1, double k2, double a1, double b1, double c1, double d1, double d );

double BMSDH_Integral_z_11_case_3( double AS, double RS, double k1, double k2, double a1, double b1, double c1, double d1, double d );

// dzH SZ
double BMSDH_Integral_z_14_case_1( double ASP, double RSP, double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d );
#define BMSDH_Integral_z_41_case_1    BMSDH_Integral_z_14_case_1

double BMSDH_Integral_z_14_case_2_sub_1( double ASP, double RSP, double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d );
#define BMSDH_Integral_z_41_case_2_sub_1    BMSDH_Integral_z_14_case_2_sub_1

double BMSDH_Integral_z_14_case_2_sub_2( double ASP, double RSP, double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d );
#define BMSDH_Integral_z_41_case_2_sub_2    BMSDH_Integral_z_14_case_2_sub_2

double BMSDH_Integral_z_14_case_3( double ASP, double RSP, double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d );
#define BMSDH_Integral_z_41_case_3    BMSDH_Integral_z_14_case_3

// dzH XX YY
double BMSDH_Integral_z_2233_case_1( double AP, double RP, double k1, double k2, double a2, double b2, double c2, double d2, double d );

double BMSDH_Integral_z_2233_case_2_sub_1( double AP, double RP, double k1, double k2, double a2, double b2, double c2, double d2, double d );

double BMSDH_Integral_z_2233_case_2_sub_2( double AP, double RP, double k1, double k2, double a2, double b2, double c2, double d2, double d );

double BMSDH_Integral_z_2233_case_3( double AP, double RP, double k1, double k2, double a2, double b2, double c2, double d2, double d );

// dzH ZZ
double BMSDH_Integral_z_44_case_1( double AP, double RP, double k1, double k2, double a2, double b2, double c2, double d2, double d );

double BMSDH_Integral_z_44_case_2_sub_1( double AP, double RP, double k1, double k2, double a2, double b2, double c2, double d2, double d );

double BMSDH_Integral_z_44_case_2_sub_2( double AP, double RP, double k1, double k2, double a2, double b2, double c2, double d2, double d );

double BMSDH_Integral_z_44_case_3( double AP, double RP, double k1, double k2, double a2, double b2, double c2, double d2, double d );



#endif
