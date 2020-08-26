#ifndef __CDH_INTEGRAL__
    #define __CDH_INTEGRAL__

// dxH SX == dyH SY
double CDH_Integral_x_12_case_1( double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d );
#define CDH_Integral_x_21_case_1    CDH_Integral_x_12_case_1
#define CDH_Integral_y_13_case_1    CDH_Integral_x_12_case_1
#define CDH_Integral_y_31_case_1    CDH_Integral_x_12_case_1

double CDH_Integral_x_12_case_2_sub_1( double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d );
#define CDH_Integral_x_21_case_2_sub_1  CDH_Integral_x_12_case_2_sub_1
#define CDH_Integral_y_13_case_2_sub_1  CDH_Integral_x_12_case_2_sub_1
#define CDH_Integral_y_31_case_2_sub_1  CDH_Integral_x_12_case_2_sub_1

double CDH_Integral_x_12_case_2_sub_2( double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d );
#define CDH_Integral_x_21_case_2_sub_2  CDH_Integral_x_12_case_2_sub_2
#define CDH_Integral_y_13_case_2_sub_2  CDH_Integral_x_12_case_2_sub_2
#define CDH_Integral_y_31_case_2_sub_2  CDH_Integral_x_12_case_2_sub_2

double CDH_Integral_x_12_case_3( double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d );
#define CDH_Integral_x_21_case_3    CDH_Integral_x_12_case_3
#define CDH_Integral_y_13_case_3    CDH_Integral_x_12_case_3
#define CDH_Integral_y_31_case_3    CDH_Integral_x_12_case_3


// dxH XZ == dyH YZ
double CDH_Integral_x_24_case_1( double k1, double k2, double a2, double b2, double c2, double d2, double d );
#define CDH_Integral_x_42_case_1    CDH_Integral_x_24_case_1
#define CDH_Integral_y_34_case_1    CDH_Integral_x_24_case_1
#define CDH_Integral_y_43_case_1    CDH_Integral_x_24_case_1

double CDH_Integral_x_24_case_2_sub_1( double k1, double k2, double a2, double b2, double c2, double d2, double d );
#define CDH_Integral_x_42_case_2_sub_1  CDH_Integral_x_24_case_2_sub_1
#define CDH_Integral_y_34_case_2_sub_1  CDH_Integral_x_24_case_2_sub_1
#define CDH_Integral_y_43_case_2_sub_1  CDH_Integral_x_24_case_2_sub_1

double CDH_Integral_x_24_case_2_sub_2( double k1, double k2, double a2, double b2, double c2, double d2, double d );
#define CDH_Integral_x_42_case_2_sub_2  CDH_Integral_x_24_case_2_sub_2
#define CDH_Integral_y_34_case_2_sub_2  CDH_Integral_x_24_case_2_sub_2
#define CDH_Integral_y_43_case_2_sub_2  CDH_Integral_x_24_case_2_sub_2

double CDH_Integral_x_24_case_3( double k1, double k2, double a2, double b2, double c2, double d2, double d );
#define CDH_Integral_x_42_case_3    CDH_Integral_x_24_case_3
#define CDH_Integral_y_34_case_3    CDH_Integral_x_24_case_3
#define CDH_Integral_y_43_case_3    CDH_Integral_x_24_case_3


// dzH SS
double CDH_Integral_z_11_case_1( double k1, double k2, double a1, double b1, double c1, double d1, double d );

double CDH_Integral_z_11_case_2_sub_1( double k1, double k2, double a1, double b1, double c1, double d1, double d );

double CDH_Integral_z_11_case_2_sub_2( double k1, double k2, double a1, double b1, double c1, double d1, double d );

double CDH_Integral_z_11_case_3( double k1, double k2, double a1, double b1, double c1, double d1, double d );

// dzH SZ
double CDH_Integral_z_14_case_1( double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d );
#define CDH_Integral_z_41_case_1    CDH_Integral_z_14_case_1

double CDH_Integral_z_14_case_2_sub_1( double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d );
#define CDH_Integral_z_41_case_2_sub_1    CDH_Integral_z_14_case_2_sub_1

double CDH_Integral_z_14_case_2_sub_2( double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d );
#define CDH_Integral_z_41_case_2_sub_2    CDH_Integral_z_14_case_2_sub_2

double CDH_Integral_z_14_case_3( double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d );
#define CDH_Integral_z_41_case_3    CDH_Integral_z_14_case_3

// dzH XX YY
double CDH_Integral_z_2233_case_1( double k1, double k2, double a2, double b2, double c2, double d2, double d );

double CDH_Integral_z_2233_case_2_sub_1( double k1, double k2, double a2, double b2, double c2, double d2, double d );

double CDH_Integral_z_2233_case_2_sub_2( double k1, double k2, double a2, double b2, double c2, double d2, double d );

double CDH_Integral_z_2233_case_3( double k1, double k2, double a2, double b2, double c2, double d2, double d );

// dzH ZZ
double CDH_Integral_z_44_case_1( double k1, double k2, double a2, double b2, double c2, double d2, double d );

double CDH_Integral_z_44_case_2_sub_1( double k1, double k2, double a2, double b2, double c2, double d2, double d );

double CDH_Integral_z_44_case_2_sub_2( double k1, double k2, double a2, double b2, double c2, double d2, double d );

double CDH_Integral_z_44_case_3( double k1, double k2, double a2, double b2, double c2, double d2, double d );



#endif
