#include<stdio.h>
#include<gsl/gsl_sf_expint.h>
#include<gsl/gsl_math.h>

// 111 _ xxx
//
//
// 111 - xxx SX

double CDDDH_Integral_xxx_12_case_1( double k1, double k2, double a1, double b1, double c1, double d1,
double a2, double b2, double c2, double d2, double d )
{   double res =3*pow(3,0.5)*pow(d,-5)*(-(d1*d2*pow(k1,4))/4. - ((c2*d1 + c1*d2)*pow(k1,5))/5. - ((c1*c2 + b2*d1 + b1*d2)*pow(k1,6))/6. - 
((b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(k1,7))/7. - ((b1*b2 + a2*c1 + a1*c2)*pow(k1,8))/8. - 
((a2*b1 + a1*b2)*pow(k1,9))/9. - (a1*a2*pow(k1,10))/10. + (d1*d2*pow(k2,4))/4. + ((c2*d1 + c1*d2)*pow(k2,5))/5. + 
((c1*c2 + b2*d1 + b1*d2)*pow(k2,6))/6. + ((b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(k2,7))/7. + 
((b1*b2 + a2*c1 + a1*c2)*pow(k2,8))/8. + ((a2*b1 + a1*b2)*pow(k2,9))/9. + (a1*a2*pow(k2,10))/10.);

return res;
}
double CDDDH_Integral_xxx_12_case_2( double k1, double k2, double a1, double b1, double c1, double d1,
double a2, double b2, double c2, double d2, double d )
{   double res =0.;

return res;
}

double CDDDH_Integral_xxx_12_case_3( double k1, double k2, double a1, double b1, double c1, double d1,
double a2, double b2, double c2, double d2, double d )
{   double res =(pow(3,-0.5)*pow(d,-5)*(pow(d,4)*(6*d2*(84*c1*d + 105*d1 + 70*b1*pow(d,2) + 60*a1*pow(d,3)) + 
d*(420*c1*c2*d + 504*c2*d1 + 60*d*(7*b2 + 6*a2*d)*d1 + 360*(b2*c1 + b1*c2)*pow(d,2) + 
315*(b1*b2 + a2*c1 + a1*c2)*pow(d,3) + 280*(a2*b1 + a1*b2)*pow(d,4) + 252*a1*a2*pow(d,5))) - 
630*d1*d2*pow(k1,4) - 504*(c2*d1 + c1*d2)*pow(k1,5) - 420*(c1*c2 + b2*d1 + b1*d2)*pow(k1,6) - 
360*(b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(k1,7) - 315*(b1*b2 + a2*c1 + a1*c2)*pow(k1,8) - 
280*(a2*b1 + a1*b2)*pow(k1,9) - 252*a1*a2*pow(k1,10)))/280.;

return res;
}

double CDDDH_Integral_xxx_12_case_4( double k1, double k2, double a1, double b1, double c1, double d1,
double a2, double b2, double c2, double d2, double d )
{   double res =0.;

return res;
}

// 111 - xxx SX Done
//
//
// 111 - xxx XZ
double CDDDH_Integral_xxx_24_case_1( double k1, double k2, double a2, double b2, double c2, double d2, double d)
{   double res =-(pow(d,-6)*(2772*pow(d2,2)*(pow(k1,5) - pow(k2,5)) + 1980*pow(c2,2)*(pow(k1,7) - pow(k2,7)) + 
165*d2*(28*c2*(pow(k1,6) - pow(k2,6)) + 24*b2*(pow(k1,7) - pow(k2,7)) + 21*a2*(pow(k1,8) - pow(k2,8))) + 
385*c2*(9*b2*(pow(k1,8) - pow(k2,8)) + 8*a2*(pow(k1,9) - pow(k2,9))) + 
28*(55*pow(b2,2)*(pow(k1,9) - pow(k2,9)) + 99*a2*b2*(pow(k1,10) - pow(k2,10)) + 
45*pow(a2,2)*(pow(k1,11) - pow(k2,11)))))/1540.;

return res;
}

double CDDDH_Integral_xxx_24_case_2( double k1, double k2, double a2, double b2, double c2, double d2, double d)
{   double res =0.;

return res;
}

double CDDDH_Integral_xxx_24_case_3( double k1, double k2, double a2, double b2, double c2, double d2, double d)
{   double res =(pow(d,-6)*(1980*pow(c2,2)*(pow(d,7) - pow(k1,7)) + 
385*c2*(12*d2*(pow(d,6) - pow(k1,6)) + 9*b2*(pow(d,8) - pow(k1,8)) + 8*a2*(pow(d,9) - pow(k1,9))) + 
1540*pow(b2,2)*(pow(d,9) - pow(k1,9)) + 396*b2*(10*d2*(pow(d,7) - pow(k1,7)) + 7*a2*(pow(d,10) - pow(k1,10))) + 
63*(44*pow(d2,2)*(pow(d,5) - pow(k1,5)) + 55*a2*d2*(pow(d,8) - pow(k1,8)) + 20*pow(a2,2)*(pow(d,11) - pow(k1,11)))))/
1540.;

return res;
}

double CDDDH_Integral_xxx_24_case_4( double k1, double k2, double a2, double b2, double c2, double d2, double d)
{   double res =0.;

return res;
}

// 111 - xxx XZ Done
//
//
// 112 - xxy
//
//
// 112 - xxy SY

double CDDDH_Integral_xxy_13_case_1( double k1, double k2, double a1, double b1, double c1, double d1,
double a2, double b2, double c2, double d2, double d)
{   double res =pow(3,0.5)*pow(d,-5)*(-(d1*d2*pow(k1,4))/4. - ((c2*d1 + c1*d2)*pow(k1,5))/5. - ((c1*c2 + b2*d1 + b1*d2)*pow(k1,6))/6. - 
((b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(k1,7))/7. - ((b1*b2 + a2*c1 + a1*c2)*pow(k1,8))/8. - 
((a2*b1 + a1*b2)*pow(k1,9))/9. - (a1*a2*pow(k1,10))/10. + (d1*d2*pow(k2,4))/4. + ((c2*d1 + c1*d2)*pow(k2,5))/5. + 
((c1*c2 + b2*d1 + b1*d2)*pow(k2,6))/6. + ((b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(k2,7))/7. + 
((b1*b2 + a2*c1 + a1*c2)*pow(k2,8))/8. + ((a2*b1 + a1*b2)*pow(k2,9))/9. + (a1*a2*pow(k2,10))/10.);

return res;
}

double CDDDH_Integral_xxy_13_case_2( double k1, double k2, double a1, double b1, double c1, double d1,
double a2, double b2, double c2, double d2, double d)
{   double res =0.;

return res;
}

double CDDDH_Integral_xxy_13_case_3( double k1, double k2, double a1, double b1, double c1, double d1,
double a2, double b2, double c2, double d2, double d)
{   double res =(pow(3,-0.5)*pow(d,-5)*(pow(d,4)*(6*d2*(84*c1*d + 105*d1 + 70*b1*pow(d,2) + 60*a1*pow(d,3)) + 
d*(420*c1*c2*d + 504*c2*d1 + 60*d*(7*b2 + 6*a2*d)*d1 + 360*(b2*c1 + b1*c2)*pow(d,2) + 
315*(b1*b2 + a2*c1 + a1*c2)*pow(d,3) + 280*(a2*b1 + a1*b2)*pow(d,4) + 252*a1*a2*pow(d,5))) - 
630*d1*d2*pow(k1,4) - 504*(c2*d1 + c1*d2)*pow(k1,5) - 420*(c1*c2 + b2*d1 + b1*d2)*pow(k1,6) - 
360*(b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(k1,7) - 315*(b1*b2 + a2*c1 + a1*c2)*pow(k1,8) - 
280*(a2*b1 + a1*b2)*pow(k1,9) - 252*a1*a2*pow(k1,10)))/840.;

return res;
}

double CDDDH_Integral_xxy_13_case_4( double k1, double k2, double a1, double b1, double c1, double d1,
double a2, double b2, double c2, double d2, double d)
{   double res =0.;

return res;
}

// 112 - xxy SY Done
//
//
// 112 - xxy YZ

double CDDDH_Integral_xxy_34_case_1( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{   double res =3*pow(d,-6)*((pow(d2,2)*(-pow(k1,5) + pow(k2,5)))/5. + (c2*d2*(-pow(k1,6) + pow(k2,6)))/3. + 
(2*b2*d2*(-pow(k1,7) + pow(k2,7)))/7. + (pow(c2,2)*(-pow(k1,7) + pow(k2,7)))/7. + 
(b2*c2*(-pow(k1,8) + pow(k2,8)))/4. + (a2*d2*(-pow(k1,8) + pow(k2,8)))/4. + (2*a2*c2*(-pow(k1,9) + pow(k2,9)))/9. + 
(pow(b2,2)*(-pow(k1,9) + pow(k2,9)))/9. + (a2*b2*(-pow(k1,10) + pow(k2,10)))/5. + 
(pow(a2,2)*(-pow(k1,11) + pow(k2,11)))/11.);

return res;
}

double CDDDH_Integral_xxy_34_case_2( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{   double res =0.;

return res;
}

double CDDDH_Integral_xxy_34_case_3( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{   double res =(pow(d,-6)*(1980*pow(c2,2)*(pow(d,7) - pow(k1,7)) + 
385*c2*(12*d2*(pow(d,6) - pow(k1,6)) + 9*b2*(pow(d,8) - pow(k1,8)) + 8*a2*(pow(d,9) - pow(k1,9))) + 
1540*pow(b2,2)*(pow(d,9) - pow(k1,9)) + 396*b2*(10*d2*(pow(d,7) - pow(k1,7)) + 7*a2*(pow(d,10) - pow(k1,10))) + 
63*(44*pow(d2,2)*(pow(d,5) - pow(k1,5)) + 55*a2*d2*(pow(d,8) - pow(k1,8)) + 20*pow(a2,2)*(pow(d,11) - pow(k1,11)))))/
4620.;

return res;
}

double CDDDH_Integral_xxy_34_case_4( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{   double res =0.;

return res;
}

// 112 - xxy YZ Done
//
//
// 113 - xxz SS

double CDDDH_Integral_xxz_11_case_1( double k1, double k2, double a1, double b1, double c1, double d1, double d )
{   double res =(pow(d,-4)*(420*pow(d1,2)*(pow(k1,3) - pow(k2,3)) + 252*pow(c1,2)*(pow(k1,5) - pow(k2,5)) + 
42*d1*(15*c1*(pow(k1,4) - pow(k2,4)) + 12*b1*(pow(k1,5) - pow(k2,5)) + 10*a1*(pow(k1,6) - pow(k2,6))) + 
60*c1*(7*b1*(pow(k1,6) - pow(k2,6)) + 6*a1*(pow(k1,7) - pow(k2,7))) + 
5*(36*pow(b1,2)*(pow(k1,7) - pow(k2,7)) + 63*a1*b1*(pow(k1,8) - pow(k2,8)) + 28*pow(a1,2)*(pow(k1,9) - pow(k2,9)))))/
420.;

return res;
}

double CDDDH_Integral_xxz_11_case_2( double k1, double k2, double a1, double b1, double c1, double d1, double d )
{   double res =0.;

return res;
}

double CDDDH_Integral_xxz_11_case_3( double k1, double k2, double a1, double b1, double c1, double d1, double d )
{   double res =-(pow(d,-4)*(252*pow(c1,2)*(pow(d,5) - pow(k1,5)) + 
30*c1*(21*d1*(pow(d,4) - pow(k1,4)) + 14*b1*(pow(d,6) - pow(k1,6)) + 12*a1*(pow(d,7) - pow(k1,7))) + 
180*pow(b1,2)*(pow(d,7) - pow(k1,7)) + 63*b1*(8*d1*(pow(d,5) - pow(k1,5)) + 5*a1*(pow(d,8) - pow(k1,8))) + 
140*(3*pow(d1,2)*(pow(d,3) - pow(k1,3)) + 3*a1*d1*(pow(d,6) - pow(k1,6)) + pow(a1,2)*(pow(d,9) - pow(k1,9)))))/420.;

return res;
}

double CDDDH_Integral_xxz_11_case_4( double k1, double k2, double a1, double b1, double c1, double d1, double d )
{   double res =0.;

return res;
}

// 113 - xxz SS Done
//
//
// 113 - xxz SZ

double CDDDH_Integral_xxz_14_case_1( double k1, double k2, double a1, double b1, double c1, double d1,
double a2, double b2, double c2, double d2, double d )
{   double res =-4*pow(3,0.5)*pow(d,-5)*(-(d1*d2*pow(k1,4))/4. - ((c2*d1 + c1*d2)*pow(k1,5))/5. - ((c1*c2 + b2*d1 + b1*d2)*pow(k1,6))/6. - 
((b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(k1,7))/7. - ((b1*b2 + a2*c1 + a1*c2)*pow(k1,8))/8. - 
((a2*b1 + a1*b2)*pow(k1,9))/9. - (a1*a2*pow(k1,10))/10. + (d1*d2*pow(k2,4))/4. + ((c2*d1 + c1*d2)*pow(k2,5))/5. + 
((c1*c2 + b2*d1 + b1*d2)*pow(k2,6))/6. + ((b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(k2,7))/7. + 
((b1*b2 + a2*c1 + a1*c2)*pow(k2,8))/8. + ((a2*b1 + a1*b2)*pow(k2,9))/9. + (a1*a2*pow(k2,10))/10.);

return res;
}

double CDDDH_Integral_xxz_14_case_2( double k1, double k2, double a1, double b1, double c1, double d1,
double a2, double b2, double c2, double d2, double d )
{   double res =0.;

return res;
}

double CDDDH_Integral_xxz_14_case_3( double k1, double k2, double a1, double b1, double c1, double d1,
double a2, double b2, double c2, double d2, double d )
{   double res =(pow(3,-0.5)*pow(d,-5)*(pow(d,4)*(-6*d2*(84*c1*d + 105*d1 + 70*b1*pow(d,2) + 60*a1*pow(d,3)) + 
d*(-12*(42*c2 + 5*d*(7*b2 + 6*a2*d))*d1 + 
d*(-420*c1*c2 - 360*(b2*c1 + b1*c2)*d - 315*(b1*b2 + a2*c1 + a1*c2)*pow(d,2) - 280*(a2*b1 + a1*b2)*pow(d,3) - 
252*a1*a2*pow(d,4)))) + 630*d1*d2*pow(k1,4) + 504*(c2*d1 + c1*d2)*pow(k1,5) + 
420*(c1*c2 + b2*d1 + b1*d2)*pow(k1,6) + 360*(b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(k1,7) + 
315*(b1*b2 + a2*c1 + a1*c2)*pow(k1,8) + 280*(a2*b1 + a1*b2)*pow(k1,9) + 252*a1*a2*pow(k1,10)))/210.;

return res;
}

double CDDDH_Integral_xxz_14_case_4( double k1, double k2, double a1, double b1, double c1, double d1,
double a2, double b2, double c2, double d2, double d )
{   double res =0.;

return res;
}

// 113 - xxz SZ Done
//
//
// 113 - xxz XX

double CDDDH_Integral_xxz_22_case_1( double k1, double k2, double a2, double b2, double c2, double d2, double d ) 
{   double res =(pow(d,-6)*(11*pow(d,2)*(420*pow(d2,2)*(pow(k1,3) - pow(k2,3)) + 252*pow(c2,2)*(pow(k1,5) - pow(k2,5)) + 
42*d2*(15*c2*(pow(k1,4) - pow(k2,4)) + 12*b2*(pow(k1,5) - pow(k2,5)) + 10*a2*(pow(k1,6) - pow(k2,6))) + 
60*c2*(7*b2*(pow(k1,6) - pow(k2,6)) + 6*a2*(pow(k1,7) - pow(k2,7))) + 
5*(36*pow(b2,2)*(pow(k1,7) - pow(k2,7)) + 63*a2*b2*(pow(k1,8) - pow(k2,8)) + 28*pow(a2,2)*(pow(k1,9) - pow(k2,9)))
) - 3*(2772*pow(d2,2)*(pow(k1,5) - pow(k2,5)) + 1980*pow(c2,2)*(pow(k1,7) - pow(k2,7)) + 
165*d2*(28*c2*(pow(k1,6) - pow(k2,6)) + 24*b2*(pow(k1,7) - pow(k2,7)) + 21*a2*(pow(k1,8) - pow(k2,8))) + 
385*c2*(9*b2*(pow(k1,8) - pow(k2,8)) + 8*a2*(pow(k1,9) - pow(k2,9))) + 
28*(55*pow(b2,2)*(pow(k1,9) - pow(k2,9)) + 99*a2*b2*(pow(k1,10) - pow(k2,10)) + 
45*pow(a2,2)*(pow(k1,11) - pow(k2,11))))))/4620.;

return res;
}

double CDDDH_Integral_xxz_22_case_2( double k1, double k2, double a2, double b2, double c2, double d2, double d ) 
{   double res =0.;

return res;
}

double CDDDH_Integral_xxz_22_case_3( double k1, double k2, double a2, double b2, double c2, double d2, double d ) 
{   double res =(pow(d,-6)*(396*pow(c2,2)*(8*pow(d,7) + 7*pow(d,2)*pow(k1,5) - 15*pow(k1,7)) + 
165*c2*(42*d2*(pow(d,6) + pow(d,2)*pow(k1,4) - 2*pow(k1,6)) + 
7*b2*(5*pow(d,8) + 4*pow(d,2)*pow(k1,6) - 9*pow(k1,8)) + 8*a2*(4*pow(d,9) + 3*pow(d,2)*pow(k1,7) - 7*pow(k1,9)))\
+ 660*pow(b2,2)*(4*pow(d,9) + 3*pow(d,2)*pow(k1,7) - 7*pow(k1,9)) + 
99*b2*(8*d2*(8*pow(d,7) + 7*pow(d,2)*pow(k1,5) - 15*pow(k1,7)) + 
7*a2*(7*pow(d,10) + 5*pow(d,2)*pow(k1,8) - 12*pow(k1,10))) + 
7*(132*pow(d2,2)*(4*pow(d,5) + 5*pow(d,2)*pow(k1,3) - 9*pow(k1,5)) + 
165*a2*d2*(5*pow(d,8) + 4*pow(d,2)*pow(k1,6) - 9*pow(k1,8)) + 
20*pow(a2,2)*(16*pow(d,11) + 11*pow(d,2)*pow(k1,9) - 27*pow(k1,11)))))/4620.;

return res;
}

double CDDDH_Integral_xxz_22_case_4( double k1, double k2, double a2, double b2, double c2, double d2, double d ) 
{   double res =0.;

return res;
}

// 113 - xxz XX Done
//
//
// 113 - xxz YY

double CDDDH_Integral_xxz_33_case_1( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{   double res =(pow(d,-6)*(-1980*pow(c2,2)*pow(k1,7) - 3465*b2*c2*pow(k1,8) - 3080*a2*c2*pow(k1,9) - 1540*pow(b2,2)*pow(k1,9) - 
2772*a2*b2*pow(k1,10) - 1260*pow(a2,2)*pow(k1,11) - 2772*pow(d2,2)*(pow(k1,5) - pow(k2,5)) + 
1980*pow(c2,2)*pow(k2,7) - 165*d2*(28*c2*(pow(k1,6) - pow(k2,6)) + 24*b2*(pow(k1,7) - pow(k2,7)) + 
21*a2*(pow(k1,8) - pow(k2,8))) + 3465*b2*c2*pow(k2,8) + 
11*pow(d,2)*(420*pow(d2,2)*(pow(k1,3) - pow(k2,3)) + 252*pow(c2,2)*(pow(k1,5) - pow(k2,5)) + 
42*d2*(15*c2*(pow(k1,4) - pow(k2,4)) + 12*b2*(pow(k1,5) - pow(k2,5)) + 10*a2*(pow(k1,6) - pow(k2,6))) + 
60*c2*(7*b2*(pow(k1,6) - pow(k2,6)) + 6*a2*(pow(k1,7) - pow(k2,7))) + 
5*(36*pow(b2,2)*(pow(k1,7) - pow(k2,7)) + 63*a2*b2*(pow(k1,8) - pow(k2,8)) + 28*pow(a2,2)*(pow(k1,9) - pow(k2,9)))
) + 3080*a2*c2*pow(k2,9) + 1540*pow(b2,2)*pow(k2,9) + 2772*a2*b2*pow(k2,10) + 1260*pow(a2,2)*pow(k2,11)))/4620.;

return res;
}

double CDDDH_Integral_xxz_33_case_2( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{   double res =0.;

return res;
}

double CDDDH_Integral_xxz_33_case_3( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{   double res =-(pow(d,-6)*(396*pow(c2,2)*(2*pow(d,7) - 7*pow(d,2)*pow(k1,5) + 5*pow(k1,7)) + 
220*pow(b2,2)*(2*pow(d,9) - 9*pow(d,2)*pow(k1,7) + 7*pow(k1,9)) + 
55*c2*(42*d2*(pow(d,6) - 3*pow(d,2)*pow(k1,4) + 2*pow(k1,6)) + 
21*b2*(pow(d,8) - 4*pow(d,2)*pow(k1,6) + 3*pow(k1,8)) + 8*a2*(2*pow(d,9) - 9*pow(d,2)*pow(k1,7) + 7*pow(k1,9)))\
+ 99*b2*(8*d2*(2*pow(d,7) - 7*pow(d,2)*pow(k1,5) + 5*pow(k1,7)) + 
7*a2*(pow(d,10) - 5*pow(d,2)*pow(k1,8) + 4*pow(k1,10))) + 
7*(132*pow(d2,2)*(2*pow(d,5) - 5*pow(d,2)*pow(k1,3) + 3*pow(k1,5)) + 
165*a2*d2*(pow(d,8) - 4*pow(d,2)*pow(k1,6) + 3*pow(k1,8)) + 
20*pow(a2,2)*(2*pow(d,11) - 11*pow(d,2)*pow(k1,9) + 9*pow(k1,11)))))/4620.;

return res;
}

double CDDDH_Integral_xxz_33_case_4( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{   double res =0.;

return res;
}

// 113 - xxz YY Done
//
//
// 113 - xxz ZZ

double CDDDH_Integral_xxz_44_case_1( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{   double res =(pow(d,-6)*(11*pow(d,2)*(420*pow(d2,2)*(pow(k1,3) - pow(k2,3)) + 252*pow(c2,2)*(pow(k1,5) - pow(k2,5)) + 
42*d2*(15*c2*(pow(k1,4) - pow(k2,4)) + 12*b2*(pow(k1,5) - pow(k2,5)) + 10*a2*(pow(k1,6) - pow(k2,6))) + 
60*c2*(7*b2*(pow(k1,6) - pow(k2,6)) + 6*a2*(pow(k1,7) - pow(k2,7))) + 
5*(36*pow(b2,2)*(pow(k1,7) - pow(k2,7)) + 63*a2*b2*(pow(k1,8) - pow(k2,8)) + 28*pow(a2,2)*(pow(k1,9) - pow(k2,9)))
) + 4*(2772*pow(d2,2)*(pow(k1,5) - pow(k2,5)) + 1980*pow(c2,2)*(pow(k1,7) - pow(k2,7)) + 
165*d2*(28*c2*(pow(k1,6) - pow(k2,6)) + 24*b2*(pow(k1,7) - pow(k2,7)) + 21*a2*(pow(k1,8) - pow(k2,8))) + 
385*c2*(9*b2*(pow(k1,8) - pow(k2,8)) + 8*a2*(pow(k1,9) - pow(k2,9))) + 
28*(55*pow(b2,2)*(pow(k1,9) - pow(k2,9)) + 99*a2*b2*(pow(k1,10) - pow(k2,10)) + 
45*pow(a2,2)*(pow(k1,11) - pow(k2,11))))))/4620.;

return res;
}

double CDDDH_Integral_xxz_44_case_2( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{   double res =0.;

return res;
}

double CDDDH_Integral_xxz_44_case_3( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{   double res =-(pow(d,-6)*(396*pow(c2,2)*(27*pow(d,7) - 7*pow(d,2)*pow(k1,5) - 20*pow(k1,7)) + 
110*c2*(21*d2*(11*pow(d,6) - 3*pow(d,2)*pow(k1,4) - 8*pow(k1,6)) + 
42*b2*(4*pow(d,8) - pow(d,2)*pow(k1,6) - 3*pow(k1,8)) + 4*a2*(37*pow(d,9) - 9*pow(d,2)*pow(k1,7) - 28*pow(k1,9)))
+ 220*pow(b2,2)*(37*pow(d,9) - 9*pow(d,2)*pow(k1,7) - 28*pow(k1,9)) + 
99*b2*(8*d2*(27*pow(d,7) - 7*pow(d,2)*pow(k1,5) - 20*pow(k1,7)) + 
7*a2*(21*pow(d,10) - 5*pow(d,2)*pow(k1,8) - 16*pow(k1,10))) + 
28*(33*pow(d2,2)*(17*pow(d,5) - 5*pow(d,2)*pow(k1,3) - 12*pow(k1,5)) + 
165*a2*d2*(4*pow(d,8) - pow(d,2)*pow(k1,6) - 3*pow(k1,8)) + 
5*pow(a2,2)*(47*pow(d,11) - 11*pow(d,2)*pow(k1,9) - 36*pow(k1,11)))))/4620.;

return res;
}

double CDDDH_Integral_xxz_44_case_4( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{   double res =0.;

return res;
}

// 113 - xxz ZZ Done
//
//
// 123 - xyz XY

double CDDDH_Integral_xyz_23_case_1( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{   double res =3*pow(d,-6)*((pow(d2,2)*(-pow(k1,5) + pow(k2,5)))/5. + (c2*d2*(-pow(k1,6) + pow(k2,6)))/3. + 
(2*b2*d2*(-pow(k1,7) + pow(k2,7)))/7. + (pow(c2,2)*(-pow(k1,7) + pow(k2,7)))/7. + 
(b2*c2*(-pow(k1,8) + pow(k2,8)))/4. + (a2*d2*(-pow(k1,8) + pow(k2,8)))/4. + (2*a2*c2*(-pow(k1,9) + pow(k2,9)))/9. + 
(pow(b2,2)*(-pow(k1,9) + pow(k2,9)))/9. + (a2*b2*(-pow(k1,10) + pow(k2,10)))/5. + 
(pow(a2,2)*(-pow(k1,11) + pow(k2,11)))/11.);

return res;
}

double CDDDH_Integral_xyz_23_case_2( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{   double res =0.;

return res;
}

double CDDDH_Integral_xyz_23_case_3( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{   double res =(pow(d,-6)*(1980*pow(c2,2)*(pow(d,7) - pow(k1,7)) + 
385*c2*(12*d2*(pow(d,6) - pow(k1,6)) + 9*b2*(pow(d,8) - pow(k1,8)) + 8*a2*(pow(d,9) - pow(k1,9))) + 
1540*pow(b2,2)*(pow(d,9) - pow(k1,9)) + 396*b2*(10*d2*(pow(d,7) - pow(k1,7)) + 7*a2*(pow(d,10) - pow(k1,10))) + 
63*(44*pow(d2,2)*(pow(d,5) - pow(k1,5)) + 55*a2*d2*(pow(d,8) - pow(k1,8)) + 20*pow(a2,2)*(pow(d,11) - pow(k1,11)))))/
4620.;

return res;
}

double CDDDH_Integral_xyz_23_case_4( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{   double res =0.;

return res;
}

// 123 - xyz XY Done
//
//
// 133 - xzz SX

double CDDDH_Integral_xzz_12_case_1( double k1, double k2, double a1, double b1, double c1, double d1,
double a2, double b2, double c2, double d2, double d )
{   double res =-4*pow(3,0.5)*pow(d,-5)*(-(d1*d2*pow(k1,4))/4. - ((c2*d1 + c1*d2)*pow(k1,5))/5. - ((c1*c2 + b2*d1 + b1*d2)*pow(k1,6))/6. - 
((b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(k1,7))/7. - ((b1*b2 + a2*c1 + a1*c2)*pow(k1,8))/8. - 
((a2*b1 + a1*b2)*pow(k1,9))/9. - (a1*a2*pow(k1,10))/10. + (d1*d2*pow(k2,4))/4. + ((c2*d1 + c1*d2)*pow(k2,5))/5. + 
((c1*c2 + b2*d1 + b1*d2)*pow(k2,6))/6. + ((b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(k2,7))/7. + 
((b1*b2 + a2*c1 + a1*c2)*pow(k2,8))/8. + ((a2*b1 + a1*b2)*pow(k2,9))/9. + (a1*a2*pow(k2,10))/10.);

return res;
}

double CDDDH_Integral_xzz_12_case_2( double k1, double k2, double a1, double b1, double c1, double d1,
double a2, double b2, double c2, double d2, double d )
{   double res =0.;

return res;
}

double CDDDH_Integral_xzz_12_case_3( double k1, double k2, double a1, double b1, double c1, double d1,
double a2, double b2, double c2, double d2, double d )
{   double res =(pow(3,-0.5)*pow(d,-5)*(pow(d,4)*(-6*d2*(84*c1*d + 105*d1 + 70*b1*pow(d,2) + 60*a1*pow(d,3)) + 
d*(-12*(42*c2 + 5*d*(7*b2 + 6*a2*d))*d1 + 
d*(-420*c1*c2 - 360*(b2*c1 + b1*c2)*d - 315*(b1*b2 + a2*c1 + a1*c2)*pow(d,2) - 280*(a2*b1 + a1*b2)*pow(d,3) - 
252*a1*a2*pow(d,4)))) + 630*d1*d2*pow(k1,4) + 504*(c2*d1 + c1*d2)*pow(k1,5) + 
420*(c1*c2 + b2*d1 + b1*d2)*pow(k1,6) + 360*(b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(k1,7) + 
315*(b1*b2 + a2*c1 + a1*c2)*pow(k1,8) + 280*(a2*b1 + a1*b2)*pow(k1,9) + 252*a1*a2*pow(k1,10)))/210.;

return res;
}

double CDDDH_Integral_xzz_12_case_4( double k1, double k2, double a1, double b1, double c1, double d1,
double a2, double b2, double c2, double d2, double d )
{   double res =0.;

return res;
}

// 113 - xzz SX Done
//
//
// 133 - xzz XZ

double CDDDH_Integral_xzz_24_case_1( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{   double res =(pow(d,-6)*(2772*pow(d2,2)*(pow(k1,5) - pow(k2,5)) + 1980*pow(c2,2)*(pow(k1,7) - pow(k2,7)) + 
165*d2*(28*c2*(pow(k1,6) - pow(k2,6)) + 24*b2*(pow(k1,7) - pow(k2,7)) + 21*a2*(pow(k1,8) - pow(k2,8))) + 
385*c2*(9*b2*(pow(k1,8) - pow(k2,8)) + 8*a2*(pow(k1,9) - pow(k2,9))) + 
28*(55*pow(b2,2)*(pow(k1,9) - pow(k2,9)) + 99*a2*b2*(pow(k1,10) - pow(k2,10)) + 
45*pow(a2,2)*(pow(k1,11) - pow(k2,11)))))/1155.;

return res;
}

double CDDDH_Integral_xzz_24_case_2( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{   double res =0.;

return res;
}

double CDDDH_Integral_xzz_24_case_3( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{   double res =-(pow(d,-6)*(1980*pow(c2,2)*(pow(d,7) - pow(k1,7)) + 
385*c2*(12*d2*(pow(d,6) - pow(k1,6)) + 9*b2*(pow(d,8) - pow(k1,8)) + 8*a2*(pow(d,9) - pow(k1,9))) + 
1540*pow(b2,2)*(pow(d,9) - pow(k1,9)) + 396*b2*(10*d2*(pow(d,7) - pow(k1,7)) + 7*a2*(pow(d,10) - pow(k1,10))) + 
63*(44*pow(d2,2)*(pow(d,5) - pow(k1,5)) + 55*a2*d2*(pow(d,8) - pow(k1,8)) + 20*pow(a2,2)*(pow(d,11) - pow(k1,11)))))
/1155.;

return res;
}

double CDDDH_Integral_xzz_24_case_4( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{   double res =0.;

return res;
}

// 133 - xzz XZ Done
//
//
// 333 - zzz SS

double CDDDH_Integral_zzz_11_case_1( double k1, double k2, double a1, double b1, double c1, double d1, double d )
{   double res =-(pow(d,-4)*(420*pow(d1,2)*(pow(k1,3) - pow(k2,3)) + 252*pow(c1,2)*(pow(k1,5) - pow(k2,5)) + 
42*d1*(15*c1*(pow(k1,4) - pow(k2,4)) + 12*b1*(pow(k1,5) - pow(k2,5)) + 10*a1*(pow(k1,6) - pow(k2,6))) + 
60*c1*(7*b1*(pow(k1,6) - pow(k2,6)) + 6*a1*(pow(k1,7) - pow(k2,7))) + 
5*(36*pow(b1,2)*(pow(k1,7) - pow(k2,7)) + 63*a1*b1*(pow(k1,8) - pow(k2,8)) + 28*pow(a1,2)*(pow(k1,9) - pow(k2,9)))))
/210.;

return res;
}

double CDDDH_Integral_zzz_11_case_2( double k1, double k2, double a1, double b1, double c1, double d1, double d )
{   double res =0.;

return res;
}

double CDDDH_Integral_zzz_11_case_3( double k1, double k2, double a1, double b1, double c1, double d1, double d )
{   double res =(pow(d,-4)*(252*pow(c1,2)*(pow(d,5) - pow(k1,5)) + 
30*c1*(21*d1*(pow(d,4) - pow(k1,4)) + 14*b1*(pow(d,6) - pow(k1,6)) + 12*a1*(pow(d,7) - pow(k1,7))) + 
180*pow(b1,2)*(pow(d,7) - pow(k1,7)) + 63*b1*(8*d1*(pow(d,5) - pow(k1,5)) + 5*a1*(pow(d,8) - pow(k1,8))) + 
140*(3*pow(d1,2)*(pow(d,3) - pow(k1,3)) + 3*a1*d1*(pow(d,6) - pow(k1,6)) + pow(a1,2)*(pow(d,9) - pow(k1,9)))))/210.;

return res;
}

double CDDDH_Integral_zzz_11_case_4( double k1, double k2, double a1, double b1, double c1, double d1, double d )
{   double res =0.;

return res;
}

// 333 - zzz SS Done
//
//
// 333 - zzz SZ

double CDDDH_Integral_zzz_14_case_1( double k1, double k2, double a1, double b1, double c1, double d1,
double a2, double b2, double c2, double d2, double d )
{   double res =8*pow(3,0.5)*pow(d,-5)*(-(d1*d2*pow(k1,4))/4. - ((c2*d1 + c1*d2)*pow(k1,5))/5. - ((c1*c2 + b2*d1 + b1*d2)*pow(k1,6))/6. - 
((b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(k1,7))/7. - ((b1*b2 + a2*c1 + a1*c2)*pow(k1,8))/8. - 
((a2*b1 + a1*b2)*pow(k1,9))/9. - (a1*a2*pow(k1,10))/10. + (d1*d2*pow(k2,4))/4. + ((c2*d1 + c1*d2)*pow(k2,5))/5. + 
((c1*c2 + b2*d1 + b1*d2)*pow(k2,6))/6. + ((b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(k2,7))/7. + 
((b1*b2 + a2*c1 + a1*c2)*pow(k2,8))/8. + ((a2*b1 + a1*b2)*pow(k2,9))/9. + (a1*a2*pow(k2,10))/10.);

return res;
}

double CDDDH_Integral_zzz_14_case_2( double k1, double k2, double a1, double b1, double c1, double d1,
double a2, double b2, double c2, double d2, double d )
{   double res =0.;

return res;
}

double CDDDH_Integral_zzz_14_case_3( double k1, double k2, double a1, double b1, double c1, double d1,
double a2, double b2, double c2, double d2, double d )
{   double res =(pow(3,-0.5)*pow(d,-5)*(pow(d,4)*(6*d2*(84*c1*d + 105*d1 + 70*b1*pow(d,2) + 60*a1*pow(d,3)) + 
d*(420*c1*c2*d + 504*c2*d1 + 60*d*(7*b2 + 6*a2*d)*d1 + 360*(b2*c1 + b1*c2)*pow(d,2) + 
315*(b1*b2 + a2*c1 + a1*c2)*pow(d,3) + 280*(a2*b1 + a1*b2)*pow(d,4) + 252*a1*a2*pow(d,5))) - 
630*d1*d2*pow(k1,4) - 504*(c2*d1 + c1*d2)*pow(k1,5) - 420*(c1*c2 + b2*d1 + b1*d2)*pow(k1,6) - 
360*(b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(k1,7) - 315*(b1*b2 + a2*c1 + a1*c2)*pow(k1,8) - 
280*(a2*b1 + a1*b2)*pow(k1,9) - 252*a1*a2*pow(k1,10)))/105.;

return res;
}

double CDDDH_Integral_zzz_14_case_4( double k1, double k2, double a1, double b1, double c1, double d1,
double a2, double b2, double c2, double d2, double d )
{   double res =0.;

return res;
}

// 333 - zzz SZ Done
//
//
// 333 - zzz XXYY

double CDDDH_Integral_zzz_2233_case_1( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{   double res =(pow(d,-6)*(3960*pow(c2,2)*pow(k1,7) + 6930*b2*c2*pow(k1,8) + 6160*a2*c2*pow(k1,9) + 3080*pow(b2,2)*pow(k1,9) + 
5544*a2*b2*pow(k1,10) + 2520*pow(a2,2)*pow(k1,11) + 5544*pow(d2,2)*(pow(k1,5) - pow(k2,5)) - 
3960*pow(c2,2)*pow(k2,7) + 330*d2*(28*c2*(pow(k1,6) - pow(k2,6)) + 24*b2*(pow(k1,7) - pow(k2,7)) + 
21*a2*(pow(k1,8) - pow(k2,8))) - 6930*b2*c2*pow(k2,8) - 
11*pow(d,2)*(420*pow(d2,2)*(pow(k1,3) - pow(k2,3)) + 252*pow(c2,2)*(pow(k1,5) - pow(k2,5)) + 
42*d2*(15*c2*(pow(k1,4) - pow(k2,4)) + 12*b2*(pow(k1,5) - pow(k2,5)) + 10*a2*(pow(k1,6) - pow(k2,6))) + 
60*c2*(7*b2*(pow(k1,6) - pow(k2,6)) + 6*a2*(pow(k1,7) - pow(k2,7))) + 
5*(36*pow(b2,2)*(pow(k1,7) - pow(k2,7)) + 63*a2*b2*(pow(k1,8) - pow(k2,8)) + 28*pow(a2,2)*(pow(k1,9) - pow(k2,9)))
) - 6160*a2*c2*pow(k2,9) - 3080*pow(b2,2)*pow(k2,9) - 5544*a2*b2*pow(k2,10) - 2520*pow(a2,2)*pow(k2,11)))/2310.;

return res;
}

double CDDDH_Integral_zzz_2233_case_2( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{   double res =0.;

return res;
}

double CDDDH_Integral_zzz_2233_case_3( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{   double res =-(pow(d,-6)*(396*pow(c2,2)*(3*pow(d,7) + 7*pow(d,2)*pow(k1,5) - 10*pow(k1,7)) + 
110*c2*(21*d2*(pow(d,6) + 3*pow(d,2)*pow(k1,4) - 4*pow(k1,6)) + 
21*b2*(pow(d,8) + 2*pow(d,2)*pow(k1,6) - 3*pow(k1,8)) + 4*a2*(5*pow(d,9) + 9*pow(d,2)*pow(k1,7) - 14*pow(k1,9)))\
+ 220*pow(b2,2)*(5*pow(d,9) + 9*pow(d,2)*pow(k1,7) - 14*pow(k1,9)) + 
99*b2*(8*d2*(3*pow(d,7) + 7*pow(d,2)*pow(k1,5) - 10*pow(k1,7)) + 
7*a2*(3*pow(d,10) + 5*pow(d,2)*pow(k1,8) - 8*pow(k1,10))) + 
14*(66*pow(d2,2)*(pow(d,5) + 5*pow(d,2)*pow(k1,3) - 6*pow(k1,5)) + 
165*a2*d2*(pow(d,8) + 2*pow(d,2)*pow(k1,6) - 3*pow(k1,8)) + 
10*pow(a2,2)*(7*pow(d,11) + 11*pow(d,2)*pow(k1,9) - 18*pow(k1,11)))))/2310.;

return res;
}

double CDDDH_Integral_zzz_2233_case_4( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{   double res =0.;

return res;
}

// 333 - zzz XXYY Done
//
//
// 333 - zzz ZZ

double CDDDH_Integral_zzz_44_case_1( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{   double res =-(pow(d,-6)*(11*pow(d,2)*(420*pow(d2,2)*(pow(k1,3) - pow(k2,3)) + 252*pow(c2,2)*(pow(k1,5) - pow(k2,5)) + 
42*d2*(15*c2*(pow(k1,4) - pow(k2,4)) + 12*b2*(pow(k1,5) - pow(k2,5)) + 10*a2*(pow(k1,6) - pow(k2,6))) + 
60*c2*(7*b2*(pow(k1,6) - pow(k2,6)) + 6*a2*(pow(k1,7) - pow(k2,7))) + 
5*(36*pow(b2,2)*(pow(k1,7) - pow(k2,7)) + 63*a2*b2*(pow(k1,8) - pow(k2,8)) + 
28*pow(a2,2)*(pow(k1,9) - pow(k2,9)))) + 
4*(2772*pow(d2,2)*(pow(k1,5) - pow(k2,5)) + 1980*pow(c2,2)*(pow(k1,7) - pow(k2,7)) + 
165*d2*(28*c2*(pow(k1,6) - pow(k2,6)) + 24*b2*(pow(k1,7) - pow(k2,7)) + 21*a2*(pow(k1,8) - pow(k2,8))) + 
385*c2*(9*b2*(pow(k1,8) - pow(k2,8)) + 8*a2*(pow(k1,9) - pow(k2,9))) + 
28*(55*pow(b2,2)*(pow(k1,9) - pow(k2,9)) + 99*a2*b2*(pow(k1,10) - pow(k2,10)) + 
45*pow(a2,2)*(pow(k1,11) - pow(k2,11))))))/2310.;

return res;
}

double CDDDH_Integral_zzz_44_case_2( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{   double res =0.;

return res;
}

double CDDDH_Integral_zzz_44_case_3( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{   double res =(pow(d,-6)*(396*pow(c2,2)*(27*pow(d,7) - 7*pow(d,2)*pow(k1,5) - 20*pow(k1,7)) + 
110*c2*(21*d2*(11*pow(d,6) - 3*pow(d,2)*pow(k1,4) - 8*pow(k1,6)) + 
42*b2*(4*pow(d,8) - pow(d,2)*pow(k1,6) - 3*pow(k1,8)) + 4*a2*(37*pow(d,9) - 9*pow(d,2)*pow(k1,7) - 28*pow(k1,9)))\
+ 220*pow(b2,2)*(37*pow(d,9) - 9*pow(d,2)*pow(k1,7) - 28*pow(k1,9)) + 
99*b2*(8*d2*(27*pow(d,7) - 7*pow(d,2)*pow(k1,5) - 20*pow(k1,7)) + 
7*a2*(21*pow(d,10) - 5*pow(d,2)*pow(k1,8) - 16*pow(k1,10))) + 
28*(33*pow(d2,2)*(17*pow(d,5) - 5*pow(d,2)*pow(k1,3) - 12*pow(k1,5)) + 
165*a2*d2*(4*pow(d,8) - pow(d,2)*pow(k1,6) - 3*pow(k1,8)) + 
5*pow(a2,2)*(47*pow(d,11) - 11*pow(d,2)*pow(k1,9) - 36*pow(k1,11)))))/2310.;

return res;
}

double CDDDH_Integral_zzz_44_case_4( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{   double res =0.;

return res;
}

// 333 - zzz ZZ Done


/* SYMMETRY
*
* 111  SX, XZ                              <=> 222 SY, YZ
*
* 112 - 121 - 211  SY, YZ                  <=> 122 - 212 - 221 SX, XZ
*
* 113 - 131 - 311  SS, SZ,*XX,*YY, ZZ      <=> 223 - 232 - 322 SS, SZ,*YY,*XX, ZZ
*
* 123 - 132 - 213 - 231 - 312 - 321    XY
*
* 133 - 313 - 331  SX, XZ                  <=> 233 - 323 - 332 SY, YZ
*
* 333  SS, SZ, XXYY, ZZ
*
*
* Total 27 - elem_cnt of rank three tensor
*
*/
