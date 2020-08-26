#include<stdio.h>
#include<gsl/gsl_sf_expint.h>
#include<gsl/gsl_math.h>

/* < S | dxH | X >  ==  < S | dyH | Y > */
double CDH_Integral_x_12_case_1( double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d )
{
double res = -((-(d1*d2*pow(4,-1)*pow(k1,4)) - (c2*d1 + c1*d2)*pow(5,-1)*pow(k1,5) - (c1*c2 + b2*d1 + b1*d2)*pow(6,-1)*pow(k1,6) - 
(b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(7,-1)*pow(k1,7) - (b1*b2 + a2*c1 + a1*c2)*pow(8,-1)*pow(k1,8) - (a2*b1 + a1*b2)*pow(9,-1)*pow(k1,9) - 
a1*a2*pow(10,-1)*pow(k1,10) + d1*d2*pow(4,-1)*pow(k2,4) + (c2*d1 + c1*d2)*pow(5,-1)*pow(k2,5) + 
(c1*c2 + b2*d1 + b1*d2)*pow(6,-1)*pow(k2,6) + (b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(7,-1)*pow(k2,7) + 
(b1*b2 + a2*c1 + a1*c2)*pow(8,-1)*pow(k2,8) + (a2*b1 + a1*b2)*pow(9,-1)*pow(k2,9) + a1*a2*pow(10,-1)*pow(k2,10))*
pow(pow(3,0.5)*pow(d,3),-1));

return res;
}

double CDH_Integral_x_12_case_2_sub_1( double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d )
{
double res = (-420*d1*d2*k2 + d*(35*d2*(6*c1*d + 12*d1 + 4*b1*pow(d,2) + 3*a1*pow(d,3)) + 
d*(35*(6*c2 + d*(4*b2 + 3*a2*d))*d1 + d*(140*c1*c2 + 105*(b2*c1 + b1*c2)*d + 84*(b1*b2 + a2*c1 + a1*c2)*pow(d,2) + 
70*(a2*b1 + a1*b2)*pow(d,3) + 60*a1*a2*pow(d,4)))) - 210*(c2*d1 + c1*d2)*pow(k2,2) - 140*(c1*c2 + b2*d1 + b1*d2)*pow(k2,3) - 
105*(b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(k2,4) - 84*(b1*b2 + a2*c1 + a1*c2)*pow(k2,5) - 70*(a2*b1 + a1*b2)*pow(k2,6) - 60*a1*a2*pow(k2,7))*
pow(420*pow(3,0.5),-1);

return res;
}

double CDH_Integral_x_12_case_2_sub_2( double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d )
{
double res = (pow(d,4)*(-6*d2*(84*c1*d + 105*d1 + 70*b1*pow(d,2) + 60*a1*pow(d,3)) + 
d*(-12*(42*c2 + 5*d*(7*b2 + 6*a2*d))*d1 + d*(-420*c1*c2 - 360*(b2*c1 + b1*c2)*d - 315*(b1*b2 + a2*c1 + a1*c2)*pow(d,2) - 
280*(a2*b1 + a1*b2)*pow(d,3) - 252*a1*a2*pow(d,4)))) + 630*d1*d2*pow(k1,4) + 504*(c2*d1 + c1*d2)*pow(k1,5) + 
420*(c1*c2 + b2*d1 + b1*d2)*pow(k1,6) + 360*(b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(k1,7) + 315*(b1*b2 + a2*c1 + a1*c2)*pow(k1,8) + 
280*(a2*b1 + a1*b2)*pow(k1,9) + 252*a1*a2*pow(k1,10))*pow(2520*pow(3,0.5)*pow(d,3),-1);

return res;
}

double CDH_Integral_x_12_case_3( double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d )
{
double res = (d1*d2*k1 - d1*d2*k2 + (c2*d1 + c1*d2)*pow(2,-1)*pow(k1,2) + (c1*c2 + b2*d1 + b1*d2)*pow(3,-1)*pow(k1,3) + 
(b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(4,-1)*pow(k1,4) + (b1*b2 + a2*c1 + a1*c2)*pow(5,-1)*pow(k1,5) + (a2*b1 + a1*b2)*pow(6,-1)*pow(k1,6) + 
a1*a2*pow(7,-1)*pow(k1,7) - (c2*d1 + c1*d2)*pow(2,-1)*pow(k2,2) - (c1*c2 + b2*d1 + b1*d2)*pow(3,-1)*pow(k2,3) - 
(b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(4,-1)*pow(k2,4) - (b1*b2 + a2*c1 + a1*c2)*pow(5,-1)*pow(k2,5) - (a2*b1 + a1*b2)*pow(6,-1)*pow(k2,6) - 
a1*a2*pow(7,-1)*pow(k2,7))*pow(pow(3,0.5),-1);

return res;
}

/* < S | dxH | X >  END */




/* < X | dxH | Z >  ==  < Y | dyH | Z > */
double CDH_Integral_x_24_case_1( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{
double res = -3*(pow(5,-1)*pow(d2,2)*(-pow(k1,5) + pow(k2,5)) + c2*d2*pow(3,-1)*(-pow(k1,6) + pow(k2,6)) + 
2*b2*d2*pow(7,-1)*(-pow(k1,7) + pow(k2,7)) + pow(7,-1)*pow(c2,2)*(-pow(k1,7) + pow(k2,7)) + 
b2*c2*pow(4,-1)*(-pow(k1,8) + pow(k2,8)) + a2*d2*pow(4,-1)*(-pow(k1,8) + pow(k2,8)) + 
2*a2*c2*pow(9,-1)*(-pow(k1,9) + pow(k2,9)) + pow(9,-1)*pow(b2,2)*(-pow(k1,9) + pow(k2,9)) + 
a2*b2*pow(5,-1)*(-pow(k1,10) + pow(k2,10)) + pow(11,-1)*pow(a2,2)*(-pow(k1,11) + pow(k2,11)))*pow(5*pow(d,4),-1);

return res;
}
double CDH_Integral_x_24_case_2_sub_1( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{
double res = d*pow(100,-1)*(30*(d - k2)*(d + k2)*pow(c2,2) + 60*log(d*pow(k2,-1))*pow(d2,2) + 
10*a2*(pow(d,3) - pow(k2,3))*(4*d2 + a2*(pow(d,3) + pow(k2,3))) + 
10*c2*(12*d2*(d - k2) + 4*b2*(pow(d,3) - pow(k2,3)) + 3*a2*(pow(d,4) - pow(k2,4))) + 
15*pow(b2,2)*(pow(d,4) - pow(k2,4)) + 12*b2*(5*d2*(d - k2)*(d + k2) + 2*a2*(pow(d,5) - pow(k2,5))));

return res;
}
double CDH_Integral_x_24_case_2_sub_2( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{
double res = -((1980*pow(c2,2)*(pow(d,7) - pow(k1,7)) + 385*c2*
(12*d2*(pow(d,6) - pow(k1,6)) + 9*b2*(pow(d,8) - pow(k1,8)) + 8*a2*(pow(d,9) - pow(k1,9))) + 
1540*pow(b2,2)*(pow(d,9) - pow(k1,9)) + 396*b2*(10*d2*(pow(d,7) - pow(k1,7)) + 7*a2*(pow(d,10) - pow(k1,10))) + 
63*(44*pow(d2,2)*(pow(d,5) - pow(k1,5)) + 55*a2*d2*(pow(d,8) - pow(k1,8)) + 20*pow(a2,2)*(pow(d,11) - pow(k1,11))))*
pow(23100*pow(d,4),-1));

return res;
}
double CDH_Integral_x_24_case_3( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{
double res = d*pow(100,-1)*(30*(k1 - k2)*(k1 + k2)*pow(c2,2) + 60*log(k1*pow(k2,-1))*pow(d2,2) + 
10*a2*(pow(k1,3) - pow(k2,3))*(4*d2 + a2*(pow(k1,3) + pow(k2,3))) + 
10*c2*(12*d2*(k1 - k2) + 4*b2*(pow(k1,3) - pow(k2,3)) + 3*a2*(pow(k1,4) - pow(k2,4))) + 
15*pow(b2,2)*(pow(k1,4) - pow(k2,4)) + 12*b2*(5*d2*(k1 - k2)*(k1 + k2) + 2*a2*(pow(k1,5) - pow(k2,5))));

return res;
}





/* < S | dzH | S > */

double CDH_Integral_z_11_case_1( double k1, double k2, double a1, double b1, double c1, double d1, double d )
{
double res = (pow(3,-1)*pow(d1,2)*(-pow(k1,3) + pow(k2,3)) + c1*d1*pow(2,-1)*(-pow(k1,4) + pow(k2,4)) + 
2*b1*d1*pow(5,-1)*(-pow(k1,5) + pow(k2,5)) + pow(5,-1)*pow(c1,2)*(-pow(k1,5) + pow(k2,5)) + 
b1*c1*pow(3,-1)*(-pow(k1,6) + pow(k2,6)) + a1*d1*pow(3,-1)*(-pow(k1,6) + pow(k2,6)) + 
2*a1*c1*pow(7,-1)*(-pow(k1,7) + pow(k2,7)) + pow(7,-1)*pow(b1,2)*(-pow(k1,7) + pow(k2,7)) + 
a1*b1*pow(4,-1)*(-pow(k1,8) + pow(k2,8)) + pow(9,-1)*pow(a1,2)*(-pow(k1,9) + pow(k2,9)))*pow(pow(d,2),-1);

return res;
}

/* < S | dzH | S > */
double CDH_Integral_z_11_case_2_sub_1( double k1, double k2, double a1, double b1, double c1, double d1, double d )
{
double res = 0.;

return res;
}

double CDH_Integral_z_11_case_2_sub_2( double k1, double k2, double a1, double b1, double c1, double d1, double d )
{
double res = (252*pow(c1,2)*(pow(d,5) - pow(k1,5)) + 30*c1*(21*d1*(pow(d,4) - pow(k1,4)) + 14*b1*(pow(d,6) - pow(k1,6)) + 
12*a1*(pow(d,7) - pow(k1,7))) + 180*pow(b1,2)*(pow(d,7) - pow(k1,7)) + 
63*b1*(8*d1*(pow(d,5) - pow(k1,5)) + 5*a1*(pow(d,8) - pow(k1,8))) + 
140*(3*pow(d1,2)*(pow(d,3) - pow(k1,3)) + 3*a1*d1*(pow(d,6) - pow(k1,6)) + pow(a1,2)*(pow(d,9) - pow(k1,9))))*
pow(1260*pow(d,2),-1);

return res;
}

double CDH_Integral_z_11_case_3( double k1, double k2, double a1, double b1, double c1, double d1, double d )
{
double res = 0.;

return res;
}

/* < S | dzH | S > done */




/* < S | dzH | Z > */

double CDH_Integral_z_14_case_1( double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d )
{
double res = 2*(-(d1*d2*pow(4,-1)*pow(k1,4)) - (c2*d1 + c1*d2)*pow(5,-1)*pow(k1,5) - (c1*c2 + b2*d1 + b1*d2)*pow(6,-1)*pow(k1,6) - 
(b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(7,-1)*pow(k1,7) - (b1*b2 + a2*c1 + a1*c2)*pow(8,-1)*pow(k1,8) - 
(a2*b1 + a1*b2)*pow(9,-1)*pow(k1,9) - a1*a2*pow(10,-1)*pow(k1,10) + d1*d2*pow(4,-1)*pow(k2,4) + 
(c2*d1 + c1*d2)*pow(5,-1)*pow(k2,5) + (c1*c2 + b2*d1 + b1*d2)*pow(6,-1)*pow(k2,6) + 
(b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(7,-1)*pow(k2,7) + (b1*b2 + a2*c1 + a1*c2)*pow(8,-1)*pow(k2,8) + 
(a2*b1 + a1*b2)*pow(9,-1)*pow(k2,9) + a1*a2*pow(10,-1)*pow(k2,10))*pow(pow(3,0.5)*pow(d,3),-1);

return res;
}


double CDH_Integral_z_14_case_2_sub_1( double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d )
{
double res = (-420*d1*d2*k2 + d*(35*d2*(6*c1*d + 12*d1 + 4*b1*pow(d,2) + 3*a1*pow(d,3)) + 
d*(35*(6*c2 + d*(4*b2 + 3*a2*d))*d1 + d*(140*c1*c2 + 105*(b2*c1 + b1*c2)*d + 84*(b1*b2 + a2*c1 + a1*c2)*pow(d,2) + 
70*(a2*b1 + a1*b2)*pow(d,3) + 60*a1*a2*pow(d,4)))) - 210*(c2*d1 + c1*d2)*pow(k2,2) - 
140*(c1*c2 + b2*d1 + b1*d2)*pow(k2,3) - 105*(b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(k2,4) - 
84*(b1*b2 + a2*c1 + a1*c2)*pow(k2,5) - 70*(a2*b1 + a1*b2)*pow(k2,6) - 60*a1*a2*pow(k2,7))*pow(420*pow(3,0.5),-1);

return res;
}


double CDH_Integral_z_14_case_2_sub_2( double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d )
{
double res = (pow(d,4)*(6*d2*(84*c1*d + 105*d1 + 70*b1*pow(d,2) + 60*a1*pow(d,3)) + 
d*(420*c1*c2*d + 504*c2*d1 + 60*d*(7*b2 + 6*a2*d)*d1 + 360*(b2*c1 + b1*c2)*pow(d,2) + 
315*(b1*b2 + a2*c1 + a1*c2)*pow(d,3) + 280*(a2*b1 + a1*b2)*pow(d,4) + 252*a1*a2*pow(d,5))) - 
630*d1*d2*pow(k1,4) - 504*(c2*d1 + c1*d2)*pow(k1,5) - 420*(c1*c2 + b2*d1 + b1*d2)*pow(k1,6) - 
360*(b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(k1,7) - 315*(b1*b2 + a2*c1 + a1*c2)*pow(k1,8) - 
280*(a2*b1 + a1*b2)*pow(k1,9) - 252*a1*a2*pow(k1,10))*pow(1260*pow(3,0.5)*pow(d,3),-1);

return res;
}


double CDH_Integral_z_14_case_3( double k1, double k2, double a1, double b1, double c1, double d1, double a2, double b2, double c2, double d2, double d )
{
double res = (d1*d2*k1 - d1*d2*k2 + (c2*d1 + c1*d2)*pow(2,-1)*pow(k1,2) + (c1*c2 + b2*d1 + b1*d2)*pow(3,-1)*pow(k1,3) + 
(b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(4,-1)*pow(k1,4) + (b1*b2 + a2*c1 + a1*c2)*pow(5,-1)*pow(k1,5) + 
(a2*b1 + a1*b2)*pow(6,-1)*pow(k1,6) + a1*a2*pow(7,-1)*pow(k1,7) - (c2*d1 + c1*d2)*pow(2,-1)*pow(k2,2) - 
(c1*c2 + b2*d1 + b1*d2)*pow(3,-1)*pow(k2,3) - (b2*c1 + b1*c2 + a2*d1 + a1*d2)*pow(4,-1)*pow(k2,4) - 
(b1*b2 + a2*c1 + a1*c2)*pow(5,-1)*pow(k2,5) - (a2*b1 + a1*b2)*pow(6,-1)*pow(k2,6) - a1*a2*pow(7,-1)*pow(k2,7))*
pow(pow(3,0.5),-1);

return res;
}

/* < S | dzH | Z > done */




/* < X or Y | dzH | X or Y > */
double CDH_Integral_z_2233_case_1( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{
double res = (8316*pow(d2,2)*(pow(k1,5) - pow(k2,5)) + 5940*pow(c2,2)*(pow(k1,7) - pow(k2,7)) + 
495*d2*(28*c2*(pow(k1,6) - pow(k2,6)) + 24*b2*(pow(k1,7) - pow(k2,7)) + 21*a2*(pow(k1,8) - pow(k2,8))) - 
55*pow(d,2)*(420*pow(d2,2)*(pow(k1,3) - pow(k2,3)) + 252*pow(c2,2)*(pow(k1,5) - pow(k2,5)) + 
42*d2*(15*c2*(pow(k1,4) - pow(k2,4)) + 12*b2*(pow(k1,5) - pow(k2,5)) + 10*a2*(pow(k1,6) - pow(k2,6))) + 
60*c2*(7*b2*(pow(k1,6) - pow(k2,6)) + 6*a2*(pow(k1,7) - pow(k2,7))) + 
5*(36*pow(b2,2)*(pow(k1,7) - pow(k2,7)) + 63*a2*b2*(pow(k1,8) - pow(k2,8)) + 28*pow(a2,2)*(pow(k1,9) - pow(k2,9)))) + 
1155*c2*(9*b2*(pow(k1,8) - pow(k2,8)) + 8*a2*(pow(k1,9) - pow(k2,9))) + 
84*(55*pow(b2,2)*(pow(k1,9) - pow(k2,9)) + 99*a2*b2*(pow(k1,10) - pow(k2,10)) + 45*pow(a2,2)*(pow(k1,11) - pow(k2,11))))*pow(69300*pow(d,4),-1);

return res;
}

double CDH_Integral_z_2233_case_2_sub_1( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{
double res = 2*d*pow(5,-1)*(log(k2*pow(d,-1))*pow(d2,2) + b2*d2*(-pow(d,2) + pow(k2,2)) + pow(2,-1)*pow(c2,2)*(-pow(d,2) + pow(k2,2)) - 
a2*pow(6,-1)*(pow(d,3) - pow(k2,3))*(4*d2 + a2*(pow(d,3) + pow(k2,3))) + pow(4,-1)*pow(b2,2)*(-pow(d,4) + pow(k2,4)) + 
c2*pow(6,-1)*(12*d2*(-d + k2) + 4*b2*(-pow(d,3) + pow(k2,3)) + 3*a2*(-pow(d,4) + pow(k2,4))) + 
2*a2*b2*pow(5,-1)*(-pow(d,5) + pow(k2,5)));

return res;
}

double CDH_Integral_z_2233_case_2_sub_2( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{   
double res = (1980*pow(c2,2)*(4*pow(d,7) - 7*pow(d,2)*pow(k1,5) + 3*pow(k1,7)) + 660*pow(b2,2)*(8*pow(d,9) - 15*pow(d,2)*pow(k1,7) + 7*pow(k1,9)) + 
165*c2*(42*d2*(3*pow(d,6) - 5*pow(d,2)*pow(k1,4) + 2*pow(k1,6)) + 7*b2*(11*pow(d,8) - 20*pow(d,2)*pow(k1,6) + 9*pow(k1,8)) + 
8*a2*(8*pow(d,9) - 15*pow(d,2)*pow(k1,7) + 7*pow(k1,9))) + 99*b2*
(40*d2*(4*pow(d,7) - 7*pow(d,2)*pow(k1,5) + 3*pow(k1,7)) + 7*a2*(13*pow(d,10) - 25*pow(d,2)*pow(k1,8) + 12*pow(k1,10))) + 
7*(132*pow(d2,2)*(16*pow(d,5) - 25*pow(d,2)*pow(k1,3) + 9*pow(k1,5)) + 165*a2*d2*(11*pow(d,8) - 20*pow(d,2)*pow(k1,6) + 9*pow(k1,8)) + 
20*pow(a2,2)*(28*pow(d,11) - 55*pow(d,2)*pow(k1,9) + 27*pow(k1,11))))*pow(69300*pow(d,4),-1);

return res;
}   

double CDH_Integral_z_2233_case_3( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{
double res = 2*d*pow(5,-1)*(log(k2*pow(k1,-1))*pow(d2,2) + b2*d2*(-pow(k1,2) + pow(k2,2)) + 
pow(2,-1)*pow(c2,2)*(-pow(k1,2) + pow(k2,2)) - 
a2*pow(6,-1)*(pow(k1,3) - pow(k2,3))*(4*d2 + a2*(pow(k1,3) + pow(k2,3))) + 
pow(4,-1)*pow(b2,2)*(-pow(k1,4) + pow(k2,4)) + 
c2*pow(6,-1)*(12*d2*(-k1 + k2) + 4*b2*(-pow(k1,3) + pow(k2,3)) + 3*a2*(-pow(k1,4) + pow(k2,4))) + 
2*a2*b2*pow(5,-1)*(-pow(k1,5) + pow(k2,5)));

return res;
}

/* < X or Y | dzH | X or Y >  done */





/* < Z | dzH | Z >  */
double CDH_Integral_z_44_case_1( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{
double res = -((55*pow(d,2)*(420*pow(d2,2)*(pow(k1,3) - pow(k2,3)) + 252*pow(c2,2)*(pow(k1,5) - pow(k2,5)) + 
42*d2*(15*c2*(pow(k1,4) - pow(k2,4)) + 12*b2*(pow(k1,5) - pow(k2,5)) + 10*a2*(pow(k1,6) - pow(k2,6))) + 
60*c2*(7*b2*(pow(k1,6) - pow(k2,6)) + 6*a2*(pow(k1,7) - pow(k2,7))) + 
5*(36*pow(b2,2)*(pow(k1,7) - pow(k2,7)) + 63*a2*b2*(pow(k1,8) - pow(k2,8)) + 28*pow(a2,2)*(pow(k1,9) - pow(k2,9)))
) + 6*(2772*pow(d2,2)*(pow(k1,5) - pow(k2,5)) + 1980*pow(c2,2)*(pow(k1,7) - pow(k2,7)) + 
165*d2*(28*c2*(pow(k1,6) - pow(k2,6)) + 24*b2*(pow(k1,7) - pow(k2,7)) + 21*a2*(pow(k1,8) - pow(k2,8))) + 
385*c2*(9*b2*(pow(k1,8) - pow(k2,8)) + 8*a2*(pow(k1,9) - pow(k2,9))) + 
28*(55*pow(b2,2)*(pow(k1,9) - pow(k2,9)) + 99*a2*b2*(pow(k1,10) - pow(k2,10)) + 
45*pow(a2,2)*(pow(k1,11) - pow(k2,11)))))*pow(69300*pow(d,4),-1));

return res;
}

double CDH_Integral_z_44_case_2_sub_1( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{
double res = d*pow(75,-1)*(30*(d - k2)*(d + k2)*pow(c2,2) + 60*log(d*pow(k2,-1))*pow(d2,2) + 
10*a2*(pow(d,3) - pow(k2,3))*(4*d2 + a2*(pow(d,3) + pow(k2,3))) + 
10*c2*(12*d2*(d - k2) + 4*b2*(pow(d,3) - pow(k2,3)) + 3*a2*(pow(d,4) - pow(k2,4))) + 
15*pow(b2,2)*(pow(d,4) - pow(k2,4)) + 12*b2*(5*d2*(d - k2)*(d + k2) + 2*a2*(pow(d,5) - pow(k2,5))));

return res;
}

double CDH_Integral_z_44_case_2_sub_2( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{   
double res = (1980*pow(c2,2)*(13*pow(d,7) - 7*pow(d,2)*pow(k1,5) - 6*pow(k1,7)) + 
330*c2*(21*d2*(9*pow(d,6) - 5*pow(d,2)*pow(k1,4) - 4*pow(k1,6)) + 
7*b2*(19*pow(d,8) - 10*pow(d,2)*pow(k1,6) - 9*pow(k1,8)) + 4*a2*(29*pow(d,9) - 15*pow(d,2)*pow(k1,7) - 14*pow(k1,9))
) + 660*pow(b2,2)*(29*pow(d,9) - 15*pow(d,2)*pow(k1,7) - 14*pow(k1,9)) + 
99*b2*(40*d2*(13*pow(d,7) - 7*pow(d,2)*pow(k1,5) - 6*pow(k1,7)) + 
7*a2*(49*pow(d,10) - 25*pow(d,2)*pow(k1,8) - 24*pow(k1,10))) + 
14*(66*pow(d2,2)*(43*pow(d,5) - 25*pow(d,2)*pow(k1,3) - 18*pow(k1,5)) + 
165*a2*d2*(19*pow(d,8) - 10*pow(d,2)*pow(k1,6) - 9*pow(k1,8)) + 
10*pow(a2,2)*(109*pow(d,11) - 55*pow(d,2)*pow(k1,9) - 54*pow(k1,11))))*pow(69300*pow(d,4),-1);

return res;
}

double CDH_Integral_z_44_case_3( double k1, double k2, double a2, double b2, double c2, double d2, double d )
{
double res = d*pow(75,-1)*(30*(k1 - k2)*(k1 + k2)*pow(c2,2) + 60*log(k1*pow(k2,-1))*pow(d2,2) + 
10*a2*(pow(k1,3) - pow(k2,3))*(4*d2 + a2*(pow(k1,3) + pow(k2,3))) + 
10*c2*(12*d2*(k1 - k2) + 4*b2*(pow(k1,3) - pow(k2,3)) + 3*a2*(pow(k1,4) - pow(k2,4))) + 
15*pow(b2,2)*(pow(k1,4) - pow(k2,4)) + 12*b2*(5*d2*(k1 - k2)*(k1 + k2) + 2*a2*(pow(k1,5) - pow(k2,5))));

return res;
}




