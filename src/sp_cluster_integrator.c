/**
 * (c) Author:    Woongkyu Jee, woong.jee.16@ucl.ac.uk, wldndrb1@gmail.com
 * Created:   02.06.2019 ~
 * 	
 * University College London, Department of Chemistry
 **/
#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_math.h>

#include"sp_cluster_integrator.h"

//#define INT_DEBUG

int sp_cluster_integrator_lut_b_search( double dist, int knot_stride, const double* integral_knot )
{
    int le = 0;
    int re = knot_stride - 1;

    while(1)
    {   if( integral_knot[(le+re)/2] < dist )
            le = (le+re)/2;
        else // dist <= integral_knot[(le+re)/2]
            re = (le+re)/2;
        if( (le+1) == re )
            break;
    }

    return le;
}


// LONG-RANGE INTEGRAL **************************************************************************************

// Calculate SS Integral
double sp_cluster_integrator_get_ch_11_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   
    double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CH_Integral_11_case_1( sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp->charge_shell*ion->charge_core*CH_Integral_11_case_2_sub_1( dist, sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp->charge_shell*ion->charge_core*CH_Integral_11_case_2_sub_2( sp->knot[i], dist,
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CH_Integral_11_case_3( sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3], dist );
        }
    }
    return Return*HA_TO_EV_UNIT;
}

// Calculate SZ Integral
double sp_cluster_integrator_get_ch_14_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CH_Integral_14_case_1( sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp->charge_shell*ion->charge_core*CH_Integral_14_case_2_sub_1( dist, sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp->charge_shell*ion->charge_core*CH_Integral_14_case_2_sub_2( sp->knot[i], dist,
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CH_Integral_14_case_3( sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*HA_TO_EV_UNIT;
}

// Calculate XX or YY Integral
double sp_cluster_integrator_get_ch_2233_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CH_Integral_2233_case_1( sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp->charge_shell*ion->charge_core*CH_Integral_2233_case_2_sub_1( dist, sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp->charge_shell*ion->charge_core*CH_Integral_2233_case_2_sub_2( sp->knot[i], dist,
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CH_Integral_2233_case_3( sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*HA_TO_EV_UNIT;
}

// Calculate ZZ Integral
double sp_cluster_integrator_get_ch_44_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CH_Integral_44_case_1( sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp->charge_shell*ion->charge_core*CH_Integral_44_case_2_sub_1( dist, sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp->charge_shell*ion->charge_core*CH_Integral_44_case_2_sub_2( sp->knot[i], dist,
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CH_Integral_44_case_3( sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*HA_TO_EV_UNIT;
}






// SHORT-RANGE INTEGRAL **************************************************************************************

// Calculate SS Integral
double sp_cluster_integrator_get_sh_11_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp->short_range_a_s)+(ion->short_range_a))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp->short_range_a_p)+(ion->short_range_a))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp->short_range_r_s*ion->short_range_r,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp->short_range_r_p*ion->short_range_r,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;
    
    int knot, ion_type;   double d, c, b, a;
    if( sp_sys->integral_lut == SP_INTEGRAL_TRUE )  // if LUT exist
    {
        dist = dist*TO_BOHR_RADII;  // back to Ang Unit
        knot = sp_cluster_integrator_lut_b_search( dist , sp_sys->knot_stride, (const double*)sp_sys->integral_knot ); 
 
        if( ion->charge_core < 0 ) // if it is anion
            ion_type = 0;
        else
            ion_type = 1;
    
        /* NOTE THAT ...
         *
         * ion_type 0 -> always aion    ( e.g., O(-1) or F(-2) )
         *
         * ion_type 1 -> always cation  ( e.g., Ba(2+)         )
         *
         */


        d = sp_sys->integral_vs_cla_s_ss[ion_type][knot][0];
        c = sp_sys->integral_vs_cla_s_ss[ion_type][knot][1];
        b = sp_sys->integral_vs_cla_s_ss[ion_type][knot][2];
        a = sp_sys->integral_vs_cla_s_ss[ion_type][knot][3];

        Return = d*pow(dist,3.) + c*pow(dist,2.) + b*dist + a;

        #ifdef INT_DEBUG
        printf("knot: %d, %lf\t%lf\t%lf\t%lf\n",knot,d,c,b,a);
        printf("ss:%lf\n",Return);
        #endif
        return Return;
    }


    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSH_Integral_11_case_1( mix_a_s, mix_r_s, sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSH_Integral_11_case_2_sub_1( mix_a_s, mix_r_s, dist, sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSH_Integral_11_case_2_sub_2( mix_a_s, mix_r_s, sp->knot[i], dist,
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSH_Integral_11_case_3( mix_a_s, mix_r_s, sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3], dist );
        }
    }
    return Return*HA_TO_EV_UNIT;
}


// Calculate SZ Integral
double sp_cluster_integrator_get_sh_14_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp->short_range_a_s)+(ion->short_range_a))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp->short_range_a_p)+(ion->short_range_a))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp->short_range_r_s*ion->short_range_r,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp->short_range_r_p*ion->short_range_r,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;

    int knot, ion_type;   double d, c, b, a;
    if( sp_sys->integral_lut == SP_INTEGRAL_TRUE )  // if LUT exist
    {
        dist = dist*TO_BOHR_RADII;  // back to Ang Unit
        knot = sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot ); 
 
        if( ion->charge_core < 0 ) // if it is anion
            ion_type = 0;
        else
            ion_type = 1;
        
        d = sp_sys->integral_vs_cla_s_sz[ion_type][knot][0];
        c = sp_sys->integral_vs_cla_s_sz[ion_type][knot][1];
        b = sp_sys->integral_vs_cla_s_sz[ion_type][knot][2];
        a = sp_sys->integral_vs_cla_s_sz[ion_type][knot][3];

        Return = d*pow(dist,3.) + c*pow(dist,2.) + b*dist + a;

        return Return;
    }


    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSH_Integral_14_case_1( (mix_a_s+mix_a_p)/2., pow(mix_r_s*mix_r_p,0.5), sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSH_Integral_14_case_2_sub_1( (mix_a_s+mix_a_p)/2., pow(mix_r_s*mix_r_p,0.5), dist, sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSH_Integral_14_case_2_sub_2( (mix_a_s+mix_a_p)/2., pow(mix_r_s*mix_r_p,0.5), sp->knot[i], dist,
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSH_Integral_14_case_3( (mix_a_s+mix_a_p)/2., pow(mix_r_s*mix_r_p,0.5), sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*HA_TO_EV_UNIT;
}


// Calculate XX or YY Integral
double sp_cluster_integrator_get_sh_2233_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp->short_range_a_s)+(ion->short_range_a))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp->short_range_a_p)+(ion->short_range_a))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp->short_range_r_s*ion->short_range_r,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp->short_range_r_p*ion->short_range_r,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;

    int knot, ion_type;   double d, c, b, a;
    if( sp_sys->integral_lut == SP_INTEGRAL_TRUE )  // if LUT exist
    {
        dist = dist*TO_BOHR_RADII;  // back to Ang Unit
        knot = sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot ); 
 
        if( ion->charge_core < 0 ) // if it is anion
            ion_type = 0;
        else
            ion_type = 1;
        
        d = sp_sys->integral_vs_cla_s_xxyy[ion_type][knot][0];
        c = sp_sys->integral_vs_cla_s_xxyy[ion_type][knot][1];
        b = sp_sys->integral_vs_cla_s_xxyy[ion_type][knot][2];
        a = sp_sys->integral_vs_cla_s_xxyy[ion_type][knot][3];

        Return = d*pow(dist,3.) + c*pow(dist,2.) + b*dist + a;

        return Return;
    }

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSH_Integral_2233_case_1( mix_a_p, mix_r_p, sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSH_Integral_2233_case_2_sub_1( mix_a_p, mix_r_p, dist, sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSH_Integral_2233_case_2_sub_2( mix_a_p, mix_r_p, sp->knot[i], dist,
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSH_Integral_2233_case_3( mix_a_p, mix_r_p, sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*HA_TO_EV_UNIT;
}


// Calculate ZZ Integral
double sp_cluster_integrator_get_sh_44_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp->short_range_a_s)+(ion->short_range_a))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp->short_range_a_p)+(ion->short_range_a))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp->short_range_r_s*ion->short_range_r,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp->short_range_r_p*ion->short_range_r,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;


    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;

    int knot, ion_type;   double d, c, b, a;
    if( sp_sys->integral_lut == SP_INTEGRAL_TRUE )  // if LUT exist
    {
        dist = dist*TO_BOHR_RADII;  // back to Ang Unit
        knot = sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot ); 
 
        if( ion->charge_core < 0 ) // if it is anion
            ion_type = 0;
        else
            ion_type = 1;
        
        d = sp_sys->integral_vs_cla_s_zz[ion_type][knot][0];
        c = sp_sys->integral_vs_cla_s_zz[ion_type][knot][1];
        b = sp_sys->integral_vs_cla_s_zz[ion_type][knot][2];
        a = sp_sys->integral_vs_cla_s_zz[ion_type][knot][3];

        Return = d*pow(dist,3.) + c*pow(dist,2.) + b*dist + a;

        return Return;
    }

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSH_Integral_44_case_1( mix_a_p, mix_r_p, sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSH_Integral_44_case_2_sub_1( mix_a_p, mix_r_p, dist, sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSH_Integral_44_case_2_sub_2( mix_a_p, mix_r_p, sp->knot[i], dist,
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSH_Integral_44_case_3( mix_a_p, mix_r_p, sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*HA_TO_EV_UNIT;
}










/* ************************************** first derivative methods ****************************************** */


// Calc CDX SX == CDY SY
double sp_cluster_integrator_get_ch_x_12_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   
    double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDH_Integral_x_12_case_1( sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp->charge_shell*ion->charge_core*CDH_Integral_x_12_case_2_sub_1( dist, sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp->charge_shell*ion->charge_core*CDH_Integral_x_12_case_2_sub_2( sp->knot[i], dist,
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDH_Integral_x_12_case_3( sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FHA_TO_FEV_UNIT;
}


// Calc CDX XZ == CDY YZ
double sp_cluster_integrator_get_ch_x_24_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   
    double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDH_Integral_x_24_case_1( sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp->charge_shell*ion->charge_core*CDH_Integral_x_24_case_2_sub_1( dist, sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp->charge_shell*ion->charge_core*CDH_Integral_x_24_case_2_sub_2( sp->knot[i], dist,
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDH_Integral_x_24_case_3( sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FHA_TO_FEV_UNIT;
}



// Calc CDZ SS
double sp_cluster_integrator_get_ch_z_11_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   
    double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDH_Integral_z_11_case_1( sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp->charge_shell*ion->charge_core*CDH_Integral_z_11_case_2_sub_1( dist, sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp->charge_shell*ion->charge_core*CDH_Integral_z_11_case_2_sub_2( sp->knot[i], dist,
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDH_Integral_z_11_case_3( sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3], dist );
        }
    }
    return Return*FHA_TO_FEV_UNIT;
}


// Calc CDZ SZ
double sp_cluster_integrator_get_ch_z_14_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   
    double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDH_Integral_z_14_case_1( sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp->charge_shell*ion->charge_core*CDH_Integral_z_14_case_2_sub_1( dist, sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp->charge_shell*ion->charge_core*CDH_Integral_z_14_case_2_sub_2( sp->knot[i], dist,
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDH_Integral_z_14_case_3( sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FHA_TO_FEV_UNIT;
}


// Calc CDZ XX YY
double sp_cluster_integrator_get_ch_z_2233_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   
    double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDH_Integral_z_2233_case_1( sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp->charge_shell*ion->charge_core*CDH_Integral_z_2233_case_2_sub_1( dist, sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp->charge_shell*ion->charge_core*CDH_Integral_z_2233_case_2_sub_2( sp->knot[i], dist,
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDH_Integral_z_2233_case_3( sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FHA_TO_FEV_UNIT;
}


// Calc CDZ ZZ
double sp_cluster_integrator_get_ch_z_44_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   
    double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDH_Integral_z_44_case_1( sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp->charge_shell*ion->charge_core*CDH_Integral_z_44_case_2_sub_1( dist, sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp->charge_shell*ion->charge_core*CDH_Integral_z_44_case_2_sub_2( sp->knot[i], dist,
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDH_Integral_z_44_case_3( sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FHA_TO_FEV_UNIT;
}


// Short Range Derivatives  !!!!!!!!


// Calc SDX SX == SDY SY
double sp_cluster_integrator_get_sh_x_12_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp->short_range_a_s)+(ion->short_range_a))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp->short_range_a_p)+(ion->short_range_a))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp->short_range_r_s*ion->short_range_r,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp->short_range_r_p*ion->short_range_r,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;

    int knot, ion_type;   double d, c, b, a;
    if( sp_sys->integral_lut == SP_INTEGRAL_TRUE )  // if LUT exist
    {
        dist = dist*TO_BOHR_RADII;  // back to Ang Unit
        knot = sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot ); 
 
        if( ion->charge_core < 0 ) // if it is anion
            ion_type = 0;
        else
            ion_type = 1;
        
        d = sp_sys->integral_vs_cla_s_x_sx[ion_type][knot][0];
        c = sp_sys->integral_vs_cla_s_x_sx[ion_type][knot][1];
        b = sp_sys->integral_vs_cla_s_x_sx[ion_type][knot][2];
        a = sp_sys->integral_vs_cla_s_x_sx[ion_type][knot][3];

        Return = d*pow(dist,3.) + c*pow(dist,2.) + b*dist + a;

        return Return;
    }

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDH_Integral_x_12_case_1( (mix_a_s+mix_a_p)/2., pow(mix_r_s*mix_r_p,0.5), sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDH_Integral_x_12_case_2_sub_1( (mix_a_s+mix_a_p)/2., pow(mix_r_s*mix_r_p,0.5), dist, sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDH_Integral_x_12_case_2_sub_2( (mix_a_s+mix_a_p)/2., pow(mix_r_s*mix_r_p,0.5), sp->knot[i], dist,
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDH_Integral_x_12_case_3( (mix_a_s+mix_a_p)/2., pow(mix_r_s*mix_r_p,0.5), sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FHA_TO_FEV_UNIT;
}



// Calc SDX XZ == SDY YZ
double sp_cluster_integrator_get_sh_x_24_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp->short_range_a_s)+(ion->short_range_a))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp->short_range_a_p)+(ion->short_range_a))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp->short_range_r_s*ion->short_range_r,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp->short_range_r_p*ion->short_range_r,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;

    int knot, ion_type;   double d, c, b, a;
    if( sp_sys->integral_lut == SP_INTEGRAL_TRUE )  // if LUT exist
    {
        dist = dist*TO_BOHR_RADII;  // back to Ang Unit
        knot = sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot ); 
 
        if( ion->charge_core < 0 ) // if it is anion
            ion_type = 0;
        else
            ion_type = 1;
        
        d = sp_sys->integral_vs_cla_s_x_xz[ion_type][knot][0];
        c = sp_sys->integral_vs_cla_s_x_xz[ion_type][knot][1];
        b = sp_sys->integral_vs_cla_s_x_xz[ion_type][knot][2];
        a = sp_sys->integral_vs_cla_s_x_xz[ion_type][knot][3];

        Return = d*pow(dist,3.) + c*pow(dist,2.) + b*dist + a;

        return Return;
    }

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDH_Integral_x_24_case_1( mix_a_p, mix_r_p, sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDH_Integral_x_24_case_2_sub_1( mix_a_p, mix_r_p, dist, sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDH_Integral_x_24_case_2_sub_2( mix_a_p, mix_r_p, sp->knot[i], dist,
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDH_Integral_x_24_case_3( mix_a_p, mix_r_p, sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FHA_TO_FEV_UNIT;
}


// Calc SDZ SS
double sp_cluster_integrator_get_sh_z_11_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp->short_range_a_s)+(ion->short_range_a))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp->short_range_a_p)+(ion->short_range_a))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp->short_range_r_s*ion->short_range_r,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp->short_range_r_p*ion->short_range_r,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;

    int knot, ion_type;   double d, c, b, a;
    if( sp_sys->integral_lut == SP_INTEGRAL_TRUE )  // if LUT exist
    {
        dist = dist*TO_BOHR_RADII;  // back to Ang Unit
        knot = sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot ); 
 
        if( ion->charge_core < 0 ) // if it is anion
            ion_type = 0;
        else
            ion_type = 1;
        
        d = sp_sys->integral_vs_cla_s_z_ss[ion_type][knot][0];
        c = sp_sys->integral_vs_cla_s_z_ss[ion_type][knot][1];
        b = sp_sys->integral_vs_cla_s_z_ss[ion_type][knot][2];
        a = sp_sys->integral_vs_cla_s_z_ss[ion_type][knot][3];

        Return = d*pow(dist,3.) + c*pow(dist,2.) + b*dist + a;

        return Return;
    }

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDH_Integral_z_11_case_1( mix_a_s, mix_r_s, sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDH_Integral_z_11_case_2_sub_1( mix_a_s, mix_r_s, dist, sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDH_Integral_z_11_case_2_sub_2( mix_a_s, mix_r_s, sp->knot[i], dist,
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDH_Integral_z_11_case_3( mix_a_s, mix_r_s, sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3], dist );
        }
    }
    return Return*FHA_TO_FEV_UNIT;
}



// Calc SDZ SZ
double sp_cluster_integrator_get_sh_z_14_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp->short_range_a_s)+(ion->short_range_a))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp->short_range_a_p)+(ion->short_range_a))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp->short_range_r_s*ion->short_range_r,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp->short_range_r_p*ion->short_range_r,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;

    int knot, ion_type;   double d, c, b, a;
    if( sp_sys->integral_lut == SP_INTEGRAL_TRUE )  // if LUT exist
    {
        dist = dist*TO_BOHR_RADII;  // back to Ang Unit
        knot = sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot ); 
 
        if( ion->charge_core < 0 ) // if it is anion
            ion_type = 0;
        else
            ion_type = 1;
        
        d = sp_sys->integral_vs_cla_s_z_sz[ion_type][knot][0];
        c = sp_sys->integral_vs_cla_s_z_sz[ion_type][knot][1];
        b = sp_sys->integral_vs_cla_s_z_sz[ion_type][knot][2];
        a = sp_sys->integral_vs_cla_s_z_sz[ion_type][knot][3];

        Return = d*pow(dist,3.) + c*pow(dist,2.) + b*dist + a;

        return Return;
    }

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDH_Integral_z_14_case_1( (mix_a_s+mix_a_p)/2., pow(mix_r_s*mix_r_p,0.5), sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDH_Integral_z_14_case_2_sub_1( (mix_a_s+mix_a_p)/2., pow(mix_r_s*mix_r_p,0.5), dist, sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDH_Integral_z_14_case_2_sub_2( (mix_a_s+mix_a_p)/2., pow(mix_r_s*mix_r_p,0.5), sp->knot[i], dist,
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDH_Integral_z_14_case_3( (mix_a_s+mix_a_p)/2., pow(mix_r_s*mix_r_p,0.5), sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FHA_TO_FEV_UNIT;
}


// Calc SDZ XXYY
double sp_cluster_integrator_get_sh_z_2233_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp->short_range_a_s)+(ion->short_range_a))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp->short_range_a_p)+(ion->short_range_a))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp->short_range_r_s*ion->short_range_r,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp->short_range_r_p*ion->short_range_r,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;

    int knot, ion_type;   double d, c, b, a;
    if( sp_sys->integral_lut == SP_INTEGRAL_TRUE )  // if LUT exist
    {
        dist = dist*TO_BOHR_RADII;  // back to Ang Unit
        knot = sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot ); 
 
        if( ion->charge_core < 0 ) // if it is anion
            ion_type = 0;
        else
            ion_type = 1;
        
        d = sp_sys->integral_vs_cla_s_z_xxyy[ion_type][knot][0];
        c = sp_sys->integral_vs_cla_s_z_xxyy[ion_type][knot][1];
        b = sp_sys->integral_vs_cla_s_z_xxyy[ion_type][knot][2];
        a = sp_sys->integral_vs_cla_s_z_xxyy[ion_type][knot][3];

        Return = d*pow(dist,3.) + c*pow(dist,2.) + b*dist + a;

        return Return;
    }

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDH_Integral_z_2233_case_1( mix_a_p, mix_r_p, sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDH_Integral_z_2233_case_2_sub_1( mix_a_p, mix_r_p, dist, sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDH_Integral_z_2233_case_2_sub_2( mix_a_p, mix_r_p, sp->knot[i], dist,
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDH_Integral_z_2233_case_3( mix_a_p, mix_r_p, sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FHA_TO_FEV_UNIT;
}


// Calc SDZ ZZ
double sp_cluster_integrator_get_sh_z_44_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp->short_range_a_s)+(ion->short_range_a))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp->short_range_a_p)+(ion->short_range_a))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp->short_range_r_s*ion->short_range_r,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp->short_range_r_p*ion->short_range_r,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;

    int knot, ion_type;   double d, c, b, a;
    if( sp_sys->integral_lut == SP_INTEGRAL_TRUE )  // if LUT exist
    {
        dist = dist*TO_BOHR_RADII;  // back to Ang Unit
        knot = sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot ); 
 
        if( ion->charge_core < 0 ) // if it is anion
            ion_type = 0;
        else
            ion_type = 1;
        
        d = sp_sys->integral_vs_cla_s_z_zz[ion_type][knot][0];
        c = sp_sys->integral_vs_cla_s_z_zz[ion_type][knot][1];
        b = sp_sys->integral_vs_cla_s_z_zz[ion_type][knot][2];
        a = sp_sys->integral_vs_cla_s_z_zz[ion_type][knot][3];

        Return = d*pow(dist,3.) + c*pow(dist,2.) + b*dist + a;

        return Return;
    }

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDH_Integral_z_44_case_1( mix_a_p, mix_r_p, sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDH_Integral_z_44_case_2_sub_1( mix_a_p, mix_r_p, dist, sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDH_Integral_z_44_case_2_sub_2( mix_a_p, mix_r_p, sp->knot[i], dist,
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDH_Integral_z_44_case_3( mix_a_p, mix_r_p, sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FHA_TO_FEV_UNIT;
}


 ///     ///     ///      ///     ///     ///      ///     ///     ///      ///     ///     ///      ///     ///     ///      ///     ///     ///      ///     ///     ///      ///     ///     ///      ///     ///     ///     
 ///     ///     ///      ///     ///     ///      ///     ///     ///      ///     ///     ///      ///     ///     ///      ///     ///     ///      ///     ///     ///      ///     ///     ///      ///     ///     ///     
 ///     ///     ///      ///     ///     ///      ///     ///     ///      ///     ///     ///      ///     ///     ///      ///     ///     ///      ///     ///     ///      ///     ///     ///      ///     ///     ///     
 ///     ///     ///      ///     ///     ///      ///     ///     ///      ///     ///     ///      ///     ///     ///      ///     ///     ///      ///     ///     ///      ///     ///     ///      ///     ///     ///     
 ///     ///     ///      ///     ///     ///      ///     ///     ///      ///     ///     ///      ///     ///     ///      ///     ///     ///      ///     ///     ///      ///     ///     ///      ///     ///     ///     








/* ************************************* 2nd Derivative Methods ************************************* */

// SECOND_DERIVATIVE_COULOMB_LONG_RANGE

// CDDH_XX

// SS   ==  CDDH_YY_SS 
double sp_cluster_integrator_get_ch_xx_11_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   
    double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_xx_11_case_1( sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_xx_11_case_2_sub_1( dist, sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_xx_11_case_2_sub_2( sp->knot[i], dist,
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_xx_11_case_3( sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}

// SZ   == CDDH_YY_SZ
double sp_cluster_integrator_get_ch_xx_14_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_xx_14_case_1( sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_xx_14_case_2_sub_1( dist, sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_xx_14_case_2_sub_2( sp->knot[i], dist,
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_xx_14_case_3( sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}




// XX   == CDDH_YY_YY
double sp_cluster_integrator_get_ch_xx_22_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_xx_22_case_1( sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_xx_22_case_2_sub_1( dist, sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_xx_22_case_2_sub_2( sp->knot[i], dist,
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_xx_22_case_3( sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}


// YY   == CDDH_YY_XX
double sp_cluster_integrator_get_ch_xx_33_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_xx_33_case_1( sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_xx_33_case_2_sub_1( dist, sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_xx_33_case_2_sub_2( sp->knot[i], dist,
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_xx_33_case_3( sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}


// ZZ   == CDDH_YY_ZZ
double sp_cluster_integrator_get_ch_xx_44_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_xx_44_case_1( sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_xx_44_case_2_sub_1( dist, sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_xx_44_case_2_sub_2( sp->knot[i], dist,
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_xx_44_case_3( sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}






// CDDH_XY

// XY
double sp_cluster_integrator_get_ch_xy_23_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_xy_23_case_1( sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_xy_23_case_2_sub_1( dist, sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_xy_23_case_2_sub_2( sp->knot[i], dist,
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_xy_23_case_3( sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}







// CDDH_XZ

// SX   == CDDH_YZ_SY
double sp_cluster_integrator_get_ch_xz_12_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_xz_12_case_1( sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_xz_12_case_2_sub_1( dist, sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_xz_12_case_2_sub_2( sp->knot[i], dist,
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_xz_12_case_3( sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}


// XZ   == CDDH_YZ_YZ
double sp_cluster_integrator_get_ch_xz_24_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_xz_24_case_1( sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_xz_24_case_2_sub_1( dist, sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_xz_24_case_2_sub_2( sp->knot[i], dist,
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_xz_24_case_3( sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}









// CDDH_ZZ

// SS 
double sp_cluster_integrator_get_ch_zz_11_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   
    double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_zz_11_case_1( sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_zz_11_case_2_sub_1( dist, sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_zz_11_case_2_sub_2( sp->knot[i], dist,
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_zz_11_case_3( sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}

// SZ
double sp_cluster_integrator_get_ch_zz_14_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_zz_14_case_1( sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_zz_14_case_2_sub_1( dist, sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_zz_14_case_2_sub_2( sp->knot[i], dist,
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_zz_14_case_3( sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}




// XXYY
double sp_cluster_integrator_get_ch_zz_2233_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_zz_2233_case_1( sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_zz_2233_case_2_sub_1( dist, sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_zz_2233_case_2_sub_2( sp->knot[i], dist,
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_zz_2233_case_3( sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}




// ZZ
double sp_cluster_integrator_get_ch_zz_44_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_zz_44_case_1( sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_zz_44_case_2_sub_1( dist, sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_zz_44_case_2_sub_2( sp->knot[i], dist,
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp->charge_shell*ion->charge_core*CDDH_Integral_zz_44_case_3( sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}



//// Second Derivative
//
//  Short-Range Repulsion



// SDDXX SS
double sp_cluster_integrator_get_sh_xx_11_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp->short_range_a_s)+(ion->short_range_a))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp->short_range_a_p)+(ion->short_range_a))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp->short_range_r_s*ion->short_range_r,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp->short_range_r_p*ion->short_range_r,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDDH_Integral_xx_11_case_1( mix_a_s, mix_r_s, sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDDH_Integral_xx_11_case_2( mix_a_s, mix_r_s, dist, sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDDH_Integral_xx_11_case_3( mix_a_s, mix_r_s, sp->knot[i], dist,
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDDH_Integral_xx_11_case_4( mix_a_s, mix_r_s, sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}



// SDDXX SZ
double sp_cluster_integrator_get_sh_xx_14_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp->short_range_a_s)+(ion->short_range_a))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp->short_range_a_p)+(ion->short_range_a))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp->short_range_r_s*ion->short_range_r,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp->short_range_r_p*ion->short_range_r,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDDH_Integral_xx_14_case_1( 0.5*(mix_a_s+mix_a_p), pow(mix_r_s*mix_r_p,0.5), sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDDH_Integral_xx_14_case_2( 0.5*(mix_a_s+mix_a_p), pow(mix_r_s*mix_r_p,0.5), dist, sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDDH_Integral_xx_14_case_3( 0.5*(mix_a_s+mix_a_p), pow(mix_r_s*mix_r_p,0.5), sp->knot[i], dist,
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDDH_Integral_xx_14_case_4( 0.5*(mix_a_s+mix_a_p), pow(mix_r_s*mix_r_p,0.5), sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}



// SDDXX XX == SDDYY YY
double sp_cluster_integrator_get_sh_xx_22_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp->short_range_a_s)+(ion->short_range_a))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp->short_range_a_p)+(ion->short_range_a))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp->short_range_r_s*ion->short_range_r,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp->short_range_r_p*ion->short_range_r,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDDH_Integral_xx_22_case_1( mix_a_p, mix_r_p, sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDDH_Integral_xx_22_case_2( mix_a_p, mix_r_p, dist, sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDDH_Integral_xx_22_case_3( mix_a_p, mix_r_p, sp->knot[i], dist,
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDDH_Integral_xx_22_case_4( mix_a_p, mix_r_p, sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}



// SDDXX YY == SDDYY XX
double sp_cluster_integrator_get_sh_xx_33_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp->short_range_a_s)+(ion->short_range_a))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp->short_range_a_p)+(ion->short_range_a))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp->short_range_r_s*ion->short_range_r,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp->short_range_r_p*ion->short_range_r,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDDH_Integral_xx_33_case_1( mix_a_p, mix_r_p, sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDDH_Integral_xx_33_case_2( mix_a_p, mix_r_p, dist, sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDDH_Integral_xx_33_case_3( mix_a_p, mix_r_p, sp->knot[i], dist,
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDDH_Integral_xx_33_case_4( mix_a_p, mix_r_p, sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}




// SDDXX ZZ == SDDYY ZZ
double sp_cluster_integrator_get_sh_xx_44_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp->short_range_a_s)+(ion->short_range_a))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp->short_range_a_p)+(ion->short_range_a))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp->short_range_r_s*ion->short_range_r,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp->short_range_r_p*ion->short_range_r,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDDH_Integral_xx_44_case_1( mix_a_p, mix_r_p, sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDDH_Integral_xx_44_case_2( mix_a_p, mix_r_p, dist, sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDDH_Integral_xx_44_case_3( mix_a_p, mix_r_p, sp->knot[i], dist,
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDDH_Integral_xx_44_case_4( mix_a_p, mix_r_p, sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}





// SDDXY XY == SDDXY XY
double sp_cluster_integrator_get_sh_xy_23_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp->short_range_a_s)+(ion->short_range_a))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp->short_range_a_p)+(ion->short_range_a))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp->short_range_r_s*ion->short_range_r,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp->short_range_r_p*ion->short_range_r,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDDH_Integral_xy_23_case_1( mix_a_p, mix_r_p, sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDDH_Integral_xy_23_case_2( mix_a_p, mix_r_p, dist, sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDDH_Integral_xy_23_case_3( mix_a_p, mix_r_p, sp->knot[i], dist,
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDDH_Integral_xy_23_case_4( mix_a_p, mix_r_p, sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}






// SDDXZ SX == SDDYZ SY
double sp_cluster_integrator_get_sh_xz_12_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp->short_range_a_s)+(ion->short_range_a))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp->short_range_a_p)+(ion->short_range_a))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp->short_range_r_s*ion->short_range_r,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp->short_range_r_p*ion->short_range_r,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDDH_Integral_xz_12_case_1( 0.5*(mix_a_s+mix_a_p), pow(mix_r_s*mix_r_p,0.5), sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDDH_Integral_xz_12_case_2( 0.5*(mix_a_s+mix_a_p), pow(mix_r_s*mix_r_p,0.5), dist, sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDDH_Integral_xz_12_case_3( 0.5*(mix_a_s+mix_a_p), pow(mix_r_s*mix_r_p,0.5), sp->knot[i], dist,
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDDH_Integral_xz_12_case_4( 0.5*(mix_a_s+mix_a_p), pow(mix_r_s*mix_r_p,0.5), sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}



// SDDXZ XZ == SDDYZ YZ
double sp_cluster_integrator_get_sh_xz_24_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp->short_range_a_s)+(ion->short_range_a))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp->short_range_a_p)+(ion->short_range_a))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp->short_range_r_s*ion->short_range_r,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp->short_range_r_p*ion->short_range_r,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDDH_Integral_xz_24_case_1( mix_a_p, mix_r_p, sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDDH_Integral_xz_24_case_2( mix_a_p, mix_r_p, dist, sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDDH_Integral_xz_24_case_3( mix_a_p, mix_r_p, sp->knot[i], dist,
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDDH_Integral_xz_24_case_4( mix_a_p, mix_r_p, sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}



// SDDZZ SS
double sp_cluster_integrator_get_sh_zz_11_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp->short_range_a_s)+(ion->short_range_a))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp->short_range_a_p)+(ion->short_range_a))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp->short_range_r_s*ion->short_range_r,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp->short_range_r_p*ion->short_range_r,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDDH_Integral_zz_11_case_1( mix_a_s, mix_r_s, sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDDH_Integral_zz_11_case_2( mix_a_s, mix_r_s, dist, sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDDH_Integral_zz_11_case_3( mix_a_s, mix_r_s, sp->knot[i], dist,
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDDH_Integral_zz_11_case_4( mix_a_s, mix_r_s, sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}




// SDDZZ SZ
double sp_cluster_integrator_get_sh_zz_14_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp->short_range_a_s)+(ion->short_range_a))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp->short_range_a_p)+(ion->short_range_a))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp->short_range_r_s*ion->short_range_r,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp->short_range_r_p*ion->short_range_r,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDDH_Integral_zz_14_case_1( 0.5*(mix_a_s+mix_a_p), pow(mix_r_s*mix_r_p,0.5), sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDDH_Integral_zz_14_case_2( 0.5*(mix_a_s+mix_a_p), pow(mix_r_s*mix_r_p,0.5), dist, sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDDH_Integral_zz_14_case_3( 0.5*(mix_a_s+mix_a_p), pow(mix_r_s*mix_r_p,0.5), sp->knot[i], dist,
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDDH_Integral_zz_14_case_4( 0.5*(mix_a_s+mix_a_p), pow(mix_r_s*mix_r_p,0.5), sp->knot[i], sp->knot[i+1],
                    sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1],
                    sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}




// SDDZZ XXYY
double sp_cluster_integrator_get_sh_zz_2233_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp->short_range_a_s)+(ion->short_range_a))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp->short_range_a_p)+(ion->short_range_a))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp->short_range_r_s*ion->short_range_r,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp->short_range_r_p*ion->short_range_r,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDDH_Integral_zz_2233_case_1( mix_a_p, mix_r_p, sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDDH_Integral_zz_2233_case_2( mix_a_p, mix_r_p, dist, sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDDH_Integral_zz_2233_case_3( mix_a_p, mix_r_p, sp->knot[i], dist,
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDDH_Integral_zz_2233_case_4( mix_a_p, mix_r_p, sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}




// SDDZZ zz
double sp_cluster_integrator_get_sh_zz_44_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp->core_position,0)-gsl_vector_get(ion->core_position,0),2.)
                + pow(gsl_vector_get(sp->core_position,1)-gsl_vector_get(ion->core_position,1),2.)
                + pow(gsl_vector_get(sp->core_position,2)-gsl_vector_get(ion->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp->short_range_a_s)+(ion->short_range_a))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp->short_range_a_p)+(ion->short_range_a))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp->short_range_r_s*ion->short_range_r,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp->short_range_r_p*ion->short_range_r,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp->number_of_knot-1;i++)
    {
        if( sp->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDDH_Integral_zz_44_case_1( mix_a_p, mix_r_p, sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( sp->knot[i] <= dist && dist < sp->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDDH_Integral_zz_44_case_2( mix_a_p, mix_r_p, dist, sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDDH_Integral_zz_44_case_3( mix_a_p, mix_r_p, sp->knot[i], dist,
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDDH_Integral_zz_44_case_4( mix_a_p, mix_r_p, sp->knot[i], sp->knot[i+1],
                    sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1],
                    sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}














///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     

//      x y and z Integrals : note that the integral observable is in Bohr Unit

double sp_cluster_integrator_get_x_12( sp_cluster_type_sp_ion* sp )
{   double Return = 0.;
    for(int i=0;i<sp->number_of_knot-1;i++)
    {   Return += Integral_x_12( sp->knot[i], sp->knot[i+1], sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1], sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1], sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3] );
    }
    return Return*TO_BOHR_RADII;    // into Ang Unit Bhor*0.529
}

double sp_cluster_integrator_get_y_13( sp_cluster_type_sp_ion* sp )
{   double Return = 0.;
    for(int i=0;i<sp->number_of_knot-1;i++)
    {   Return += Integral_y_13( sp->knot[i], sp->knot[i+1], sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1], sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1], sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3] );
    }
    return Return*TO_BOHR_RADII;
}

double sp_cluster_integrator_get_z_14( sp_cluster_type_sp_ion* sp )
{   double Return = 0.;
    for(int i=0;i<sp->number_of_knot-1;i++)
    {   Return += Integral_z_14( sp->knot[i], sp->knot[i+1], sp->radial_s_coefficient[i][0], sp->radial_s_coefficient[i][1], sp->radial_s_coefficient[i][2], sp->radial_s_coefficient[i][3],
                sp->radial_p_coefficient[i][0], sp->radial_p_coefficient[i][1], sp->radial_p_coefficient[i][2], sp->radial_p_coefficient[i][3] );
    }
    return Return*TO_BOHR_RADII;
}


// N Centre Integral Functions
//
// Naming convention ...spsp...


// Mono-pole related
//
// Coulomb
double sp_cluster_spsp_mono_integrator_get_ch_11_element( sp_cluster_type_sp_ion* sp1, sp_cluster_type_sp_ion* sp2 )
{
    // the integrals are carried over sp2
    // i.e., sp2 centre is pointing sp1 centre ..... sp2 -> sp1



    double Return = 0.;
    // get core(sp1) - core(sp2) distance
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {   // estimating charge of core & monopole shell together (w.r.t sp1 ion)
            Return += sp2->charge_shell*(sp1->charge_shell)*CH_Integral_11_case_1( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp2->charge_shell*(sp1->charge_shell)*CH_Integral_11_case_2_sub_1( dist, sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp2->charge_shell*(sp1->charge_shell)*CH_Integral_11_case_2_sub_2( sp2->knot[i], dist,
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp2->charge_shell*(sp1->charge_shell)*CH_Integral_11_case_3( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], dist );
        }
    }
    return Return*HA_TO_EV_UNIT;
}




double sp_cluster_spsp_mono_integrator_get_ch_14_element( sp_cluster_type_sp_ion* sp1, sp_cluster_type_sp_ion* sp2 )
{
    double Return = 0.;
    // get core(sp1) - core(sp2) distance
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {   // estimating charge of core & monopole shell together (w.r.t sp1 ion)
            Return += sp2->charge_shell*(sp1->charge_shell)*CH_Integral_14_case_1( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], 
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp2->charge_shell*(sp1->charge_shell)*CH_Integral_14_case_2_sub_1( dist, sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], 
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp2->charge_shell*(sp1->charge_shell)*CH_Integral_14_case_2_sub_2( sp2->knot[i], dist,
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], 
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp2->charge_shell*(sp1->charge_shell)*CH_Integral_14_case_3( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], 
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*HA_TO_EV_UNIT;
}




double sp_cluster_spsp_mono_integrator_get_ch_2233_element( sp_cluster_type_sp_ion* sp1, sp_cluster_type_sp_ion* sp2 )
{
    double Return = 0.;
    // get core(sp1) - core(sp2) distance
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {   // estimating charge of core & monopole shell together (w.r.t sp1 ion)
            Return += sp2->charge_shell*(sp1->charge_shell)*CH_Integral_2233_case_1( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp2->charge_shell*(sp1->charge_shell)*CH_Integral_2233_case_2_sub_1( dist, sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp2->charge_shell*(sp1->charge_shell)*CH_Integral_2233_case_2_sub_2( sp2->knot[i], dist,
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp2->charge_shell*(sp1->charge_shell)*CH_Integral_2233_case_3( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*HA_TO_EV_UNIT;
}




double sp_cluster_spsp_mono_integrator_get_ch_44_element( sp_cluster_type_sp_ion* sp1, sp_cluster_type_sp_ion* sp2 )
{
    double Return = 0.;
    // get core(sp1) - core(sp2) distance
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {   // estimating charge of core & monopole shell together (w.r.t sp1 ion)
            Return += sp2->charge_shell*(sp1->charge_shell)*CH_Integral_44_case_1( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp2->charge_shell*(sp1->charge_shell)*CH_Integral_44_case_2_sub_1( dist, sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp2->charge_shell*(sp1->charge_shell)*CH_Integral_44_case_2_sub_2( sp2->knot[i], dist,
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp2->charge_shell*(sp1->charge_shell)*CH_Integral_44_case_3( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*HA_TO_EV_UNIT;
}




// BM
double sp_cluster_spsp_mono_integrator_get_sh_11_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, sp_cluster_type_sp_ion* sp2 )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp1->short_range_a_s)+(sp2->short_range_a_s))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp1->short_range_a_p)+(sp2->short_range_a_p))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp1->short_range_r_s*sp2->short_range_r_s,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp1->short_range_r_p*sp2->short_range_r_p,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;

    int knot, ion_type;   double d, c, b, a;
    if( sp_sys->integral_lut == SP_INTEGRAL_TRUE )  // if LUT exist
    {
        dist = dist*TO_BOHR_RADII;  // back to Ang Unit
        knot = sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot ); 
 
        d = sp_sys->integral_vs_sp_s_ss[knot][0];
        c = sp_sys->integral_vs_sp_s_ss[knot][1];
        b = sp_sys->integral_vs_sp_s_ss[knot][2];
        a = sp_sys->integral_vs_sp_s_ss[knot][3];

        Return = d*pow(dist,3.) + c*pow(dist,2.) + b*dist + a;

        return Return;
    }

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp1->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSH_Integral_11_case_1( mix_a_s, mix_r_s, sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSH_Integral_11_case_2_sub_1( mix_a_s, mix_r_s, dist, sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSH_Integral_11_case_2_sub_2( mix_a_s, mix_r_s, sp2->knot[i], dist,
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSH_Integral_11_case_3( mix_a_s, mix_r_s, sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], dist );
        }
    }
    return Return*HA_TO_EV_UNIT;
}




double sp_cluster_spsp_mono_integrator_get_sh_14_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, sp_cluster_type_sp_ion* sp2 )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp1->short_range_a_s)+(sp2->short_range_a_s))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp1->short_range_a_p)+(sp2->short_range_a_p))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp1->short_range_r_s*sp2->short_range_r_s,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp1->short_range_r_p*sp2->short_range_r_p,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;

    int knot, ion_type;   double d, c, b, a;
    if( sp_sys->integral_lut == SP_INTEGRAL_TRUE )  // if LUT exist
    {
        dist = dist*TO_BOHR_RADII;  // back to Ang Unit
        knot = sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot ); 
 
        d = sp_sys->integral_vs_sp_s_sz[knot][0];
        c = sp_sys->integral_vs_sp_s_sz[knot][1];
        b = sp_sys->integral_vs_sp_s_sz[knot][2];
        a = sp_sys->integral_vs_sp_s_sz[knot][3];

        Return = d*pow(dist,3.) + c*pow(dist,2.) + b*dist + a;

        return Return;
    }

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp1->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSH_Integral_14_case_1( 0.5*(mix_a_s+mix_a_p) , pow(mix_r_s*mix_r_p,0.5) , sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], 
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSH_Integral_14_case_2_sub_1( 0.5*(mix_a_s+mix_a_p) , pow(mix_r_s*mix_r_p,0.5) , dist, sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], 
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSH_Integral_14_case_2_sub_2( 0.5*(mix_a_s+mix_a_p) , pow(mix_r_s*mix_r_p,0.5) , sp2->knot[i], dist,
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], 
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSH_Integral_14_case_3( 0.5*(mix_a_s+mix_a_p) , pow(mix_r_s*mix_r_p,0.5) , sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], 
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*HA_TO_EV_UNIT;
}




double sp_cluster_spsp_mono_integrator_get_sh_2233_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, sp_cluster_type_sp_ion* sp2 )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp1->short_range_a_s)+(sp2->short_range_a_s))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp1->short_range_a_p)+(sp2->short_range_a_p))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp1->short_range_r_s*sp2->short_range_r_s,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp1->short_range_r_p*sp2->short_range_r_p,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;

    int knot, ion_type;   double d, c, b, a;
    if( sp_sys->integral_lut == SP_INTEGRAL_TRUE )  // if LUT exist
    {
        dist = dist*TO_BOHR_RADII;  // back to Ang Unit
        knot = sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot ); 
 
        d = sp_sys->integral_vs_sp_s_xxyy[knot][0];
        c = sp_sys->integral_vs_sp_s_xxyy[knot][1];
        b = sp_sys->integral_vs_sp_s_xxyy[knot][2];
        a = sp_sys->integral_vs_sp_s_xxyy[knot][3];

        Return = d*pow(dist,3.) + c*pow(dist,2.) + b*dist + a;

        return Return;
    }

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp1->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSH_Integral_2233_case_1( mix_a_s, mix_r_s, sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSH_Integral_2233_case_2_sub_1( mix_a_s, mix_r_s, dist, sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSH_Integral_2233_case_2_sub_2( mix_a_s, mix_r_s, sp2->knot[i], dist,
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSH_Integral_2233_case_3( mix_a_s, mix_r_s, sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*HA_TO_EV_UNIT;
}




double sp_cluster_spsp_mono_integrator_get_sh_44_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, sp_cluster_type_sp_ion* sp2 )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp1->short_range_a_s)+(sp2->short_range_a_s))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp1->short_range_a_p)+(sp2->short_range_a_p))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp1->short_range_r_s*sp2->short_range_r_s,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp1->short_range_r_p*sp2->short_range_r_p,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;

    int knot, ion_type;   double d, c, b, a;
    if( sp_sys->integral_lut == SP_INTEGRAL_TRUE )  // if LUT exist
    {
        dist = dist*TO_BOHR_RADII;  // back to Ang Unit
        knot = sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot ); 
 
        d = sp_sys->integral_vs_sp_s_zz[knot][0];
        c = sp_sys->integral_vs_sp_s_zz[knot][1];
        b = sp_sys->integral_vs_sp_s_zz[knot][2];
        a = sp_sys->integral_vs_sp_s_zz[knot][3];

        Return = d*pow(dist,3.) + c*pow(dist,2.) + b*dist + a;

        return Return;
    }

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp1->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSH_Integral_44_case_1( mix_a_s, mix_r_s, sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSH_Integral_44_case_2_sub_1( mix_a_s, mix_r_s, dist, sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSH_Integral_44_case_2_sub_2( mix_a_s, mix_r_s, sp2->knot[i], dist,
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSH_Integral_44_case_3( mix_a_s, mix_r_s, sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*HA_TO_EV_UNIT;
}
//
// Mono-pole relaated Done












// Di-pole related
//
// Coulomb
double sp_cluster_spsp_di_integrator_get_ch_x_12_element( sp_cluster_type_sp_ion* sp1, sp_cluster_type_sp_ion* sp2 )
{   
    double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDH_Integral_x_12_case_1( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp2->charge_shell*sp1->charge_shell*CDH_Integral_x_12_case_2_sub_1( dist, sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp2->charge_shell*sp1->charge_shell*CDH_Integral_x_12_case_2_sub_2( sp2->knot[i], dist,
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDH_Integral_x_12_case_3( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FHA_TO_FEV_UNIT;
}



double sp_cluster_spsp_di_integrator_get_ch_x_24_element( sp_cluster_type_sp_ion* sp1, sp_cluster_type_sp_ion* sp2 )
{   
    double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDH_Integral_x_24_case_1( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp2->charge_shell*sp1->charge_shell*CDH_Integral_x_24_case_2_sub_1( dist, sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp2->charge_shell*sp1->charge_shell*CDH_Integral_x_24_case_2_sub_2( sp2->knot[i], dist,
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDH_Integral_x_24_case_3( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FHA_TO_FEV_UNIT;
}





// ch_x_12 == ch_y_13 & ch_x_24 == ch_y_34
double sp_cluster_spsp_di_integrator_get_ch_z_11_element( sp_cluster_type_sp_ion* sp1, sp_cluster_type_sp_ion* sp2 )
{   
    double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDH_Integral_z_11_case_1( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp2->charge_shell*sp1->charge_shell*CDH_Integral_z_11_case_2_sub_1( dist, sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp2->charge_shell*sp1->charge_shell*CDH_Integral_z_11_case_2_sub_2( sp2->knot[i], dist,
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDH_Integral_z_11_case_3( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], dist );
        }
    }
    return Return*FHA_TO_FEV_UNIT;
}




double sp_cluster_spsp_di_integrator_get_ch_z_14_element( sp_cluster_type_sp_ion* sp1, sp_cluster_type_sp_ion* sp2 )
{   
    double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDH_Integral_z_14_case_1( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp2->charge_shell*sp1->charge_shell*CDH_Integral_z_14_case_2_sub_1( dist, sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp2->charge_shell*sp1->charge_shell*CDH_Integral_z_14_case_2_sub_2( sp2->knot[i], dist,
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDH_Integral_z_14_case_3( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FHA_TO_FEV_UNIT;
}




double sp_cluster_spsp_di_integrator_get_ch_z_2233_element( sp_cluster_type_sp_ion* sp1, sp_cluster_type_sp_ion* sp2 )
{   
    double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDH_Integral_z_2233_case_1( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp2->charge_shell*sp1->charge_shell*CDH_Integral_z_2233_case_2_sub_1( dist, sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp2->charge_shell*sp1->charge_shell*CDH_Integral_z_2233_case_2_sub_2( sp2->knot[i], dist,
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDH_Integral_z_2233_case_3( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FHA_TO_FEV_UNIT;
}




double sp_cluster_spsp_di_integrator_get_ch_z_44_element( sp_cluster_type_sp_ion* sp1, sp_cluster_type_sp_ion* sp2 )
{   
    double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDH_Integral_z_44_case_1( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp2->charge_shell*sp1->charge_shell*CDH_Integral_z_44_case_2_sub_1( dist, sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp2->charge_shell*sp1->charge_shell*CDH_Integral_z_44_case_2_sub_2( sp2->knot[i], dist,
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDH_Integral_z_44_case_3( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FHA_TO_FEV_UNIT;
}







// BM
double sp_cluster_spsp_di_integrator_get_sh_x_12_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, sp_cluster_type_sp_ion* sp2 )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp1->short_range_a_s)+(sp2->short_range_a_s))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp1->short_range_a_p)+(sp2->short_range_a_p))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp1->short_range_r_s*sp2->short_range_r_s,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp1->short_range_r_p*sp2->short_range_r_p,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;

    int knot, ion_type;   double d, c, b, a;
    if( sp_sys->integral_lut == SP_INTEGRAL_TRUE )  // if LUT exist
    {
        dist = dist*TO_BOHR_RADII;  // back to Ang Unit
        knot = sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot ); 
 
        d = sp_sys->integral_vs_sp_s_x_sx[knot][0];
        c = sp_sys->integral_vs_sp_s_x_sx[knot][1];
        b = sp_sys->integral_vs_sp_s_x_sx[knot][2];
        a = sp_sys->integral_vs_sp_s_x_sx[knot][3];

        Return = d*pow(dist,3.) + c*pow(dist,2.) + b*dist + a;

        return Return;
    }


    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDH_Integral_x_12_case_1( (mix_a_s+mix_a_p)/2., pow(mix_r_s*mix_r_p,0.5), sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDH_Integral_x_12_case_2_sub_1( (mix_a_s+mix_a_p)/2., pow(mix_r_s*mix_r_p,0.5), dist, sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDH_Integral_x_12_case_2_sub_2( (mix_a_s+mix_a_p)/2., pow(mix_r_s*mix_r_p,0.5), sp2->knot[i], dist,
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDH_Integral_x_12_case_3( (mix_a_s+mix_a_p)/2., pow(mix_r_s*mix_r_p,0.5), sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FHA_TO_FEV_UNIT;
}




double sp_cluster_spsp_di_integrator_get_sh_x_24_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, sp_cluster_type_sp_ion* sp2 )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp1->short_range_a_s)+(sp2->short_range_a_s))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp1->short_range_a_p)+(sp2->short_range_a_p))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp1->short_range_r_s*sp2->short_range_r_s,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp1->short_range_r_p*sp2->short_range_r_p,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;

    int knot, ion_type;   double d, c, b, a;
    if( sp_sys->integral_lut == SP_INTEGRAL_TRUE )  // if LUT exist
    {
        dist = dist*TO_BOHR_RADII;  // back to Ang Unit
        knot = sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot ); 
 
        d = sp_sys->integral_vs_sp_s_x_xz[knot][0];
        c = sp_sys->integral_vs_sp_s_x_xz[knot][1];
        b = sp_sys->integral_vs_sp_s_x_xz[knot][2];
        a = sp_sys->integral_vs_sp_s_x_xz[knot][3];

        Return = d*pow(dist,3.) + c*pow(dist,2.) + b*dist + a;

        return Return;
    }

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDH_Integral_x_24_case_1( mix_a_p, mix_r_p, sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDH_Integral_x_24_case_2_sub_1( mix_a_p, mix_r_p, dist, sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDH_Integral_x_24_case_2_sub_2( mix_a_p, mix_r_p, sp2->knot[i], dist,
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDH_Integral_x_24_case_3( mix_a_p, mix_r_p, sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FHA_TO_FEV_UNIT;
}




// sh_x_12 == sh_y_13 & sh_x_24 == sh_y_34
double sp_cluster_spsp_di_integrator_get_sh_z_11_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, sp_cluster_type_sp_ion* sp2 )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp1->short_range_a_s)+(sp2->short_range_a_s))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp1->short_range_a_p)+(sp2->short_range_a_p))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp1->short_range_r_s*sp2->short_range_r_s,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp1->short_range_r_p*sp2->short_range_r_p,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;

    int knot, ion_type;   double d, c, b, a;
    if( sp_sys->integral_lut == SP_INTEGRAL_TRUE )  // if LUT exist
    {
        dist = dist*TO_BOHR_RADII;  // back to Ang Unit
        knot = sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot ); 
 
        d = sp_sys->integral_vs_sp_s_z_ss[knot][0];
        c = sp_sys->integral_vs_sp_s_z_ss[knot][1];
        b = sp_sys->integral_vs_sp_s_z_ss[knot][2];
        a = sp_sys->integral_vs_sp_s_z_ss[knot][3];

        Return = d*pow(dist,3.) + c*pow(dist,2.) + b*dist + a;

        return Return;
    }

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDH_Integral_z_11_case_1( mix_a_s, mix_r_s, sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDH_Integral_z_11_case_2_sub_1( mix_a_s, mix_r_s, dist, sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDH_Integral_z_11_case_2_sub_2( mix_a_s, mix_r_s, sp2->knot[i], dist,
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDH_Integral_z_11_case_3( mix_a_s, mix_r_s, sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], dist );
        }
    }
    return Return*FHA_TO_FEV_UNIT;
}




double sp_cluster_spsp_di_integrator_get_sh_z_14_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, sp_cluster_type_sp_ion* sp2 )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp1->short_range_a_s)+(sp2->short_range_a_s))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp1->short_range_a_p)+(sp2->short_range_a_p))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp1->short_range_r_s*sp2->short_range_r_s,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp1->short_range_r_p*sp2->short_range_r_p,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;
  
    int knot, ion_type;   double d, c, b, a;
    if( sp_sys->integral_lut == SP_INTEGRAL_TRUE )  // if LUT exist
    {
        dist = dist*TO_BOHR_RADII;  // back to Ang Unit
        knot = sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot ); 
 
        d = sp_sys->integral_vs_sp_s_z_sz[knot][0];
        c = sp_sys->integral_vs_sp_s_z_sz[knot][1];
        b = sp_sys->integral_vs_sp_s_z_sz[knot][2];
        a = sp_sys->integral_vs_sp_s_z_sz[knot][3];

        Return = d*pow(dist,3.) + c*pow(dist,2.) + b*dist + a;

        return Return;
    }

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDH_Integral_z_14_case_1( (mix_a_s+mix_a_p)/2., pow(mix_r_s*mix_r_p,0.5), sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDH_Integral_z_14_case_2_sub_1( (mix_a_s+mix_a_p)/2., pow(mix_r_s*mix_r_p,0.5), dist, sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDH_Integral_z_14_case_2_sub_2( (mix_a_s+mix_a_p)/2., pow(mix_r_s*mix_r_p,0.5), sp2->knot[i], dist,
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDH_Integral_z_14_case_3( (mix_a_s+mix_a_p)/2., pow(mix_r_s*mix_r_p,0.5), sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FHA_TO_FEV_UNIT;
}




double sp_cluster_spsp_di_integrator_get_sh_z_2233_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, sp_cluster_type_sp_ion* sp2 )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp1->short_range_a_s)+(sp2->short_range_a_s))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp1->short_range_a_p)+(sp2->short_range_a_p))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp1->short_range_r_s*sp2->short_range_r_s,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp1->short_range_r_p*sp2->short_range_r_p,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;

    int knot, ion_type;   double d, c, b, a;
    if( sp_sys->integral_lut == SP_INTEGRAL_TRUE )  // if LUT exist
    {
        dist = dist*TO_BOHR_RADII;  // back to Ang Unit
        knot = sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot ); 
 
        d = sp_sys->integral_vs_sp_s_z_xxyy[knot][0];
        c = sp_sys->integral_vs_sp_s_z_xxyy[knot][1];
        b = sp_sys->integral_vs_sp_s_z_xxyy[knot][2];
        a = sp_sys->integral_vs_sp_s_z_xxyy[knot][3];

        Return = d*pow(dist,3.) + c*pow(dist,2.) + b*dist + a;

        return Return;
    }

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDH_Integral_z_2233_case_1( mix_a_p, mix_r_p, sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDH_Integral_z_2233_case_2_sub_1( mix_a_p, mix_r_p, dist, sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDH_Integral_z_2233_case_2_sub_2( mix_a_p, mix_r_p, sp2->knot[i], dist,
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDH_Integral_z_2233_case_3( mix_a_p, mix_r_p, sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FHA_TO_FEV_UNIT;
}



double sp_cluster_spsp_di_integrator_get_sh_z_44_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, sp_cluster_type_sp_ion* sp2 )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp1->short_range_a_s)+(sp2->short_range_a_s))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp1->short_range_a_p)+(sp2->short_range_a_p))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp1->short_range_r_s*sp2->short_range_r_s,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp1->short_range_r_p*sp2->short_range_r_p,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;

    int knot, ion_type;   double d, c, b, a;
    if( sp_sys->integral_lut == SP_INTEGRAL_TRUE )  // if LUT exist
    {
        dist = dist*TO_BOHR_RADII;  // back to Ang Unit
        knot = sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot ); 
 
        d = sp_sys->integral_vs_sp_s_z_zz[knot][0];
        c = sp_sys->integral_vs_sp_s_z_zz[knot][1];
        b = sp_sys->integral_vs_sp_s_z_zz[knot][2];
        a = sp_sys->integral_vs_sp_s_z_zz[knot][3];

        Return = d*pow(dist,3.) + c*pow(dist,2.) + b*dist + a;

        return Return;
    }

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDH_Integral_z_44_case_1( mix_a_p, mix_r_p, sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDH_Integral_z_44_case_2_sub_1( mix_a_p, mix_r_p, dist, sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDH_Integral_z_44_case_2_sub_2( mix_a_p, mix_r_p, sp2->knot[i], dist,
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDH_Integral_z_44_case_3( mix_a_p, mix_r_p, sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FHA_TO_FEV_UNIT;
}
//
// Di-pole related Done




///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///
// Quadru-pole related


/* ************************************* 2nd Derivative Methods ************************************* */

// SECOND_DERIVATIVE_COULOMB_LONG_RANGE

// CDDH_XX

// SS   ==  CDDH_YY_SS

double sp_cluster_spsp_quadru_integrator_get_ch_xx_11_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2)
{   
    double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_xx_11_case_1( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_xx_11_case_2_sub_1( dist, sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_xx_11_case_2_sub_2( sp2->knot[i], dist,
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_xx_11_case_3( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}

// SZ   == CDDH_YY_SZ
double sp_cluster_spsp_quadru_integrator_get_ch_xx_14_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_xx_14_case_1( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_xx_14_case_2_sub_1( dist, sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_xx_14_case_2_sub_2( sp2->knot[i], dist,
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_xx_14_case_3( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}




// XX   == CDDH_YY_YY
double sp_cluster_spsp_quadru_integrator_get_ch_xx_22_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_xx_22_case_1( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_xx_22_case_2_sub_1( dist, sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_xx_22_case_2_sub_2( sp2->knot[i], dist,
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_xx_22_case_3( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}


// YY   == CDDH_YY_XX
double sp_cluster_spsp_quadru_integrator_get_ch_xx_33_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_xx_33_case_1( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_xx_33_case_2_sub_1( dist, sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_xx_33_case_2_sub_2( sp2->knot[i], dist,
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_xx_33_case_3( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}


// ZZ   == CDDH_YY_ZZ
double sp_cluster_spsp_quadru_integrator_get_ch_xx_44_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_xx_44_case_1( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_xx_44_case_2_sub_1( dist, sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_xx_44_case_2_sub_2( sp2->knot[i], dist,
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_xx_44_case_3( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}






// CDDH_XY

// XY
double sp_cluster_spsp_quadru_integrator_get_ch_xy_23_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_xy_23_case_1( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_xy_23_case_2_sub_1( dist, sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_xy_23_case_2_sub_2( sp2->knot[i], dist,
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_xy_23_case_3( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}







// CDDH_XZ

// SX   == CDDH_YZ_SY
double sp_cluster_spsp_quadru_integrator_get_ch_xz_12_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_xz_12_case_1( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_xz_12_case_2_sub_1( dist, sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_xz_12_case_2_sub_2( sp2->knot[i], dist,
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_xz_12_case_3( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}


// XZ   == CDDH_YZ_YZ
double sp_cluster_spsp_quadru_integrator_get_ch_xz_24_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_xz_24_case_1( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_xz_24_case_2_sub_1( dist, sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_xz_24_case_2_sub_2( sp2->knot[i], dist,
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_xz_24_case_3( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}









// CDDH_ZZ

// SS 
double sp_cluster_spsp_quadru_integrator_get_ch_zz_11_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 )
{   
    double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_zz_11_case_1( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_zz_11_case_2_sub_1( dist, sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_zz_11_case_2_sub_2( sp2->knot[i], dist,
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_zz_11_case_3( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}

// SZ
double sp_cluster_spsp_quadru_integrator_get_ch_zz_14_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_zz_14_case_1( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_zz_14_case_2_sub_1( dist, sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_zz_14_case_2_sub_2( sp2->knot[i], dist,
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_zz_14_case_3( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}




// XXYY
double sp_cluster_spsp_quadru_integrator_get_ch_zz_2233_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_zz_2233_case_1( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_zz_2233_case_2_sub_1( dist, sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_zz_2233_case_2_sub_2( sp2->knot[i], dist,
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_zz_2233_case_3( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}




// ZZ
double sp_cluster_spsp_quadru_integrator_get_ch_zz_44_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=0;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_zz_44_case_1( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_zz_44_case_2_sub_1( dist, sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_zz_44_case_2_sub_2( sp2->knot[i], dist,
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += sp2->charge_shell*sp1->charge_shell*CDDH_Integral_zz_44_case_3( sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}



//// Second Derivative
//
//  Short-Range Repulsion



// SDDXX SS
double sp_cluster_spsp_quadru_integrator_get_sh_xx_11_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp1->short_range_a_s)+(sp2->short_range_a_s))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp1->short_range_a_p)+(sp2->short_range_a_p))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp1->short_range_r_s*sp2->short_range_r_s,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp1->short_range_r_p*sp2->short_range_r_p,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;

    int knot, ion_type;   double d, c, b, a;
    if( sp_sys->integral_lut == SP_INTEGRAL_TRUE )  // if LUT exist
    {
        dist = dist*TO_BOHR_RADII;  // back to Ang Unit
        knot = sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot ); 
 
        d = sp_sys->integral_vs_sp_s_xx_ss[knot][0];
        c = sp_sys->integral_vs_sp_s_xx_ss[knot][1];
        b = sp_sys->integral_vs_sp_s_xx_ss[knot][2];
        a = sp_sys->integral_vs_sp_s_xx_ss[knot][3];

        Return = d*pow(dist,3.) + c*pow(dist,2.) + b*dist + a;

        return Return;
    }

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDDH_Integral_xx_11_case_1( mix_a_s, mix_r_s, sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDDH_Integral_xx_11_case_2( mix_a_s, mix_r_s, dist, sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDDH_Integral_xx_11_case_3( mix_a_s, mix_r_s, sp2->knot[i], dist,
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDDH_Integral_xx_11_case_4( mix_a_s, mix_r_s, sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}



// SDDXX SZ
double sp_cluster_spsp_quadru_integrator_get_sh_xx_14_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp1->short_range_a_s)+(sp2->short_range_a_s))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp1->short_range_a_p)+(sp2->short_range_a_p))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp1->short_range_r_s*sp2->short_range_r_s,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp1->short_range_r_p*sp2->short_range_r_p,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;

    int knot, ion_type;   double d, c, b, a;
    if( sp_sys->integral_lut == SP_INTEGRAL_TRUE )  // if LUT exist
    {
        dist = dist*TO_BOHR_RADII;  // back to Ang Unit
        knot = sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot ); 
 
        d = sp_sys->integral_vs_sp_s_xx_sz[knot][0];
        c = sp_sys->integral_vs_sp_s_xx_sz[knot][1];
        b = sp_sys->integral_vs_sp_s_xx_sz[knot][2];
        a = sp_sys->integral_vs_sp_s_xx_sz[knot][3];

        Return = d*pow(dist,3.) + c*pow(dist,2.) + b*dist + a;

        return Return;
    }

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDDH_Integral_xx_14_case_1( 0.5*(mix_a_s+mix_a_p), pow(mix_r_s*mix_r_p,0.5), sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDDH_Integral_xx_14_case_2( 0.5*(mix_a_s+mix_a_p), pow(mix_r_s*mix_r_p,0.5), dist, sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDDH_Integral_xx_14_case_3( 0.5*(mix_a_s+mix_a_p), pow(mix_r_s*mix_r_p,0.5), sp2->knot[i], dist,
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDDH_Integral_xx_14_case_4( 0.5*(mix_a_s+mix_a_p), pow(mix_r_s*mix_r_p,0.5), sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}



// SDDXX XX == SDDYY YY
double sp_cluster_spsp_quadru_integrator_get_sh_xx_22_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp1->short_range_a_s)+(sp2->short_range_a_s))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp1->short_range_a_p)+(sp2->short_range_a_p))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp1->short_range_r_s*sp2->short_range_r_s,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp1->short_range_r_p*sp2->short_range_r_p,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;

    int knot, ion_type;   double d, c, b, a;
    if( sp_sys->integral_lut == SP_INTEGRAL_TRUE )  // if LUT exist
    {
        dist = dist*TO_BOHR_RADII;  // back to Ang Unit
        knot = sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot ); 
 
        d = sp_sys->integral_vs_sp_s_xx_xx[knot][0];
        c = sp_sys->integral_vs_sp_s_xx_xx[knot][1];
        b = sp_sys->integral_vs_sp_s_xx_xx[knot][2];
        a = sp_sys->integral_vs_sp_s_xx_xx[knot][3];

        Return = d*pow(dist,3.) + c*pow(dist,2.) + b*dist + a;

        return Return;
    }

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDDH_Integral_xx_22_case_1( mix_a_p, mix_r_p, sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDDH_Integral_xx_22_case_2( mix_a_p, mix_r_p, dist, sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDDH_Integral_xx_22_case_3( mix_a_p, mix_r_p, sp2->knot[i], dist,
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDDH_Integral_xx_22_case_4( mix_a_p, mix_r_p, sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}



// SDDXX YY == SDDYY XX
double sp_cluster_spsp_quadru_integrator_get_sh_xx_33_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp1->short_range_a_s)+(sp2->short_range_a_s))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp1->short_range_a_p)+(sp2->short_range_a_p))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp1->short_range_r_s*sp2->short_range_r_s,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp1->short_range_r_p*sp2->short_range_r_p,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;

    int knot, ion_type;   double d, c, b, a;
    if( sp_sys->integral_lut == SP_INTEGRAL_TRUE )  // if LUT exist
    {
        dist = dist*TO_BOHR_RADII;  // back to Ang Unit
        knot = sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot ); 
 
        d = sp_sys->integral_vs_sp_s_xx_yy[knot][0];
        c = sp_sys->integral_vs_sp_s_xx_yy[knot][1];
        b = sp_sys->integral_vs_sp_s_xx_yy[knot][2];
        a = sp_sys->integral_vs_sp_s_xx_yy[knot][3];

        Return = d*pow(dist,3.) + c*pow(dist,2.) + b*dist + a;

        return Return;
    }

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDDH_Integral_xx_33_case_1( mix_a_p, mix_r_p, sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDDH_Integral_xx_33_case_2( mix_a_p, mix_r_p, dist, sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDDH_Integral_xx_33_case_3( mix_a_p, mix_r_p, sp2->knot[i], dist,
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDDH_Integral_xx_33_case_4( mix_a_p, mix_r_p, sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}




// SDDXX ZZ == SDDYY ZZ
double sp_cluster_spsp_quadru_integrator_get_sh_xx_44_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp1->short_range_a_s)+(sp2->short_range_a_s))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp1->short_range_a_p)+(sp2->short_range_a_p))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp1->short_range_r_s*sp2->short_range_r_s,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp1->short_range_r_p*sp2->short_range_r_p,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;
 
    int knot, ion_type;   double d, c, b, a;
    if( sp_sys->integral_lut == SP_INTEGRAL_TRUE )  // if LUT exist
    {
        dist = dist*TO_BOHR_RADII;  // back to Ang Unit
        knot = sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot ); 
 
        d = sp_sys->integral_vs_sp_s_xx_zz[knot][0];
        c = sp_sys->integral_vs_sp_s_xx_zz[knot][1];
        b = sp_sys->integral_vs_sp_s_xx_zz[knot][2];
        a = sp_sys->integral_vs_sp_s_xx_zz[knot][3];

        Return = d*pow(dist,3.) + c*pow(dist,2.) + b*dist + a;

        return Return;
    }

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDDH_Integral_xx_44_case_1( mix_a_p, mix_r_p, sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDDH_Integral_xx_44_case_2( mix_a_p, mix_r_p, dist, sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDDH_Integral_xx_44_case_3( mix_a_p, mix_r_p, sp2->knot[i], dist,
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDDH_Integral_xx_44_case_4( mix_a_p, mix_r_p, sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}





// SDDXY XY == SDDXY XY
double sp_cluster_spsp_quadru_integrator_get_sh_xy_23_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, sp_cluster_type_sp_ion* sp2 )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp1->short_range_a_s)+(sp2->short_range_a_s))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp1->short_range_a_p)+(sp2->short_range_a_p))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp1->short_range_r_s*sp2->short_range_r_s,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp1->short_range_r_p*sp2->short_range_r_p,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;

    int knot, ion_type;   double d, c, b, a;
    if( sp_sys->integral_lut == SP_INTEGRAL_TRUE )  // if LUT exist
    {
        dist = dist*TO_BOHR_RADII;  // back to Ang Unit
        knot = sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot ); 
 
        d = sp_sys->integral_vs_sp_s_xy_xy[knot][0];
        c = sp_sys->integral_vs_sp_s_xy_xy[knot][1];
        b = sp_sys->integral_vs_sp_s_xy_xy[knot][2];
        a = sp_sys->integral_vs_sp_s_xy_xy[knot][3];

        Return = d*pow(dist,3.) + c*pow(dist,2.) + b*dist + a;

        return Return;
    }

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDDH_Integral_xy_23_case_1( mix_a_p, mix_r_p, sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDDH_Integral_xy_23_case_2( mix_a_p, mix_r_p, dist, sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDDH_Integral_xy_23_case_3( mix_a_p, mix_r_p, sp2->knot[i], dist,
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDDH_Integral_xy_23_case_4( mix_a_p, mix_r_p, sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}






// SDDXZ SX == SDDYZ SY
double sp_cluster_spsp_quadru_integrator_get_sh_xz_12_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp1->short_range_a_s)+(sp2->short_range_a_s))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp1->short_range_a_p)+(sp2->short_range_a_p))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp1->short_range_r_s*sp2->short_range_r_s,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp1->short_range_r_p*sp2->short_range_r_p,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;

    int knot, ion_type;   double d, c, b, a;
    if( sp_sys->integral_lut == SP_INTEGRAL_TRUE )  // if LUT exist
    {
        dist = dist*TO_BOHR_RADII;  // back to Ang Unit
        knot = sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot ); 
 
        d = sp_sys->integral_vs_sp_s_xz_sx[knot][0];
        c = sp_sys->integral_vs_sp_s_xz_sx[knot][1];
        b = sp_sys->integral_vs_sp_s_xz_sx[knot][2];
        a = sp_sys->integral_vs_sp_s_xz_sx[knot][3];

        Return = d*pow(dist,3.) + c*pow(dist,2.) + b*dist + a;

        return Return;
    }

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDDH_Integral_xz_12_case_1( 0.5*(mix_a_s+mix_a_p), pow(mix_r_s*mix_r_p,0.5), sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDDH_Integral_xz_12_case_2( 0.5*(mix_a_s+mix_a_p), pow(mix_r_s*mix_r_p,0.5), dist, sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDDH_Integral_xz_12_case_3( 0.5*(mix_a_s+mix_a_p), pow(mix_r_s*mix_r_p,0.5), sp2->knot[i], dist,
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDDH_Integral_xz_12_case_4( 0.5*(mix_a_s+mix_a_p), pow(mix_r_s*mix_r_p,0.5), sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}



// SDDXZ XZ == SDDYZ YZ
double sp_cluster_spsp_quadru_integrator_get_sh_xz_24_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp1->short_range_a_s)+(sp2->short_range_a_s))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp1->short_range_a_p)+(sp2->short_range_a_p))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp1->short_range_r_s*sp2->short_range_r_s,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp1->short_range_r_p*sp2->short_range_r_p,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;

    int knot, ion_type;   double d, c, b, a;
    if( sp_sys->integral_lut == SP_INTEGRAL_TRUE )  // if LUT exist
    {
        dist = dist*TO_BOHR_RADII;  // back to Ang Unit
        knot = sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot ); 
 
        d = sp_sys->integral_vs_sp_s_xz_xz[knot][0];
        c = sp_sys->integral_vs_sp_s_xz_xz[knot][1];
        b = sp_sys->integral_vs_sp_s_xz_xz[knot][2];
        a = sp_sys->integral_vs_sp_s_xz_xz[knot][3];

        Return = d*pow(dist,3.) + c*pow(dist,2.) + b*dist + a;

        return Return;
    }

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDDH_Integral_xz_24_case_1( mix_a_p, mix_r_p, sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDDH_Integral_xz_24_case_2( mix_a_p, mix_r_p, dist, sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDDH_Integral_xz_24_case_3( mix_a_p, mix_r_p, sp2->knot[i], dist,
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDDH_Integral_xz_24_case_4( mix_a_p, mix_r_p, sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}



// SDDZZ SS
double sp_cluster_spsp_quadru_integrator_get_sh_zz_11_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, sp_cluster_type_sp_ion* sp2 )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp1->short_range_a_s)+(sp2->short_range_a_s))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp1->short_range_a_p)+(sp2->short_range_a_p))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp1->short_range_r_s*sp2->short_range_r_s,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp1->short_range_r_p*sp2->short_range_r_p,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;

    int knot, ion_type;   double d, c, b, a;
    if( sp_sys->integral_lut == SP_INTEGRAL_TRUE )  // if LUT exist
    {
        dist = dist*TO_BOHR_RADII;  // back to Ang Unit
        knot = sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot ); 
 
        d = sp_sys->integral_vs_sp_s_zz_ss[knot][0];
        c = sp_sys->integral_vs_sp_s_zz_ss[knot][1];
        b = sp_sys->integral_vs_sp_s_zz_ss[knot][2];
        a = sp_sys->integral_vs_sp_s_zz_ss[knot][3];

        Return = d*pow(dist,3.) + c*pow(dist,2.) + b*dist + a;

        return Return;
    }

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDDH_Integral_zz_11_case_1( mix_a_s, mix_r_s, sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDDH_Integral_zz_11_case_2( mix_a_s, mix_r_s, dist, sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDDH_Integral_zz_11_case_3( mix_a_s, mix_r_s, sp2->knot[i], dist,
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDDH_Integral_zz_11_case_4( mix_a_s, mix_r_s, sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}




// SDDZZ SZ
double sp_cluster_spsp_quadru_integrator_get_sh_zz_14_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, sp_cluster_type_sp_ion* sp2 )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp1->short_range_a_s)+(sp2->short_range_a_s))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp1->short_range_a_p)+(sp2->short_range_a_p))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp1->short_range_r_s*sp2->short_range_r_s,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp1->short_range_r_p*sp2->short_range_r_p,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;

    int knot, ion_type;   double d, c, b, a;
    if( sp_sys->integral_lut == SP_INTEGRAL_TRUE )  // if LUT exist
    {
        dist = dist*TO_BOHR_RADII;  // back to Ang Unit
        knot = sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot ); 
 
        d = sp_sys->integral_vs_sp_s_zz_sz[knot][0];
        c = sp_sys->integral_vs_sp_s_zz_sz[knot][1];
        b = sp_sys->integral_vs_sp_s_zz_sz[knot][2];
        a = sp_sys->integral_vs_sp_s_zz_sz[knot][3];

        Return = d*pow(dist,3.) + c*pow(dist,2.) + b*dist + a;

        return Return;
    }

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDDH_Integral_zz_14_case_1( 0.5*(mix_a_s+mix_a_p), pow(mix_r_s*mix_r_p,0.5), sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDDH_Integral_zz_14_case_2( 0.5*(mix_a_s+mix_a_p), pow(mix_r_s*mix_r_p,0.5), dist, sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDDH_Integral_zz_14_case_3( 0.5*(mix_a_s+mix_a_p), pow(mix_r_s*mix_r_p,0.5), sp2->knot[i], dist,
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDDH_Integral_zz_14_case_4( 0.5*(mix_a_s+mix_a_p), pow(mix_r_s*mix_r_p,0.5), sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_s_coefficient[i][0], sp2->radial_s_coefficient[i][1],
                    sp2->radial_s_coefficient[i][2], sp2->radial_s_coefficient[i][3],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}




// SDDZZ XXYY
double sp_cluster_spsp_quadru_integrator_get_sh_zz_2233_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp1->short_range_a_s)+(sp2->short_range_a_s))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp1->short_range_a_p)+(sp2->short_range_a_p))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp1->short_range_r_s*sp2->short_range_r_s,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp1->short_range_r_p*sp2->short_range_r_p,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;

    int knot, ion_type;   double d, c, b, a;
    if( sp_sys->integral_lut == SP_INTEGRAL_TRUE )  // if LUT exist
    {
        dist = dist*TO_BOHR_RADII;  // back to Ang Unit
        knot = sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot ); 
 
        d = sp_sys->integral_vs_sp_s_zz_xxyy[knot][0];
        c = sp_sys->integral_vs_sp_s_zz_xxyy[knot][1];
        b = sp_sys->integral_vs_sp_s_zz_xxyy[knot][2];
        a = sp_sys->integral_vs_sp_s_zz_xxyy[knot][3];

        Return = d*pow(dist,3.) + c*pow(dist,2.) + b*dist + a;

        return Return;
    }

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDDH_Integral_zz_2233_case_1( mix_a_p, mix_r_p, sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDDH_Integral_zz_2233_case_2( mix_a_p, mix_r_p, dist, sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDDH_Integral_zz_2233_case_3( mix_a_p, mix_r_p, sp2->knot[i], dist,
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDDH_Integral_zz_2233_case_4( mix_a_p, mix_r_p, sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}




// SDDZZ zz
double sp_cluster_spsp_quadru_integrator_get_sh_zz_44_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 )
{   double Return = 0.;
    // distance calc
    double dist = pow(pow(gsl_vector_get(sp1->core_position,0)-gsl_vector_get(sp2->core_position,0),2.)
                + pow(gsl_vector_get(sp1->core_position,1)-gsl_vector_get(sp2->core_position,1),2.)
                + pow(gsl_vector_get(sp1->core_position,2)-gsl_vector_get(sp2->core_position,2),2.),0.5)/TO_BOHR_RADII;
    // Short Range Parameter Following Standard Mixing Rule
    double mix_a_s = ((sp1->short_range_a_s)+(sp2->short_range_a_s))/2./HA_TO_EV_UNIT;
    double mix_a_p = ((sp1->short_range_a_p)+(sp2->short_range_a_p))/2./HA_TO_EV_UNIT;   // A
    double mix_r_s = pow(sp1->short_range_r_s*sp2->short_range_r_s,0.5)/TO_BOHR_RADII;
    double mix_r_p = pow(sp1->short_range_r_p*sp2->short_range_r_p,0.5)/TO_BOHR_RADII;   // Rho
    // Integral Tolerance
    int int_tol = SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE;

    if( dist > SP_INTEGRAL_SHORT_RANGE_CUTOFF/TO_BOHR_RADII )
        return 0.;

    int knot, ion_type;   double d, c, b, a;
    if( sp_sys->integral_lut == SP_INTEGRAL_TRUE )  // if LUT exist
    {
        dist = dist*TO_BOHR_RADII;  // back to Ang Unit
        knot = sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot ); 
 
        d = sp_sys->integral_vs_sp_s_zz_zz[knot][0];
        c = sp_sys->integral_vs_sp_s_zz_zz[knot][1];
        b = sp_sys->integral_vs_sp_s_zz_zz[knot][2];
        a = sp_sys->integral_vs_sp_s_zz_zz[knot][3];

        Return = d*pow(dist,3.) + c*pow(dist,2.) + b*dist + a;

        return Return;
    }

    // Calc partial integrals from the integral libray by the knot vs dist criteria
    for(int i=int_tol;i<sp2->number_of_knot-1;i++)
    {
        if( sp2->knot[i+1] <= dist )  // case 1: k1 < k2 < d
        {
            Return += BMSDDH_Integral_zz_44_case_1( mix_a_p, mix_r_p, sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( sp2->knot[i] <= dist && dist < sp2->knot[i+1] ) // case 2: k1 < d < k2
        {
            // sub case 2-1: d < r < k2
            Return += BMSDDH_Integral_zz_44_case_2( mix_a_p, mix_r_p, dist, sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
            // sub case 2-2: k1 < r < d
            Return += BMSDDH_Integral_zz_44_case_3( mix_a_p, mix_r_p, sp2->knot[i], dist,
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
        else if( dist < sp2->knot[i] ) // case 3: k1 < k2 < d
        {
            Return += BMSDDH_Integral_zz_44_case_4( mix_a_p, mix_r_p, sp2->knot[i], sp2->knot[i+1],
                    sp2->radial_p_coefficient[i][0], sp2->radial_p_coefficient[i][1],
                    sp2->radial_p_coefficient[i][2], sp2->radial_p_coefficient[i][3], dist );
        }
    }
    return Return*FFHA_TO_FFEV_UNIT;
}






































