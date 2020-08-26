/**
 * (c) Author:    Woongkyu Jee, woong.jee.16@ucl.ac.uk, wldndrb1@gmail.com
 * Created:   02.06.2019 ~
 * 	
 * University College London, Department of Chemistry
 **/
#ifndef __SP_CLUSTER_INTEGRATOR__
    #define __SP_CLUSTER_INTEGRATOR__

#include"sp_cluster_integral_lib/CH_Integrals.h"
//#include"sp_cluster_integral_lib/SH_Integrals.h"
#include"sp_cluster_integral_lib/BM_SH_Integrals.h"
#include"sp_cluster_integral_lib/CDH_Integrals.h"
//#include"sp_cluster_integral_lib/SDH_Integrals.h"
#include"sp_cluster_integral_lib/BM_SDH_Integrals.h"
#include"sp_cluster_integral_lib/CDDH_Integrals.h"
#include"sp_cluster_integral_lib/BM_SDDH_Integrals.h"
#include"sp_cluster_integral_lib/CDDDH_Integrals.h"
#include"sp_cluster_integral_lib/BM_SDDDH_Integrals.h"
#include"sp_cluster_type.h"

#define SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE 145	// tight integral setup
//#define SP_CLUSTER_INTEGRATOR_INTEGRAL_TOLERANCE 125
#define EV_UNIT 14.39964390675221758120                 // mks -> eV

/*
#define TO_BOHR_RADII 0.52917721067                     // 1bohr == 0.529...Angs
#define HA_TO_EV_UNIT 27.21138582                       // q1q2/r (atomic unit /Hartree(Ha)) -> eV
#define FHA_TO_FEV_UNIT 51.42208629083232               // Force Ha/Bohr -> eV/Angs
#define FFHA_TO_FFEV_UNIT 97.1553005                    // Ha/Bohr^2 -> eV/Angs^2       ... second derivatives
*/
#define _RADIAL_FUNCTION_SCALE_ 1.0

#define TO_BOHR_RADII (0.52917721067*_RADIAL_FUNCTION_SCALE_)
#define HA_TO_EV_UNIT (EV_UNIT/TO_BOHR_RADII)
#define FHA_TO_FEV_UNIT (EV_UNIT/TO_BOHR_RADII/TO_BOHR_RADII)
#define FFHA_TO_FFEV_UNIT (EV_UNIT/TO_BOHR_RADII/TO_BOHR_RADII/TO_BOHR_RADII)



#define IN
#define OUT

#define SP_INTEGRAL_SHORT_RANGE_CUTOFF 9.50
#define SP_INTEGRAL_SHORT_RANGE_START 1.00
#define SP_INTEGRAL_SHORT_RANGE_STEP  1.0123
// Beware !! These Are In Angs Unit

#define SP_INTEGRAL_TRUE 1
#define SP_INTEGRAL_FALSE -1

int sp_cluster_integrator_lut_b_search( double dist, int knot_stride, const double* integral_knot );

// Long Range Integral Methods

// SS
double sp_cluster_integrator_get_ch_11_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );
// SZ
double sp_cluster_integrator_get_ch_14_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );
#define sp_cluster_integrator_get_ch_41_element sp_cluster_integrator_get_ch_14_elemetn
// XX or YY
double sp_cluster_integrator_get_ch_2233_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );
// ZZ
double sp_cluster_integrator_get_ch_44_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );


// Short Range Integral Methods

// SS
double sp_cluster_integrator_get_sh_11_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );
// SZ
double sp_cluster_integrator_get_sh_14_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );
#define sp_cluster_integrator_get_sh_41_element sp_cluster_integrator_get_sh_14_element
// XX or YY
double sp_cluster_integrator_get_sh_2233_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );
// ZZ
double sp_cluster_integrator_get_sh_44_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );


// 1st Derivative methods

// LongRange

// CDX SX == CDY SY
double sp_cluster_integrator_get_ch_x_12_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );
#define sp_cluster_integrator_get_ch_x_21_element   sp_cluster_integrator_get_ch_x_12_element
#define sp_cluster_integrator_get_ch_y_13_element   sp_cluster_integrator_get_ch_x_12_element
#define sp_cluster_integrator_get_ch_y_31_element   sp_cluster_integrator_get_ch_x_12_element

// CDX XZ == CDY YZ
double sp_cluster_integrator_get_ch_x_24_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );
#define sp_cluster_integrator_get_ch_x_42_element   sp_cluster_integrator_get_ch_x_24_element
#define sp_cluster_integrator_get_ch_y_34_element   sp_cluster_integrator_get_ch_x_24_element
#define sp_cluster_integrator_get_ch_y_43_element   sp_cluster_integrator_get_ch_x_24_element

// CDZ SS
double sp_cluster_integrator_get_ch_z_11_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );
// CDZ SZ
double sp_cluster_integrator_get_ch_z_14_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );
// CDZ XX YY
double sp_cluster_integrator_get_ch_z_2233_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );
// CDZ ZZ
double sp_cluster_integrator_get_ch_z_44_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );

// ShortRange

// SDX SX == SDY SY
double sp_cluster_integrator_get_sh_x_12_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );
#define sp_cluster_integrator_get_sh_x_21_element   sp_cluster_integrator_get_sh_x_12_element
#define sp_cluster_integrator_get_sh_y_13_element   sp_cluster_integrator_get_sh_x_12_element
#define sp_cluster_integrator_get_sh_y_31_element   sp_cluster_integrator_get_sh_x_12_element

// SDX XZ == SDY YZ
double sp_cluster_integrator_get_sh_x_24_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );
#define sp_cluster_integrator_get_sh_x_42_element   sp_cluster_integrator_get_sh_x_24_element
#define sp_cluster_integrator_get_sh_y_34_element   sp_cluster_integrator_get_sh_x_24_element
#define sp_cluster_integrator_get_sh_y_43_element   sp_cluster_integrator_get_sh_x_24_element

// SDZ SS
double sp_cluster_integrator_get_sh_z_11_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );
// SDZ SZ
double sp_cluster_integrator_get_sh_z_14_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );
// SDZ XX YY
double sp_cluster_integrator_get_sh_z_2233_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );
// SDZ ZZ
double sp_cluster_integrator_get_sh_z_44_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );


// 2nd Derivative methods

// LongRange

// CDDXX SS == CDDYY SS
double sp_cluster_integrator_get_ch_xx_11_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );
// CDDXX SZ == CDDYY SZ
double sp_cluster_integrator_get_ch_xx_14_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );
// CDDXX XX == CDDYY YY
double sp_cluster_integrator_get_ch_xx_22_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );
// CDDXX YY == CDDYY XX
double sp_cluster_integrator_get_ch_xx_33_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );
// CDDXX ZZ == CDDYY ZZ
double sp_cluster_integrator_get_ch_xx_44_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );


// CDDXY XY
double sp_cluster_integrator_get_ch_xy_23_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );


// CDDXZ SX == CDDYZ SY
double sp_cluster_integrator_get_ch_xz_12_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );
// CDDXZ XZ == CDDYZ YZ
double sp_cluster_integrator_get_ch_xz_24_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );


// CDDZZ SS
double sp_cluster_integrator_get_ch_zz_11_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );
// CDDZZ SZ
double sp_cluster_integrator_get_ch_zz_14_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );
// CDDZZ XXYY
double sp_cluster_integrator_get_ch_zz_2233_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );
// CDDZZ ZZ
double sp_cluster_integrator_get_ch_zz_44_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );


// ShortRange

// SDDXX SS == SDDYY SS
double sp_cluster_integrator_get_sh_xx_11_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );
// SDDXX SZ == SDDYY SZ
double sp_cluster_integrator_get_sh_xx_14_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );
// SDDXX XX == SDDYY YY
double sp_cluster_integrator_get_sh_xx_22_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );
// SDDXX YY == SDDYY XX
double sp_cluster_integrator_get_sh_xx_33_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );
// SDDXX ZZ == SDDYY ZZ
double sp_cluster_integrator_get_sh_xx_44_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );


// SDDXY XY == SDDXY XY
double sp_cluster_integrator_get_sh_xy_23_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );


// SDDXZ SX == SDDYZ SY
double sp_cluster_integrator_get_sh_xz_12_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );
// SDDXZ XZ == SDDYZ YZ
double sp_cluster_integrator_get_sh_xz_24_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );


// SDDZZ SS
double sp_cluster_integrator_get_sh_zz_11_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );
// SDDZZ SZ
double sp_cluster_integrator_get_sh_zz_14_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );
// SDDZZ XXYY
double sp_cluster_integrator_get_sh_zz_2233_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );
// SDDZZ zz
double sp_cluster_integrator_get_sh_zz_44_element( sp_cluster_type_sp_ion* sp, sp_cluster_type_classic_ion* ion );

















// Direction Integral ... x, y and z
double sp_cluster_integrator_get_x_12( sp_cluster_type_sp_ion* sp );
double sp_cluster_integrator_get_y_13( sp_cluster_type_sp_ion* sp );
double sp_cluster_integrator_get_z_14( sp_cluster_type_sp_ion* sp );


// N Centre Integral Functions
//
// Naming convention ...spsp...

// For a Given State 'sp1' calculating 'sp2' Integrals !!!

// Mono-pole related
//
// Coulomb
double sp_cluster_spsp_mono_integrator_get_ch_11_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );
double sp_cluster_spsp_mono_integrator_get_ch_14_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );
double sp_cluster_spsp_mono_integrator_get_ch_2233_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );
double sp_cluster_spsp_mono_integrator_get_ch_44_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );
// BM
double sp_cluster_spsp_mono_integrator_get_sh_11_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );
double sp_cluster_spsp_mono_integrator_get_sh_14_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );
double sp_cluster_spsp_mono_integrator_get_sh_2233_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );
double sp_cluster_spsp_mono_integrator_get_sh_44_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );
//
// Mono-pole relaated Done

// Di-pole related
//
// Coulomb
double sp_cluster_spsp_di_integrator_get_ch_x_12_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );
double sp_cluster_spsp_di_integrator_get_ch_x_24_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );  // wrt x
// ch_x_12 == ch_y_13 & ch_x_24 == ch_y_34
double sp_cluster_spsp_di_integrator_get_ch_z_11_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );
double sp_cluster_spsp_di_integrator_get_ch_z_14_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );
double sp_cluster_spsp_di_integrator_get_ch_z_2233_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );
double sp_cluster_spsp_di_integrator_get_ch_z_44_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );
// BM
double sp_cluster_spsp_di_integrator_get_sh_x_12_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );
double sp_cluster_spsp_di_integrator_get_sh_x_24_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );
// sh_x_12 == sh_y_13 & sh_x_24 == sh_y_34
double sp_cluster_spsp_di_integrator_get_sh_z_11_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );
double sp_cluster_spsp_di_integrator_get_sh_z_14_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );
double sp_cluster_spsp_di_integrator_get_sh_z_2233_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );
double sp_cluster_spsp_di_integrator_get_sh_z_44_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );
//
// Di-pole related Done

// Quadru-pole related
//
// Coulomb
double sp_cluster_spsp_quadru_integrator_get_ch_xx_11_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );
double sp_cluster_spsp_quadru_integrator_get_ch_xx_14_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );
double sp_cluster_spsp_quadru_integrator_get_ch_xx_22_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );
double sp_cluster_spsp_quadru_integrator_get_ch_xx_33_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );
double sp_cluster_spsp_quadru_integrator_get_ch_xx_44_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );      //xx

double sp_cluster_spsp_quadru_integrator_get_ch_xy_23_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );      //xy

double sp_cluster_spsp_quadru_integrator_get_ch_xz_12_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );      //xz
double sp_cluster_spsp_quadru_integrator_get_ch_xz_24_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );

double sp_cluster_spsp_quadru_integrator_get_ch_zz_11_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );      //zz
double sp_cluster_spsp_quadru_integrator_get_ch_zz_14_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );      
double sp_cluster_spsp_quadru_integrator_get_ch_zz_2233_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );      
double sp_cluster_spsp_quadru_integrator_get_ch_zz_44_element( sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );      
// BM
double sp_cluster_spsp_quadru_integrator_get_sh_xx_11_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );
double sp_cluster_spsp_quadru_integrator_get_sh_xx_14_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );
double sp_cluster_spsp_quadru_integrator_get_sh_xx_22_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );
double sp_cluster_spsp_quadru_integrator_get_sh_xx_33_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );
double sp_cluster_spsp_quadru_integrator_get_sh_xx_44_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );      //xx

double sp_cluster_spsp_quadru_integrator_get_sh_xy_23_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );      //xy

double sp_cluster_spsp_quadru_integrator_get_sh_xz_12_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );      //xz
double sp_cluster_spsp_quadru_integrator_get_sh_xz_24_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );

double sp_cluster_spsp_quadru_integrator_get_sh_zz_11_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );      //zz
double sp_cluster_spsp_quadru_integrator_get_sh_zz_14_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );      
double sp_cluster_spsp_quadru_integrator_get_sh_zz_2233_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );      
double sp_cluster_spsp_quadru_integrator_get_sh_zz_44_element( sp_cluster_system* sp_sys, sp_cluster_type_sp_ion* sp1, OUT sp_cluster_type_sp_ion* sp2 );      
// Quadru End











#endif
