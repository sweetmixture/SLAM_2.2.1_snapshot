/**
 * (c) Author:    Woongkyu Jee, woong.jee.16@ucl.ac.uk, wldndrb1@gmail.com
 * Created:   02.06.2019 ~
 * 	
 * University College London, Department of Chemistry
 **/
#ifndef __SP_CLUSTER_TYPE__
    #define __SP_CLUSTER_TYPE__

#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include"sp_cluster_atom_type.h"

enum _Element{H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, Cr, Mn, Fe, Co, Ni,
	Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe, Cs,
	Ba, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb,
	Bi, Po, At, Rn};

typedef struct _sp_cluster_type_sp_ion
{	
	/**/
	char atom_name[3];
	int _fix_flag;
	enum _Element atom_number;
	/**/

    double charge_core;
    double charge_shell;        // charge of sp-lone pair electrons
    double esp;                 // energy difference between s and p states

    // initiall noting is saved when the system is initialised
    gsl_vector* eigen_value;    // eigen_value  type gsl_vector
    gsl_matrix* eigen_vector;   // eigen_vector type gsl_matrix
    
    gsl_matrix* h_matrix;       // Hamiltonian matrix
    gsl_matrix** dh_matrix;     // First Derivatives of Hamiltonian matrix elements
    gsl_matrix*** ddh_matrix;   // Second Derivatives
    gsl_matrix**** dddh_matrix; // Third Derivatives

    // coordinate of sp-lone pair core
    gsl_vector* core_position;

    // saving spline information of sp-lone model basis set
    int number_of_knot;
    double* knot;
    double** radial_s_coefficient;
    double** radial_p_coefficient;

    // short range parameter
    double short_range_a_s; double short_range_r_s; // ifthere is core buckingham then use this parameters !!
    double short_range_a_p; double short_range_r_p;

    // physical properties
    double force[3];            // force acting on sp core by all the ions vs sp_electron density .. to be derivative sign must be inversed
    double classic_force[3];    // derivative w.r.t sp core by all the ions ... to be force sign must be inversed

    /// OLD ENERGY VAL ... NEED TO BE CHANGED ! 07252019
    double classic_energy; 

    
    // Sp Elec Props
    double sp_energy;           // Energy by sp-els interactions ... / vs point charges/ vs sp_ions / vs sp_onsite
    double force_by_sp_elec[3];
    
    double elec_force_by_sp[3];
    double elec_force_by_ion[3];


    // Classical Props
    double classic_energy_by_sp_core;      // Coulomb Energy by sp core vs external ions //+ maybe + sp_core vs external ion Born-Mayer term
    double classic_energy_by_ion_core;
    double force_by_sp_core[3];
    double force_by_ion_core[3];
    // derivative evec wrt xyz ... only total
    

    double cent_short_range_a;
    double cent_short_range_r;
    double cent_short_range_vdw;
    

}sp_cluster_type_sp_ion;

typedef struct _sp_cluster_type_classic_ion
{
	/**/
	char atom_name[3];
	int _fix_flag;
	enum _Element atom_number;
	/**/


    double charge_core;             // charge of classical ion
    gsl_vector* core_position;      // coordinate of classical ion core
    // short range parameter
    // UPDATE INPUT FILE + INPUT LOAD METHOD IN 'sp_cluster_system_load ...'
    double short_range_a;   double short_range_r;   double short_range_c;

    double force_by_sp_shell[3];            // force on the ion core by sp_electron density               .. to be derivative, sign must be inversed
    double classic_force[3];                // classical derivative of ion core w.r.t  all the other ions .. to be force, sign must be inversed

    // OLD ENERGY VAL ... NEED TO BE CHANGED ! 07252019
    double classic_energy;

    // Sp Elec Props
    double elec_force_by_sp[3];

    // Classical Props
    double classic_energy_by_sp_core;       // Buckingham interaction on sp-core vs external ion 'shell or core'
    double classic_energy_by_ion_core;
    double force_by_sp_core[3];             // classical derivative of ion core w.r.t  sp core            .. to be force, signe must be inversed
	double force_by_ion_core[3];			// Buckingham interaction by sp_core (Coulomb + Born-Mayer)

    // CORE PART SUMMARY ... SHELL IS NOT ON TOP OF IT YET !!
    double total_classic_energy;
    double total_force[3];                  // total derivative classic_force + force_by_sp_core - force_by_sp_shell    // THIS IS FORCE ACTING ON 'Core'
    
    ////    ////    ////    ABOVE ALL ABOUT CORES

    ////    ////    WHEN ADDING SHELL FEATURES BELOW, KEEP IT HIGHLY INDEPENDENT

    /// /// /// /// /// /// /// /// /// /// /// /// /// /// /// /// /// /// /// /// ///
    /* SHELL PROPS */
    int if_shell;                        // Flag ... if this species has shell -> then set 'True'
    int if_has_shell;			 // Flag ... if this is core and has shell ?
    //double charge_shell;                    // CHARGE OF SHELL PART
    double k2_const;                        // k2 SPRING CONST ... UNIT eV/Angs/Angs
    double k4_const;
    //double short_range_a_shell;   
    //double short_range_r_shell;             // Shell BuckingHam Params  
    //double short_range_c_shell; 
    double elastic_energy;                  // shell elastic energy
    double elastic_force_shell[3];          // elastic force on shell
    double elastic_force_core[3];           // elastic force on core ... this is 'counter force'
    //gsl_vector* shell_position;             // SHELL_POSITION

    // classic energies & forces on the 'shell'
    //double classic_energy_by_sp_core_shell;     // sp-core vs shell
    //double classic_energy_by_ion_core_shell;    // core    vs shell
    //double classic_energy_by_ion_shell_shell;   // shell   vs shell
    //double force_by_sp_core_shell[3];
    //double force_by_ion_core_shell[3];
    //double force_by_ion_shell_shell[3];

    // sp elec density relates ...
    //double elec_force_by_sp_shell[3];           // note that the corresponding elec energy is already counted by LP elec density

    // summary
    //double total_classic_energy_shell;          // total energy of the shell by MM species
    //double total_force_shell[3];                // total force acting on the shell

    // RELATED FEATURES ARE NOT IMPLEMENTED YET ! 28.02.2020
    
    /* LIST OF THE UPCOMING FEATURES TO ADD
     *
     *		# CLASSICAL MODULES
     *	
     *	1. SHELL - SHELL INTERACTION
     *	
     *		- ENERGY, FORCE
     *
     *	2. SHELL - CORE  INTERACTION
     *
     *		- ENERGY, FORCE
     *
     *	3. SELF ENERGY
     *
     *		- SHELL-CORE ENERGY
     *		- SHELL-CORE FORCE
     *
     *		#	4. CORE  - CORE  INTERACTION - DONE
     *
     * 	///	///	///	///	///	///
     *
     * 		# QM MODULES
     *
     *	1. SHELL - SP INTERACTION
     *
     * 		- ENERGY, EVEC DERIVATIVES, FORCE
     *
     *		#	2. CORE - SP INTERACTION - DONE
     *	
     *
     *	if there is a shell... then sp density - core might not include short_range interaction
     *
     *	the short range only in sp density - shell.
     *
     */

    /// /// /// /// /// /// /// /// /// /// /// /// /// /// /// /// /// /// /// /// ///

}sp_cluster_type_classic_ion;
    

typedef struct _sp_cluster_system_general_info
{

    int if_first_scf_trial;

    int number_of_classic_ion;                  // saves number of classical ion
    int number_of_sp_ion;                       // saves number of sp ion
    int number_of_species;                      // saves number of species

    sp_cluster_type_classic_ion* classic_ion;   // saves classic ion info
    sp_cluster_type_sp_ion* sp_ion;                 // saves sp-lone pair ion info

    double gnorm;


    // SCF WORKSPACE
    int scf_cnt;        // count number of scf cycle performed.
    int is_scf_done;    // Flag for checking if it is scf converged e.g., if '1' converged; else if '0' not converged.

    gsl_matrix** scf_h_matrix_vs_classic_ion;
    gsl_matrix** scf_h_matrix_vs_sp_ion_monopole;   // This term includes spe - spe Coulomb + Short .. + spe - spc Coulomb
    gsl_matrix** scf_h_matrix_vs_sp_ion_onsite;

    gsl_matrix** scf_h_matrix_vs_sp_core;       // This term includes spe - spc Coulomb

    gsl_matrix*** scf_dh_x_matrix_vs_sp_ion_dipole;
    gsl_matrix*** scf_dh_y_matrix_vs_sp_ion_dipole;
    gsl_matrix*** scf_dh_z_matrix_vs_sp_ion_dipole;

    // FORCE CALCULATION WORK SPACE
    
    /// vs cla ion
    
    gsl_matrix*** cla_dh_matrix_workspace_x;  gsl_matrix*** cla_dh_matrix_workspace_y;  gsl_matrix*** cla_dh_matrix_workspace_z;
    gsl_matrix*** cla_dh_matrix_x;  gsl_matrix*** cla_dh_matrix_y;  gsl_matrix*** cla_dh_matrix_z;

    /// vs sp ion

    // sp-sp mono + core Buffer
    gsl_matrix*** dh_matrix_workspace_x; gsl_matrix*** dh_matrix_workspace_y; gsl_matrix*** dh_matrix_workspace_z;
    gsl_matrix*** dh_matrix_x; gsl_matrix*** dh_matrix_y; gsl_matrix*** dh_matrix_z;
    // sp-sp Dipolar Buffer
    gsl_matrix*** ddh_matrix_workspace_xx;  gsl_matrix*** ddh_matrix_workspace_xy;  gsl_matrix*** ddh_matrix_workspace_xz;  
    gsl_matrix*** ddh_matrix_workspace_yx;  gsl_matrix*** ddh_matrix_workspace_yy;  gsl_matrix*** ddh_matrix_workspace_yz;
    gsl_matrix*** ddh_matrix_workspace_zx;  gsl_matrix*** ddh_matrix_workspace_zy;  gsl_matrix*** ddh_matrix_workspace_zz; 
    gsl_matrix*** ddh_matrix_xx;    gsl_matrix*** ddh_matrix_xy;    gsl_matrix*** ddh_matrix_xz;    
    gsl_matrix*** ddh_matrix_yx;    gsl_matrix*** ddh_matrix_yy;    gsl_matrix*** ddh_matrix_yz;
    gsl_matrix*** ddh_matrix_zx;    gsl_matrix*** ddh_matrix_zy;    gsl_matrix*** ddh_matrix_zz;


    // SPLINE INTEGRALS
    int integral_lut;
    int knot_stride;

    double* integral_knot;

    // ShortRange ... Depending on the number of species ... e.g. if Ba, O and Sn -> need three 3C2 types
    
    // LUT DATATPYES
    double*** integral_vs_cla_s_ss;
    double*** integral_vs_cla_s_sz;
    double*** integral_vs_cla_s_xxyy;
    double*** integral_vs_cla_s_zz;

    // ShortRange - First Derivative
    double*** integral_vs_cla_s_x_sx;
    double*** integral_vs_cla_s_x_xz;

    double*** integral_vs_cla_s_z_ss;
    double*** integral_vs_cla_s_z_sz;
    double*** integral_vs_cla_s_z_xxyy;
    double*** integral_vs_cla_s_z_zz;

    // ShortRange - Second Derivative ... only doing with sp-sp interaction

    double** integral_vs_sp_s_ss;
    double** integral_vs_sp_s_sz;
    double** integral_vs_sp_s_xxyy;
    double** integral_vs_sp_s_zz;       //mono

    double** integral_vs_sp_s_x_sx;
    double** integral_vs_sp_s_x_xz;

    double** integral_vs_sp_s_z_ss;
    double** integral_vs_sp_s_z_sz;
    double** integral_vs_sp_s_z_xxyy;
    double** integral_vs_sp_s_z_zz;     //di

    double** integral_vs_sp_s_xx_ss;
    double** integral_vs_sp_s_xx_sz;
    double** integral_vs_sp_s_xx_xx;
    double** integral_vs_sp_s_xx_yy;
    double** integral_vs_sp_s_xx_zz;

    double** integral_vs_sp_s_xy_xy;

    double** integral_vs_sp_s_xz_sx;
    double** integral_vs_sp_s_xz_xz;

    double** integral_vs_sp_s_zz_ss;
    double** integral_vs_sp_s_zz_sz;
    double** integral_vs_sp_s_zz_xxyy;
    double** integral_vs_sp_s_zz_zz;          //quadru
    // ShortRange End


    // eigenvector derivatives

    // w.r.t sp
    double**** deriv_evec_sp;
    // w.r.t cla
    double**** deriv_evec_cla;
    // convention
    /* 
     * [differentiate with this parameter][differentiation target 'm'th eigenvector]['m'th eigen vector, e.g. 0 -> c_0][direction x,y and z ... 012];
     *
     */

	/* TYPE OF QM INTERACTIONS */
	int number_of_qm_interaction_bm;
	char interaction_qm_type_bm[100][16];
	double interaction_qm_AR_bm[100][4];

	/* TYPE OF BUCKINGHAM INTERACTIONS */
	int number_of_mm_interaction_buck;
	char interaction_type_buck[100][4][16];			
	// e.g., 'interaction_type[0][..][..]' keeps '0'th interaction b/t type 'interaction_type[0][0][..] <-> interaction_type[0][1][..]'
	// 	  interaction_type[..][..][16]' keeps type of interacting body, for instance, 'SHELL' <-> 'CORE' ...	
	double interaction_ARC_buck[100][3];	
	// keeps interaction parameters ..
	// A RHO C are saved by same orders with 'char interaction_type[..][..][..]' specified right above
	// Again, here ARC are parameters in Buckingham potential

	/* RECOMMENDED FORMAT OF INPUT FILE 
	 *
	 * QM_POT #
	 * .
	 * .
	 * .
	 * TYPE #
	 * RIM 	 TYPE 	CHARGE_C
	 * SHELL TYPE 	CHARGE_C CHARGE_S SPRING
	 * SP_TYPE 	CHARGE_C CHARGE_S ESP
	 * .
	 * .
	 * MM_POT #
	 * SHELL/CORE TYPE SHELL/CORE TYPE	A	RHO	C
	 * SHELL/CORE TYPE SHELL/CORE TYPE	A	RHO	C
	 * SP_TYPE(CORE) TYPE SEHLL/CORE 	A	RHO	C
	 * SP_TYPE(CORE) TYPE SHELL/CORE	A	RHO	C
	 * .
	 * .
	 */
/*
	char element_name[100][4] = {"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", 
	"Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "Cr", "Mn", "Fe", 
	"Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", 
	"Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", 
	"La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy","Ho", "Er", "Tm", "Yb", "Lu", 
	"Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn"};
*/
}sp_cluster_system;

#endif
