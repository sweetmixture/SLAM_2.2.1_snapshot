/**
 * (c) Author:    Woongkyu Jee, woong.jee.16@ucl.ac.uk, wldndrb1@gmail.com
 * Created:   02.06.2019 ~
 * 	
 * University College London, Department of Chemistry
 **/
#ifndef __SP_CLUSTER_SYSTEM__
    #define __SP_CLUSTER_SYSTEM__

#include"sp_cluster_type.h"
#include"sp_cluster_support.h"

// LOAD SP_LONE_PAIR PARAMETERS ON STRUCT 'sp_cluster_system'
// SUB METHOD USED IN 'sp_cluster_system_init()'
sp_cluster_type_sp_ion* sp_cluster_load_sp_parameter();

// LOAD CLASSICAL IONS PARAMETERS ON STRUCT 'sp_cluster_system'
// SUB METHOD USED IN 'sp_cluset_system_init()'
sp_cluster_type_classic_ion* sp_cluster_load_classic_parameter();

// LOAD SP_CLUSTER_SYSTEM
sp_cluster_system* sp_cluster_system_init();
// DETACH SP_CLUSTER_SYSTEM ... memory release
void sp_cluster_system_detach( sp_cluster_system* ptr );

void sp_cluster_system_get_spline_integral( sp_cluster_system* sp_sys );
void sp_cluster_system_free_spline_integral( sp_cluster_system* sp_sys );


// testing ... 2019 2 14
void sp_cluster_system_get_h_matrix_mpi( sp_cluster_system* sp_sys, int rank, int numtasks );

void sp_cluster_system_get_eigensystem( sp_cluster_system* sp_sys );

int sp_cluster_system_set_is_scf_done( sp_cluster_system* sp_sys, int t_or_f );

int sp_cluster_system_scf_mpi( sp_cluster_system* sp_sys, const int rank, const int numtasks );


// force calc ( on sp core )
void sp_cluster_system_get_force( sp_cluster_system* sp_sys );
void sp_cluster_system_get_force_mpi( sp_cluster_system* sp_sys, int rank, int numtasks );

// classic energy
void sp_cluster_system_get_classic_energy( sp_cluster_system* sp_sys );

// force classic calc pure coulom w.r.t classic cores vs sp core
void sp_cluster_system_get_classic_force( sp_cluster_system* sp_sys );

// Get total energy after SCF & Classic Energy are estimated
double sp_cluster_system_get_cluster_energy( sp_cluster_system* sp_sys );

void sp_cluster_system_get_density( sp_cluster_system* sp_sys, const double voxel_in, FILE* fp );

void sp_cluster_system_get_next_config_mpi( sp_cluster_system* sp_sys, FILE* fp );


void sp_cluster_system_bfgs_support_get_alpha_mpi_ls( sp_cluster_system* sp_sys, gsl_vector* p_k, double stepmx, int rank, int numtasks);
void sp_cluster_system_bfgs_support_load_x( sp_cluster_system* sp_sys, gsl_vector* x );
void sp_cluster_system_bfgs_support_load_g( sp_cluster_system* sp_sys, gsl_vector* g );
int sp_cluster_system_call_bfgs_algorithm_mpi( sp_cluster_system* sp_sys, double stepmx /* maximum stepsize */, const int cyclemx, int rank, int numtasks );







void sp_cluster_system_write_xyz( sp_cluster_system* sp_sys, FILE* fp, char* fn );


// supportive methods
double sp_cluster_system_get_gnorm( sp_cluster_system* sp_sys );

#endif
