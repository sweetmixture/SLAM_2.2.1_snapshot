/**
 * (c) Author:    Woongkyu Jee, woong.jee.16@ucl.ac.uk, wldndrb1@gmail.com
 * Created:   02.06.2019 ~
 * 	
 * University College London, Department of Chemistry
 **/
#ifndef __SP_CLUSTER_SUPPORT__
    #define __SP_CLUSTER_SUPPORT__

// Get the lowest energy index
int sp_cluster_support_get_lowest_state( gsl_vector* v );


// TENSOR VISUALISATION MODULES + INTERNAL PARAMETER VIEWERS + FIND TRANSFORMATION

// Matrix Viewer
void sp_cluster_support_matrix_view( const gsl_matrix* m );
void sp_cluster_support_matrix_view_f( FILE* fp, const gsl_matrix* m );

// Vector Viwer
void sp_cluster_support_vector_view( const gsl_vector* v );
void sp_cluster_support_vector_view_f( FILE* fp, const gsl_vector* v );

// Get Transformation Matrix for 'Global (Original)' -> 'Local (vector on z' axis)' Symmetry
gsl_matrix* sp_cluster_support_get_transformation_matrix( const gsl_vector* v );

// kronecker_delta
double sp_cluster_support_kronecker_delta( int a, int b );

// get norm of a vector 3D
double sp_cluster_support_get_norm( double v1, double v2, double v3 );

double** sp_cluster_support_get_spline( const double** data, const int knot_stride );


void sp_cluster_support_print_xyz( void* sp_sys_void /* in "sp_cluster_support.c" cast the type with "sp_cluster_system*" */, 
		const double cur_energy, const int rank, const int numtasks );

#endif
