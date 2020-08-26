/**
 * (c) Author:    Woongkyu Jee, woong.jee.16@ucl.ac.uk, wldndrb1@gmail.com
 * Created:   02.06.2019 ~
 * 	
 * University College London, Department of Chemistry
 **/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<gsl/gsl_eigen.h>
#include"sp_cluster_system.h"
#include<mpi.h>
#include<assert.h>

#include<time.h>

#define MAIN_FALSE -1
#define MAIN_TRUE   1

void MatViewf( FILE* fp, gsl_matrix* m )
{
    for(int i=0;i<m->size1;i++)
    {   for(int j=0;j<m->size2;j++)
            fprintf(fp,"%s%lf\t",gsl_matrix_get(m,i,j)>0.?"+":"",gsl_matrix_get(m,i,j));
        fprintf(fp,"\n");
    }

    return;
}

void VecViewf( FILE* fp, gsl_vector* v )
{
    for(int i=0;i<v->size;i++)
        fprintf(fp,"%lf\t",gsl_vector_get(v,i));
    fprintf(fp,"\n");
    return;
}

void liner_sharp()
{	for(int i=0;i<77;i++)	printf("#");
	printf("\n");
	return;
}
void liner_line()
{	for(int i=0;i<77;i++)	printf("-");
	printf("\n");
	return;
}


int main( int argc, char* argv[] )
{
    FILE* xyz;
    FILE* next_config;
    FILE* cube;

    int if_success = MAIN_FALSE;

    double elapsed_time;
    double opti_elapsed_time;
    /// /// /// /// MPI Variables
    int numtasks, rank;
    /// /// /// /// 

    sp_cluster_system* sp_sys = NULL;
    double grad_x, grad_y, grad_z;

    /// Creat MPI Environment
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

    sp_sys = sp_cluster_system_init();  // init sp_lone pair system single centre

    elapsed_time = -MPI_Wtime();
    if(rank == 0)
    {   
	liner_sharp();
        printf("\n");
        printf(" Sp Lone pair involved Atomistic Model ( S L A M ) \n\n");
        printf(" Author (c)     :   Woongkyu Jee\n");
        printf(" Affiliation    :   Universitiy College London, Department of Chemistry\n");
        printf(" Contact        :   woong.jee.16@ucl.ac.uk / wldndrb1@gmail.com\n");
        printf(" Version        :   slam_v_2.2  ( Lastest Update : 04. 28. 2020 SNAPSHOT )\n");
        printf("\n");
        liner_sharp();
        printf("\n");
        printf(" Initiating SLAM Calculation\n\n");
        printf(" Number of CPUs Requested   : %d\n\n",numtasks);
        printf(" Making Integral Table ( Using Cubic Spline Interpolation Algorithm )\n\n");

    }
    MPI_Barrier(MPI_COMM_WORLD);
    sp_cluster_system_get_spline_integral( sp_sys );
    elapsed_time += MPI_Wtime();
    if(rank == 0)
    {
        printf(" Table Set Up Wtime : %.6lf s\n\n",elapsed_time);
    	liner_sharp();
        /* HERE PUT GENERAL INFO ABOUT THE CALCULATION */
        printf("\n");
        printf(" General Input Info\n");
        printf("\n");	
	liner_line();
	printf(" QM - Born-Mayer potential on sp-lone pair density\n");
	liner_line();
	for(int i=0;i<sp_sys->number_of_qm_interaction_bm;i++)
	{
		if( i < sp_sys->number_of_qm_interaction_bm - 1 )
			printf("%4s%18.6lf%18.6lf\n",sp_sys->interaction_qm_type_bm[i],sp_sys->interaction_qm_AR_bm[i][0],sp_sys->interaction_qm_AR_bm[i][1]);
		else
			printf("%4s%18.6lf%18.6lf%18.6lf%18.6lf\n",sp_sys->interaction_qm_type_bm[i],sp_sys->interaction_qm_AR_bm[i][0],sp_sys->interaction_qm_AR_bm[i][1],
				sp_sys->interaction_qm_AR_bm[i][2],sp_sys->interaction_qm_AR_bm[i][3]);
	}
	liner_line();
	printf(" MM - Buckingham potential is used as default\n");
	liner_line();
	for(int i=0;i<sp_sys->number_of_mm_interaction_buck;i++)
	{
		printf(" %6.5s%6.4s%8.5s%8.4s",sp_sys->interaction_type_buck[i][0],sp_sys->interaction_type_buck[i][1],sp_sys->interaction_type_buck[i][2],sp_sys->interaction_type_buck[i][3]);
		printf("%18.6lf%12.6lf%12.6lf\n",sp_sys->interaction_ARC_buck[i][0],sp_sys->interaction_ARC_buck[i][1],sp_sys->interaction_ARC_buck[i][2]);
	}
	liner_line();
	printf("\n");
	printf(" MM atoms/ions\n");
	printf("\n");
	liner_line();
	printf(" %s%8s%12s%12s%16s%19s\n","Species.","x","y","z","Charge","Spring(k2/k4)");
	liner_line();
	for(int i=0;i<sp_sys->number_of_classic_ion;i++)
	{
		printf("%3s%2s%15.6lf%12.6lf%12.6lf%15.6lf",gsl_vector_get(sp_sys->classic_ion[i].core_position,0),gsl_vector_get(sp_sys->classic_ion[i].core_position,1),gsl_vector_get(sp_sys->classic_ion[i].core_position,2),
			sp_sys->classic_ion[i].atom_name,(sp_sys->classic_ion[i].if_shell==MAIN_TRUE)?"s":"c",sp_sys->classic_ion[i].charge_core);
		if( sp_sys->classic_ion[i].if_shell == MAIN_TRUE )
			printf("%11.4lf/%.4lf\n",sp_sys->classic_ion[i].k2_const,sp_sys->classic_ion[i].k4_const);
		else
			printf("\n");
	}
	liner_line();
	printf("\n");
	printf(" SP atoms/ions\n");
	printf("\n");
	liner_line();
	printf(" %s%8s%12s%12s%16s%14s\n","Species.","x","y","z","q_c/q_sp","Esp");
	liner_line();
	for(int i=0;i<sp_sys->number_of_sp_ion;i++)
	{
		printf(" %s%17.6lf%12.6lf%12.6lf%10.3lf/%.3lf%13.6lf\n",sp_sys->sp_ion[i].atom_name,
			gsl_vector_get(sp_sys->sp_ion[i].core_position,0),gsl_vector_get(sp_sys->sp_ion[i].core_position,1),gsl_vector_get(sp_sys->sp_ion[i].core_position,2),
			sp_sys->sp_ion[i].charge_core,sp_sys->sp_ion[i].charge_shell,sp_sys->sp_ion[i].esp);
	}

	//liner_line();
	//printf("\n");
	//printf(" Position Integral Reference: %.12lf\n",sp_cluster_integrator_get_x_12( sp_sys->sp_ion[0] ));
	//printf("\n");

	liner_line();
	liner_sharp();
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //  END OF INITIATION
    


    int opti_res = MAIN_FALSE;

    opti_elapsed_time = -MPI_Wtime();
    if_success = sp_cluster_system_call_bfgs_algorithm_mpi( sp_sys, 0.20 /*in Angs*/, atoi(argv[2]) , rank,numtasks);
    opti_elapsed_time+= MPI_Wtime();

    if( rank == 0 ) liner_sharp();
    if( if_success == MAIN_FALSE )
    {	if( rank == 0 )
	{	printf("\n");
		printf(" Optimiser Failed ...\n");
		printf("\n");
		printf(" SLAM may need more optimisation steps ...\n");
	}
    }

    if(rank == 0) 
    {
        xyz = fopen(argv[1],"w");
        next_config = fopen("geo.txt.next","w");


        sp_cluster_system_write_xyz(sp_sys,xyz,argv[1]);
        sp_cluster_system_get_next_config_mpi(sp_sys,next_config);
        fclose(xyz);
        fclose(next_config);

        // density test
        // if argv[3] contains the name of the cube file, then write cube file
        if( argv[3] != NULL )
        {   cube = fopen(argv[3],"w");
            sp_cluster_system_get_density( sp_sys, 0.188973, cube );
            fclose(cube);
        }
        // density test
        printf("\n");
        printf(" Computation Wtime :    %.6lf s\n", opti_elapsed_time ); 
	printf("\n");
	if( if_success != MAIN_FALSE )
	{	
		time_t t = time(NULL);
 		struct tm tm = *localtime(&t);
  		printf(" Date		   :	%d-%02d-%02d %02d:%02d:%02d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
		printf("\n");
		printf(" output xyz 	   :	%s\n",argv[1]);
		if( argv[3] != NULL )
			printf(" output cube	   :	%s\n",argv[3]);
		printf(" next config	   :	%s\n","geo.txt.next");
		printf("\n");

	}
	liner_sharp();
    }   

    MPI_Barrier(MPI_COMM_WORLD);

    sp_cluster_system_free_spline_integral( sp_sys );
    sp_cluster_system_detach(sp_sys);

    MPI_Finalize();

    return 0;
}

