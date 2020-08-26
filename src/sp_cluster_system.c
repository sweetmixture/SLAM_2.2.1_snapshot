/**
 * (c) Author:    Woongkyu Jee, woong.jee.16@ucl.ac.uk, wldndrb1@gmail.com
 * Created:   02.06.2019 ~
 * 	
 * University College London, Department of Chemistry
 **/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include<gsl/gsl_eigen.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_math.h>
#include<mpi.h>

#include"sp_cluster_support.h"
#include"sp_cluster_integrator.h"

#define SP_CLUSTER_SUPPORT_KNOT_PATH "sp_cluster_parameter_src/sp_cluster_knot.txt"
#define SP_CLUSTER_SUPPORT_RADIAL_S_PATH "sp_cluster_parameter_src/sp_cluster_radial_s.txt"
#define SP_CLUSTER_SUPPORT_RADIAL_P_PATH "sp_cluster_parameter_src/sp_cluster_radial_p.txt"
#define SP_CLUSTER_SUPPORT_GEO_PATH "geo.txt"

#define SP_CLUSTER_SUPPORT_SPECIES_PATH "sp_cluster_species.txt"

#define vIN
#define vOUT
#define vINOUT

#define MIN(a,b)        ((a)>=(b)?(b):(a))
#define ARI_MEAN(a,b)   (((a)+(b))/2.)
#define GEO_MEAN(a,b)   sqrt((a)*(b))

#define MAX_SCF_CYCLE 1200

#define SP_SYSTEM_TRUE 1
#define SP_SYSTEM_FALSE -1

#define SP_SYSTEM_SCF_TOL          10E-12
#define SP_SYSTEM_EVEC_TOL         10E-12
//#define SP_SYSTEM_GNORM_TOL        0.0000099999
#define SP_SYSTEM_GNORM_TOL        0.0000499999


void liner()
{   printf("\n");
    for(int f=0;f<77;f++)  printf("#");
    printf("\n");
    return;
}


//#define DEBUG

sp_cluster_type_sp_ion* sp_cluster_load_sp_parameter()
{
    sp_cluster_type_sp_ion* pReturn = NULL;
    int line_number;    char dummy[128];    double tmp;
    int NumberOfClassicIons;    int NumberOfSpIons;
    FILE* fp_knot = NULL; FILE* fp_radial_s = NULL; FILE* fp_radial_p = NULL;
    FILE* fp_geo = NULL;
 
	char _if_fix[8];	int fp_cur;
   
    double tmp_knot;
    double tmp_s_coefficient[4];
    double tmp_p_coefficient[4];

    fp_knot = fopen(SP_CLUSTER_SUPPORT_KNOT_PATH,"r");
    fp_radial_s = fopen(SP_CLUSTER_SUPPORT_RADIAL_S_PATH,"r");
    fp_radial_p = fopen(SP_CLUSTER_SUPPORT_RADIAL_P_PATH,"r");
    fp_geo = fopen(SP_CLUSTER_SUPPORT_GEO_PATH,"r");

    fgets(dummy,sizeof(dummy),fp_geo);      // read out first line in "geo.txt"
    fscanf(fp_geo,"%d",&NumberOfClassicIons);   // read Number Of Classic ions
    fscanf(fp_geo,"%d",&NumberOfSpIons);        // read Number Of sp - ions

    pReturn = (sp_cluster_type_sp_ion*)malloc(NumberOfSpIons*sizeof(sp_cluster_type_sp_ion));
    //
    // Read Radial Distribution Functions
    fscanf(fp_knot,"%d",&line_number);      // # of knots
    for(int i=0;i<NumberOfSpIons;i++)
        pReturn[i].number_of_knot = line_number;    // read # of knots in ith sp - ion
    fscanf(fp_radial_s,"%*d");  fscanf(fp_radial_p,"%*d");
    for(int i=0;i<NumberOfSpIons;i++)
    {
        pReturn[i].knot = (double*)calloc(line_number,sizeof(double));        // knot type: (double*)
        pReturn[i].radial_s_coefficient = (double**)malloc((line_number-1)*sizeof(double));
        pReturn[i].radial_p_coefficient = (double**)malloc((line_number-1)*sizeof(double));   // type (double**)

        for(int j=0;j<line_number-1;j++)
        {   pReturn[i].radial_s_coefficient[j] = (double*)calloc(4,sizeof(double));
            pReturn[i].radial_p_coefficient[j] = (double*)calloc(4,sizeof(double));   }

    }// Memory allocation for sp_ion radial functions
    for(int i=0;i<line_number;i++)
    {   
        fscanf(fp_knot,"%lf", &tmp_knot);
        if( i < line_number-1 )
        {   fscanf(fp_radial_s,"%lf",&tmp_s_coefficient[0]);
            fscanf(fp_radial_s,"%lf",&tmp_s_coefficient[1]);
            fscanf(fp_radial_s,"%lf",&tmp_s_coefficient[2]);
            fscanf(fp_radial_s,"%lf",&tmp_s_coefficient[3]);

            fscanf(fp_radial_p,"%lf",&tmp_p_coefficient[0]);
            fscanf(fp_radial_p,"%lf",&tmp_p_coefficient[1]);
            fscanf(fp_radial_p,"%lf",&tmp_p_coefficient[2]);
            fscanf(fp_radial_p,"%lf",&tmp_p_coefficient[3]);
        }

        for(int j=0;j<NumberOfSpIons;j++)
        {
            pReturn[j].knot[i] = tmp_knot;      // save knot

            if( i < line_number-1 )
            {   pReturn[j].radial_s_coefficient[i][0] = tmp_s_coefficient[0];   pReturn[j].radial_s_coefficient[i][1] = tmp_s_coefficient[1];
                pReturn[j].radial_s_coefficient[i][2] = tmp_s_coefficient[2];   pReturn[j].radial_s_coefficient[i][3] = tmp_s_coefficient[3];   // save s function coefficients

                pReturn[j].radial_p_coefficient[i][0] = tmp_p_coefficient[0];   pReturn[j].radial_p_coefficient[i][1] = tmp_p_coefficient[1];
                pReturn[j].radial_p_coefficient[i][2] = tmp_p_coefficient[2];   pReturn[j].radial_p_coefficient[i][3] = tmp_p_coefficient[3];   // save p function coefficients
            }
        }
    }
    fclose(fp_knot);    fclose(fp_radial_s);    fclose(fp_radial_p);
    // radial input done
    
	// READ "geo.txt"	!!!	!!!	!!!	!!!	!!!

    // Read Dummy Classic ion Info
    for(int i=0;i<NumberOfClassicIons;i++)
    { 
	fscanf(fp_geo,"%*s");	// atom_name
	fscanf(fp_geo,"%*s");	// c/s types
	fscanf(fp_geo,"%*lf");  // x coord
	fscanf(fp_geo,"%*lf");  // y coord
	fscanf(fp_geo,"%*lf");  // z coord


	fp_cur = ftell(fp_geo);
	fscanf(fp_geo,"%s",_if_fix);
	if( strcmp(_if_fix,"_fix_") == 0 )
		continue;
	else
		fseek(fp_geo,fp_cur,0);


	/*
	fscanf(fp_geo,"%*lf");  // charge
	fscanf(fp_geo,"%*lf");  // short A
	fscanf(fp_geo,"%*lf");  // short Rho
	fscanf(fp_geo,"%*lf");  // short C
	*/
        // read through classic ions
    }
//printf("sp read .. after read classic\n");
//printf("noa sp : %d\n",NumberOfSpIons);
    // READ sp ions
    for(int n=0;n<NumberOfSpIons;n++)
    {
//printf("sp alloc %d \n",n);
        pReturn[n].eigen_value = gsl_vector_calloc(4);
        //  Initialise eval with dummy lowest state by giving dummy val
        gsl_vector_set(pReturn[n].eigen_value,0,-100.);
        pReturn[n].eigen_vector = gsl_matrix_calloc(4,4);
        //   Initialise eigenvectors with pure s state
        for(int i=0;i<4;i++)    gsl_matrix_set(pReturn[n].eigen_vector,0,i,1.);
        pReturn[n].core_position = gsl_vector_calloc(3);
        pReturn[n].h_matrix = gsl_matrix_calloc(4,4);                                 // Zeroth Derivatives H
        pReturn[n].dh_matrix = (gsl_matrix**)malloc(3*sizeof(gsl_matrix*));           // First  Derivatives H
        for(int i=0;i<3;i++)    pReturn[n].dh_matrix[i] = gsl_matrix_calloc(4,4);

        pReturn[n].ddh_matrix = (gsl_matrix***)malloc(3*sizeof(gsl_matrix**));
        for(int i=0;i<3;i++)    
            pReturn[n].ddh_matrix[i] = (gsl_matrix**)malloc(3*sizeof(gsl_matrix*));
        for(int i=0;i<3;i++)
        {   for(int j=0;j<3;j++)
                pReturn[n].ddh_matrix[i][j] = gsl_matrix_calloc(4,4);     }           // Second Deriavtives H 

	/* READ NAME + CORE_COORD */
	// READ ATOM NAME
	fscanf(fp_geo,"%s",&pReturn[n].atom_name[0]);
	// READ COORD 
        fscanf(fp_geo,"%lf",&tmp);  gsl_vector_set(pReturn[n].core_position,0,tmp);
        fscanf(fp_geo,"%lf",&tmp);  gsl_vector_set(pReturn[n].core_position,1,tmp);
        fscanf(fp_geo,"%lf",&tmp);  gsl_vector_set(pReturn[n].core_position,2,tmp);
//printf("%s\t%12.6lf\t%12.6lf\t%12.6lf\n",pReturn[n].atom_name, gsl_vector_get(pReturn[n].core_position,0),gsl_vector_get(pReturn[n].core_position,1),gsl_vector_get(pReturn[n].core_position,2));

	fp_cur = ftell(fp_geo);			// RECORD CUR FILE PTR
	fscanf(fp_geo,"%s",_if_fix);
	if( strcmp(_if_fix,"_fix_") == 0 )	// CHECK IF FLAG EXISTS
		pReturn[n]._fix_flag = SP_SYSTEM_TRUE;
	else
	{	pReturn[n]._fix_flag = SP_SYSTEM_FALSE;
		fseek(fp_geo,fp_cur,0);
	}

	/* THIS SECTION, READING SP_PARAM + MM_POT WILL BE MOVED TO "sp_cluster_init()" */
	/*
        fscanf(fp_geo,"%lf",&pReturn[n].charge_core); 
        fscanf(fp_geo,"%lf",&pReturn[n].charge_shell);    // core and shell (sp-els) charges
        fscanf(fp_geo,"%lf",&pReturn[n].esp);
        // order of short-range parameter of sp-lone pair ... ## IMPORTANT CONVENTION AS -> AP ->  RS -> RP ORDER
        fscanf(fp_geo,"%lf",&pReturn[n].cent_short_range_a);
        fscanf(fp_geo,"%lf",&pReturn[n].cent_short_range_r);
	*/
    }
    fclose(fp_geo);
    return pReturn;
}


sp_cluster_type_classic_ion* sp_cluster_load_classic_parameter()
{
    sp_cluster_type_classic_ion* pReturn = NULL;
    int NumberOfClassicIons;    int NumberOfSpIons;
    char dummy[128];    double tmp;
    FILE* fp_geo = NULL;
	
	char _if_fix[8];	int fp_cur;

    char core_shell_checker[3];    

    fp_geo = fopen(SP_CLUSTER_SUPPORT_GEO_PATH,"r");
    fgets(dummy,sizeof(dummy),fp_geo);

    fscanf(fp_geo,"%d",&NumberOfClassicIons);
    fscanf(fp_geo,"%d",&NumberOfSpIons);

    pReturn = (sp_cluster_type_classic_ion*)malloc(NumberOfClassicIons*sizeof(sp_cluster_type_classic_ion));
	/// alocation done: Here "pReturn" is classic_ion type struct
	
	/// READ DATA
    for(int i=0;i<NumberOfClassicIons;i++)
    {   
	pReturn[i].core_position = gsl_vector_calloc(3);    // allocate pos var
	// load atom type
	fscanf(fp_geo,"%s",&pReturn[i].atom_name[0]);
	// check if it is core or shell
	memset(&core_shell_checker[0],0,3);
	fscanf(fp_geo,"%s",&core_shell_checker);	

	if( strcmp( core_shell_checker, "c" ) == 0 )			// check if its core
		pReturn[i].if_shell =  SP_SYSTEM_FALSE;			// set flag false
	else if( strcmp( core_shell_checker, "s" ) == 0 )		// check if its shell
		pReturn[i].if_shell =  SP_SYSTEM_TRUE;			// check flag true
	else
	{	printf("geo read shell/core checker error !\n");
		exit(1);
	}
        // load coordinates
        fscanf(fp_geo,"%lf",&tmp);  gsl_vector_set(pReturn[i].core_position,0,tmp);
        fscanf(fp_geo,"%lf",&tmp);  gsl_vector_set(pReturn[i].core_position,1,tmp);
        fscanf(fp_geo,"%lf",&tmp);  gsl_vector_set(pReturn[i].core_position,2,tmp);
	// READ:	ATOM_NAME	COORD_X		COORD_Y		COORD_Z


	fp_cur = ftell(fp_geo);
	fscanf(fp_geo,"%s",_if_fix);
	if( strcmp(_if_fix,"_fix_") == 0 )
		pReturn[i]._fix_flag = SP_SYSTEM_TRUE;
	else
	{	pReturn[i]._fix_flag = SP_SYSTEM_FALSE;
		fseek(fp_geo,fp_cur,0);
	}

	/* READ 	CHARGE		MM_SHORT_A	MM_SHORT_R	MM_SHORT_C	... BUCKINGHAM POT */
	// SECTION WILL BE TRANSFERED TO "sp_cluster_init()"
    }
    fclose(fp_geo);
    MPI_Barrier(MPI_COMM_WORLD);
    return pReturn;
}

///	///	///	///     ///     ///	///     ///     ///	///     ///     ///	///     ///     ///	///     ///     ///


sp_cluster_system* sp_cluster_system_init()
{   
    int number_of_sp_ion, number_of_classic_ion;
    sp_cluster_system* pReturn = NULL;
    FILE* fp_geo = NULL;    char dummy[128];
    fp_geo = fopen(SP_CLUSTER_SUPPORT_GEO_PATH,"r");
    fgets(dummy,sizeof(dummy),fp_geo);
    pReturn = (sp_cluster_system*)malloc(sizeof(sp_cluster_system));
    pReturn->sp_ion = sp_cluster_load_sp_parameter();               // Init sp-ions & parameteres
    pReturn->classic_ion = sp_cluster_load_classic_parameter();     // Init classic ions
    fscanf(fp_geo,"%d",&pReturn->number_of_classic_ion);
    fscanf(fp_geo,"%d",&pReturn->number_of_sp_ion);
    fclose(fp_geo);

    number_of_sp_ion = pReturn->number_of_sp_ion;   number_of_classic_ion = pReturn->number_of_classic_ion;

    // SCF INDEX
    pReturn->if_first_scf_trial = SP_SYSTEM_TRUE;

    // SCF WorkSpace
    pReturn->scf_h_matrix_vs_classic_ion     = (gsl_matrix**)malloc(pReturn->number_of_sp_ion*sizeof(gsl_matrix*));     // saving a sp_ion vs classical ions iteraction
    pReturn->scf_h_matrix_vs_sp_ion_monopole = (gsl_matrix**)malloc(pReturn->number_of_sp_ion*sizeof(gsl_matrix*));     // saving a sp_ion vs sp ions monopole iteraction
    pReturn->scf_h_matrix_vs_sp_ion_onsite   = (gsl_matrix**)malloc(pReturn->number_of_sp_ion*sizeof(gsl_matrix*));     // saving a sp_ion onsite term
    pReturn->scf_h_matrix_vs_sp_core         = (gsl_matrix**)malloc(pReturn->number_of_sp_ion*sizeof(gsl_matrix*));     // saving a sp_ion vs sp_ion_core Coulomb
    for(int i=0;i<pReturn->number_of_sp_ion;i++)
    {   pReturn->scf_h_matrix_vs_classic_ion[i]     = gsl_matrix_calloc(4,4);
        pReturn->scf_h_matrix_vs_sp_ion_monopole[i] = gsl_matrix_calloc(4,4);
        pReturn->scf_h_matrix_vs_sp_ion_onsite[i]   = gsl_matrix_calloc(4,4);        
        pReturn->scf_h_matrix_vs_sp_core[i]         = gsl_matrix_calloc(4,4);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    pReturn->scf_dh_x_matrix_vs_sp_ion_dipole = (gsl_matrix***)malloc(pReturn->number_of_sp_ion*sizeof(gsl_matrix**));
    pReturn->scf_dh_y_matrix_vs_sp_ion_dipole = (gsl_matrix***)malloc(pReturn->number_of_sp_ion*sizeof(gsl_matrix**));
    pReturn->scf_dh_z_matrix_vs_sp_ion_dipole = (gsl_matrix***)malloc(pReturn->number_of_sp_ion*sizeof(gsl_matrix**));
    for(int i=0;i<pReturn->number_of_sp_ion;i++)
    {   pReturn->scf_dh_x_matrix_vs_sp_ion_dipole[i] = (gsl_matrix**)malloc(pReturn->number_of_sp_ion*sizeof(gsl_matrix*));
        pReturn->scf_dh_y_matrix_vs_sp_ion_dipole[i] = (gsl_matrix**)malloc(pReturn->number_of_sp_ion*sizeof(gsl_matrix*));
        pReturn->scf_dh_z_matrix_vs_sp_ion_dipole[i] = (gsl_matrix**)malloc(pReturn->number_of_sp_ion*sizeof(gsl_matrix*));     }
    for(int i=0;i<pReturn->number_of_sp_ion;i++)
    {   for(int j=0;j<pReturn->number_of_sp_ion;j++)
        {   pReturn->scf_dh_x_matrix_vs_sp_ion_dipole[i][j] = gsl_matrix_calloc(4,4);
            pReturn->scf_dh_y_matrix_vs_sp_ion_dipole[i][j] = gsl_matrix_calloc(4,4);
            pReturn->scf_dh_z_matrix_vs_sp_ion_dipole[i][j] = gsl_matrix_calloc(4,4);       }}
    MPI_Barrier(MPI_COMM_WORLD);

    // FORCE WORKSPACE

    // vs cla - ion
    pReturn->cla_dh_matrix_workspace_x = (gsl_matrix***)malloc(number_of_sp_ion*sizeof(gsl_matrix**));    pReturn->cla_dh_matrix_workspace_y = (gsl_matrix***)malloc(number_of_sp_ion*sizeof(gsl_matrix**));
    pReturn->cla_dh_matrix_workspace_z = (gsl_matrix***)malloc(number_of_sp_ion*sizeof(gsl_matrix**));
    pReturn->cla_dh_matrix_x = (gsl_matrix***)malloc(number_of_sp_ion*sizeof(gsl_matrix**)); pReturn->cla_dh_matrix_y = (gsl_matrix***)malloc(number_of_sp_ion*sizeof(gsl_matrix**));
    pReturn->cla_dh_matrix_z = (gsl_matrix***)malloc(number_of_sp_ion*sizeof(gsl_matrix**));
    MPI_Barrier(MPI_COMM_WORLD);
    for(int n=0;n<number_of_sp_ion;n++)
    {   pReturn->cla_dh_matrix_workspace_x[n] = (gsl_matrix**)malloc(number_of_classic_ion*sizeof(gsl_matrix*));  pReturn->cla_dh_matrix_x[n] = (gsl_matrix**)malloc(number_of_classic_ion*sizeof(gsl_matrix*));
        pReturn->cla_dh_matrix_workspace_y[n] = (gsl_matrix**)malloc(number_of_classic_ion*sizeof(gsl_matrix*));  pReturn->cla_dh_matrix_y[n] = (gsl_matrix**)malloc(number_of_classic_ion*sizeof(gsl_matrix*));
        pReturn->cla_dh_matrix_workspace_z[n] = (gsl_matrix**)malloc(number_of_classic_ion*sizeof(gsl_matrix*));  pReturn->cla_dh_matrix_z[n] = (gsl_matrix**)malloc(number_of_classic_ion*sizeof(gsl_matrix*));    }
    MPI_Barrier(MPI_COMM_WORLD);
    for(int n=0;n<number_of_sp_ion;n++)
    {   for(int m=0;m<number_of_classic_ion;m++)
        {   pReturn->cla_dh_matrix_workspace_x[n][m] = gsl_matrix_calloc(4,4);   pReturn->cla_dh_matrix_workspace_y[n][m] = gsl_matrix_calloc(4,4);   pReturn->cla_dh_matrix_workspace_z[n][m] = gsl_matrix_calloc(4,4);
            pReturn->cla_dh_matrix_x[n][m] = gsl_matrix_calloc(4,4); pReturn->cla_dh_matrix_y[n][m] = gsl_matrix_calloc(4,4); pReturn->cla_dh_matrix_z[n][m] = gsl_matrix_calloc(4,4);                                         }}
    MPI_Barrier(MPI_COMM_WORLD);


    /// vs sp - ion

    // Memory allocation    ... sp-sp mono + core buffer
    pReturn->dh_matrix_workspace_x = (gsl_matrix***)malloc(number_of_sp_ion*sizeof(gsl_matrix**));    pReturn->dh_matrix_workspace_y = (gsl_matrix***)malloc(number_of_sp_ion*sizeof(gsl_matrix**));
    pReturn->dh_matrix_workspace_z = (gsl_matrix***)malloc(number_of_sp_ion*sizeof(gsl_matrix**));
    pReturn->dh_matrix_x = (gsl_matrix***)malloc(number_of_sp_ion*sizeof(gsl_matrix**)); pReturn->dh_matrix_y = (gsl_matrix***)malloc(number_of_sp_ion*sizeof(gsl_matrix**));
    pReturn->dh_matrix_z = (gsl_matrix***)malloc(number_of_sp_ion*sizeof(gsl_matrix**));
    MPI_Barrier(MPI_COMM_WORLD);
    for(int n=0;n<number_of_sp_ion;n++)
    {   pReturn->dh_matrix_workspace_x[n] = (gsl_matrix**)malloc(number_of_sp_ion*sizeof(gsl_matrix*));  pReturn->dh_matrix_x[n] = (gsl_matrix**)malloc(number_of_sp_ion*sizeof(gsl_matrix*));
        pReturn->dh_matrix_workspace_y[n] = (gsl_matrix**)malloc(number_of_sp_ion*sizeof(gsl_matrix*));  pReturn->dh_matrix_y[n] = (gsl_matrix**)malloc(number_of_sp_ion*sizeof(gsl_matrix*));
        pReturn->dh_matrix_workspace_z[n] = (gsl_matrix**)malloc(number_of_sp_ion*sizeof(gsl_matrix*));  pReturn->dh_matrix_z[n] = (gsl_matrix**)malloc(number_of_sp_ion*sizeof(gsl_matrix*));    }
    for(int n=0;n<number_of_sp_ion;n++)
    {   for(int m=0;m<number_of_sp_ion;m++)
        {   pReturn->dh_matrix_workspace_x[n][m] = gsl_matrix_calloc(4,4);   pReturn->dh_matrix_workspace_y[n][m] = gsl_matrix_calloc(4,4);   pReturn->dh_matrix_workspace_z[n][m] = gsl_matrix_calloc(4,4);
            pReturn->dh_matrix_x[n][m] = gsl_matrix_calloc(4,4); pReturn->dh_matrix_y[n][m] = gsl_matrix_calloc(4,4); pReturn->dh_matrix_z[n][m] = gsl_matrix_calloc(4,4);                                         }}
    MPI_Barrier(MPI_COMM_WORLD);

    // Memory allocation    ... sp-sp Dipolar Buffer
    pReturn->ddh_matrix_workspace_xx = (gsl_matrix***)malloc(number_of_sp_ion*sizeof(gsl_matrix**)); pReturn->ddh_matrix_workspace_xy = (gsl_matrix***)malloc(number_of_sp_ion*sizeof(gsl_matrix**));
    pReturn->ddh_matrix_workspace_xz = (gsl_matrix***)malloc(number_of_sp_ion*sizeof(gsl_matrix**)); pReturn->ddh_matrix_workspace_yx = (gsl_matrix***)malloc(number_of_sp_ion*sizeof(gsl_matrix**));
    pReturn->ddh_matrix_workspace_yy = (gsl_matrix***)malloc(number_of_sp_ion*sizeof(gsl_matrix**)); pReturn->ddh_matrix_workspace_yz = (gsl_matrix***)malloc(number_of_sp_ion*sizeof(gsl_matrix**));
    pReturn->ddh_matrix_workspace_zx = (gsl_matrix***)malloc(number_of_sp_ion*sizeof(gsl_matrix**)); pReturn->ddh_matrix_workspace_zy = (gsl_matrix***)malloc(number_of_sp_ion*sizeof(gsl_matrix**));
    pReturn->ddh_matrix_workspace_zz = (gsl_matrix***)malloc(number_of_sp_ion*sizeof(gsl_matrix**));

    pReturn->ddh_matrix_xx = (gsl_matrix***)malloc(number_of_sp_ion*sizeof(gsl_matrix**));   pReturn->ddh_matrix_xy = (gsl_matrix***)malloc(number_of_sp_ion*sizeof(gsl_matrix**));
    pReturn->ddh_matrix_xz = (gsl_matrix***)malloc(number_of_sp_ion*sizeof(gsl_matrix**));   pReturn->ddh_matrix_yx = (gsl_matrix***)malloc(number_of_sp_ion*sizeof(gsl_matrix**));
    pReturn->ddh_matrix_yy = (gsl_matrix***)malloc(number_of_sp_ion*sizeof(gsl_matrix**));   pReturn->ddh_matrix_yz = (gsl_matrix***)malloc(number_of_sp_ion*sizeof(gsl_matrix**));
    pReturn->ddh_matrix_zx = (gsl_matrix***)malloc(number_of_sp_ion*sizeof(gsl_matrix**));   pReturn->ddh_matrix_zy = (gsl_matrix***)malloc(number_of_sp_ion*sizeof(gsl_matrix**));
    pReturn->ddh_matrix_zz = (gsl_matrix***)malloc(number_of_sp_ion*sizeof(gsl_matrix**));
    MPI_Barrier(MPI_COMM_WORLD);

    for(int n=0;n<number_of_sp_ion;n++)
    {   
        pReturn->ddh_matrix_workspace_xx[n] = (gsl_matrix**)malloc(number_of_sp_ion*sizeof(gsl_matrix*));    pReturn->ddh_matrix_workspace_xy[n] = (gsl_matrix**)malloc(number_of_sp_ion*sizeof(gsl_matrix*));
        pReturn->ddh_matrix_workspace_xz[n] = (gsl_matrix**)malloc(number_of_sp_ion*sizeof(gsl_matrix*));    pReturn->ddh_matrix_workspace_yx[n] = (gsl_matrix**)malloc(number_of_sp_ion*sizeof(gsl_matrix*));
        pReturn->ddh_matrix_workspace_yy[n] = (gsl_matrix**)malloc(number_of_sp_ion*sizeof(gsl_matrix*));    pReturn->ddh_matrix_workspace_yz[n] = (gsl_matrix**)malloc(number_of_sp_ion*sizeof(gsl_matrix*));
        pReturn->ddh_matrix_workspace_zx[n] = (gsl_matrix**)malloc(number_of_sp_ion*sizeof(gsl_matrix*));    pReturn->ddh_matrix_workspace_zy[n] = (gsl_matrix**)malloc(number_of_sp_ion*sizeof(gsl_matrix*));
        pReturn->ddh_matrix_workspace_zz[n] = (gsl_matrix**)malloc(number_of_sp_ion*sizeof(gsl_matrix*));

        pReturn->ddh_matrix_xx[n] = (gsl_matrix**)malloc(number_of_sp_ion*sizeof(gsl_matrix*));    pReturn->ddh_matrix_xy[n] = (gsl_matrix**)malloc(number_of_sp_ion*sizeof(gsl_matrix*));
        pReturn->ddh_matrix_xz[n] = (gsl_matrix**)malloc(number_of_sp_ion*sizeof(gsl_matrix*));    pReturn->ddh_matrix_yx[n] = (gsl_matrix**)malloc(number_of_sp_ion*sizeof(gsl_matrix*));                        
        pReturn->ddh_matrix_yy[n] = (gsl_matrix**)malloc(number_of_sp_ion*sizeof(gsl_matrix*));    pReturn->ddh_matrix_yz[n] = (gsl_matrix**)malloc(number_of_sp_ion*sizeof(gsl_matrix*));                        
        pReturn->ddh_matrix_zx[n] = (gsl_matrix**)malloc(number_of_sp_ion*sizeof(gsl_matrix*));    pReturn->ddh_matrix_zy[n] = (gsl_matrix**)malloc(number_of_sp_ion*sizeof(gsl_matrix*));                        
        pReturn->ddh_matrix_zz[n] = (gsl_matrix**)malloc(number_of_sp_ion*sizeof(gsl_matrix*));                        
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for(int n=0;n<number_of_sp_ion;n++)
    {   for(int m=0;m<number_of_sp_ion;m++)
        {   
            pReturn->ddh_matrix_workspace_xx[n][m] = gsl_matrix_calloc(4,4);    pReturn->ddh_matrix_workspace_xy[n][m] = gsl_matrix_calloc(4,4);
            pReturn->ddh_matrix_workspace_xz[n][m] = gsl_matrix_calloc(4,4);    pReturn->ddh_matrix_workspace_yx[n][m] = gsl_matrix_calloc(4,4);
            pReturn->ddh_matrix_workspace_yy[n][m] = gsl_matrix_calloc(4,4);    pReturn->ddh_matrix_workspace_yz[n][m] = gsl_matrix_calloc(4,4);
            pReturn->ddh_matrix_workspace_zx[n][m] = gsl_matrix_calloc(4,4);    pReturn->ddh_matrix_workspace_zy[n][m] = gsl_matrix_calloc(4,4);
            pReturn->ddh_matrix_workspace_zz[n][m] = gsl_matrix_calloc(4,4);

            pReturn->ddh_matrix_xx[n][m] = gsl_matrix_calloc(4,4);    pReturn->ddh_matrix_xy[n][m] = gsl_matrix_calloc(4,4);
            pReturn->ddh_matrix_xz[n][m] = gsl_matrix_calloc(4,4);    pReturn->ddh_matrix_yx[n][m] = gsl_matrix_calloc(4,4);                   
            pReturn->ddh_matrix_yy[n][m] = gsl_matrix_calloc(4,4);    pReturn->ddh_matrix_yz[n][m] = gsl_matrix_calloc(4,4);                   
            pReturn->ddh_matrix_zx[n][m] = gsl_matrix_calloc(4,4);    pReturn->ddh_matrix_zy[n][m] = gsl_matrix_calloc(4,4);                   
            pReturn->ddh_matrix_zz[n][m] = gsl_matrix_calloc(4,4);                  
        }
    }

    // for evec deriv w.r.t. sp ion moves
    pReturn->deriv_evec_sp = (double****)malloc(number_of_sp_ion*sizeof(double***));
    for(int n=0;n<number_of_sp_ion;n++)
        pReturn->deriv_evec_sp[n] = (double***)malloc(number_of_sp_ion*sizeof(double**));
    for(int n=0;n<number_of_sp_ion;n++)
    {   for(int m=0;m<number_of_sp_ion;m++)
            pReturn->deriv_evec_sp[n][m] = (double**)malloc(4*sizeof(double*));
    }
    for(int n=0;n<number_of_sp_ion;n++)
    {   for(int m=0;m<number_of_sp_ion;m++)
        {   for(int o=0;o<4;o++)
                pReturn->deriv_evec_sp[n][m][o] = (double*)calloc(3,sizeof(double));
        }
    }
    // for evec deriv w.r.t. classic ion moves
    pReturn->deriv_evec_cla = (double****)malloc(number_of_classic_ion*sizeof(double***));
    for(int n=0;n<number_of_classic_ion;n++)
        pReturn->deriv_evec_cla[n] = (double***)malloc(number_of_sp_ion*sizeof(double**));
    for(int n=0;n<number_of_classic_ion;n++)
    {   for(int m=0;m<number_of_sp_ion;m++)
            pReturn->deriv_evec_cla[n][m] = (double**)malloc(4*sizeof(double*));
    }
    for(int n=0;n<number_of_classic_ion;n++)
    {   for(int m=0;m<number_of_sp_ion;m++)
        {   for(int o=0;o<4;o++)
                pReturn->deriv_evec_cla[n][m][o] = (double*)calloc(3,sizeof(double));
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    // Note that the force convention when 'n' sp core sees -> 'm' sp core

    return pReturn;
}


void sp_cluster_system_detach( sp_cluster_system* ptr )
{   
    const int number_of_sp_ion = ptr->number_of_sp_ion;
    const int number_of_classic_ion = ptr->number_of_classic_ion;
    
    // detach classic ion
    for(int n=0;n<ptr->number_of_classic_ion;n++)
        gsl_vector_free( ptr->classic_ion[n].core_position );
    // detach sp-lone pair ion
    // ptr->sp_ion
    for(int n=0;n<ptr->number_of_sp_ion;n++)
    {   gsl_vector_free(ptr->sp_ion[n].eigen_value);              // detach eval
        gsl_matrix_free(ptr->sp_ion[n].eigen_vector);             // detach evec
        gsl_matrix_free(ptr->sp_ion[n].h_matrix);                 // detach Hamiltonial matrix

        for(int i=0;i<3;i++)    
            gsl_matrix_free(ptr->sp_ion[n].dh_matrix[i]);
        free(ptr->sp_ion[n].dh_matrix);                           // detach first derivatives of Hamiltonian matrix

        for(int i=0;i<3;i++)
        {   for(int j=0;j<3;j++)    
                gsl_matrix_free(ptr->sp_ion[n].ddh_matrix[i][j]); 
        }
        for(int i=0;i<3;i++)        
            free(ptr->sp_ion[n].ddh_matrix[i]);
        free(ptr->sp_ion[n].ddh_matrix);                          // detach second dertivative of Hamiltonian matrix

        gsl_vector_free(ptr->sp_ion[n].core_position);            // detach core position

        free(ptr->sp_ion[n].knot);
        for(int i=0;i<ptr->sp_ion[n].number_of_knot-1;i++)
        {   free(ptr->sp_ion[n].radial_s_coefficient[i]);
            free(ptr->sp_ion[n].radial_p_coefficient[i]);
        }
        free( ptr->sp_ion[n].radial_s_coefficient );
        free( ptr->sp_ion[n].radial_p_coefficient );
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // detach scf workspace
    
    // detach vs_classic_ion/ vs_sp_ion_monopole/ vs_sp_ion_onsite
    for(int i=0;i<ptr->number_of_sp_ion;i++)
    {   gsl_matrix_free(ptr->scf_h_matrix_vs_classic_ion[i]);
        gsl_matrix_free(ptr->scf_h_matrix_vs_sp_ion_monopole[i]);
        gsl_matrix_free(ptr->scf_h_matrix_vs_sp_ion_onsite[i]);        
        gsl_matrix_free(ptr->scf_h_matrix_vs_sp_core[i]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    free(ptr->scf_h_matrix_vs_classic_ion);
    free(ptr->scf_h_matrix_vs_sp_ion_monopole);
    free(ptr->scf_h_matrix_vs_sp_ion_onsite);
    free(ptr->scf_h_matrix_vs_sp_core);
    MPI_Barrier(MPI_COMM_WORLD);
    // detach vs sp_dipole x,y and z terms
    for(int i=0;i<ptr->number_of_sp_ion;i++)
    {   for(int j=0;j<ptr->number_of_sp_ion;j++)
        {   gsl_matrix_free(ptr->scf_dh_x_matrix_vs_sp_ion_dipole[i][j]);
            gsl_matrix_free(ptr->scf_dh_y_matrix_vs_sp_ion_dipole[i][j]);
            gsl_matrix_free(ptr->scf_dh_z_matrix_vs_sp_ion_dipole[i][j]);       }}
    MPI_Barrier(MPI_COMM_WORLD);
    for(int i=0;i<ptr->number_of_sp_ion;i++)
    {   free(ptr->scf_dh_x_matrix_vs_sp_ion_dipole[i]);
        free(ptr->scf_dh_y_matrix_vs_sp_ion_dipole[i]);
        free(ptr->scf_dh_z_matrix_vs_sp_ion_dipole[i]);     }
    free(ptr->scf_dh_x_matrix_vs_sp_ion_dipole);
    free(ptr->scf_dh_y_matrix_vs_sp_ion_dipole);
    free(ptr->scf_dh_z_matrix_vs_sp_ion_dipole);
    MPI_Barrier(MPI_COMM_WORLD);

    //printf("Flag\n");
    // detach force workspace

    /// vs cla workspace
    for(int n=0;n<number_of_sp_ion;n++)
    {   for(int m=0;m<number_of_classic_ion;m++)
        {   gsl_matrix_free(ptr->cla_dh_matrix_workspace_x[n][m]);  gsl_matrix_free(ptr->cla_dh_matrix_x[n][m]);
            gsl_matrix_free(ptr->cla_dh_matrix_workspace_y[n][m]);  gsl_matrix_free(ptr->cla_dh_matrix_y[n][m]);
            gsl_matrix_free(ptr->cla_dh_matrix_workspace_z[n][m]);  gsl_matrix_free(ptr->cla_dh_matrix_z[n][m]);    }}
    MPI_Barrier(MPI_COMM_WORLD);
    for(int n=0;n<number_of_sp_ion;n++)
    {   free(ptr->cla_dh_matrix_workspace_x[n]); free(ptr->cla_dh_matrix_x[n]);
        free(ptr->cla_dh_matrix_workspace_y[n]); free(ptr->cla_dh_matrix_y[n]);
        free(ptr->cla_dh_matrix_workspace_z[n]); free(ptr->cla_dh_matrix_z[n]);   }
    free(ptr->cla_dh_matrix_workspace_x);    free(ptr->cla_dh_matrix_x);
    free(ptr->cla_dh_matrix_workspace_y);    free(ptr->cla_dh_matrix_y);
    free(ptr->cla_dh_matrix_workspace_z);    free(ptr->cla_dh_matrix_z);
    MPI_Barrier(MPI_COMM_WORLD);
    /// vs sp workspace

    // Force Memory Detach
    for(int n=0;n<number_of_sp_ion;n++)
    {   for(int m=0;m<number_of_sp_ion;m++)
        {   // monopole
            gsl_matrix_free(ptr->dh_matrix_workspace_x[n][m]);   gsl_matrix_free(ptr->dh_matrix_x[n][m]);
            gsl_matrix_free(ptr->dh_matrix_workspace_y[n][m]);   gsl_matrix_free(ptr->dh_matrix_y[n][m]);
            gsl_matrix_free(ptr->dh_matrix_workspace_z[n][m]);   gsl_matrix_free(ptr->dh_matrix_z[n][m]);
            // dipole
            gsl_matrix_free(ptr->ddh_matrix_workspace_xx[n][m]); gsl_matrix_free(ptr->ddh_matrix_workspace_xy[n][m]);
            gsl_matrix_free(ptr->ddh_matrix_workspace_xz[n][m]); gsl_matrix_free(ptr->ddh_matrix_workspace_yx[n][m]);
            gsl_matrix_free(ptr->ddh_matrix_workspace_yy[n][m]); gsl_matrix_free(ptr->ddh_matrix_workspace_yz[n][m]);
            gsl_matrix_free(ptr->ddh_matrix_workspace_zx[n][m]); gsl_matrix_free(ptr->ddh_matrix_workspace_zy[n][m]);
            gsl_matrix_free(ptr->ddh_matrix_workspace_zz[n][m]);

            gsl_matrix_free(ptr->ddh_matrix_xx[n][m]);   gsl_matrix_free(ptr->ddh_matrix_xy[n][m]);
            gsl_matrix_free(ptr->ddh_matrix_xz[n][m]);   gsl_matrix_free(ptr->ddh_matrix_yx[n][m]);                   
            gsl_matrix_free(ptr->ddh_matrix_yy[n][m]);   gsl_matrix_free(ptr->ddh_matrix_yz[n][m]);                   
            gsl_matrix_free(ptr->ddh_matrix_zx[n][m]);   gsl_matrix_free(ptr->ddh_matrix_zy[n][m]);                   
            gsl_matrix_free(ptr->ddh_matrix_zz[n][m]);                   
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for(int n=0;n<number_of_sp_ion;n++)
    {   // monopole
        free(ptr->dh_matrix_workspace_x[n]);   free(ptr->dh_matrix_workspace_y[n]);    free(ptr->dh_matrix_workspace_z[n]);
        free(ptr->dh_matrix_x[n]);   free(ptr->dh_matrix_y[n]);   free(ptr->dh_matrix_z[n]);
        // dipole
        free(ptr->ddh_matrix_workspace_xx[n]);   free(ptr->ddh_matrix_workspace_xy[n]);   free(ptr->ddh_matrix_workspace_xz[n]);   
        free(ptr->ddh_matrix_workspace_yx[n]);   free(ptr->ddh_matrix_workspace_yy[n]);   free(ptr->ddh_matrix_workspace_yz[n]);   
        free(ptr->ddh_matrix_workspace_zx[n]);   free(ptr->ddh_matrix_workspace_zy[n]);   free(ptr->ddh_matrix_workspace_zz[n]);   

        free(ptr->ddh_matrix_xx[n]); free(ptr->ddh_matrix_xy[n]); free(ptr->ddh_matrix_xz[n]);                                                 
        free(ptr->ddh_matrix_yx[n]); free(ptr->ddh_matrix_yy[n]); free(ptr->ddh_matrix_yz[n]);                                                 
        free(ptr->ddh_matrix_zx[n]); free(ptr->ddh_matrix_zy[n]); free(ptr->ddh_matrix_zz[n]);                                                 
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // monopole
    free(ptr->dh_matrix_workspace_x);    free(ptr->dh_matrix_workspace_y);    free(ptr->dh_matrix_workspace_z);
    free(ptr->dh_matrix_x);  free(ptr->dh_matrix_y);  free(ptr->dh_matrix_z);
    // dipole
    free(ptr->ddh_matrix_workspace_xx);  free(ptr->ddh_matrix_workspace_xy);  free(ptr->ddh_matrix_workspace_xz);
    free(ptr->ddh_matrix_workspace_yx);  free(ptr->ddh_matrix_workspace_yy);  free(ptr->ddh_matrix_workspace_yz);
    free(ptr->ddh_matrix_workspace_zx);  free(ptr->ddh_matrix_workspace_zy);  free(ptr->ddh_matrix_workspace_zz);

    free(ptr->ddh_matrix_xx);    free(ptr->ddh_matrix_xy);    free(ptr->ddh_matrix_xz);
    free(ptr->ddh_matrix_yx);    free(ptr->ddh_matrix_yy);    free(ptr->ddh_matrix_yz);
    free(ptr->ddh_matrix_zx);    free(ptr->ddh_matrix_zy);    free(ptr->ddh_matrix_zz);
    // End of Force Memory Detach
    

    // evec_deriv
    for(int n=0;n<number_of_sp_ion;n++)
    {   for(int m=0;m<number_of_sp_ion;m++)
        {   for(int o=0;o<4;o++)
            {   free(ptr->deriv_evec_sp[n][m][o]);      }}}
    for(int n=0;n<number_of_sp_ion;n++)
    {   for(int m=0;m<number_of_sp_ion;m++)
        {   free(ptr->deriv_evec_sp[n][m]);     }}
    for(int n=0;n<number_of_sp_ion;n++)
    {   free(ptr->deriv_evec_sp[n]);            }
    free(ptr->deriv_evec_sp);

    for(int n=0;n<number_of_classic_ion;n++)
    {   for(int m=0;m<number_of_sp_ion;m++)
        {   for(int o=0;o<4;o++)
            {   free(ptr->deriv_evec_cla[n][m][o]);      }}}
    for(int n=0;n<number_of_classic_ion;n++)
    {   for(int m=0;m<number_of_sp_ion;m++)
        {   free(ptr->deriv_evec_cla[n][m]);     }}
    for(int n=0;n<number_of_classic_ion;n++)
    {   free(ptr->deriv_evec_cla[n]);            }
    free(ptr->deriv_evec_cla);


    free(ptr);
    return;
}

void sp_cluster_system_get_spline_integral( sp_cluster_system* sp_sys )         // 0709 ... short-range integral only..
{   
	char atom_name_cmp[3];
	int sp_cnt = 0;	int classic_cnt = 0;
	char MM_type[10];	
	char atom_type[10];

    sp_sys->integral_lut = SP_SYSTEM_FALSE;

    FILE* fp = fopen(SP_CLUSTER_SUPPORT_SPECIES_PATH,"r");
    fscanf(fp,"%*s");	// read dummy single word comment 
    fscanf(fp,"%d",&sp_sys->number_of_species);
    sp_sys->number_of_qm_interaction_bm = sp_sys->number_of_species;

    const int number_of_classic_ion = sp_sys->number_of_species - 1;
    const int number_of_sp_ion = sp_sys->number_of_species - number_of_classic_ion + 1;
    double file_read_dummy;


    FILE* fp_knot = NULL; FILE* fp_radial_s = NULL; FILE* fp_radial_p = NULL;
    fp_knot = fopen(SP_CLUSTER_SUPPORT_KNOT_PATH,"r");
    fp_radial_s = fopen(SP_CLUSTER_SUPPORT_RADIAL_S_PATH,"r");
    fp_radial_p = fopen(SP_CLUSTER_SUPPORT_RADIAL_P_PATH,"r");


    sp_cluster_type_classic_ion* cla = (sp_cluster_type_classic_ion*)malloc(number_of_classic_ion*sizeof(sp_cluster_type_classic_ion));
    sp_cluster_type_sp_ion* sp       = (sp_cluster_type_sp_ion*)malloc(number_of_sp_ion*sizeof(sp_cluster_type_sp_ion));
    // allocate sp / cla pos vectors
    for(int i=0;i<number_of_classic_ion;i++)
        cla[i].core_position = gsl_vector_calloc(3);
    for(int i=0;i<number_of_sp_ion;i++)
        sp[i].core_position = gsl_vector_calloc(3);
    // at this point, dont need to consider the positions, they will be set just before the calculations
    
    
    // temporal saving space for potential parameters;
    double short_range_potential_a_vs_classic_ion[number_of_classic_ion];
    double short_range_potential_r_vs_classic_ion[number_of_classic_ion];
    double short_range_potential_c_vs_classic_ion[number_of_classic_ion];

    double short_range_potential_a_s_vs_sp_ion;
    double short_range_potential_r_s_vs_sp_ion;
    double short_range_potential_a_p_vs_sp_ion;
    double short_range_potential_r_p_vs_sp_ion;	// END OF QM_POT

	char if_shell[10];	memset(if_shell,0,10);	// if use shell model

	int number_of_MM_interaction;
	int number_of_atom_type;
	double MM_Q;	double MM_A;	double MM_R;	double MM_C;
	double MM_SP_QC;	double MM_SP_QS;	double MM_SP_ESP;
	double MM_SP_A;		double MM_SP_R;	// END OF MM_POT

	double MM_QS;	double MM_SPRING; /* k2 */	double MM_SPRING_K4;

    double MM_SP_VDW;
    memset(short_range_potential_a_vs_classic_ion,0.,number_of_classic_ion*sizeof(double));
    memset(short_range_potential_r_vs_classic_ion,0.,number_of_classic_ion*sizeof(double));

    // read short-range potential vs classic_ion
    for(int i=0;i<sp_sys->number_of_qm_interaction_bm-1;i++)
    {   
	//fscanf(fp,"%*s");	// read atom name
	fscanf(fp,"%s",sp_sys->interaction_qm_type_bm[i]);
        fscanf(fp,"%lf",&short_range_potential_a_vs_classic_ion[i]);
	sp_sys->interaction_qm_AR_bm[i][0] = short_range_potential_a_vs_classic_ion[i];
        fscanf(fp,"%lf",&short_range_potential_r_vs_classic_ion[i]);    
	sp_sys->interaction_qm_AR_bm[i][1] = short_range_potential_r_vs_classic_ion[i];
    }
    // read short-range potential vs sp_ion
    //fscanf(fp,"%*s");		// read atom name
    fscanf(fp,"%s",sp_sys->interaction_qm_type_bm[sp_sys->number_of_qm_interaction_bm-1]);
    fscanf(fp,"%lf",&short_range_potential_a_s_vs_sp_ion);
    sp_sys->interaction_qm_AR_bm[sp_sys->number_of_qm_interaction_bm-1][0] = short_range_potential_a_s_vs_sp_ion;
    fscanf(fp,"%lf",&short_range_potential_a_p_vs_sp_ion);
    sp_sys->interaction_qm_AR_bm[sp_sys->number_of_qm_interaction_bm-1][1] = short_range_potential_a_p_vs_sp_ion;
    fscanf(fp,"%lf",&short_range_potential_r_s_vs_sp_ion);
    sp_sys->interaction_qm_AR_bm[sp_sys->number_of_qm_interaction_bm-1][2] = short_range_potential_r_s_vs_sp_ion;
    fscanf(fp,"%lf",&short_range_potential_r_p_vs_sp_ion);
    sp_sys->interaction_qm_AR_bm[sp_sys->number_of_qm_interaction_bm-1][3] = short_range_potential_r_p_vs_sp_ion;
	// QM_SP_POT READ END

	///		///		///		///		///		///		///		///
	// HERE, NEED TO SET POT TABLE BY DEPENDING ON TYPES
	fscanf(fp,"%*s");			// read dummy single word comment 
   	fscanf(fp,"%d",&number_of_atom_type);	// read dummy number_of_species
	
	// READ ATOM TYPES ... 
	for(int i=0;i<number_of_atom_type;i++)
	{
		memset(&atom_type[0],0,10);	// refresh atom type
		fscanf(fp,"%s",&atom_type[0]);  // read atom type

		if ( strcmp(atom_type,"SHELL") == 0 )
		{
			memset(atom_name_cmp,0,sizeof(atom_name_cmp));
			fscanf(fp,"%s",&atom_name_cmp[0]);		// READ ATOM NAME
			fscanf(fp,"%lf",&MM_Q);				// READ CORE  CHARGE
			fscanf(fp,"%lf",&MM_QS);			// READ SHELL CHARGE
			fscanf(fp,"%lf",&MM_SPRING);			// READ SPRING CONSTANT
			fscanf(fp,"%lf",&MM_SPRING_K4);			// READ k4 spring const

			for(int j=0;j<sp_sys->number_of_classic_ion;j++)
			{
				if( ( strcmp( atom_name_cmp, sp_sys->classic_ion[j].atom_name ) == 0 ) /* if atom name is same */ && ( sp_sys->classic_ion[j].if_shell == SP_SYSTEM_TRUE ) /* if atom is shell */ )
				{
					sp_sys->classic_ion[j].charge_core = MM_QS;	// save shell charge
					sp_sys->classic_ion[j].k2_const    = MM_SPRING;	// save spring constant
					sp_sys->classic_ion[j].k4_const	   = MM_SPRING_K4;	// save anharmonic spring const

					classic_cnt++;
				}
				else if( strcmp( atom_name_cmp, sp_sys->classic_ion[j].atom_name ) == 0 /* if atom name is same */ && sp_sys->classic_ion[j].if_shell == SP_SYSTEM_FALSE /*  if atom is core (not shell) */ )
				{
					sp_sys->classic_ion[j].if_has_shell = SP_SYSTEM_TRUE;	// set shell flag ... this core has shell
					sp_sys->classic_ion[j].charge_core = MM_Q;	// save core charge
					sp_sys->classic_ion[j].k2_const    = MM_SPRING;	// save spring constant
					sp_sys->classic_ion[j].k4_const	   = MM_SPRING_K4;	// save anharmonic spring const

					classic_cnt++;
				}
				else if( strcmp( atom_name_cmp, sp_sys->classic_ion[j].atom_name ) != 0 )	// if atom name is not same, then continue (pass)
					continue;	
			}
		}
		else if ( strcmp(atom_type,"RIM") == 0 )		// Rigid Ion Model (RIM)
		{
			memset(atom_name_cmp,0,sizeof(atom_name_cmp));
			fscanf(fp,"%s",&atom_name_cmp[0]);		// READ ATOM NAME
			fscanf(fp,"%lf",&MM_Q);				// READ CORE  CHARGE

			for(int j=0;j<sp_sys->number_of_classic_ion;j++)
			{
				if( strcmp( atom_name_cmp, sp_sys->classic_ion[j].atom_name ) == 0 /* if atom name is same */ && sp_sys->classic_ion[j].if_shell == SP_SYSTEM_FALSE /* if atom is not shell (core) */ )
				{
					sp_sys->classic_ion[j].if_has_shell = SP_SYSTEM_FALSE;	// set shell flag false .. this is RIM
					sp_sys->classic_ion[j].charge_core = MM_Q;		// save core charge

					classic_cnt++;
				}
				else if( strcmp( atom_name_cmp, sp_sys->classic_ion[j].atom_name ) != 0 )	// if atom name is not in match, then continue (pass)
					continue;								
			}
		}
		else if ( strcmp(atom_type,"SP") == 0 )			// SP (SLAM)
		{
			memset(atom_name_cmp,0,sizeof(atom_name_cmp));
			fscanf(fp,"%s",&atom_name_cmp[0]);		// READ ATOM NAME
			fscanf(fp,"%lf",&MM_SP_QC);			// READ SP_CORE CHARGE
			fscanf(fp,"%lf",&MM_SP_QS);			// READ SP_LP	CHARGE
			fscanf(fp,"%lf",&MM_SP_ESP);			// READ SP_ESP	(ENERGY DIFFERENCE)

			for(int j=0;j<sp_sys->number_of_sp_ion;j++)
			{
				if( strcmp( atom_name_cmp, sp_sys->sp_ion[j].atom_name ) == 0 /* if atom name is same */ )
				{
					sp_sys->sp_ion[j].charge_core = MM_SP_QC;	// SAVE SP_CORE CHARGE
					sp_sys->sp_ion[j].charge_shell= MM_SP_QS;	// SAVE SP_LP	CHARGE
					sp_sys->sp_ion[j].esp	      = MM_SP_ESP;	// SAVE SP_ESP  (ENERGY_DIFFERENCE)

					sp_cnt++;
				}
				else
					continue;
			}
		}
		else
		{	// ERROR HANDELER
			printf("InputReadErr, atom type in 'sp_cluster_species.txt' has to be either one of SHELL, RIM and SP\n");
			exit(1);
		}
	}
	if( classic_cnt != sp_sys->number_of_classic_ion || sp_cnt != sp_sys->number_of_sp_ion )
	{	
		printf("InputReadErr, number of specified atoms in 'geo.txt' does not match with 'sp_cluster_species.txt'\n");
		exit(1);
	}


	/// READ POTENTIAL (INTERATOMIC MM)
	fscanf(fp,"%*s");				// read dummy single word comment 
   	fscanf(fp,"%d",&number_of_MM_interaction);	// read dummy number_of_MM_interaction

	sp_sys->number_of_mm_interaction_buck = number_of_MM_interaction;

	for(int i=0;i<number_of_MM_interaction;i++)
	{		
		for(int j=0;j<4;j++)
			memset(&sp_sys->interaction_type_buck[i][j][0],0,16*sizeof(char));
		memset(&sp_sys->interaction_ARC_buck[i][0],0.,3*sizeof(double));

		fscanf(fp,"%s",&sp_sys->interaction_type_buck[i][0][0]);	// READ INTERACTION BODY ... e.g., 'SHELL', 'RIM', 'SP' etc.	
		fscanf(fp,"%s",&sp_sys->interaction_type_buck[i][1][0]);	// READ BODY's TYPE	 ... e.g., 'O', 'Ba', 'Sn' etc.
		fscanf(fp,"%s",&sp_sys->interaction_type_buck[i][2][0]);	// READ INTERACTION BODY ... e.g., SAME WITH ABOVE
		fscanf(fp,"%s",&sp_sys->interaction_type_buck[i][3][0]);	// READ BODY's TYPE	 ... e.g., SAME WITH ABOVE

		fscanf(fp,"%lf",&sp_sys->interaction_ARC_buck[i][0]);		// A
		fscanf(fp,"%lf",&sp_sys->interaction_ARC_buck[i][1]);		// RHO
		fscanf(fp,"%lf",&sp_sys->interaction_ARC_buck[i][2]);		// C
	}

	//test
	//printf("%s\t%s\t%s\t%s\n",sp_sys->interaction_type_buck[0][0],sp_sys->interaction_type_buck[0][1],sp_sys->interaction_type_buck[0][2],sp_sys->interaction_type_buck[0][3]);
	//printf("%s\t%s\t%s\t%s\n",sp_sys->interaction_type_buck[1][0],sp_sys->interaction_type_buck[1][1],sp_sys->interaction_type_buck[1][2],sp_sys->interaction_type_buck[1][3]);
	//printf("test end \n");
	
	//exit(1);


    fclose(fp);
    ///     ///     ///     ///     ///     ///     /// END OF "sp_cluster_species.txt" READING


    ///     Get Radial Function
    int line_number;    double tmp_knot;    double tmp_s_coefficient[4];    double tmp_p_coefficient[4];

    fscanf(fp_knot,"%d",&line_number);      // # of knots
    for(int i=0;i<number_of_sp_ion;i++)
        sp[i].number_of_knot = line_number;    // read # of knots in ith sp - ion
    fscanf(fp_radial_s,"%*d");  fscanf(fp_radial_p,"%*d");
    for(int i=0;i<number_of_sp_ion;i++)
    {
        sp[i].knot = (double*)calloc(line_number,sizeof(double));        // knot type: (double*)
        sp[i].radial_s_coefficient = (double**)malloc((line_number-1)*sizeof(double));
        sp[i].radial_p_coefficient = (double**)malloc((line_number-1)*sizeof(double));   // type (double**)

        for(int j=0;j<line_number-1;j++)
        {   sp[i].radial_s_coefficient[j] = (double*)calloc(4,sizeof(double));
            sp[i].radial_p_coefficient[j] = (double*)calloc(4,sizeof(double));   }

    }// Memory allocation for sp_ion radial functions
    for(int i=0;i<line_number;i++)
    {   
        fscanf(fp_knot,"%lf", &tmp_knot);
        if( i < line_number-1 )
        {   fscanf(fp_radial_s,"%lf",&tmp_s_coefficient[0]);
            fscanf(fp_radial_s,"%lf",&tmp_s_coefficient[1]);
            fscanf(fp_radial_s,"%lf",&tmp_s_coefficient[2]);
            fscanf(fp_radial_s,"%lf",&tmp_s_coefficient[3]);

            fscanf(fp_radial_p,"%lf",&tmp_p_coefficient[0]);
            fscanf(fp_radial_p,"%lf",&tmp_p_coefficient[1]);
            fscanf(fp_radial_p,"%lf",&tmp_p_coefficient[2]);
            fscanf(fp_radial_p,"%lf",&tmp_p_coefficient[3]);
        }

        for(int j=0;j<number_of_sp_ion;j++)
        {
            sp[j].knot[i] = tmp_knot;      // save knot

            if( i < line_number-1 )
            {   sp[j].radial_s_coefficient[i][0] = tmp_s_coefficient[0];   sp[j].radial_s_coefficient[i][1] = tmp_s_coefficient[1];
                sp[j].radial_s_coefficient[i][2] = tmp_s_coefficient[2];   sp[j].radial_s_coefficient[i][3] = tmp_s_coefficient[3];   // save s function coefficients

                sp[j].radial_p_coefficient[i][0] = tmp_p_coefficient[0];   sp[j].radial_p_coefficient[i][1] = tmp_p_coefficient[1];
                sp[j].radial_p_coefficient[i][2] = tmp_p_coefficient[2];   sp[j].radial_p_coefficient[i][3] = tmp_p_coefficient[3];   // save p function coefficients
            }
        }
    }
    fclose(fp_knot);    fclose(fp_radial_s);    fclose(fp_radial_p);
    // radial input done

    ///


    double dist;
    double** data_ss, **data_sz, **data_xxyy, **data_zz;
    double** data_x_sx, **data_x_xz, **data_z_ss, **data_z_sz, **data_z_xxyy, **data_z_zz;
    double** data_xx_ss, **data_xx_sz, **data_xx_xx, **data_xx_yy, **data_xx_zz, **data_xy_xy, **data_xz_sx, **data_xz_xz;
    double** data_zz_ss, **data_zz_sz, **data_zz_xxyy, **data_zz_zz;
    // Data Samping buffers

    // Count Knot Stride and Create Knot
    dist = SP_INTEGRAL_SHORT_RANGE_START;
    sp_sys->knot_stride = 1;
    while( dist < SP_INTEGRAL_SHORT_RANGE_CUTOFF )
    {   dist = dist*SP_INTEGRAL_SHORT_RANGE_STEP;
        sp_sys->knot_stride++;                              }
    sp_sys->integral_knot = (double*)calloc(sp_sys->knot_stride,sizeof(double));

    sp_sys->integral_knot[0] = SP_INTEGRAL_SHORT_RANGE_START;
    dist = SP_INTEGRAL_SHORT_RANGE_START;
    for(int i=1;i<sp_sys->knot_stride;i++)
    {   dist = dist*SP_INTEGRAL_SHORT_RANGE_STEP;
        sp_sys->integral_knot[i] = dist;                    }
    //
    // sp_sys->knot_stride      ... saves knot length
    // sp_sys->integral_knot[i] ... saves knot points

    // Create tempolar data
    data_ss     = (double**)malloc(sp_sys->knot_stride*sizeof(double*));    data_sz    = (double**)malloc(sp_sys->knot_stride*sizeof(double*));
    data_xxyy   = (double**)malloc(sp_sys->knot_stride*sizeof(double*));    data_zz    = (double**)malloc(sp_sys->knot_stride*sizeof(double*));

    data_x_sx   = (double**)malloc(sp_sys->knot_stride*sizeof(double*));    data_x_xz  = (double**)malloc(sp_sys->knot_stride*sizeof(double*));

    data_z_ss   = (double**)malloc(sp_sys->knot_stride*sizeof(double*));    data_z_sz  = (double**)malloc(sp_sys->knot_stride*sizeof(double*));
    data_z_xxyy = (double**)malloc(sp_sys->knot_stride*sizeof(double*));    data_z_zz  = (double**)malloc(sp_sys->knot_stride*sizeof(double*));

    data_xx_ss  = (double**)malloc(sp_sys->knot_stride*sizeof(double*));    data_xx_sz = (double**)malloc(sp_sys->knot_stride*sizeof(double*));
    data_xx_xx  = (double**)malloc(sp_sys->knot_stride*sizeof(double*));    data_xx_yy = (double**)malloc(sp_sys->knot_stride*sizeof(double*));
    data_xx_zz  = (double**)malloc(sp_sys->knot_stride*sizeof(double*)); 
    
    data_xy_xy  = (double**)malloc(sp_sys->knot_stride*sizeof(double*));

    data_xz_sx  = (double**)malloc(sp_sys->knot_stride*sizeof(double*));    data_xz_xz = (double**)malloc(sp_sys->knot_stride*sizeof(double*));

    data_zz_ss  = (double**)malloc(sp_sys->knot_stride*sizeof(double*));    data_zz_sz = (double**)malloc(sp_sys->knot_stride*sizeof(double*));
    data_zz_xxyy= (double**)malloc(sp_sys->knot_stride*sizeof(double*));    data_zz_zz = (double**)malloc(sp_sys->knot_stride*sizeof(double*));


    for(int i=0;i<sp_sys->knot_stride;i++)
    {
        data_ss[i]      = (double*)calloc(2,sizeof(double));        data_sz[i]      =   (double*)calloc(2,sizeof(double));
        data_xxyy[i]    = (double*)calloc(2,sizeof(double));        data_zz[i]      =   (double*)calloc(2,sizeof(double));

        data_x_sx[i]    = (double*)calloc(2,sizeof(double));        data_x_xz[i]    =   (double*)calloc(2,sizeof(double));
        
        data_z_ss[i]    = (double*)calloc(2,sizeof(double));        data_z_sz[i]    =   (double*)calloc(2,sizeof(double));
        data_z_xxyy[i]  = (double*)calloc(2,sizeof(double));        data_z_zz[i]    =   (double*)calloc(2,sizeof(double));

        data_xx_ss[i]   = (double*)calloc(2,sizeof(double));        data_xx_sz[i]   =   (double*)calloc(2,sizeof(double));
        data_xx_xx[i]   = (double*)calloc(2,sizeof(double));        data_xx_yy[i]   =   (double*)calloc(2,sizeof(double));
        data_xx_zz[i]   = (double*)calloc(2,sizeof(double));

        data_xy_xy[i]   = (double*)calloc(2,sizeof(double));

        data_xz_sx[i]   = (double*)calloc(2,sizeof(double));        data_xz_xz[i]   =   (double*)calloc(2,sizeof(double));

        data_zz_ss[i]   = (double*)calloc(2,sizeof(double));        data_zz_sz[i]   =   (double*)calloc(2,sizeof(double));
        data_zz_xxyy[i] = (double*)calloc(2,sizeof(double));        data_zz_zz[i]   =   (double*)calloc(2,sizeof(double));
    }

    
    // allocation pre-calculated short-range integrals

    // sp vs classic ions
    //
    // Convention ... [0] -> anion // [1] -> addtional cation
    
    sp_sys->integral_vs_cla_s_ss   = (double***)malloc((number_of_classic_ion)*sizeof(double**));
    sp_sys->integral_vs_cla_s_sz   = (double***)malloc((number_of_classic_ion)*sizeof(double**));
    sp_sys->integral_vs_cla_s_xxyy = (double***)malloc((number_of_classic_ion)*sizeof(double**));
    sp_sys->integral_vs_cla_s_zz   = (double***)malloc((number_of_classic_ion)*sizeof(double**));

    sp_sys->integral_vs_cla_s_x_sx = (double***)malloc((number_of_classic_ion)*sizeof(double**));
    sp_sys->integral_vs_cla_s_x_xz = (double***)malloc((number_of_classic_ion)*sizeof(double**));

    sp_sys->integral_vs_cla_s_z_ss   = (double***)malloc((number_of_classic_ion)*sizeof(double**));
    sp_sys->integral_vs_cla_s_z_sz   = (double***)malloc((number_of_classic_ion)*sizeof(double**));
    sp_sys->integral_vs_cla_s_z_xxyy = (double***)malloc((number_of_classic_ion)*sizeof(double**));
    sp_sys->integral_vs_cla_s_z_zz   = (double***)malloc((number_of_classic_ion)*sizeof(double**));
    
    MPI_Barrier(MPI_COMM_WORLD);

    //sp_cluster_integrator_get_sh_11_element( &sp_sys->sp_ion[idx_sp1], &sp_sys->classic_ion[idx_cla] );
    //sp_cluster_integrator_get_sh_x_12_element(&sp_sys->sp_ion[idxn],&sp_sys->classic_ion[idxm])

    /* 02/11/2019 
     *
     * prbly in this section, the short-ranage parameters for " sp-cation vs classic-ion " can be set
     *
     *
     * e.g.,
     *
     * #    set short-range parameter
     *
     * FOR i < NUMBER_OF_CLASSIC_IONS, i++
     *
     *  DO  
     *      SHORT_RANGE POTENTIAL PARAMETERS
     *          CLA[i] vs SP-CATION
     *  DONE
     *
     * ENDFOR
     *
     */

    // Phase 1 : classic ion vs sp ion LUT calc
    //
    // sp ion at origin vs classic ion on z-axis
    gsl_vector_set(sp[0].core_position,0,0.);   gsl_vector_set(sp[0].core_position,1,0.);   gsl_vector_set(sp[0].core_position,2,0.);

    for(int i=0;i<number_of_classic_ion;i++)
    {
        gsl_vector_set(cla[i].core_position,0,0.);  gsl_vector_set(cla[i].core_position,1,0.);  // set 'i'th cla ion x=y=0

        // new feature 02/11/2019 ... set potentail

        cla[i].short_range_a = short_range_potential_a_vs_classic_ion[i];   cla[i].short_range_r = short_range_potential_r_vs_classic_ion[i];
        sp[0].short_range_a_s = short_range_potential_a_vs_classic_ion[i];  sp[0].short_range_a_p = short_range_potential_a_vs_classic_ion[i];
        sp[0].short_range_r_s = short_range_potential_r_vs_classic_ion[i];  sp[0].short_range_r_p = short_range_potential_r_vs_classic_ion[i];
        // sp-e SR pot vs ith classic ion;

        // potential set done

        for(int j=0;j<sp_sys->knot_stride;j++)
        {
            gsl_vector_set(cla[i].core_position,2,sp_sys->integral_knot[j]);        // set 'i'th cal ion z = knot[j]
            
            data_ss[j][0]   = sp_sys->integral_knot[j];         data_ss[j][1]   = sp_cluster_integrator_get_sh_11_element( sp_sys, &sp[0], &cla[i] ); 
            data_sz[j][0]   = sp_sys->integral_knot[j];         data_sz[j][1]   = sp_cluster_integrator_get_sh_14_element( sp_sys, &sp[0], &cla[i] );
            data_xxyy[j][0] = sp_sys->integral_knot[j];         data_xxyy[j][1] = sp_cluster_integrator_get_sh_2233_element( sp_sys, &sp[0], &cla[i] );
            data_zz[j][0]   = sp_sys->integral_knot[j];         data_zz[j][1]   = sp_cluster_integrator_get_sh_44_element( sp_sys, &sp[0], &cla[i] );
            
            data_x_sx[j][0] = sp_sys->integral_knot[j];         data_x_sx[j][1] = sp_cluster_integrator_get_sh_x_12_element( sp_sys, &sp[0], &cla[i] );
            data_x_xz[j][0] = sp_sys->integral_knot[j];         data_x_xz[j][1] = sp_cluster_integrator_get_sh_x_24_element( sp_sys, &sp[0], &cla[i] );

            data_z_ss[j][0]   = sp_sys->integral_knot[j];         data_z_ss[j][1]   = sp_cluster_integrator_get_sh_z_11_element( sp_sys, &sp[0], &cla[i] );
            data_z_sz[j][0]   = sp_sys->integral_knot[j];         data_z_sz[j][1]   = sp_cluster_integrator_get_sh_z_14_element( sp_sys, &sp[0], &cla[i] );
            data_z_xxyy[j][0] = sp_sys->integral_knot[j];         data_z_xxyy[j][1] = sp_cluster_integrator_get_sh_z_2233_element( sp_sys, &sp[0], &cla[i] );
            data_z_zz[j][0]   = sp_sys->integral_knot[j];         data_z_zz[j][1]   = sp_cluster_integrator_get_sh_z_44_element( sp_sys, &sp[0], &cla[i] );
    
            #ifdef DEBUG
            fprintf(stdout,"%lf\t\t%e\t%e\t%e\t%e\t\t%e\t%e\t%e\t%e\n",sp_sys->integral_knot[j],
                    data_ss[j][1],data_sz[j][1],data_xxyy[j][1],data_zz[j][1],
                    data_z_ss[j][1],data_z_sz[j][1],data_z_xxyy[j][1],data_z_zz[j][1]);
            #endif
        }
        // CALL C-SPLINE FUNCTION
        //
        //
        //
        //  input ... data_ss[j][..]  ... sp_sys->knot_stride ... 
        //  output ... sp_sys->integral_vs_cla_s_z_zz[i][j]  
        //  or Return;
        sp_sys->integral_vs_cla_s_ss[i]   = sp_cluster_support_get_spline( (const double**)data_ss, sp_sys->knot_stride );
        sp_sys->integral_vs_cla_s_sz[i]   = sp_cluster_support_get_spline( (const double**)data_sz, sp_sys->knot_stride );
        sp_sys->integral_vs_cla_s_xxyy[i] = sp_cluster_support_get_spline( (const double**)data_xxyy, sp_sys->knot_stride );
        sp_sys->integral_vs_cla_s_zz[i]   = sp_cluster_support_get_spline( (const double**)data_zz, sp_sys->knot_stride );

        sp_sys->integral_vs_cla_s_x_sx[i] = sp_cluster_support_get_spline( (const double**)data_x_sx, sp_sys->knot_stride );
        sp_sys->integral_vs_cla_s_x_xz[i] = sp_cluster_support_get_spline( (const double**)data_x_xz, sp_sys->knot_stride );

        sp_sys->integral_vs_cla_s_z_ss[i]   = sp_cluster_support_get_spline( (const double**)data_z_ss, sp_sys->knot_stride );
        sp_sys->integral_vs_cla_s_z_sz[i]   = sp_cluster_support_get_spline( (const double**)data_z_sz, sp_sys->knot_stride );
        sp_sys->integral_vs_cla_s_z_xxyy[i] = sp_cluster_support_get_spline( (const double**)data_z_xxyy, sp_sys->knot_stride );
        sp_sys->integral_vs_cla_s_z_zz[i]   = sp_cluster_support_get_spline( (const double**)data_z_zz, sp_sys->knot_stride );

    }

    #ifdef DEBUG
    fp = fopen("dummy1.txt","w");
    for(int i=0;i<sp_sys->knot_stride;i++)
    {
        fprintf(fp,"%lf\t\t%e\t%e\t%e\t%e\t\t%e\t%e\t%e\t%e\n",sp_sys->integral_knot[i],
                data_ss[i][1],data_sz[i][1],data_xxyy[i][1],data_zz[i][1],
                data_z_ss[i][1],data_z_sz[i][1],data_z_xxyy[i][1],data_z_zz[i][1]);
    }
    fclose(fp);
    #endif

    /*  02/11/2019
     *
     *  SET " SP-ELEC VS SP-ELEC " SHORT_RANGE POTENTIAL AGAIN BEFORE MAKING SPLINE FUNCTION
     *
     */
    for(int i=0;i<2;i++)
    {   sp[i].short_range_a_s = short_range_potential_a_s_vs_sp_ion;    sp[i].short_range_a_p = short_range_potential_a_p_vs_sp_ion;
        sp[i].short_range_r_s = short_range_potential_r_s_vs_sp_ion;    sp[i].short_range_r_p = short_range_potential_r_p_vs_sp_ion;    }
    // resetting sp-e vs sp-e shortrange potential

    // Phase 2 : sp ion vs sp ion LUT calc
    //
    // 1 sp ion at origin vs 2 sp ion on z-axis
    gsl_vector_set(sp[0].core_position,0,0.);   gsl_vector_set(sp[0].core_position,1,0.);   gsl_vector_set(sp[0].core_position,2,0.);
    gsl_vector_set(sp[1].core_position,0,0.);   gsl_vector_set(sp[1].core_position,1,0.);   gsl_vector_set(sp[1].core_position,2,0.);

    for(int i=0;i<sp_sys->knot_stride;i++)
    {   
        gsl_vector_set(sp[1].core_position,2,sp_sys->integral_knot[i]);        // set '2'th sp ion z = knot[j], i.e., moving sp[1] on z-axis
        
        // monopole
        data_ss[i][0]   = sp_sys->integral_knot[i];         data_ss[i][1]   = sp_cluster_spsp_mono_integrator_get_sh_11_element( sp_sys, &sp[1], &sp[0] );                  
        data_sz[i][0]   = sp_sys->integral_knot[i];         data_sz[i][1]   = sp_cluster_spsp_mono_integrator_get_sh_14_element( sp_sys, &sp[1], &sp[0] );
        data_xxyy[i][0] = sp_sys->integral_knot[i];         data_xxyy[i][1] = sp_cluster_spsp_mono_integrator_get_sh_2233_element( sp_sys, &sp[1], &sp[0] );
        data_zz[i][0]   = sp_sys->integral_knot[i];         data_zz[i][1]   = sp_cluster_spsp_mono_integrator_get_sh_44_element( sp_sys, &sp[1], &sp[0] );

        // dipole
        data_x_sx[i][0] = sp_sys->integral_knot[i];         data_x_sx[i][1] = sp_cluster_spsp_di_integrator_get_sh_x_12_element( sp_sys, &sp[1], &sp[0] );
        data_x_xz[i][0] = sp_sys->integral_knot[i];         data_x_xz[i][1] = sp_cluster_spsp_di_integrator_get_sh_x_24_element( sp_sys, &sp[1], &sp[0] );

        data_z_ss[i][0]   = sp_sys->integral_knot[i];       data_z_ss[i][1]   = sp_cluster_spsp_di_integrator_get_sh_z_11_element( sp_sys, &sp[1], &sp[0] );
        data_z_sz[i][0]   = sp_sys->integral_knot[i];       data_z_sz[i][1]   = sp_cluster_spsp_di_integrator_get_sh_z_14_element( sp_sys, &sp[1], &sp[0] );
        data_z_xxyy[i][0] = sp_sys->integral_knot[i];       data_z_xxyy[i][1] = sp_cluster_spsp_di_integrator_get_sh_z_2233_element( sp_sys, &sp[1], &sp[0] );
        data_z_zz[i][0]   = sp_sys->integral_knot[i];       data_z_zz[i][1]   = sp_cluster_spsp_di_integrator_get_sh_z_44_element( sp_sys, &sp[1], &sp[0] );

        // quadrupole
        data_xx_ss[i][0]   = sp_sys->integral_knot[i];      data_xx_ss[i][1]   = sp_cluster_spsp_quadru_integrator_get_sh_xx_11_element( sp_sys, &sp[1], &sp[0] );
        data_xx_sz[i][0]   = sp_sys->integral_knot[i];      data_xx_sz[i][1]   = sp_cluster_spsp_quadru_integrator_get_sh_xx_14_element( sp_sys, &sp[1], &sp[0] );
        data_xx_xx[i][0]   = sp_sys->integral_knot[i];      data_xx_xx[i][1]   = sp_cluster_spsp_quadru_integrator_get_sh_xx_22_element( sp_sys, &sp[1], &sp[0] );
        data_xx_yy[i][0]   = sp_sys->integral_knot[i];      data_xx_yy[i][1]   = sp_cluster_spsp_quadru_integrator_get_sh_xx_33_element( sp_sys, &sp[1], &sp[0] );
        data_xx_zz[i][0]   = sp_sys->integral_knot[i];      data_xx_zz[i][1]   = sp_cluster_spsp_quadru_integrator_get_sh_xx_44_element( sp_sys, &sp[1], &sp[0] );

        data_xy_xy[i][0]   = sp_sys->integral_knot[i];      data_xy_xy[i][1]   = sp_cluster_spsp_quadru_integrator_get_sh_xy_23_element( sp_sys, &sp[1], &sp[0] );

        data_xz_sx[i][0]   = sp_sys->integral_knot[i];      data_xz_sx[i][1]   = sp_cluster_spsp_quadru_integrator_get_sh_xz_12_element( sp_sys, &sp[1], &sp[0] );
        data_xz_xz[i][0]   = sp_sys->integral_knot[i];      data_xz_xz[i][1]   = sp_cluster_spsp_quadru_integrator_get_sh_xz_24_element( sp_sys, &sp[1], &sp[0] );
    
        data_zz_ss[i][0]   = sp_sys->integral_knot[i];      data_zz_ss[i][1]   = sp_cluster_spsp_quadru_integrator_get_sh_zz_11_element( sp_sys, &sp[1], &sp[0] );
        data_zz_sz[i][0]   = sp_sys->integral_knot[i];      data_zz_sz[i][1]   = sp_cluster_spsp_quadru_integrator_get_sh_zz_14_element( sp_sys, &sp[1], &sp[0] );
        data_zz_xxyy[i][0] = sp_sys->integral_knot[i];      data_zz_xxyy[i][1] = sp_cluster_spsp_quadru_integrator_get_sh_zz_2233_element( sp_sys, &sp[1], &sp[0] );
        data_zz_zz[i][0]   = sp_sys->integral_knot[i];      data_zz_zz[i][1]   = sp_cluster_spsp_quadru_integrator_get_sh_zz_44_element( sp_sys, &sp[1], &sp[0] );

        #ifdef DEBUG
        printf("spsp # .. %.3d\n",i+1);
        #endif

    }
    // CALL C-SPLINE FUNCTION
    sp_sys->integral_vs_sp_s_ss   = sp_cluster_support_get_spline( (const double**)data_ss, sp_sys->knot_stride );
    sp_sys->integral_vs_sp_s_sz   = sp_cluster_support_get_spline( (const double**)data_sz, sp_sys->knot_stride );
    sp_sys->integral_vs_sp_s_xxyy = sp_cluster_support_get_spline( (const double**)data_xxyy, sp_sys->knot_stride );
    sp_sys->integral_vs_sp_s_zz   = sp_cluster_support_get_spline( (const double**)data_zz, sp_sys->knot_stride );

    sp_sys->integral_vs_sp_s_x_sx = sp_cluster_support_get_spline( (const double**)data_x_sx, sp_sys->knot_stride );
    sp_sys->integral_vs_sp_s_x_xz = sp_cluster_support_get_spline( (const double**)data_x_xz, sp_sys->knot_stride );

    sp_sys->integral_vs_sp_s_z_ss   = sp_cluster_support_get_spline( (const double**)data_z_ss, sp_sys->knot_stride );
    sp_sys->integral_vs_sp_s_z_sz   = sp_cluster_support_get_spline( (const double**)data_z_sz, sp_sys->knot_stride );
    sp_sys->integral_vs_sp_s_z_xxyy = sp_cluster_support_get_spline( (const double**)data_z_xxyy, sp_sys->knot_stride );
    sp_sys->integral_vs_sp_s_z_zz   = sp_cluster_support_get_spline( (const double**)data_z_zz, sp_sys->knot_stride );

    sp_sys->integral_vs_sp_s_xx_ss    = sp_cluster_support_get_spline( (const double**)data_xx_ss, sp_sys->knot_stride );
    sp_sys->integral_vs_sp_s_xx_sz    = sp_cluster_support_get_spline( (const double**)data_xx_sz, sp_sys->knot_stride );
    sp_sys->integral_vs_sp_s_xx_xx    = sp_cluster_support_get_spline( (const double**)data_xx_xx, sp_sys->knot_stride );
    sp_sys->integral_vs_sp_s_xx_yy    = sp_cluster_support_get_spline( (const double**)data_xx_yy, sp_sys->knot_stride );
    sp_sys->integral_vs_sp_s_xx_zz    = sp_cluster_support_get_spline( (const double**)data_xx_zz, sp_sys->knot_stride );

    sp_sys->integral_vs_sp_s_xy_xy    = sp_cluster_support_get_spline( (const double**)data_xy_xy, sp_sys->knot_stride );

    sp_sys->integral_vs_sp_s_xz_sx    = sp_cluster_support_get_spline( (const double**)data_xz_sx, sp_sys->knot_stride );
    sp_sys->integral_vs_sp_s_xz_xz    = sp_cluster_support_get_spline( (const double**)data_xz_xz, sp_sys->knot_stride );

    sp_sys->integral_vs_sp_s_zz_ss   = sp_cluster_support_get_spline( (const double**)data_zz_ss, sp_sys->knot_stride );
    sp_sys->integral_vs_sp_s_zz_sz   = sp_cluster_support_get_spline( (const double**)data_zz_sz, sp_sys->knot_stride );
    sp_sys->integral_vs_sp_s_zz_xxyy = sp_cluster_support_get_spline( (const double**)data_zz_xxyy, sp_sys->knot_stride );
    sp_sys->integral_vs_sp_s_zz_zz   = sp_cluster_support_get_spline( (const double**)data_zz_zz, sp_sys->knot_stride );
    // END OF SPLINE

    #ifdef DEBUG
    fp = fopen("dummy2.txt","w");
    for(int i=0;i<sp_sys->knot_stride;i++)
    {
        fprintf(fp,"%lf\t\t%e\t%e\t%e\t%e\t\t%e\t%e\t%e\t%e\t\t%e\t%e\t%e\t%e\n",sp_sys->integral_knot[i],
                data_ss[i][1],data_sz[i][1],data_xxyy[i][1],data_zz[i][1],
                data_z_ss[i][1],data_z_sz[i][1],data_z_xxyy[i][1],data_z_zz[i][1],
                data_zz_ss[i][1],data_zz_sz[i][1],data_zz_xxyy[i][1],data_zz_zz[i][1]);
    }
    fclose(fp);

    //sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot );

    dist = 5.46;
    printf("!!!!!!dist : %lf // knot: %d\n", dist, sp_cluster_integrator_lut_b_search( dist, sp_sys->knot_stride, (const double*)sp_sys->integral_knot ));

    #endif


    // SET LUT DONE
    sp_sys->integral_lut = SP_SYSTEM_TRUE;

    /* ******************** FINALISE WORKSPACE DEALLOCATION ******************** */

    // dealloc
    for(int i=0;i<sp_sys->knot_stride;i++)
    {   free(data_ss[i]);       free(data_sz[i]);       free(data_xxyy[i]);     free(data_zz[i]);
        free(data_x_sx[i]);     free(data_x_xz[i]);
        free(data_z_ss[i]);     free(data_z_sz[i]);     free(data_z_xxyy[i]);   free(data_z_zz[i]);
        free(data_xx_ss[i]);    free(data_xx_sz[i]);    free(data_xx_xx[i]);    free(data_xx_yy[i]);    free(data_xx_zz[i]);
        free(data_xy_xy[i]);
        free(data_xz_sx[i]);    free(data_xz_xz[i]);
        free(data_zz_ss[i]);    free(data_zz_sz[i]);    free(data_zz_xxyy[i]);  free(data_zz_zz[i]);
    }
    free(data_ss);      free(data_sz);      free(data_xxyy);        free(data_zz);
    free(data_x_sx);    free(data_x_xz);
    free(data_z_ss);    free(data_z_sz);    free(data_z_xxyy);      free(data_z_zz);
    free(data_xx_ss);   free(data_xx_sz);   free(data_xx_xx);       free(data_xx_yy);   free(data_xx_zz);
    free(data_xy_xy);
    free(data_xz_sx);   free(data_xz_xz);
    free(data_zz_ss);   free(data_zz_sz);   free(data_zz_xxyy);     free(data_zz_zz);


    for(int i=0;i<number_of_classic_ion;i++)
        gsl_vector_free(cla[i].core_position);
    free(cla);

    for(int n=0;n<number_of_sp_ion;n++)
    {   
        gsl_vector_free(sp[n].core_position);            // detach core position

        free(sp[n].knot);
        for(int i=0;i<sp[n].number_of_knot-1;i++)
        {   free(sp[n].radial_s_coefficient[i]);
            free(sp[n].radial_p_coefficient[i]);
        }
        free( sp[n].radial_s_coefficient );
        free( sp[n].radial_p_coefficient );
    }
    free(sp);

    MPI_Barrier(MPI_COMM_WORLD);
        
    return;
}

void sp_cluster_system_free_spline_integral( sp_cluster_system* sp_sys )
{
    // free cla
    for(int i=0;i<sp_sys->number_of_species-1;i++)
    {
        for(int j=0;j<sp_sys->knot_stride-1;j++)
        {   free( sp_sys->integral_vs_cla_s_ss[i][j] );       free( sp_sys->integral_vs_cla_s_sz[i][j] );       free( sp_sys->integral_vs_cla_s_xxyy[i][j] );       free( sp_sys->integral_vs_cla_s_zz[i][j] );
            free( sp_sys->integral_vs_cla_s_x_sx[i][j] );     free( sp_sys->integral_vs_cla_s_x_xz[i][j] );
            free( sp_sys->integral_vs_cla_s_z_ss[i][j] );     free( sp_sys->integral_vs_cla_s_z_sz[i][j] );     free( sp_sys->integral_vs_cla_s_z_xxyy[i][j] );     free( sp_sys->integral_vs_cla_s_z_zz[i][j] );
        }
        free( sp_sys->integral_vs_cla_s_ss[i] );       free( sp_sys->integral_vs_cla_s_sz[i] );       free( sp_sys->integral_vs_cla_s_xxyy[i] );       free( sp_sys->integral_vs_cla_s_zz[i] );
        free( sp_sys->integral_vs_cla_s_x_sx[i] );     free( sp_sys->integral_vs_cla_s_x_xz[i] );
        free( sp_sys->integral_vs_cla_s_z_ss[i] );     free( sp_sys->integral_vs_cla_s_z_sz[i] );     free( sp_sys->integral_vs_cla_s_z_xxyy[i] );     free( sp_sys->integral_vs_cla_s_z_zz[i] );
    }
    free( sp_sys->integral_vs_cla_s_ss );       free( sp_sys->integral_vs_cla_s_sz );       free( sp_sys->integral_vs_cla_s_xxyy );       free( sp_sys->integral_vs_cla_s_zz );
    free( sp_sys->integral_vs_cla_s_x_sx );     free( sp_sys->integral_vs_cla_s_x_xz );
    free( sp_sys->integral_vs_cla_s_z_ss );     free( sp_sys->integral_vs_cla_s_z_sz );     free( sp_sys->integral_vs_cla_s_z_xxyy );     free( sp_sys->integral_vs_cla_s_z_zz );

    // free sp
    #ifdef DEBUG
    printf(" CHK POINT \n");
    #endif


    for(int i=0;i<sp_sys->knot_stride-1;i++)
    {
        free( sp_sys->integral_vs_sp_s_ss[i] );       free( sp_sys->integral_vs_sp_s_sz[i] );       free( sp_sys->integral_vs_sp_s_xxyy[i] );       free( sp_sys->integral_vs_sp_s_zz[i] );
        free( sp_sys->integral_vs_sp_s_x_sx[i] );     free( sp_sys->integral_vs_sp_s_x_xz[i] );
        free( sp_sys->integral_vs_sp_s_z_ss[i] );     free( sp_sys->integral_vs_sp_s_z_sz[i] );     free( sp_sys->integral_vs_sp_s_z_xxyy[i] );     free( sp_sys->integral_vs_sp_s_z_zz[i] );

        free( sp_sys->integral_vs_sp_s_xx_ss[i] );     free( sp_sys->integral_vs_sp_s_xx_sz[i] );     free( sp_sys->integral_vs_sp_s_xx_xx[i] );    free( sp_sys->integral_vs_sp_s_xx_yy[i] );   free( sp_sys->integral_vs_sp_s_xx_zz[i] );
        free( sp_sys->integral_vs_sp_s_xy_xy[i] );     
        free( sp_sys->integral_vs_sp_s_xz_sx[i] );      free( sp_sys->integral_vs_sp_s_xz_xz[i] );
        free( sp_sys->integral_vs_sp_s_zz_ss[i] );     free( sp_sys->integral_vs_sp_s_zz_sz[i] );     free( sp_sys->integral_vs_sp_s_zz_xxyy[i] );     free( sp_sys->integral_vs_sp_s_zz_zz[i] );
    }
    free( sp_sys->integral_vs_sp_s_ss );       free( sp_sys->integral_vs_sp_s_sz );       free( sp_sys->integral_vs_sp_s_xxyy );       free( sp_sys->integral_vs_sp_s_zz );
    free( sp_sys->integral_vs_sp_s_x_sx );     free( sp_sys->integral_vs_sp_s_x_xz );
    free( sp_sys->integral_vs_sp_s_z_ss );     free( sp_sys->integral_vs_sp_s_z_sz );     free( sp_sys->integral_vs_sp_s_z_xxyy );     free( sp_sys->integral_vs_sp_s_z_zz );

    free( sp_sys->integral_vs_sp_s_xx_ss );      free( sp_sys->integral_vs_sp_s_xx_sz );      free( sp_sys->integral_vs_sp_s_xx_xx );      free( sp_sys->integral_vs_sp_s_xx_yy );      free( sp_sys->integral_vs_sp_s_xx_zz );
    free( sp_sys->integral_vs_sp_s_xy_xy );
    free( sp_sys->integral_vs_sp_s_xz_sx );      free( sp_sys->integral_vs_sp_s_xz_xz );
    free( sp_sys->integral_vs_sp_s_zz_ss );      free( sp_sys->integral_vs_sp_s_zz_sz );      free( sp_sys->integral_vs_sp_s_zz_xxyy );      free( sp_sys->integral_vs_sp_s_zz_zz );


    #ifdef DEBUG
    printf(" CHK POINT - 2 \n");
    #endif

    // free knot
    free(sp_sys->integral_knot);


    #ifdef DEBUG
    printf(" CHK POINT - 3\n");
    #endif

    sp_sys->integral_lut = SP_SYSTEM_FALSE;

    MPI_Barrier(MPI_COMM_WORLD);

    return;
}



/// GA COST FUNCTION Calcualte Gnorm of the structure
/// This Method is assumed to be only called when all derivative calculation is done
double sp_cluster_system_get_gnorm( sp_cluster_system* sp_sys )
{
    // Gnorm is derined w.r.t. the last sp-lone pair ion centre
    double Return = 0.;
    double grad_x=0.;
    double grad_y=0.;
    double grad_z=0.;

    for(int i=0;i<sp_sys->number_of_classic_ion;i++)
    {
		// force acting on classic ions
		grad_x = pow(-sp_sys->classic_ion[i].elec_force_by_sp[0]-sp_sys->classic_ion[i].force_by_sp_core[0]-sp_sys->classic_ion[i].force_by_ion_core[0],2.);
		grad_y = pow(-sp_sys->classic_ion[i].elec_force_by_sp[1]-sp_sys->classic_ion[i].force_by_sp_core[1]-sp_sys->classic_ion[i].force_by_ion_core[1],2.);
		grad_z = pow(-sp_sys->classic_ion[i].elec_force_by_sp[2]-sp_sys->classic_ion[i].force_by_sp_core[2]-sp_sys->classic_ion[i].force_by_ion_core[2],2.);
		Return += (grad_x+grad_y+grad_z);
    }
    for(int i=0;i<sp_sys->number_of_sp_ion;i++)
    {
		// force acting on sp ions
		grad_x = pow(-sp_sys->sp_ion[i].elec_force_by_sp[0]-sp_sys->sp_ion[i].elec_force_by_ion[0]-sp_sys->sp_ion[i].force_by_sp_core[0]-sp_sys->sp_ion[i].force_by_sp_core[0],2.);
		grad_y = pow(-sp_sys->sp_ion[i].elec_force_by_sp[1]-sp_sys->sp_ion[i].elec_force_by_ion[1]-sp_sys->sp_ion[i].force_by_sp_core[1]-sp_sys->sp_ion[i].force_by_sp_core[1],2.);
		grad_z = pow(-sp_sys->sp_ion[i].elec_force_by_sp[2]-sp_sys->sp_ion[i].elec_force_by_ion[2]-sp_sys->sp_ion[i].force_by_sp_core[2]-sp_sys->sp_ion[i].force_by_sp_core[2],2.);
		Return += (grad_x+grad_y+grad_z);
    }
    Return = pow(Return,0.5)/((double)(sp_sys->number_of_classic_ion+sp_sys->number_of_sp_ion))/3.;
    return Return;
}

// CALCULATE EIGEN_SYSTEM OF SP_LONE PAIR SYSTEM
// Just Before the diagonalisation !
void sp_cluster_system_get_h_matrix_mpi( sp_cluster_system* sp_sys, int rank, int numtasks )
{   
    // MPI TREATMENT VARIABLES
    int L,R;                       // round-robin scheme index 1.
    int ista,iend, len;                 // round-robin scheme index 2.
    // H matrix temporal saving workspace
    gsl_matrix** h_matrix_tmp = NULL;
    gsl_matrix** h_matrix_tmp2= NULL;   // sp-e vs sp-core Coulomb buffer workspace
    // Monopolar Approximation WorkSpace
    gsl_vector* r = NULL;
    gsl_matrix* h_matrix_local = NULL;
    gsl_matrix* h_matrix_global_tmp = NULL;
    gsl_matrix* trans_matrix = NULL;                    // workspace for saving transformation matrix
    double tmp = 0.;    double charge_scaler;

    // Dipolar Approximation WorkSpace
    gsl_matrix* fx_matrix = NULL; gsl_matrix* fy_matrix = NULL; gsl_matrix* fz_matrix = NULL;
    double dxh = 0.; double dyh = 0.; double dzh = 0.;
    gsl_vector* global_dh = NULL;   gsl_vector* local_dh = NULL;

    gsl_matrix*** dummy_dh_matrix_x = NULL; gsl_matrix*** dummy_dh_matrix_y = NULL; gsl_matrix*** dummy_dh_matrix_z = NULL;

    double cs, cx, cy, cz;  // evec elem
    double dipx, dipy, dipz;
    int low_idx;
    int idx_sp1, idx_sp2, idx_cla;


    if( sp_sys->scf_cnt == 0 )
    {
        // memory allocation
        dummy_dh_matrix_x = (gsl_matrix***)malloc(sp_sys->number_of_sp_ion*sizeof(gsl_matrix**));
        dummy_dh_matrix_y = (gsl_matrix***)malloc(sp_sys->number_of_sp_ion*sizeof(gsl_matrix**));
        dummy_dh_matrix_z = (gsl_matrix***)malloc(sp_sys->number_of_sp_ion*sizeof(gsl_matrix**));

        for(int i=0;i<sp_sys->number_of_sp_ion;i++)
        {   dummy_dh_matrix_x[i] = (gsl_matrix**)malloc(sp_sys->number_of_sp_ion*sizeof(gsl_matrix*));
            dummy_dh_matrix_y[i] = (gsl_matrix**)malloc(sp_sys->number_of_sp_ion*sizeof(gsl_matrix*));
            dummy_dh_matrix_z[i] = (gsl_matrix**)malloc(sp_sys->number_of_sp_ion*sizeof(gsl_matrix*));  }
        for(int i=0;i<sp_sys->number_of_sp_ion;i++)
        {   for(int j=0;j<sp_sys->number_of_sp_ion;j++)
            {   dummy_dh_matrix_x[i][j] = gsl_matrix_calloc(4,4);
                dummy_dh_matrix_y[i][j] = gsl_matrix_calloc(4,4);
                dummy_dh_matrix_z[i][j] = gsl_matrix_calloc(4,4);       }}

        h_matrix_tmp  = (gsl_matrix**)malloc(sp_sys->number_of_sp_ion*sizeof(gsl_matrix*));
        h_matrix_tmp2 = (gsl_matrix**)malloc(sp_sys->number_of_sp_ion*sizeof(gsl_matrix*));
        for(int n=0;n<sp_sys->number_of_sp_ion;n++) 
        {   h_matrix_tmp[n]  = gsl_matrix_calloc(4,4);
            h_matrix_tmp2[n] = gsl_matrix_calloc(4,4);  }

        r = gsl_vector_calloc(3);               // workspace for a vector: sp_core -> external point charge (classical ion)
        h_matrix_local = gsl_matrix_calloc(4,4);  // workspace for temporal saving pairwise h_matrix w.r.t a classical ion in local symm
        h_matrix_global_tmp = gsl_matrix_calloc(4,4);   // workspace for temporal saving pairwise h_matrix in global symm

        fx_matrix = gsl_matrix_calloc(4,4);    fy_matrix = gsl_matrix_calloc(4,4);    fz_matrix = gsl_matrix_calloc(4,4); 
        global_dh = gsl_vector_calloc(4);
        local_dh  = gsl_vector_calloc(4);
        // memory allocation end


        ///     sp - point charge 
        ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///
        
        len = sp_sys->number_of_sp_ion*sp_sys->number_of_classic_ion;    // total number of for loop cycles
        // round robin algorithm
        L = len/numtasks;
        R = len%numtasks;
        ista = L*rank + MIN(R,rank);
        iend = ista + L - 1;
        if( R > rank )
            iend++;
        // round robin algorithm
        
        for(int n=ista;n<iend+1;n++)    /// decomposition by the number of sp-ions
        {   
            if( sp_sys->number_of_sp_ion == sp_sys->number_of_classic_ion )
            {   idx_sp1 = n/sp_sys->number_of_sp_ion;
                idx_cla = n%sp_sys->number_of_classic_ion;                                      } 
            else
            {   idx_sp1 = n/sp_sys->number_of_classic_ion;
                idx_cla = n-(n/sp_sys->number_of_classic_ion)*sp_sys->number_of_classic_ion;    }

            // r_classic_core - r_sp_core = r_sp_core_classic_core ... calc vector
            // core of ith classic ion: (gsl_vector*) sp_sys->classic_ion[i].core_position
            // core of sp-lone pair   : (gsl_vector*) sp_sys->sp_ion[n].core_position
            // using gsl support blas ... gsl_vector_sub( gsl_vector* a, gsl_vector* b ):  a <- a - b (NB. a is destroyed after operation)
            gsl_vector_memcpy( /*dest*/r,/*src*/sp_sys->classic_ion[idx_cla].core_position );
            gsl_vector_sub( r, sp_sys->sp_ion[idx_sp1].core_position ); // r = r_classic_core - r_sp_core(of nth core)
         
            trans_matrix = sp_cluster_support_get_transformation_matrix( r );   // get transformation matrix

            //if( sp_sys->classic_ion[idx_cla].charge_core < 0. ) // if the ion is anion
            
		// FIX this bit, shell -> then BM or RIM -> then BM ( if it is core and doesn't have shell )
            if( sp_sys->classic_ion[idx_cla].if_shell == SP_SYSTEM_TRUE \
		|| (sp_sys->classic_ion[idx_cla].if_shell == SP_SYSTEM_FALSE && sp_sys->classic_ion[idx_cla].if_has_shell == SP_SYSTEM_FALSE) )	// IF IT IS SHELL MODEL, THEN INCLUDE SHORT_RANGE INTERACTION 
            {   // Calculate H matrix in the transformed symmetry and save the result in 'h_matrix_local'
                gsl_matrix_set(h_matrix_local,0,0, sp_cluster_integrator_get_ch_11_element( &sp_sys->sp_ion[idx_sp1], &sp_sys->classic_ion[idx_cla] )
                        + sp_cluster_integrator_get_sh_11_element( sp_sys, &sp_sys->sp_ion[idx_sp1], &sp_sys->classic_ion[idx_cla] ));   // SS ... Long + Short
                gsl_matrix_set(h_matrix_local,0,3, sp_cluster_integrator_get_ch_14_element( &sp_sys->sp_ion[idx_sp1], &sp_sys->classic_ion[idx_cla] )
                        + sp_cluster_integrator_get_sh_14_element( sp_sys, &sp_sys->sp_ion[idx_sp1], &sp_sys->classic_ion[idx_cla] ));   // SZ ... Long + Short
                gsl_matrix_set(h_matrix_local,3,0, gsl_matrix_get(h_matrix_local,0,3)); // set m_03 == m_30
                gsl_matrix_set(h_matrix_local,1,1, sp_cluster_integrator_get_ch_2233_element( &sp_sys->sp_ion[idx_sp1], &sp_sys->classic_ion[idx_cla] )
                        + sp_cluster_integrator_get_sh_2233_element( sp_sys, &sp_sys->sp_ion[idx_sp1], &sp_sys->classic_ion[idx_cla] )); // XX ... Long + Short
                gsl_matrix_set(h_matrix_local,2,2, gsl_matrix_get(h_matrix_local,1,1)); // set m_11 == m_22
                gsl_matrix_set(h_matrix_local,3,3, sp_cluster_integrator_get_ch_44_element( &sp_sys->sp_ion[idx_sp1], &sp_sys->classic_ion[idx_cla] )
                        + sp_cluster_integrator_get_sh_44_element( sp_sys, &sp_sys->sp_ion[idx_sp1], &sp_sys->classic_ion[idx_cla] ));    // ZZ ... Long + Short
            }
            //else    // if it is cation 
            else if( sp_sys->classic_ion[idx_cla].if_shell == SP_SYSTEM_FALSE \
		&& sp_sys->classic_ion[idx_cla].if_has_shell == SP_SYSTEM_TRUE )	// IF IT IS RIM, THEN INCLUDE ONLY ELECTROSTATIC
            {   // Calculate H matrix in the transformed symmetry and save the result in 'h_matrix_local'
                gsl_matrix_set(h_matrix_local,0,0, sp_cluster_integrator_get_ch_11_element( &sp_sys->sp_ion[idx_sp1], &sp_sys->classic_ion[idx_cla] )
                        );   // SS ... Long
                gsl_matrix_set(h_matrix_local,0,3, sp_cluster_integrator_get_ch_14_element( &sp_sys->sp_ion[idx_sp1], &sp_sys->classic_ion[idx_cla] )
                        );   // SZ ... Long
                gsl_matrix_set(h_matrix_local,3,0, gsl_matrix_get(h_matrix_local,0,3)); // set m_03 == m_30
                gsl_matrix_set(h_matrix_local,1,1, sp_cluster_integrator_get_ch_2233_element( &sp_sys->sp_ion[idx_sp1], &sp_sys->classic_ion[idx_cla] )
                        ); // XX ... Long
                gsl_matrix_set(h_matrix_local,2,2, gsl_matrix_get(h_matrix_local,1,1)); // set m_11 == m_22
                gsl_matrix_set(h_matrix_local,3,3, sp_cluster_integrator_get_ch_44_element( &sp_sys->sp_ion[idx_sp1], &sp_sys->classic_ion[idx_cla] )
                        );    // ZZ ... Long + Short
            }
            // Exclude cat - cat short-range interaction

            // Transform the partial h matrix contribution in local -> global symm
            // 
            //  H_ij = T_ai T_bj H'_ab
            //
            //  recall:  r'_a = T_ai r_i   &&  r_i = T_ai r'_a  (Inverse Transformation of rank 1 tensor)
            //
            for(int alpha=0;alpha<4;alpha++)
            {   for(int beta=0;beta<4;beta++)
                {   for(int A=0;A<4;A++)
                        for(int B=0;B<4;B++)
                            tmp += gsl_matrix_get(trans_matrix,A,alpha)*gsl_matrix_get(trans_matrix,B,beta)*gsl_matrix_get(h_matrix_local,A,B);
                    gsl_matrix_set(h_matrix_global_tmp,alpha,beta,tmp);
                    tmp = 0.;
                }
            }

            gsl_matrix_add( h_matrix_tmp[idx_sp1], h_matrix_global_tmp );

        
            // refresh the workspace
            gsl_matrix_set_zero( h_matrix_local );
            gsl_matrix_set_zero( h_matrix_global_tmp );
            gsl_vector_set_zero( r );
            gsl_matrix_free( trans_matrix );

        }
        
        // Saving each sp_ions ('n'th) h_matrix w.r.t. classic ions -> scf_workspace (scf_h_matrix_vs_classic_ion[n])
        for(int n=0;n<sp_sys->number_of_sp_ion;n++)
            MPI_Allreduce( h_matrix_tmp[n]->data, sp_sys->scf_h_matrix_vs_classic_ion[n]->data, 16, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

        // Refreshing 'h_matrix_tmp[n] workspace'
        for(int n=0;n<sp_sys->number_of_sp_ion;n++)
            gsl_matrix_set_zero(h_matrix_tmp[n]);

        ///     sp - point charge - DONE 
        ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///




        // For loop work decomposition for sp - monopole interactions
        len = sp_sys->number_of_sp_ion*sp_sys->number_of_sp_ion;    // total number of for loop cycles
        // round robin algorithm
        L = len/numtasks;
        R = len%numtasks;
        ista = L*rank + MIN(R,rank);
        iend = ista + L - 1;
        if( R > rank )
            iend++;

        ///     sp - sp : Mono + Dipolar Approximation
        
        ///     ///     ///     ///     ///     ///    Monpolar Term Start
        for(int n=ista;n<iend+1;n++)
        {   
            idx_sp1 = n/sp_sys->number_of_sp_ion;   // for 'n'
            idx_sp2 = n%sp_sys->number_of_sp_ion;   // for 'm'    // this is beacuse 'n == m' always
                
            if( idx_sp1 != idx_sp2 ) // Excluding self-interaction
            {   // Monopolar Term Start
                gsl_vector_memcpy( r, sp_sys->sp_ion[idx_sp2].core_position );
                gsl_vector_sub( r, sp_sys->sp_ion[idx_sp1].core_position ); // r(sp1->sp2) = r_sp_core(of sp1th core) - r_sp_core(of sp2th core)
                trans_matrix = sp_cluster_support_get_transformation_matrix( r );   // get transformation matrix

                charge_scaler = (sp_sys->sp_ion[idx_sp2].charge_core+0.5*sp_sys->sp_ion[idx_sp2].charge_shell)/(sp_sys->sp_ion[idx_sp2].charge_shell);
                // Matrix Elem Estimation
                // For a Given state ('sp2'th sp-ion) determine a new state of 'sp1'th sp-ion
                // in  'sp_cluster_spsp_mono_integrator_get_ch_11_element(sp_k, spl)' like function the second arg 'spl' is the cenre
                // i.e., spl -> spk
                // This is beacuse, the transformation matrix is defined above w.r.t. sp1 ion centre
                gsl_matrix_set(h_matrix_local,0,0, charge_scaler*sp_cluster_spsp_mono_integrator_get_ch_11_element( &sp_sys->sp_ion[idx_sp2], &sp_sys->sp_ion[idx_sp1] )
                      + 0.5*sp_cluster_spsp_mono_integrator_get_sh_11_element( sp_sys, &sp_sys->sp_ion[idx_sp2], &sp_sys->sp_ion[idx_sp1] ));   // SS ... Long + Short
                gsl_matrix_set(h_matrix_local,0,3, charge_scaler*sp_cluster_spsp_mono_integrator_get_ch_14_element( &sp_sys->sp_ion[idx_sp2], &sp_sys->sp_ion[idx_sp1] )
                      + 0.5*sp_cluster_spsp_mono_integrator_get_sh_14_element( sp_sys, &sp_sys->sp_ion[idx_sp2], &sp_sys->sp_ion[idx_sp1] ));   // SZ ... Long + Short
                gsl_matrix_set(h_matrix_local,3,0, gsl_matrix_get(h_matrix_local,0,3)); // set m_03 == m_30
                gsl_matrix_set(h_matrix_local,1,1, charge_scaler*sp_cluster_spsp_mono_integrator_get_ch_2233_element( &sp_sys->sp_ion[idx_sp2], &sp_sys->sp_ion[idx_sp1] )
                      + 0.5*sp_cluster_spsp_mono_integrator_get_sh_2233_element( sp_sys, &sp_sys->sp_ion[idx_sp2], &sp_sys->sp_ion[idx_sp1] )); // XX ... Long + Short
                gsl_matrix_set(h_matrix_local,2,2, gsl_matrix_get(h_matrix_local,1,1)); // set m_11 == m_22
                gsl_matrix_set(h_matrix_local,3,3, charge_scaler*sp_cluster_spsp_mono_integrator_get_ch_44_element( &sp_sys->sp_ion[idx_sp2], &sp_sys->sp_ion[idx_sp1] )
                      + 0.5*sp_cluster_spsp_mono_integrator_get_sh_44_element( sp_sys, &sp_sys->sp_ion[idx_sp2], &sp_sys->sp_ion[idx_sp1] ));    // ZZ ... Long + Short
                // Transform the partial h matrix contribution in local -> global symm
                // 
                //  H_ij = T_ai T_bj H'_ab
                //
                //  recall:  r'_a = T_ai r_i   &&  r_i = T_ai r'_a  (Inverse Transformation of rank 1 tensor)
                //
                for(int alpha=0;alpha<4;alpha++)
                {   for(int beta=0;beta<4;beta++)
                    {   for(int A=0;A<4;A++)
                            for(int B=0;B<4;B++)
                                tmp += gsl_matrix_get(trans_matrix,A,alpha)*gsl_matrix_get(trans_matrix,B,beta)*gsl_matrix_get(h_matrix_local,A,B);
                        gsl_matrix_set(h_matrix_global_tmp,alpha,beta,tmp);
                        tmp = 0.;              
                    }
                }
                // Refresh Monopolar related workspace
                //gsl_matrix_add( sp_sys->sp_ion[n].h_matrix, h_matrix_global_tmp );    // save the contribution to h_matrix (in sp_ion object)
                gsl_matrix_add( h_matrix_tmp[idx_sp1] , h_matrix_global_tmp );    // save the contribution to h_matrix (in sp_ion object)
                /*
                 * At the start of this if( idx_sp1 != idx_sp2 ) statement right above
                 * the symmetry is following " idx_sp1 -> idx_sp2 " the standard in symm is sp1
                 * 
                 * thats why 'h_mat_global_tmp' is added to h_matrix_tm[idx_sp1] but Not ..[idx_sp2]
                 */

                gsl_matrix_set_zero(h_matrix_global_tmp);
                gsl_matrix_set_zero(h_matrix_local);
                // Monopolar Term End

                charge_scaler = sp_sys->sp_ion[idx_sp2].charge_core/(sp_sys->sp_ion[idx_sp2].charge_shell);
                // Scaling factor to make sp_core_charge*sp_shell_charge


                /// Estimating sp-e vs sp-core interaction (only Coulomb);
                gsl_matrix_set(h_matrix_local,0,0, charge_scaler*sp_cluster_spsp_mono_integrator_get_ch_11_element( &sp_sys->sp_ion[idx_sp2], &sp_sys->sp_ion[idx_sp1] ));
                gsl_matrix_set(h_matrix_local,0,3, charge_scaler*sp_cluster_spsp_mono_integrator_get_ch_14_element( &sp_sys->sp_ion[idx_sp2], &sp_sys->sp_ion[idx_sp1] ));
                gsl_matrix_set(h_matrix_local,3,0, gsl_matrix_get(h_matrix_local,0,3)); // set m_03 == m_30
                gsl_matrix_set(h_matrix_local,1,1, charge_scaler*sp_cluster_spsp_mono_integrator_get_ch_2233_element( &sp_sys->sp_ion[idx_sp2], &sp_sys->sp_ion[idx_sp1] ));
                gsl_matrix_set(h_matrix_local,2,2, gsl_matrix_get(h_matrix_local,1,1)); // set m_11 == m_22
                gsl_matrix_set(h_matrix_local,3,3, charge_scaler*sp_cluster_spsp_mono_integrator_get_ch_44_element( &sp_sys->sp_ion[idx_sp2], &sp_sys->sp_ion[idx_sp1] ));

                for(int alpha=0;alpha<4;alpha++)
                {   for(int beta=0;beta<4;beta++)
                    {   for(int A=0;A<4;A++)
                            for(int B=0;B<4;B++)
                                tmp += gsl_matrix_get(trans_matrix,A,alpha)*gsl_matrix_get(trans_matrix,B,beta)*gsl_matrix_get(h_matrix_local,A,B);
                        gsl_matrix_set(h_matrix_global_tmp,alpha,beta,tmp);
                        tmp = 0.;              
                    }
                }
                gsl_matrix_add( h_matrix_tmp2[idx_sp1] , h_matrix_global_tmp );    // save the contribution to h_matrix (in sp_ion object)
                gsl_matrix_set_zero(h_matrix_global_tmp);
                gsl_matrix_set_zero(h_matrix_local);
                // sp-e vs sp-core term END

                // Refresh Transformation workspace
                gsl_vector_set_zero(r);     gsl_matrix_free(trans_matrix);
            }

        }
        // End of sp-sp monopole h_matrix estimation    // End of sp-e vs sp-c Coulomb h_matrix estimation
        
        // Save the 'n'th sp-ions monopole interactions into scf workspace ( scf_h_matrix_vs_sp_ion_monopole[n] );
        for(int n=0;n<sp_sys->number_of_sp_ion;n++)
        {   MPI_Allreduce( h_matrix_tmp[n]->data, sp_sys->scf_h_matrix_vs_sp_ion_monopole[n]->data, 16, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce( h_matrix_tmp2[n]->data, sp_sys->scf_h_matrix_vs_sp_core[n]->data,16,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);               }
 
        // Refreshing h_matrix_tmp workspace
        for(int n=0;n<sp_sys->number_of_sp_ion;n++)
            gsl_matrix_set_zero( h_matrix_tmp[n] );

        // Setting 'n'th sp_ion onsite term
        for(int n=0;n<sp_sys->number_of_sp_ion;n++)
        {   gsl_matrix_set(sp_sys->scf_h_matrix_vs_sp_ion_onsite[n],3,3,sp_sys->sp_ion[n].esp);
            gsl_matrix_set(sp_sys->scf_h_matrix_vs_sp_ion_onsite[n],2,2,sp_sys->sp_ion[n].esp);
            gsl_matrix_set(sp_sys->scf_h_matrix_vs_sp_ion_onsite[n],1,1,sp_sys->sp_ion[n].esp);     }

        /// End Of vs Classic/ vs Sp-monopole / vs Sp-onsite


        ///     ///     ///     ///     ///     ///    Dipolar Term Start
        for(int n=ista;n<iend+1;n++)
        {   
            idx_sp1 = n/sp_sys->number_of_sp_ion;   // for 'n'
            idx_sp2 = n%sp_sys->number_of_sp_ion;   // for 'm'    // this is beacuse 'n == m' always
                
            if( idx_sp1 != idx_sp2 ) // Excluding self-interaction
            {   
                gsl_vector_memcpy( r, sp_sys->sp_ion[idx_sp2].core_position );
                gsl_vector_sub( r, sp_sys->sp_ion[idx_sp1].core_position ); // r(sp1->sp2) = r_sp2_core - r_sp1_core
                trans_matrix = sp_cluster_support_get_transformation_matrix( r );   // get transformation matrix w.r.t the above defined geometry
                // Likewise sp1 is seeing sp2 ... i.e., sp1 -> sp2

                // !!! set fx matrix
                gsl_matrix_set(fx_matrix,0,1,0.5*(sp_cluster_spsp_di_integrator_get_ch_x_12_element(&sp_sys->sp_ion[idx_sp2],&sp_sys->sp_ion[idx_sp1]) 
                        + sp_cluster_spsp_di_integrator_get_sh_x_12_element( sp_sys, &sp_sys->sp_ion[idx_sp2],&sp_sys->sp_ion[idx_sp1])) );       // set fx 12 element
                gsl_matrix_set(fx_matrix,1,3,0.5*(sp_cluster_spsp_di_integrator_get_ch_x_24_element(&sp_sys->sp_ion[idx_sp2],&sp_sys->sp_ion[idx_sp1])
                        + sp_cluster_spsp_di_integrator_get_sh_x_24_element( sp_sys, &sp_sys->sp_ion[idx_sp2],&sp_sys->sp_ion[idx_sp1])) );       // set fx 24 element
                gsl_matrix_set(fx_matrix,1,0,gsl_matrix_get(fx_matrix,0,1)); gsl_matrix_set(fx_matrix,3,1,gsl_matrix_get(fx_matrix,1,3));           // SIGN IS NOT INVERSED YET
                // !!! set fy matrix
                gsl_matrix_set(fy_matrix,0,2,gsl_matrix_get(fx_matrix,0,1));    // note that DxH SX == DyH SY
                gsl_matrix_set(fy_matrix,2,3,gsl_matrix_get(fx_matrix,1,3));    // note that DxH XZ == DyH YZ
                gsl_matrix_set(fy_matrix,2,0,gsl_matrix_get(fy_matrix,0,2)); gsl_matrix_set(fy_matrix,3,2,gsl_matrix_get(fy_matrix,2,3));
                // !!! set fz matrix
                gsl_matrix_set(fz_matrix,0,0,0.5*(sp_cluster_spsp_di_integrator_get_ch_z_11_element(&sp_sys->sp_ion[idx_sp2],&sp_sys->sp_ion[idx_sp1])
                        + sp_cluster_spsp_di_integrator_get_sh_z_11_element( sp_sys, &sp_sys->sp_ion[idx_sp2],&sp_sys->sp_ion[idx_sp1])) );       // set fz 11 element
                gsl_matrix_set(fz_matrix,0,3,0.5*(sp_cluster_spsp_di_integrator_get_ch_z_14_element(&sp_sys->sp_ion[idx_sp2],&sp_sys->sp_ion[idx_sp1])
                        + sp_cluster_spsp_di_integrator_get_sh_z_14_element( sp_sys, &sp_sys->sp_ion[idx_sp2],&sp_sys->sp_ion[idx_sp1])) );       // set fz 14 element
                gsl_matrix_set(fz_matrix,3,0,gsl_matrix_get(fz_matrix,0,3));                                        // sz == zs
                gsl_matrix_set(fz_matrix,1,1,0.5*(sp_cluster_spsp_di_integrator_get_ch_z_2233_element(&sp_sys->sp_ion[idx_sp2],&sp_sys->sp_ion[idx_sp1])
                        + sp_cluster_spsp_di_integrator_get_sh_z_2233_element( sp_sys, &sp_sys->sp_ion[idx_sp2],&sp_sys->sp_ion[idx_sp1])) );     // set fz 22 element
                gsl_matrix_set(fz_matrix,2,2,gsl_matrix_get(fz_matrix,1,1));                                        // xx == yy
                gsl_matrix_set(fz_matrix,3,3,0.5*(sp_cluster_spsp_di_integrator_get_ch_z_44_element(&sp_sys->sp_ion[idx_sp2],&sp_sys->sp_ion[idx_sp1])
                        + sp_cluster_spsp_di_integrator_get_sh_z_44_element( sp_sys, &sp_sys->sp_ion[idx_sp2],&sp_sys->sp_ion[idx_sp1])) );       // set fz 44 element

                // Caculating global fxyz matrices   wrt  ith classic ion
                for(int m=0;m<4;m++)
                {   for(int n=0;n<4;n++)
                    {   for(int h=0;h<4;h++)
                        {   for(int l=0;l<4;l++)
                            {   dxh += gsl_matrix_get(trans_matrix,h,m)*gsl_matrix_get(trans_matrix,l,n)*gsl_matrix_get(fx_matrix,h,l); 
                                dyh += gsl_matrix_get(trans_matrix,h,m)*gsl_matrix_get(trans_matrix,l,n)*gsl_matrix_get(fy_matrix,h,l); 
                                dzh += gsl_matrix_get(trans_matrix,h,m)*gsl_matrix_get(trans_matrix,l,n)*gsl_matrix_get(fz_matrix,h,l); 
                            }   // Calculate local dH matrix element
                                //
                                //  < m | d H | n > == T_hm T_ln < h' | d'H | l' >
                                //
                        }       // Get local expression for < m | d H | n >
                        gsl_vector_set(local_dh,1,dxh); gsl_vector_set(local_dh,2,dyh); gsl_vector_set(local_dh,3,dzh);
                        dxh = 0.; dyh = 0.; dzh = 0.;
                        // Inverse transformation to global symm ...  
                        gsl_blas_dgemv(CblasTrans,1.,trans_matrix,local_dh,0.,global_dh);
                    
                        gsl_matrix_set(dummy_dh_matrix_x[idx_sp1][idx_sp2],m,n,gsl_vector_get(global_dh,1));
                        gsl_matrix_set(dummy_dh_matrix_y[idx_sp1][idx_sp2],m,n,gsl_vector_get(global_dh,2));
                        gsl_matrix_set(dummy_dh_matrix_z[idx_sp1][idx_sp2],m,n,gsl_vector_get(global_dh,3));      // Global dH matrix is saved here 'n'th sp ion w.r.t 'm'th sp ion
        
                        /* 
                         * The data structure of dummy_dh_matrix_x,y and z 
                         *
                         * 'sp1' is at row and 'sp2' is at col
                         *
                         * e.g., row2, col4 -> (2,4) means second sp ion is pointing forth sp ion
                         *
                         * i.e., sp2 -> sp4
                         *
                         * Likewise (4,2) means sp4 -> sp2
                         */

                        gsl_vector_set_zero(local_dh); gsl_vector_set_zero(global_dh);
                    }
                }   // The section Right Below can be improved by the obtained global < m | d H | n >

                gsl_matrix_set_zero(fx_matrix); gsl_matrix_set_zero(fy_matrix); gsl_matrix_set_zero(fz_matrix);
                gsl_vector_set_zero(r);     gsl_matrix_free(trans_matrix);
            }
        }
        //MPI_Barrier(MPI_COMM_WORLD);
        // End of scf_dh_x,y and z_matrix_vs_sp_dipole end 
        // Note that here the workspace buffers are dummy_dh_matrix_x,y and z
        for(int n=0;n<sp_sys->number_of_sp_ion;n++)
        {   for(int m=0;m<sp_sys->number_of_sp_ion;m++)
            {   if( n != m )
                {   MPI_Allreduce( dummy_dh_matrix_x[n][m]->data, sp_sys->scf_dh_x_matrix_vs_sp_ion_dipole[n][m]->data, 16, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
                    MPI_Allreduce( dummy_dh_matrix_y[n][m]->data, sp_sys->scf_dh_y_matrix_vs_sp_ion_dipole[n][m]->data, 16, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
                    MPI_Allreduce( dummy_dh_matrix_z[n][m]->data, sp_sys->scf_dh_z_matrix_vs_sp_ion_dipole[n][m]->data, 16, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );  }}}
        // Save into scf workspace ... 'scf_dh_x,y and z_matrix_vs_sp_ion_dipole
        //MPI_Barrier(MPI_COMM_WORLD);
        
        // Memory Detach
        gsl_vector_free(global_dh); gsl_vector_free(local_dh);
        gsl_matrix_free(fx_matrix);         gsl_matrix_free(fy_matrix);         gsl_matrix_free(fz_matrix);
        gsl_matrix_free(h_matrix_global_tmp);
        gsl_matrix_free(h_matrix_local);
        gsl_vector_free(r);
        for(int n=0;n<sp_sys->number_of_sp_ion;n++) 
        {   gsl_matrix_free(h_matrix_tmp[n]);
            gsl_matrix_free(h_matrix_tmp2[n]);      }
        free(h_matrix_tmp);     // tmp matrix mem-detach
        free(h_matrix_tmp2);
        //free dh matrices
        
        for(int n=0;n<sp_sys->number_of_sp_ion;n++)
        {   for(int m=0;m<sp_sys->number_of_sp_ion;m++)
            {   gsl_matrix_free( dummy_dh_matrix_x[n][m] );
                gsl_matrix_free( dummy_dh_matrix_y[n][m] );
                gsl_matrix_free( dummy_dh_matrix_z[n][m] );     }}
        for(int n=0;n<sp_sys->number_of_sp_ion;n++)
        {   free( dummy_dh_matrix_x[n] );
            free( dummy_dh_matrix_y[n] );
            free( dummy_dh_matrix_z[n] );   }
        free(dummy_dh_matrix_x);    free(dummy_dh_matrix_y);    free(dummy_dh_matrix_z);
        // Memory Detach End
   


        ///     ///     ///     ///     ///     ///     ///     ///     ///
        // H_MAT ESTIMATION
        ///     ///     ///     ///     ///     ///     ///     ///     ///
        //
        //
        for(int n=0;n<sp_sys->number_of_sp_ion;n++)
        {  
            // refresh sp_ion[n].h_matrix into zero
            gsl_matrix_set_zero(sp_sys->sp_ion[n].h_matrix);
            // add h contribution by classic ions
            gsl_matrix_add(sp_sys->sp_ion[n].h_matrix,sp_sys->scf_h_matrix_vs_classic_ion[n]);
            // add h contribution by onsite
            gsl_matrix_add(sp_sys->sp_ion[n].h_matrix,sp_sys->scf_h_matrix_vs_sp_ion_onsite[n]);
            // add h contribution by sp ions monopole
            gsl_matrix_add(sp_sys->sp_ion[n].h_matrix,sp_sys->scf_h_matrix_vs_sp_ion_monopole[n]);


            for(int m=0;m<sp_sys->number_of_sp_ion;m++)
            {
                if( n != m )
                {
                    // Estimate Dipolar H Matrix 
                    /*
                     *       < alpha_n| { 2Cm1Cm2< Sm | rx | pmx > (-dxh) + 2Cm1Cm3< Sm | ry | Pmy > (-dyh) + 2Cm1Cm3< Sm | rz | Pmz > (-dzh) } | beta_n >
                     */
                    
                    // Get evec of 'm'th sp-ion
                    // This is because, at the moment we are estimating 'n' ions h_matrix to determine its eigensystem w.r.t 'm' ions density(eigenstate)
                    low_idx = sp_cluster_support_get_lowest_state( sp_sys->sp_ion[m].eigen_value ); // get lowest eval index
                    cs=gsl_matrix_get(sp_sys->sp_ion[m].eigen_vector,0,low_idx);    cx=gsl_matrix_get(sp_sys->sp_ion[m].eigen_vector,1,low_idx);
                    cy=gsl_matrix_get(sp_sys->sp_ion[m].eigen_vector,2,low_idx);    cz=gsl_matrix_get(sp_sys->sp_ion[m].eigen_vector,3,low_idx);
                    // Pre-Calc-Leading-Terms 2*CS*CX<S|rx|PX> ....
                    dipx = 2.*cs*cx*sp_cluster_integrator_get_x_12(&sp_sys->sp_ion[m]); dipy = 2.*cs*cy*sp_cluster_integrator_get_y_13(&sp_sys->sp_ion[m]); dipz = 2.*cs*cz*sp_cluster_integrator_get_z_14(&sp_sys->sp_ion[m]);

                    // Note that Charge*Charge is alrdy pre-calc-ed when fx fy fz matriced estimation
                    
                    // Dipolar Term addition operation
                    // When 'n'th ion -> 'm'th ion
                    for(int ii=0;ii<4;ii++)
                    {   for(int jj=ii;jj<4;jj++)
                        {   
                            gsl_matrix_set(sp_sys->sp_ion[n].h_matrix,ii,jj,gsl_matrix_get(sp_sys->sp_ion[n].h_matrix,ii,jj)
                                    - dipx*gsl_matrix_get(sp_sys->scf_dh_x_matrix_vs_sp_ion_dipole[n][m],ii,jj) 
                                    - dipy*gsl_matrix_get(sp_sys->scf_dh_y_matrix_vs_sp_ion_dipole[n][m],ii,jj) 
                                    - dipz*gsl_matrix_get(sp_sys->scf_dh_z_matrix_vs_sp_ion_dipole[n][m],ii,jj) );
                            if( ii != jj )
                                gsl_matrix_set(sp_sys->sp_ion[n].h_matrix,jj,ii, gsl_matrix_get(sp_sys->sp_ion[n].h_matrix,ii,jj));
                        }
                    }
                }// End of If

            }
        }
    }
    else    // if its not the first scf cycle ... keep using the pre-estimated data ( from the first scf iteration )
    {
        ///     ///     ///     ///     ///     ///     ///     ///     ///
        // H_MAT ESTIMATION
        ///     ///     ///     ///     ///     ///     ///     ///     ///
        //
        //
        for(int n=0;n<sp_sys->number_of_sp_ion;n++)
        {  
            // refresh sp_ion[n].h_matrix into zero
            gsl_matrix_set_zero(sp_sys->sp_ion[n].h_matrix);
            // add h contribution by classic ions
            gsl_matrix_add(sp_sys->sp_ion[n].h_matrix,sp_sys->scf_h_matrix_vs_classic_ion[n]);
            // add h contribution by onsite
            gsl_matrix_add(sp_sys->sp_ion[n].h_matrix,sp_sys->scf_h_matrix_vs_sp_ion_onsite[n]);
            // add h contribution by sp ions monopole
            gsl_matrix_add(sp_sys->sp_ion[n].h_matrix,sp_sys->scf_h_matrix_vs_sp_ion_monopole[n]);
    
            for(int m=0;m<sp_sys->number_of_sp_ion;m++)
            {
                if( n != m )
                {
                    // Estimate Dipolar H Matrix 
                    /*
                     *       < alpha_n| { 2Cm1Cm2< Sm | rx | pmx > (-dxh) + 2Cm1Cm3< Sm | ry | Pmy > (-dyh) + 2Cm1Cm3< Sm | rz | Pmz > (-dzh) } | beta_n >
                     */
                    
                    // Get evec of 'm'th sp-ion
                    // This is because, at the moment we are estimating 'n' ions h_matrix to determine its eigensystem w.r.t 'm' ions density(eigenstate)
                    low_idx = sp_cluster_support_get_lowest_state( sp_sys->sp_ion[m].eigen_value ); // get lowest eval index
                    cs=gsl_matrix_get(sp_sys->sp_ion[m].eigen_vector,0,low_idx);    cx=gsl_matrix_get(sp_sys->sp_ion[m].eigen_vector,1,low_idx);
                    cy=gsl_matrix_get(sp_sys->sp_ion[m].eigen_vector,2,low_idx);    cz=gsl_matrix_get(sp_sys->sp_ion[m].eigen_vector,3,low_idx);
                    // Pre-Calc-Leading-Terms 2*CS*CX<S|rx|PX> ....
                    dipx = 2.*cs*cx*sp_cluster_integrator_get_x_12(&sp_sys->sp_ion[m]); dipy = 2.*cs*cy*sp_cluster_integrator_get_y_13(&sp_sys->sp_ion[m]); dipz = 2.*cs*cz*sp_cluster_integrator_get_z_14(&sp_sys->sp_ion[m]);

                    //printf("dip r\t\t:\t %lf\t%lf\t%lf\n",dipx/2./cs/cx,dipy/2./cs/cy,dipz/2./cs/cz);

                    // Note that Charge*Charge is alrdy pre-calc-ed when fx fy fz matriced estimation
                    
                    // Dipolar Term addition operation
                    // When 'n'th ion -> 'm'th ion
                    for(int ii=0;ii<4;ii++)
                    {   for(int jj=ii;jj<4;jj++)
                        {   
                            gsl_matrix_set(sp_sys->sp_ion[n].h_matrix,ii,jj,gsl_matrix_get(sp_sys->sp_ion[n].h_matrix,ii,jj)
                                    - dipx*gsl_matrix_get(sp_sys->scf_dh_x_matrix_vs_sp_ion_dipole[n][m],ii,jj) 
                                    - dipy*gsl_matrix_get(sp_sys->scf_dh_y_matrix_vs_sp_ion_dipole[n][m],ii,jj) 
                                    - dipz*gsl_matrix_get(sp_sys->scf_dh_z_matrix_vs_sp_ion_dipole[n][m],ii,jj) );

                            if( ii != jj )
                                gsl_matrix_set(sp_sys->sp_ion[n].h_matrix,jj,ii, gsl_matrix_get(sp_sys->sp_ion[n].h_matrix,ii,jj));
                        }
                    }
                }// End of If

            }
        }

    }// END OF IF & ELSE


    return;
}


// DIALGONALISATION MODULE
// Diagonalise the sp-lone pair h matrix of the 'labeled' ion
void sp_cluster_system_get_eigensystem( sp_cluster_system* sp_sys )
{   //USING QR DECOMPO BASED GSL SUPPORT ALGORITHM
    gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(4);     // 4x4 h matrix ( real symmetric )
    gsl_matrix* h_tmp = gsl_matrix_calloc(4,4);                 // dummy h_matrix workspace, take it up for the matrix destruction
    // CALL DIAGONALISATION ROUTINE
    // Syntax: gsl_eigen_symmv( gsl_matrix* A, gsl_vector* eval, gsl_matrix* evec, gsl_eigen_symmv_workspace* w )
    // A: Target matrix to diagonalise
    // eval: type gsl_vector saves eigenvalues unordered
    // evec: type gsl_matrix saves corresponding eigenvector (mutually orthogonal and normalised) as a column vector
    for(int n=0;n<sp_sys->number_of_sp_ion;n++)
    {   gsl_matrix_memcpy( /*dest*/ h_tmp, sp_sys->sp_ion[n].h_matrix );
        gsl_eigen_symmv( h_tmp, sp_sys->sp_ion[n].eigen_value, sp_sys->sp_ion[n].eigen_vector, w );     }
    gsl_matrix_free(h_tmp);
    gsl_eigen_symmv_free(w);
    return;
}

// This is a function for updating sp-derived energy
// must be called after scf convergence achieved !
void sp_cluster_system_sp_energy_update( sp_cluster_system* sp_sys )
{   
    const int number_of_sp_ion = sp_sys->number_of_sp_ion;

    int low_idx;
    double evec[4] = {0.,0.,0.,0.};
    double sp_energy_tmp = 0.;

    for(int n=0;n<number_of_sp_ion;n++)
    {   
        low_idx = sp_cluster_support_get_lowest_state(sp_sys->sp_ion[n].eigen_value);

        evec[0] = gsl_matrix_get(sp_sys->sp_ion[n].eigen_vector,0,low_idx);
        evec[1] = gsl_matrix_get(sp_sys->sp_ion[n].eigen_vector,1,low_idx);
        evec[2] = gsl_matrix_get(sp_sys->sp_ion[n].eigen_vector,2,low_idx);
        evec[3] = gsl_matrix_get(sp_sys->sp_ion[n].eigen_vector,3,low_idx);     // get evec

        // Calculate Energy
        for(int k=0;k<4;k++)
        {   for(int l=0;l<4;l++)
                sp_energy_tmp += evec[k]*evec[l]*gsl_matrix_get(sp_sys->sp_ion[n].h_matrix,k,l);
        }

        sp_sys->sp_ion[n].sp_energy = sp_energy_tmp;
        sp_energy_tmp = 0.;
    }

    return;
}

int sp_cluster_system_set_is_scf_done( sp_cluster_system* sp_sys, int t_or_f )
{
    sp_sys->is_scf_done = t_or_f;
    return t_or_f;
}

// SCF CALL ADDED 07192019
int sp_cluster_system_scf_mpi( sp_cluster_system* sp_sys, const int rank, const int numtasks )
{
    double prev_eval[sp_sys->number_of_sp_ion];
    double convergence_checker;
    int ret = SP_SYSTEM_FALSE;
    const int max_scf = MAX_SCF_CYCLE;

	// damping work space variables
	int damping_low_state_index;
	double damping_prev_evec[sp_sys->number_of_sp_ion][4];
	double damping_norm;
	for(int i=0;i<sp_sys->number_of_sp_ion;i++)	memset(damping_prev_evec[i],0.,4*sizeof(double));
	// damping work space variables - end

	// recording initial evec info
	for(int n=0;n<sp_sys->number_of_sp_ion;n++)
	{
		damping_low_state_index = sp_cluster_support_get_lowest_state( sp_sys->sp_ion[n].eigen_value );
		for(int k=0;k<4;k++)
			damping_prev_evec[n][k] = gsl_matrix_get( sp_sys->sp_ion[n].eigen_vector, k, damping_low_state_index );
	}
	// end of recording initial evec info

    // Refresh SCF WorkSpace    
    for(int i=0;i<sp_sys->number_of_sp_ion;i++) 
    {   gsl_matrix_set_zero(sp_sys->scf_h_matrix_vs_classic_ion[i]);
        gsl_matrix_set_zero(sp_sys->scf_h_matrix_vs_sp_ion_onsite[i]);
        gsl_matrix_set_zero(sp_sys->scf_h_matrix_vs_sp_ion_monopole[i]);
        gsl_matrix_set_zero(sp_sys->scf_h_matrix_vs_sp_core[i]);
    }
    for(int i=0;i<sp_sys->number_of_sp_ion;i++)
    {   for(int j=0;j<sp_sys->number_of_sp_ion;j++)
        {   gsl_matrix_set_zero(sp_sys->scf_dh_x_matrix_vs_sp_ion_dipole[i][j]);
            gsl_matrix_set_zero(sp_sys->scf_dh_x_matrix_vs_sp_ion_dipole[i][j]);
            gsl_matrix_set_zero(sp_sys->scf_dh_x_matrix_vs_sp_ion_dipole[i][j]); }}
    sp_sys->scf_cnt = 0;


    if( sp_sys->is_scf_done == SP_SYSTEM_FALSE )
    {
        // Reset SP_ION_EVEC
        if( sp_sys->if_first_scf_trial == SP_SYSTEM_TRUE ) 
        {   
		/*
            if( rank == 0 )
            {   
		printf("\n");
                printf(" Initialising First SCF Cycle ! \n");
            }
		*/
            for(int i=0;i<sp_sys->number_of_sp_ion;i++)
            {
                gsl_matrix_set_zero(sp_sys->sp_ion[i].eigen_vector);
                gsl_vector_set_zero(sp_sys->sp_ion[i].eigen_value);

                for(int j=0;j<4;j++)
                    gsl_matrix_set(sp_sys->sp_ion[i].eigen_vector,0,j,1.);   // Set All with pure s state
                gsl_vector_set(sp_sys->sp_ion[i].eigen_value,0,-1.);
            }
        
            sp_sys->if_first_scf_trial = SP_SYSTEM_FALSE;
        }
	// Reset done ///	///	///	///	///	///	///	///	///	///	///	///	///	///
	
        for(int k=0;k<max_scf;k++)
        {   // Refresh Convergence Checker
            convergence_checker = 0.;
            // Get old eigen_values
            for(int n=0;n<sp_sys->number_of_sp_ion;n++)
                prev_eval[n] = gsl_vector_get(sp_sys->sp_ion[n].eigen_value, sp_cluster_support_get_lowest_state(sp_sys->sp_ion[n].eigen_value));
            // Estimate Eval/Evec
     
            sp_cluster_system_get_h_matrix_mpi(sp_sys,rank,numtasks);
            sp_cluster_system_get_eigensystem(sp_sys);

            // If the number of sp ion is one !!!
            if( sp_sys->number_of_sp_ion == 1 )
            {   ret = SP_SYSTEM_TRUE;
                sp_sys->is_scf_done = SP_SYSTEM_FALSE;   // Turn Flag On; SCF IS DONE               
                sp_cluster_system_sp_energy_update(sp_sys);

                MPI_Barrier(MPI_COMM_WORLD);
                return ret;                     
            }
            // Single sp ion case

            sp_sys->scf_cnt++;

            //if(rank==0) printf("SCF CYCLE %d:\t\t",k+1);

            // SCF Convergence Test by evaluating NewEval - PrevEval
            for(int n=0;n<sp_sys->number_of_sp_ion;n++)
            {    
                convergence_checker += pow(( gsl_vector_get(sp_sys->sp_ion[n].eigen_value,sp_cluster_support_get_lowest_state(sp_sys->sp_ion[n].eigen_value)) - prev_eval[n] ), 2.);
                
                //TEMPOLAR PT.
                //if(rank==0) printf("%.8lf\t", gsl_vector_get(sp_sys->sp_ion[n].eigen_value, sp_cluster_support_get_lowest_state(sp_sys->sp_ion[n].eigen_value) ) );

            }
            //if(rank==0) printf("\n");

            convergence_checker = sqrt(convergence_checker)/sp_sys->number_of_sp_ion;

            if( convergence_checker < SP_SYSTEM_SCF_TOL )
            {   ret = SP_SYSTEM_TRUE;
                sp_sys->is_scf_done = SP_SYSTEM_TRUE;
                // ONCE SCF IS ACHIEVED ENERGY UPDATE NEEDS TO BE CARRIED OUT !
                sp_cluster_system_sp_energy_update(sp_sys);

                MPI_Barrier(MPI_COMM_WORLD);

                return ret;
            }
            else
                continue;
        }
	// if the line below is read -> SCF is failed .... implementation of damping method is necessary
    if( rank == 0 ) 
	{	//printf("SCF CONVERGENCE FAILED ... ------> ");
		//printf("DAMPING AUXILLARY ROUTINE IS CALLED ....\n");
	}
	
	// start from initial evec 
	for(int n=0;n<sp_sys->number_of_sp_ion;n++)
	{
		damping_low_state_index = sp_cluster_support_get_lowest_state( sp_sys->sp_ion[n].eigen_value );
		
		for(int k=0;k<4;k++)
			gsl_matrix_set( sp_sys->sp_ion[n].eigen_vector, k, damping_low_state_index, damping_prev_evec[n][k] );
	}

	/// DAMPING ROUTINE	///	///	///	///	///	///	///	///	///	///	///	///	///	
        sp_sys->scf_cnt = 0;

        for(int k=0;k<max_scf+1200;k++)
        {   
	    
	    // Refresh Convergence Checker
            convergence_checker = 0.;
            // Get old eigen_values
            for(int n=0;n<sp_sys->number_of_sp_ion;n++)
                prev_eval[n] = gsl_vector_get(sp_sys->sp_ion[n].eigen_value, sp_cluster_support_get_lowest_state(sp_sys->sp_ion[n].eigen_value));
            // Estimate Eval/Evec
     
            sp_cluster_system_get_h_matrix_mpi(sp_sys,rank,numtasks);
            sp_cluster_system_get_eigensystem(sp_sys);

            // If the number of sp ion is one !!!
            if( sp_sys->number_of_sp_ion == 1 )
            {   ret = SP_SYSTEM_TRUE;
                sp_sys->is_scf_done = SP_SYSTEM_FALSE;   // Turn Flag On; SCF IS DONE               
                sp_cluster_system_sp_energy_update(sp_sys);

                MPI_Barrier(MPI_COMM_WORLD);
                return ret;                     
            }
            // Single sp ion case

            sp_sys->scf_cnt++;

            //if(rank==0) 
                //printf("DAMPING SCF CYCLE %d:\t\t",k+1);
            // SCF Convergence Test by evaluating NewEval - PrevEval
            for(int n=0;n<sp_sys->number_of_sp_ion;n++)
            {    
                convergence_checker += pow(( gsl_vector_get(sp_sys->sp_ion[n].eigen_value,sp_cluster_support_get_lowest_state(sp_sys->sp_ion[n].eigen_value)) - prev_eval[n] ), 2.);
                
                //TEMPOLAR PT.
                //if(rank==0) printf("%.8lf\t", gsl_vector_get(sp_sys->sp_ion[n].eigen_value, sp_cluster_support_get_lowest_state(sp_sys->sp_ion[n].eigen_value) ) );

            }
            //if(rank==0) printf("\n");

            convergence_checker = sqrt(convergence_checker)/sp_sys->number_of_sp_ion;

            if( convergence_checker < SP_SYSTEM_SCF_TOL )
            {   ret = SP_SYSTEM_TRUE;
                sp_sys->is_scf_done = SP_SYSTEM_TRUE;
                // ONCE SCF IS ACHIEVED ENERGY UPDATE NEEDS TO BE CARRIED OUT !
                sp_cluster_system_sp_energy_update(sp_sys);

                MPI_Barrier(MPI_COMM_WORLD);

                return ret;
            }
            else
	    {
		for(int n=0;n<sp_sys->number_of_sp_ion;n++)
		{
			damping_low_state_index = sp_cluster_support_get_lowest_state( sp_sys->sp_ion[n].eigen_value );
			
			// check if s_k-1 and s_k has same sign
			if( gsl_matrix_get(sp_sys->sp_ion[n].eigen_vector,0,damping_low_state_index)*damping_prev_evec[n][0] < 0. )	// if sign different -> result (-)
			{	for(int ii=0;ii<4;ii++)
					damping_prev_evec[n][ii] = damping_prev_evec[n][ii]*-1.;
			}
			// else if sign is same
			// sign check end !!!
			damping_norm = 0.;
			// mixing
			for(int ii=0;ii<4;ii++)
			{
				damping_prev_evec[n][ii] = (0.9*damping_prev_evec[n][ii]+0.1*gsl_matrix_get(sp_sys->sp_ion[n].eigen_vector,ii,damping_low_state_index));
				damping_norm += damping_prev_evec[n][ii]*damping_prev_evec[n][ii];
			}
			damping_norm = sqrt(damping_norm);

			for(int ii=0;ii<4;ii++)
				gsl_matrix_set(sp_sys->sp_ion[n].eigen_vector,ii,damping_low_state_index,damping_prev_evec[n][ii]/damping_norm);
			// mixing end !!!
		}

		// prev info record
		for(int n=0;n<sp_sys->number_of_sp_ion;n++)
		{	damping_low_state_index = sp_cluster_support_get_lowest_state(sp_sys->sp_ion[n].eigen_value);	// get low index
			damping_prev_evec[n][0] = gsl_matrix_get(sp_sys->sp_ion[n].eigen_vector,0,damping_low_state_index);	// get s evec
			damping_prev_evec[n][1] = gsl_matrix_get(sp_sys->sp_ion[n].eigen_vector,1,damping_low_state_index);	// get s evec
			damping_prev_evec[n][2] = gsl_matrix_get(sp_sys->sp_ion[n].eigen_vector,2,damping_low_state_index);	// get s evec
			damping_prev_evec[n][3] = gsl_matrix_get(sp_sys->sp_ion[n].eigen_vector,3,damping_low_state_index);	// get s evec
		}
		// prev info record done
		
	    }
        }
	/// DAMPING ROUTINE END ///	///	///	///	///	///	///	///	///	///	///	///	///

    }
    else if( sp_sys->is_scf_done == SP_SYSTEM_TRUE )
    {
        // update only h_matriced without diagonalisation, i.e., keep the eigenvectors
        sp_cluster_system_get_h_matrix_mpi(sp_sys,rank,numtasks);
        // Based on the recalculated h_matrices .. Calculate sp-energies
        sp_cluster_system_sp_energy_update(sp_sys);

    }
    MPI_Barrier(MPI_COMM_WORLD);
    return ret;
}





// This is a function for updating sp-derived energy
// must be called after scf convergence achieved !
double sp_cluster_system_get_sp_energy( sp_cluster_type_sp_ion* sp, const int k /*which state?*/)
{     
    double Return = 0.;
    double evec[4];
    evec[0] = gsl_matrix_get(sp->eigen_vector,0,k);
    evec[1] = gsl_matrix_get(sp->eigen_vector,1,k);
    evec[2] = gsl_matrix_get(sp->eigen_vector,2,k);
    evec[3] = gsl_matrix_get(sp->eigen_vector,3,k);     // get evec
    // Calculate Energy
    for(int k=0;k<4;k++)
    {   for(int l=0;l<4;l++)
            Return += evec[k]*evec[l]*gsl_matrix_get(sp->h_matrix,k,l);
    }
    return Return;
}


// CALL D_EVEC CALCULATOR; case when i=j
// 
//      input: sp_cluster_system* sp_sys, double*** dh_matrix_tmp, int i, int j, const int type (=0 vs_sp//=1 vs_cal) 
//
// dHj;uv/di -> derive_evec_sp[i][j][..][..]; 
// 2019/09/25
void test__test_support( sp_cluster_system* sp_sys, double*** dh_matrix, const int i /*diff with*/, const int j/*diff tar*/, const int type, double**** out )
{
    double e[4];
    double de[3] = {0,0,0};
    double tmp[3] = {0,0,0};
    const int low_idx = sp_cluster_support_get_lowest_state(sp_sys->sp_ion[j].eigen_value);

    // get lonepair eigenvalues ... 0th 1st 2nd 3rd
    e[0] = sp_cluster_system_get_sp_energy( &sp_sys->sp_ion[j], 0 );
    e[1] = sp_cluster_system_get_sp_energy( &sp_sys->sp_ion[j], 1 );
    e[2] = sp_cluster_system_get_sp_energy( &sp_sys->sp_ion[j], 2 );
    e[3] = sp_cluster_system_get_sp_energy( &sp_sys->sp_ion[j], 3 );

    // get derivatives of eigenvalue of the lowest eigenvalue ... at this point the dh_matrix is calculated based on the 'n'th iterations dh_matrix !!!
    for(int u=0;u<4;u++)
    {   for(int v=0;v<4;v++)
        {   de[0] += gsl_matrix_get(sp_sys->sp_ion[j].eigen_vector,u,low_idx)*
                        gsl_matrix_get(sp_sys->sp_ion[j].eigen_vector,v,low_idx)*dh_matrix[u][v][0];
            de[1] += gsl_matrix_get(sp_sys->sp_ion[j].eigen_vector,u,low_idx)*
                        gsl_matrix_get(sp_sys->sp_ion[j].eigen_vector,v,low_idx)*dh_matrix[u][v][1];
            de[2] += gsl_matrix_get(sp_sys->sp_ion[j].eigen_vector,u,low_idx)*
                        gsl_matrix_get(sp_sys->sp_ion[j].eigen_vector,v,low_idx)*dh_matrix[u][v][2];    }}

    if( type == 0 ) // vs_sp
    {
        for(int c=0;c<4;c++)     // loop through c (... c1 ,c2 ,c3 and c4 )
        {
            for(int l=0;l<4;l++)
            {
                if( l != low_idx ) // l => j, low_idx => k
                {
                    for(int u=0;u<4;u++)
                    {   for(int v=0;v<4;v++)
                        {
                            tmp[0] += ( gsl_matrix_get(sp_sys->sp_ion[j].eigen_vector,c,l)*gsl_matrix_get(sp_sys->sp_ion[j].eigen_vector,u,l)/(e[l]-e[low_idx])
                                     *( gsl_matrix_get(sp_sys->sp_ion[j].eigen_vector,v,low_idx)*de[0]*sp_cluster_support_kronecker_delta(u,v) 
                                      - gsl_matrix_get(sp_sys->sp_ion[j].eigen_vector,v,low_idx)*dh_matrix[u][v][0] ));

                            tmp[1] += ( gsl_matrix_get(sp_sys->sp_ion[j].eigen_vector,c,l)*gsl_matrix_get(sp_sys->sp_ion[j].eigen_vector,u,l)/(e[l]-e[low_idx])
                                     *( gsl_matrix_get(sp_sys->sp_ion[j].eigen_vector,v,low_idx)*de[1]*sp_cluster_support_kronecker_delta(u,v) 
                                      - gsl_matrix_get(sp_sys->sp_ion[j].eigen_vector,v,low_idx)*dh_matrix[u][v][1] ));

                            tmp[2] += ( gsl_matrix_get(sp_sys->sp_ion[j].eigen_vector,c,l)*gsl_matrix_get(sp_sys->sp_ion[j].eigen_vector,u,l)/(e[l]-e[low_idx])
                                     *( gsl_matrix_get(sp_sys->sp_ion[j].eigen_vector,v,low_idx)*de[2]*sp_cluster_support_kronecker_delta(u,v) 
                                      - gsl_matrix_get(sp_sys->sp_ion[j].eigen_vector,v,low_idx)*dh_matrix[u][v][2] ));
                        }//v
                    }//u
                }//if( l != low_idx(k) )
            }//l
            out[i][j][c][0] = tmp[0]; out[i][j][c][1] = tmp[1]; out[i][j][c][2] = tmp[2];
            memset(tmp,0,3*sizeof(double));
        }//c
    }//if type
    else if( type == 1 ) // vs_cla
    {
        for(int c=0;c<4;c++)     // loop through c (... c1 ,c2 ,c3 and c4 )
        {
            for(int l=0;l<4;l++)
            {
                if( l != low_idx ) // l => j, low_idx => k
                {
                    for(int u=0;u<4;u++)
                    {   for(int v=0;v<4;v++)
                        {
                            tmp[0] += ( gsl_matrix_get(sp_sys->sp_ion[j].eigen_vector,c,l)*gsl_matrix_get(sp_sys->sp_ion[j].eigen_vector,u,l)/(e[l]-e[low_idx])
                                     *( gsl_matrix_get(sp_sys->sp_ion[j].eigen_vector,v,low_idx)*de[0]*sp_cluster_support_kronecker_delta(u,v) 
                                      - gsl_matrix_get(sp_sys->sp_ion[j].eigen_vector,v,low_idx)*dh_matrix[u][v][0] ));

                            tmp[1] += ( gsl_matrix_get(sp_sys->sp_ion[j].eigen_vector,c,l)*gsl_matrix_get(sp_sys->sp_ion[j].eigen_vector,u,l)/(e[l]-e[low_idx])
                                     *( gsl_matrix_get(sp_sys->sp_ion[j].eigen_vector,v,low_idx)*de[1]*sp_cluster_support_kronecker_delta(u,v) 
                                      - gsl_matrix_get(sp_sys->sp_ion[j].eigen_vector,v,low_idx)*dh_matrix[u][v][1] ));

                            tmp[2] += ( gsl_matrix_get(sp_sys->sp_ion[j].eigen_vector,c,l)*gsl_matrix_get(sp_sys->sp_ion[j].eigen_vector,u,l)/(e[l]-e[low_idx])
                                     *( gsl_matrix_get(sp_sys->sp_ion[j].eigen_vector,v,low_idx)*de[2]*sp_cluster_support_kronecker_delta(u,v) 
                                      - gsl_matrix_get(sp_sys->sp_ion[j].eigen_vector,v,low_idx)*dh_matrix[u][v][2] ));
                        }//v
                    }//u
                }//if( l != low_idx(k) )
            }//l
            out[i][j][c][0] = tmp[0]; out[i][j][c][1] = tmp[1]; out[i][j][c][2] = tmp[2];
            memset(tmp,0,3*sizeof(double));
        }//c
    }//if type

    return;
}

//2019/09/23 test
void test__test( sp_cluster_system* sp_sys, const int rank, const int numtasks )
{   
    const int cycmx = 25;
    const int number_of_sp_ion = sp_sys->number_of_sp_ion;
    const int number_of_classic_ion = sp_sys->number_of_classic_ion;
    double ssqr = 0.;
    // This Routine estimates the derivatives of eigenvectors in SLAM system by an iterative method
    // double dh_matrix_tmp[4][4][3];
    double*** dh_matrix_tmp = (double***)malloc(4*sizeof(double**));
    for(int i=0;i<4;i++)
        dh_matrix_tmp[i] = (double**)malloc(4*sizeof(double*));
    for(int i=0;i<4;i++)
    {   for(int j=0;j<4;j++)
        {   dh_matrix_tmp[i][j] = (double*)calloc(3,sizeof(double));    }}

    double**** deriv_evec_sp_tmp = (double****)malloc(number_of_sp_ion*sizeof(double***));
    for(int i=0;i<number_of_sp_ion;i++)
        deriv_evec_sp_tmp[i] = (double***)malloc(number_of_sp_ion*sizeof(double**));
    for(int i=0;i<number_of_sp_ion;i++)
    {   for(int j=0;j<number_of_sp_ion;j++)
        {   deriv_evec_sp_tmp[i][j] = (double**)malloc(4*sizeof(double*));      }}
    for(int i=0;i<number_of_sp_ion;i++)
    {   for(int j=0;j<number_of_sp_ion;j++)
        {   for(int k=0;k<4;k++)
            {   deriv_evec_sp_tmp[i][j][k] = (double*)calloc(3,sizeof(double)); }}}

    double**** deriv_evec_cla_tmp = (double****)malloc(number_of_classic_ion*sizeof(double***));
    for(int i=0;i<number_of_classic_ion;i++)
        deriv_evec_cla_tmp[i] = (double***)malloc(number_of_sp_ion*sizeof(double**));
    for(int i=0;i<number_of_classic_ion;i++)
    {   for(int j=0;j<number_of_sp_ion;j++)
        {   deriv_evec_cla_tmp[i][j] = (double**)malloc(4*sizeof(double*));      }}
    for(int i=0;i<number_of_classic_ion;i++)
    {   for(int j=0;j<number_of_sp_ion;j++)
        {   for(int k=0;k<4;k++)
            {   deriv_evec_cla_tmp[i][j][k] = (double*)calloc(3,sizeof(double)); }}}


    double integral_x, integral_y, integral_z;
    double evk[4];

    // memset evec derivatives
    for(int n=0;n<number_of_sp_ion;n++)
    {   for(int m=0;m<number_of_sp_ion;m++)
        {   for(int o=0;o<4;o++)
            {   memset(sp_sys->deriv_evec_sp[n][m][o],0.,3*sizeof(double));     }}}
    for(int n=0;n<number_of_classic_ion;n++)
    {   for(int m=0;m<number_of_sp_ion;m++)
        {   for(int o=0;o<4;o++)
            {   memset(sp_sys->deriv_evec_cla[n][m][o],0.,3*sizeof(double));    }}}
    // memset dh_matrix_tmp
    for(int i=0;i<4;i++)
    {   for(int j=0;j<4;j++)
        {   memset(dh_matrix_tmp[i][j],0.,3*sizeof(double));        }}

    /* PHASE 1 - Calculate derivatives of H matrix w.r.t. sp ion moves
     *
     *      calculate evec derivatives
     *
     * PHASE 2 - Calculate derivatives of H matrix w.r.t. claion moves
     *
     *      calculate evec derivatives
     */


for(int cycle=0;cycle<cycmx;cycle++)
{
    // PHASE 1 : vs sp_ion
    for(int i=0;i<number_of_sp_ion;i++)
    {
        // loop through sp ions   
        for(int j=0;j<number_of_sp_ion;j++)
        {
            if( i == j ) // 
            {
                for(int k=0;k<number_of_classic_ion;k++)
                {
                    for(int u=0;u<4;u++)
                    {   for(int v=0;v<4;v++)
                        {
                            dh_matrix_tmp[u][v][0] += gsl_matrix_get(sp_sys->cla_dh_matrix_x[j][k],u,v);
                            dh_matrix_tmp[u][v][1] += gsl_matrix_get(sp_sys->cla_dh_matrix_y[j][k],u,v);
                            dh_matrix_tmp[u][v][2] += gsl_matrix_get(sp_sys->cla_dh_matrix_z[j][k],u,v);    
                        }
                    }
                }//vs classic ions ... for 'j'th sp-ion with respect to 'k'th classic ion

                for(int k=0;k<number_of_sp_ion;k++)
                {
                    if( j != k )
                    {   
                        integral_x = sp_cluster_integrator_get_x_12(&sp_sys->sp_ion[k]);    
                        integral_y = integral_x;    integral_z = integral_x;
                        //integral_y = sp_cluster_integrator_get_y_13(&sp_sys->sp_ion[k]);    
                        //integral_z = sp_cluster_integrator_get_z_14(&sp_sys->sp_ion[k]);        
                        evk[0] = gsl_matrix_get(sp_sys->sp_ion[k].eigen_vector,0,sp_cluster_support_get_lowest_state(sp_sys->sp_ion[k].eigen_value));
                        evk[1] = gsl_matrix_get(sp_sys->sp_ion[k].eigen_vector,1,sp_cluster_support_get_lowest_state(sp_sys->sp_ion[k].eigen_value));
                        evk[2] = gsl_matrix_get(sp_sys->sp_ion[k].eigen_vector,2,sp_cluster_support_get_lowest_state(sp_sys->sp_ion[k].eigen_value));
                        evk[3] = gsl_matrix_get(sp_sys->sp_ion[k].eigen_vector,3,sp_cluster_support_get_lowest_state(sp_sys->sp_ion[k].eigen_value));

                        for(int u=0;u<4;u++)
                        {   for(int v=0;v<4;v++)
                            {
                                dh_matrix_tmp[u][v][0] += gsl_matrix_get(sp_sys->dh_matrix_x[j][k],u,v);
                                dh_matrix_tmp[u][v][1] += gsl_matrix_get(sp_sys->dh_matrix_y[j][k],u,v);
                                dh_matrix_tmp[u][v][2] += gsl_matrix_get(sp_sys->dh_matrix_z[j][k],u,v);    
                                //vs sp ion monopole

                                dh_matrix_tmp[u][v][0] += ( -2.*evk[0]*evk[1]*integral_x*gsl_matrix_get(sp_sys->ddh_matrix_xx[j][k],u,v)
                                                            -2.*evk[0]*evk[2]*integral_y*gsl_matrix_get(sp_sys->ddh_matrix_yx[j][k],u,v)
                                                            -2.*evk[0]*evk[3]*integral_z*gsl_matrix_get(sp_sys->ddh_matrix_zx[j][k],u,v) ); 

                                dh_matrix_tmp[u][v][1] += ( -2.*evk[0]*evk[1]*integral_x*gsl_matrix_get(sp_sys->ddh_matrix_xy[j][k],u,v)
                                                            -2.*evk[0]*evk[2]*integral_y*gsl_matrix_get(sp_sys->ddh_matrix_yy[j][k],u,v)
                                                            -2.*evk[0]*evk[3]*integral_z*gsl_matrix_get(sp_sys->ddh_matrix_zy[j][k],u,v) ); 

                                dh_matrix_tmp[u][v][2] += ( -2.*evk[0]*evk[1]*integral_x*gsl_matrix_get(sp_sys->ddh_matrix_xz[j][k],u,v)
                                                            -2.*evk[0]*evk[2]*integral_y*gsl_matrix_get(sp_sys->ddh_matrix_yz[j][k],u,v)
                                                            -2.*evk[0]*evk[3]*integral_z*gsl_matrix_get(sp_sys->ddh_matrix_zz[j][k],u,v) ); 
                                //vs sp ion dipole      
                                
                                    // contribution by evec derivatives
                                dh_matrix_tmp[u][v][0] += ( -2.*integral_x*gsl_matrix_get(sp_sys->scf_dh_x_matrix_vs_sp_ion_dipole[j][k],u,v)*(sp_sys->deriv_evec_sp[i][k][0][0]*evk[1]+evk[0]*sp_sys->deriv_evec_sp[i][k][1][0])
                                                            -2.*integral_y*gsl_matrix_get(sp_sys->scf_dh_y_matrix_vs_sp_ion_dipole[j][k],u,v)*(sp_sys->deriv_evec_sp[i][k][0][0]*evk[2]+evk[0]*sp_sys->deriv_evec_sp[i][k][2][0])
                                                            -2.*integral_z*gsl_matrix_get(sp_sys->scf_dh_z_matrix_vs_sp_ion_dipole[j][k],u,v)*(sp_sys->deriv_evec_sp[i][k][0][0]*evk[3]+evk[0]*sp_sys->deriv_evec_sp[i][k][3][0]) );

                                dh_matrix_tmp[u][v][1] += ( -2.*integral_x*gsl_matrix_get(sp_sys->scf_dh_x_matrix_vs_sp_ion_dipole[j][k],u,v)*(sp_sys->deriv_evec_sp[i][k][0][1]*evk[1]+evk[0]*sp_sys->deriv_evec_sp[i][k][1][1])
                                                            -2.*integral_y*gsl_matrix_get(sp_sys->scf_dh_y_matrix_vs_sp_ion_dipole[j][k],u,v)*(sp_sys->deriv_evec_sp[i][k][0][1]*evk[2]+evk[0]*sp_sys->deriv_evec_sp[i][k][2][1])
                                                            -2.*integral_z*gsl_matrix_get(sp_sys->scf_dh_z_matrix_vs_sp_ion_dipole[j][k],u,v)*(sp_sys->deriv_evec_sp[i][k][0][1]*evk[3]+evk[0]*sp_sys->deriv_evec_sp[i][k][3][1]) );

                                dh_matrix_tmp[u][v][2] += ( -2.*integral_x*gsl_matrix_get(sp_sys->scf_dh_x_matrix_vs_sp_ion_dipole[j][k],u,v)*(sp_sys->deriv_evec_sp[i][k][0][2]*evk[1]+evk[0]*sp_sys->deriv_evec_sp[i][k][1][2])
                                                            -2.*integral_y*gsl_matrix_get(sp_sys->scf_dh_y_matrix_vs_sp_ion_dipole[j][k],u,v)*(sp_sys->deriv_evec_sp[i][k][0][2]*evk[2]+evk[0]*sp_sys->deriv_evec_sp[i][k][2][2])
                                                            -2.*integral_z*gsl_matrix_get(sp_sys->scf_dh_z_matrix_vs_sp_ion_dipole[j][k],u,v)*(sp_sys->deriv_evec_sp[i][k][0][2]*evk[3]+evk[0]*sp_sys->deriv_evec_sp[i][k][3][2]) );
                            }//v
                        }//u

                    }// end if
                }//vs sp ions

                // CALL D_EVEC CALCULATOR; case when i=j
                // 
                //      input: sp_cluster_system* sp_sys, double*** dh_matrix_tmp, int j, int i, const int type (=0 vs_sp//=1 vs_cal) 
                //
                // dHj;uv/di -> derive_evec_sp[i][j][..][..]; 
                //test__test_support( sp_cluster_system* sp_sys, double*** dh_matrix, const int i /*diff with*/, const int j/*diff tar*/, const int type )
                test__test_support( sp_sys, dh_matrix_tmp, i /*diff with*/, j/*diff tar*/, 0 , deriv_evec_sp_tmp );

                // memset dh_matrix_tmp
                for(int u=0;u<4;u++)
                {   for(int v=0;v<4;v++)
                    {   memset(dh_matrix_tmp[u][v],0.,3*sizeof(double));    }}
            }// end if( i == j )
            else // if( i != j )// when derivative w.r.t other sp-ion centres -> MM integral terms vanish
            {
                for(int k=0;k<number_of_sp_ion;k++)
                {   
                    if( j != k )
                    {
                        integral_x = sp_cluster_integrator_get_x_12(&sp_sys->sp_ion[k]);    
                        integral_y = integral_x;    integral_z = integral_x;
                        //integral_y = sp_cluster_integrator_get_y_13(&sp_sys->sp_ion[j]);    
                        //integral_z = sp_cluster_integrator_get_z_14(&sp_sys->sp_ion[j]);        
                        evk[0] = gsl_matrix_get(sp_sys->sp_ion[k].eigen_vector,0,sp_cluster_support_get_lowest_state(sp_sys->sp_ion[k].eigen_value));
                        evk[1] = gsl_matrix_get(sp_sys->sp_ion[k].eigen_vector,1,sp_cluster_support_get_lowest_state(sp_sys->sp_ion[k].eigen_value));
                        evk[2] = gsl_matrix_get(sp_sys->sp_ion[k].eigen_vector,2,sp_cluster_support_get_lowest_state(sp_sys->sp_ion[k].eigen_value));
                        evk[3] = gsl_matrix_get(sp_sys->sp_ion[k].eigen_vector,3,sp_cluster_support_get_lowest_state(sp_sys->sp_ion[k].eigen_value));

                        for(int u=0;u<4;u++)
                        {   for(int v=0;v<4;v++)
                            {
                                if( i == k )// when 'x' == 'beta'
                                {
                                    dh_matrix_tmp[u][v][0] += -gsl_matrix_get(sp_sys->dh_matrix_x[j][k],u,v);
                                    dh_matrix_tmp[u][v][1] += -gsl_matrix_get(sp_sys->dh_matrix_y[j][k],u,v);
                                    dh_matrix_tmp[u][v][2] += -gsl_matrix_get(sp_sys->dh_matrix_z[j][k],u,v);    
                                    //vs sp ion monopole

                                    dh_matrix_tmp[u][v][0] += -( -2.*evk[0]*evk[1]*integral_x*gsl_matrix_get(sp_sys->ddh_matrix_xx[j][k],u,v)
                                                                 -2.*evk[0]*evk[2]*integral_y*gsl_matrix_get(sp_sys->ddh_matrix_yx[j][k],u,v)
                                                                 -2.*evk[0]*evk[3]*integral_z*gsl_matrix_get(sp_sys->ddh_matrix_zx[j][k],u,v) ); 

                                    dh_matrix_tmp[u][v][1] += -( -2.*evk[0]*evk[1]*integral_x*gsl_matrix_get(sp_sys->ddh_matrix_xy[j][k],u,v)
                                                                 -2.*evk[0]*evk[2]*integral_y*gsl_matrix_get(sp_sys->ddh_matrix_yy[j][k],u,v)
                                                                 -2.*evk[0]*evk[3]*integral_z*gsl_matrix_get(sp_sys->ddh_matrix_zy[j][k],u,v) ); 

                                    dh_matrix_tmp[u][v][2] += -( -2.*evk[0]*evk[1]*integral_x*gsl_matrix_get(sp_sys->ddh_matrix_xz[j][k],u,v)
                                                                 -2.*evk[0]*evk[2]*integral_y*gsl_matrix_get(sp_sys->ddh_matrix_yz[j][k],u,v)
                                                                 -2.*evk[0]*evk[3]*integral_z*gsl_matrix_get(sp_sys->ddh_matrix_zz[j][k],u,v) ); 
                                    //vs sp ion dipole      
                                }// end if( i == k )


                                    // contribution by evec derivatives
                                dh_matrix_tmp[u][v][0] += ( -2.*integral_x*gsl_matrix_get(sp_sys->scf_dh_x_matrix_vs_sp_ion_dipole[j][k],u,v)*(sp_sys->deriv_evec_sp[i][k][0][0]*evk[1]+evk[0]*sp_sys->deriv_evec_sp[i][k][1][0])
                                                            -2.*integral_y*gsl_matrix_get(sp_sys->scf_dh_y_matrix_vs_sp_ion_dipole[j][k],u,v)*(sp_sys->deriv_evec_sp[i][k][0][0]*evk[2]+evk[0]*sp_sys->deriv_evec_sp[i][k][2][0])
                                                            -2.*integral_z*gsl_matrix_get(sp_sys->scf_dh_z_matrix_vs_sp_ion_dipole[j][k],u,v)*(sp_sys->deriv_evec_sp[i][k][0][0]*evk[3]+evk[0]*sp_sys->deriv_evec_sp[i][k][3][0]) );

                                dh_matrix_tmp[u][v][1] += ( -2.*integral_x*gsl_matrix_get(sp_sys->scf_dh_x_matrix_vs_sp_ion_dipole[j][k],u,v)*(sp_sys->deriv_evec_sp[i][k][0][1]*evk[1]+evk[0]*sp_sys->deriv_evec_sp[i][k][1][1])
                                                            -2.*integral_y*gsl_matrix_get(sp_sys->scf_dh_y_matrix_vs_sp_ion_dipole[j][k],u,v)*(sp_sys->deriv_evec_sp[i][k][0][1]*evk[2]+evk[0]*sp_sys->deriv_evec_sp[i][k][2][1])
                                                            -2.*integral_z*gsl_matrix_get(sp_sys->scf_dh_z_matrix_vs_sp_ion_dipole[j][k],u,v)*(sp_sys->deriv_evec_sp[i][k][0][1]*evk[3]+evk[0]*sp_sys->deriv_evec_sp[i][k][3][1]) );

                                dh_matrix_tmp[u][v][2] += ( -2.*integral_x*gsl_matrix_get(sp_sys->scf_dh_x_matrix_vs_sp_ion_dipole[j][k],u,v)*(sp_sys->deriv_evec_sp[i][k][0][2]*evk[1]+evk[0]*sp_sys->deriv_evec_sp[i][k][1][2])
                                                            -2.*integral_y*gsl_matrix_get(sp_sys->scf_dh_y_matrix_vs_sp_ion_dipole[j][k],u,v)*(sp_sys->deriv_evec_sp[i][k][0][2]*evk[2]+evk[0]*sp_sys->deriv_evec_sp[i][k][2][2])
                                                            -2.*integral_z*gsl_matrix_get(sp_sys->scf_dh_z_matrix_vs_sp_ion_dipole[j][k],u,v)*(sp_sys->deriv_evec_sp[i][k][0][2]*evk[3]+evk[0]*sp_sys->deriv_evec_sp[i][k][3][2]) );
                            }//v
                        }//u

                    }// if( k !=j )
                }//end k

                // CALL D_EVEC CALCULATOR; case when i!=j
                // dHj;uv/di -> derive_evec_sp[i][j][..][..]; 
                //test__test_support( sp_cluster_system* sp_sys, double*** dh_matrix, const int i /*diff with*/, const int j/*diff tar*/, const int type )
                test__test_support( sp_sys, dh_matrix_tmp, i /*diff with*/, j/*diff tar*/, 0, deriv_evec_sp_tmp );
                
                // memset dh_matrix_tmp
                for(int u=0;u<4;u++)
                {   for(int v=0;v<4;v++)
                    {   memset(dh_matrix_tmp[u][v],0.,3*sizeof(double));    }}
            }// end of if 

        }// end of sp ion j(alpha)
    }// i loop('x')

    // PHASE 2: vs classic ion
    for(int i=0;i<number_of_classic_ion;i++)
    {
        // loop through sp ions   
        for(int j=0;j<number_of_sp_ion;j++)
        {
            for(int k=0;k<number_of_classic_ion;k++)
            {
                if( i == k )
                {
                    for(int u=0;u<4;u++)
                    {   for(int v=0;v<4;v++)
                        {   dh_matrix_tmp[u][v][0] += -gsl_matrix_get(sp_sys->cla_dh_matrix_x[j][k],u,v);
                            dh_matrix_tmp[u][v][1] += -gsl_matrix_get(sp_sys->cla_dh_matrix_y[j][k],u,v);
                            dh_matrix_tmp[u][v][2] += -gsl_matrix_get(sp_sys->cla_dh_matrix_z[j][k],u,v);    }} 
                    //////// THIS BIT ERROR -> [i] goes over the limit index !!!!!!!!!!!!!!!!!!!!!!!! 2019/09/25
                }
            }//vs classic ions
            

            
            // EVEC DERIVATIVE CONTRIBUTION w.r.t classic ions
            for(int k=0;k<number_of_sp_ion;k++)
            {
                if( j != k )
                {   
                    integral_x = sp_cluster_integrator_get_x_12(&sp_sys->sp_ion[k]);    
                    integral_y = integral_x;    integral_z = integral_x;
                    //integral_y = sp_cluster_integrator_get_y_13(&sp_sys->sp_ion[k]);    
                    //integral_z = sp_cluster_integrator_get_z_14(&sp_sys->sp_ion[k]);        
                    evk[0] = gsl_matrix_get(sp_sys->sp_ion[k].eigen_vector,0,sp_cluster_support_get_lowest_state(sp_sys->sp_ion[k].eigen_value));
                    evk[1] = gsl_matrix_get(sp_sys->sp_ion[k].eigen_vector,1,sp_cluster_support_get_lowest_state(sp_sys->sp_ion[k].eigen_value));
                    evk[2] = gsl_matrix_get(sp_sys->sp_ion[k].eigen_vector,2,sp_cluster_support_get_lowest_state(sp_sys->sp_ion[k].eigen_value));
                    evk[3] = gsl_matrix_get(sp_sys->sp_ion[k].eigen_vector,3,sp_cluster_support_get_lowest_state(sp_sys->sp_ion[k].eigen_value));

                    for(int u=0;u<4;u++)
                    {   for(int v=0;v<4;v++)
                        {
                                // contribution by evec derivatives
                            dh_matrix_tmp[u][v][0] += ( -2.*integral_x*gsl_matrix_get(sp_sys->scf_dh_x_matrix_vs_sp_ion_dipole[j][k],u,v)*(sp_sys->deriv_evec_cla[i][k][0][0]*evk[1]+evk[0]*sp_sys->deriv_evec_cla[i][k][1][0])
                                                        -2.*integral_y*gsl_matrix_get(sp_sys->scf_dh_y_matrix_vs_sp_ion_dipole[j][k],u,v)*(sp_sys->deriv_evec_cla[i][k][0][0]*evk[2]+evk[0]*sp_sys->deriv_evec_cla[i][k][2][0])
                                                        -2.*integral_z*gsl_matrix_get(sp_sys->scf_dh_z_matrix_vs_sp_ion_dipole[j][k],u,v)*(sp_sys->deriv_evec_cla[i][k][0][0]*evk[3]+evk[0]*sp_sys->deriv_evec_cla[i][k][3][0]) );

                            dh_matrix_tmp[u][v][1] += ( -2.*integral_x*gsl_matrix_get(sp_sys->scf_dh_x_matrix_vs_sp_ion_dipole[j][k],u,v)*(sp_sys->deriv_evec_cla[i][k][0][1]*evk[1]+evk[0]*sp_sys->deriv_evec_cla[i][k][1][1])
                                                        -2.*integral_y*gsl_matrix_get(sp_sys->scf_dh_y_matrix_vs_sp_ion_dipole[j][k],u,v)*(sp_sys->deriv_evec_cla[i][k][0][1]*evk[2]+evk[0]*sp_sys->deriv_evec_cla[i][k][2][1])
                                                        -2.*integral_z*gsl_matrix_get(sp_sys->scf_dh_z_matrix_vs_sp_ion_dipole[j][k],u,v)*(sp_sys->deriv_evec_cla[i][k][0][1]*evk[3]+evk[0]*sp_sys->deriv_evec_cla[i][k][3][1]) );

                            dh_matrix_tmp[u][v][2] += ( -2.*integral_x*gsl_matrix_get(sp_sys->scf_dh_x_matrix_vs_sp_ion_dipole[j][k],u,v)*(sp_sys->deriv_evec_cla[i][k][0][2]*evk[1]+evk[0]*sp_sys->deriv_evec_cla[i][k][1][2])
                                                        -2.*integral_y*gsl_matrix_get(sp_sys->scf_dh_y_matrix_vs_sp_ion_dipole[j][k],u,v)*(sp_sys->deriv_evec_cla[i][k][0][2]*evk[2]+evk[0]*sp_sys->deriv_evec_cla[i][k][2][2])
                                                        -2.*integral_z*gsl_matrix_get(sp_sys->scf_dh_z_matrix_vs_sp_ion_dipole[j][k],u,v)*(sp_sys->deriv_evec_cla[i][k][0][2]*evk[3]+evk[0]*sp_sys->deriv_evec_cla[i][k][3][2]) );
                        }//v
                    }//u

                }// end if
            }//vs sp ions

            // CALL D_EVEC CALCULATOR; case when i=j
            // dHj;uv/di -> derive_evec_cla[i][j][..][..]; 
            //test__test_support( sp_cluster_system* sp_sys, double*** dh_matrix, const int i /*diff with*/, const int j/*diff tar*/, const int type )
            test__test_support( sp_sys, dh_matrix_tmp, i /*diff with*/, j/*diff tar*/, 1, deriv_evec_cla_tmp );

            // memset dh_matrix_tmp
            for(int u=0;u<4;u++)
            {   for(int v=0;v<4;v++)
                    {   memset(dh_matrix_tmp[u][v],0.,3*sizeof(double));    }}
        }// end of sp ion j(alpha)

    }// i loop('x')

    // update evec deriv
    //
    //
    // Indices how to read ?
    //
    // [i][j] (first two) for 'j' w.r.t 'i' 
    //
    // e.g., deriv_evec_sp[3][1][..][..] 
    //
    //  for '1' sp ion w.r.t '3' sp ion moves !!!
    //
    // [k][l] (last  two) for 'k'th eigenvector coefficient w.r.t 'l'th -> (x,y,z)
    //
    // e.g., [2][2]
    //
    // for '2' eigenvector element diff w.r.t y
    //
    for(int i=0;i<number_of_sp_ion;i++)
    {   for(int j=0;j<number_of_sp_ion;j++)
        {   for(int k=0;k<4;k++)
            {   for(int l=0;l<3;l++)
                {   ssqr += pow(sp_sys->deriv_evec_sp[i][j][k][l] - deriv_evec_sp_tmp[i][j][k][l],2.);
                    sp_sys->deriv_evec_sp[i][j][k][l] = deriv_evec_sp_tmp[i][j][k][l];  
                    //printf("%lf\n",sp_sys->deriv_evec_sp[i][j][k][l]);
                }}}}
    for(int i=0;i<number_of_classic_ion;i++)
    {   for(int j=0;j<number_of_sp_ion;j++)
        {   for(int k=0;k<4;k++)
            {   for(int l=0;l<3;l++)
                {   ssqr += pow(sp_sys->deriv_evec_cla[i][j][k][l] - deriv_evec_cla_tmp[i][j][k][l],2.);
                    sp_sys->deriv_evec_cla[i][j][k][l] = deriv_evec_cla_tmp[i][j][k][l];  
                }}}}

    ssqr = sqrt(ssqr)/(3*(number_of_sp_ion*number_of_sp_ion*4+number_of_classic_ion*number_of_sp_ion*4));
    
    /* ALGORITHM DESCRIPTION 
     *
     * 
     *
     *
     *
     *
     *
     */

    if( ssqr < SP_SYSTEM_EVEC_TOL )
        break;
    ssqr = 0.;

}//cycle

    // Check Convergence (?)

    ///     ///     ///     ///     ///     ///     ///     ///     ///     FINALISE    ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///
    ///     ///     ///     ///     ///     ///     ///     ///     ///     FINALISE    ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///
    ///     ///     ///     ///     ///     ///     ///     ///     ///     FINALISE    ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///
    ///     ///     ///     ///     ///     ///     ///     ///     ///     FINALISE    ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///
    ///     ///     ///     ///     ///     ///     ///     ///     ///     FINALISE    ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///

    // workspace dealloc
    for(int i=0;i<number_of_sp_ion;i++)
    {   for(int j=0;j<number_of_sp_ion;j++)
        {   for(int k=0;k<4;k++)
            {   free(deriv_evec_sp_tmp[i][j][k]);  }}}
    for(int i=0;i<number_of_sp_ion;i++)
    {   for(int j=0;j<number_of_sp_ion;j++)
        {   free(deriv_evec_sp_tmp[i][j]);  }}
    for(int i=0;i<number_of_sp_ion;i++)
    {   free(deriv_evec_sp_tmp[i]);  }
    free(deriv_evec_sp_tmp);

    for(int i=0;i<number_of_classic_ion;i++)
    {   for(int j=0;j<number_of_sp_ion;j++)
        {   for(int k=0;k<4;k++)
            {   free(deriv_evec_cla_tmp[i][j][k]);  }}}
    for(int i=0;i<number_of_classic_ion;i++)
    {   for(int j=0;j<number_of_sp_ion;j++)
        {   free(deriv_evec_cla_tmp[i][j]);  }}
    for(int i=0;i<number_of_classic_ion;i++)
    {   free(deriv_evec_cla_tmp[i]);  }
    free(deriv_evec_cla_tmp);

    for(int i=0;i<4;i++)
    {   for(int j=0;j<4;j++)
        {   free(dh_matrix_tmp[i][j]);  }}
    for(int i=0;i<4;i++)    
        free(dh_matrix_tmp[i]);
    free(dh_matrix_tmp);

    return;
}









// FORCE CALCULATION MODULE
void sp_cluster_system_get_force_mpi( sp_cluster_system* sp_sys, int rank, int numtasks )
{
    /* FUNCTION OUTLINE...
     *
     *  Phase 1 - Estimation of sp - sp interaction force
     *  
     *          1. sp - sp derivative of monopole + sp core interaction
     *
     *          2. sp - sp derivative of dipolar interaction
     *
     *          e.g., 'n'th sp ion vs 'm'th sp ion;
     *
     *  Phase 2 - Estimation of sp - classic ion interaction force
     *
     *          e.g., 'n'th sp ion vs 'm'th classic ion  
     *
     */

    // MPI Variables
    int idxm,idxn;
    int L,R;                        // round-robin scheme index 1.
    int ista,iend;                  // round-robin scheme index 2.
    int len;                        // total number of for loop cycles
    // Common Variables
    const int number_of_classic_ion = sp_sys->number_of_classic_ion;
    const int number_of_sp_ion      = sp_sys->number_of_sp_ion;
    // vs Classic Ion Variables

    // Common Variables - sp - sp mono + dipole / sp vs classic_ion 
    gsl_matrix* trans_matrix = NULL;
    gsl_vector* r = gsl_vector_calloc(3);       
    gsl_vector* global_dh = gsl_vector_calloc(4);
    gsl_vector* local_dh  = gsl_vector_calloc(4);

    gsl_matrix* fx_matrix = gsl_matrix_calloc(4,4);
    gsl_matrix* fy_matrix = gsl_matrix_calloc(4,4);
    gsl_matrix* fz_matrix = gsl_matrix_calloc(4,4);     // workspace for tmp saving 1st derivative matrix in the local symm

    double charge_scaler;

    double dxx, dxy, dxz, dyx, dyy, dyz, dzx, dzy, dzz;
    dxx = 0.; dxy = 0.; dxz = 0.; dyx = 0.; dyy = 0.; dyz = 0.; dzx = 0.; dzy = 0.; dzz = 0.;

    gsl_matrix* fxx_matrix = gsl_matrix_calloc(4,4); gsl_matrix* fxy_matrix = gsl_matrix_calloc(4,4); gsl_matrix* fxz_matrix = gsl_matrix_calloc(4,4);
    gsl_matrix* fyx_matrix = gsl_matrix_calloc(4,4); gsl_matrix* fyy_matrix = gsl_matrix_calloc(4,4); gsl_matrix* fyz_matrix = gsl_matrix_calloc(4,4);
    gsl_matrix* fzx_matrix = gsl_matrix_calloc(4,4); gsl_matrix* fzy_matrix = gsl_matrix_calloc(4,4); gsl_matrix* fzz_matrix = gsl_matrix_calloc(4,4);

    gsl_matrix* global_ddh = gsl_matrix_calloc(4,4); gsl_matrix* local_ddh = gsl_matrix_calloc(4,4); gsl_matrix* trans_tmp = gsl_matrix_calloc(4,4);

    // vs Sp Ion Variables
    double evn[4], evm[4];
    double low_index_n, low_index_m;
    double tmp_deriv_x, tmp_deriv_y, tmp_deriv_z;   tmp_deriv_x = 0.; tmp_deriv_y = 0.; tmp_deriv_z = 0.;
    double integral_x, integral_y, integral_z;

    // Refresh Force Saving Buffer;
    for(int n=0;n<number_of_sp_ion;n++) 
    {   memset(&sp_sys->sp_ion[n].elec_force_by_sp[0] ,0.,3*sizeof(double));
        memset(&sp_sys->sp_ion[n].elec_force_by_ion[0],0.,3*sizeof(double));     }
    for(int n=0;n<number_of_classic_ion;n++)    memset(&sp_sys->classic_ion[n].elec_force_by_sp[0],0.,3*sizeof(double));
    
    // Refresh Force WorkSpace
    for(int n=0;n<number_of_sp_ion;n++)
    {   for(int m=0;m<number_of_sp_ion;m++)
        {   gsl_matrix_set_zero(sp_sys->dh_matrix_x[n][m]); gsl_matrix_set_zero(sp_sys->dh_matrix_y[n][m]);  gsl_matrix_set_zero(sp_sys->dh_matrix_z[n][m]);

            gsl_matrix_set_zero(sp_sys->ddh_matrix_xx[n][m]);   gsl_matrix_set_zero(sp_sys->ddh_matrix_xy[n][m]);   gsl_matrix_set_zero(sp_sys->ddh_matrix_xz[n][m]);
            gsl_matrix_set_zero(sp_sys->ddh_matrix_yx[n][m]);   gsl_matrix_set_zero(sp_sys->ddh_matrix_yy[n][m]);   gsl_matrix_set_zero(sp_sys->ddh_matrix_yz[n][m]);
            gsl_matrix_set_zero(sp_sys->ddh_matrix_zx[n][m]);   gsl_matrix_set_zero(sp_sys->ddh_matrix_zy[n][m]);   gsl_matrix_set_zero(sp_sys->ddh_matrix_zz[n][m]);   }
        for(int m=0;m<number_of_classic_ion;m++)
        {   gsl_matrix_set_zero(sp_sys->cla_dh_matrix_x[n][m]); gsl_matrix_set_zero(sp_sys->cla_dh_matrix_y[n][m]); gsl_matrix_set_zero(sp_sys->cla_dh_matrix_z[n][m]); }
    }

            
    ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///
    // Phase 1 - force estimation: 'n'th sp vs 'm'th sp
    len = number_of_sp_ion*number_of_sp_ion;    // total number of for loop cycles
    L = len/numtasks;   R = len%numtasks;
    ista = L*rank + MIN(R,rank);    iend = ista + L - 1;
    if( R > rank )  iend++;
    // loop decomposition with round robin fashion
    
    for(int global_n=ista;global_n<iend+1;global_n++)
    {   idxn = global_n/sp_sys->number_of_sp_ion;   // for 'n'
        idxm = global_n%sp_sys->number_of_sp_ion;   // for 'm'    // cuz 'n == m' always

        if( idxn != idxm )
        {
            gsl_vector_memcpy( r, sp_sys->sp_ion[idxm].core_position );
            gsl_vector_sub( r, sp_sys->sp_ion[idxn].core_position ); // r(sp1->sp2) = r_sp2_core - r_sp1_core
            trans_matrix = sp_cluster_support_get_transformation_matrix( r );   // get transformation matrix w.r.t the above defined geometry
            // Likewise sp1 is seeing sp2 ... i.e., sp1 -> sp2

            // electrostatic charge_scaler ... sp - sp mono was including (sp-core + sp-elec), and sp-sp di calculates sp_elec*sp_elec
            // therefore (sp_elec*sp_elec) should be -> sp_elec*(sp_elec+sp_core)
            charge_scaler = (0.5*sp_sys->sp_ion[idxm].charge_shell+sp_sys->sp_ion[idxm].charge_core)/sp_sys->sp_ion[idxm].charge_shell;

            // !!! set fx matrix
            gsl_matrix_set(fx_matrix,0,1, charge_scaler*sp_cluster_spsp_di_integrator_get_ch_x_12_element(&sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn]) 
                    + 0.5*sp_cluster_spsp_di_integrator_get_sh_x_12_element( sp_sys, &sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn]) );       // set fx 12 element
            gsl_matrix_set(fx_matrix,1,3, charge_scaler*sp_cluster_spsp_di_integrator_get_ch_x_24_element(&sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn])
                    + 0.5*sp_cluster_spsp_di_integrator_get_sh_x_24_element( sp_sys, &sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn]) );       // set fx 24 element
            gsl_matrix_set(fx_matrix,1,0,gsl_matrix_get(fx_matrix,0,1)); gsl_matrix_set(fx_matrix,3,1,gsl_matrix_get(fx_matrix,1,3));           // SIGN IS NOT INVERSED YET
            // !!! set fy matrix
            gsl_matrix_set(fy_matrix,0,2,gsl_matrix_get(fx_matrix,0,1));    // note that DxH SX == DyH SY
            gsl_matrix_set(fy_matrix,2,3,gsl_matrix_get(fx_matrix,1,3));    // note that DxH XZ == DyH YZ
            gsl_matrix_set(fy_matrix,2,0,gsl_matrix_get(fy_matrix,0,2)); gsl_matrix_set(fy_matrix,3,2,gsl_matrix_get(fy_matrix,2,3));
            // !!! set fz matrix
            gsl_matrix_set(fz_matrix,0,0, charge_scaler*sp_cluster_spsp_di_integrator_get_ch_z_11_element(&sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn])
                    + 0.5*sp_cluster_spsp_di_integrator_get_sh_z_11_element( sp_sys, &sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn]) );       // set fz 11 element
            gsl_matrix_set(fz_matrix,0,3, charge_scaler*sp_cluster_spsp_di_integrator_get_ch_z_14_element(&sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn])
                    + 0.5*sp_cluster_spsp_di_integrator_get_sh_z_14_element( sp_sys, &sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn]) );       // set fz 14 element
            gsl_matrix_set(fz_matrix,3,0,gsl_matrix_get(fz_matrix,0,3));                                        // sz == zs
            gsl_matrix_set(fz_matrix,1,1, charge_scaler*sp_cluster_spsp_di_integrator_get_ch_z_2233_element(&sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn])
                    + 0.5*sp_cluster_spsp_di_integrator_get_sh_z_2233_element( sp_sys, &sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn]) );     // set fz 22 element
            gsl_matrix_set(fz_matrix,2,2,gsl_matrix_get(fz_matrix,1,1));                                        // xx == yy
            gsl_matrix_set(fz_matrix,3,3, charge_scaler*sp_cluster_spsp_di_integrator_get_ch_z_44_element(&sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn])
                    + 0.5*sp_cluster_spsp_di_integrator_get_sh_z_44_element( sp_sys, &sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn]) );       // set fz 44 element

            // Caculating global fxyz matrices   w.r.t 'idxm' sp-ion
            for(int m=0;m<4;m++)
            {   for(int n=0;n<4;n++)
                {   for(int h=0;h<4;h++)
                    {   for(int l=0;l<4;l++)
                        {   tmp_deriv_x += gsl_matrix_get(trans_matrix,h,m)*gsl_matrix_get(trans_matrix,l,n)*gsl_matrix_get(fx_matrix,h,l); 
                            tmp_deriv_y += gsl_matrix_get(trans_matrix,h,m)*gsl_matrix_get(trans_matrix,l,n)*gsl_matrix_get(fy_matrix,h,l); 
                            tmp_deriv_z += gsl_matrix_get(trans_matrix,h,m)*gsl_matrix_get(trans_matrix,l,n)*gsl_matrix_get(fz_matrix,h,l); 
                        }   // Calculate local dH matrix element
                            //
                            //  < m | d H | n > == T_hm T_ln < h' | d'H | l' >
                            //
                    }       // Get local expression for < m | d H | n >
                    gsl_vector_set(local_dh,1,tmp_deriv_x); gsl_vector_set(local_dh,2,tmp_deriv_y); gsl_vector_set(local_dh,3,tmp_deriv_z);
                    tmp_deriv_x = 0.; tmp_deriv_y = 0.; tmp_deriv_z = 0.;
                    // Inverse transformation to global symm ...  
                    gsl_blas_dgemv(CblasTrans,1.,trans_matrix,local_dh,0.,global_dh);
                
                    gsl_matrix_set(sp_sys->dh_matrix_workspace_x[idxn][idxm],m,n,gsl_vector_get(global_dh,1));
                    gsl_matrix_set(sp_sys->dh_matrix_workspace_y[idxn][idxm],m,n,gsl_vector_get(global_dh,2));
                    gsl_matrix_set(sp_sys->dh_matrix_workspace_z[idxn][idxm],m,n,gsl_vector_get(global_dh,3));      // Global dH matrix is saved here 'n'th sp ion w.r.t 'm'th sp ion
    
                    /* 
                     * The data structure of dummy_dh_matrix_x,y and z 
                     *
                     * 'sp1' is at row and 'sp2' is at col
                     *
                     * e.g., row2, col4 -> (2,4) means second sp ion is pointing forth sp ion
                     *
                     * i.e., sp2 -> sp4
                     *
                     * Likewise (4,2) means sp4 -> sp2
                     */

                    gsl_vector_set_zero(local_dh); gsl_vector_set_zero(global_dh);
                }
            }   // The section Right Below can be improved by the obtained global < m | d H | n >

            gsl_matrix_set_zero(fx_matrix); gsl_matrix_set_zero(fy_matrix); gsl_matrix_set_zero(fz_matrix);
            gsl_vector_set_zero(r);     gsl_matrix_free(trans_matrix);

        }// end if if( idxn != idxm )
    }// end of sp - sp mono + core
    MPI_Barrier(MPI_COMM_WORLD);

    for(int n=0;n<sp_sys->number_of_sp_ion;n++)
    {   for(int m=0;m<sp_sys->number_of_sp_ion;m++)
        {   if( n != m )
            {   MPI_Allreduce( sp_sys->dh_matrix_workspace_x[n][m]->data, sp_sys->dh_matrix_x[n][m]->data, 16, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
                MPI_Allreduce( sp_sys->dh_matrix_workspace_y[n][m]->data, sp_sys->dh_matrix_y[n][m]->data, 16, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
                MPI_Allreduce( sp_sys->dh_matrix_workspace_z[n][m]->data, sp_sys->dh_matrix_z[n][m]->data, 16, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );   }}}

    // END OF Phase 1 - force estimation: 'n'th sp vs 'm'th sp
    ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///



    ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///
    // Phase 1 - force estimation: 'n'th sp vs 'm'th sp Dipolar Term

    for(int global_n=ista;global_n<iend+1;global_n++)
    {   idxn = global_n/sp_sys->number_of_sp_ion;   // for 'n'
        idxm = global_n%sp_sys->number_of_sp_ion;   // for 'm'    // cuz 'n == m' always

        if( idxn != idxm )
        {
            gsl_vector_memcpy( r, sp_sys->sp_ion[idxm].core_position );
            gsl_vector_sub( r, sp_sys->sp_ion[idxn].core_position ); // r(sp1->sp2) = r_sp2_core - r_sp1_core
            trans_matrix = sp_cluster_support_get_transformation_matrix( r );   // get transformation matrix w.r.t the above defined geometry

            // Estimate local matriced for the derivatives of dipolar terms
            gsl_matrix_set(fxx_matrix,0,0,0.5*(sp_cluster_spsp_quadru_integrator_get_ch_xx_11_element(&sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn]) 
                    + sp_cluster_spsp_quadru_integrator_get_sh_xx_11_element( sp_sys, &sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn]))  );   //xx11
            gsl_matrix_set(fyy_matrix,0,0,gsl_matrix_get(fxx_matrix,0,0));  // yy11 == xx11
            gsl_matrix_set(fxx_matrix,0,3,0.5*(sp_cluster_spsp_quadru_integrator_get_ch_xx_14_element(&sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn]) 
                    + sp_cluster_spsp_quadru_integrator_get_sh_xx_14_element( sp_sys, &sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn]))  );   //xx14
            gsl_matrix_set(fxx_matrix,3,0,gsl_matrix_get(fxx_matrix,0,3));  // symm xx14 == xx41
            gsl_matrix_set(fyy_matrix,0,3,gsl_matrix_get(fxx_matrix,0,3));  // yy14 == xx14
            gsl_matrix_set(fyy_matrix,3,0,gsl_matrix_get(fyy_matrix,0,3));  // symm yy14 == yy41
            gsl_matrix_set(fxx_matrix,1,1,0.5*(sp_cluster_spsp_quadru_integrator_get_ch_xx_22_element(&sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn]) 
                    + sp_cluster_spsp_quadru_integrator_get_sh_xx_22_element( sp_sys, &sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn]))  );   //xx22
            gsl_matrix_set(fyy_matrix,2,2,gsl_matrix_get(fxx_matrix,1,1));  // yy33 == xx22
            gsl_matrix_set(fxx_matrix,2,2,0.5*(sp_cluster_spsp_quadru_integrator_get_ch_xx_33_element(&sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn])  
                    + sp_cluster_spsp_quadru_integrator_get_sh_xx_33_element( sp_sys, &sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn]))  );   //xx33
            gsl_matrix_set(fyy_matrix,1,1,gsl_matrix_get(fxx_matrix,2,2));  // yy22 == xx33
            gsl_matrix_set(fxx_matrix,3,3,0.5*(sp_cluster_spsp_quadru_integrator_get_ch_xx_44_element(&sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn]) 
                    + sp_cluster_spsp_quadru_integrator_get_sh_xx_44_element( sp_sys, &sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn]))  );   //xx44
            gsl_matrix_set(fyy_matrix,3,3,gsl_matrix_get(fxx_matrix,3,3));  // yy44 == xx44                                                                                                                                     fxx fyy !    

            gsl_matrix_set(fxy_matrix,1,2,0.5*(sp_cluster_spsp_quadru_integrator_get_ch_xy_23_element(&sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn]) 
                    + sp_cluster_spsp_quadru_integrator_get_sh_xy_23_element( sp_sys, &sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn]))  );   //xy23
            gsl_matrix_set(fxy_matrix,2,1,gsl_matrix_get(fxy_matrix,1,2));  // symm xy23 == xy32
            gsl_matrix_memcpy(fyx_matrix,fxy_matrix); // fxy == fyx; memcpy(*dest,*src)                                                                                                                                         fxy fyx !   

            gsl_matrix_set(fxz_matrix,0,1,0.5*(sp_cluster_spsp_quadru_integrator_get_ch_xz_12_element(&sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn]) 
                    + sp_cluster_spsp_quadru_integrator_get_sh_xz_12_element( sp_sys, &sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn]))  );   //xz12
            gsl_matrix_set(fxz_matrix,1,0,gsl_matrix_get(fxz_matrix,0,1));  // symm xz12 == xz21
            gsl_matrix_set(fxz_matrix,1,3,0.5*(sp_cluster_spsp_quadru_integrator_get_ch_xz_24_element(&sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn]) 
                    + sp_cluster_spsp_quadru_integrator_get_sh_xz_24_element( sp_sys, &sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn]))  );   //xz24
            gsl_matrix_set(fxz_matrix,3,1,gsl_matrix_get(fxz_matrix,1,3));   // symm xz24 == xz42
            gsl_matrix_memcpy(fzx_matrix,fxz_matrix); // fxz == fzx; memcpy(*dest,*src)                                                                                                                                         fxz fzx !

            gsl_matrix_set(fyz_matrix,0,2,gsl_matrix_get(fxz_matrix,0,1));  // xz12 == yz13
            gsl_matrix_set(fyz_matrix,2,0,gsl_matrix_get(fyz_matrix,0,2));  // symm yz13 == yz31
            gsl_matrix_set(fyz_matrix,2,3,gsl_matrix_get(fxz_matrix,1,3));  // yz34 == xz24
            gsl_matrix_set(fyz_matrix,3,2,gsl_matrix_get(fyz_matrix,2,3));  // symm yz34 == yz43
            gsl_matrix_memcpy(fzy_matrix,fyz_matrix); // fyz == fzy; memcpy(*dest,*src)                                                                                                                                         fyz fzy !

            gsl_matrix_set(fzz_matrix,0,0,0.5*(sp_cluster_spsp_quadru_integrator_get_ch_zz_11_element(&sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn]) 
                    + sp_cluster_spsp_quadru_integrator_get_sh_zz_11_element( sp_sys, &sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn]))  );   //zz11
            gsl_matrix_set(fzz_matrix,0,3,0.5*(sp_cluster_spsp_quadru_integrator_get_ch_zz_14_element(&sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn]) 
                    + sp_cluster_spsp_quadru_integrator_get_sh_zz_14_element( sp_sys, &sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn]))  );   //zz14
            gsl_matrix_set(fzz_matrix,3,0,gsl_matrix_get(fzz_matrix,0,3));  // symm zz41 == zz14
            gsl_matrix_set(fzz_matrix,1,1,0.5*(sp_cluster_spsp_quadru_integrator_get_ch_zz_2233_element(&sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn]) 
                    + sp_cluster_spsp_quadru_integrator_get_sh_zz_2233_element( sp_sys, &sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn]))  );   //zz2233
            gsl_matrix_set(fzz_matrix,2,2,gsl_matrix_get(fzz_matrix,1,1));  // zz22 == zz33
            gsl_matrix_set(fzz_matrix,3,3,0.5*(sp_cluster_spsp_quadru_integrator_get_ch_zz_44_element(&sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn]) 
                    + sp_cluster_spsp_quadru_integrator_get_sh_zz_44_element( sp_sys, &sp_sys->sp_ion[idxm],&sp_sys->sp_ion[idxn]))  );   //zz44      fzz !

            
            // second derivative matrices in local coordinate system
            for(int m=0;m<4;m++)
            {   for(int n=0;n<4;n++)
                {   for(int h=0;h<4;h++)
                    {   for(int l=0;l<4;l++)
                        {   dxx += gsl_matrix_get(trans_matrix,h,m)*gsl_matrix_get(trans_matrix,l,n)*gsl_matrix_get(fxx_matrix,h,l);
                            dxy += gsl_matrix_get(trans_matrix,h,m)*gsl_matrix_get(trans_matrix,l,n)*gsl_matrix_get(fxy_matrix,h,l);
                            dxz += gsl_matrix_get(trans_matrix,h,m)*gsl_matrix_get(trans_matrix,l,n)*gsl_matrix_get(fxz_matrix,h,l);

                            dyx += gsl_matrix_get(trans_matrix,h,m)*gsl_matrix_get(trans_matrix,l,n)*gsl_matrix_get(fyx_matrix,h,l);
                            dyy += gsl_matrix_get(trans_matrix,h,m)*gsl_matrix_get(trans_matrix,l,n)*gsl_matrix_get(fyy_matrix,h,l);
                            dyz += gsl_matrix_get(trans_matrix,h,m)*gsl_matrix_get(trans_matrix,l,n)*gsl_matrix_get(fyz_matrix,h,l);

                            dzx += gsl_matrix_get(trans_matrix,h,m)*gsl_matrix_get(trans_matrix,l,n)*gsl_matrix_get(fzx_matrix,h,l);
                            dzy += gsl_matrix_get(trans_matrix,h,m)*gsl_matrix_get(trans_matrix,l,n)*gsl_matrix_get(fzy_matrix,h,l);
                            dzz += gsl_matrix_get(trans_matrix,h,m)*gsl_matrix_get(trans_matrix,l,n)*gsl_matrix_get(fzz_matrix,h,l);
                        }
                    }
                    gsl_matrix_set(local_ddh,1,1,dxx);  gsl_matrix_set(local_ddh,1,2,dxy);  gsl_matrix_set(local_ddh,1,3,dxz);
                    gsl_matrix_set(local_ddh,2,1,dyx);  gsl_matrix_set(local_ddh,2,2,dyy);  gsl_matrix_set(local_ddh,2,3,dyz);
                    gsl_matrix_set(local_ddh,3,1,dzx);  gsl_matrix_set(local_ddh,3,2,dzy);  gsl_matrix_set(local_ddh,3,3,dzz);
                    // Refresh variables
                    dxx = 0.; dxy = 0.; dxz = 0.; dyx = 0.; dyy = 0.; dyz = 0.; dzx = 0.; dzy = 0.; dzz = 0.;
                    // Inverse to global
                    //
                    //  M = TM'T'
                    //
                    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.,trans_matrix,local_ddh,0.,trans_tmp);  // T^-1 M'
                    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,trans_tmp,trans_matrix,0.,global_ddh);   // (T^-1 M')T

                    // convention ... when 'n'sp sees 'm'sp
                    gsl_matrix_set(sp_sys->ddh_matrix_workspace_xx[idxn][idxm],m,n,gsl_matrix_get(global_ddh,1,1)); gsl_matrix_set(sp_sys->ddh_matrix_workspace_xy[idxn][idxm],m,n,gsl_matrix_get(global_ddh,1,2));
                    gsl_matrix_set(sp_sys->ddh_matrix_workspace_xz[idxn][idxm],m,n,gsl_matrix_get(global_ddh,1,3)); gsl_matrix_set(sp_sys->ddh_matrix_workspace_yx[idxn][idxm],m,n,gsl_matrix_get(global_ddh,2,1));
                    gsl_matrix_set(sp_sys->ddh_matrix_workspace_yy[idxn][idxm],m,n,gsl_matrix_get(global_ddh,2,2)); gsl_matrix_set(sp_sys->ddh_matrix_workspace_yz[idxn][idxm],m,n,gsl_matrix_get(global_ddh,2,3));
                    gsl_matrix_set(sp_sys->ddh_matrix_workspace_zx[idxn][idxm],m,n,gsl_matrix_get(global_ddh,3,1)); gsl_matrix_set(sp_sys->ddh_matrix_workspace_zy[idxn][idxm],m,n,gsl_matrix_get(global_ddh,3,2));
                    gsl_matrix_set(sp_sys->ddh_matrix_workspace_zz[idxn][idxm],m,n,gsl_matrix_get(global_ddh,3,3));

                    gsl_matrix_set_zero(local_ddh);     gsl_matrix_set_zero(global_ddh);
                }
            }
            gsl_matrix_free(trans_matrix);  gsl_vector_set_zero(r);

            gsl_matrix_set_zero(fxx_matrix);    gsl_matrix_set_zero(fxy_matrix);    gsl_matrix_set_zero(fxz_matrix);
            gsl_matrix_set_zero(fyx_matrix);    gsl_matrix_set_zero(fyy_matrix);    gsl_matrix_set_zero(fyz_matrix);
            gsl_matrix_set_zero(fzx_matrix);    gsl_matrix_set_zero(fzy_matrix);    gsl_matrix_set_zero(fzz_matrix);

        }// End of if( idxn != idxm )
    }
    MPI_Barrier(MPI_COMM_WORLD);

    for(int n=0;n<number_of_sp_ion;n++)
    {   for(int m=0;m<number_of_sp_ion;m++)
        {   if( n != m )
            {   MPI_Allreduce(sp_sys->ddh_matrix_workspace_xx[n][m]->data,sp_sys->ddh_matrix_xx[n][m]->data,16,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                MPI_Allreduce(sp_sys->ddh_matrix_workspace_xy[n][m]->data,sp_sys->ddh_matrix_xy[n][m]->data,16,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                MPI_Allreduce(sp_sys->ddh_matrix_workspace_xz[n][m]->data,sp_sys->ddh_matrix_xz[n][m]->data,16,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

                MPI_Allreduce(sp_sys->ddh_matrix_workspace_yx[n][m]->data,sp_sys->ddh_matrix_yx[n][m]->data,16,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                MPI_Allreduce(sp_sys->ddh_matrix_workspace_yy[n][m]->data,sp_sys->ddh_matrix_yy[n][m]->data,16,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                MPI_Allreduce(sp_sys->ddh_matrix_workspace_yz[n][m]->data,sp_sys->ddh_matrix_yz[n][m]->data,16,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

                MPI_Allreduce(sp_sys->ddh_matrix_workspace_zx[n][m]->data,sp_sys->ddh_matrix_zx[n][m]->data,16,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                MPI_Allreduce(sp_sys->ddh_matrix_workspace_zy[n][m]->data,sp_sys->ddh_matrix_zy[n][m]->data,16,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                MPI_Allreduce(sp_sys->ddh_matrix_workspace_zz[n][m]->data,sp_sys->ddh_matrix_zz[n][m]->data,16,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            }
        }
    }


    ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///
    // Phase 2 - force estimation: 'n'th sp vs 'm'th classic ion 
    len = number_of_sp_ion*number_of_classic_ion;    // total number of for loop cycles
    L = len/numtasks;   R = len%numtasks;
    ista = L*rank + MIN(R,rank);    iend = ista + L - 1;
    if( R > rank )  iend++;
    // loop decomposition with round robin fashion
    
    for(int global_n=ista;global_n<iend+1;global_n++)
    {   
        if( number_of_sp_ion == number_of_classic_ion )
        {   idxn = global_n/sp_sys->number_of_sp_ion;
            idxm = global_n%sp_sys->number_of_classic_ion;                                      } 
        else
        {
            idxn = global_n/sp_sys->number_of_classic_ion;
            idxm = global_n - (global_n/sp_sys->number_of_classic_ion)*sp_sys->number_of_classic_ion;   }
        // idxn for sp-ion index ... idxm for classic-ion index
        
        gsl_vector_memcpy( r, sp_sys->classic_ion[idxm].core_position );
        gsl_vector_sub( r, sp_sys->sp_ion[idxn].core_position ); // r(sp1->sp2) = r_sp2_core - r_sp1_core
        trans_matrix = sp_cluster_support_get_transformation_matrix( r );   // get transformation matrix w.r.t the above defined geometry

        // vs classic ion ... if the charge is positive then excluding short_range interaction
        //if( sp_sys->classic_ion[idxm].charge_core < 0. ) // SR -- only when the charge is negative
        //if( sp_sys->classic_ion[idxm].if_shell == SP_SYSTEM_TRUE )	// IF THE ATOM IS SHELL, THEN USE SHORT_RANGE INTERACTION 
	if( sp_sys->classic_ion[idxm].if_shell == SP_SYSTEM_TRUE \
		|| (sp_sys->classic_ion[idxm].if_shell == SP_SYSTEM_FALSE && sp_sys->classic_ion[idxm].if_has_shell == SP_SYSTEM_FALSE) )	// IF IT IS SHELL MODEL, THEN INCLUDE SHORT_RANGE INTERACTION 
        {
            // !!! set fx matrix
            gsl_matrix_set(fx_matrix,0,1,sp_cluster_integrator_get_ch_x_12_element(&sp_sys->sp_ion[idxn],&sp_sys->classic_ion[idxm]) 
                    + sp_cluster_integrator_get_sh_x_12_element( sp_sys, &sp_sys->sp_ion[idxn],&sp_sys->classic_ion[idxm]));       // set fx 12 element
            gsl_matrix_set(fx_matrix,1,3,sp_cluster_integrator_get_ch_x_24_element(&sp_sys->sp_ion[idxn],&sp_sys->classic_ion[idxm])
                    + sp_cluster_integrator_get_sh_x_24_element( sp_sys, &sp_sys->sp_ion[idxn],&sp_sys->classic_ion[idxm]));       // set fx 24 element
            // anyway its hermitian matrix
            gsl_matrix_set(fx_matrix,1,0,gsl_matrix_get(fx_matrix,0,1)); gsl_matrix_set(fx_matrix,3,1,gsl_matrix_get(fx_matrix,1,3));

            // !!! set fy matrix
            gsl_matrix_set(fy_matrix,0,2,gsl_matrix_get(fx_matrix,0,1));    // note that DxH SX == DyH SY
            gsl_matrix_set(fy_matrix,2,3,gsl_matrix_get(fx_matrix,1,3));    // note that DxH XZ == DyH YZ
            // anayway its hermitian matrix
            gsl_matrix_set(fy_matrix,2,0,gsl_matrix_get(fy_matrix,0,2)); gsl_matrix_set(fy_matrix,3,2,gsl_matrix_get(fy_matrix,2,3));

            // !!! set fz matrix
            gsl_matrix_set(fz_matrix,0,0,sp_cluster_integrator_get_ch_z_11_element(&sp_sys->sp_ion[idxn],&sp_sys->classic_ion[idxm])
                    + sp_cluster_integrator_get_sh_z_11_element( sp_sys, &sp_sys->sp_ion[idxn],&sp_sys->classic_ion[idxm]));       // set fz 11 element
            gsl_matrix_set(fz_matrix,0,3,sp_cluster_integrator_get_ch_z_14_element(&sp_sys->sp_ion[idxn],&sp_sys->classic_ion[idxm])
                    + sp_cluster_integrator_get_sh_z_14_element( sp_sys, &sp_sys->sp_ion[idxn],&sp_sys->classic_ion[idxm]));       // set fz 14 element
            gsl_matrix_set(fz_matrix,3,0,gsl_matrix_get(fz_matrix,0,3));                                        // sz == zs
            gsl_matrix_set(fz_matrix,1,1,sp_cluster_integrator_get_ch_z_2233_element(&sp_sys->sp_ion[idxn],&sp_sys->classic_ion[idxm])
                    + sp_cluster_integrator_get_sh_z_2233_element( sp_sys, &sp_sys->sp_ion[idxn],&sp_sys->classic_ion[idxm]));     // set fz 22 element
            gsl_matrix_set(fz_matrix,2,2,gsl_matrix_get(fz_matrix,1,1));                                        // xx == yy
            gsl_matrix_set(fz_matrix,3,3,sp_cluster_integrator_get_ch_z_44_element(&sp_sys->sp_ion[idxn],&sp_sys->classic_ion[idxm])
                    + sp_cluster_integrator_get_sh_z_44_element( sp_sys, &sp_sys->sp_ion[idxn],&sp_sys->classic_ion[idxm]));       // set fz 44 element
        }
        else if( sp_sys->classic_ion[idxm].if_shell == SP_SYSTEM_FALSE \
		&& sp_sys->classic_ion[idxm].if_has_shell == SP_SYSTEM_TRUE )	// CORE and HAS SHELL
        {
            // !!! set fx matrix
            gsl_matrix_set(fx_matrix,0,1,sp_cluster_integrator_get_ch_x_12_element(&sp_sys->sp_ion[idxn],&sp_sys->classic_ion[idxm]));       // set fx 12 element
            gsl_matrix_set(fx_matrix,1,3,sp_cluster_integrator_get_ch_x_24_element(&sp_sys->sp_ion[idxn],&sp_sys->classic_ion[idxm]));       // set fx 24 element
            // anyway its hermitian matrix
            gsl_matrix_set(fx_matrix,1,0,gsl_matrix_get(fx_matrix,0,1)); gsl_matrix_set(fx_matrix,3,1,gsl_matrix_get(fx_matrix,1,3));

            // !!! set fy matrix
            gsl_matrix_set(fy_matrix,0,2,gsl_matrix_get(fx_matrix,0,1));    // note that DxH SX == DyH SY
            gsl_matrix_set(fy_matrix,2,3,gsl_matrix_get(fx_matrix,1,3));    // note that DxH XZ == DyH YZ
            // anayway its hermitian matrix
            gsl_matrix_set(fy_matrix,2,0,gsl_matrix_get(fy_matrix,0,2)); gsl_matrix_set(fy_matrix,3,2,gsl_matrix_get(fy_matrix,2,3));

            // !!! set fz matrix
            gsl_matrix_set(fz_matrix,0,0,sp_cluster_integrator_get_ch_z_11_element(&sp_sys->sp_ion[idxn],&sp_sys->classic_ion[idxm]));       // set fz 11 element
            gsl_matrix_set(fz_matrix,0,3,sp_cluster_integrator_get_ch_z_14_element(&sp_sys->sp_ion[idxn],&sp_sys->classic_ion[idxm]));       // set fz 14 element
            gsl_matrix_set(fz_matrix,3,0,gsl_matrix_get(fz_matrix,0,3));                                        // sz == zs
            gsl_matrix_set(fz_matrix,1,1,sp_cluster_integrator_get_ch_z_2233_element(&sp_sys->sp_ion[idxn],&sp_sys->classic_ion[idxm]));     // set fz 22 element
            gsl_matrix_set(fz_matrix,2,2,gsl_matrix_get(fz_matrix,1,1));                                        // xx == yy
            gsl_matrix_set(fz_matrix,3,3,sp_cluster_integrator_get_ch_z_44_element(&sp_sys->sp_ion[idxn],&sp_sys->classic_ion[idxm]));       // set fz 44 element
        }

        // Caculating global fxyz matrices   wrt 'idxm'th classic ion
        for(int m=0;m<4;m++)
        {   for(int n=0;n<4;n++)
            {   for(int h=0;h<4;h++)
                {   for(int l=0;l<4;l++)
                    {   tmp_deriv_x += gsl_matrix_get(trans_matrix,h,m)*gsl_matrix_get(trans_matrix,l,n)*gsl_matrix_get(fx_matrix,h,l); 
                        tmp_deriv_y += gsl_matrix_get(trans_matrix,h,m)*gsl_matrix_get(trans_matrix,l,n)*gsl_matrix_get(fy_matrix,h,l); 
                        tmp_deriv_z += gsl_matrix_get(trans_matrix,h,m)*gsl_matrix_get(trans_matrix,l,n)*gsl_matrix_get(fz_matrix,h,l); 
                    }   // Calculate local dH matrix element
                        //
                        //  < m | d H | n > == T_hm T_ln < h' | d'H | l' >
                        //
                }       // Get local expression for < m | d H | n >
                gsl_vector_set(local_dh,1,tmp_deriv_x); gsl_vector_set(local_dh,2,tmp_deriv_y); gsl_vector_set(local_dh,3,tmp_deriv_z);
                tmp_deriv_x = 0.; tmp_deriv_y = 0.; tmp_deriv_z = 0.;
                // Inverse transformation to global symm ...  
                gsl_blas_dgemv(CblasTrans,1.,trans_matrix,local_dh,0.,global_dh);
            
                gsl_matrix_set(sp_sys->cla_dh_matrix_workspace_x[idxn][idxm],m,n,gsl_vector_get(global_dh,1));
                gsl_matrix_set(sp_sys->cla_dh_matrix_workspace_y[idxn][idxm],m,n,gsl_vector_get(global_dh,2));
                gsl_matrix_set(sp_sys->cla_dh_matrix_workspace_z[idxn][idxm],m,n,gsl_vector_get(global_dh,3));      // Global dH matrix is saved here 'n'th sp ion w.r.t 'm'th classic ion

                /* 
                 * The data structure of dummy_dh_matrix_x,y and z 
                 *
                 * 'sp1' is at row and 'sp2' is at col
                 *
                 * e.g., row2, col4 -> (2,4) means second sp ion is pointing forth sp ion
                 *
                 * i.e., sp2 -> sp4
                 *
                 * Likewise (4,2) means sp4 -> sp2
                 */

                gsl_vector_set_zero(local_dh); gsl_vector_set_zero(global_dh);
            }
        }   // The section Right Below can be improved by the obtained global < m | d H | n >

        gsl_matrix_set_zero(fx_matrix); gsl_matrix_set_zero(fy_matrix); gsl_matrix_set_zero(fz_matrix);
        gsl_vector_set_zero(r); gsl_matrix_free(trans_matrix);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Allrudce for sp-ion vs cla-ion
    for(int n=0;n<sp_sys->number_of_sp_ion;n++)
    {   for(int m=0;m<sp_sys->number_of_classic_ion;m++)
        {   MPI_Allreduce( sp_sys->cla_dh_matrix_workspace_x[n][m]->data, sp_sys->cla_dh_matrix_x[n][m]->data, 16, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
            MPI_Allreduce( sp_sys->cla_dh_matrix_workspace_y[n][m]->data, sp_sys->cla_dh_matrix_y[n][m]->data, 16, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
            MPI_Allreduce( sp_sys->cla_dh_matrix_workspace_z[n][m]->data, sp_sys->cla_dh_matrix_z[n][m]->data, 16, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );   }}

    MPI_Barrier(MPI_COMM_WORLD);
    ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     End Of Force Calc Prep .. all matriced estimated
    
    ///     Derivative Estimation Confirmed 31072019

    ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     End Of Force Calc Prep .. all matriced estimated

    // FORCE CALCULATION

    //double dip[number_of_sp_ion][3]; for(int i=0;i<number_of_sp_ion;i++)    memset(dip[i],0.,3*sizeof(double));

    for(int n=0;n<number_of_sp_ion;n++)
    {
        low_index_n = sp_cluster_support_get_lowest_state(sp_sys->sp_ion[n].eigen_value);    // Get Evec Low Index of 'n' sp
        evn[0] = gsl_matrix_get(sp_sys->sp_ion[n].eigen_vector,0,low_index_n);  evn[1] = gsl_matrix_get(sp_sys->sp_ion[n].eigen_vector,1,low_index_n);
        evn[2] = gsl_matrix_get(sp_sys->sp_ion[n].eigen_vector,2,low_index_n);  evn[3] = gsl_matrix_get(sp_sys->sp_ion[n].eigen_vector,3,low_index_n);      // Get Evec of 'n' sp

        for(int m=0;m<number_of_sp_ion;m++)               // sp - sp Force Calculation ... Convention: Derivative of 'n' sp w.r.t 'm' sp 
        {
            if( n != m )
            {

            low_index_m = sp_cluster_support_get_lowest_state(sp_sys->sp_ion[m].eigen_value);    // Get Evec Low Index of 'm' sp
            evm[0] = gsl_matrix_get(sp_sys->sp_ion[m].eigen_vector,0,low_index_m);  evm[1] = gsl_matrix_get(sp_sys->sp_ion[m].eigen_vector,1,low_index_m);
            evm[2] = gsl_matrix_get(sp_sys->sp_ion[m].eigen_vector,2,low_index_m);  evm[3] = gsl_matrix_get(sp_sys->sp_ion[m].eigen_vector,3,low_index_m);      // Get Evec of 'm' sp

            integral_x = sp_cluster_integrator_get_x_12(&sp_sys->sp_ion[m]);    integral_y = sp_cluster_integrator_get_y_13(&sp_sys->sp_ion[m]);    integral_z = sp_cluster_integrator_get_z_14(&sp_sys->sp_ion[m]);        // Get position operator integrals

            for(int k=0;k<4;k++)
            {
                for(int l=0;l<4;l++)
                {   
                    // Monopole Contribution
                    tmp_deriv_x += evn[k]*evn[l]*gsl_matrix_get(sp_sys->dh_matrix_x[n][m],k,l);   
                    tmp_deriv_y += evn[k]*evn[l]*gsl_matrix_get(sp_sys->dh_matrix_y[n][m],k,l);   
                    tmp_deriv_z += evn[k]*evn[l]*gsl_matrix_get(sp_sys->dh_matrix_z[n][m],k,l);
                    
                    // Dipole Contribution

                    tmp_deriv_x += evn[k]*evn[l]*(
                                    -2.*evm[0]*evm[1]*integral_x*gsl_matrix_get(sp_sys->ddh_matrix_xx[n][m],k,l)
                                    -2.*evm[0]*evm[2]*integral_y*gsl_matrix_get(sp_sys->ddh_matrix_xy[n][m],k,l)
                                    -2.*evm[0]*evm[3]*integral_z*gsl_matrix_get(sp_sys->ddh_matrix_xz[n][m],k,l) ); 

                    tmp_deriv_y += evn[k]*evn[l]*(
                                    -2.*evm[0]*evm[1]*integral_x*gsl_matrix_get(sp_sys->ddh_matrix_yx[n][m],k,l)
                                    -2.*evm[0]*evm[2]*integral_y*gsl_matrix_get(sp_sys->ddh_matrix_yy[n][m],k,l)
                                    -2.*evm[0]*evm[3]*integral_z*gsl_matrix_get(sp_sys->ddh_matrix_yz[n][m],k,l) ); 

                    tmp_deriv_z += evn[k]*evn[l]*(
                                    -2.*evm[0]*evm[1]*integral_x*gsl_matrix_get(sp_sys->ddh_matrix_zx[n][m],k,l)
                                    -2.*evm[0]*evm[2]*integral_y*gsl_matrix_get(sp_sys->ddh_matrix_zy[n][m],k,l)
                                    -2.*evm[0]*evm[3]*integral_z*gsl_matrix_get(sp_sys->ddh_matrix_zz[n][m],k,l) ); 
                }//l       
            }//k

            // Force Update ... on 'n' sp by 'm' sp
            sp_sys->sp_ion[n].elec_force_by_sp[0] -= tmp_deriv_x;   sp_sys->sp_ion[m].elec_force_by_sp[0] += tmp_deriv_x; 
            sp_sys->sp_ion[n].elec_force_by_sp[1] -= tmp_deriv_y;   sp_sys->sp_ion[m].elec_force_by_sp[1] += tmp_deriv_y; 
            sp_sys->sp_ion[n].elec_force_by_sp[2] -= tmp_deriv_z;   sp_sys->sp_ion[m].elec_force_by_sp[2] += tmp_deriv_z;

            // Refresg WS
            tmp_deriv_x = 0.;   tmp_deriv_y = 0.; tmp_deriv_z = 0.;

            }// End of If( n!=m )

        }// End of sp 'm' loop



        for(int m=0;m<number_of_classic_ion;m++)            // sp - classic_ion Force Calculation ... Convention: Derivative of 'n' sp w.r.t 'm' classic ion
        {
            for(int k=0;k<4;k++)
            {   for(int l=0;l<4;l++)
                {   tmp_deriv_x += evn[k]*evn[l]*gsl_matrix_get(sp_sys->cla_dh_matrix_x[n][m],k,l);
                    tmp_deriv_y += evn[k]*evn[l]*gsl_matrix_get(sp_sys->cla_dh_matrix_y[n][m],k,l);
                    tmp_deriv_z += evn[k]*evn[l]*gsl_matrix_get(sp_sys->cla_dh_matrix_z[n][m],k,l);     }
            }
            
            //  Force Update ... on 'n' sp & on 'm' classic ion
            sp_sys->sp_ion[n].elec_force_by_ion[0] -= tmp_deriv_x;  sp_sys->classic_ion[m].elec_force_by_sp[0] += tmp_deriv_x;
            sp_sys->sp_ion[n].elec_force_by_ion[1] -= tmp_deriv_y;  sp_sys->classic_ion[m].elec_force_by_sp[1] += tmp_deriv_y;
            sp_sys->sp_ion[n].elec_force_by_ion[2] -= tmp_deriv_z;  sp_sys->classic_ion[m].elec_force_by_sp[2] += tmp_deriv_z;

            // Refresh WS
            tmp_deriv_x = 0.;   tmp_deriv_y = 0.; tmp_deriv_z = 0.;

        }// End of ion 'm' loop

    }// End of sp 'n' loop
 
    // UP TO HERE ... THE FORCE CONTRIBUTION COMPUTED ARE THE TERMS, INDEPENDENT ON THE DERIVATIVES OF EVECS 

    
    // EVEC DERIVATIVE CORRECTION TERM !!!!!! 2019/09/26  VER 1.2 

    test__test( sp_sys, rank, numtasks );

    // EVEC DERIVATIVE CORRECTION TERM - VERY IMPORTANT !!!

    integral_x = sp_cluster_integrator_get_x_12(&sp_sys->sp_ion[0]);    integral_y = sp_cluster_integrator_get_y_13(&sp_sys->sp_ion[0]);    integral_z = sp_cluster_integrator_get_z_14(&sp_sys->sp_ion[0]);        

    // evec contri - sp part
    for(int l=0;l<number_of_sp_ion;l++)
    {
        for(int n=0;n<number_of_sp_ion;n++)
        {
            low_index_n = sp_cluster_support_get_lowest_state(sp_sys->sp_ion[n].eigen_value);    // Get Evec Low Index of 'n' sp
            evn[0] = gsl_matrix_get(sp_sys->sp_ion[n].eigen_vector,0,low_index_n);  evn[1] = gsl_matrix_get(sp_sys->sp_ion[n].eigen_vector,1,low_index_n);
            evn[2] = gsl_matrix_get(sp_sys->sp_ion[n].eigen_vector,2,low_index_n);  evn[3] = gsl_matrix_get(sp_sys->sp_ion[n].eigen_vector,3,low_index_n);      // Get Evec of 'n' sp

            for(int m=0;m<number_of_sp_ion;m++)               // sp - sp Force Calculation ... Convention: Derivative of 'n' sp w.r.t 'm' sp 
            {
                if( n != m )
                {
                    low_index_m = sp_cluster_support_get_lowest_state(sp_sys->sp_ion[m].eigen_value);    // Get Evec Low Index of 'm' sp
                    evm[0] = gsl_matrix_get(sp_sys->sp_ion[m].eigen_vector,0,low_index_m);  evm[1] = gsl_matrix_get(sp_sys->sp_ion[m].eigen_vector,1,low_index_m);
                    evm[2] = gsl_matrix_get(sp_sys->sp_ion[m].eigen_vector,2,low_index_m);  evm[3] = gsl_matrix_get(sp_sys->sp_ion[m].eigen_vector,3,low_index_m);      // Get Evec of 'm' sp

                    //integral_x = sp_cluster_integrator_get_x_12(&sp_sys->sp_ion[m]);    integral_y = sp_cluster_integrator_get_y_13(&sp_sys->sp_ion[m]);    integral_z = sp_cluster_integrator_get_z_14(&sp_sys->sp_ion[m]);        
                    // Get position operator integrals

                    for(int u=0;u<4;u++)
                    {
                        for(int v=0;v<4;v++)
                        {
                            tmp_deriv_x += ( -2.*evn[u]*evn[v]*integral_x*gsl_matrix_get(sp_sys->scf_dh_x_matrix_vs_sp_ion_dipole[n][m],u,v)*(sp_sys->deriv_evec_sp[l][m][0][0]*evm[1]+evm[0]*sp_sys->deriv_evec_sp[l][m][1][0])
                                             -2.*evn[u]*evn[v]*integral_y*gsl_matrix_get(sp_sys->scf_dh_y_matrix_vs_sp_ion_dipole[n][m],u,v)*(sp_sys->deriv_evec_sp[l][m][0][0]*evm[2]+evm[0]*sp_sys->deriv_evec_sp[l][m][2][0])
                                             -2.*evn[u]*evn[v]*integral_z*gsl_matrix_get(sp_sys->scf_dh_z_matrix_vs_sp_ion_dipole[n][m],u,v)*(sp_sys->deriv_evec_sp[l][m][0][0]*evm[3]+evm[0]*sp_sys->deriv_evec_sp[l][m][3][0]) );

                            tmp_deriv_y += ( -2.*evn[u]*evn[v]*integral_x*gsl_matrix_get(sp_sys->scf_dh_x_matrix_vs_sp_ion_dipole[n][m],u,v)*(sp_sys->deriv_evec_sp[l][m][0][1]*evm[1]+evm[0]*sp_sys->deriv_evec_sp[l][m][1][1])
                                             -2.*evn[u]*evn[v]*integral_y*gsl_matrix_get(sp_sys->scf_dh_y_matrix_vs_sp_ion_dipole[n][m],u,v)*(sp_sys->deriv_evec_sp[l][m][0][1]*evm[2]+evm[0]*sp_sys->deriv_evec_sp[l][m][2][1])
                                             -2.*evn[u]*evn[v]*integral_z*gsl_matrix_get(sp_sys->scf_dh_z_matrix_vs_sp_ion_dipole[n][m],u,v)*(sp_sys->deriv_evec_sp[l][m][0][1]*evm[3]+evm[0]*sp_sys->deriv_evec_sp[l][m][3][1]) );

                            tmp_deriv_z += ( -2.*evn[u]*evn[v]*integral_x*gsl_matrix_get(sp_sys->scf_dh_x_matrix_vs_sp_ion_dipole[n][m],u,v)*(sp_sys->deriv_evec_sp[l][m][0][2]*evm[1]+evm[0]*sp_sys->deriv_evec_sp[l][m][1][2])
                                             -2.*evn[u]*evn[v]*integral_y*gsl_matrix_get(sp_sys->scf_dh_y_matrix_vs_sp_ion_dipole[n][m],u,v)*(sp_sys->deriv_evec_sp[l][m][0][2]*evm[2]+evm[0]*sp_sys->deriv_evec_sp[l][m][2][2])
                                             -2.*evn[u]*evn[v]*integral_z*gsl_matrix_get(sp_sys->scf_dh_z_matrix_vs_sp_ion_dipole[n][m],u,v)*(sp_sys->deriv_evec_sp[l][m][0][2]*evm[3]+evm[0]*sp_sys->deriv_evec_sp[l][m][3][2]) );
                        }//v
                    }//u
                }//end if
            }//m
        }//n

        sp_sys->sp_ion[l].elec_force_by_sp[0] += -tmp_deriv_x;
        sp_sys->sp_ion[l].elec_force_by_sp[1] += -tmp_deriv_y;
        sp_sys->sp_ion[l].elec_force_by_sp[2] += -tmp_deriv_z;

        // Refresh WS
        tmp_deriv_x = 0.;   tmp_deriv_y = 0.; tmp_deriv_z = 0.;

    }//l


    // evec contri - classic part
    // derivative contribution by eigenvector derivatives.
    for(int l=0;l<number_of_classic_ion;l++)
    {
        for(int n=0;n<number_of_sp_ion;n++)
        {
            low_index_n = sp_cluster_support_get_lowest_state(sp_sys->sp_ion[n].eigen_value);    // Get Evec Low Index of 'n' sp
            evn[0] = gsl_matrix_get(sp_sys->sp_ion[n].eigen_vector,0,low_index_n);  evn[1] = gsl_matrix_get(sp_sys->sp_ion[n].eigen_vector,1,low_index_n);
            evn[2] = gsl_matrix_get(sp_sys->sp_ion[n].eigen_vector,2,low_index_n);  evn[3] = gsl_matrix_get(sp_sys->sp_ion[n].eigen_vector,3,low_index_n);      // Get Evec of 'n' sp

            for(int m=0;m<number_of_sp_ion;m++)               // sp - sp Force Calculation ... Convention: Derivative of 'n' sp w.r.t 'm' sp 
            {
                if( n != m )
                {
                    low_index_m = sp_cluster_support_get_lowest_state(sp_sys->sp_ion[m].eigen_value);    // Get Evec Low Index of 'm' sp
                    evm[0] = gsl_matrix_get(sp_sys->sp_ion[m].eigen_vector,0,low_index_m);  evm[1] = gsl_matrix_get(sp_sys->sp_ion[m].eigen_vector,1,low_index_m);
                    evm[2] = gsl_matrix_get(sp_sys->sp_ion[m].eigen_vector,2,low_index_m);  evm[3] = gsl_matrix_get(sp_sys->sp_ion[m].eigen_vector,3,low_index_m);      // Get Evec of 'm' sp

                    //integral_x = sp_cluster_integrator_get_x_12(&sp_sys->sp_ion[m]);    integral_y = sp_cluster_integrator_get_y_13(&sp_sys->sp_ion[m]);    integral_z = sp_cluster_integrator_get_z_14(&sp_sys->sp_ion[m]);        
                    // Get position operator integrals

                    for(int u=0;u<4;u++)
                    {
                        for(int v=0;v<4;v++)
                        {
                            // name convention move when move 'l' then the changes in evec of 'm' than energy change?
                            tmp_deriv_x += ( -2.*evn[u]*evn[v]*integral_x*gsl_matrix_get(sp_sys->scf_dh_x_matrix_vs_sp_ion_dipole[n][m],u,v)*(sp_sys->deriv_evec_cla[l][m][0][0]*evm[1]+evm[0]*sp_sys->deriv_evec_cla[l][m][1][0])
                                             -2.*evn[u]*evn[v]*integral_y*gsl_matrix_get(sp_sys->scf_dh_y_matrix_vs_sp_ion_dipole[n][m],u,v)*(sp_sys->deriv_evec_cla[l][m][0][0]*evm[2]+evm[0]*sp_sys->deriv_evec_cla[l][m][2][0])
                                             -2.*evn[u]*evn[v]*integral_z*gsl_matrix_get(sp_sys->scf_dh_z_matrix_vs_sp_ion_dipole[n][m],u,v)*(sp_sys->deriv_evec_cla[l][m][0][0]*evm[3]+evm[0]*sp_sys->deriv_evec_cla[l][m][3][0]) );

                            tmp_deriv_y += ( -2.*evn[u]*evn[v]*integral_x*gsl_matrix_get(sp_sys->scf_dh_x_matrix_vs_sp_ion_dipole[n][m],u,v)*(sp_sys->deriv_evec_cla[l][m][0][1]*evm[1]+evm[0]*sp_sys->deriv_evec_cla[l][m][1][1])
                                             -2.*evn[u]*evn[v]*integral_y*gsl_matrix_get(sp_sys->scf_dh_y_matrix_vs_sp_ion_dipole[n][m],u,v)*(sp_sys->deriv_evec_cla[l][m][0][1]*evm[2]+evm[0]*sp_sys->deriv_evec_cla[l][m][2][1])
                                             -2.*evn[u]*evn[v]*integral_z*gsl_matrix_get(sp_sys->scf_dh_z_matrix_vs_sp_ion_dipole[n][m],u,v)*(sp_sys->deriv_evec_cla[l][m][0][1]*evm[3]+evm[0]*sp_sys->deriv_evec_cla[l][m][3][1]) );

                            tmp_deriv_z += ( -2.*evn[u]*evn[v]*integral_x*gsl_matrix_get(sp_sys->scf_dh_x_matrix_vs_sp_ion_dipole[n][m],u,v)*(sp_sys->deriv_evec_cla[l][m][0][2]*evm[1]+evm[0]*sp_sys->deriv_evec_cla[l][m][1][2])
                                             -2.*evn[u]*evn[v]*integral_y*gsl_matrix_get(sp_sys->scf_dh_y_matrix_vs_sp_ion_dipole[n][m],u,v)*(sp_sys->deriv_evec_cla[l][m][0][2]*evm[2]+evm[0]*sp_sys->deriv_evec_cla[l][m][2][2])
                                             -2.*evn[u]*evn[v]*integral_z*gsl_matrix_get(sp_sys->scf_dh_z_matrix_vs_sp_ion_dipole[n][m],u,v)*(sp_sys->deriv_evec_cla[l][m][0][2]*evm[3]+evm[0]*sp_sys->deriv_evec_cla[l][m][3][2]) );
                        }//v
                    }//u
                }//end if
            }//m
        }//n

        //printf("%lf\t\t%lf\t\t%lf\n",tmp_deriv_x,tmp_deriv_y,tmp_deriv_z);
        sp_sys->classic_ion[l].elec_force_by_sp[0] += -tmp_deriv_x;
        sp_sys->classic_ion[l].elec_force_by_sp[1] += -tmp_deriv_y;
        sp_sys->classic_ion[l].elec_force_by_sp[2] += -tmp_deriv_z;

        // Refresh WS
        tmp_deriv_x = 0.;   tmp_deriv_y = 0.; tmp_deriv_z = 0.;

    }//l
 

    ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///

    ///     END OF FORCE COMPUTATION

    // detach memory
    gsl_vector_free(local_dh); gsl_vector_free(global_dh);

    gsl_matrix_free(fx_matrix); gsl_matrix_free(fy_matrix); gsl_matrix_free(fz_matrix);
    gsl_vector_free(r);



    // 2nd derivative test variables detach
    free(fxx_matrix);   free(fxy_matrix);   free(fxz_matrix);
    free(fyx_matrix);   free(fyy_matrix);   free(fyz_matrix);
    free(fzx_matrix);   free(fzy_matrix);   free(fzz_matrix);
    
    gsl_matrix_free(global_ddh);    gsl_matrix_free(local_ddh);

    gsl_matrix_free(trans_tmp);
    
    MPI_Barrier(MPI_COMM_WORLD);
    return;
}






/************************************* Classical Treatments ************************************************/
// Tag:classic_energy
void sp_cluster_system_get_classic_energy( sp_cluster_system* sp_sys )
{
	// 2019 - 03 - 03 Sunday Repulsive interaction on sp_core vs external ions
	//
	// in this module Born-Mayer Interaction will be added on sp-core vs external ions
	// and this result will be saved on 'classic_ion -> classic energy' & 'sp_ion -> classic energy'

    // 2019 - 07 - 25 Update sp_vs_sp // sp_vs_ion // ion_vs_ion

    /* energy equation is given by
     *
     *                     q_sp_core*q_classic_ion
     *   ------------------------------------------------------------
     *    [(x_ion - x_sp)^2 + (y_ion - y_sp)^2 +(z_ion - z_sp)^2]^1/2
     *
     *
     *	 A_ion_sp * Exp( -[(x_ion - x_sp)^2 + (y_ion - y_sp)^2 +(z_ion - z_sp)^2]^1/2 / R_ion_sp );
     *
     *
     *   Aij Exp[ - [(xi - xj)^2 + (yi - yj)^2 +(zi - zj)^2]^1/2 / rhoij ] - Cij/([(xi - xj)^2 + (yi - yj)^2 +(zi - zj)^2]^1/2)^6
     *
     */

    // COMMON VARIABLES
    double x1,y1,z1, x2,y2,z2;  // coordinate tmp
    double short_range_a1, short_range_a2, short_range_r1, short_range_r2; // BM-like potential short_range A rho param tmp
    double vdw1, vdw2;  // van der waals param for (1/r^6) term among classical ions
    double charge1, charge2;    // charge tmp
    double dist; // interionic distance tmp
    double energy_tmp = 0.;

    /// Overview ...
    /*
     *  Prep    - Init all configs ( set into zero )
     *
     *  Phase 1 - sp_ion_core vs sp_ion_core
     *
     *  Phase 2 - classic_ion_core vs classic_ion_core
     *
     *  Phase 3 - sp_ion_core vs classic_ion_core
     */

    /// Prep 
    for(int n=0;n<sp_sys->number_of_sp_ion;n++)
    {   sp_sys->sp_ion[n].classic_energy_by_sp_core = 0.;
        sp_sys->sp_ion[n].classic_energy_by_ion_core = 0.;  }
    for(int n=0;n<sp_sys->number_of_classic_ion;n++)
    {   sp_sys->classic_ion[n].classic_energy_by_sp_core = 0.;
        sp_sys->classic_ion[n].classic_energy_by_ion_core = 0.; }
    // refreshing energy 

    /// Phase 1: sp_core(n) vs sp_core(m) interaction

	int SP_SP_BUCK = SP_SYSTEM_FALSE;
	int SP_SP_BUCK_INDEX;
	double SP_SP_A;	double SP_SP_RHO; double SP_SP_C;
	for(int i=0;i<sp_sys->number_of_mm_interaction_buck;i++)
	{
		if( strcmp( sp_sys->interaction_type_buck[i][0], "SP" ) == 0 && strcmp( sp_sys->interaction_type_buck[i][2], "SP" ) == 0 )
		{	SP_SP_BUCK_INDEX = i;
			SP_SP_BUCK = SP_SYSTEM_TRUE;
			
			SP_SP_A   = sp_sys->interaction_ARC_buck[i][0];
			SP_SP_RHO = sp_sys->interaction_ARC_buck[i][1];
			SP_SP_C   = sp_sys->interaction_ARC_buck[i][2];
			
			break;
		}
		//printf("%s\t%s\t%s\t%s\t\t%lf\t%lf\t%lf\n",sp_sys->interaction_type_buck[i][0],sp_sys->interaction_type_buck[i][1],sp_sys->interaction_type_buck[i][2],sp_sys->interaction_type_buck[i][3],
		//sp_sys->interaction_ARC_buck[i][0],sp_sys->interaction_ARC_buck[i][1],sp_sys->interaction_ARC_buck[i][2]);	

	}	// if there is SP-SP Buck then set flag 'True'

    for(int n=0;n<sp_sys->number_of_sp_ion;n++)
    {   for(int m=n+1;m<sp_sys->number_of_sp_ion;m++)   // excluding double cnting & self-interaction   SUM{ N<M }
        {   
            // init pair parameters
            // of 'n'th sp - core
            x1 = gsl_vector_get(sp_sys->sp_ion[n].core_position,0); y1 = gsl_vector_get(sp_sys->sp_ion[n].core_position,1); z1 = gsl_vector_get(sp_sys->sp_ion[n].core_position,2);
            charge1 = sp_sys->sp_ion[n].charge_core;    
            // of 'm'th sp - core
            x2 = gsl_vector_get(sp_sys->sp_ion[m].core_position,0); y2 = gsl_vector_get(sp_sys->sp_ion[m].core_position,1); z2 = gsl_vector_get(sp_sys->sp_ion[m].core_position,2);
            charge2 = sp_sys->sp_ion[m].charge_core;    
            // get dist
            dist = sqrt( pow(x2-x1,2.) + pow( y2-y1,2. ) + pow( z2-z1,2. ) );
            
            ///     ///     ///     ///     ///     ///     ///     ///     ///
            // Get 'n'sp-core vs 'm'sp-core Coulomb Interaction Energy
            
            energy_tmp = EV_UNIT*(charge1*charge2/dist);

	    if( SP_SP_BUCK == SP_SYSTEM_TRUE )
	    {
		energy_tmp += SP_SP_A*pow(M_E,-dist/SP_SP_RHO);		// SHORT_RANGE INTERACTION SP-SP BM
		energy_tmp += SP_SP_C*pow(dist,-6.);			// SHORT_RANGE INTERACTION SP-SP VDW
	    }

            // vs sp_core_classic_energy update
            sp_sys->sp_ion[n].classic_energy_by_sp_core += energy_tmp/2.;
            sp_sys->sp_ion[m].classic_energy_by_sp_core += energy_tmp/2.;

//printf("sp-sp %d - %d e : %lf\n",n,m,energy_tmp);

        }
    }// End of Phase 1: sp_core(n) vs sp_core(m) interaction

    // Phase 2: ion_core(n) vs ion_core(m) interaction
    
	int MM_BUCK = SP_SYSTEM_FALSE;
	int MM_BUCK_INDEX;
	double MM_A;	
	double MM_RHO;	
	double MM_C;

	int n_ion_shell = SP_SYSTEM_FALSE;
	int n_ion_core  = SP_SYSTEM_FALSE;
	int m_ion_shell = SP_SYSTEM_FALSE;
	int m_ion_core  = SP_SYSTEM_FALSE;

    for(int n=0;n<sp_sys->number_of_classic_ion;n++)
    {   for(int m=n+1;m<sp_sys->number_of_classic_ion;m++)   // excluding double cnting & self-interaction   SUM{ N<M }
        {	
	    	if( (sp_sys->classic_ion[n].if_shell == SP_SYSTEM_FALSE && sp_sys->classic_ion[n].if_has_shell == SP_SYSTEM_TRUE ) && ( (m-n) == 1 ) )
			continue;	// if 'n' is core and has shell, then skip 'm' if( m = n+1 ) ... avoiding self interaction

		if( sp_sys->classic_ion[n].if_shell == SP_SYSTEM_TRUE )			// if 'n' th ion is SHELL
			n_ion_shell = SP_SYSTEM_TRUE;
		else if ( sp_sys->classic_ion[n].if_shell == SP_SYSTEM_FALSE )		// if 'n' th ion is CORE	.. it doesnt matter if the core has shell or not
			n_ion_core  = SP_SYSTEM_TRUE;
		
		if( sp_sys->classic_ion[m].if_shell == SP_SYSTEM_TRUE )			// if 'm' th ion is SHELL
			m_ion_shell = SP_SYSTEM_TRUE;
		else if ( sp_sys->classic_ion[m].if_shell == SP_SYSTEM_FALSE )		// if 'm' th ion is CORE
			m_ion_core  = SP_SYSTEM_TRUE;

		// get buck pair potential - check if there is corresponding interaction
int cc;
		for(int i=0;i<sp_sys->number_of_mm_interaction_buck;i++)
		{
			if ( (strcmp( sp_sys->classic_ion[n].atom_name, sp_sys->interaction_type_buck[i][1] ) == 0 && strcmp( sp_sys->classic_ion[m].atom_name, sp_sys->interaction_type_buck[i][3] ) == 0) )	// IF ATOM TYPE IN MATCH 
			{
				if( (strcmp( sp_sys->interaction_type_buck[i][0], "SHELL" ) == 0 && n_ion_shell == SP_SYSTEM_TRUE) && (strcmp( sp_sys->interaction_type_buck[i][2], "SHELL") == 0 && m_ion_shell == SP_SYSTEM_TRUE ) )
				{
					MM_BUCK_INDEX = i;	MM_A = sp_sys->interaction_ARC_buck[i][0];	MM_RHO = sp_sys->interaction_ARC_buck[i][1];	MM_C = sp_sys->interaction_ARC_buck[i][2];
					MM_BUCK = SP_SYSTEM_TRUE;
//cc = 1; printf("cc %d\t",cc);
					break;
				}
				else if( (strcmp( sp_sys->interaction_type_buck[i][0], "SHELL" ) == 0 && n_ion_shell == SP_SYSTEM_TRUE) && (strcmp( sp_sys->interaction_type_buck[i][2],"RIM" ) == 0 && m_ion_core == SP_SYSTEM_TRUE ) )
				{
					MM_BUCK_INDEX = i;	MM_A = sp_sys->interaction_ARC_buck[i][0];	MM_RHO = sp_sys->interaction_ARC_buck[i][1];	MM_C = sp_sys->interaction_ARC_buck[i][2];
					MM_BUCK = SP_SYSTEM_TRUE;
//cc = 2; printf("cc %d\t",cc);
					break;
				}
				else if( (strcmp( sp_sys->interaction_type_buck[i][0], "RIM" ) == 0 && n_ion_core == SP_SYSTEM_TRUE) && (strcmp( sp_sys->interaction_type_buck[i][2],"SHELL" ) == 0 && m_ion_shell == SP_SYSTEM_TRUE ) )
				{
					MM_BUCK_INDEX = i;	MM_A = sp_sys->interaction_ARC_buck[i][0];	MM_RHO = sp_sys->interaction_ARC_buck[i][1];	MM_C = sp_sys->interaction_ARC_buck[i][2];
					MM_BUCK = SP_SYSTEM_TRUE;
//cc = 3; printf("cc %d\t",cc);
					break;
				}
				else if( (strcmp( sp_sys->interaction_type_buck[i][0], "RIM" ) == 0 && n_ion_core == SP_SYSTEM_TRUE) && (strcmp( sp_sys->interaction_type_buck[i][2],"RIM" ) == 0 && m_ion_core == SP_SYSTEM_TRUE ) )
				{
					MM_BUCK_INDEX = i;	MM_A = sp_sys->interaction_ARC_buck[i][0];	MM_RHO = sp_sys->interaction_ARC_buck[i][1];	MM_C = sp_sys->interaction_ARC_buck[i][2];
					MM_BUCK = SP_SYSTEM_TRUE;
//cc = 4; printf("cc %d\t",cc);
					break;
				}
			}
			else if ( (strcmp( sp_sys->classic_ion[n].atom_name, sp_sys->interaction_type_buck[i][3] ) == 0 && strcmp( sp_sys->classic_ion[m].atom_name, sp_sys->interaction_type_buck[i][1] ) == 0) )	// IF ATOM TYPE IN MATCH 
			{
				if( (strcmp( sp_sys->interaction_type_buck[i][2], "SHELL" ) == 0 && n_ion_shell == SP_SYSTEM_TRUE) && (strcmp( sp_sys->interaction_type_buck[i][0], "SHELL") == 0 && m_ion_shell == SP_SYSTEM_TRUE ) )
				{
					MM_BUCK_INDEX = i;	MM_A = sp_sys->interaction_ARC_buck[i][0];	MM_RHO = sp_sys->interaction_ARC_buck[i][1];	MM_C = sp_sys->interaction_ARC_buck[i][2];
					MM_BUCK = SP_SYSTEM_TRUE;
//cc = 5; printf("cc %d\t",cc);
					break;
				}
				else if( (strcmp( sp_sys->interaction_type_buck[i][2], "SHELL" ) == 0 && n_ion_shell == SP_SYSTEM_TRUE) && (strcmp( sp_sys->interaction_type_buck[i][0],"RIM" ) == 0 && m_ion_core == SP_SYSTEM_TRUE ) )
				{
					MM_BUCK_INDEX = i;	MM_A = sp_sys->interaction_ARC_buck[i][0];	MM_RHO = sp_sys->interaction_ARC_buck[i][1];	MM_C = sp_sys->interaction_ARC_buck[i][2];
					MM_BUCK = SP_SYSTEM_TRUE;
//cc = 6; printf("cc %d\t",cc);
					break;
				}
				else if( (strcmp( sp_sys->interaction_type_buck[i][2], "RIM" ) == 0 && n_ion_core == SP_SYSTEM_TRUE) && (strcmp( sp_sys->interaction_type_buck[i][0],"SHELL" ) == 0 && m_ion_shell == SP_SYSTEM_TRUE ) )
				{
					MM_BUCK_INDEX = i;	MM_A = sp_sys->interaction_ARC_buck[i][0];	MM_RHO = sp_sys->interaction_ARC_buck[i][1];	MM_C = sp_sys->interaction_ARC_buck[i][2];
					MM_BUCK = SP_SYSTEM_TRUE;
//cc = 7; printf("cc %d\t",cc);
					break;
				}
				else if( (strcmp( sp_sys->interaction_type_buck[i][2], "RIM" ) == 0 && n_ion_core == SP_SYSTEM_TRUE) && (strcmp( sp_sys->interaction_type_buck[i][0],"RIM" ) == 0 && m_ion_core == SP_SYSTEM_TRUE ) )
				{
					MM_BUCK_INDEX = i;	MM_A = sp_sys->interaction_ARC_buck[i][0];	MM_RHO = sp_sys->interaction_ARC_buck[i][1];	MM_C = sp_sys->interaction_ARC_buck[i][2];
					MM_BUCK = SP_SYSTEM_TRUE;
//cc = 8; printf("cc %d\t",cc);
					break;
				}
			}
   		}	
		n_ion_shell = SP_SYSTEM_FALSE;	n_ion_core = SP_SYSTEM_FALSE;	m_ion_shell = SP_SYSTEM_FALSE;	m_ion_core = SP_SYSTEM_FALSE;	// reset flags

		// SET PAIR POTENTIAL DONE (PP)

	    // init pair parameters
	    // of 'n'th species
	    x1 = gsl_vector_get(sp_sys->classic_ion[n].core_position,0); y1 = gsl_vector_get(sp_sys->classic_ion[n].core_position,1); z1 = gsl_vector_get(sp_sys->classic_ion[n].core_position,2);
	    charge1 = sp_sys->classic_ion[n].charge_core;
	    // of 'm'th species
	    x2 = gsl_vector_get(sp_sys->classic_ion[m].core_position,0); y2 = gsl_vector_get(sp_sys->classic_ion[m].core_position,1); z2 = gsl_vector_get(sp_sys->classic_ion[m].core_position,2);
	    charge2 = sp_sys->classic_ion[m].charge_core;
	    // get dist
	    dist = sqrt( pow(x2-x1,2.) + pow( y2-y1,2. ) + pow( z2-z1,2. ) );

	    // reset 'energy_tmp'
	    // energy_tmp = 0.;
		energy_tmp = EV_UNIT*(charge1*charge2/dist);			// COULOMB INTERACTION
	
		if ( MM_BUCK == SP_SYSTEM_TRUE )
		{	
		//printf("flag\n");
		//printf("%lf\t%lf\t%lf\n",MM_A,MM_RHO,MM_C);
		//printf("%d\t%d\n",sp_sys->classic_ion[n].if_shell,sp_sys->classic_ion[m].if_shell);
			energy_tmp += MM_A*pow(M_E,-dist/MM_RHO);
			energy_tmp += MM_C*pow(dist,-6.); 
		//printf("%lf\t%d\n",energy_tmp,cnt);
		}								// SHORT_RANGE INTERACTION ... IF EXISTS
		MM_BUCK = SP_SYSTEM_FALSE;
//printf("mm-mm %d - %d e : %lf\n",n,m,energy_tmp);
	    // vs sp_core_classic_energy update
            sp_sys->classic_ion[n].classic_energy_by_ion_core += energy_tmp/2.;
	    sp_sys->classic_ion[m].classic_energy_by_ion_core += energy_tmp/2.;
        }
    }// End of Phase 2: ion_core(n) vs ion_core(m) interaction

    // Phase 2 - 1 : shell model elastic energy
    for(int n=0;n<sp_sys->number_of_classic_ion;n++)
    {
	if( sp_sys->classic_ion[n].if_shell == SP_SYSTEM_FALSE && sp_sys->classic_ion[n].if_has_shell == SP_SYSTEM_TRUE )	// it is core and has shell ... then
	{
		x1 = gsl_vector_get(sp_sys->classic_ion[n].core_position,0); y1 = gsl_vector_get(sp_sys->classic_ion[n].core_position,1); 
		z1 = gsl_vector_get(sp_sys->classic_ion[n].core_position,2);							// core position
		x2 = gsl_vector_get(sp_sys->classic_ion[n+1].core_position,0); y2 = gsl_vector_get(sp_sys->classic_ion[n+1].core_position,1); 
		z2 = gsl_vector_get(sp_sys->classic_ion[n+1].core_position,2);							// shell position
		dist = sqrt( pow(x2-x1,2.) + pow( y2-y1,2. ) + pow( z2-z1,2. ) );
		// Elastic energy
		//energy_tmp = 1./2.*sp_sys->classic_ion[n].k2_const*pow(dist,2.);

		// k2 + k4 elastic effect
		energy_tmp = 1./2.*sp_sys->classic_ion[n].k2_const*pow(dist,2.) + 1./24.*sp_sys->classic_ion[n].k4_const*pow(dist,4.);

		// update elastic energy
		sp_sys->classic_ion[n].classic_energy_by_ion_core += energy_tmp/2.;
		sp_sys->classic_ion[n+1].classic_energy_by_ion_core += energy_tmp/2.;
	}
    }


    ///     ///     ///     ///     ///     ///     ///     ///     ///

    // Phase 3: sp_core(n) vs ion_core(m) interaction

    int SP_MM_BUCK = SP_SYSTEM_FALSE;
    int SP_MM_BUCK_INDEX;	
    // int m_ion_core/shell -> reuse
    double SP_MM_A;	
    double SP_MM_RHO;	
    double SP_MM_C;


    for(int n=0;n<sp_sys->number_of_sp_ion;n++)
    {   for(int m=0;m<sp_sys->number_of_classic_ion;m++)
        {
		if( sp_sys->classic_ion[m].if_shell == SP_SYSTEM_TRUE )
			m_ion_shell = SP_SYSTEM_TRUE;
		else if( sp_sys->classic_ion[m].if_shell == SP_SYSTEM_FALSE )
			m_ion_core  = SP_SYSTEM_TRUE;

		// get buck pair potential - check if there is corresponding interaction
		for(int i=0;i<sp_sys->number_of_mm_interaction_buck;i++)
		{
			if( strcmp(sp_sys->sp_ion[n].atom_name,sp_sys->interaction_type_buck[i][1]) == 0 && strcmp(sp_sys->classic_ion[m].atom_name,sp_sys->interaction_type_buck[i][3]) == 0 )
			{
				if( (strcmp(sp_sys->interaction_type_buck[i][0],"SP") == 0) && ( strcmp(sp_sys->interaction_type_buck[i][2],"SHELL") == 0 && m_ion_shell == SP_SYSTEM_TRUE ) )
				{
					SP_MM_BUCK_INDEX = i;	SP_MM_A = sp_sys->interaction_ARC_buck[i][0];	SP_MM_RHO = sp_sys->interaction_ARC_buck[i][1];	SP_MM_C = sp_sys->interaction_ARC_buck[i][2]; 
					SP_MM_BUCK = SP_SYSTEM_TRUE;
					break;
				} 
				else if( (strcmp(sp_sys->interaction_type_buck[i][0],"SP") == 0) && ( strcmp(sp_sys->interaction_type_buck[i][2],"RIM") == 0 && m_ion_core == SP_SYSTEM_TRUE ) )
				{
					SP_MM_BUCK_INDEX = i;	SP_MM_A = sp_sys->interaction_ARC_buck[i][0];	SP_MM_RHO = sp_sys->interaction_ARC_buck[i][1];	SP_MM_C = sp_sys->interaction_ARC_buck[i][2]; 
					SP_MM_BUCK = SP_SYSTEM_TRUE;
					break;
				}
			}
			else if( strcmp(sp_sys->sp_ion[n].atom_name,sp_sys->interaction_type_buck[i][3]) == 0 && strcmp(sp_sys->classic_ion[m].atom_name,sp_sys->interaction_type_buck[i][1]) == 0 )
			{
				if( (strcmp(sp_sys->interaction_type_buck[i][2],"SP") == 0) && ( strcmp(sp_sys->interaction_type_buck[i][0],"SHELL") == 0 && m_ion_shell == SP_SYSTEM_TRUE ) )
				{
					SP_MM_BUCK_INDEX = i;	SP_MM_A = sp_sys->interaction_ARC_buck[i][0];	SP_MM_RHO = sp_sys->interaction_ARC_buck[i][1];	SP_MM_C = sp_sys->interaction_ARC_buck[i][2]; 
					SP_MM_BUCK = SP_SYSTEM_TRUE;
					break;
				} 
				else if( (strcmp(sp_sys->interaction_type_buck[i][2],"SP") == 0) && ( strcmp(sp_sys->interaction_type_buck[i][0],"RIM") == 0 && m_ion_core == SP_SYSTEM_TRUE ) )
				{
					SP_MM_BUCK_INDEX = i;	SP_MM_A = sp_sys->interaction_ARC_buck[i][0];	SP_MM_RHO = sp_sys->interaction_ARC_buck[i][1];	SP_MM_C = sp_sys->interaction_ARC_buck[i][2]; 
					SP_MM_BUCK = SP_SYSTEM_TRUE;
					break;
				}
			}
		}
		m_ion_shell = SP_SYSTEM_FALSE;	m_ion_core = SP_SYSTEM_FALSE;

            // init pair parameters
            
            // of 'n'th sp - core
            x1 = gsl_vector_get(sp_sys->sp_ion[n].core_position,0); y1 = gsl_vector_get(sp_sys->sp_ion[n].core_position,1); z1 = gsl_vector_get(sp_sys->sp_ion[n].core_position,2);
            charge1 = sp_sys->sp_ion[n].charge_core;    vdw1 = sp_sys->sp_ion[n].cent_short_range_vdw;
            // of 'm'th ion - core
            x2 = gsl_vector_get(sp_sys->classic_ion[m].core_position,0); y2 = gsl_vector_get(sp_sys->classic_ion[m].core_position,1); z2 = gsl_vector_get(sp_sys->classic_ion[m].core_position,2);
            charge2 = sp_sys->classic_ion[m].charge_core;   vdw2 = sp_sys->classic_ion[m].short_range_c;
            // get dist
            dist = sqrt( pow(x2-x1,2.) + pow( y2-y1,2. ) + pow( z2-z1,2. ) );

	    // reset 'energy_tmp'
	    energy_tmp = 0.;
	    energy_tmp = EV_UNIT*(charge1*charge2/dist);
	    
	    if( SP_MM_BUCK == SP_SYSTEM_TRUE )
	    {
		energy_tmp += SP_MM_A*pow(M_E,-dist/SP_MM_RHO);
		energy_tmp += SP_MM_C*pow(dist,-6.);
	    }
	    SP_MM_BUCK = SP_SYSTEM_FALSE;
//printf("sp[%d] - mm[%d] e : %lf\n",n+1,m+1,energy_tmp);
            // vs_core_classic_energy_update
            sp_sys->sp_ion[n].classic_energy_by_ion_core += energy_tmp/2.;
            sp_sys->classic_ion[n].classic_energy_by_sp_core += energy_tmp/2.;
        }
    }// END OF Phase 3;

    ///     ///     ///     ///     ///     ///     ///     ///     ///

    // README
    /*
     * In the update 07252019 (Thur)
     *
     * sp(n)-sps  core interaction classic energy is in 'sp_sys->sp_ion[n].classic_energy_by_sp_core'
     *
     * sp(n)-ions core interaction classic energy is in 'sp_sys->sp_ion[n].classic_energy_by_ion_core'
     *
     *
     * ion(n)-sps  core interaction classic energy is in 'sp-sys->classic_ion[n].classic_energy_by_sp_core'
     *
     * ion(n)-ions core interaction classic energy is in 'sp_sys->classic_ion[n].classic energy_by_ion_core'
     *
     * None of changes made yet for further treatments ...
     *
     */


    return;
}






// Tag:classic_force
void sp_cluster_system_get_classic_force( sp_cluster_system* sp_sys )
{
    // 25072019 update

    // COMMON VARIABLES
    double xn,yn,zn, xm,ym,zm;  // coordinate tmp
    double short_range_an, short_range_am, short_range_rn, short_range_rm; // BM-like potential short_range A rho param tmp
    double vdwn, vdwm;  // van der waals param for (1/r^6) term among classical ions
    double chargen, chargem;    // charge tmp
    double dist; // interionic distance tmp
    double deriv_tmp[3] = {0.,0.,0.};

    /* force equation is given by ( Direction Convention )
     *
     *  when foo -> f( r_m - r_n ) ... i.e., r_n is pointing (or watching) r_m
     *
     *       d                          q_m*q_n
     *  ------------  ------------------------------------------------------------
     *     d x_n      [(x_m - x_n)^2 + (y_m - y_n)^2 +(z_m - z_n)^2]^1/2
     *
     *
     *                    q_m*q_n*(x_m - x_n)
     *     =    ------------------------------------------------------------
     *           [(x_m - x_n)^2 + (y_m - y_n)^2 +(z_m - z_m)^2]^3/2
     *
     *
     * i.e., it calculates derivative w.r.t atom 'n', and its inversely signed value will tell the derivative
     * w.r.t atom 'm'. if one multiply '-1' at each, then it will tell the direction where the force acts on.
     *
     *       d                           
     *  ------------  Amn Exp[ - [(x_m - x_n)^2 + (y_m - y_n)^2 +(z_m - z_n)^2]^1/2 / rho_mn ] - Cmn/([(x_m - x_n)^2 + (y_m - y_n)^2 +(z_m - z_n)^2]^1/2)^6
     *     d x_n     
     *
     *
     *
     *      =   (xm - xn)( 6*Cmn/[(xm - xn)^2 + (ym - yn)^2 +(zm - zn)^2]^4 + Amn Exp[ - [(xm - xn)^2 + (ym - yn)^2 +(zm - zn)^2]^1/2 / rhomn ]
     *                                                                        -----------------------------------------------------------------
     *                                                                              [(xm - xn)^2 + (ym - yn)^2 +(zm - zn)^2]^1/2 * rhomn
     *
     */

    /*
     *  Prep    - Init all configs ( set into zero )
     *
     *  Phase 1 - sp_ion_core vs sp_ion_core
     *
     *  Phase 2 - classic_ion_core vs classic_ion_core
     *
     *  Phase 3 - sp_ion_core vs classic_ion_core
     */
    
    // Prep ... Refreshing memory
    for(int n=0;n<sp_sys->number_of_sp_ion;n++)
    {   memset(sp_sys->sp_ion[n].force_by_sp_core,0.,3*sizeof(double));
        memset(sp_sys->sp_ion[n].force_by_ion_core,0.,3*sizeof(double));    }
    for(int n=0;n<sp_sys->number_of_classic_ion;n++)
    {   memset(sp_sys->classic_ion[n].force_by_sp_core,0.,3*sizeof(double));
        memset(sp_sys->classic_ion[n].force_by_ion_core,0.,3*sizeof(double));   }

    // Phase 1 : sp_core(n) vs sp_core(m) force estimation
    
	int SP_SP_BUCK = SP_SYSTEM_FALSE;
	int SP_SP_BUCK_INDEX;
	double SP_SP_A;	double SP_SP_RHO; double SP_SP_C;

	for(int i=0;i<sp_sys->number_of_mm_interaction_buck;i++)
	{
		if( strcmp( sp_sys->interaction_type_buck[i][0], "SP" ) == 0 && strcmp( sp_sys->interaction_type_buck[i][2], "SP" ) == 0 )
		{	SP_SP_BUCK_INDEX = i;
			SP_SP_BUCK = SP_SYSTEM_TRUE;
			
			SP_SP_A   = sp_sys->interaction_ARC_buck[i][0];
			SP_SP_RHO = sp_sys->interaction_ARC_buck[i][1];
			SP_SP_C   = sp_sys->interaction_ARC_buck[i][2];
			break;
		}
	}	// if there is SP-SP Buck then set flag 'True'

    for(int n=0;n<sp_sys->number_of_sp_ion;n++)
    {   for(int m=n+1;m<sp_sys->number_of_sp_ion;m++)
        {
            // init pair parameters (m,n) == (2,1) index rule
            // of 'n'th sp - core
            xn = gsl_vector_get(sp_sys->sp_ion[n].core_position,0); yn = gsl_vector_get(sp_sys->sp_ion[n].core_position,1); zn = gsl_vector_get(sp_sys->sp_ion[n].core_position,2);
            chargen = sp_sys->sp_ion[n].charge_core;    
            // of 'm'th sp - core
            xm = gsl_vector_get(sp_sys->sp_ion[m].core_position,0); ym = gsl_vector_get(sp_sys->sp_ion[m].core_position,1); zm = gsl_vector_get(sp_sys->sp_ion[m].core_position,2);
            chargem = sp_sys->sp_ion[m].charge_core;    
            // get dist
            dist = sqrt( pow(xm-xn,2.) + pow( ym-yn,2. ) + pow( zm-zn,2. ) );   // 'r_mn'
            
            ///     ///     ///     ///     ///     ///     ///     ///     ///
            // Get derivative of 'n'sp-core w.r.t 'm'sp-core Coulomb Force 
            deriv_tmp[0] = EV_UNIT*(chargen*chargem*(xm-xn)/pow(dist,3.));       
            deriv_tmp[1] = EV_UNIT*(chargen*chargem*(ym-yn)/pow(dist,3.));       
            deriv_tmp[2] = EV_UNIT*(chargen*chargem*(zm-zn)/pow(dist,3.));       
            
		if( SP_SP_BUCK == SP_SYSTEM_TRUE )
		{
			deriv_tmp[0] += SP_SP_A/SP_SP_RHO*pow(M_E,-dist/SP_SP_RHO)/dist*(xm-xn);
			deriv_tmp[1] += SP_SP_A/SP_SP_RHO*pow(M_E,-dist/SP_SP_RHO)/dist*(ym-yn);
			deriv_tmp[2] += SP_SP_A/SP_SP_RHO*pow(M_E,-dist/SP_SP_RHO)/dist*(zm-zn);			// BM FORCE
		
			deriv_tmp[0] += SP_SP_C*6./pow(dist,8.)*(xm-xn);
			deriv_tmp[1] += SP_SP_C*6./pow(dist,8.)*(ym-yn);
			deriv_tmp[2] += SP_SP_C*6./pow(dist,8.)*(zm-zn);					// VDW FORCE
		}

            // force update ... again note that the above calculation estimates derivative of 'n'th sp-core w.r.t 'm'th sp-core
            sp_sys->sp_ion[n].force_by_sp_core[0] -= deriv_tmp[0];  sp_sys->sp_ion[n].force_by_sp_core[1] -= deriv_tmp[1];  sp_sys->sp_ion[n].force_by_sp_core[2] -= deriv_tmp[2];
            sp_sys->sp_ion[m].force_by_sp_core[0] += deriv_tmp[0];  sp_sys->sp_ion[m].force_by_sp_core[1] += deriv_tmp[1];  sp_sys->sp_ion[m].force_by_sp_core[2] += deriv_tmp[2];
        }
    }// End of Phase 1

    // Phase 2 : ion_core(n) vs ion_core(m) force estimation
    
	int MM_BUCK = SP_SYSTEM_FALSE;
	int MM_BUCK_INDEX;
	double MM_A;	
	double MM_RHO;	
	double MM_C;

	int n_ion_shell = SP_SYSTEM_FALSE;
	int n_ion_core  = SP_SYSTEM_FALSE;
	int m_ion_shell = SP_SYSTEM_FALSE;
	int m_ion_core  = SP_SYSTEM_FALSE;	// identifier

    for(int n=0;n<sp_sys->number_of_classic_ion;n++)
    {   for(int m=n+1;m<sp_sys->number_of_classic_ion;m++)
        {
	    	if( (sp_sys->classic_ion[n].if_shell == SP_SYSTEM_FALSE && sp_sys->classic_ion[n].if_has_shell == SP_SYSTEM_TRUE ) && ( (m-n) == 1 ) )
			continue;	// if 'n' is core and has shell, then skip 'm' if( m = n+1 ) ... avoiding self interaction

		if( sp_sys->classic_ion[n].if_shell == SP_SYSTEM_TRUE )			// if 'n' th ion is SHELL
			n_ion_shell = SP_SYSTEM_TRUE;
		else if ( sp_sys->classic_ion[n].if_shell == SP_SYSTEM_FALSE )		// if 'n' th ion is CORE	.. it doesnt matter if the core has shell or not
			n_ion_core  = SP_SYSTEM_TRUE;
		
		if( sp_sys->classic_ion[m].if_shell == SP_SYSTEM_TRUE )			// if 'm' th ion is SHELL
			m_ion_shell = SP_SYSTEM_TRUE;
		else if ( sp_sys->classic_ion[m].if_shell == SP_SYSTEM_FALSE )		// if 'm' th ion is CORE
			m_ion_core  = SP_SYSTEM_TRUE;

		// get buck pair potential - check if there is corresponding interaction
		for(int i=0;i<sp_sys->number_of_mm_interaction_buck;i++)
		{
			if ( (strcmp( sp_sys->classic_ion[n].atom_name, sp_sys->interaction_type_buck[i][1] ) == 0 && strcmp( sp_sys->classic_ion[m].atom_name, sp_sys->interaction_type_buck[i][3] ) == 0) )	// IF ATOM TYPE IN MATCH 
			{
				if( (strcmp( sp_sys->interaction_type_buck[i][0], "SHELL" ) == 0 && n_ion_shell == SP_SYSTEM_TRUE) && (strcmp( sp_sys->interaction_type_buck[i][2], "SHELL") == 0 && m_ion_shell == SP_SYSTEM_TRUE ) )
				{
					MM_BUCK_INDEX = i;	MM_A = sp_sys->interaction_ARC_buck[i][0];	MM_RHO = sp_sys->interaction_ARC_buck[i][1];	MM_C = sp_sys->interaction_ARC_buck[i][2];
					MM_BUCK = SP_SYSTEM_TRUE;
					break;
				}
				else if( (strcmp( sp_sys->interaction_type_buck[i][0], "SHELL" ) == 0 && n_ion_shell == SP_SYSTEM_TRUE) && (strcmp( sp_sys->interaction_type_buck[i][2],"RIM" ) == 0 && m_ion_core == SP_SYSTEM_TRUE ) )
				{
					MM_BUCK_INDEX = i;	MM_A = sp_sys->interaction_ARC_buck[i][0];	MM_RHO = sp_sys->interaction_ARC_buck[i][1];	MM_C = sp_sys->interaction_ARC_buck[i][2];
					MM_BUCK = SP_SYSTEM_TRUE;
					break;
				}
				else if( (strcmp( sp_sys->interaction_type_buck[i][0], "RIM" ) == 0 && n_ion_core == SP_SYSTEM_TRUE) && (strcmp( sp_sys->interaction_type_buck[i][2],"SHELL" ) == 0 && m_ion_shell == SP_SYSTEM_TRUE ) )
				{
					MM_BUCK_INDEX = i;	MM_A = sp_sys->interaction_ARC_buck[i][0];	MM_RHO = sp_sys->interaction_ARC_buck[i][1];	MM_C = sp_sys->interaction_ARC_buck[i][2];
					MM_BUCK = SP_SYSTEM_TRUE;
					break;
				}
				else if( (strcmp( sp_sys->interaction_type_buck[i][0], "RIM" ) == 0 && n_ion_core == SP_SYSTEM_TRUE) && (strcmp( sp_sys->interaction_type_buck[i][2],"RIM" ) == 0 && m_ion_core == SP_SYSTEM_TRUE ) )
				{
					MM_BUCK_INDEX = i;	MM_A = sp_sys->interaction_ARC_buck[i][0];	MM_RHO = sp_sys->interaction_ARC_buck[i][1];	MM_C = sp_sys->interaction_ARC_buck[i][2];
					MM_BUCK = SP_SYSTEM_TRUE;
					break;
				}
			}
			else if ( (strcmp( sp_sys->classic_ion[n].atom_name, sp_sys->interaction_type_buck[i][3] ) == 0 && strcmp( sp_sys->classic_ion[m].atom_name, sp_sys->interaction_type_buck[i][1] ) == 0) )	// IF ATOM TYPE IN MATCH 
			{
				if( (strcmp( sp_sys->interaction_type_buck[i][2], "SHELL" ) == 0 && n_ion_shell == SP_SYSTEM_TRUE) && (strcmp( sp_sys->interaction_type_buck[i][0], "SHELL") == 0 && m_ion_shell == SP_SYSTEM_TRUE ) )
				{
					MM_BUCK_INDEX = i;	MM_A = sp_sys->interaction_ARC_buck[i][0];	MM_RHO = sp_sys->interaction_ARC_buck[i][1];	MM_C = sp_sys->interaction_ARC_buck[i][2];
					MM_BUCK = SP_SYSTEM_TRUE;
					break;
				}
				else if( (strcmp( sp_sys->interaction_type_buck[i][2], "SHELL" ) == 0 && n_ion_shell == SP_SYSTEM_TRUE) && (strcmp( sp_sys->interaction_type_buck[i][0],"RIM" ) == 0 && m_ion_core == SP_SYSTEM_TRUE ) )
				{
					MM_BUCK_INDEX = i;	MM_A = sp_sys->interaction_ARC_buck[i][0];	MM_RHO = sp_sys->interaction_ARC_buck[i][1];	MM_C = sp_sys->interaction_ARC_buck[i][2];
					MM_BUCK = SP_SYSTEM_TRUE;
					break;
				}
				else if( (strcmp( sp_sys->interaction_type_buck[i][2], "RIM" ) == 0 && n_ion_core == SP_SYSTEM_TRUE) && (strcmp( sp_sys->interaction_type_buck[i][0],"SHELL" ) == 0 && m_ion_shell == SP_SYSTEM_TRUE ) )
				{
					MM_BUCK_INDEX = i;	MM_A = sp_sys->interaction_ARC_buck[i][0];	MM_RHO = sp_sys->interaction_ARC_buck[i][1];	MM_C = sp_sys->interaction_ARC_buck[i][2];
					MM_BUCK = SP_SYSTEM_TRUE;
					break;
				}
				else if( (strcmp( sp_sys->interaction_type_buck[i][2], "RIM" ) == 0 && n_ion_core == SP_SYSTEM_TRUE) && (strcmp( sp_sys->interaction_type_buck[i][0],"RIM" ) == 0 && m_ion_core == SP_SYSTEM_TRUE ) )
				{
					MM_BUCK_INDEX = i;	MM_A = sp_sys->interaction_ARC_buck[i][0];	MM_RHO = sp_sys->interaction_ARC_buck[i][1];	MM_C = sp_sys->interaction_ARC_buck[i][2];
					MM_BUCK = SP_SYSTEM_TRUE;
					break;
				}
			}
   		}	
		n_ion_shell = SP_SYSTEM_FALSE;	n_ion_core = SP_SYSTEM_FALSE;	m_ion_shell = SP_SYSTEM_FALSE;	m_ion_core = SP_SYSTEM_FALSE;	// reset flags

		// SET PAIR POTENTIAL DONE (PP)

            // init params
            //
            // of 'n'th classic ion
            xn = gsl_vector_get(sp_sys->classic_ion[n].core_position,0); yn = gsl_vector_get(sp_sys->classic_ion[n].core_position,1); zn = gsl_vector_get(sp_sys->classic_ion[n].core_position,2);
            short_range_an = sp_sys->classic_ion[n].short_range_a;  short_range_rn = sp_sys->classic_ion[n].short_range_r;  vdwn = sp_sys->classic_ion[n].short_range_c;
            chargen = sp_sys->classic_ion[n].charge_core;
            // of 'm'th classic ion
            xm = gsl_vector_get(sp_sys->classic_ion[m].core_position,0); ym = gsl_vector_get(sp_sys->classic_ion[m].core_position,1); zm = gsl_vector_get(sp_sys->classic_ion[m].core_position,2);
            short_range_am = sp_sys->classic_ion[m].short_range_a;  short_range_rm = sp_sys->classic_ion[m].short_range_r;  vdwm = sp_sys->classic_ion[m].short_range_c;
            chargem = sp_sys->classic_ion[m].charge_core;
            // get dist
            dist = sqrt( pow(xm-xn,2.) + pow( ym-yn,2. ) + pow( zm-zn,2. ) );   // 'r_mn'

		deriv_tmp[0] = EV_UNIT*(chargen*chargem)*(xm-xn)/pow(dist,3.);
		deriv_tmp[1] = EV_UNIT*(chargen*chargem)*(ym-yn)/pow(dist,3.);
		deriv_tmp[2] = EV_UNIT*(chargen*chargem)*(zm-zn)/pow(dist,3.);

		if( MM_BUCK == SP_SYSTEM_TRUE )
		{
			deriv_tmp[0] += MM_A/MM_RHO*pow(M_E,-dist/MM_RHO)/dist*(xm-xn);
			deriv_tmp[1] += MM_A/MM_RHO*pow(M_E,-dist/MM_RHO)/dist*(ym-yn);
			deriv_tmp[2] += MM_A/MM_RHO*pow(M_E,-dist/MM_RHO)/dist*(zm-zn);

			deriv_tmp[0] += MM_C*6./pow(dist,8.)*(xm-xn);
			deriv_tmp[1] += MM_C*6./pow(dist,8.)*(ym-yn);
			deriv_tmp[2] += MM_C*6./pow(dist,8.)*(zm-zn);
		}
		MM_BUCK = SP_SYSTEM_FALSE;

            // force update ... again note that the above calculation estimates derivative of 'n'th sp-core w.r.t 'm'th sp-core
            sp_sys->classic_ion[n].force_by_ion_core[0] -= deriv_tmp[0];  sp_sys->classic_ion[n].force_by_ion_core[1] -= deriv_tmp[1];  sp_sys->classic_ion[n].force_by_ion_core[2] -= deriv_tmp[2];
            sp_sys->classic_ion[m].force_by_ion_core[0] += deriv_tmp[0];  sp_sys->classic_ion[m].force_by_ion_core[1] += deriv_tmp[1];  sp_sys->classic_ion[m].force_by_ion_core[2] += deriv_tmp[2];
        }
    }// End of Phase 2

	// Phase 2 - 1 : shell model elastic force
	for(int n=0;n<sp_sys->number_of_classic_ion;n++)
	{
		if( sp_sys->classic_ion[n].if_shell == SP_SYSTEM_FALSE && sp_sys->classic_ion[n].if_has_shell == SP_SYSTEM_TRUE )	// it is core and has shell ... then
		{
			xn = gsl_vector_get(sp_sys->classic_ion[n].core_position,0); yn = gsl_vector_get(sp_sys->classic_ion[n].core_position,1); 
			zn = gsl_vector_get(sp_sys->classic_ion[n].core_position,2);							// core position
			xm = gsl_vector_get(sp_sys->classic_ion[n+1].core_position,0); ym = gsl_vector_get(sp_sys->classic_ion[n+1].core_position,1); 
			zm = gsl_vector_get(sp_sys->classic_ion[n+1].core_position,2);							// shell position
			// elastic force ... (derivative w.r.t. r1 or core)
			/*
			deriv_tmp[0] = sp_sys->classic_ion[n].k2_const*(xn-xm);
			deriv_tmp[1] = sp_sys->classic_ion[n].k2_const*(yn-ym);
			deriv_tmp[2] = sp_sys->classic_ion[n].k2_const*(zn-zm);
			*/
			// k4 elastic const
			deriv_tmp[0] = sp_sys->classic_ion[n].k2_const*(xn-xm) + 1./6.*sp_sys->classic_ion[n].k4_const*(xn-xm)*(pow(xn-xm,2.)+pow(yn-ym,2.)+pow(zn-zm,2.));
			deriv_tmp[1] = sp_sys->classic_ion[n].k2_const*(yn-ym) + 1./6.*sp_sys->classic_ion[n].k4_const*(yn-ym)*(pow(xn-xm,2.)+pow(yn-ym,2.)+pow(zn-zm,2.));
			deriv_tmp[2] = sp_sys->classic_ion[n].k2_const*(zn-zm) + 1./6.*sp_sys->classic_ion[n].k4_const*(zn-zm)*(pow(xn-xm,2.)+pow(yn-ym,2.)+pow(zn-zm,2.));


			// update elastic force
			sp_sys->classic_ion[n].force_by_ion_core[0] -= deriv_tmp[0];  sp_sys->classic_ion[n].force_by_ion_core[1] -= deriv_tmp[1];  sp_sys->classic_ion[n].force_by_ion_core[2] -= deriv_tmp[2];
			sp_sys->classic_ion[n+1].force_by_ion_core[0] += deriv_tmp[0];  sp_sys->classic_ion[n+1].force_by_ion_core[1] += deriv_tmp[1];  sp_sys->classic_ion[n+1].force_by_ion_core[2] += deriv_tmp[2];
		}
	}

    // Phase 3 : sp-core(n) vs ion-core(m)  force calculation
 
    int SP_MM_BUCK = SP_SYSTEM_FALSE;
    int SP_MM_BUCK_INDEX;	
    // int m_ion_core/shell -> reuse
    double SP_MM_A;	
    double SP_MM_RHO;	
    double SP_MM_C;
   
    for(int n=0;n<sp_sys->number_of_sp_ion;n++)
    {   for(int m=0;m<sp_sys->number_of_classic_ion;m++)
        {
		if( sp_sys->classic_ion[m].if_shell == SP_SYSTEM_TRUE )
			m_ion_shell = SP_SYSTEM_TRUE;
		else if( sp_sys->classic_ion[m].if_shell == SP_SYSTEM_FALSE )
			m_ion_core  = SP_SYSTEM_TRUE;

		// get buck pair potential - check if there is corresponding interaction
		for(int i=0;i<sp_sys->number_of_mm_interaction_buck;i++)
		{
			if( strcmp(sp_sys->sp_ion[n].atom_name,sp_sys->interaction_type_buck[i][1]) == 0 && strcmp(sp_sys->classic_ion[m].atom_name,sp_sys->interaction_type_buck[i][3]) == 0 )
			{
				if( (strcmp(sp_sys->interaction_type_buck[i][0],"SP") == 0) && ( strcmp(sp_sys->interaction_type_buck[i][2],"SHELL") == 0 && m_ion_shell == SP_SYSTEM_TRUE ) )
				{
					SP_MM_BUCK_INDEX = i;	SP_MM_A = sp_sys->interaction_ARC_buck[i][0];	SP_MM_RHO = sp_sys->interaction_ARC_buck[i][1];	SP_MM_C = sp_sys->interaction_ARC_buck[i][2]; 
					SP_MM_BUCK = SP_SYSTEM_TRUE;
					break;
				} 
				else if( (strcmp(sp_sys->interaction_type_buck[i][0],"SP") == 0) && ( strcmp(sp_sys->interaction_type_buck[i][2],"RIM") == 0 && m_ion_core == SP_SYSTEM_TRUE ) )
				{
					SP_MM_BUCK_INDEX = i;	SP_MM_A = sp_sys->interaction_ARC_buck[i][0];	SP_MM_RHO = sp_sys->interaction_ARC_buck[i][1];	SP_MM_C = sp_sys->interaction_ARC_buck[i][2]; 
					SP_MM_BUCK = SP_SYSTEM_TRUE;
					break;
				}
			}
			else if( strcmp(sp_sys->sp_ion[n].atom_name,sp_sys->interaction_type_buck[i][3]) == 0 && strcmp(sp_sys->classic_ion[m].atom_name,sp_sys->interaction_type_buck[i][1]) == 0 )
			{
				if( (strcmp(sp_sys->interaction_type_buck[i][2],"SP") == 0) && ( strcmp(sp_sys->interaction_type_buck[i][0],"SHELL") == 0 && m_ion_shell == SP_SYSTEM_TRUE ) )
				{
					SP_MM_BUCK_INDEX = i;	SP_MM_A = sp_sys->interaction_ARC_buck[i][0];	SP_MM_RHO = sp_sys->interaction_ARC_buck[i][1];	SP_MM_C = sp_sys->interaction_ARC_buck[i][2]; 
					SP_MM_BUCK = SP_SYSTEM_TRUE;
					break;
				} 
				else if( (strcmp(sp_sys->interaction_type_buck[i][2],"SP") == 0) && ( strcmp(sp_sys->interaction_type_buck[i][0],"RIM") == 0 && m_ion_core == SP_SYSTEM_TRUE ) )
				{
					SP_MM_BUCK_INDEX = i;	SP_MM_A = sp_sys->interaction_ARC_buck[i][0];	SP_MM_RHO = sp_sys->interaction_ARC_buck[i][1];	SP_MM_C = sp_sys->interaction_ARC_buck[i][2]; 
					SP_MM_BUCK = SP_SYSTEM_TRUE;
					break;
				}
			}
		}
		m_ion_shell = SP_SYSTEM_FALSE;	m_ion_core = SP_SYSTEM_FALSE;

            // init pair parameters


            // init pair parameters (m,n) == (2,1) index rule
            //
            // of 'n'th sp - core
            xn = gsl_vector_get(sp_sys->sp_ion[n].core_position,0); yn = gsl_vector_get(sp_sys->sp_ion[n].core_position,1); zn = gsl_vector_get(sp_sys->sp_ion[n].core_position,2);
            chargen = sp_sys->sp_ion[n].charge_core;        vdwn = sp_sys->sp_ion[n].cent_short_range_vdw;
            // of 'm'th classic ion
            xm = gsl_vector_get(sp_sys->classic_ion[m].core_position,0); ym = gsl_vector_get(sp_sys->classic_ion[m].core_position,1); zm = gsl_vector_get(sp_sys->classic_ion[m].core_position,2);
            chargem = sp_sys->classic_ion[m].charge_core;   vdwm = sp_sys->classic_ion[m].short_range_c;
            // get dist
            dist = sqrt( pow(xm-xn,2.) + pow( ym-yn,2. ) + pow( zm-zn,2. ) );   // 'r_mn'

		deriv_tmp[0] = EV_UNIT*(chargen*chargem*(xm-xn)/pow(dist,3.));
		deriv_tmp[1] = EV_UNIT*(chargen*chargem*(ym-yn)/pow(dist,3.));
		deriv_tmp[2] = EV_UNIT*(chargen*chargem*(zm-zn)/pow(dist,3.));

		if( SP_MM_BUCK == SP_SYSTEM_TRUE )
		{
			deriv_tmp[0] += SP_MM_A/SP_MM_RHO*pow(M_E,-dist/SP_MM_RHO)/dist*(xm-xn);	
			deriv_tmp[1] += SP_MM_A/SP_MM_RHO*pow(M_E,-dist/SP_MM_RHO)/dist*(ym-yn);	
			deriv_tmp[2] += SP_MM_A/SP_MM_RHO*pow(M_E,-dist/SP_MM_RHO)/dist*(zm-zn);	

			deriv_tmp[0] += SP_MM_C*6./pow(dist,8.)*(xm-xn);
			deriv_tmp[1] += SP_MM_C*6./pow(dist,8.)*(ym-yn);
			deriv_tmp[2] += SP_MM_C*6./pow(dist,8.)*(zm-zn);
		}
		SP_MM_BUCK = SP_SYSTEM_FALSE;
//printf("sp[%d] - mm[%d] f : %lf\t%lf\t%lf\n",n+1,m+1,deriv_tmp[0],deriv_tmp[1],deriv_tmp[2]);
		// force update ... again note that the above calculation estimates derivative of 'n'th sp-core w.r.t 'm'th ion-core
		sp_sys->sp_ion[n].force_by_ion_core[0] -= deriv_tmp[0];  sp_sys->sp_ion[n].force_by_ion_core[1] -= deriv_tmp[1];  sp_sys->sp_ion[n].force_by_ion_core[2] -= deriv_tmp[2];
		sp_sys->classic_ion[m].force_by_sp_core[0] += deriv_tmp[0];  sp_sys->classic_ion[m].force_by_sp_core[1] += deriv_tmp[1];  sp_sys->classic_ion[m].force_by_sp_core[2] += deriv_tmp[2];
        }
    }// End of Phase 3

    return;
}







void sp_cluster_system_write_xyz( sp_cluster_system* sp_sys, FILE* fp, char* fn )
{   
    //printf("xyz file is written : %s\n",fn);
    int offset;
    double origin_energy = 0.;

    // WIRTE NUMBER OF ATOMS IN A SYSTEM
    int atom_number = 0;

    for(int i=0;i<sp_sys->number_of_classic_ion;i++)
    {	if( sp_sys->classic_ion[i].if_shell == SP_SYSTEM_FALSE )
		atom_number++;
	// COUNT ONLY WHEN CORE IS THERE
    }

    //fprintf(fp,"\t%d\n",sp_sys->number_of_sp_ion+sp_sys->number_of_classic_ion);
    fprintf(fp,"\t%d\n",sp_sys->number_of_sp_ion+atom_number);

    // GET Energy of the structure
    for(int i=0;i<sp_sys->number_of_classic_ion;i++)
    {   origin_energy += sp_sys->classic_ion[i].classic_energy_by_sp_core;      // energy contri by sp_core
        origin_energy += sp_sys->classic_ion[i].classic_energy_by_ion_core;     // energy contri by ion_core
    }    
    for(int i=0;i<sp_sys->number_of_sp_ion;i++)
    {   origin_energy += sp_sys->sp_ion[i].classic_energy_by_sp_core;           // energy contri by sp_core
        origin_energy += sp_sys->sp_ion[i].classic_energy_by_ion_core;          // energy contri by ion_core

        origin_energy += sp_sys->sp_ion[i].sp_energy;

        //origin_energy += gsl_vector_get(sp_sys->sp_ion[i].eigen_value,sp_cluster_support_get_lowest_state(sp_sys->sp_ion[i].eigen_value));  // energy contri by sp_elec
    }    
    // this is in eV unit
    fprintf(fp,"SCF DONE %12.9lf\n",origin_energy);

    for(int n=0;n<sp_sys->number_of_sp_ion+sp_sys->number_of_classic_ion;n++)
    {
        offset = n-sp_sys->number_of_classic_ion;
        if( n < sp_sys->number_of_classic_ion )
        {	
		if( sp_sys->classic_ion[n].if_shell == SP_SYSTEM_FALSE )	// if it is core position
		{
			fprintf(fp,"%2s%12.6lf%12.6lf%12.6lf\n", sp_sys->classic_ion[n].atom_name,
			gsl_vector_get(sp_sys->classic_ion[n].core_position,0), gsl_vector_get(sp_sys->classic_ion[n].core_position,1),	gsl_vector_get(sp_sys->classic_ion[n].core_position,2));
		}
        }
	else // this is for printing sp-ions
	{
		offset = n - sp_sys->number_of_classic_ion;
		fprintf(fp,"%2s%12.6lf%12.6lf%12.6lf\n", sp_sys->sp_ion[offset].atom_name,
			gsl_vector_get(sp_sys->sp_ion[offset].core_position,0), gsl_vector_get(sp_sys->sp_ion[offset].core_position,1), gsl_vector_get(sp_sys->sp_ion[offset].core_position,2));
	}
    }
    
    return;
}





void sp_cluster_system_get_next_config_mpi( sp_cluster_system* sp_sys, FILE* fp )
{
    fprintf(fp,"#configuration for the next generation\n");
    fprintf(fp,"%d\t%d\n",sp_sys->number_of_classic_ion,sp_sys->number_of_sp_ion);
    for(int i=0;i<sp_sys->number_of_classic_ion;i++)
    {   //classic ion print coordinate
	
	if( sp_sys->classic_ion[i].if_shell == SP_SYSTEM_FALSE )	// if it is core
	{
		fprintf(fp,"%2s%2s%14.6lf%12.6lf%12.6lf\n", sp_sys->classic_ion[i].atom_name,"c",
		gsl_vector_get(sp_sys->classic_ion[i].core_position,0),gsl_vector_get(sp_sys->classic_ion[i].core_position,1),gsl_vector_get(sp_sys->classic_ion[i].core_position,2));
	}
	else if( sp_sys->classic_ion[i].if_shell == SP_SYSTEM_TRUE )	// if it is shell
	{
		fprintf(fp,"%2s%2s%14.6lf%12.6lf%12.6lf\n", sp_sys->classic_ion[i].atom_name,"s",
		gsl_vector_get(sp_sys->classic_ion[i].core_position,0),gsl_vector_get(sp_sys->classic_ion[i].core_position,1),gsl_vector_get(sp_sys->classic_ion[i].core_position,2));
	}
    }
    for(int i=0;i<sp_sys->number_of_sp_ion;i++)
    {
	fprintf(fp,"%2s%16.6lf%12.6lf%12.6lf\n", sp_sys->sp_ion[i].atom_name,
		gsl_vector_get(sp_sys->sp_ion[i].core_position,0),gsl_vector_get(sp_sys->sp_ion[i].core_position,1),gsl_vector_get(sp_sys->sp_ion[i].core_position,2));
    }

    return;
}


int sp_cluster_system_get_density_support_sub_knot_b_search( double dist, int knot_stride, const double* integral_knot )
{   // dist is in a_0 unit
    int le = 0; 
    int re = knot_stride - 1; 

    if( dist > integral_knot[knot_stride-1] )
        return SP_SYSTEM_FALSE;

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

double sp_cluster_system_get_density_support_sub_radxrad( double dist, double* rad1, double* rad2 )
{   // here dist has to be in bohr unit
    double Return;
    Return = (rad1[0]*pow(dist,3.)+rad1[1]*pow(dist,2.)+rad1[2]*dist+rad1[3])
                *(rad2[0]*pow(dist,3.)+rad2[1]*pow(dist,2.)+rad2[2]*dist+rad2[3]);
    return Return;
}

double sp_cluster_system_get_density_support( sp_cluster_type_sp_ion* sp, double* r_e )
{
    double Return = 0.;
    
    const double dist = sqrt(r_e[0]*r_e[0]+r_e[1]*r_e[1]+r_e[2]*r_e[2]);  //in atomic unit .. bohr
    const double q    = sp->charge_shell;
    const double theta= acos( r_e[2]/dist )*(180./M_PI); // in degree unit ... for calc cos -> cos( deg/(180/M_PI) );
    //const double phi  = acos( r_e[1]/dist/sin( theta/(180./M_PI) ))*(180./M_PI); // in degree unit
    const double phi = atan2( r_e[1], r_e[0] )*(180./M_PI); // HERE !! atan2 function is used instead of using "atan(y/x)" to tell the correct quadrant 
    double ev[4];
    double rs[4];
    double rp[4];
    const int knot = sp_cluster_system_get_density_support_sub_knot_b_search( dist, sp->number_of_knot, sp->knot );

    if( knot == SP_SYSTEM_FALSE )
        return 0.;

    ev[0] = gsl_matrix_get(sp->eigen_vector,0,sp_cluster_support_get_lowest_state(sp->eigen_value));
    ev[1] = gsl_matrix_get(sp->eigen_vector,1,sp_cluster_support_get_lowest_state(sp->eigen_value));
    ev[2] = gsl_matrix_get(sp->eigen_vector,2,sp_cluster_support_get_lowest_state(sp->eigen_value));
    ev[3] = gsl_matrix_get(sp->eigen_vector,3,sp_cluster_support_get_lowest_state(sp->eigen_value));

    rs[0] = sp->radial_s_coefficient[knot][0];  rs[1] = sp->radial_s_coefficient[knot][1];
    rs[2] = sp->radial_s_coefficient[knot][2];  rs[3] = sp->radial_s_coefficient[knot][3];

    rp[0] = sp->radial_p_coefficient[knot][0];  rp[1] = sp->radial_p_coefficient[knot][1];
    rp[2] = sp->radial_p_coefficient[knot][2];  rp[3] = sp->radial_p_coefficient[knot][3];

    Return = ev[0]*ev[0]*sp_cluster_system_get_density_support_sub_radxrad(dist,rs,rs)
           + ev[1]*ev[1]*sp_cluster_system_get_density_support_sub_radxrad(dist,rp,rp)*pow( sin(theta/(180./M_PI))*cos(phi/(180./M_PI)), 2. )
           + ev[2]*ev[2]*sp_cluster_system_get_density_support_sub_radxrad(dist,rp,rp)*pow( sin(theta/(180./M_PI))*sin(phi/(180./M_PI)), 2. )
           + ev[3]*ev[3]*sp_cluster_system_get_density_support_sub_radxrad(dist,rp,rp)*pow( cos(theta/(180./M_PI)), 2. )
        +2*( ev[0]*ev[1]*sp_cluster_system_get_density_support_sub_radxrad(dist,rs,rp)*( sin(theta/(180./M_PI))*cos(phi/(180./M_PI)) )
           + ev[0]*ev[2]*sp_cluster_system_get_density_support_sub_radxrad(dist,rs,rp)*( sin(theta/(180./M_PI))*sin(phi/(180./M_PI)) )
           + ev[0]*ev[3]*sp_cluster_system_get_density_support_sub_radxrad(dist,rs,rp)*( cos(theta/(180./M_PI)) )
           + ev[1]*ev[2]*sp_cluster_system_get_density_support_sub_radxrad(dist,rp,rp)*( sin(theta/(180./M_PI))*cos(phi/(180./M_PI)) )*( sin(theta/(180./M_PI))*sin(phi/(180./M_PI)) )
           + ev[1]*ev[3]*sp_cluster_system_get_density_support_sub_radxrad(dist,rp,rp)*( sin(theta/(180./M_PI))*cos(phi/(180./M_PI)) )*( cos(theta/(180./M_PI)) )
           + ev[2]*ev[3]*sp_cluster_system_get_density_support_sub_radxrad(dist,rp,rp)*( sin(theta/(180./M_PI))*sin(phi/(180./M_PI)) )*( cos(theta/(180./M_PI)) ) );

    return Return*q;
}

void sp_cluster_system_get_density_support_sub_boundry( sp_cluster_system* sp_sys, double* r_cnt, double* pos_mx, double* neg_mx )
{
    double pos_mx_cmp[3];
    double neg_mx_cmp[3];

    const int nsp = sp_sys->number_of_sp_ion;
    const int ncla= sp_sys->number_of_classic_ion;

    for(int i=0;i<nsp;i++)
    {
        if( i == 0 )
        {
            pos_mx[0] = gsl_vector_get(sp_sys->sp_ion[i].core_position,0)/TO_BOHR_RADII - r_cnt[0];
            pos_mx[1] = gsl_vector_get(sp_sys->sp_ion[i].core_position,1)/TO_BOHR_RADII - r_cnt[1];
            pos_mx[2] = gsl_vector_get(sp_sys->sp_ion[i].core_position,2)/TO_BOHR_RADII - r_cnt[2];

            neg_mx[0] = gsl_vector_get(sp_sys->sp_ion[i].core_position,0)/TO_BOHR_RADII - r_cnt[0];
            neg_mx[1] = gsl_vector_get(sp_sys->sp_ion[i].core_position,1)/TO_BOHR_RADII - r_cnt[1];
            neg_mx[2] = gsl_vector_get(sp_sys->sp_ion[i].core_position,2)/TO_BOHR_RADII - r_cnt[2];
        }

        pos_mx_cmp[0] = gsl_vector_get(sp_sys->sp_ion[i].core_position,0)/TO_BOHR_RADII - r_cnt[0];
        pos_mx_cmp[1] = gsl_vector_get(sp_sys->sp_ion[i].core_position,1)/TO_BOHR_RADII - r_cnt[1];
        pos_mx_cmp[2] = gsl_vector_get(sp_sys->sp_ion[i].core_position,2)/TO_BOHR_RADII - r_cnt[2];

        neg_mx_cmp[0] = gsl_vector_get(sp_sys->sp_ion[i].core_position,0)/TO_BOHR_RADII - r_cnt[0];
        neg_mx_cmp[1] = gsl_vector_get(sp_sys->sp_ion[i].core_position,1)/TO_BOHR_RADII - r_cnt[1];
        neg_mx_cmp[2] = gsl_vector_get(sp_sys->sp_ion[i].core_position,2)/TO_BOHR_RADII - r_cnt[2];

        if( pos_mx_cmp[0] > pos_mx[0] ) pos_mx[0] = pos_mx_cmp[0];
        if( pos_mx_cmp[1] > pos_mx[1] ) pos_mx[1] = pos_mx_cmp[1];
        if( pos_mx_cmp[2] > pos_mx[2] ) pos_mx[2] = pos_mx_cmp[2];

        if( neg_mx_cmp[0] < neg_mx[0] ) neg_mx[0] = neg_mx_cmp[0];
        if( neg_mx_cmp[1] < neg_mx[1] ) neg_mx[1] = neg_mx_cmp[1];
        if( neg_mx_cmp[2] < neg_mx[2] ) neg_mx[2] = neg_mx_cmp[2];
    }
    for(int i=0;i<ncla;i++)
    {
        pos_mx_cmp[0] = gsl_vector_get(sp_sys->classic_ion[i].core_position,0)/TO_BOHR_RADII - r_cnt[0];
        pos_mx_cmp[1] = gsl_vector_get(sp_sys->classic_ion[i].core_position,1)/TO_BOHR_RADII - r_cnt[1];
        pos_mx_cmp[2] = gsl_vector_get(sp_sys->classic_ion[i].core_position,2)/TO_BOHR_RADII - r_cnt[2];

        neg_mx_cmp[0] = gsl_vector_get(sp_sys->classic_ion[i].core_position,0)/TO_BOHR_RADII - r_cnt[0];
        neg_mx_cmp[1] = gsl_vector_get(sp_sys->classic_ion[i].core_position,1)/TO_BOHR_RADII - r_cnt[1];
        neg_mx_cmp[2] = gsl_vector_get(sp_sys->classic_ion[i].core_position,2)/TO_BOHR_RADII - r_cnt[2];

        if( pos_mx_cmp[0] > pos_mx[0] ) pos_mx[0] = pos_mx_cmp[0];
        if( pos_mx_cmp[1] > pos_mx[1] ) pos_mx[1] = pos_mx_cmp[1];
        if( pos_mx_cmp[2] > pos_mx[2] ) pos_mx[2] = pos_mx_cmp[2];

        if( neg_mx_cmp[0] < neg_mx[0] ) neg_mx[0] = neg_mx_cmp[0];
        if( neg_mx_cmp[1] < neg_mx[1] ) neg_mx[1] = neg_mx_cmp[1];
        if( neg_mx_cmp[2] < neg_mx[2] ) neg_mx[2] = neg_mx_cmp[2];
    }
    // return pos & neg mx
    return;
}

void sp_cluster_system_get_density( sp_cluster_system* sp_sys, const double voxel_in, FILE* fp )
{
    // Note that the type of output is same with Gaussian cube file ... and all the values are in atomic unit

    const double vx = voxel_in; // may be the default = 0.188973 bohr ;
    int grid[3];
    double data;

    double r_obs[3];
    double r_cnt[3] = {0,0,0};
    double r_e[3];

    double pos_mx[3];
    double neg_mx[3];

    int atom_number;

	char element_name[100][4] = {"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
	"Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe",
	"Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb",
	"Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba",
	"La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy","Ho", "Er", "Tm", "Yb", "Lu",
	"Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn"};

    // get cent coord (in bohr)
    for(int i=0;i<sp_sys->number_of_sp_ion;i++)
    {   r_cnt[0] += gsl_vector_get(sp_sys->sp_ion[i].core_position,0)/TO_BOHR_RADII;
        r_cnt[1] += gsl_vector_get(sp_sys->sp_ion[i].core_position,1)/TO_BOHR_RADII;
        r_cnt[2] += gsl_vector_get(sp_sys->sp_ion[i].core_position,2)/TO_BOHR_RADII;      }
    for(int i=0;i<sp_sys->number_of_classic_ion;i++)
    {   r_cnt[0] += gsl_vector_get(sp_sys->classic_ion[i].core_position,0)/TO_BOHR_RADII;
        r_cnt[1] += gsl_vector_get(sp_sys->classic_ion[i].core_position,1)/TO_BOHR_RADII;
        r_cnt[2] += gsl_vector_get(sp_sys->classic_ion[i].core_position,2)/TO_BOHR_RADII;      }

    r_cnt[0] = r_cnt[0]/(sp_sys->number_of_sp_ion+sp_sys->number_of_classic_ion);
    r_cnt[1] = r_cnt[1]/(sp_sys->number_of_sp_ion+sp_sys->number_of_classic_ion);
    r_cnt[2] = r_cnt[2]/(sp_sys->number_of_sp_ion+sp_sys->number_of_classic_ion);

    memset(r_cnt,0.,3*sizeof(double));

    sp_cluster_system_get_density_support_sub_boundry( sp_sys vIN , r_cnt vIN , pos_mx vOUT, neg_mx vOUT ); // get pos neg boundary
    // calculating box   add 3 bohr dist from each edges of the box
    pos_mx[0] += 3./TO_BOHR_RADII;    pos_mx[1] += 3./TO_BOHR_RADII;    pos_mx[2] += 3./TO_BOHR_RADII;
    neg_mx[0] -= 3./TO_BOHR_RADII;    neg_mx[1] -= 3./TO_BOHR_RADII;    neg_mx[2] -= 3./TO_BOHR_RADII;
    // get grid number
    grid[0] = (int)(( pos_mx[0] - neg_mx[0] )/vx) + 1; // get number of grids on x direction
    grid[1] = (int)(( pos_mx[1] - neg_mx[1] )/vx) + 1; // get number of grids on y direction
    grid[2] = (int)(( pos_mx[2] - neg_mx[2] )/vx) + 1; // get number of grids on z direction

    // reset org
    r_cnt[0] = neg_mx[0];   r_cnt[1] = neg_mx[1];   r_cnt[2] = neg_mx[2];

    // file printing
    fprintf(fp,"SLAM CUBE FILE.\n");
    fprintf(fp,"COMMENT\n");
    fprintf(fp,"%4.d%12.6lf%12.6lf%12.6lf\n",
            sp_sys->number_of_classic_ion+sp_sys->number_of_sp_ion,r_cnt[0],r_cnt[1],r_cnt[2]); // Num of atoms , origin of the box x,y,z
    fprintf(fp,"%4.d%12.6lf%12.6lf%12.6lf\n",grid[0],vx,0.,0.); // Num of grids on x direction
    fprintf(fp,"%4.d%12.6lf%12.6lf%12.6lf\n",grid[1],0.,vx,0.); // Num of grids on y direction
    fprintf(fp,"%4.d%12.6lf%12.6lf%12.6lf\n",grid[2],0.,0.,vx); // Num of grids on z direction
 
    // printing classic ion coord (/bohr)
    for(int i=0;i<sp_sys->number_of_classic_ion;i++)
    {   
	for(int j=0;j<100;j++)
	{	if( strcmp( element_name[j], sp_sys->classic_ion[i].atom_name ) == 0 )
		{	atom_number = j+1;
			break;
		}
	}

	fprintf(fp,"%4.d%12.6lf%12.6lf%12.6lf%12.6lf\n",atom_number,0.,gsl_vector_get(sp_sys->classic_ion[i].core_position,0)/TO_BOHR_RADII,
		gsl_vector_get(sp_sys->classic_ion[i].core_position,1)/TO_BOHR_RADII,gsl_vector_get(sp_sys->classic_ion[i].core_position,2)/TO_BOHR_RADII); 
    }
    // printing sp ion coord
    for(int i=0;i<sp_sys->number_of_sp_ion;i++)
    {   
	for(int j=0;j<100;j++)
	{	if( strcmp( element_name[j], sp_sys->sp_ion[i].atom_name ) == 0 )
		{	atom_number = j+1;
			break;
		}
	}

        fprintf(fp,"%4.d%12.6lf%12.6lf%12.6lf%12.6lf\n",atom_number,0.,gsl_vector_get(sp_sys->sp_ion[i].core_position,0)/TO_BOHR_RADII,
		gsl_vector_get(sp_sys->sp_ion[i].core_position,1)/TO_BOHR_RADII,gsl_vector_get(sp_sys->sp_ion[i].core_position,2)/TO_BOHR_RADII);  
    }

    // printing electronic potential by sp - lone pairs
    for(int ix=0;ix<grid[0];ix++)
    {   for(int iy=0;iy<grid[1];iy++)
        {   for(int iz=0;iz<grid[2];iz++)
            {
                data = 0.;  // refresh
                // Get observation point
                r_obs[0] = r_cnt[0] + vx*(double)ix;
                r_obs[1] = r_cnt[1] + vx*(double)iy;
                r_obs[2] = r_cnt[2] + vx*(double)iz;
    
                for(int n=0;n<sp_sys->number_of_sp_ion;n++)
                {
                    r_e[0] = r_obs[0] - gsl_vector_get(sp_sys->sp_ion[n].core_position,0)/TO_BOHR_RADII; 
                    r_e[1] = r_obs[1] - gsl_vector_get(sp_sys->sp_ion[n].core_position,1)/TO_BOHR_RADII; 
                    r_e[2] = r_obs[2] - gsl_vector_get(sp_sys->sp_ion[n].core_position,2)/TO_BOHR_RADII; 
                    data += sp_cluster_system_get_density_support( &sp_sys->sp_ion[n], &r_e[0] );
                }
                fprintf(fp,"%e  ",data);
                if( iz%6 == 5 )
                    fprintf(fp,"\n");
            }
            fprintf(fp,"\n");
        }
    }

    return;
}




// Optimiser // Optimiser // Optimiser // Optimiser // Optimiser // Optimiser // Optimiser // Optimiser // Optimiser // Optimiser // Optimiser // Optimiser 
// Optimiser // Optimiser // Optimiser // Optimiser // Optimiser // Optimiser // Optimiser // Optimiser // Optimiser // Optimiser // Optimiser // Optimiser 
// Optimiser // Optimiser // Optimiser // Optimiser // Optimiser // Optimiser // Optimiser // Optimiser // Optimiser // Optimiser // Optimiser // Optimiser 
// Optimiser // Optimiser // Optimiser // Optimiser // Optimiser // Optimiser // Optimiser // Optimiser // Optimiser // Optimiser // Optimiser // Optimiser













// Get Energy
double sp_cluster_system_get_cluster_energy( sp_cluster_system* sp_sys )
{
    double origin_energy = 0.;

    for(int i=0;i<sp_sys->number_of_classic_ion;i++)
    {   origin_energy += sp_sys->classic_ion[i].classic_energy_by_sp_core;      // energy contri by sp_core
        origin_energy += sp_sys->classic_ion[i].classic_energy_by_ion_core;     // energy contri by ion_core
    }    
    for(int i=0;i<sp_sys->number_of_sp_ion;i++)
    {   origin_energy += sp_sys->sp_ion[i].classic_energy_by_sp_core;           // energy contri by sp_core
        origin_energy += sp_sys->sp_ion[i].classic_energy_by_ion_core;          // energy contri by ion_core

        origin_energy += sp_sys->sp_ion[i].sp_energy;

        //origin_energy += gsl_vector_get(sp_sys->sp_ion[i].eigen_value,sp_cluster_support_get_lowest_state(sp_sys->sp_ion[i].eigen_value));  // energy contri by sp_elec
    }    
    
    return origin_energy;
}




    // ALGORITHM DESCRIPTION - SIMPLE BACK-TRACKING
    /*
     *  Define objective function := f(x), and at kth iteration f_k := f(x_k)
     *
     *  INIT
     *      
     *      SET     t = t_init      // normally t_init = 1.
     *              
     *              const alpha = 0.5
     *
     *              const beta  = 0.95
     *
     *  Step 1.
     *  
     *      t = t_init
     *
     *      IF  THE INITAL grad[f(x_k)] is too huge
     *
     *          e.g.,   || t*grad[f(x_k)] || > 0.5 (0.5 Angstrom)
     *
     *      THEN    t = beta*t
     *
     *  Step 2.
     *
     *      IF  f( x_k - t*grad[f(x_k)] )    >   f(x_k) - alpha*t*||grad[f(x_k)]||^2 
     *
     *      THEN    t = beta*t
     *
     *      ELSE    UPDATE 'x_k' BY
     *          
     *              x_k+1 = x - t*grad[f(x_k)]
     *
     */


void sp_cluster_system_bfgs_support_load_x( sp_cluster_system* sp_sys, gsl_vector* x )
{
    int offset;
    const int number_of_classic_ion = sp_sys->number_of_classic_ion; const int number_of_sp_ion = sp_sys->number_of_sp_ion;

    for(int i=0;i<number_of_classic_ion+number_of_sp_ion;i++)
    {
        if( i < number_of_classic_ion )
        {
            gsl_vector_set(x,i*3+0, gsl_vector_get(sp_sys->classic_ion[i].core_position,0));
            gsl_vector_set(x,i*3+1, gsl_vector_get(sp_sys->classic_ion[i].core_position,1));
            gsl_vector_set(x,i*3+2, gsl_vector_get(sp_sys->classic_ion[i].core_position,2));
        }
        else
        {   offset = i-number_of_classic_ion;
            gsl_vector_set(x,i*3+0, gsl_vector_get(sp_sys->sp_ion[offset].core_position,0));
            gsl_vector_set(x,i*3+1, gsl_vector_get(sp_sys->sp_ion[offset].core_position,1));
            gsl_vector_set(x,i*3+2, gsl_vector_get(sp_sys->sp_ion[offset].core_position,2));
        }
    }
    return;
}


void sp_cluster_system_bfgs_support_load_g( sp_cluster_system* sp_sys, gsl_vector* g )
{   
    int offset;
    const int number_of_classic_ion = sp_sys->number_of_classic_ion;    const int number_of_sp_ion = sp_sys->number_of_sp_ion;

    for(int i=0;i<number_of_classic_ion+number_of_sp_ion;i++)
    {   
        if( i< number_of_classic_ion )
        {   

		if( sp_sys->classic_ion[i]._fix_flag == SP_SYSTEM_FALSE )	// if atom is not fixed
		{
		    gsl_vector_set(g,i*3+0,-sp_sys->classic_ion[i].elec_force_by_sp[0]-sp_sys->classic_ion[i].force_by_sp_core[0]-sp_sys->classic_ion[i].force_by_ion_core[0]);
		    gsl_vector_set(g,i*3+1,-sp_sys->classic_ion[i].elec_force_by_sp[1]-sp_sys->classic_ion[i].force_by_sp_core[1]-sp_sys->classic_ion[i].force_by_ion_core[1]);
		    gsl_vector_set(g,i*3+2,-sp_sys->classic_ion[i].elec_force_by_sp[2]-sp_sys->classic_ion[i].force_by_sp_core[2]-sp_sys->classic_ion[i].force_by_ion_core[2]);
		}
        }   // gradient on ion core
        else
        {   offset = i-number_of_classic_ion;
		
		if( sp_sys->sp_ion[offset]._fix_flag == SP_SYSTEM_FALSE )	// if atom is not fixed
		{
		    gsl_vector_set(g,i*3+0,-sp_sys->sp_ion[offset].elec_force_by_sp[0]-sp_sys->sp_ion[offset].elec_force_by_ion[0]-sp_sys->sp_ion[offset].force_by_sp_core[0]-sp_sys->sp_ion[offset].force_by_ion_core[0]);
		    gsl_vector_set(g,i*3+1,-sp_sys->sp_ion[offset].elec_force_by_sp[1]-sp_sys->sp_ion[offset].elec_force_by_ion[1]-sp_sys->sp_ion[offset].force_by_sp_core[1]-sp_sys->sp_ion[offset].force_by_ion_core[1]);
		    gsl_vector_set(g,i*3+2,-sp_sys->sp_ion[offset].elec_force_by_sp[2]-sp_sys->sp_ion[offset].elec_force_by_ion[2]-sp_sys->sp_ion[offset].force_by_sp_core[2]-sp_sys->sp_ion[offset].force_by_ion_core[2]);
		}
        }   // gradient on sp core
    }
    return;
}



double sp_cluster_system_bfgs_support_get_alpha_mpi_ls( sp_cluster_system* sp_sys, gsl_vector* p_k, double stepmx, int rank, int numtasks )
{
    int offset;
    const int _tmx_ = 120;                          // Internal trial cnt max
    // Backup memory space to save the original coordinates;
    const int number_of_classic_ion = sp_sys->number_of_classic_ion;
    const int number_of_sp_ion      = sp_sys->number_of_sp_ion;
    const int BFGS_Stride = 3*(number_of_sp_ion+number_of_classic_ion);
    const double c2 = 0.9;
    double curvature_LT, curvature_RT;

    gsl_vector* g_trial = gsl_vector_calloc(BFGS_Stride);
    gsl_vector* g_org   = gsl_vector_calloc(BFGS_Stride);
    // Load G_ORG
    sp_cluster_system_bfgs_support_load_g( sp_sys, g_org );

    double** origin_xyz = (double**)malloc((number_of_classic_ion+number_of_sp_ion)*sizeof(double*));    // classic ions + 1 sp-lone pair ion
    double** next_xyz   = (double**)malloc((number_of_classic_ion+number_of_sp_ion)*sizeof(double*));    // saving next step config
    double** grad       = (double**)malloc((number_of_classic_ion+number_of_sp_ion)*sizeof(double*));    // saving gradient on each atom
    for(int i=0;i<number_of_classic_ion+number_of_sp_ion;i++)
    {   origin_xyz[i] = (double*)calloc(3,sizeof(double));
        next_xyz[i]   = (double*)calloc(3,sizeof(double));
        grad[i]       = (double*)calloc(3,sizeof(double));
    }

    double origin_energy, trial_energy;
    double total_norm_square; // ||grad f(x)||^2
    double backtrack_std;   // bactracking cost function

    const double alpha = 0.5;   const double beta = 0.85;   double t;
    int trial_counter = 0;  // this is the internal counter
    int ls_counter = 0;
    /// END OF VARIABLE DERCLARATION

    //gsl_blas_ddot(p_k,g_org,&total_norm_square);

    // Backup the inital ion configuration into 'origin_xyz'
    // 1. back up coordinates
    for(int i=0;i<(number_of_classic_ion+number_of_sp_ion);i++)
    {   if( i < number_of_classic_ion )
        {   origin_xyz[i][0] = gsl_vector_get(sp_sys->classic_ion[i].core_position,0);
            origin_xyz[i][1] = gsl_vector_get(sp_sys->classic_ion[i].core_position,1);
            origin_xyz[i][2] = gsl_vector_get(sp_sys->classic_ion[i].core_position,2);
        }
        else
        {   offset = i-number_of_classic_ion;
            origin_xyz[i][0] = gsl_vector_get(sp_sys->sp_ion[offset].core_position,0);
            origin_xyz[i][1] = gsl_vector_get(sp_sys->sp_ion[offset].core_position,1);
            origin_xyz[i][2] = gsl_vector_get(sp_sys->sp_ion[offset].core_position,2);
        }
    }
    // End Backup

    // Get original struct energy
    origin_energy = sp_cluster_system_get_cluster_energy(sp_sys);
    // Get Energy

    // Initialising Derivative(-Force);
    // get geometrical derivatives acting on atoms
    for(int i=0;i<number_of_classic_ion+number_of_sp_ion;i++)
    {   grad[i][0] = gsl_vector_get(p_k,3*i+0);
        grad[i][1] = gsl_vector_get(p_k,3*i+1);
        grad[i][2] = gsl_vector_get(p_k,3*i+2);
    }
    // A Comment here for avoiding reader's confusion: the original vector 'p_k' is 
    // already 'force' calculated by the equation
    //
    // p_k = B^-1 * ( - grad f(x_k) );
    //
    // therefore here 'grad' isn't really grad but force
    
    // Set Initial t - value
    t = 1.;

    // trial move
    while(1)
    {   trial_counter++;
        // step size check
        while(1)
        {   
            int boolean = 1;    // if step size smaller then 'stepmx'? if false then 'boolean' = 0;
            
            for(int i=0;i<number_of_classic_ion+number_of_sp_ion;i++)
            {   
                if( t*sp_cluster_support_get_norm(grad[i][0],grad[i][1],grad[i][2]) > stepmx )
                {   t = beta*t; // if the step size is too huge then rescaling
                    boolean = 0;
                }
            }
            if( boolean == 1 )  // all step sizes less then 'stepmx'
            {   
                break;          // if all t*gnorm less then 'stepmx' then break
            }
        }
        // step size check done

        total_norm_square=0.;
        // trial (next) move coordinates
        //
        // Set the trial xyz coordinates "EXCEPT THE LAST SP_ION CENTRE": -> OPTIMISING REFERENCE COORDINATE!!

	// MOVING ATOMS
        for(int i=0;i<number_of_classic_ion+number_of_sp_ion;i++)
        {   
            // next xyz  =   xk + alpha(t)*p_k;
            // s_k       =   alpha(t)*p_k

            if( i < number_of_classic_ion )
            {
		if( sp_sys->classic_ion[i]._fix_flag == SP_SYSTEM_FALSE )	// if atom is not fixed
		{
			next_xyz[i][0] = origin_xyz[i][0] + t*grad[i][0];
			next_xyz[i][1] = origin_xyz[i][1] + t*grad[i][1];
			next_xyz[i][2] = origin_xyz[i][2] + t*grad[i][2];                       // calc tiral move classic ion coordinates
			gsl_vector_set(sp_sys->classic_ion[i].core_position,0,next_xyz[i][0]);
			gsl_vector_set(sp_sys->classic_ion[i].core_position,1,next_xyz[i][1]);
			gsl_vector_set(sp_sys->classic_ion[i].core_position,2,next_xyz[i][2]);  // put calculated coordinates into sp_sys class
			total_norm_square += (grad[i][0]*grad[i][0] + grad[i][1]*grad[i][1] + grad[i][2]*grad[i][2]);
		}
            }
            else
            {   
		offset = i-number_of_classic_ion;

		if( sp_sys->sp_ion[offset]._fix_flag == SP_SYSTEM_FALSE )	// if atom is not fixed
		{
			next_xyz[i][0] = origin_xyz[i][0] + t*grad[i][0];
			next_xyz[i][1] = origin_xyz[i][1] + t*grad[i][1];
			next_xyz[i][2] = origin_xyz[i][2] + t*grad[i][2];                  // calc tiral move sp ion coordinates
			gsl_vector_set(sp_sys->sp_ion[offset].core_position,0,next_xyz[i][0]);
			gsl_vector_set(sp_sys->sp_ion[offset].core_position,1,next_xyz[i][1]);
			gsl_vector_set(sp_sys->sp_ion[offset].core_position,2,next_xyz[i][2]);  // put calculated coordinates into sp_sys class
			total_norm_square += (grad[i][0]*grad[i][0] + grad[i][1]*grad[i][1] + grad[i][2]*grad[i][2]);
		}
            }
            // total norm sqr
            //total_norm_square += (grad[i][0]*grad[i][0] + grad[i][1]*grad[i][1] + grad[i][2]*grad[i][2]);
        }
        MPI_Barrier(MPI_COMM_WORLD);    // SET SYNC
        
        /// CALCULATE TRIAL ENERGY
        //TEST
        sp_cluster_system_set_is_scf_done(sp_sys,SP_SYSTEM_FALSE);  // FLAG -> FALSE means to recalibrate the SCF w.r.t the new trial position 
        //
        sp_cluster_system_scf_mpi(sp_sys,rank,numtasks);    // Get SCF Achieved sp set
        // At This point, Molucular Orbital EigenVectors are fixed
        sp_cluster_system_get_classic_energy(sp_sys);       // Get Classical energy
        MPI_Barrier(MPI_COMM_WORLD);   
        trial_energy = sp_cluster_system_get_cluster_energy(sp_sys);   // Get 'trial energy'
        
        // Calculate:  f(r_origin) - alpha*t*||(grad_x,grad_y,grad_z)||^2 -> I will denote this with 'backtrack_std'
        backtrack_std = origin_energy - alpha*t*total_norm_square;
        //backtrack_std = origin_energy + alpha*t*total_norm_square;
        MPI_Barrier(MPI_COMM_WORLD);   

        // Backtracking Criteria
        if( trial_energy > backtrack_std )
            t = beta*t;
        else //if( trial_energy <= backtrack_std || trial_energy < origin_energy ) // trial move accepted or whatever the criteria is if the trial energy is lower then the original energy
        {   
            //if(rank==0) printf("Trial_counter: %d \t origin_energy/trial_energy/backtrack_std/t: %18.12lf%18.12lf%18.12lf%18.9e\n",trial_counter,origin_energy,trial_energy,backtrack_std,t);
            
            ls_counter++;     // this is for using LS with Wolfe - condition
            //ls_counter = 5; // this is for only using armijo line-search ... to use it remove the commenting
            // WOLFE - CONDITION ... Curvature Check

            // Derivative Update
            sp_cluster_system_get_force_mpi(sp_sys,rank,numtasks);
            sp_cluster_system_get_classic_force(sp_sys);
            MPI_Barrier(MPI_COMM_WORLD);   
	    // Derivative Calculation Done

            // only check Wolfe Condition when MULTI CENTRE ROUTINE IS USED ...
            if( number_of_sp_ion > 1 )
            {
                curvature_LT = 0.; curvature_RT = 0.;
                // Load G Trial
                sp_cluster_system_bfgs_support_load_g( sp_sys, g_trial );
                gsl_blas_ddot(p_k,g_trial,&curvature_LT);
                gsl_blas_ddot(p_k,g_org  ,&curvature_RT);

                //if( -1.*curvature_LT <= -1.*c2*curvature_RT )
                if( -1.*curvature_LT <= -1.*c2*curvature_RT || ls_counter > 15 )
                {
                    if(rank==0) 
                    {   
                        if( ls_counter > 15 )
                        {   
                            printf("\n");
                            printf(" Armijo-Condition is Achieved !\n\n");
                            printf(" Trial Count            :    %d\n",trial_counter);
                            printf(" Original Energy (eV)   :   %.12lf\n",origin_energy);
                            printf(" Trial    Energy (eV)   :   %.12lf\n",trial_energy);
                            printf(" Diff     Energy (eV)   :   %.12lf\n\n",trial_energy-origin_energy);
                        }
                        else
                        {   printf("\n");
                            printf(" Wolfe-Condition is Achieved !\n\n");
                            printf(" Trial Count            :    %d\n",trial_counter);
                            printf(" Original Energy (eV)   :   %.12lf\n",origin_energy);
                            printf(" Trial    Energy (eV)   :   %.12lf\n",trial_energy);
                            printf(" Diff     Energy (eV)   :   %.12lf\n\n",trial_energy-origin_energy);
                        }
                    }
                    MPI_Barrier(MPI_COMM_WORLD);
                    break;  // Break; Internal Backtracking Search While Loop
                }

                t = beta*t; // if the wolfe-condition is not achieved then t = beta*t
                //if(rank==0) printf("Trial_counter: %d \t trial_energy/backtrack_std/t: %lf \t %lf \t %.9e\n",trial_counter,trial_energy,backtrack_std,t);
                MPI_Barrier(MPI_COMM_WORLD);   
            }
            else    // if number_of_sp_ion == 1     // Armijo Backtracking Search
            {
                if(rank==0) 
                {   
                    printf("\n");
                    printf(" Armijo-Condition is Achieved !\n\n");
                    printf(" Trial Count            :    %d\n",trial_counter);
                    printf(" Original Energy (eV)   :   %.12lf\n",origin_energy);
                    printf(" Trial    Energy (eV)   :   %.12lf\n",trial_energy);
                    printf(" Diff     Energy (eV)   :   %.12lf\n\n",trial_energy-origin_energy);
                }
                //t = SP_SYSTEM_FALSE;    // Armijo condition doesnt guarantee positive definite BFGS update
                MPI_Barrier(MPI_COMM_WORLD);
                break;  // Break; Internal Backtracking Search While Loop
            }

        } // END OF TRIAL MOVE ACCEPTED IF SENTENCE

        if( t < 10E-9  )  //0.025?
        {   
            if( rank == 0 )
            {   
                printf("\n");
                printf(" Line Search Failed (too small step size) : %.6e ...\n\n",t);
                printf(" Trial Count    :    %d\n",trial_counter);
                printf(" Original Energy:   %.8lf\n",origin_energy);
                printf(" Trial    Energy:   %.8lf\n\n",trial_energy);
                printf("\n");
            }
            // Maybe this move is in the scale of 10E-9 ... might be just negligible

            for(int i=0;i<(number_of_classic_ion+number_of_sp_ion);i++)
            {   if( i < number_of_classic_ion )
                {   
                    gsl_vector_set(sp_sys->classic_ion[i].core_position,0,origin_xyz[i][0]);
                    gsl_vector_set(sp_sys->classic_ion[i].core_position,1,origin_xyz[i][1]);
                    gsl_vector_set(sp_sys->classic_ion[i].core_position,2,origin_xyz[i][2]);
                }
                else
                {   offset = i-number_of_classic_ion;
                    gsl_vector_set(sp_sys->sp_ion[offset].core_position,0,origin_xyz[i][0]);
                    gsl_vector_set(sp_sys->sp_ion[offset].core_position,1,origin_xyz[i][1]);
                    gsl_vector_set(sp_sys->sp_ion[offset].core_position,2,origin_xyz[i][2]);
                }
            }
            // Restore H mat
            sp_cluster_system_scf_mpi(sp_sys,rank,numtasks);    // Get SCF Achieved sp set
            MPI_Barrier(MPI_COMM_WORLD);
            t = SP_SYSTEM_FALSE; // FAIL FLAG
            break;  // Break; Internal Backtracking Search While Loop
        }

    } // trial move done
    MPI_Barrier(MPI_COMM_WORLD);   
    // MEMORY DETACH
    for(int i=0;i<number_of_classic_ion+number_of_sp_ion;i++)
    {   free(origin_xyz[i]);
        free(next_xyz[i]);
        free(grad[i]);
    }
    free(origin_xyz);
    free(next_xyz);
    free(grad);

    gsl_vector_free(g_trial);   gsl_vector_free(g_org);

    MPI_Barrier(MPI_COMM_WORLD);   
    return t;
}


// MPI_SUPPORT added 8/12/2019
// BFGS Algorithm
int sp_cluster_system_call_bfgs_algorithm_mpi( sp_cluster_system* sp_sys, double stepmx /* maximum stepsize */, const int cyclemx ,int rank, int numtasks )
{   int ret = SP_SYSTEM_FALSE; 
    double gnorm_prev, gnorm_next;
    int same_gnorm_cnt = 0;
    double e_prev, e_next;
    int same_e_cnt = 0;

    const double _gnorm_tol_ = SP_SYSTEM_GNORM_TOL;
    const double _sum_sqr_tol_ = 0.05;
    const double _force_tol_ = 0.;
    const int _max_bfgs_step_ = cyclemx;
    const int _max_trial_ = 500;
        
    const int number_of_classic_ion = sp_sys->number_of_classic_ion;
    const int number_of_sp_ion      = sp_sys->number_of_sp_ion;
	
	int deg_free = 0;
	for(int i=0;i<number_of_classic_ion;i++)
	{	if( sp_sys->classic_ion[i]._fix_flag == SP_SYSTEM_FALSE )
			deg_free++;
	}
	for(int i=0;i<number_of_sp_ion;i++)
	{	if( sp_sys->sp_ion[i]._fix_flag == SP_SYSTEM_FALSE )
			deg_free++;
	}
	deg_free = deg_free*3;


    const int BFGS_Stride = 3*(number_of_sp_ion+number_of_classic_ion);
    int offset;
    // BFGS WORKSPACE
    int if_bfgs = SP_SYSTEM_FALSE;
    double alpha=0.;
    double grad_x,grad_y,grad_z;
    double s_tmp1, s_tmp2;
    double sum_sqr;
    gsl_matrix* m_tmp1 = gsl_matrix_calloc(BFGS_Stride,BFGS_Stride);
    gsl_matrix* m_tmp2 = gsl_matrix_calloc(BFGS_Stride,BFGS_Stride);
    gsl_matrix* m_tmp3 = gsl_matrix_calloc(BFGS_Stride,BFGS_Stride);
    gsl_matrix* m_tmp4 = gsl_matrix_calloc(BFGS_Stride,BFGS_Stride);
    gsl_vector* v_tmp1 = gsl_vector_calloc(BFGS_Stride);

    gsl_matrix* Inv_B_k      = gsl_matrix_calloc(BFGS_Stride,BFGS_Stride);         // B_k^-1 Matrix (Approximated Hessian)
    gsl_vector* p_k          = gsl_vector_calloc(BFGS_Stride);                                                    // p_k    Vector
    gsl_vector* g_k_next     = gsl_vector_calloc(BFGS_Stride);                                                    // g_k+1
    gsl_vector* g_k_prev     = gsl_vector_calloc(BFGS_Stride);                                                    // g_k
    gsl_vector* s_k          = gsl_vector_calloc(BFGS_Stride);
    gsl_vector* x_k_next     = gsl_vector_calloc(BFGS_Stride);
    gsl_vector* x_k_prev     = gsl_vector_calloc(BFGS_Stride);
    gsl_vector* y_k          = gsl_vector_calloc(BFGS_Stride);
    
    double wtime;

    // Variables for tmp saving eigenvector derivatives
    double evec_tmp[number_of_sp_ion][4];
    double deriv_evec_sp[number_of_sp_ion][number_of_sp_ion][4][3];
    double deriv_evec_cla[number_of_classic_ion][number_of_sp_ion][4][3];
    double deriv_en_sp[number_of_sp_ion][3];
    double deriv_en_cla[number_of_classic_ion][3];
    double dl = 0.00025;
    double ds,dx,dy,dz;
    double prev_en, next_en;
    //sp_cluster_support_get_lowest_state(sp_sys->sp_ion[kk].eigen_value);
    // [i][j][k] ... i=spion, j=evec_elem, k=x,y,z

    MPI_Barrier(MPI_COMM_WORLD);

    // INITIALISING B_k Matrix ... for the first iteration, the matrix is I.
    for(int i=0;i<BFGS_Stride;i++)  gsl_matrix_set(Inv_B_k,i,i,1.);

    // INITIALISING FIRST CONFIGURATION
    sp_cluster_system_set_is_scf_done(sp_sys,SP_SYSTEM_FALSE);
    /// CALCULATE TRIAL ENERGY & FORCE
    sp_cluster_system_scf_mpi(sp_sys,rank,numtasks);    // Get SCF Achieved sp set
    sp_cluster_system_get_classic_energy(sp_sys);
    sp_cluster_system_get_force_mpi(sp_sys,rank,numtasks);
    sp_cluster_system_get_classic_force(sp_sys);

    // Load x_k_prev (x_k)
    sp_cluster_system_bfgs_support_load_x( sp_sys, x_k_prev );
    // Load g_k_prev (g_k)
    sp_cluster_system_bfgs_support_load_g( sp_sys, g_k_prev );

    gnorm_prev = gsl_blas_dnrm2(g_k_prev)/(double)BFGS_Stride;
    gnorm_next = 1.;
    e_prev = sp_cluster_system_get_cluster_energy(sp_sys);
    e_next = 1.;

	
	// if the max cycle requested with '0' take it as a single point calculation !!!
	if( _max_bfgs_step_ == 0 )
	{
		if(rank==0)
		{
			printf("\n");
			printf(" Single point calculation is requested ... \n");
			printf("\n");
			printf(" Position Integral Reference : %.12lf\n",sp_cluster_integrator_get_x_12( &sp_sys->sp_ion[0] ));
			printf("\n");
			printf(" Single Point Result\n");
			printf("\n");
			printf(" Geometric Derivatives ( eV / Angstrom )\n");
			printf("\n");
			printf("-----------------------------------------------------------------------------\n");
			printf(" Species.        x           y           z          |r| \n");
			printf("-----------------------------------------------------------------------------\n");

			for(int i=0;i<number_of_classic_ion+number_of_sp_ion;i++)
			{   
				if( i< number_of_classic_ion )
				{   
				    grad_x = -sp_sys->classic_ion[i].elec_force_by_sp[0]-sp_sys->classic_ion[i].force_by_sp_core[0]-sp_sys->classic_ion[i].force_by_ion_core[0];
				    grad_y = -sp_sys->classic_ion[i].elec_force_by_sp[1]-sp_sys->classic_ion[i].force_by_sp_core[1]-sp_sys->classic_ion[i].force_by_ion_core[1];
				    grad_z = -sp_sys->classic_ion[i].elec_force_by_sp[2]-sp_sys->classic_ion[i].force_by_sp_core[2]-sp_sys->classic_ion[i].force_by_ion_core[2];
				    if( sp_sys->classic_ion[i].if_shell == SP_SYSTEM_TRUE )		// if the MM type is shell
					printf("%3s%3s%15.6e%15.6e%15.6e%15.6e\n",sp_sys->classic_ion[i].atom_name,"s",grad_x,grad_y,grad_z,sqrt(grad_x*grad_x+grad_y*grad_y+grad_z*grad_z));
				    else if( sp_sys->classic_ion[i].if_shell == SP_SYSTEM_FALSE )	// if the MM type is core
					printf("%3s%3s%15.6e%15.6e%15.6e%15.6e\n",sp_sys->classic_ion[i].atom_name,"c",grad_x,grad_y,grad_z,sqrt(grad_x*grad_x+grad_y*grad_y+grad_z*grad_z));

				    //Backing Up - SinglePoint if
				    deriv_en_cla[i][0] = grad_x;    deriv_en_cla[i][1] = grad_y;    deriv_en_cla[i][2] = grad_z;

				}   // gradient on ion core
				else
				{   offset = i-number_of_classic_ion;
				    grad_x = -sp_sys->sp_ion[offset].elec_force_by_sp[0]-sp_sys->sp_ion[offset].elec_force_by_ion[0]-sp_sys->sp_ion[offset].force_by_sp_core[0]-sp_sys->sp_ion[offset].force_by_ion_core[0];
				    grad_y = -sp_sys->sp_ion[offset].elec_force_by_sp[1]-sp_sys->sp_ion[offset].elec_force_by_ion[1]-sp_sys->sp_ion[offset].force_by_sp_core[1]-sp_sys->sp_ion[offset].force_by_ion_core[1];
				    grad_z = -sp_sys->sp_ion[offset].elec_force_by_sp[2]-sp_sys->sp_ion[offset].elec_force_by_ion[2]-sp_sys->sp_ion[offset].force_by_sp_core[2]-sp_sys->sp_ion[offset].force_by_ion_core[2];
				    printf("%3s%18.6e%15.6e%15.6e%15.6e\n",sp_sys->sp_ion[offset].atom_name,grad_x,grad_y,grad_z,sqrt(grad_x*grad_x+grad_y*grad_y+grad_z*grad_z));
                    //Backing Up - SinglePoint if
                    deriv_en_sp[offset][0] = grad_x;    deriv_en_sp[offset][1] = grad_y;    deriv_en_sp[offset][2] = grad_z;
				}   // gradient on sp core
			}
			printf("-----------------------------------------------------------------------------\n");
			printf("\n");
			printf(" Gnorm  (eV/Angs)   :   %.6lf\n",gsl_blas_dnrm2(g_k_prev)/(double)BFGS_Stride);
			printf(" Energy (eV)        :  %.12lf\n",sp_cluster_system_get_cluster_energy(sp_sys));
			printf("\n");
			printf("-----------------------------------------------------------------------------\n");
			printf("\n");
			printf(" Lone Pair Molecular Orbital Info ( Lowest EigenValue / EigenVector Set )\n");
			printf("\n");
			printf("-----------------------------------------------------------------------------\n");
			printf(" Species.  Energy(eV)        s           px           py           pz \n");
			printf("-----------------------------------------------------------------------------\n");

			for(int kk=0;kk<number_of_sp_ion;kk++)
			{   int low_stat = sp_cluster_support_get_lowest_state(sp_sys->sp_ion[kk].eigen_value);
				printf("%3s%18.8lf%12.8lf%13.8lf%13.8lf%13.8lf\n",sp_sys->sp_ion[kk].atom_name,gsl_vector_get(sp_sys->sp_ion[kk].eigen_value,low_stat),
					gsl_matrix_get(sp_sys->sp_ion[kk].eigen_vector,0,low_stat),gsl_matrix_get(sp_sys->sp_ion[kk].eigen_vector,1,low_stat),
					gsl_matrix_get(sp_sys->sp_ion[kk].eigen_vector,2,low_stat),gsl_matrix_get(sp_sys->sp_ion[kk].eigen_vector,3,low_stat));
			}
            printf("-----------------------------------------------------------------------------\n");
            // CONVENTION ... sp_sys->deriv_evec_sp[0][1] ...  is derivatives of [1] atom when [0] is moved
            // BELOW NEW
            printf("\n");
            printf(" (1) EigenVector Derivatives ( 1/Angs )\n");
            printf("\n");
            printf(" ( Changes in MO Coefficients of (i) versus the changes in (xyz) of (j) )\n");
            printf("\n");
            printf(" Format\n");
            printf("\n");
            printf("-----------------------------------------------------------------------------\n");
			printf(" sp(i)-sp(j)   Type.           s           px           py           pz \n");
            printf("-----------------------------------------------------------------------------\n");
            printf("\n");
            printf(" (2) Total Energy Geometric Derivatives ( eV/Angs )\n");
            printf("\n");
            printf(" ( Changes in Total Energy versus the changes in (xyz) of (i) )\n");
            printf("\n");
            printf(" Format\n");
            printf("\n");
            printf("-----------------------------------------------------------------------------\n");
			printf(" ion (i)       Type.           Species.\n");
            printf("-----------------------------------------------------------------------------\n");
            printf("\n");
		}
		MPI_Barrier(MPI_COMM_WORLD);

//////////
        // BCAST ORG GEOMETRIC DERIVATIVES
        for(int i=0;i<number_of_sp_ion;i++)
            MPI_Bcast(&deriv_en_sp[i][0],3,MPI_DOUBLE,0,MPI_COMM_WORLD);
        for(int i=0;i<number_of_classic_ion;i++)
            MPI_Bcast(&deriv_en_cla[i][0],3,MPI_DOUBLE,0,MPI_COMM_WORLD);
        // Backing Up the original cluster energy
        prev_en = sp_cluster_system_get_cluster_energy(sp_sys);
        // Backing Up the original MO info
        for(int i=0;i<number_of_sp_ion;i++)
        {   for(int j=0;j<4;j++)
                evec_tmp[i][j] = gsl_matrix_get(sp_sys->sp_ion[i].eigen_vector,j,sp_cluster_support_get_lowest_state(sp_sys->sp_ion[i].eigen_value));
        }
        // Backing up the original MO derivative info
        for(int i=0;i<number_of_sp_ion;i++)
        {   for(int j=0;j<number_of_sp_ion;j++)
            {   for(int u=0;u<4;u++)
                {   for(int v=0;v<3;v++)
                        deriv_evec_sp[i][j][u][v] = sp_sys->deriv_evec_sp[i][j][u][v];
                }
            }
        }
        for(int i=0;i<number_of_classic_ion;i++)
        {   for(int j=0;j<number_of_sp_ion;j++)
            {   for(int u=0;u<4;u++)
                {   for(int v=0;v<3;v++)
                        deriv_evec_cla[i][j][u][v] = sp_sys->deriv_evec_cla[i][j][u][v];
                }
            }
        }
        // FDM CHECKER
        // EVEC COEFF DERIVATIVES               W.R.T SP_ION MOVES 
        // ENERGY GEOMETRIC DERIVATIVES         W.R.T SP_ION MOVES
        for(int i=0;i<number_of_sp_ion;i++)
        {   for(int j=0;j<3;j++)
            {   gsl_vector_set( sp_sys->sp_ion[i].core_position, j, gsl_vector_get(sp_sys->sp_ion[i].core_position,j) + dl );
                sp_cluster_system_set_is_scf_done(sp_sys,SP_SYSTEM_FALSE);
                sp_cluster_system_scf_mpi(sp_sys,rank,numtasks);    // Get SCF Achieved sp set
                sp_cluster_system_get_classic_energy(sp_sys);
                sp_cluster_system_get_force_mpi(sp_sys,rank,numtasks);
                sp_cluster_system_get_classic_force(sp_sys);
         
                // ADD FORCE FDM
                next_en = sp_cluster_system_get_cluster_energy(sp_sys);
                if(rank==0)
                {   
                    printf("-----------------------------------------------------------------------------\n");
                    printf("\n");
                    printf(" Geometric Derivative ( Total Energy )\n");
                    printf("\n");
                    if( j == 0 )
                    {   
			printf(" sp(%d)            d/dx          %s\n\n",i+1,sp_sys->sp_ion[i].atom_name);
                        printf(" \t   Numerical(FDM) %9.6lf\n",(next_en-prev_en)/dl);
                        printf(" \t   Analytical     %9.6lf\n",deriv_en_sp[i][j]);
                        printf(" \t   Error (%%)      %9.6lf%s\n",(deriv_en_sp[i][j]-(next_en-prev_en)/dl)/deriv_en_sp[i][j]*100.," %");
                        printf("\n");
                    }   else if( j == 1 )
                    {   printf(" sp(%d)            d/dy          %s\n\n",i+1,sp_sys->sp_ion[i].atom_name);
                        printf(" \t   Numerical(FDM) %9.6lf\n",(next_en-prev_en)/dl);
                        printf(" \t   Analytical     %9.6lf\n",deriv_en_sp[i][j]);
                        printf(" \t   Error (%%)      %9.6lf%s\n",(deriv_en_sp[i][j]-(next_en-prev_en)/dl)/deriv_en_sp[i][j]*100.," %");
                        printf("\n");
                    }   else
                    {   printf(" sp(%d)            d/dz          %s\n\n",i+1,sp_sys->sp_ion[i].atom_name);
                        printf(" \t   Numerical(FDM) %9.6lf\n",(next_en-prev_en)/dl);
                        printf(" \t   Analytical     %9.6lf\n",deriv_en_sp[i][j]);
                        printf(" \t   Error (%%)      %9.6lf%s\n",(deriv_en_sp[i][j]-(next_en-prev_en)/dl)/deriv_en_sp[i][j]*100.," %");
                        printf("\n");
                    }
                    printf("-----------------------------------------------------------------------------\n");
                }

                for(int k=0;k<number_of_sp_ion;k++)
                {
                    if (rank == 0 )
                    {   
                        printf("-----------------------------------------------------------------------------\n");
                        printf("\n");
                        printf(" MO Coefficient Derivative \n");
                        printf("\n");
                        
                        if( j == 0 )
                            printf(" sp(%d) - sp(%d)   d/dx \n",k+1,i+1);
                        else if ( j == 1 )
                            printf(" sp(%d) - sp(%d)   d/dy \n",k+1,i+1);
                        else
                            printf(" sp(%d) - sp(%d)   d/dz \n",k+1,i+1);
                        printf("\n");
                    }

                    ds = gsl_matrix_get(sp_sys->sp_ion[k].eigen_vector,0,sp_cluster_support_get_lowest_state(sp_sys->sp_ion[k].eigen_value)) - evec_tmp[k][0];
                    dx = gsl_matrix_get(sp_sys->sp_ion[k].eigen_vector,1,sp_cluster_support_get_lowest_state(sp_sys->sp_ion[k].eigen_value)) - evec_tmp[k][1];
                    dy = gsl_matrix_get(sp_sys->sp_ion[k].eigen_vector,2,sp_cluster_support_get_lowest_state(sp_sys->sp_ion[k].eigen_value)) - evec_tmp[k][2];
                    dz = gsl_matrix_get(sp_sys->sp_ion[k].eigen_vector,3,sp_cluster_support_get_lowest_state(sp_sys->sp_ion[k].eigen_value)) - evec_tmp[k][3];

                    ds = ds/dl; dx = dx/dl; dy = dy/dl; dz = dz/dl;
            
                    if( rank == 0 )
                    {
                        //printf(" sp(i)-sp(j)   Type.           s           px           py           pz \n");
                        printf(" \t   Numerical(FDM) %9.6lf%12.6lf%12.6lf%12.6lf\n",ds,dx,dy,dz);
                        printf(" \t   Analytical     %9.6lf%12.6lf%12.6lf%12.6lf\n",deriv_evec_sp[i][k][0][j],deriv_evec_sp[i][k][1][j],deriv_evec_sp[i][k][2][j],deriv_evec_sp[i][k][3][j]);
                        printf(" \t   Error (%%)      %9.6lf%s%10.6lf%s%10.6lf%s%10.6lf%s\n",(deriv_evec_sp[i][k][0][j]-ds)/deriv_evec_sp[i][k][0][j]*100.," %",                                 \
                                (deriv_evec_sp[i][k][1][j]-dx)/deriv_evec_sp[i][k][1][j]*100.," %",(deriv_evec_sp[i][k][2][j]-dy)/deriv_evec_sp[i][k][2][j]*100.," %",      \
                                (deriv_evec_sp[i][k][3][j]-dz)/deriv_evec_sp[i][k][3][j]*100.," %");
                        printf("\n");
                    }
                }
                gsl_vector_set( sp_sys->sp_ion[i].core_position, j, gsl_vector_get(sp_sys->sp_ion[i].core_position,j) - dl );
            }
        }

        // FDM CHECKER
        // EVEC COEFF DERIVATIVES               W.R.T SP_ION MOVES 
        // ENERGY GEOMETRIC DERIVATIVES         W.R.T SP_ION MOVES
        for(int i=0;i<number_of_classic_ion;i++)
        {   for(int j=0;j<3;j++)
            {   gsl_vector_set( sp_sys->classic_ion[i].core_position, j, gsl_vector_get(sp_sys->classic_ion[i].core_position,j) + dl );
                sp_cluster_system_set_is_scf_done(sp_sys,SP_SYSTEM_FALSE);
                sp_cluster_system_scf_mpi(sp_sys,rank,numtasks);    // Get SCF Achieved sp set
                sp_cluster_system_get_classic_energy(sp_sys);
                sp_cluster_system_get_force_mpi(sp_sys,rank,numtasks);
                sp_cluster_system_get_classic_force(sp_sys);
         
                // ADD FORCE FDM
                next_en = sp_cluster_system_get_cluster_energy(sp_sys);
                if(rank==0)
                {   
                    printf("-----------------------------------------------------------------------------\n");
                    printf("\n");
                    printf(" Geometric Derivative ( Total Energy )\n");
                    printf("\n");
                    if( j == 0 )
                    {   printf(" cla(%d)            d/dx          %s\n\n",i+1,sp_sys->classic_ion[i].atom_name);
                        printf(" \t   Numerical(FDM) %9.6lf\n",(next_en-prev_en)/dl);
                        printf(" \t   Analytical     %9.6lf\n",deriv_en_cla[i][j]);
                        printf(" \t   Error (%%)      %9.6lf%s\n",(deriv_en_cla[i][j]-(next_en-prev_en)/dl)/deriv_en_cla[i][j]*100.," %");
                        printf("\n");
                    }   else if( j == 1 )
                    {   printf(" cla(%d)            d/dy          %s\n\n",i+1,sp_sys->classic_ion[i].atom_name);
                        printf(" \t   Numerical(FDM) %9.6lf\n",(next_en-prev_en)/dl);
                        printf(" \t   Analytical     %9.6lf\n",deriv_en_cla[i][j]);
                        printf(" \t   Error (%%)      %9.6lf%s\n",(deriv_en_cla[i][j]-(next_en-prev_en)/dl)/deriv_en_cla[i][j]*100.," %");
                        printf("\n");
                    }   else
                    {   printf(" cla(%d)            d/dz          %s\n\n",i+1,sp_sys->classic_ion[i].atom_name);
                        printf(" \t   Numerical(FDM) %9.6lf\n",(next_en-prev_en)/dl);
                        printf(" \t   Analytical     %9.6lf\n",deriv_en_cla[i][j]);
                        printf(" \t   Error (%%)      %9.6lf%s\n",(deriv_en_cla[i][j]-(next_en-prev_en)/dl)/deriv_en_cla[i][j]*100.," %");
                        printf("\n");
                    }
                    printf("-----------------------------------------------------------------------------\n");
                }

                for(int k=0;k<number_of_sp_ion;k++)
                {
                    if (rank == 0 )
                    {   
                        printf("-----------------------------------------------------------------------------\n");
                        printf("\n");
                        printf(" MO Coefficient Derivative \n");
                        printf("\n");
                        
                        if( j == 0 )
                            printf(" sp(%d) - cla(%d)   d/dx \n",k+1,i+1);
                        else if ( j == 1 )
                            printf(" sp(%d) - cla(%d)   d/dy \n",k+1,i+1);
                        else
                            printf(" sp(%d) - cla(%d)   d/dz \n",k+1,i+1);
                        printf("\n");
                    }

                    ds = gsl_matrix_get(sp_sys->sp_ion[k].eigen_vector,0,sp_cluster_support_get_lowest_state(sp_sys->sp_ion[k].eigen_value)) - evec_tmp[k][0];
                    dx = gsl_matrix_get(sp_sys->sp_ion[k].eigen_vector,1,sp_cluster_support_get_lowest_state(sp_sys->sp_ion[k].eigen_value)) - evec_tmp[k][1];
                    dy = gsl_matrix_get(sp_sys->sp_ion[k].eigen_vector,2,sp_cluster_support_get_lowest_state(sp_sys->sp_ion[k].eigen_value)) - evec_tmp[k][2];
                    dz = gsl_matrix_get(sp_sys->sp_ion[k].eigen_vector,3,sp_cluster_support_get_lowest_state(sp_sys->sp_ion[k].eigen_value)) - evec_tmp[k][3];

                    ds = ds/dl; dx = dx/dl; dy = dy/dl; dz = dz/dl;
            
                    if( rank == 0 )
                    {
                        //printf(" sp(i)-sp(j)   Type.           s           px           py           pz \n");
                        printf(" \t   Numerical(FDM) %9.6lf%12.6lf%12.6lf%12.6lf\n",ds,dx,dy,dz);
                        printf(" \t   Analytical     %9.6lf%12.6lf%12.6lf%12.6lf\n",deriv_evec_cla[i][k][0][j],deriv_evec_cla[i][k][1][j],deriv_evec_cla[i][k][2][j],deriv_evec_cla[i][k][3][j]);
                        printf(" \t   Error (%%)      %9.6lf%s%10.6lf%s%10.6lf%s%10.6lf%s\n",(deriv_evec_cla[i][k][0][j]-ds)/deriv_evec_cla[i][k][0][j]*100.," %",                                 \
                                (deriv_evec_cla[i][k][1][j]-dx)/deriv_evec_cla[i][k][1][j]*100.," %",(deriv_evec_cla[i][k][2][j]-dy)/deriv_evec_cla[i][k][2][j]*100.," %",      \
                                (deriv_evec_cla[i][k][3][j]-dz)/deriv_evec_cla[i][k][3][j]*100.," %");
                        printf("\n");
                    }
                }
                gsl_vector_set( sp_sys->classic_ion[i].core_position, j, gsl_vector_get(sp_sys->classic_ion[i].core_position,j) - dl );
            }
        }

		// Free BFGS WorkSpace
		gsl_vector_free(v_tmp1);
		gsl_matrix_free(m_tmp1);    gsl_matrix_free(m_tmp3);
		gsl_matrix_free(m_tmp2);    gsl_matrix_free(m_tmp4);
		gsl_matrix_free(Inv_B_k);
		gsl_vector_free(p_k);
		gsl_vector_free(g_k_next);  gsl_vector_free(g_k_prev);
		gsl_vector_free(s_k);
		gsl_vector_free(x_k_next);  gsl_vector_free(x_k_prev);
		gsl_vector_free(y_k);
		MPI_Barrier(MPI_COMM_WORLD);
		ret = SP_SYSTEM_TRUE;
		return ret;
	}

	///	///	///	///	///	///	///	///	///	///	///	///	///	END SINGLE POINT
    if( rank == 0 )
    {   
        printf("\n");
        printf(" Initialising Gradient Descent Optimiser ... \n");
        printf("\n");
	printf(" Position Integral Reference : %.12lf\n",sp_cluster_integrator_get_x_12( &sp_sys->sp_ion[0] ));
	printf("\n");
        printf(" Internal Tolerances\n");	printf("\n");
        printf(" SCF Max Trials		   (-)		:	%d\n",MAX_SCF_CYCLE);
        printf(" SCF Eigenvalue  Tolerance (eV)		:	%.6e\n",SP_SYSTEM_SCF_TOL);
        printf(" SCF Eigenvector Tolerance (-)		:	%.6e\n",SP_SYSTEM_EVEC_TOL);
        printf(" Gnorm		 Tolerance (eV/Angs)	:	%.6e\n",SP_SYSTEM_GNORM_TOL);
    }
	///
	
    for(int bfgs_step = 0; bfgs_step < _max_bfgs_step_ ; bfgs_step++)
    {   
        /*
        gsl_matrix_set_zero(Inv_B_k);
        for(int i=0;i<BFGS_Stride;i++)  
            gsl_matrix_set(Inv_B_k,i,i,1.);
        */
        // if the above comment is removed, the optimiser will be carried out like a simple gradient descent method
    
        // Load x_k_prev (x_k_prev)
        sp_cluster_system_bfgs_support_load_x( sp_sys, x_k_prev );
        // Load g_k_prev (g_k)
        sp_cluster_system_bfgs_support_load_g( sp_sys, g_k_prev );
        // Get p_k
        gsl_blas_dgemv( CblasNoTrans, -1., Inv_B_k, g_k_prev, 0., p_k );   // This operation is ... p_k = -1.*Inv_B_k*g_k_prev + 0.*p_k ... remind if 'bfgs_step = 0' Inv_B_k = I;
        //MPI_Barrier(MPI_COMM_WORLD);

        ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     
        ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     
        if(rank == 0)
        {   liner();

            printf("\n");
            printf(" Configuration : %d \n",bfgs_step+1);
            printf("\n");
            printf(" Geometric Derivatives ( eV / Angstrom )\n");
            printf("\n");
            printf("-----------------------------------------------------------------------------\n");
            printf(" Species.        x           y           z          |r| \n");
            printf("-----------------------------------------------------------------------------\n");

            for(int i=0;i<number_of_classic_ion+number_of_sp_ion;i++)
            {   
                if( i< number_of_classic_ion )
                {   
                    grad_x = -sp_sys->classic_ion[i].elec_force_by_sp[0]-sp_sys->classic_ion[i].force_by_sp_core[0]-sp_sys->classic_ion[i].force_by_ion_core[0];
                    grad_y = -sp_sys->classic_ion[i].elec_force_by_sp[1]-sp_sys->classic_ion[i].force_by_sp_core[1]-sp_sys->classic_ion[i].force_by_ion_core[1];
                    grad_z = -sp_sys->classic_ion[i].elec_force_by_sp[2]-sp_sys->classic_ion[i].force_by_sp_core[2]-sp_sys->classic_ion[i].force_by_ion_core[2];
		    if( sp_sys->classic_ion[i].if_shell == SP_SYSTEM_TRUE )		// if the MM type is shell
			printf("%3s%3s%15.6lf%12.6lf%12.6lf%12.6lf\n",sp_sys->classic_ion[i].atom_name,"s",grad_x,grad_y,grad_z,sqrt(grad_x*grad_x+grad_y*grad_y+grad_z*grad_z));
		    else if( sp_sys->classic_ion[i].if_shell == SP_SYSTEM_FALSE )	// if the MM type is core
			printf("%3s%3s%15.6lf%12.6lf%12.6lf%12.6lf\n",sp_sys->classic_ion[i].atom_name,"c",grad_x,grad_y,grad_z,sqrt(grad_x*grad_x+grad_y*grad_y+grad_z*grad_z));
                }   // gradient on ion core
                else
                {   offset = i-number_of_classic_ion;
                    grad_x = -sp_sys->sp_ion[offset].elec_force_by_sp[0]-sp_sys->sp_ion[offset].elec_force_by_ion[0]-sp_sys->sp_ion[offset].force_by_sp_core[0]-sp_sys->sp_ion[offset].force_by_ion_core[0];
                    grad_y = -sp_sys->sp_ion[offset].elec_force_by_sp[1]-sp_sys->sp_ion[offset].elec_force_by_ion[1]-sp_sys->sp_ion[offset].force_by_sp_core[1]-sp_sys->sp_ion[offset].force_by_ion_core[1];
                    grad_z = -sp_sys->sp_ion[offset].elec_force_by_sp[2]-sp_sys->sp_ion[offset].elec_force_by_ion[2]-sp_sys->sp_ion[offset].force_by_sp_core[2]-sp_sys->sp_ion[offset].force_by_ion_core[2];
                    printf("%3s%18.6lf%12.6lf%12.6lf%12.6lf\n",sp_sys->sp_ion[offset].atom_name,grad_x,grad_y,grad_z,sqrt(grad_x*grad_x+grad_y*grad_y+grad_z*grad_z));
                }   // gradient on sp core
            }
            printf("-----------------------------------------------------------------------------\n");
            printf("\n");
            printf(" Cycle              :   %d\n",bfgs_step+1);
            //printf(" Gnorm  (eV/Angs)   :   %.6lf\n",gsl_blas_dnrm2(g_k_prev)/(double)BFGS_Stride);
            printf(" Gnorm  (eV/Angs)   :   %.6lf\n",gsl_blas_dnrm2(g_k_prev)/(double)deg_free);
            printf(" Energy (eV)        :  %.12lf\n",sp_cluster_system_get_cluster_energy(sp_sys));
            printf("\n");
            printf("-----------------------------------------------------------------------------\n");
            printf("\n");
            printf(" Lone Pair Molecular Orbital Info ( Lowest EigenValue / EigenVector Set )\n");
            printf("\n");
            printf("-----------------------------------------------------------------------------\n");
            printf(" Species.   Energy(eV)       s           px           py          pz \n");
            printf("-----------------------------------------------------------------------------\n");
            
            for(int kk=0;kk<number_of_sp_ion;kk++)
            {   int low_stat = sp_cluster_support_get_lowest_state(sp_sys->sp_ion[kk].eigen_value);
                printf("%3s%18.6lf%12.6lf%13.6lf%13.6lf%13.6lf\n",sp_sys->sp_ion[kk].atom_name,gsl_vector_get(sp_sys->sp_ion[kk].eigen_value,low_stat),
                        gsl_matrix_get(sp_sys->sp_ion[kk].eigen_vector,0,low_stat),gsl_matrix_get(sp_sys->sp_ion[kk].eigen_vector,1,low_stat),
                        gsl_matrix_get(sp_sys->sp_ion[kk].eigen_vector,2,low_stat),gsl_matrix_get(sp_sys->sp_ion[kk].eigen_vector,3,low_stat));
            }
            printf("-----------------------------------------------------------------------------\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);

	
	// ADD TEMPORAL RESULT  ... write xyz
	if( rank == 0 )
	{	   
        printf("\n");
        printf(" CONFIGURATION_XYZ_INFO ( step / number of atoms )   :   %d  /  %d\n",bfgs_step+1,sp_sys->number_of_classic_ion+sp_sys->number_of_sp_ion);
		sp_cluster_support_print_xyz( sp_sys, sp_cluster_system_get_cluster_energy(sp_sys), rank, numtasks );
        printf("\n");
        printf("-----------------------------------------------------------------------------\n");
	}
	


        ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     
        ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     

        // Find alpha ... Inexact LineSearch (loosely) by Backtracking (or Armijo linesearch)
        if( rank == 0 )
        {   printf("\n");
            printf(" Initiate Line Search ( Armijo-Goldstein-Wolfe )\n");
        }
	wtime = -MPI_Wtime();
        alpha = sp_cluster_system_bfgs_support_get_alpha_mpi_ls( sp_sys, p_k, stepmx, rank, numtasks );
	// AT THIS POINT THE FORCE IS ALREADY CALCULATED !
	wtime += MPI_Wtime();
        if( rank == 0 )
	{   printf(" Line Search Wtime	:  %.6lf s\n",wtime);	
	    printf("\n");
            printf("-----------------------------------------------------------------------------\n");
	}

        if( alpha == SP_SYSTEM_FALSE )  // if simple ls failed ... then try Damping
        {   
            if( rank == 0 )
                printf(" \nLINE SEARCH FAILED !!! TERMINATING THE PROGRAM ...\n");
            MPI_Barrier(MPI_COMM_WORLD);
            break;
        }
        else    // if bfgs ls succeeds
        {   
            // Load g_k_next (g_k+1)
            sp_cluster_system_bfgs_support_load_g( sp_sys, g_k_next );

            // POSITIVE DEFINITENESS CHECKER
            
            // Get s_k ... at this point you dont need p_k anymore
            gsl_vector_memcpy(s_k,p_k);
            gsl_vector_scale(s_k,alpha);            // SET : s_k = alpha*p_k
            // Get y_k = g_k_next - g_k_prev
            gsl_vector_memcpy(y_k,g_k_next);
            gsl_vector_sub(y_k,g_k_prev);           // SET : y_k = g_k+1 - g_k
            // s_k(T)*y_k
            gsl_blas_ddot(s_k,y_k,&s_tmp1);  

            //if( s_tmp1 > 0. && gsl_blas_dnrm2(g_k_prev)/(double)BFGS_Stride < SP_SYSTEM_GNORM_TOL /*10E-5*/ ) // if positive definiteness is satisfied
            if( s_tmp1 > 0. && gsl_blas_dnrm2(g_k_prev)/(double)deg_free < SP_SYSTEM_GNORM_TOL /*10E-5*/ ) // if positive definiteness is satisfied
            {
                if_bfgs = SP_SYSTEM_TRUE;
                if( rank == 0 )
                {   printf("\n");
                    printf(" Positive Definiteness ( TRUE )  :%12.6lf\n",s_tmp1);
                    printf("\n");
                }
            }
            else             // if positive deriniteness isn't satisfied
            {   
                if_bfgs = SP_SYSTEM_FALSE;
                if( rank == 0 )
                {   
                    printf("\n");
                    if( s_tmp1 > 0. )
                        printf(" Positive Definiteness ( TRUE )  :%12.6lf\n",s_tmp1);
                    else
                        printf(" Positive Definiteness ( FALSE ) :%12.6lf\n",s_tmp1);
		            //printf("B_Inv Matrix Reset -> I\n");                
                }
                // reset B_Inv Matrix into 'I'
                gsl_matrix_set_zero(Inv_B_k);
                for(int i=0;i<BFGS_Stride;i++)  
                    gsl_matrix_set(Inv_B_k,i,i,1.);
            }

        }
        // At this point already a step is taken
        // CHECK SUM_SQR
        sp_cluster_system_bfgs_support_load_x( sp_sys, x_k_next );
        gsl_vector_sub(x_k_prev,x_k_next);
        sum_sqr = gsl_blas_dnrm2(x_k_prev)/(double)BFGS_Stride;

        if( rank == 0 )
        {   printf("\n");
            printf(" Config RMS ( current_step - previous_step ) : %.6lf\n\n",sum_sqr);
            printf(" Finalising Optimisation Step ... \n");
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // CHECK SUM_SQR END
        // END OF ECKART & SUM_SQR TEST ///     ///     ///     ///     ///     ///     ///     ///     ///     ///     ///

    
        // if you wanna do a simple gradient descent ... use the line below
        //if_bfgs = SP_SYSTEM_FALSE;

        //else if( number_of_sp_ion > 1 && if_bfgs == SP_SYSTEM_TRUE )
        if( if_bfgs == SP_SYSTEM_TRUE )
        {   
            // reset if_bfgs flag
            if_bfgs = SP_SYSTEM_FALSE;

            // Get s_k ... at this point you dont need p_k anymore
            gsl_vector_memcpy(s_k,p_k);
            gsl_vector_scale(s_k,alpha);            // SET : s_k = alpha*p_k
            // Get y_k = g_k_next - g_k_prev
            gsl_vector_memcpy(y_k,g_k_next);
            gsl_vector_sub(y_k,g_k_prev);           // SET : y_k = g_k+1 - g_k

            // Now We Have ... InvB, y_k, and s_k
            
            // s_k(T)*y_k
            gsl_blas_ddot(s_k,y_k,&s_tmp1);  
            // y_k(T)*Inv_B_K*y_k
            gsl_blas_dgemv(CblasNoTrans,1.,Inv_B_k,y_k,0.,v_tmp1);
            gsl_blas_ddot(y_k,v_tmp1,&s_tmp2);
            // s_k*s_k(T)
            for(int i=0;i<BFGS_Stride;i++)
            {   for(int j=0;j<BFGS_Stride;j++)
                    gsl_matrix_set(m_tmp1,i,j,gsl_vector_get(s_k,i)*gsl_vector_get(s_k,j)*(s_tmp1+s_tmp2)/(s_tmp1*s_tmp1));
            }
            // 1st Matrix Calculation Done ... m_tmp1
            // s_tmp1 , m_tmp1 are in use
            
            // Inv_B_k*y_k*s_k(T)
            //gsl_blas_dgemv(CblasNoTrans,1.,Inv_B_k,y_k,0 .... is alrdy in v_tmp1
            for(int i=0;i<BFGS_Stride;i++)
            {   for(int j=0;j<BFGS_Stride;j++)
                    gsl_matrix_set(m_tmp2,i,j,gsl_vector_get(v_tmp1,i)*gsl_vector_get(s_k,j));
            }
            // s_k*y_k(T)*Inv_B_k
            for(int i=0;i<BFGS_Stride;i++)
            {   for(int j=0;j<BFGS_Stride;j++)
                    gsl_matrix_set(m_tmp3,i,j,gsl_vector_get(s_k,i)*gsl_vector_get(y_k,j));
            }
            gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,m_tmp3,Inv_B_k,0.,m_tmp4);

            gsl_matrix_set_zero(m_tmp3);
            // ( m_tmp2 + m_tmp4 )/s_tmp1
            gsl_matrix_add(m_tmp2,m_tmp4);  
            gsl_matrix_scale(m_tmp2,1./s_tmp1);
            gsl_matrix_set_zero(m_tmp4);

            // 2nd Matrix Calculation Done ... m_tmp2

            // GET Inv_B_k+1 = Inv_B_k + m_tmp1 - m_tmp2

            gsl_matrix_add(Inv_B_k,m_tmp1); // Inv_B_k = Inv_B_k + m_tmp1;
            gsl_matrix_sub(Inv_B_k,m_tmp2); // Inv_B_k = Inv_B_k - m_tmp2;
            
            MPI_Barrier(MPI_COMM_WORLD);
            // Inv_B_k+1 update Done
        }


        // CONVERGENCE CONDITION
        
        //gnorm_next = gsl_blas_dnrm2(g_k_next)/(double)BFGS_Stride;
        gnorm_next = gsl_blas_dnrm2(g_k_next)/(double)deg_free;
        e_next = sp_cluster_system_get_cluster_energy(sp_sys);

        //if( gsl_blas_dnrm2(g_k_next)/(double)BFGS_Stride < _gnorm_tol_ )
        if( gsl_blas_dnrm2(g_k_next)/(double)deg_free < _gnorm_tol_ )
        {	
	    ret = SP_SYSTEM_TRUE;
            // INITIALISING FIRST CONFIGURATION
            sp_cluster_system_set_is_scf_done(sp_sys,SP_SYSTEM_FALSE);
            MPI_Barrier(MPI_COMM_WORLD);
            /// CALCULATE TRIAL ENERGY & FORCE
            sp_cluster_system_scf_mpi(sp_sys,rank,numtasks);    // Get SCF Achieved sp set
            MPI_Barrier(MPI_COMM_WORLD);
            sp_cluster_system_get_force_mpi(sp_sys,rank,numtasks);
            sp_cluster_system_get_classic_force(sp_sys);
            MPI_Barrier(MPI_COMM_WORLD);
            
            sp_cluster_system_bfgs_support_load_g( sp_sys, g_k_next );
            //gnorm_next = gsl_blas_dnrm2(g_k_next)/(double)BFGS_Stride;
            gnorm_next = gsl_blas_dnrm2(g_k_next)/(double)deg_free;
            
            if( gnorm_next < _gnorm_tol_ )
            {
                if( rank == 0 ) 
                {   liner();
                    printf("\n");
                    printf(" Optimisation Meets Termination Condition, Final Configuration is\n");
                    printf("\n");
                    printf(" Geometric Derivatives ( eV / Angstrom )\n");
                    printf("\n");
                    printf("-----------------------------------------------------------------------------\n");
                    printf(" Species.        x           y           z          |r| \n");
                    printf("-----------------------------------------------------------------------------\n");


                    for(int i=0;i<number_of_classic_ion+number_of_sp_ion;i++)
                    {   
                        if( i< number_of_classic_ion )
                        {   
                            grad_x = -sp_sys->classic_ion[i].elec_force_by_sp[0]-sp_sys->classic_ion[i].force_by_sp_core[0]-sp_sys->classic_ion[i].force_by_ion_core[0];
                            grad_y = -sp_sys->classic_ion[i].elec_force_by_sp[1]-sp_sys->classic_ion[i].force_by_sp_core[1]-sp_sys->classic_ion[i].force_by_ion_core[1];
                            grad_z = -sp_sys->classic_ion[i].elec_force_by_sp[2]-sp_sys->classic_ion[i].force_by_sp_core[2]-sp_sys->classic_ion[i].force_by_ion_core[2];
			    if( sp_sys->classic_ion[i].if_shell == SP_SYSTEM_TRUE )		// if the MM type is shell
				printf("%3s%3s%15.6lf%12.6lf%12.6lf%12.6lf\n",sp_sys->classic_ion[i].atom_name,"s",grad_x,grad_y,grad_z,sqrt(grad_x*grad_x+grad_y*grad_y+grad_z*grad_z));
			    else if( sp_sys->classic_ion[i].if_shell == SP_SYSTEM_FALSE )	// if the MM type is core
				printf("%3s%3s%15.6lf%12.6lf%12.6lf%12.6lf\n",sp_sys->classic_ion[i].atom_name,"c",grad_x,grad_y,grad_z,sqrt(grad_x*grad_x+grad_y*grad_y+grad_z*grad_z));
                        }   // gradient on ion core
                        else
                        {   offset = i-number_of_classic_ion;
                            grad_x = -sp_sys->sp_ion[offset].elec_force_by_sp[0]-sp_sys->sp_ion[offset].elec_force_by_ion[0]-sp_sys->sp_ion[offset].force_by_sp_core[0]-sp_sys->sp_ion[offset].force_by_ion_core[0];
                            grad_y = -sp_sys->sp_ion[offset].elec_force_by_sp[1]-sp_sys->sp_ion[offset].elec_force_by_ion[1]-sp_sys->sp_ion[offset].force_by_sp_core[1]-sp_sys->sp_ion[offset].force_by_ion_core[1];
                            grad_z = -sp_sys->sp_ion[offset].elec_force_by_sp[2]-sp_sys->sp_ion[offset].elec_force_by_ion[2]-sp_sys->sp_ion[offset].force_by_sp_core[2]-sp_sys->sp_ion[offset].force_by_ion_core[2];
                            printf("%3s%18.6lf%12.6lf%12.6lf%12.6lf\n",sp_sys->sp_ion[offset].atom_name,grad_x,grad_y,grad_z,sqrt(grad_x*grad_x+grad_y*grad_y+grad_z*grad_z));
                        }   // gradient on sp core
                    }
                    printf("-----------------------------------------------------------------------------\n");
                    printf("\n");
                    printf(" Cycle              :   %d\n",bfgs_step+1);
                    //printf(" Gnorm  (eV/Angs)   :   %.6lf\n",gsl_blas_dnrm2(g_k_next)/(double)BFGS_Stride);
                    printf(" Gnorm  (eV/Angs)   :   %.6lf\n",gsl_blas_dnrm2(g_k_next)/(double)deg_free);
                    printf(" Energy (eV)        :  %.12lf\n",sp_cluster_system_get_cluster_energy(sp_sys));
                    printf("\n");
                    printf("-----------------------------------------------------------------------------\n");
                    printf("\n");
                    printf(" Lone Pair Molecular Orbital ( Lowest Eigenvalue / EigenVector )\n");
                    printf("\n");
                    printf("-----------------------------------------------------------------------------\n");
                    printf(" Species.   Energy(eV)       s           px           py          pz \n");
                    printf("-----------------------------------------------------------------------------\n");
                    
                    for(int kk=0;kk<number_of_sp_ion;kk++)
                    {   int low_stat = sp_cluster_support_get_lowest_state(sp_sys->sp_ion[kk].eigen_value);
                        printf("%3s%18.6lf%12.6lf%13.6lf%13.6lf%13.6lf\n",sp_sys->sp_ion[kk].atom_name,gsl_vector_get(sp_sys->sp_ion[kk].eigen_value,low_stat),
                                gsl_matrix_get(sp_sys->sp_ion[kk].eigen_vector,0,low_stat),gsl_matrix_get(sp_sys->sp_ion[kk].eigen_vector,1,low_stat),
                                gsl_matrix_get(sp_sys->sp_ion[kk].eigen_vector,2,low_stat),gsl_matrix_get(sp_sys->sp_ion[kk].eigen_vector,3,low_stat));
                    }
                    printf("-----------------------------------------------------------------------------\n");
                    /*
                    */
                }
                MPI_Barrier(MPI_COMM_WORLD);
                break;
            }
            else
            {   gsl_matrix_set_zero(Inv_B_k);
                for(int i=0;i<BFGS_Stride;i++)  gsl_matrix_set(Inv_B_k,i,i,1.);
            }
        }
        // BREAK;


    }

    // ADD TEMPORAL RESULT  ... write xyz
    if( rank == 0 )
    {	
        printf("\n");
        printf(" CONFIGURATION_XYZ_INFO ( final / number of atoms )   :   %d\n",sp_sys->number_of_classic_ion+sp_sys->number_of_sp_ion);
		sp_cluster_support_print_xyz( sp_sys, sp_cluster_system_get_cluster_energy(sp_sys), rank, numtasks );
        printf("\n");
        printf("-----------------------------------------------------------------------------\n");
    }


    // Free BFGS WorkSpace
    gsl_vector_free(v_tmp1);
    gsl_matrix_free(m_tmp1);    gsl_matrix_free(m_tmp3);
    gsl_matrix_free(m_tmp2);    gsl_matrix_free(m_tmp4);
    gsl_matrix_free(Inv_B_k);
    gsl_vector_free(p_k);
    gsl_vector_free(g_k_next);  gsl_vector_free(g_k_prev);
    gsl_vector_free(s_k);
    gsl_vector_free(x_k_next);  gsl_vector_free(x_k_prev);
    gsl_vector_free(y_k);
    MPI_Barrier(MPI_COMM_WORLD);

    return ret;
}   

