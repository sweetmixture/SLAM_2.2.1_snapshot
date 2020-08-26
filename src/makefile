#
# (c) Author:    Woongkyu Jee, woong.jee.16@ucl.ac.uk, wldndrb1@gmail.com
# Created:   02.06.2019 ~
# 	
# University College London, Department of Chemistry
#
MCC=mpiicc -std=c11 -mcmodel=large
CC=icc -std=c11 -mcmodel=large
LD=-lgsl -lgslcblas -lm
ILP=sp_cluster_integral_lib/



LIB= $(ILP)CH_Integrals.o $(ILP)CDH_Integrals.o $(ILP)CDDH_Integrals.o \
	 $(ILP)BM_SH_Integrals.o $(ILP)BM_SDH_Integrals.o $(ILP)CDDDH_Integrals.o $(ILP)BM_SDDH_Integrals.o \
	 $(ILP)BM_SDDDH_Integrals.o


#LIB= $(ILP)CH_Integrals.o $(ILP)SH_Integrals.o $(ILP)CDH_Integrals.o $(ILP)SDH_Integrals.o $(ILP)CDDH_Integrals.o \
#	 $(ILP)BM_SH_Integrals.o $(ILP)BM_SDH_Integrals.o $(ILP)CDDDH_Integrals.o $(ILP)BM_SDDH_Integrals.o \
#	 $(ILP)BM_SDDDH_Integrals.o


slam_v_2.2_opti.mpi.x: main_opti.o sp_cluster_system.o sp_cluster_support.o sp_cluster_integrator.o $(LIB)
	$(MCC) -o $@ $^ $(LD)

main_opti.o: main_opti.c
	$(MCC) -o $@ -c $^

main_fdm.o: main_fdm.c
	$(MCC) -o $@ -c $^

sp_cluster_system.o: sp_cluster_system.c sp_cluster_system.h sp_cluster_type.h sp_cluster_atom_type.h
	$(MCC) -o $@ -c sp_cluster_system.c

sp_cluster_support.o: sp_cluster_support.c sp_cluster_support.h
	$(CC) -o $@ -c sp_cluster_support.c


# Compile Integral Library

CH_Integrals.o: sp_cluster_integral_lib/CH_Integrals.c sp_cluster_integral_lib/CH_Integrals.h
	$(CC) -o $@ -c sp_cluster_integral_lib/CH_Integrals.c

#SH_Integrals.o: sp_cluster_integral_lib/SH_Integrals.c sp_cluster_integral_lib/SH_Integrals.h
#	$(CC) -o $@ -c sp_cluster_integral_lib/SH_Integrals.c

BM_SH_Integrals.o: sp_cluster_integral_lib/BM_SH_Integrals.c sp_cluster_integral_lib/BM_SH_Integrals.h
	$(CC) -o $@ -c sp_cluster_integral_lib/BM_SH_Integrals.c

CDH_Integrals.o: sp_cluster_integral_lib/CDH_Integrals.c sp_cluster_integral_lib/CDH_Integrals.h
	$(CC) -o $@ -c sp_cluster_integral_lib/CDH_Integrals.c

#SDH_Integrals.o: sp_cluster_integral_lib/SDH_Integrals.c sp_cluster_integral_lib/SDH_Integrals.h
#	$(CC) -o $@ -c sp_cluster_integral_lib/SDH_Integrals.c

BM_SDH_Integrals.o: sp_cluster_integral_lib/BM_SDH_Integrals.c sp_cluster_integral_lib/BM_SDH_Integrals.h
	$(CC) -o $@ -c sp_cluster_integral_lib/BM_SDH_Integrals.c

CDDH_Integrals.o: sp_cluster_integral_lib/CDDH_Integrals.c sp_cluster_integral_lib/CDDH_Integrals.h
	$(CC) -o $@ -c sp_cluster_integral_lib/CDDH_Integrals.c

BM_SDDH_Integrals.o: sp_cluster_integral_lib/BM_SDDH_Integrals.c sp_cluster_integral_lib/BM_SDDH_Integrals.h
	$(CC) -o $@ -c sp_cluster_integral_lib/BM_SDDH_Integrals.c

CDDDH_Integrals.o: sp_cluster_Integral_lib/CDDDH_Integrals.c sp_cluster_integral_lib/CDDDH_Integrals.h
	$(CC) -o $@ -c sp_cluster_integral_lib/CDDDH_Integrals.c

BM_SDDDH_Integrals.o: sp_cluster_integral_lib/BM_SDDDH_Integrals.c sp_cluster_integral_lib/BM_SDDDH_Integrals.h
	$(CC) -o $@ -c sp_cluster_integral_lib/BM_SDDDH_Integrals.c

sp_cluster_integrator.o: sp_cluster_integrator.c sp_cluster_integrator.h sp_cluster_type.h
	$(CC) -o $@ -c sp_cluster_integrator.c



clean:
	rm *.o
