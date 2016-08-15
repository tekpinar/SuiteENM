#This is my makefile to produce libtekpinar.a
#If the library is changed, this make file will reproduce 
#the library.

#C Compiler
CC=gcc

#Fortran Compiler
FC=gfortran

OPTFLAGS=-O3 -funroll-loops

CFLAGS=-Wall ${OPTFLAGS}

#Library path
MYLIBDIR=/Users/mt/Documents/ccr_home/mylib
INC_DIR=${MYLIBDIR}/include
SRC_DIR=${MYLIBDIR}/src

OBJ=${MYLIBDIR}/enm_functions.o ${MYLIBDIR}/pair_functions.o ${MYLIBDIR}/rmsd_functions.o ${MYLIBDIR}/distance_functions.o ${MYLIBDIR}/aa_functions.o ${MYLIBDIR}/pdb_io.o ${MYLIBDIR}/form_factors.o ${MYLIBDIR}/stat_functions.o ${MYLIBDIR}/waxs_functions.o ${MYLIBDIR}/rc_functions.o ${MYLIBDIR}/blast_functions.o ${MYLIBDIR}/lapack_blas_functions.o ${MYLIBDIR}/aa_form_factors.o ${MYLIBDIR}/time_functions.o ${MYLIBDIR}/nd_functions.o 

#I removed waxs_functions.o from the upper OBJ variable since I didn't have gsl library and it was causing compilation problem. Added by MT on July 24, 2013.
#OBJ=${MYLIBDIR}/enm_functions.o ${MYLIBDIR}/pair_functions.o ${MYLIBDIR}/rmsd_functions.o ${MYLIBDIR}/distance_functions.o ${MYLIBDIR}/aa_functions.o ${MYLIBDIR}/pdb_io.o ${MYLIBDIR}/form_factors.o ${MYLIBDIR}/stat_functions.o  ${MYLIBDIR}/rc_functions.o ${MYLIBDIR}/blast_functions.o ${MYLIBDIR}/lapack_blas_functions.o ${MYLIBDIR}/aa_form_factors.o ${MYLIBDIR}/time_functions.o ${MYLIBDIR}/nd_functions.o 

default:
	(make -f pair_functions.mk)
	(make -f enm_functions.mk)
	(make -f rmsd_functions.mk)
	(make -f distance_functions.mk)
	(make -f aa_functions.mk)
	(make -f pdb_io.mk)
	(make -f stat_functions.mk)
	(make -f form_factors.mk)
	(make -f waxs_functions.mk)
	(make -f rc_functions.mk)
	(make -f blast_functions.mk)
	(make -f lapack_blas_functions.mk)
	(make -f aa_form_factors.mk)
	(make -f time_functions.mk)
	(make -f nd_functions.mk)
	(${FC} -c ${SRC_DIR}/u3b.f)
	(rm -f libtekpinar.a)
	(ar rcs libtekpinar.a ${OBJ} ${MYLIBDIR}/u3b.o)

#all: libtekpinar.a 

#libtekpinar.a: ${OBJ} ${MYLIBDIR}/u3b.o
#	rm -f $@
#	ar rcs $@  ${OBJ} ${MYLIBDIR}/u3b.o

#%.o: ./src/%.c
#	$(CC) $(CFLAGS) -I${INC_DIR} -c $<

#u3b.o: ${SRC_DIR}/u3b.f
#	g77 -c ${SRC_DIR}/u3b.f

clean: 
	rm -rf ${OBJ} ${MYLIBDIR}/u3b.o
