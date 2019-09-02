#This is Makefile for cpu_enm_nma program
CC=icc

CFLAGS= -Wall -pg

#Optimization flags
OPTFLAG=-O3 -funroll-loops #-mtune=generic

#Include files directory
INC_DIR=./include

#Source files directory
SRC_DIR=./src

#Library files directory
LIB_DIR=./lib

#Library files directory
BIN_DIR=./bin

OBJ=cpu_enm_nma.o

# cpu_enm_nma.exe: ${LIB_DIR}/${OBJ}
# 	${CC} ${CFLAGS} ${OBJ} -lm -framework Accelerate -o cpu_enm_nma.exe

# cpu_enm_nma.o: ${SRC_DIR}/cpu_enm_nma.c
# 	${CC}  ${CFLAGS} -c ${SRC_DIR}/cpu_enm_nma.c -o ${LIB_DIR}/cpu_enm_nma.o

MKLROOT=/opt/intel/mkl

cpu_enm_nma.exe: ${SRC_DIR}/cpu_enm_nma.c
	icc ${SRC_DIR}/cpu_enm_nma.c -O3 -axCORE-AVX2,AVX,SSE4.2 -o ${BIN_DIR}/cpu_enm_nma.exe -L ${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -Wl,-rpath=$GSL_LIBDIR -Wl,-rpath=${MKLROOT}/lib -fopenmp

clean:
	rm -rf ${LIB_DIR}/cpu_enm_nma.o cpu_enm_nma.exe
