#This is my Makefile for lapack_blas_functions.c
CC=gcc

CFLAG=-Wall

#Optimization flags
OPTFLAG=-O3 -funroll-loops #-mtune=generic

#Include files directory
INC_DIR=./include

#Source files directory
SRC_DIR=./src

lapack_blas_functions.o: ${SRC_DIR}/lapack_blas_functions.c ${INC_DIR}/lapack_blas_functions.h
	${CC}  ${CFLAG} ${OPTFLAG} -c -I${INC_DIR} ${SRC_DIR}/lapack_blas_functions.c -o lapack_blas_functions.o