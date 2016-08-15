#This is my Makefile for form_factors.c
CC=gcc

CFLAG=-Wall

#Optimization flags
OPTFLAG=-O3 -funroll-loops #-mtune=generic

#Include files directory
INC_DIR=./include

#Source files directory
SRC_DIR=./src

form_factors.o: ${SRC_DIR}/form_factors.c ${INC_DIR}/structures.h ${INC_DIR}/pdb_io.h ${INC_DIR}/form_factors.h
	${CC}  ${CFLAG} ${OPTFLAG} -c  -I${INC_DIR}/ ${SRC_DIR}/form_factors.c -o form_factors.o

