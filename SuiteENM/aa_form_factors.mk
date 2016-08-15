#This is my Makefile for aa_form_factors.c
CC=gcc

CFLAG=-Wall

#Optimization flags
OPTFLAG=-O3 -funroll-loops #-mtune=generic

#Include files directory
INC_DIR=./include

#Source files directory
SRC_DIR=./src

aa_form_factors.o: ${SRC_DIR}/aa_form_factors.c ${INC_DIR}/structures.h ${INC_DIR}/form_factors.h ${INC_DIR}/distance_functions.h ${INC_DIR}/aa_form_factors.h
	${CC}  ${CFLAG} ${OPTFLAG} -c -I${INC_DIR} ${SRC_DIR}/aa_form_factors.c -o aa_form_factors.o