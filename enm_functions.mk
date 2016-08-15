#This is my Makefile for enm_functions.c
CC=gcc

CFLAG=-Wall

#Optimization flags
OPTFLAG=-O3 -funroll-loops #-mtune=generic

#Include files directory
INC_DIR=./include

#Source files directory
SRC_DIR=./src

enm_functions.o: ${SRC_DIR}/enm_functions.c ${INC_DIR}/structures.h ${INC_DIR}/distance_functions.h ${INC_DIR}/enm_functions.h
	${CC}  ${CFLAG} ${OPTFLAG} -c -I${INC_DIR} ${SRC_DIR}/enm_functions.c -o enm_functions.o