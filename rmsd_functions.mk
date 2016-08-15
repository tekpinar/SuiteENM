#This is my Makefile for rmsd_functions.c
CC=gcc

CFLAG=-Wall

#Optimization flags
OPTFLAG=-O3 -funroll-loops #-mtune=generic

#Include files directory
INC_DIR=./include

#Source files directory
SRC_DIR=./src

rmsd_functions.o: ${SRC_DIR}/rmsd_functions.c ${INC_DIR}/defines.h ${INC_DIR}/structures.h ${INC_DIR}/rmsd_functions.h
	${CC}  ${CFLAG} ${OPTFLAG} -c -I${INC_DIR} ${SRC_DIR}/rmsd_functions.c -o rmsd_functions.o