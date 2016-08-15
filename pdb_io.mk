#This is my Makefile for pdb_io.c
CC=gcc

CFLAG=-Wall

#Include files directory
INC_DIR=./include

#Source files directory
SRC_DIR=./src

pdb_io.o: ${SRC_DIR}/pdb_io.c ${INC_DIR}/structures.h ${INC_DIR}/rmsd_functions.h ${INC_DIR}/aa_functions.h ${INC_DIR}/pdb_io.h
	${CC}  ${CFLAG} -c -I${INC_DIR} ${SRC_DIR}/pdb_io.c -o pdb_io.o

