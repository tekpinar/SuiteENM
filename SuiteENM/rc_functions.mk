#This is my Makefile for rc_functions.c
CC=gcc

CFLAG=-Wall

#Include files directory
INC_DIR=./include

#Source files directory
SRC_DIR=./src

rc_functions.o: ${SRC_DIR}/rc_functions.c ${INC_DIR}/structures.h ${INC_DIR}/pdb_io.h ${INC_DIR}/rmsd_functions.h ${INC_DIR}/rc_functions.h
	${CC}  ${CFLAG} -c -I${INC_DIR} ${SRC_DIR}/rc_functions.c -o rc_functions.o
