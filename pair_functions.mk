# Makefile for pair_functions.c  
#Set compiler
CC=gcc

#Warning flags
CFLAGS=-Wall

#Debugging flag 
DBGFLAG=-g

#Profiling flag
PRFLFLAG=-pg

#Optimization flags
OPTFLAG=-O3 -funroll-loops #-mtune=generic

#Include files directory
INC_DIR=./include

#Source files directory
SRC_DIR=./src

#Library files directory
LIB_DIR=./lib

pair_functions.o: ${SRC_DIR}/pair_functions.c ${INC_DIR}/defines.h ${INC_DIR}/structures.h ${INC_DIR}/distance_functions.h ${INC_DIR}/pair_functions.h
	${CC} ${CFLAGS} -c -std=c99  ${SRC_DIR}/pair_functions.c -I${INC_DIR} -o ${LIB_DIR}/pair_functions.o
