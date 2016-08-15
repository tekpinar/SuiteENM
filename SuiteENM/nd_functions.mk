#This is my Makefile for nd_functions.c
CC=gcc

CFLAG=-Wall

#Optimization flags
OPTFLAG=-O3 -funroll-loops #-mtune=generic

#Include files directory
INC_DIR=./include

#Source files directory
SRC_DIR=./src

nd_functions.o: ${SRC_DIR}/nd_functions.c ${INC_DIR}/nd_functions.h
	${CC}  ${CFLAG} ${OPTFLAG} -c -I${INC_DIR} ${SRC_DIR}/nd_functions.c -o nd_functions.o