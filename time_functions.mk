#This is my Makefile for time_functions.c
CC=gcc

CFLAG=-Wall

#Optimization flags
OPTFLAG=-O3 -funroll-loops #-mtune=generic

#Include files directory
INC_DIR=./include

#Source files directory
SRC_DIR=./src

time_functions.o: ${SRC_DIR}/time_functions.c ${INC_DIR}/time_functions.h
	${CC}  ${CFLAG} ${OPTFLAG} -c -I${INC_DIR} ${SRC_DIR}/time_functions.c -o time_functions.o