#This is my Makefile for enm_functions.c
CC=gcc

CFLAG=-Wall

#Optimization flags
OPTFLAG=-O3 -funroll-loops #-mtune=generic

#Include files directory
INC_DIR=./include 

#Source files directory
SRC_DIR=./src

#Library files directory                                                           |
LIB_DIR=./lib

waxs_functions.o: ${SRC_DIR}/waxs_functions.c
	${CC}  ${CFLAG} ${OPTFLAG} -c -I/opt/local/include -I${INC_DIR} ${SRC_DIR}/waxs_functions.c -o ${LIB_DIR}/waxs_functions.o
