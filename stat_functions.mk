#This is my Makefile for stat_functions.c
CC=gcc

CFLAG=-Wall

#Include files directory
INC_DIR=./include

#Source files directory
SRC_DIR=./src

#Library files directory                                                           |
LIB_DIR=./lib

stat_functions.o: ${SRC_DIR}/stat_functions.c ${INC_DIR}/stat_functions.h
	${CC}  ${CFLAG} -c  -I${INC_DIR} ${SRC_DIR}/stat_functions.c -o ${LIB_DIR}/stat_functions.o
