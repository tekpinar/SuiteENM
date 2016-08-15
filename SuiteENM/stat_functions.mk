#This is my Makefile for stat_functions.c
CC=gcc

CFLAG=-Wall

#Include files directory
INC_DIR=./include

#Source files directory
SRC_DIR=./src

stat_functions.o: ${SRC_DIR}/stat_functions.c ${INC_DIR}/stat_functions.h
	${CC}  ${CFLAG} -c  -I${INC_DIR} ${SRC_DIR}/stat_functions.c -o stat_functions.o