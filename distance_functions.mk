#This is my Makefile for distance_functions.c
CC=gcc

CFLAGS=-Wall

#Include files directory
INC_DIR=./include

#Source files directory
SRC_DIR=./src

distance_functions.o: ${SRC_DIR}/distance_functions.c ${INC_DIR}/structures.h ${INC_DIR}/distance_functions.h
	${CC} ${CFLAGS} -c ${SRC_DIR}/distance_functions.c -I${INC_DIR}/ -o ${LIB_DIR}/distance_functions.o
