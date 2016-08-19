#This is my Makefile for blast_functions.c

CC=gcc

CFLAG=-Wall

#Include files directory
INC_DIR=./include

#Source files directory
SRC_DIR=./src

#Library files directory
LIB_DIR=./lib

blast_functions.o: ${SRC_DIR}/blast_functions.c ${INC_DIR}/defines.h ${INC_DIR}/structures.h ${INC_DIR}/blast_functions.h
	${CC}  ${CFLAG} -c -I${INC_DIR} ${SRC_DIR}/blast_functions.c -o ${LIB_DIR}/blast_functions.o
