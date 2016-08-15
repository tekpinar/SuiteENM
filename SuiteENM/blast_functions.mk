#This is my Makefile for blast_functions.c

CC=gcc

CFLAG=-Wall

#Include files directory
INC_DIR=./include

#Source files directory
SRC_DIR=./src

blast_functions.o: ${SRC_DIR}/blast_functions.c ${INC_DIR}/defines.h ${INC_DIR}/structures.h ${INC_DIR}/blast_functions.h
	${CC}  ${CFLAG} -c -I${INC_DIR} ${SRC_DIR}/blast_functions.c -o blast_functions.o