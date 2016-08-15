#This is my Makefile for enm_functions.c
CC=gcc

CFLAG=-Wall

#Include files directory
INC_DIR=./include

#Source files directory
SRC_DIR=./src

aa_functions.o: ${SRC_DIR}/aa_functions.c ${INC_DIR}/aa_functions.h
	${CC}  ${CFLAG} -c  -I${INC_DIR} ${SRC_DIR}/aa_functions.c -o aa_functions.o