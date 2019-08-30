#This is Makefile for cpu_enm_nma program
CC=gcc

CFLAGS= -Wall -pg

#Optimization flags
OPTFLAG=-O3 -funroll-loops #-mtune=generic

#Include files directory
INC_DIR=./include

#Source files directory
SRC_DIR=./src

#Library files directory
LIB_DIR=./lib

OBJ=cpu_enm_nma.o

cpu_enm_nma.exe: ${LIB_DIR}/${OBJ}
	${CC} ${CFLAGS} ${OBJ} -lm -framework Accelerate -o cpu_enm_nma.exe

cpu_enm_nma.o: ${SRC_DIR}/cpu_enm_nma.c
	${CC}  ${CFLAGS} -c ${SRC_DIR}/cpu_enm_nma.c -o ${LIB_DIR}/cpu_enm_nma.o

clean:
	rm -rf ${LIB_DIR}/cpu_enm_nma.o ${SRC_DIR}/cpu_enm_nma.exe