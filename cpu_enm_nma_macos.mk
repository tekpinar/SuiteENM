#This is Makefile for cpu_enm_nma program

CC=gcc
CFLAGS= -Wall -pg

OBJ=cpu_enm_nma.o

cpu_enm_nma.exe: ${OBJ}
	${CC} ${CFLAGS} ${OBJ} -lm -framework Accelerate -o cpu_enm_nma.exe

cpu_enm_nma.o: cpu_enm_nma.c
	${CC}  ${CFLAGS} -c cpu_enm_nma.c -o cpu_enm_nma.o

clean:
	rm -rf cpu_enm_nma.o cpu_enm_nma.exe