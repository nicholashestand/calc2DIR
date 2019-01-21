# Makefile

CC      = gcc

CC_FLAGS= -lm -O3

all: stochastic

stochastic: stochastic.o randomlib.o
	$(CC) randomlib.o stochastic.o -o stochastic $(CC_FLAGS) 

stochastic.o: stochastic.c stochastic.h

randomlib.o: randomlib.h randomlib.c
