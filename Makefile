src     = calc2DIR.cpp
exes    = calc2DIR.exe
CC      = g++
LIBS    = -lm -lfftw3

all: ${exes}

${exes}: ${src} calc2DIR.h
	$(CC) $(src) -o $(exes) $(LIBS) -std=c++11 -fmax-errors=10 -O3

clean:
	rm calc2DIR.exe
