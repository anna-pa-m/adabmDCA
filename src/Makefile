CFLAGS=-Wall -O3 -g -fopenmp
CC=g++
SRC= adabmDCA.cpp BMmc.cpp BMaux.cpp BMlib.cpp
OBJ=$(SRC:.cpp = -o)



adabmDCA: $(OBJ)
	$(CC) ${CFLAGS} -o adabmDCA $(OBJ)
clean:
	rm -f adabmDCA *.o
