CFLAGS=-Wall -O3 -g 
CC=g++
SRC= adabmDCA.cpp BMaux.cpp BMmc.cpp BMlib.cpp
DEPS= 
OBJ=$(SRC:.cpp = -o)

adabmDCA: $(OBJ) $(DEPS)
	$(CC) ${CFLAGS} -o adabmDCA $(OBJ)
clean:
	rm -f adabmDCA *.o
