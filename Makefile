CFLAGS=-Wall -O3 -g


adabmDCA: adabmDCA.c
	gcc ${CFLAGS} adabmDCA.c -lm -o $@

clean:
	rm -f adabmDCA
