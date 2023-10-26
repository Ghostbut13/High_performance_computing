
opt = O2
CC = gcc
CFLAG = -Wall -g -$(opt) -lm -lpthread

SRC=newton

clean:
	clear
	touch ./$(SRC)
	rm ./$(SRC)
	touch *.ppm
	rm ./*.ppm
	touch *.png
	rm ./*.png

$(SRC): $(SRC).c
	$(CC) $(CFLAG) $< -o $@ 

