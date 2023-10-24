
opt = O2
CC = gcc
#CFLAG = -Wall -g -$(opt) -fopenmp -lm

SRC=new

clean:
	clear
	touch ./$(SRC)
	rm ./$(SRC)

$(SRC): $(SRC).c
	$(CC) $< -o $@

