SOURCE_PATH =  ./src/
OUTPUT_PATH =  ./object/

#-march=native

opt = O2
CC = gcc
CFLAG = -Wall -g -$(opt) -fopenmp -lm

SRC=distance_omp

clean:
	clear

$(SRC).o: $(SOURCE_PATH)/$(SRC).c
	$(CC) $(CFLAG) $< -o $@
	mv $@ $(OUTPUT_PATH)/$@
