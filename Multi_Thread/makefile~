#SOURCE_PATH =  /home/hpcuser053/work/High_performance_computing/OpenMP/src
#OUTPUT_PATH =  /home/hpcuser053/work/High_performance_computing/OpenMP/object

#-march=native

opt = O2
CC = gcc
CFLAG = -Wall -g -$(opt) -fopenmp -lm

SRC=cell_distances

clean:
	clear
	touch ./$(SRC)
	rm ./$(SRC)

$(SRC): $(SRC).c
	$(CC) $(CFLAG) $< -o $@
	#mv $@ $(OUTPUT_PATH)/$@
