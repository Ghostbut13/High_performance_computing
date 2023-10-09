#!/bin/bash

touch /home/hpcuser053/work/High_performance_computing/OpenMP/extracted
rm -r /home/hpcuser053/work/High_performance_computing/OpenMP/extracted
tar -czf cell_distances.tar.gz ./makefile ./cell_distances.c
/home/hpc2023/cell_distances/check_submission.jl /home/hpcuser053/work/High_performance_computing/OpenMP/cell_distances.tar.gz
