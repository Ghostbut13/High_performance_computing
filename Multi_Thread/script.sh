#!/bin/bash


touch newton.tar.gz
rm newton.tar.gz
rm -r extracted
rm -r pictures


tar -czf newton.tar.gz newton.c makefile
/home/hpc2023/newton_iteration/check_submission.jl ./newton.tar.gz


