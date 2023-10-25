#!/bin/bash


tar -czf newton.tar.gz newton.c makefile

rm -r extracted
rm -r pictures

/home/hpc2023/newton_iteration/check_submission.jl ./newton.tar.gz

