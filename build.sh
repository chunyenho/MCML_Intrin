#!/bin/bash

icc -no-vec -fopenmp -O3 -vec-report3 -w -I/opt/intel/composer_xe_2013_sp1.2.144/mkl/include -Wl,--start-group /opt/intel/composer_xe_2013_sp1.2.144/mkl/lib/intel64/libmkl_intel_lp64.so /opt/intel/composer_xe_2013_sp1.2.144/mkl/lib/intel64/libmkl_intel_thread.so /opt/intel/composer_xe_2013_sp1.2.144/mkl/lib/intel64/libmkl_core.so -Wl,--end-group -L/opt/intel/composer_xe_2013_sp1.2.144/mkl/../compiler/lib/intel64 -liomp5 -lm -ldl -lpthread mcmlmain.c mcmlgo.c mcmlio.c mcmlnr.c -fopenmp -I. -lm -o mcml.out 
