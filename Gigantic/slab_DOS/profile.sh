#!/bin/bash
# clang++-3.8 PDOS.cpp -O3 -std=c++14 -DMKL_LP64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -lprofiler && CPUPROFILE=foo.prof ./a.out
g++ -g ./PDOS.cpp -O3 -std=c++11 -DMKL_LP64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -lprofiler
./a.out
go tool pprof --text ./a.out profile.log

# go tool pprof ./a.out foo.prof 

