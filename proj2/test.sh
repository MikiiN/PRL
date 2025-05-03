#!/bin/bash

if [ $# -ne 1 ]; then
  exit 1;
fi;


# check input tree len
length=${#1}
if [ $length -lt 1 ]; then
  >&2 echo "empty input tree"
  exit 1; 
fi;


mpic++ --prefix /usr/local/share/OpenMPI -o vuv vuv.cpp

proc_number=$(((2*length)-2))

# check if number of processes is not 0 (input tree has single node)
if [ $proc_number == 0 ]; then
  mpirun --oversubscribe -np 1 ./vuv $1
else
  mpirun --oversubscribe -np $proc_number ./vuv $1 
fi;

# clean up
rm -f vuv