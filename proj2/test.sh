#!/bin/bash

if [ $# -ne 1 ]; then
  exit 1;
fi;


length=${#1}
if [ $length -lt 1 ]; then
  exit 1; 
fi;

# preklad
mpic++ --prefix /usr/local/share/OpenMPI -o vuv vuv.cpp

proc_number=$(((2*length)-2))
if [ $proc_number == 0 ]; then
  mpirun --oversubscribe -np 1 ./vuv $1
else
  mpirun --oversubscribe -np $proc_number ./vuv $1 
fi;

# uklid
rm -f vuv