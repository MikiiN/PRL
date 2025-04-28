#!/bin/bash

if [ $# -ne 1 ]; then
  exit 1;
fi;

# preklad
mpic++ --prefix /usr/local/share/OpenMPI -o vuv vuv.cpp

length=${#1}
proc_number=$(((2*length)-2))
# spusteni aplikace (oversubscribe - vice procesu nez fyzicky k dispozici)
mpirun --oversubscribe -np $proc_number ./vuv $1 

# uklid
rm -f vuv