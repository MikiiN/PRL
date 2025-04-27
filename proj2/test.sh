#!/bin/bash

if [ $# -ne 1 ]; then
  exit 1;
fi;

# preklad
mpic++ --prefix /usr/local/share/OpenMPI -o vuv vuv.cpp


# spusteni aplikace (oversubscribe - vice procesu nez fyzicky k dispozici)
mpirun --oversubscribe -np 1 ./vuv $1 

# uklid
rm -f vuv