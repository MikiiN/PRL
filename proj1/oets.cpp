#include <mpi.h>
#include <iostream>
#include <fstream>

#define INPUT "numbers"
#define ERROR 1
#define MASTER_ID 0
#define TAG 0

using namespace std;

bool is_odd(int num){
    return num % 2 == 1;
}

bool is_even(int num){
    return num % 2 == 0;
}

int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);

    int procID, processesNum;
    MPI_Comm_size(MPI_COMM_WORLD, &processesNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procID);

    if(procID == MASTER_ID){
        char loadedNum; 
        int actualprocNum = 0;

        ifstream fileNum(INPUT);
        if(!fileNum){
            cerr << "Err: Can't open file numbers" << endl;
            MPI_Abort(MPI_COMM_WORLD, ERROR);
        }

        while(fileNum.get(loadedNum)){
            int num = (unsigned char) loadedNum;
            cout << num << " ";
            MPI_Send(&num, 1, MPI_INT, actualprocNum, TAG, MPI_COMM_WORLD);
            actualprocNum++;
        }
        cout << endl;
        fileNum.close();
    }

    MPI_Status stat;
    int myNum;
    MPI_Recv(&myNum, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD, &stat);

    // TODO deadlock
    int numCycles = processesNum/2 + processesNum%2;
    for(int i = 0; i < numCycles; i++){
        if(is_even(procID) && procID < (processesNum-1)){
            MPI_Send(&myNum, 1, MPI_INT, procID + 1, TAG, MPI_COMM_WORLD);
            MPI_Recv(&myNum, 1, MPI_INT, procID + 1, TAG, MPI_COMM_WORLD, &stat);
        }
        else if(is_odd(procID)){
            int neighbourNum;
            MPI_Recv(&neighbourNum, 1, MPI_INT, procID - 1, TAG, MPI_COMM_WORLD, &stat);
            if(myNum < neighbourNum){
                MPI_Send(&myNum, 1, MPI_INT, procID - 1, TAG, MPI_COMM_WORLD);
                myNum = neighbourNum;
            }
            else{
                MPI_Send(&neighbourNum, 1, MPI_INT, procID - 1, TAG, MPI_COMM_WORLD);
            }
        }

        if(is_odd(procID) && procID < (processesNum-1)){
            MPI_Send(&myNum, 1, MPI_INT, procID + 1, TAG, MPI_COMM_WORLD);
            MPI_Recv(&myNum, 1, MPI_INT, procID + 1, TAG, MPI_COMM_WORLD, &stat);
        }
        else if(is_even(procID)){
            int neighbourNum;
            MPI_Recv(&neighbourNum, 1, MPI_INT, procID - 1, TAG, MPI_COMM_WORLD, &stat);
            if(myNum < neighbourNum){
                MPI_Send(&myNum, 1, MPI_INT, procID - 1, TAG, MPI_COMM_WORLD);
                myNum = neighbourNum;
            }
            else{
                MPI_Send(&neighbourNum, 1, MPI_INT, procID - 1, TAG, MPI_COMM_WORLD);
            }
        }
    } 

    int result[processesNum];
    for(int id = 1; id < processesNum; id++){
        if(procID == id){
            MPI_Send(&myNum, 1, MPI_INT, MASTER_ID, TAG, MPI_COMM_WORLD);
        }
        if(procID == MASTER_ID){
            MPI_Recv(&(result[id]), 1, MPI_INT, id, TAG, MPI_COMM_WORLD, &stat);
        }
    }

    if(procID == MASTER_ID){
        result[MASTER_ID] = myNum;
        for(int i = 0; i < processesNum; i++){
            cout << result[i] << endl;
        }
    }

    MPI_Finalize();
}