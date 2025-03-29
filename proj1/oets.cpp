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

bool has_left_neighbour(int proc){
    return proc > 0;
}

bool has_right_neighbour(int proc, int count){
    return proc < (count-1); 
}

void load_and_send_numbers(){
    char loadedNum; 
    int actualprocNum = 0;

    ifstream fileNum(INPUT);
    if(!fileNum){
        cerr << "Err: Can't open file numbers" << endl;
        MPI_Abort(MPI_COMM_WORLD, ERROR);
    }

    // load all numbers, print them on STDOUT and send them to processes
    while(fileNum.get(loadedNum)){
        int num = (unsigned char) loadedNum;
        cout << num << " ";
        MPI_Send(&num, 1, MPI_INT, actualprocNum, TAG, MPI_COMM_WORLD);
        actualprocNum++;
    }
    cout << endl;

    fileNum.close();
}

int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);
    int procID, processesNum;

    MPI_Comm_size(MPI_COMM_WORLD, &processesNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procID);

    // distribute numbers to processes
    if(procID == MASTER_ID){
        load_and_send_numbers();
    }

    MPI_Status status;
    int myNum;
    // process receive it's number
    MPI_Recv(&myNum, 1, MPI_INT, MASTER_ID, TAG, MPI_COMM_WORLD, &status);

    // sorting
    int numCycles = processesNum/2 + processesNum%2;
    for(int i = 0; i < numCycles; i++){
        // even processes
        if(is_even(procID) && has_right_neighbour(procID, processesNum)){
            // send its number and waiting until neighbour process send number back
            MPI_Send(&myNum, 1, MPI_INT, procID + 1, TAG, MPI_COMM_WORLD);
            MPI_Recv(&myNum, 1, MPI_INT, procID + 1, TAG, MPI_COMM_WORLD, &status);
        }
        // odd processes
        else if(is_odd(procID)){
            int neighbourNum;
            // receive number from even neighbour
            MPI_Recv(&neighbourNum, 1, MPI_INT, procID - 1, TAG, MPI_COMM_WORLD, &status);
            // compare with own number
            if(myNum < neighbourNum){
                // process's number is smaller, so process keeps neighbour number
                // and send its number back (swap)
                MPI_Send(&myNum, 1, MPI_INT, procID - 1, TAG, MPI_COMM_WORLD);
                myNum = neighbourNum;
            }
            else{
                // process's number is bigger or equal (no swap)
                MPI_Send(&neighbourNum, 1, MPI_INT, procID - 1, TAG, MPI_COMM_WORLD);
            }
        }

        // same like before, just now odd processes send and even processes compare
        if(is_odd(procID) && has_right_neighbour(procID, processesNum)){
            MPI_Send(&myNum, 1, MPI_INT, procID + 1, TAG, MPI_COMM_WORLD);
            MPI_Recv(&myNum, 1, MPI_INT, procID + 1, TAG, MPI_COMM_WORLD, &status);
        }
        else if(is_even(procID) && has_left_neighbour(procID)){
            int neighbourNum;
            MPI_Recv(&neighbourNum, 1, MPI_INT, procID - 1, TAG, MPI_COMM_WORLD, &status);
            if(myNum < neighbourNum){
                MPI_Send(&myNum, 1, MPI_INT, procID - 1, TAG, MPI_COMM_WORLD);
                myNum = neighbourNum;
            }
            else{
                MPI_Send(&neighbourNum, 1, MPI_INT, procID - 1, TAG, MPI_COMM_WORLD);
            }
        }
    } 

    // master process gather all processes numbers
    int result[processesNum];
    for(int id = 1; id < processesNum; id++){
        if(procID == id){
            MPI_Send(&myNum, 1, MPI_INT, MASTER_ID, TAG, MPI_COMM_WORLD);
        }
        if(procID == MASTER_ID){
            MPI_Recv(&(result[id]), 1, MPI_INT, id, TAG, MPI_COMM_WORLD, &status);
        }
    }

    // master process print sorted numbers on STDOUT
    if(procID == MASTER_ID){
        result[MASTER_ID] = myNum;
        for(int i = 0; i < processesNum; i++){
            cout << result[i] << endl;
        }
    }

    MPI_Finalize();
}