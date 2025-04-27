/*
* Project: PRL-proj2
* Author:  Michal Žatečka xzatec02
* Date:    27.04.2025
*/

#include <mpi.h>
#include <iostream>
#include <string>

using namespace std;

typedef struct node {
    char name;
    node *left;
    node *right;
} BinNode;


int main(int argc, char *argv[]){
    string input_tree = argv[1];
    cout << input_tree << endl;
}