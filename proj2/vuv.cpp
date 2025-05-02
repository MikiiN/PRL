/*
* Project: PRL-proj2
* Author:  Michal Žatečka xzatec02
* Date:    27.04.2025
*/

#include <bits/stdc++.h>
#include <cstddef>
#include <mpi.h>
#include <iostream>
#include <string>
#include <tuple>

using namespace std;

#define BASE 2
#define FORWARD -1
#define MASTER_ID 0
#define NOT_EXIST -1
#define REVERSE 1
#define ROOT_INDEX 0
#define TAG 0


typedef struct g_edge {
    int start;
    int end;
}Edge; 

typedef struct adj_node {
    Edge forward;
    Edge reverse;
} AdjacencyListNode;

tuple<size_t,size_t> getNodeChildrenIndexes(size_t node_index, size_t tree_len){
    int firstChild = 2*node_index+1;
    int secondChild = firstChild+1;

    if(firstChild >= tree_len){
        firstChild = NOT_EXIST;
        secondChild = NOT_EXIST;
    }
    else if(secondChild >= tree_len)
        secondChild = NOT_EXIST;
    return make_tuple(firstChild, secondChild);
}

size_t getNodeParent(size_t node_index){
    if(node_index == 0)
        return NOT_EXIST;
    return (node_index-1) / 2;
}

AdjacencyListNode createNode(int node0, int node1){
    Edge edge = {node0, node1};
    Edge revert = {node1, node0};
    AdjacencyListNode node = {edge, revert};
    return node;
}

vector<AdjacencyListNode> getListForNode(size_t node_index, string tree){
    vector<AdjacencyListNode> result;
    size_t left_child, right_child;
    size_t tree_len = tree.length();
    size_t parent = getNodeParent(node_index);
    tie(left_child, right_child) = getNodeChildrenIndexes(node_index, tree_len);

    if(parent != NOT_EXIST){
        AdjacencyListNode p_node = createNode(node_index, parent);
        result.push_back(p_node);
    }
    if(left_child != NOT_EXIST){
        AdjacencyListNode l_node = createNode(node_index, left_child);
        result.push_back(l_node);
    }
    if(right_child != NOT_EXIST){
        AdjacencyListNode r_node = createNode(node_index, right_child);
        result.push_back(r_node);
    }
    return result;
}

tuple<MPI_Datatype,MPI_Datatype> getCustomTypes(){
    const int nitems=2;
    int          blocklengths[2] = {1,1};

    MPI_Datatype edge_types[2] = {MPI_INT, MPI_INT};
    MPI_Datatype mpi_edge_type;
    MPI_Aint     edge_offsets[2];
    edge_offsets[0] = offsetof(Edge, start);
    edge_offsets[1] = offsetof(Edge, end);
    MPI_Type_create_struct(nitems, blocklengths, edge_offsets, edge_types, &mpi_edge_type);
    MPI_Type_commit(&mpi_edge_type);

    MPI_Datatype adj_types[2] = {mpi_edge_type, mpi_edge_type};
    MPI_Datatype mpi_adj_list_type;
    MPI_Aint     adj_offsets[2];
    adj_offsets[0] = offsetof(AdjacencyListNode, forward);
    adj_offsets[1] = offsetof(AdjacencyListNode, reverse);
    MPI_Type_create_struct(nitems, blocklengths, adj_offsets, adj_types, &mpi_adj_list_type);
    MPI_Type_commit(&mpi_adj_list_type);

    return make_tuple(mpi_edge_type, mpi_adj_list_type);
}

Edge getNextEdge(Edge edge, vector<AdjacencyListNode> *list, int list_size){
    bool flag = false;
    vector<AdjacencyListNode> tmp = list[edge.start];
    for(size_t i = 0; i < tmp.size(); i++){
        Edge to_compare = tmp[i].forward; 
        if(flag){
            return to_compare;
        }
        if(to_compare.start == edge.start && to_compare.end == edge.end){
            flag = true;
        }
    }
    Edge not_exist = {NOT_EXIST, NOT_EXIST};
    return not_exist;
}

Edge getFirstNodeEdge(Edge edge, vector<AdjacencyListNode> *list, int list_size){
    vector<AdjacencyListNode> tmp = list[edge.start];
    return tmp[0].forward;
}

bool compareEdges(Edge edge0, Edge edge1){
    if(edge0.start != edge1.start)
        return false;
    if(edge0.end != edge1.end)
        return false;
    return true;
}

bool isForwardEdge(Edge edge){
    return edge.start < edge.end;
}

void printResult(int *level, string &tree){
    int size = tree.length();
    for(int i = 0; i < size-1; i++){
        cout << tree[i] << ":" << level[i] << ",";
    }
    cout << tree[size-1] << ":" << level[size-1] << endl;
}

void calculateAndSendEulersTourParts(vector<AdjacencyListNode> *ajd_list, int node_number, MPI_Status *status){
    MPI_Datatype mpi_edge_type, mpi_adj_list_type; 
    tie(mpi_edge_type, mpi_adj_list_type) = getCustomTypes();

    // process receive its node
    AdjacencyListNode proc_node;
    MPI_Recv(&proc_node, 1, mpi_adj_list_type, MASTER_ID, TAG, MPI_COMM_WORLD, status);

    // do Euler's tour in parallel
    Edge next = getNextEdge(proc_node.reverse, ajd_list, node_number);
    Edge edge_next;
    if(next.start == NOT_EXIST){
        edge_next = getFirstNodeEdge(proc_node.reverse, ajd_list, node_number);
    }
    else{
        edge_next = next;
    }
    MPI_Send(&edge_next, 1, mpi_edge_type, MASTER_ID, TAG, MPI_COMM_WORLD);
}

int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);
    int proc_id, processes_num;
    MPI_Comm_size(MPI_COMM_WORLD, &processes_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Status status;

    // prepare custom MPI data types
    MPI_Datatype mpi_edge_type, mpi_adj_list_type; 
    tie(mpi_edge_type, mpi_adj_list_type) = getCustomTypes();

    string input_tree = argv[1];
    size_t node_number = input_tree.length();

    // check if tree is not only root node
    if(node_number < 2){
        cout << input_tree[0] << ":0" << endl;
        MPI_Finalize();
        exit(0);
    }

    vector<AdjacencyListNode> ajd_list[node_number];
    
    // get adjacency list for specific node and send it to other processes
    if(proc_id < node_number){
        vector<AdjacencyListNode> tmp = getListForNode(proc_id, input_tree);
        int size = (int) tmp.size();
        for(int i = 0; i < processes_num; i++){
            MPI_Send(&size, 1, MPI_INT, i, TAG, MPI_COMM_WORLD);
        }
        for(int i = 0; i < size; i++){
            AdjacencyListNode el = tmp[i];
            for(int j = 0; j < processes_num; j++){
                MPI_Send(&el, 1, mpi_adj_list_type, j, TAG+1, MPI_COMM_WORLD);
            }
        }   
    }

    // processes receive parts of adjacency list and assemble them to complete list
    for(size_t i = 0; i < node_number; i++){
        int size;
        vector<AdjacencyListNode> list;
        MPI_Recv(&size, 1, MPI_INT, i, TAG, MPI_COMM_WORLD, &status);
        for(int j = 0; j < size; j++){
            AdjacencyListNode tmp;
            MPI_Recv(&tmp, 1, mpi_adj_list_type, i, TAG+1, MPI_COMM_WORLD, &status);
            list.push_back(tmp);
        }
        ajd_list[i] = list;
    }

    // array to store which edge was given to which process
    Edge given_edge[processes_num];

    // assign node to process
    if(proc_id == MASTER_ID){
        int iter = 0;
        for(int i = 0; i < node_number; i++){
            for(size_t j = 0; j < ajd_list[i].size(); j++){
                AdjacencyListNode tmp = ajd_list[i][j];
                MPI_Send(&tmp, 1, mpi_adj_list_type, iter, TAG, MPI_COMM_WORLD);
                given_edge[iter] = tmp.forward;
                iter++;
            }
        }
    }

    calculateAndSendEulersTourParts(&ajd_list[0], node_number, &status);
    
    // master process receive parts of Euler tour from processes
    Edge euler_tour[processes_num];
    if(proc_id == MASTER_ID){
        Edge received_edge[processes_num];
        for(int i = 0; i < processes_num; i++){
            Edge received;
            MPI_Recv(&received, 1, mpi_edge_type, i, TAG, MPI_COMM_WORLD, &status);
            received_edge[i] = received;
        }

        // assemble Euler's tour path
        Edge tmp = given_edge[0];
        euler_tour[0] = tmp;
        for(int i = 1; i < processes_num; i++){
            int iter = 0;
            while(iter < processes_num){
                if(compareEdges(tmp, given_edge[iter])){
                    euler_tour[i] = received_edge[iter];
                    tmp = received_edge[iter];
                    break;
                }
                iter++;
            }
        }
    }
    MPI_Bcast(euler_tour, processes_num, mpi_edge_type, MASTER_ID, MPI_COMM_WORLD);

    int weights[processes_num];
    int w;
    Edge proc_edge = euler_tour[proc_id];
    if(proc_edge.start < proc_edge.end){
        w = FORWARD;
    }
    else{
        w = REVERSE;
    }
    MPI_Send(&w, 1, MPI_INT, MASTER_ID, TAG, MPI_COMM_WORLD);

    if(proc_id == MASTER_ID){
        for(int i = 0; i < processes_num; i++){
            int tmp;
            MPI_Recv(&tmp, 1, MPI_INT, i, TAG, MPI_COMM_WORLD, &status);
            weights[i] = tmp;
        }
    }
    MPI_Bcast(weights, processes_num, MPI_INT, MASTER_ID, MPI_COMM_WORLD);
    
    // suffix sum
    int suffix_val[processes_num];
    int succ = proc_id + 1;
    int suffix;
    if(succ >= processes_num){
        suffix = weights[proc_id];
    }
    else {
        suffix = weights[proc_id];
    }

    MPI_Send(&suffix, 1, MPI_INT, MASTER_ID, TAG, MPI_COMM_WORLD);
    if(proc_id == MASTER_ID){
        for(int i = 0; i < processes_num; i++){
            int tmp;
            MPI_Recv(&tmp, 1, MPI_INT, i, TAG, MPI_COMM_WORLD, &status);
            suffix_val[i] = tmp;
        }
    }
    MPI_Bcast(suffix_val, processes_num, MPI_INT, MASTER_ID, MPI_COMM_WORLD);

    int limit = (int) (log2(processes_num)+0.5);
    int successor, successor_i;
    for(int k = 0; k <= limit; k++){
        successor_i = proc_id + (int) pow(BASE, k);
        if(successor_i >= processes_num)
            successor = 0;
        else
            successor = suffix_val[successor_i];
        suffix_val[proc_id] += successor;

        MPI_Send(&(suffix_val[proc_id]), 1, MPI_INT, MASTER_ID, TAG, MPI_COMM_WORLD);

        if(proc_id == MASTER_ID){
            for(int i = 0; i < processes_num; i++){
                int tmp;
                MPI_Recv(&tmp, 1, MPI_INT, i, TAG, MPI_COMM_WORLD, &status);
                suffix_val[i] = tmp;
            }
        }
        MPI_Bcast(suffix_val, processes_num, MPI_INT, MASTER_ID, MPI_COMM_WORLD);
    }
    
    int result[node_number];
    int node_level = 0;
    int pos = -1;
    if(isForwardEdge(euler_tour[proc_id])){
        pos = euler_tour[proc_id].end;
        node_level = suffix_val[proc_id] + 1;
    }
    MPI_Send(&pos, 1, MPI_INT, MASTER_ID, TAG, MPI_COMM_WORLD);
    MPI_Send(&node_level, 1, MPI_INT, MASTER_ID, TAG+1, MPI_COMM_WORLD);
    
    if(proc_id == MASTER_ID){
        for(int i = 0; i < processes_num; i++){
            int tmp_value, tmp_pos;
            MPI_Recv(&tmp_pos, 1, MPI_INT, i, TAG, MPI_COMM_WORLD, &status);
            MPI_Recv(&tmp_value, 1, MPI_INT, i, TAG+1, MPI_COMM_WORLD, &status);
            if(tmp_pos != NOT_EXIST)
                result[tmp_pos] = tmp_value;
        }
        result[MASTER_ID] = 0;

        printResult(result, input_tree);
    }

    MPI_Finalize();
}