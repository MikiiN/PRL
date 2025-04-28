/*
* Project: PRL-proj2
* Author:  Michal Žatečka xzatec02
* Date:    27.04.2025
*/

#include <bits/stdc++.h>
#include <mpi.h>
#include <iostream>
#include <string>
#include <tuple>

using namespace std;

#define BASE 2
#define NOT_EXIST -1
#define ROOT_INDEX 0

typedef struct g_edge {
    size_t start;
    size_t end;
}Edge;

typedef struct adj_node {
    Edge forward;
    Edge reverse;
} AdjacencyListNode;

int getNodeDepth(size_t node_index, size_t tree_len){
    int i = 0;
    size_t sum = 0;
    while(sum <= tree_len){
        sum += (int) pow(BASE, i);
        if (node_index < sum){
            break;
        }
        i++;
    }
    return i;
}

tuple<size_t,size_t> getNodeChildrenIndexes(size_t node_index, size_t tree_len){
    int depth = getNodeDepth(node_index, tree_len);
    size_t shift = (size_t) pow(BASE, depth);
    int firstChild = node_index+shift;
    int secondChild = firstChild+1;

    if(firstChild >= tree_len){
        firstChild = NOT_EXIST;
        secondChild = NOT_EXIST;
    }
    else if(secondChild >= tree_len)
        secondChild = NOT_EXIST;
    return make_tuple(firstChild, secondChild);
}

size_t getNodeParent(size_t node_index, size_t tree_len){
    int depth = getNodeDepth(node_index, tree_len);
    if(depth <= 0)
        return NOT_EXIST;
    int is_left_child = node_index % BASE;
    size_t shift = (size_t) pow(BASE, depth-1);
    if(!is_left_child)
        shift++;
    return node_index-shift;
    
}

AdjacencyListNode createNode(size_t node0, size_t node1){
    Edge edge = {node0, node1};
    Edge revert = {node1, node0};
    AdjacencyListNode node = {edge, revert};
    return node;
}

vector<AdjacencyListNode> getListForNode(size_t node_index, string tree){
    vector<AdjacencyListNode> result;
    size_t left_child, right_child;
    size_t tree_len = tree.length();
    size_t parent = getNodeParent(node_index, tree_len);
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

int main(int argc, char *argv[]){
    string input_tree = argv[1];
    cout << input_tree << endl;
}