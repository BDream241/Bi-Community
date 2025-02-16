#include<stdio.h>
#include<set> 
#include<map> 
#include<queue>
#include<tuple>
#include<vector>
#include<chrono>
#include<iostream> 
#include<algorithm>
#include<filesystem>
#include<unordered_set>
#include<unordered_map>
using namespace std;

class UnionFind {
public:
    int num;
    vector<int> father, rank;
public:
    UnionFind() {
       num = 0;
    }
    UnionFind(int n) {
        num = n;
        father.resize(n);
        for (int i = 0; i < n; i++)
            father[i] = i;
        rank.resize(n, 0);
    }
    int add() {
        father.push_back(num);
        rank.push_back(0);
        num++;
        return num - 1;
    }
    int Find(int x) {
        int f = x;
        if (x == father[x])
            return x;
        while (x != father[x])
            x = father[x];
        return father[f] = x;
    }
    int merge(int x, int y) {
        if (x != y) {
            if (rank[x] > rank[y]) {
                father[y] = x;
                rank[y]++;
            }
            else {
                father[x] = y;
                rank[x]++;
            } 
            return 1;
        }
        return 0;
    }
    //void UnionFind
    int BatchUnite(set<int>& root_nodes) {
        if (root_nodes.size() == 0) {
            father.push_back(num);
            rank.push_back(0);
            num++;
            return num - 1;
        }
        int root=0, flag = 0;
        for (int node : root_nodes) {
            if (flag == 0) {
                flag = 1;
                root = node;
            }
            if (rank[root] < rank[node])
                root = node;
        }
        for (int node : root_nodes) {
            father[node] = root;
            if (node != root)
                rank[node]++;
        }
        return root;
    }
};

#pragma once
