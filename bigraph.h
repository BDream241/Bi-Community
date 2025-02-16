#include"base.h"

class BiGraph {
public:
	//变量
	int n1, n2, m, maxAlpha, maxBeta, maxK, maxUDeg, maxVDeg;
	vector<int>uDeg, vDeg;
	vector<vector<int>> uNeighbor, vNeighbor;
	vector<vector<tuple<int, int, int>>> uNumber, vNumber, uBeta, vBeta;
	
public:
	//base
	BiGraph();
	BiGraph(string str);
	void addEdge(int u, int v);
	void deleteEdge(int u, int v);
	void loadGraph(string str);
	void init(int n1, int n2);
	void print();
	void printBiCore();

public:
	//alpha、beta分解
	int decompose();
	int coreDecompose();
	int alphaDecompose(int alpha, vector<int>& leftDeg, vector<int>& rightDeg, vector<bool>& uDelete,
		vector<bool>& vDelete, vector<vector<int>>& uDAG,
		vector<vector<int>>& vDAG, vector<vector<pair<int, bool>>>& order);
	int betaDecompose(int beta, vector<int>& leftDeg, vector<int>& rightDeg, vector<bool>& uDelete,
		vector<bool>& vDelete, vector<vector<int>>& uDAG,
		vector<vector<int>>& vDAG, vector<vector<pair<int, bool>>>& order);
	int dealDAG(pair<int, bool>ab, vector<vector<int>>& uDAG, vector<vector<int>>& vDAG, vector<vector<pair<int, bool>>>& order);
	int updateNumber(int alpha, int beta, int x, pair<int, bool>p, bool isAlpha);

};
#pragma once

