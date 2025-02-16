#include"MHG.h"

class Maintenance {
public:
	int n1, n2;
	vector<vector<int>>uNeighbor, vNeighbor;
	vector<vector<int>>uNumber, vNumber;
	vector<bool>uT, vT, uC, vC;
	vector<int>uDeg, vDeg;
	BiGraph g;
	ASG asg;
	set<int>nodeQ;
public:
	Maintenance();
	Maintenance(BiGraph& g, ASG& a);
	void init(BiGraph& g, ASG& a);
	Maintenance(string adds, ASG& myAsg, BiGraph& g);
	bool update(int u, bool isLeft, int node);
	int findOff(int u, bool isLeft, int value, bool isAlpha);
	void dealC(int a, int b, vector<pair<bool, int>>R);
	void dealC2(int a, int b, vector<pair<bool, int>>R);
public:
	void edgeInsert(int u, int v);
	void alphaIncrease(int u, int v, int alpha, int beta, vector<pair<int, int>>& H);
	void betaIncrease(int u, int v, int alpha, int beta, vector<pair<int, int>>& H);
	void alphaIncrease2(int u, int v, int alpha, int beta, vector<pair<int, int>>& H);
	void betaIncrease2(int u, int v, int alpha, int beta, vector<pair<int, int>>& H);
	void neighborBaseIns(int u, bool isLeft, int a, int b);
	void neighborBaseIns2(int u, bool isLeft, int a, int b);
	void deleteCand(int w, bool isLeft, int a, int b);
public:
	void edgeDelete(int u, int v);
	void alphaDecrease(int u, int v, int alpha, int beta, vector<pair<bool, int>>& R);
	void betaDecrease(int u, int v, int alpha, int beta, vector<pair<bool, int>>& R);
	void alphaDecrease2(int u, int v, int alpha, int beta, vector<pair<bool, int>>& R);
	void betaDecrease2(int u, int v, int alpha, int beta, vector<pair<bool, int>>& R);
	vector<pair<bool, int>> neighborBaseDel(int u, bool isLeft, int a, int b);
	void addCand(vector<pair<bool, int>>&R, int w, bool isLeft, int a, int b, queue<pair<bool, int>>& temp, vector<bool>& uS, vector<bool>& vS);

};
#pragma once
