#include"ASG.h"

typedef struct treeNode {//tree
	int x1, y1, x2, y2;
	int father;
	vector<int>children;
	vector<int>leaf;
}treeNode;
typedef struct mNode {//mhg
	int x1, y1, x2, y2;
	int begin, end;
	vector<int>slide;
	vector<int>neighbor;
}mNode;

class MHG {
public:
	vector<mNode>tf;
	vector<Node>gf;
	vector<vector<int>>stairs;
	vector<vector<int>>uIndex, vIndex;
	vector<treeNode>qtree;
	int num;
public:
	MHG();
	MHG(vector<Node>& asg, int k, int a, int b, int n1, int n2);
	void gridding(vector<Node>& asg, int c);
	void pre(int c, vector<int>& order, vector<int>&father);
	void connection(vector<int>& father);
	void Query(int alpha, int beta, int u, bool isU, vector<bool>& leftResult, vector<bool>& rightResult);
	void QueryViaMHG(int alpha, int beta, int u, bool isU, vector<bool>& leftResult, vector<bool>& rightResult);
};

#pragma once
