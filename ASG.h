#include"bigraph.h"

typedef struct Node {//³¬½Úµã
	int alpha, beta;
	vector<int>left, right;
	vector<int>neighbor;
}Node;

class ASG {
public:
	vector<Node>asg;
	vector<vector<int>>uIndex, vIndex;
	vector<int>visited;
	int rec;
	vector<vector<vector<pair<int, int>>>>block;
	vector<vector<pair<int, int>>>rank;
public:
	ASG();
	ASG(BiGraph& g);
	void addNode(BiGraph& g);
	void creatSEviaNode(BiGraph& g);
	void creatSEviaNode(BiGraph& g, set<int>& nodeQ);
	void addEdge();
	void BFSonASG(int a, int b, int w, vector<int>& wmap);
	void Query(int alpha, int beta, int u, bool isU, vector<bool>&leftResult, vector<bool>& rightResult);
};
int verifyCom(BiGraph& g, vector<bool>& leftResult, vector<bool>& rightResult, int alpha, int beta, int q, bool isLeft);

