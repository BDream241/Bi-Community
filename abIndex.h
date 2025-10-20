#pragma once
#include"bigraph.h"


typedef pair<int, vector<int>> group;

typedef struct node {
	int alpha, beta;
	vector<int>neighbor;
	vector<int>left, right;
	vector<group>leftGroup, rightGroup;
}node;

typedef struct tree {
	pair<int,int> lowLeft, topRight;
	vector<int>child;
	vector<int>leaf;
}tree;

class abIndex {
public:
	vector<tree> T;
	vector<node> S, SOpt;
	vector<vector<int>> uNumber, vNumber;
	vector<int>visited;
	int rec;
	vector<vector<set<pair<int, int>>>>block,block_2;
	vector<vector<pair<int, int>>>rank;
	vector<int>wmap;
	double ASG, ASGL, ASGT;
	vector<bool>uT, vT, uC, vC;
	vector<int>uDeg, vDeg;
	set<int>nodeQ;
	BiGraph *g;
	vector<vector<int>>uNeighbor, vNeighbor;
	int n1, n2;
public:
	abIndex();
	abIndex(BiGraph& graph);
	void creatSG(BiGraph& graph);
	void creatSG1(BiGraph& graph);
	void creatSG2(BiGraph& graph);
	void creatSGOpt(BiGraph& graph);
	void creatSGOpt_2(BiGraph& graph);
	void creatTree();
	void addEdge();
	void addEdge_2();
	void BFSonASG(int a, int b, int w, vector<int>& wmap);
	void query(int q, int alpha, int beta, vector<bool>& leftResult, vector<bool>& rightResult);
	void queryOpt(int q, int alpha, int beta, vector<bool>& leftResult, vector<bool>& rightResult, int k);
	void queryCore(int alpha, int beta, vector<int>& leftCom, vector<int>& rightCom, int k);
	void BFSonASG(int s, int alpha, int beta, vector<int>& leftCom, vector<int>& rightCom, vector<bool>& F, int k, int sum);
	void writeIndex();

	void Maintenance(BiGraph& g);
	void initMain(BiGraph& g);
	void Maintenance(string adds, BiGraph& g);
	void Maintenance(string adds, BiGraph graph, int mup);
	bool update(int u, bool isLeft, int node);
	int findOff(int u, bool isLeft, int value, bool isAlpha);
	void dealC(int a, int b, vector<pair<bool, int>>R);
	void dealC2(int a, int b, vector<pair<bool, int>>R);
	void edgeInsert(int u, int v);
	void alphaIncrease(int u, int v, int alpha, int beta, vector<pair<int, int>>& H);
	void betaIncrease(int u, int v, int alpha, int beta, vector<pair<int, int>>& H);
	void alphaIncrease2(int u, int v, int alpha, int beta, vector<pair<int, int>>& H);
	void alphaIncrease2(int u, int v, int alpha, int beta, vector<pair<int, int>>& H);
	void betaIncrease2(int u, int v, int alpha, int beta, vector<pair<int, int>>& H);
	void neighborBaseIns(int u, bool isLeft, int a, int b);
	void neighborBaseIns2(int u, bool isLeft, int a, int b);
	void deleteCand(int w, bool isLeft, int a, int b);
	void edgeDelete(int u, int v);
	void alphaDecrease(int u, int v, int alpha, int beta, vector<pair<bool, int>>& R);
	void betaDecrease(int u, int v, int alpha, int beta, vector<pair<bool, int>>& R);
	void alphaDecrease2(int u, int v, int alpha, int beta, vector<pair<bool, int>>& R);
	void betaDecrease2(int u, int v, int alpha, int beta, vector<pair<bool, int>>& R);
	vector<pair<bool, int>> neighborBaseDel(int u, bool isLeft, int a, int b);
	void addCand(vector<pair<bool, int>>& R, int w, bool isLeft, int a, int b, queue<pair<bool, int>>& temp, vector<bool>& uS, vector<bool>& vS);
	void creatSEviaVertex(set<int>& nodeQ);
	void creatSEviaNode_2(set<int>& nodeQ);
};

abIndex::abIndex(){}

abIndex::abIndex(BiGraph& graph) {
	uNumber = graph.uNumber;
	graph.uNumber.clear();
	vNumber = graph.vNumber;
	graph.vNumber.clear();
	//creatSG(graph);
	creatSGOpt_2(graph);
	//cout << S.size() << "  " << SOpt.size() << endl;
	creatTree();
	//writeIndex();
}

void abIndex::creatSG(BiGraph& graph) {
	for (auto i : graph.cluster) {
		node s;
		s.alpha = get<0>(i);
		s.beta = get<1>(i);
		S.push_back(s);
	}

	for (int i = 0; i < uNumber.size(); i++)
		for (auto j : uNumber[i])
			S[j].left.push_back(i);
	for (int i = 0; i < vNumber.size(); i++)
		for (auto j : vNumber[i])
			S[j].right.push_back(i);

	for (auto i : graph.block1) {
		if (i.first == i.second)
			continue;
		S[i.first].neighbor.push_back(i.second);
		S[i.second].neighbor.push_back(i.first);
	}
	//graph.block1.clear();
	for (auto i : graph.block2) {
		if (i.first == i.second)
			continue;
		S[i.first].neighbor.push_back(i.second);
		S[i.second].neighbor.push_back(i.first);
	}
	//graph.block2.clear();
}

void abIndex::creatSG1(BiGraph& graph) {
	for (auto i : graph.cluster) {
		node s;
		s.alpha = get<0>(i);
		s.beta = get<1>(i);
		SOpt.push_back(s);
	}
	for (int i = 0; i < uNumber.size(); i++)
		for (auto j : uNumber[i])
			SOpt[j].left.push_back(i);
	for (int i = 0; i < vNumber.size(); i++)
		for (auto j : vNumber[i])
			SOpt[j].right.push_back(i);
	for (auto i : graph.block1) {
		if (i.first == i.second)
			continue;
		SOpt[i.first].neighbor.push_back(i.second);
		SOpt[i.second].neighbor.push_back(i.first);
	}
	visited.resize(SOpt.size(), -1);
	rec = 0;
	block.resize(graph.maxAlpha + 1);
	for (auto i : graph.block2) {
		if (i.first == i.second)
			continue;
		int a = min(SOpt[i.first].alpha, SOpt[i.second].alpha);
		int b = min(SOpt[i.first].beta, SOpt[i.second].beta);
		if (block[a].size() <= b)
			block[a].resize(b + 1);
		block[a][b].insert(make_pair(i.first, i.second));
	}
	addEdge();
	S = SOpt;
	SOpt.clear();
}

void abIndex::creatSG2(BiGraph& graph) {
	for (auto i : graph.cluster) {
		node s;
		s.alpha = get<0>(i);
		s.beta = get<1>(i);
		SOpt.push_back(s);
	}
	for (auto i : graph.block1) {
		if (i.first == i.second)
			continue;
		SOpt[i.first].neighbor.push_back(i.second);
		SOpt[i.second].neighbor.push_back(i.first);
	}
	for (auto i : graph.block2) {
		if (i.first == i.second)
			continue;
		SOpt[i.first].neighbor.push_back(i.second);
		SOpt[i.second].neighbor.push_back(i.first);
	}
	vector<vector<vector<int>>>left, right;
	left.resize(SOpt.size());
	right.resize(SOpt.size());
	for (int i = 0; i < uNumber.size(); i++)
		for (int j = 0; j < uNumber[i].size(); j++) {
			int c = uNumber[i][j];
			if (SOpt[c].alpha <= graph.maxK) {
				if (j == 0) {
					if (left[c].size() == 0)
						left[c].resize(1);
					left[c][j].push_back(i);
				}
				else {
					int c0 = uNumber[i][j - 1];
					if (left[c].size() <= SOpt[c0].alpha)
						left[c].resize(SOpt[c0].alpha + 1);
					left[c][SOpt[c0].alpha].push_back(i);
				}
			}
			else {
				if (j == uNumber[i].size() - 1) {
					if (left[c].size() == 0)
						left[c].resize(1);
					left[c][0].push_back(i);
				}
				else {
					int c2 = uNumber[i][j + 1];
					if (left[c].size() <= SOpt[c2].beta)
						left[c].resize(SOpt[c2].beta + 1);
					left[c][SOpt[c2].beta].push_back(i);
				}
			}
		}
	for (int i = 0; i < vNumber.size(); i++)
		for (int j = 0; j < vNumber[i].size(); j++) {
			int c = vNumber[i][j];
			if (SOpt[c].alpha <= graph.maxK) {
				if (j == 0) {
					if (right[c].size() == 0)
						right[c].resize(1);
					right[c][j].push_back(i);
				}
				else {
					int c0 = vNumber[i][j - 1];
					if (right[c].size() <= SOpt[c0].alpha)
						right[c].resize(SOpt[c0].alpha + 1);
					right[c][SOpt[c0].alpha].push_back(i);
				}
			}
			else {
				if (j == vNumber[i].size() - 1) {
					if (right[c].size() == 0)
						right[c].resize(1);
					right[c][0].push_back(i);
				}
				else {
					int c2 = vNumber[i][j + 1];
					if (right[c].size() <= SOpt[c2].beta)
						right[c].resize(SOpt[c2].beta + 1);
					right[c][SOpt[c2].beta].push_back(i);
				}
			}
		}
	for (int i = 0; i < SOpt.size(); i++) {
		for (int j = 0; j < left[i].size(); j++)
			if (left[i][j].size() > 0)
				SOpt[i].leftGroup.push_back({ j,left[i][j] });
		for (int j = 0; j < right[i].size(); j++)
			if (right[i][j].size() > 0)
				SOpt[i].rightGroup.push_back({ j,right[i][j] });
	}
	left.clear();
	right.clear();
}

void abIndex::creatSGOpt(BiGraph& graph) {
	for (auto i : graph.cluster) {
		node s;
		s.alpha = get<0>(i);
		s.beta = get<1>(i);
		SOpt.push_back(s);
	}
	for (auto i : graph.block1) {
		if (i.first == i.second)
			continue;
		SOpt[i.first].neighbor.push_back(i.second);
		SOpt[i.second].neighbor.push_back(i.first);
	}
	graph.block1.clear();

	visited.resize(SOpt.size(), -1);
	rec = 0;
	block.resize(graph.maxAlpha + 1);
	for (auto i : graph.block2) {
		if (i.first == i.second)
			continue;
		int a = min(SOpt[i.first].alpha, SOpt[i.second].alpha);
		int b = min(SOpt[i.first].beta, SOpt[i.second].beta);
		if (block[a].size() <= b)
			block[a].resize(b + 1);
		block[a][b].insert(make_pair(i.first, i.second));
	}
	graph.block2.clear();
	addEdge();

	vector<vector<vector<int>>>left, right;
	left.resize(SOpt.size());
	right.resize(SOpt.size());
	for (int i = 0; i < uNumber.size(); i++)
		for (int j = 0; j < uNumber[i].size(); j++) {
			int c = uNumber[i][j];
			if (SOpt[c].alpha <= graph.maxK) {
				if (j == 0) {
					if (left[c].size() == 0)
						left[c].resize(1);
					left[c][j].push_back(i);
				}
				else {
					int c0 = uNumber[i][j - 1];
					if (left[c].size() <= SOpt[c0].alpha)
						left[c].resize(SOpt[c0].alpha + 1);
					left[c][SOpt[c0].alpha].push_back(i);
				}
			}
			else {
				if (j == uNumber[i].size() - 1) {
					if (left[c].size() == 0)
						left[c].resize(1);
					left[c][0].push_back(i);
				}
				else {
					int c2 = uNumber[i][j + 1];
					if (left[c].size() <= SOpt[c2].beta)
						left[c].resize(SOpt[c2].beta + 1);
					left[c][SOpt[c2].beta].push_back(i);
				}
			}
		}
	for (int i = 0; i < vNumber.size(); i++)
		for (int j = 0; j < vNumber[i].size(); j++) {
			int c = vNumber[i][j];
			if (SOpt[c].alpha <= graph.maxK) {
				if (j == 0) {
					if (right[c].size() == 0)
						right[c].resize(1);
					right[c][j].push_back(i);
				}
				else {
					int c0 = vNumber[i][j - 1];
					if (right[c].size() <= SOpt[c0].alpha)
						right[c].resize(SOpt[c0].alpha + 1);
					right[c][SOpt[c0].alpha].push_back(i);
				}
			}
			else {
				if (j == vNumber[i].size() - 1) {
					if (right[c].size() == 0)
						right[c].resize(1);
					right[c][0].push_back(i);
				}
				else {
					int c2 = vNumber[i][j + 1];
					if (right[c].size() <= SOpt[c2].beta)
						right[c].resize(SOpt[c2].beta + 1);
					right[c][SOpt[c2].beta].push_back(i);
				}
			}
		}
	for (int i = 0; i < SOpt.size(); i++) {
		for (int j = 0; j < left[i].size(); j++) 
			if (left[i][j].size() > 0) 
				SOpt[i].leftGroup.push_back({ j,left[i][j] });
		for (int j = 0; j < right[i].size(); j++) 
			if (right[i][j].size() > 0) 
				SOpt[i].rightGroup.push_back({ j,right[i][j] });
	}
	left.clear();
	right.clear();
}

void abIndex::creatSGOpt_2(BiGraph& graph) {
	for (auto i : graph.cluster) {
		node s;
		s.alpha = get<0>(i);
		s.beta = get<1>(i);
		SOpt.push_back(s);
	}

	visited.resize(SOpt.size(), -1);
	rec = 0;
	block.resize(graph.maxAlpha + 1);
	block_2.resize(graph.maxAlpha + 1);
	for (auto i : graph.block1) {
		if (i.first == i.second)
			continue;
		int a = min(SOpt[i.first].alpha, SOpt[i.second].alpha);
		int b = min(SOpt[i.first].beta, SOpt[i.second].beta);
		if (block[a].size() <= b)
			block[a].resize(b + 1);
		block[a][b].insert(make_pair(i.first, i.second));
	}
	graph.block1.clear();
	for (auto i : graph.block2) {
		if (i.first == i.second)
			continue;
		int a = min(SOpt[i.first].alpha, SOpt[i.second].alpha);
		int b = min(SOpt[i.first].beta, SOpt[i.second].beta);
		if (block_2[a].size() <= b)
			block_2[a].resize(b + 1);
		block_2[a][b].insert(make_pair(i.first, i.second));
	}
	graph.block2.clear();
	addEdge_2();

	vector<vector<vector<int>>>left, right;
	left.resize(SOpt.size());
	right.resize(SOpt.size());
	for (int i = 0; i < uNumber.size(); i++)
		for (int j = 0; j < uNumber[i].size(); j++) {
			int c = uNumber[i][j];
			if (SOpt[c].alpha <= graph.maxK) {
				if (j == 0) {
					if (left[c].size() == 0)
						left[c].resize(1);
					left[c][j].push_back(i);
				}
				else {
					int c0 = uNumber[i][j - 1];
					if (left[c].size() <= SOpt[c0].alpha)
						left[c].resize(SOpt[c0].alpha + 1);
					left[c][SOpt[c0].alpha].push_back(i);
				}
			}
			else {
				if (j == uNumber[i].size() - 1) {
					if (left[c].size() == 0)
						left[c].resize(1);
					left[c][0].push_back(i);
				}
				else {
					int c2 = uNumber[i][j + 1];
					if (left[c].size() <= SOpt[c2].beta)
						left[c].resize(SOpt[c2].beta + 1);
					left[c][SOpt[c2].beta].push_back(i);
				}
			}
		}
	for (int i = 0; i < vNumber.size(); i++)
		for (int j = 0; j < vNumber[i].size(); j++) {
			int c = vNumber[i][j];
			if (SOpt[c].alpha <= graph.maxK) {
				if (j == 0) {
					if (right[c].size() == 0)
						right[c].resize(1);
					right[c][j].push_back(i);
				}
				else {
					int c0 = vNumber[i][j - 1];
					if (right[c].size() <= SOpt[c0].alpha)
						right[c].resize(SOpt[c0].alpha + 1);
					right[c][SOpt[c0].alpha].push_back(i);
				}
			}
			else {
				if (j == vNumber[i].size() - 1) {
					if (right[c].size() == 0)
						right[c].resize(1);
					right[c][0].push_back(i);
				}
				else {
					int c2 = vNumber[i][j + 1];
					if (right[c].size() <= SOpt[c2].beta)
						right[c].resize(SOpt[c2].beta + 1);
					right[c][SOpt[c2].beta].push_back(i);
				}
			}
		}
	for (int i = 0; i < SOpt.size(); i++) {
		for (int j = 0; j < left[i].size(); j++)
			if (left[i][j].size() > 0)
				SOpt[i].leftGroup.push_back({ j,left[i][j] });
		for (int j = 0; j < right[i].size(); j++)
			if (right[i][j].size() > 0)
				SOpt[i].rightGroup.push_back({ j,right[i][j] });
	}
	left.clear();
	right.clear();
}

void abIndex::creatTree() {
	tree t;
	T.push_back(t);
	T[0].lowLeft = { 1,1 };
	int maxA = 0, maxB = 0;
	for (int i = 0; i < SOpt.size(); i++) {
		bool add = true;
		for (auto j : SOpt[i].neighbor)
			if (SOpt[i].alpha <= SOpt[j].alpha && SOpt[i].beta <= SOpt[j].beta) {
				add = false;
				break;
			}
		if (i == 207763 && add)
			cout << T[0].leaf.size() << endl;
		if (add) {
			T[0].leaf.push_back(i);
			if (SOpt[i].alpha > maxA)
				maxA = SOpt[i].alpha;
			if (SOpt[i].beta > maxB)
				maxB = SOpt[i].beta;
		}
	}
	T[0].topRight = { maxA,maxB };
	for (int i = 0; i < T.size(); i++) {
		if ((T[i].topRight.first - T[i].lowLeft.first <= 2) || (T[i].topRight.second - T[i].lowLeft.second <= 2))
			continue;
		int midA = (T[i].topRight.first + T[i].lowLeft.first) / 2;
		int midB = (T[i].topRight.second + T[i].lowLeft.second) / 2;
		tree t1, t2, t3, t4;
		for (auto i : T[i].leaf) {
			if (SOpt[i].alpha >= midA && SOpt[i].beta >= midB)
				t1.leaf.push_back(i);
			else if (SOpt[i].alpha >= midA && SOpt[i].beta <= midB)
				t2.leaf.push_back(i);
			else if (SOpt[i].alpha <= midA && SOpt[i].beta <= midB)
				t3.leaf.push_back(i);
			else
				t4.leaf.push_back(i);
		}
		T[i].leaf.clear();
		if (t1.leaf.size() > 0) {
			t1.lowLeft = { midA,midB };
			t1.topRight = T[i].topRight;
			T[i].child.push_back(T.size());
			T.push_back(t1);
		}
		if (t2.leaf.size() > 0) {
			t2.lowLeft = { midA,T[i].lowLeft.second };
			t2.topRight = { T[i].topRight.first, midB };
			T[i].child.push_back(T.size());
			T.push_back(t2);
		}
		if (t3.leaf.size() > 0) {
			t3.lowLeft = T[i].lowLeft;
			t3.topRight = { midA, midB };
			T[i].child.push_back(T.size());
			T.push_back(t3);
		}
		if (t4.leaf.size() > 0) {
			t4.lowLeft = { T[i].lowLeft.first,midB };
			t4.topRight = { midA, T[i].topRight.second };
			T[i].child.push_back(T.size());
			T.push_back(t4);
		}
	}
}

void abIndex::addEdge() {
	int max = 0;
	for (int i = 1; i < block.size(); i++)
		if (max < (i + block[i].size() - 1))
			max = i + block[i].size() - 1;
	rank.clear();
	rank.resize(max + 1);
	for (int i = 1; i < block.size(); i++)
		for (int j = 1; j < block[i].size(); j++)
			rank[i + j].push_back(make_pair(i, j));

	for (int i = max; i > 1; i--)
		for (int j = 0; j < rank[i].size(); j++) {
			int a = rank[i][j].first;
			int b = rank[i][j].second;
			UnionFind uf;
			wmap.resize(SOpt.size(), 0);
			for (auto k : block[a][b]) {
				if (wmap[k.first] == 0) {
					wmap[k.first] = uf.add();
					BFSonASG(a, b, k.first, wmap);
				}
				if (wmap[k.second] == 0) {
					wmap[k.second] = uf.add();
					BFSonASG(a, b, k.second, wmap);
				}

				if (uf.merge(uf.Find(wmap[k.first]), uf.Find(wmap[k.second]))) {
					SOpt[k.first].neighbor.push_back(k.second);
					SOpt[k.second].neighbor.push_back(k.first);
				}
			}
			wmap.clear();
			block[a][b].clear();
		}
	block.clear();
}

void abIndex::addEdge_2() {
	int max = 0;
	for (int i = 1; i < block.size(); i++)
		if (max < (i + block[i].size() - 1))
			max = i + block[i].size() - 1;
	for (int i = 1; i < block_2.size(); i++)
		if (max < (i + block_2[i].size() - 1))
			max = i + block_2[i].size() - 1;
	rank.clear();
	rank.resize(max + 1);
	for (int i = 1; i < block.size(); i++)
		for (int j = 1; j < block[i].size(); j++)
			rank[i + j].push_back(make_pair(i, j));
	for (int i = 1; i < rank.size(); i++)
		rank[i].push_back(make_pair(0, 0));

	for (int i = 1; i < block_2.size(); i++)
		for (int j = 1; j < block_2[i].size(); j++)
			rank[i + j].push_back(make_pair(i, j));
	UnionFind uf;
	wmap.resize(SOpt.size(), 0);

	for (int i = max; i > 1; i--) {
		int falg = 0;
		for (int j = 0; j < rank[i].size(); j++) {
			int a = rank[i][j].first;
			int b = rank[i][j].second;

			if (a == 0 && b == 0) 
				falg = 1;

			if (falg == 0 && block[a].size()> b) {
				for (auto k : block[a][b]) {
					if (wmap[k.first] == 0) {
						wmap[k.first] = uf.add();
						BFSonASG(a, b, k.first, wmap);
					}
					if (wmap[k.second] == 0) {
						wmap[k.second] = uf.add();
						BFSonASG(a, b, k.second, wmap);
					}
					if (uf.merge(uf.Find(wmap[k.first]), uf.Find(wmap[k.second]))) {
					}
					SOpt[k.first].neighbor.push_back(k.second);
					SOpt[k.second].neighbor.push_back(k.first);
				}
				wmap.clear();
				block[a][b].clear();
			}
			else if(falg==1 && block_2[a].size() > b){
				for (auto k : block_2[a][b]) {
					if (wmap[k.first] == 0) {
						wmap[k.first] = uf.add();
						BFSonASG(a, b, k.first, wmap);
					}
					if (wmap[k.second] == 0) {
						wmap[k.second] = uf.add();
						BFSonASG(a, b, k.second, wmap);
					}
					if (uf.merge(uf.Find(wmap[k.first]), uf.Find(wmap[k.second]))) {
						SOpt[k.first].neighbor.push_back(k.second);
						SOpt[k.second].neighbor.push_back(k.first);
					}
				}
				wmap.clear();
				block_2[a][b].clear();
			}
		}
	}
	block.clear();
	block_2.clear();
}

void abIndex::BFSonASG(int a, int b, int w, vector<int>& wmap) {
	rec++;
	queue<int>bfs;
	bfs.push(w);
	visited[w] = rec;
	while (!bfs.empty()) {
		int p = bfs.front();
		bfs.pop();
		for (auto nbr : SOpt[p].neighbor)
			if (visited[nbr] != rec) {
				visited[nbr] = rec;
				if (SOpt[nbr].alpha >= a && SOpt[nbr].beta >= b) {
					bfs.push(nbr);
					wmap[nbr] = wmap[w];
				}
			}
	}
}

void abIndex::query(int q, int alpha, int beta, vector<bool>& leftResult, vector<bool>& rightResult) {
	vector<bool> visitedCluster(S.size(), false);
	queue<int>Q;
	for (auto i : uNumber[q])
		if (S[i].alpha >= alpha && S[i].beta >= beta) {
			Q.push(i);
			visitedCluster[i] = true;
			break;
		}
	//int sum = 0;
	while (!Q.empty()) {
		int w = Q.front();
		//sum++;
		Q.pop();
		for (auto k : S[w].neighbor)
			if (!visitedCluster[k]) {
				visitedCluster[k] = true;
				if (S[k].alpha >= alpha && S[k].beta >= beta)
					Q.push(k);
			}
		for (auto k : S[w].left)
			leftResult[k] = true;
		for (auto k : S[w].right)
			rightResult[k] = true;
	}
	//cout << "baseSum: " << sum << endl;
}

void abIndex::queryOpt(int q, int alpha, int beta, vector<bool>& leftResult, vector<bool>& rightResult, int k) {
	vector<bool> visitedCluster(SOpt.size(), false);
	queue<int>Q;
	for (auto i : uNumber[q])
		if (SOpt[i].alpha >= alpha && SOpt[i].beta >= beta) {
			Q.push(i);
			visitedCluster[i] = true;
			break;
		}
	//int sum = 0;

	while (!Q.empty()) {
		int w = Q.front();
		//sum++;
		Q.pop();
		for (auto k : SOpt[w].neighbor)
			if (!visitedCluster[k]) {
				visitedCluster[k] = true;
				if (SOpt[k].alpha >= alpha && SOpt[k].beta >= beta)
					Q.push(k);
			}
		if (SOpt[w].alpha > k ) {
			for (int i = 0; i < SOpt[w].leftGroup.size(); i++) {
				if (SOpt[w].leftGroup[i].first >= beta)
					break;
				for (auto j : SOpt[w].leftGroup[i].second)
					leftResult[j] = true;
			}

			for (int i = 0; i < SOpt[w].rightGroup.size(); i++) {
				if (SOpt[w].rightGroup[i].first >= beta)
					break;
				for (auto j : SOpt[w].rightGroup[i].second)
					rightResult[j] = true;
			}
		}
		else {
			for (int i = 0; i < SOpt[w].leftGroup.size(); i++) {
				if (SOpt[w].leftGroup[i].first >= alpha)
					break;
				for (auto j : SOpt[w].leftGroup[i].second)
					leftResult[j] = true;
			}
			for (int i = 0; i < SOpt[w].rightGroup.size(); i++) {
				if (SOpt[w].rightGroup[i].first >= alpha)
					break;
				for (auto j : SOpt[w].rightGroup[i].second)
					rightResult[j] = true;
			}
		}
	}
	//cout << "optSum: " << sum << endl;
}

void abIndex::queryCore(int alpha, int beta, vector<int>& leftCom, vector<int>& rightCom, int k) {
	vector<int>Q;
	int sum = 0;
	vector<bool>F(SOpt.size(), false);
	if (T[0].topRight.first < alpha || T[0].topRight.second < beta)
		return;
	Q.push_back(0);
	for (int i = 0; i < Q.size(); i++) {
		int t = Q[i];
		if (T[t].child.size() == 0) {
			for (auto j : T[t].leaf)
				if (SOpt[j].alpha >= alpha && SOpt[j].beta >= beta && !F[j]) {
					sum++;
					BFSonASG(j, alpha, beta, leftCom, rightCom, F, k, sum);
					//cout << SOpt[j].alpha << " " << SOpt[j].beta << endl;
				}
		}
		else  {
			for (auto j : T[t].child)
				if (T[j].topRight.first >= alpha && T[j].topRight.second >= beta)
					Q.push_back(j);
		}
	}
	//cout << sum << endl;
}

void abIndex::BFSonASG(int s, int alpha, int beta, vector<int>& leftCom, vector<int>& rightCom, vector<bool>&F, int k, int sum) {
	vector<bool> visitedCluster(SOpt.size(), false);
	queue<int>Q;
	Q.push(s);
	F[s] = true;
	while (!Q.empty()) {
		int w = Q.front();
		//sum++;
		Q.pop();
		for (auto k : SOpt[w].neighbor)
			if (!visitedCluster[k]) {
				visitedCluster[k] = true;
				if (SOpt[k].alpha >= alpha && SOpt[k].beta >= beta) {
					Q.push(k);
					F[k] = true;
				}
			}
		if (SOpt[w].alpha > k) {
			for (int i = 0; i < SOpt[w].leftGroup.size(); i++) {
				if (SOpt[w].leftGroup[i].first >= beta)
					break;
				for (auto j : SOpt[w].leftGroup[i].second)
					leftCom[j] = sum;
			}

			for (int i = 0; i < SOpt[w].rightGroup.size(); i++) {
				if (SOpt[w].rightGroup[i].first >= beta)
					break;
				for (auto j : SOpt[w].rightGroup[i].second)
					rightCom[j] = sum;
			}
		}
		else {
			for (int i = 0; i < SOpt[w].leftGroup.size(); i++) {
				if (SOpt[w].leftGroup[i].first >= alpha)
					break;
				for (auto j : SOpt[w].leftGroup[i].second)
					leftCom[j] = sum;
			}
			for (int i = 0; i < SOpt[w].rightGroup.size(); i++) {
				if (SOpt[w].rightGroup[i].first >= alpha)
					break;
				for (auto j : SOpt[w].rightGroup[i].second)
					rightCom[j] = sum;
			}
		}
	}
}

void abIndex::writeIndex() {
	int num = 0, snode = 0, tnode = 0;
	for (auto i : uNumber)
		num += i.size() * 4;
	for (auto i : vNumber)
		num += i.size() * 4;
	for (auto i : SOpt) {
		snode += 2 * 4;
		snode += i.neighbor.size() * 4;
		for (auto j : i.leftGroup) {
			snode += 4;
			snode += j.second.size() * 4;
		}
	}
	for (auto i : T) {
		tnode += 4 * 4;
		tnode += i.child.size() * 4;
		tnode += i.leaf.size() * 4;
	}
	ASG = (static_cast<long long>(num) + tnode + snode) / 1048576.0;
	ASGL = (static_cast<long long>(num) + snode) / 1048576.0;
	ASGT = (static_cast<long long>(tnode) + snode) / 1048576.0;

}

void abIndex::Maintenance(BiGraph& g) {
	this->g = &g;
	n1 = g.n1;
	n2 = g.n2;
	uNeighbor = g.uNeighbor;
	vNeighbor = g.vNeighbor;
}
void abIndex::initMain(BiGraph& g) {
	this->g = &g;
	n1 = g.n1;
	n2 = g.n2;
	uNeighbor = g.uNeighbor;
	vNeighbor = g.vNeighbor;
}
void abIndex::Maintenance(string adds, BiGraph& graph) {
	vector<pair<int, int>>edges;
	FILE* fp = fopen(adds.c_str(), "r");
	int u, v;
	while (fscanf(fp, "%ld %ld", &u, &v) != EOF)
		edges.push_back({ u,v });
	fclose(fp);
	auto start = chrono::system_clock::now();
	initMain(graph);
	int sss = 0;
	for (auto i : edges) {
		//cout << i.first << " " << i.second << endl;
		edgeInsert(i.first, i.second);
		sss++;
		if (sss == 1000)
			break;
	}
	creatSEviaVertex(nodeQ);
	addEdge();
	auto end = chrono::system_clock::now();
	chrono::duration<double> time = end - start;
	cout << "insert time: " << time.count() * 1000 << "ms" << endl;
}
void abIndex::Maintenance(string adds, BiGraph graph, int mup) {
	vector<pair<int, int>>edges;
	FILE* fp = fopen(adds.c_str(), "r");
	int u, v;
	while (fscanf(fp, "%ld %ld", &u, &v) != EOF)
		edges.push_back({ u,v });
	fclose(fp);
	auto start = chrono::system_clock::now();
	initMain(graph);
	int sss = 0;
	for (auto i : edges) {
		//cout << i.first << " " << i.second << endl;
		edgeInsert(i.first, i.second);
		sss++;
		if (sss == mup)
			break;
	}
	creatSEviaVertex(nodeQ);
	addEdge();
	auto end = chrono::system_clock::now();
	chrono::duration<double> time = end - start;
	cout << "insert time: " << time.count() * 1000 << "ms" << endl;
}
void abIndex::edgeInsert(int u, int v) {
	if (u >= uNeighbor.size())
		return;
	for (int i = 0; i < uNeighbor[u].size(); i++)
		if (uNeighbor[u][i] == v)
			return;
	g->addEdge(u, v);
	uNeighbor[u].push_back(v);
	vNeighbor[v].push_back(u);

	vector<pair<int, int>>H;
	//分析需要考虑的双核数
	vector<int> ubn = uNumber[u], vbn = vNumber[v];

	for (int i = 0; i < ubn.size(); i++) {
		int alpha = SOpt[ubn[i]].alpha;
		int beta = SOpt[ubn[i]].beta;
		int Bu = findOff(u, true, alpha + 1, true);
		int Bv = findOff(v, false, alpha + 1, true);
		if (Bv != 0 && Bv >= Bu)
			alphaIncrease(u, v, alpha, beta, H);
		Bv = findOff(v, false, alpha, true);
		if (Bv != 0 && Bv > beta)
			betaIncrease(u, v, alpha, beta, H);
	}
	for (int i = 0; i < vbn.size(); i++) {
		int alpha = SOpt[vbn[i]].alpha;
		int beta = SOpt[vbn[i]].beta;
		int Bu = findOff(u, true, beta + 1, false);
		int Bv = findOff(v, false, beta + 1, false);
		if (Bv != 0 && Bu >= Bv)
			alphaIncrease2(u, v, alpha, beta, H);
		Bu = findOff(v, false, beta, false);
		if (Bv != 0 && Bv > alpha)
			betaIncrease2(u, v, alpha, beta, H);
	}
	H.clear();
}
void abIndex::alphaIncrease(int u, int v, int alpha, int beta, vector<pair<int, int>>& H) {
	//top-beta
	vector<int>topX;
	for (int i = 0; i < uNeighbor[u].size(); i++)
		topX.push_back(findOff(uNeighbor[u][i], false, alpha + 1, true));
	sort(topX.begin(), topX.end(), [](int a, int b) {return a > b; });
	int Lu;
	if (alpha + 1 < topX.size())
		Lu = topX[alpha + 1];
	else
		return;
	//update
	bool had = false;
	for (int i = 0; i < H.size(); i++)
		if (H[i].first == alpha + 1 && H[i].second == Lu + 1) {
			had = true;
			break;
		}
	if (!had) {
		//bool flag = update(u, true, alpha + 1, Lu+1);
		neighborBaseIns(u, true, alpha + 1, Lu + 1);
		H.push_back(make_pair(alpha + 1, Lu + 1));
	}
}
void abIndex::betaIncrease(int u, int v, int alpha, int beta, vector<pair<int, int>>& H) {
	//top-beta
	vector<int>topX;
	for (int i = 0; i < uNeighbor[u].size(); i++)
		topX.push_back(findOff(uNeighbor[u][i], false, alpha, true));
	sort(topX.begin(), topX.end(), [](int a, int b) {return a > b; });
	int Lu;
	if (alpha < topX.size())
		Lu = topX[alpha];
	else
		return;
	//update
	bool had = false;
	for (int i = 0; i < H.size(); i++)
		if (H[i].first == alpha && H[i].second == Lu + 1) {
			had = true;
			break;
		}
	if (!had) {
		//bool flag = update(u, true, alpha, Lu+1);
		neighborBaseIns(u, true, alpha, Lu + 1);
		H.push_back(make_pair(alpha, Lu + 1));
	}
}
void abIndex::alphaIncrease2(int u, int v, int alpha, int beta, vector<pair<int, int>>& H) {
	//top-beta
	vector<int>topX;
	for (int i = 0; i < vNeighbor[v].size(); i++)
		topX.push_back(findOff(vNeighbor[v][i], true, beta, false));
	sort(topX.begin(), topX.end(), [](int a, int b) {return a > b; });
	int Lu;
	if (beta < topX.size())
		Lu = topX[beta];
	else
		return;
	//update
	bool had = false;
	for (int i = 0; i < H.size(); i++)
		if (H[i].first == Lu + 1 && H[i].second == beta) {
			had = true;
			break;
		}
	if (!had) {
		//bool flag = update(v, false, Lu+1, beta);
		neighborBaseIns2(v, false, Lu + 1, beta);
		H.push_back(make_pair(Lu + 1, beta));
	}
}
void abIndex::betaIncrease2(int u, int v, int alpha, int beta, vector<pair<int, int>>& H) {
	//top-beta
	vector<int>topX;
	for (int i = 0; i < vNeighbor[v].size(); i++)
		topX.push_back(findOff(vNeighbor[v][i], true, beta + 1, false));
	sort(topX.begin(), topX.end(), [](int a, int b) {return a > b; });
	int Lu;
	if (beta + 1 < topX.size())
		Lu = topX[beta + 1];
	else
		return;
	//update
	bool had = false;
	for (int i = 0; i < H.size(); i++)
		if (H[i].first == Lu + 1 && H[i].second == beta + 1) {
			had = true;
			break;
		}
	if (!had) {
		//bool flag = update(v, false, Lu+1, beta+1);
		neighborBaseIns2(v, false, Lu + 1, beta + 1);
		H.push_back(make_pair(Lu + 1, beta + 1));
	}
}
void abIndex::neighborBaseIns(int u, bool isLeft, int a, int b) {
	uT.clear(); uC.clear(); uDeg.clear();
	vT.clear(); vC.clear(); vDeg.clear();
	uT.resize(n1 + 1, 0); uC.resize(n1 + 1, 0); uDeg.resize(n1 + 1, 0);
	vT.resize(n2 + 1, 0); vC.resize(n2 + 1, 0); vDeg.resize(n2 + 1, 0);
	queue<pair<bool, int>>S;
	vector<pair<bool, int>>R;
	S.push(make_pair(isLeft, u));
	while (!S.empty()) {
		bool is = S.front().first;
		int w = S.front().second;
		R.push_back({ is,w });
		S.pop();
		queue<pair<bool, int>>temp = S;
		if (is) {
			uC[w] = true;
			uT[w] = true;
			for (auto i : uNeighbor[w]) {
				for (auto j : vNumber[i]) {
					if (!vC[i] && SOpt[j].alpha >= a && SOpt[j].beta >= b - 1) {
						uDeg[w]++;
						break;
					}
					if (!vT[i] && SOpt[j].alpha == a && SOpt[j].beta == b - 1) {
						uDeg[w]++;
						temp.push({ false,i });
						break;
					}
				}
			}
			if (uDeg[w] >= a)
				swap(S, temp);
			else
				deleteCand(w, true, a, b);
		}
		else {
			vC[w] = true;
			vT[w] = true;
			for (auto i : vNeighbor[w]) {
				for (auto j : uNumber[i]) {
					if (!uC[i] && SOpt[j].alpha >= a && SOpt[j].beta >= b - 1) {
						vDeg[w]++;
						break;
					}
					if (!uT[i] && SOpt[j].alpha == a && SOpt[j].beta == b - 1) {
						vDeg[w]++;
						temp.push({ false,i });
						break;
					}
				}
			}
			if (vDeg[w] >= b)
				swap(S, temp);
			else
				deleteCand(w, false, a, b);
		}
	}
	dealC(a, b, R);
}
void abIndex::neighborBaseIns2(int u, bool isLeft, int a, int b) {
	uT.clear(); uC.clear(); uDeg.clear();
	vT.clear(); vC.clear(); vDeg.clear();
	uT.resize(n1 + 1, 0); uC.resize(n1 + 1, 0); uDeg.resize(n1 + 1, 0);
	vT.resize(n2 + 1, 0); vC.resize(n2 + 1, 0); vDeg.resize(n2 + 1, 0);
	queue<pair<bool, int>>S;
	vector<pair<bool, int>>R;
	S.push(make_pair(isLeft, u));
	while (!S.empty()) {
		bool is = S.front().first;
		int w = S.front().second;
		R.push_back({ is,w });
		S.pop();
		queue<pair<bool, int>>temp = S;
		if (is) {
			uC[w] = true;
			uT[w] = true;
			for (auto i : uNeighbor[w]) {
				for (auto j : vNumber[i]) {
					if (!vC[i] && SOpt[j].alpha >= a - 1 && SOpt[j].beta >= b) {
						uDeg[w]++;
						break;
					}
					if (!vT[i] && SOpt[j].alpha == a - 1 && SOpt[j].beta == b) {
						uDeg[w]++;
						temp.push({ false,i });
						break;
					}
				}
			}
			if (uDeg[w] >= a)
				swap(S, temp);
			else
				deleteCand(w, true, a, b);
		}
		else {
			vC[w] = true;
			vT[w] = true;
			for (auto i : vNeighbor[w]) {
				for (auto j : uNumber[i]) {
					if (!uC[i] && SOpt[j].alpha >= a && SOpt[j].beta >= b - 1) {
						vDeg[w]++;
						break;
					}
					if (!uT[i] && SOpt[j].alpha == a && SOpt[j].beta == b - 1) {
						vDeg[w]++;
						temp.push({ false,i });
						break;
					}
				}
			}
			if (vDeg[w] >= b)
				swap(S, temp);
			else
				deleteCand(w, false, a, b);
		}
	}
	dealC(a, b, R);
}
void abIndex::deleteCand(int w, bool isLeft, int a, int b) {
	if (isLeft) {
		uC[w] = false;
		for (auto i : uNeighbor[w])
			if (vC[i]) {
				vDeg[i]--;
				if (vDeg[i] < b)
					deleteCand(i, false, a, b);
			}
	}
	else {
		vC[w] = false;
		for (auto i : vNeighbor[w])
			if (uC[i]) {
				uDeg[i]--;
				if (uDeg[i] < b)
					deleteCand(i, true, a, b);
			}
	}
}
bool abIndex::update(int u, bool isLeft, int node) {
	if (isLeft) {
		for (int i = 0; i < uNumber[u].size(); i++)
			if (SOpt[uNumber[u][i]].alpha >= SOpt[node].alpha
				&& SOpt[uNumber[u][i]].beta >= SOpt[node].beta)
				return false;
		uNumber[u].push_back(node);
		return true;
	}
	else {
		for (int i = 0; i < vNumber[u].size(); i++)
			if (SOpt[vNumber[u][i]].alpha >= SOpt[node].alpha
				&& SOpt[vNumber[u][i]].beta >= SOpt[node].beta)
				return false;
		vNumber[u].push_back(node);
		return true;
	}

}
int abIndex::findOff(int u, bool isLeft, int value, bool isAlpha) {
	int off = 0;
	if (isLeft) {
		if (isAlpha) {
			for (int i = 0; i < uNumber[u].size(); i += 2)
				if (SOpt[uNumber[u][i]].alpha >= value && SOpt[uNumber[u][i]].beta > off)
					off = SOpt[uNumber[u][i]].beta;
		}
		else {
			for (int i = 0; i < uNumber[u].size(); i += 2)
				if (SOpt[uNumber[u][i]].alpha >= off && SOpt[uNumber[u][i]].beta > value)
					off = SOpt[uNumber[u][i]].alpha;
		}
	}
	else {
		if (isAlpha) {
			for (int i = 0; i < vNumber[u].size(); i += 2)
				if (SOpt[vNumber[u][i]].alpha >= value && SOpt[vNumber[u][i]].beta > off)
					off = SOpt[vNumber[u][i]].beta;
		}
		else {
			for (int i = 0; i < vNumber[u].size(); i += 2)
				if (SOpt[vNumber[u][i]].alpha >= off && SOpt[vNumber[u][i]].beta > value)
					off = SOpt[vNumber[u][i]].alpha;
		}
	}
	return off;
}
void abIndex::dealC(int a, int b, vector<pair<bool, int>>R) {
	pair<vector<int>, vector<int>>C;
	int num = 0;
	for (int i = 0; i < R.size(); i++) {
		if (R[i].first) {
			if (uC[R[i].second]) {
				C.first.push_back(R[i].second);
				num++;
			}
		}
		else {
			if (vC[R[i].second]) {
				C.second.push_back(R[i].second);
				num++;
			}
		}
	}
	//cout << num << endl;
	int nodeNum = -1;
	if (C.first.size() > 0) {
		for (auto w : uNeighbor[C.first[0]])
			for (int i = 0; i < vNumber[w].size(); i++)
				if (SOpt[vNumber[w][i]].alpha == a && SOpt[vNumber[w][i]].beta == b) {
					nodeNum = vNumber[w][i];
					break;
				}
	}
	else if (C.second.size() > 0) {
		int nodeNum = -1;
		for (auto w : vNeighbor[C.second[0]])
			for (int i = 0; i < uNumber[w].size(); i++)
				if (SOpt[uNumber[w][i]].alpha == a && SOpt[uNumber[w][i]].beta == b) {
					nodeNum = uNumber[w][i];
					break;
				}
	}
	if (nodeNum == -1) {
		nodeNum = SOpt.size();
		node* newNode = new node;
		newNode->alpha = a;
		newNode->beta = b;
		SOpt.push_back(*newNode);
	}
	nodeQ.insert(nodeNum);
	SOpt[nodeNum].neighbor.clear();
	for (int i = 0; i < C.first.size(); i++) {
		update(C.first[i], true, nodeNum);
		uNumber[C.first[i]].push_back(nodeNum);
		SOpt[nodeNum].left.push_back(C.first[i]);
		for (auto j : uNumber[C.first[i]])
			if (SOpt[j].alpha <= a && SOpt[j].beta <= b) {
				SOpt[j].left.erase(std::remove(SOpt[j].left.begin(), SOpt[j].left.end(), C.first[i]), SOpt[j].left.end());
				nodeQ.insert(j);
				SOpt[j].neighbor.clear();
			}
	}
	for (int i = 0; i < C.second.size(); i++) {
		update(C.second[i], false, nodeNum);
		vNumber[C.second[i]].push_back(nodeNum);
		SOpt[nodeNum].right.push_back(C.second[i]);
		for (auto j : vNumber[C.second[i]])
			if (SOpt[j].alpha <= a && SOpt[j].beta <= b) {
				SOpt[j].right.erase(std::remove(SOpt[j].right.begin(), SOpt[j].right.end(), C.second[i]), SOpt[j].right.end());
				nodeQ.insert(j);
				SOpt[j].neighbor.clear();
			}
	}
}
void abIndex::creatSEviaVertex(set<int>& nodeQ) {
	int sum = 0;
	visited.resize(SOpt.size(), -1);
	rec = 0;
	block.resize(g->maxAlpha + 1);
	for (auto i : nodeQ) {
		for (int j = 0; j < SOpt[i].left.size(); j++) {
			int u = SOpt[i].left[j];
			for (int k = 0; k < g->uNeighbor[u].size(); k++) {
				int nbr = g->uNeighbor[u][k];
				for (int p = 0; p < vNumber[nbr].size(); p++)
					if (vNumber[nbr][p] > i && visited[vNumber[nbr][p]] != rec) {
						visited[vNumber[nbr][p]] = rec;
						int a = min(SOpt[i].alpha, SOpt[vNumber[nbr][p]].alpha);
						int b = min(SOpt[i].beta, SOpt[vNumber[nbr][p]].beta);
						if (block[a].size() <= b)
							block[a].resize(b + 1);
						block[a][b].insert(make_pair(i, vNumber[nbr][p]));
						sum++;
					}
			}
		}
		for (int j = 0; j < SOpt[i].right.size(); j++) {
			int v = SOpt[i].right[j];
			for (int k = 0; k < g->vNeighbor[v].size(); k++) {
				int nbr = g->vNeighbor[v][k];
				for (int p = 0; p < uNumber[nbr].size(); p++)
					if (uNumber[nbr][p] > i && visited[uNumber[nbr][p]] != rec) {
						visited[uNumber[nbr][p]] = rec;
						int a = min(SOpt[i].alpha, SOpt[uNumber[nbr][p]].alpha);
						int b = min(SOpt[i].beta, SOpt[uNumber[nbr][p]].beta);
						if (a < 1 || b < 1)
							continue;
						if (block[a].size() <= b)
							block[a].resize(b + 1);
						block[a][b].insert(make_pair(i, uNumber[nbr][p]));
						sum++;
					}
			}
		}
		rec++;
	}
	//cout << sum << endl;
}
void abIndex::creatSEviaNode_2(set<int>& nodeQ) {
	visited.resize(SOpt.size(), 0);
	block.resize(g->maxAlpha + 1);
	for (auto i : nodeQ) {
		visited[i] = 1;
		vector<int>c, f, non;
		for (auto j : SOpt[i].neighbor) {
			if (SOpt[i].alpha >= SOpt[j].alpha && SOpt[i].beta >= SOpt[j].beta)
				c.push_back(j);
			else if (SOpt[i].alpha <= SOpt[j].alpha && SOpt[i].beta <= SOpt[j].beta)
				f.push_back(j);
			else
				non.push_back(j);
		}
		if (f.size() >= 0) {
			for (int j = 1; j < f.size();j++) {
				int a =SOpt[i].alpha;
				int b = SOpt[i].beta;
				if (block[a].size() <= b)
					block[a].resize(b + 1);
				block[a][b].insert(make_pair(f[0], f[j]));
			}
			for (int j = 1; j < c.size(); j++) {
				int a = SOpt[i].alpha;
				int b = SOpt[i].beta;
				if (block[a].size() <= b)
					block[a].resize(b + 1);
				block[a][b].insert(make_pair(f[0], c[j]));
			}
			for (int j = 1; j < non.size(); j++) {
				int a = SOpt[i].alpha;
				int b = SOpt[i].beta;
				if (block[a].size() <= b)
					block[a].resize(b + 1);
				block[a][b].insert(make_pair(f[0], non[j]));
			}
		}
		else {
			for (int j = 1; j < non.size(); j++) {
				for (int k = 1; k < non.size(); k++) {
					int a = min(SOpt[non[k]].alpha, SOpt[non[j]].alpha);
					int b = min(SOpt[non[k]].beta, SOpt[non[j]].beta);
					if (block[a].size() <= b)
						block[a].resize(b + 1);
					block[a][b].insert(make_pair(non[k], non[j]));
				}
			}
		}

	}
}

void abIndex::edgeDelete(int u, int v) {
	int falg = 0;
	for (int i = 0; i < uNeighbor[u].size(); i++)
		if (uNeighbor[u][i] == v)
			falg = 1;
	if (falg == 0)
		return;
	g->deleteEdge(u, v);
	uNeighbor[u].erase(std::remove(uNeighbor[u].begin(), uNeighbor[u].end(), v), uNeighbor[u].end());
	vNeighbor[v].erase(std::remove(vNeighbor[v].begin(), vNeighbor[v].end(), u), vNeighbor[v].end());

	//分析需要考虑的双核数
	vector<int>ubn = uNumber[u], vbn = vNumber[v];
	for (int i = 0; i < ubn.size(); i++) {
		int alpha = SOpt[ubn[i]].alpha;
		int beta = SOpt[ubn[i]].beta;
		int Bv = findOff(v, false, alpha, true);
		if (Bv < beta)
			continue;
		vector<pair<bool, int>> R = neighborBaseDel(u, true, alpha, beta);
		alphaDecrease(u, v, alpha, beta, R);
		betaDecrease(u, v, alpha, beta, R);
	}
	for (int i = 0; i < vbn.size(); i++) {
		int alpha = SOpt[vbn[i]].alpha;
		int beta = SOpt[vbn[i]].beta;
		int Bu = findOff(u, true, beta, false);
		if (Bu < alpha)
			continue;
		vector<pair<bool, int>> R = neighborBaseDel(v, false, alpha, beta);
		alphaDecrease2(u, v, alpha, beta, R);
		betaDecrease2(u, v, alpha, beta, R);
	}
}
void abIndex::alphaDecrease(int u, int v, int alpha, int beta, vector<pair<bool, int>>& R) {
	//top-beta
	vector<int>topX;
	for (int i = 0; i < uNeighbor[u].size(); i++)
		topX.push_back(findOff(uNeighbor[u][i], false, alpha - 1, true));
	sort(topX.begin(), topX.end(), [](int a, int b) {return a > b; });
	int Lu;
	if (alpha - 1 < topX.size())
		Lu = topX[alpha - 1];
	else
		return;
	//update
	if (R.size() == 0)
		return;
	dealC2(alpha - 1, Lu, R);
}
void abIndex::betaDecrease(int u, int v, int alpha, int beta, vector<pair<bool, int>>& R) {
	//top-beta
	vector<int>topX;
	for (int i = 0; i < uNeighbor[u].size(); i++)
		topX.push_back(findOff(uNeighbor[u][i], false, alpha, true));
	sort(topX.begin(), topX.end(), [](int a, int b) {return a > b; });
	int Lu;
	if (alpha < topX.size())
		Lu = topX[alpha];
	else
		return;

	//update
	if (R.size() == 0)
		return;
	dealC2(alpha, beta - 1, R);
}
void abIndex::alphaDecrease2(int u, int v, int alpha, int beta, vector<pair<bool, int>>& R) {
	//top-beta
	vector<int>topX;
	for (int i = 0; i < vNeighbor[v].size(); i++)
		topX.push_back(findOff(vNeighbor[v][i], true, beta, false));
	sort(topX.begin(), topX.end(), [](int a, int b) {return a > b; });
	int Lu;
	if (beta < topX.size())
		Lu = topX[beta];
	else
		return;
	//update
	if (R.size() == 0)
		return;
	dealC2(alpha - 1, beta, R);
}
void abIndex::betaDecrease2(int u, int v, int alpha, int beta, vector<pair<bool, int>>& R) {
	//top-beta
	vector<int>topX;
	for (int i = 0; i < vNeighbor[v].size(); i++)
		topX.push_back(findOff(vNeighbor[v][i], true, beta - 1, false));
	sort(topX.begin(), topX.end(), [](int a, int b) {return a > b; });
	int Lu;
	if (beta - 1 < topX.size())
		Lu = topX[beta - 1];
	else
		return;

	//update
	if (R.size() == 0)
		return;
	dealC2(Lu, beta - 1, R);
}
vector<pair<bool, int>> abIndex::neighborBaseDel(int u, bool isLeft, int a, int b) {
	uT.clear(); uC.clear(); uDeg.clear();
	vT.clear(); vC.clear(); vDeg.clear();
	uT.resize(n1 + 1, 0); uC.resize(n1 + 1, 0); uDeg.resize(n1 + 1, 0);
	vT.resize(n2 + 1, 0); vC.resize(n2 + 1, 0); vDeg.resize(n2 + 1, 0);
	queue<pair<bool, int>>S;
	vector<pair<bool, int>>R;
	S.push(make_pair(isLeft, u));
	while (!S.empty()) {
		bool is = S.front().first;
		int w = S.front().second;
		S.pop();
		queue<pair<bool, int>>temp = S;
		vector<bool>uS(n1 + 1, 0), vS(n2 + 1, 0);
		if (is) {
			uT[w] = true;
			for (auto i : uNeighbor[w])
				for (int j = 0; j < vNumber[i].size(); j++)
					if (SOpt[vNumber[i][j]].alpha >= a && SOpt[vNumber[i][j]].beta >= b) {
						uDeg[w]++;
						break;
					}
			if (uDeg[w] < a) {
				addCand(R, w, true, a, b, temp, uS, vS);
				swap(S, temp);
			}
		}
		else {
			vT[w] = true;
			for (auto i : vNeighbor[w])
				for (int j = 0; j < uNumber[i].size(); j += 2)
					if (SOpt[uNumber[i][j]].alpha >= a && SOpt[uNumber[i][j]].beta >= b) {
						vDeg[w]++;
						break;
					}
			if (vDeg[w] < b) {
				addCand(R, w, false, a, b, temp, uS, vS);
				swap(S, temp);
			}
		}
	}
	return R;
}
void abIndex::addCand(vector<pair<bool, int>>& R, int w, bool isLeft, int a, int b, queue<pair<bool, int>>& S, vector<bool>& uS, vector<bool>& vS) {
	if (isLeft) {
		uC[w] = true;
		R.push_back({ true,w });
		for (auto i : uNeighbor[w]) {
			if (!vT[i] && findOff(i, false, a, true) == b && !vS[i]) {
				S.push({ false,i });
			}
			if (vT[i] && !vC[i]) {
				vDeg[i]--;
				if (vDeg[i] < b)
					addCand(R, i, false, a, b, S, uS, vS);
			}
		}
	}
	else {
		vC[w] = true;
		R.push_back({ false,w });
		for (auto i : vNeighbor[w]) {
			if (!uT[i] && findOff(i, true, a, true) == b && !uS[i]) {
				S.push({ true,i });
			}
			if (uT[i] && !uC[i]) {
				uDeg[i]--;
				if (uDeg[i] < a)
					addCand(R, i, true, a, b, S, uS, vS);
			}
		}
	}
}
void abIndex::dealC2(int a, int b, vector<pair<bool, int>>R) {
	pair<vector<int>, vector<int>>C;
	int num = 0;
	for (int i = 0; i < R.size(); i++) {
		if (R[i].first) {
			if (uC[R[i].second]) {
				C.first.push_back(R[i].second);
				num++;
			}
		}
		else {
			if (vC[R[i].second]) {
				C.second.push_back(R[i].second);
				num++;
			}
		}
	}
	//cout << num << endl;
	int nodeNum = -1;
	if (C.first.size() > 0) {
		for (auto w : uNeighbor[C.first[0]])
			for (int i = 0; i < vNumber[w].size(); i++)
				if (SOpt[vNumber[w][i]].alpha == a && SOpt[vNumber[w][i]].beta == b) {
					nodeNum = vNumber[w][i];
					break;
				}
	}
	else if (C.second.size() > 0) {
		int nodeNum = -1;
		for (auto w : vNeighbor[C.second[0]])
			for (int i = 0; i < uNumber[w].size(); i++)
				if (SOpt[uNumber[w][i]].alpha == a && SOpt[uNumber[w][i]].beta == b) {
					nodeNum = uNumber[w][i];
					break;
				}
	}
	if (nodeNum == -1) {
		nodeNum = SOpt.size();
		node* newNode = new node;
		newNode->alpha = a;
		newNode->beta = b;
		SOpt.push_back(*newNode);
	}
	nodeQ.insert(nodeNum);
	SOpt[nodeNum].neighbor.clear();
	for (int i = 0; i < C.first.size(); i++) {
		update(C.first[i], true, nodeNum);
		uNumber[C.first[i]].push_back(nodeNum);
		SOpt[nodeNum].left.push_back(C.first[i]);
		for (auto j : uNumber[C.first[i]])
			if (SOpt[j].alpha >= a && SOpt[j].beta >= b) {
				SOpt[j].left.erase(std::remove(SOpt[j].left.begin(), SOpt[j].left.end(), C.first[i]), SOpt[j].left.end());
				nodeQ.insert(j);
				SOpt[j].neighbor.clear();
			}
	}
	for (int i = 0; i < C.second.size(); i++) {
		update(C.second[i], false, nodeNum);
		vNumber[C.second[i]].push_back(nodeNum);
		SOpt[nodeNum].right.push_back(C.second[i]);
		for (auto j : vNumber[C.second[i]])
			if (SOpt[j].alpha >= a && SOpt[j].beta >= b) {
				SOpt[j].right.erase(std::remove(SOpt[j].right.begin(), SOpt[j].right.end(), C.second[i]), SOpt[j].right.end());
				nodeQ.insert(j);
				SOpt[j].neighbor.clear();
			}
	}
}
