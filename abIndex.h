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
	vector<vector<set<pair<int, int>>>>block;
	vector<vector<pair<int, int>>>rank;
	vector<int>wmap;

	double ASG, ASGL, ASGT;
public:
	abIndex();
	abIndex(BiGraph& graph);
	void creatSG(BiGraph& graph);
	void creatSG1(BiGraph& graph);
	void creatSG2(BiGraph& graph);
	void creatSGOpt(BiGraph& graph);
	void creatTree();
	void addEdge();
	void BFSonASG(int a, int b, int w, vector<int>& wmap);
	void query(int q, int alpha, int beta, vector<bool>& leftResult, vector<bool>& rightResult);
	void queryOpt(int q, int alpha, int beta, vector<bool>& leftResult, vector<bool>& rightResult, int k);
	void queryCore(int alpha, int beta, vector<int>& leftCom, vector<int>& rightCom, int k);
	void BFSonASG(int s, int alpha, int beta, vector<int>& leftCom, vector<int>& rightCom, vector<bool>& F, int k, int sum);
	void writeIndex();
};


abIndex::abIndex() {}
abIndex::abIndex(BiGraph& graph) {
	uNumber = graph.uNumber;
	graph.uNumber.clear();
	vNumber = graph.vNumber;
	graph.vNumber.clear();
	//creatSG(graph);
	//creatSGOpt(graph);

	///cout << S.size() << "  " << SOpt.size() << endl;
	//creatTree();
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