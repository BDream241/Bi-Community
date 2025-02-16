#include"ASG.h"

ASG::ASG() {}
ASG::ASG(BiGraph& g) {
	uIndex.resize(g.n1 + 1);
	vIndex.resize(g.n2 + 1);
	addNode(g);
	visited.resize(asg.size(), -1);
	rec = 0;
	block.resize(g.maxAlpha + 1);
	creatSEviaNode(g);
	addEdge();
}
void ASG::addNode(BiGraph& g) {
	vector<vector<vector<pair<int, int>>>>si(g.maxAlpha + 1);
	for (int i = 1; i <= g.n1; i++)
		for (auto j : g.uNumber[i]) {
			int a = get<0>(j);
			int b = get<1>(j);
			int c = get<2>(j);
			if (si[a].size() <= b)
				si[a].resize(b + 1);
			int f = -1;
			for (auto k : si[a][b])
				if (k.first == c) {
					f = k.second;
					break;
				}
			if (f == -1) {
				Node node;
				node.alpha = a;
				node.beta = b;
				node.left.push_back(i);
				asg.push_back(node);
				si[a][b].push_back({ c,asg.size() - 1 });
				uIndex[i].push_back(asg.size() - 1);
			}
			else {
				asg[f].left.push_back(i);
				uIndex[i].push_back(f);
			}
		}
	for (int i = 1; i <= g.n2; i++)
		for (auto j : g.vNumber[i]) {
			int a = get<0>(j);
			int b = get<1>(j);
			int c = get<2>(j);
			if (si[a].size() <= b)
				si[a].resize(b + 1);
			int f = -1;
			for (auto k : si[a][b])
				if (k.first == c) {
					f = k.second;
					break;
				}
			if (f == -1) {
				Node node;
				node.alpha = a;
				node.beta = b;
				node.right.push_back(i);
				asg.push_back(node);
				si[a][b].push_back({ c,asg.size() - 1 });
				vIndex[i].push_back(asg.size() - 1);
			}
			else {
				asg[f].right.push_back(i);
				vIndex[i].push_back(f);
			}
		}
	si.clear();
}
void ASG::creatSEviaNode(BiGraph& g) {
	int sum = 0;
	for (int i = 0; i < asg.size(); i++) {
		for (int j = 0; j < asg[i].left.size(); j++) {
			int u = asg[i].left[j];
			for (int k = 0; k < g.uNeighbor[u].size(); k++) {
				int nbr = g.uNeighbor[u][k];
				for (int p = 0; p < vIndex[nbr].size(); p++)
					if (vIndex[nbr][p] > i && visited[vIndex[nbr][p]] != rec) {
						visited[vIndex[nbr][p]] = rec;
						int a = min(asg[i].alpha, asg[vIndex[nbr][p]].alpha);
						int b = min(asg[i].beta, asg[vIndex[nbr][p]].beta);
						if (block[a].size() <= b)
							block[a].resize(b + 1);
						block[a][b].push_back(make_pair(i, vIndex[nbr][p]));
						sum++;
					}
			}
		}
		for (int j = 0; j < asg[i].right.size(); j++) {
			int v = asg[i].right[j];
			for (int k = 0; k < g.vNeighbor[v].size(); k++) {
				int nbr = g.vNeighbor[v][k];
				for (int p = 0; p < uIndex[nbr].size(); p++)
					if (uIndex[nbr][p] > i && visited[uIndex[nbr][p]] != rec) {
						visited[uIndex[nbr][p]] = rec;
						int a = min(asg[i].alpha, asg[uIndex[nbr][p]].alpha);
						int b = min(asg[i].beta, asg[uIndex[nbr][p]].beta);
						if (block[a].size() <= b)
							block[a].resize(b + 1);
						block[a][b].push_back(make_pair(i, uIndex[nbr][p]));
						sum++;
					}
			}
		}
		rec++;
	}
	//cout << sum << endl;
}
void ASG::addEdge() {
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
			vector<int>wmap(asg.size(), 0);
			for (int k = 0; k < block[a][b].size(); k++) {
				if (wmap[block[a][b][k].first] == 0) {
					wmap[block[a][b][k].first] = uf.add();
					BFSonASG(a, b, block[a][b][k].first, wmap);
				}
				if (wmap[block[a][b][k].second] == 0) {
					wmap[block[a][b][k].second] = uf.add();
					BFSonASG(a, b, block[a][b][k].second, wmap);
				}
				
				if (uf.merge(uf.Find(wmap[block[a][b][k].first]), uf.Find(wmap[block[a][b][k].second]))) {
					asg[block[a][b][k].first].neighbor.push_back(block[a][b][k].second);
					asg[block[a][b][k].second].neighbor.push_back(block[a][b][k].first);
				}
			}
			wmap.clear();
		}
	block.clear();
}
void ASG::BFSonASG(int a,int b, int w, vector<int>&wmap) {
	rec++;
	queue<int>bfs;
	bfs.push(w);
	visited[w] = rec;
	while (!bfs.empty()) {
		int p = bfs.front();
		bfs.pop();
		for (auto nbr : asg[p].neighbor)
			if (visited[nbr] != rec) {
				visited[nbr] = rec;
				if (asg[nbr].alpha >= a && asg[nbr].beta >= b) {
					bfs.push(nbr);
					wmap[nbr] = wmap[w];
				}
			}
	}
}
void ASG::Query(int alpha, int beta, int u, bool isU, vector<bool>& leftResult, vector<bool>& rightResult) {
	vector<bool>had(asg.size(), false);
	queue<int>Q;
	if (isU) {
		for (auto i : uIndex[u]) {
			had[i] = true;
			if (asg[i].alpha >= alpha && asg[i].beta >= beta) {
				Q.push(i);
				for (auto j : asg[i].left)
					leftResult[j] = true;
				for (auto j : asg[i].right)
					rightResult[j] = true;
			}
		}
	}
	else {
		for (auto i : vIndex[u]) {
			had[i] = true;
			if (asg[i].alpha >= alpha && asg[i].beta >= beta) {
				Q.push(i);
				for (auto j : asg[i].left)
					leftResult[j] = true;
				for (auto j : asg[i].right)
					rightResult[j] = true;
			}
		}
	}
	while (!Q.empty()) {
		int w = Q.front();
		Q.pop();
		for(auto i:asg[w].neighbor)
			if (!had[i]) {
				had[i] = true;
				if (asg[i].alpha >= alpha && asg[i].beta >= beta) {
					Q.push(i);
					for (auto j : asg[i].left)
						leftResult[j] = true;
					for (auto j : asg[i].right)
						rightResult[j] = true;
				}
			}
	}
}
int verifyCom(BiGraph& g, vector<bool>& leftResult, vector<bool>& rightResult, int alpha, int beta, int q, bool isLeft) {
	vector<bool>leftcom(g.n1 + 1, false), rightcom(g.n2 + 1, false);
	vector<int>tempLeft(g.n1 + 1), tempRight(g.n2 + 1);
	for (int i = 0; i <= g.n1; i++)
		tempLeft[i] = g.uNeighbor[i].size();
	for (int i = 0; i <= g.n2; i++)
		tempRight[i] = g.vNeighbor[i].size();
	
	vector<bool>leftVerify(g.n1 + 1, true);
	vector<bool>rightVerify(g.n2 + 1, true);
	vector<int>leftQ, rightQ;
	for (int i = 1; i < tempLeft.size(); i++)
		if (tempLeft[i] < alpha) {
			leftQ.push_back(i);
			leftVerify[i] = false;
		}
	for (int i = 1; i < tempRight.size(); i++)
		if (tempRight[i] < beta) {
			rightQ.push_back(i);
			rightVerify[i] = false;
		}
	while (!leftQ.empty() || !rightQ.empty()) {
		for (int i = 0; i < leftQ.size(); i++) {
			int u = leftQ[i];
			for (int j = 0; j < g.uNeighbor[u].size(); j++) {
				int v = g.uNeighbor[u][j];
				if (rightVerify[v]) {
					tempRight[v]--;
					if (tempRight[v] < beta) {
						rightQ.push_back(v);
						rightVerify[v] = false;
					}
				}
			}
		}
		leftQ.clear();
		for (int i = 0; i < rightQ.size(); i++) {
			int v = rightQ[i];
			for (int j = 0; j < g.vNeighbor[v].size(); j++) {
				int u = g.vNeighbor[v][j];
				if (leftVerify[u]) {
					tempLeft[u]--;
					if (tempLeft[u] < alpha) {
						leftQ.push_back(u);
						leftVerify[u] = false;
					}
				}
			}
		}
		rightQ.clear();
	}
	if (isLeft && leftVerify[q]) {
		leftcom[q] = true;
		leftQ.push_back(q);
	}
	else if (!isLeft && rightVerify[q]) {
		rightcom[q] = true;
		rightQ.push_back(q);
	}
	while (!leftQ.empty() || !rightQ.empty()) {
		for (int i = 0; i < leftQ.size(); i++) {
			int u = leftQ[i];
			for (int j = 0; j < g.uNeighbor[u].size(); j++) {
				int v = g.uNeighbor[u][j];
				if (rightVerify[v] && !rightcom[v]) {
					rightQ.push_back(v);
					rightcom[v] = true;
				}
			}
		}
		leftQ.clear();
		for (int i = 0; i < rightQ.size(); i++) {
			int v = rightQ[i];
			for (int j = 0; j < g.vNeighbor[v].size(); j++) {
				int u = g.vNeighbor[v][j];
				if (leftVerify[u] && !leftcom[u]) {
					leftQ.push_back(u);
					leftcom[u] = true;
				}
			}
		}
		rightQ.clear();
	}
	int more = 0, lower = 0;
	for (int i = 1; i <= g.n1; i++)
		if (leftcom[i] != leftResult[i]) {
			leftcom[i] ? lower++ : more++;
		}
	for (int i = 1; i <= g.n2; i++)
		if (rightcom[i] != rightResult[i]) {
			rightcom[i] ? lower++ : more++;
		}
	if ((more + lower) == 0) {
		cout << "Verify: true" << endl;
		return 1;
	}
	else {
		cout << "Verify: false   more: " << more << "  lower: " << lower ;
		cout << "  Input:" << alpha << "  " << beta << "  " << q << endl;
		return 0;
	}
}

void ASG::creatSEviaNode(BiGraph& g, set<int>& nodeQ) {
	int sum = 0;
	visited.resize(asg.size(), -1);
	rec = 0;
	block.resize(g.maxAlpha + 1);
	for (auto i : nodeQ) {
		for (int j = 0; j < asg[i].left.size(); j++) {
			int u = asg[i].left[j];
			for (int k = 0; k < g.uNeighbor[u].size(); k++) {
				int nbr = g.uNeighbor[u][k];
				for (int p = 0; p < vIndex[nbr].size(); p++)
					if (vIndex[nbr][p] > i && visited[vIndex[nbr][p]] != rec) {
						visited[vIndex[nbr][p]] = rec;
						int a = min(asg[i].alpha, asg[vIndex[nbr][p]].alpha);
						int b = min(asg[i].beta, asg[vIndex[nbr][p]].beta);
						if (block[a].size() <= b)
							block[a].resize(b + 1);
						block[a][b].push_back(make_pair(i, vIndex[nbr][p]));
						sum++;
					}
			}
		}
		for (int j = 0; j < asg[i].right.size(); j++) {
			int v = asg[i].right[j];
			for (int k = 0; k < g.vNeighbor[v].size(); k++) {
				int nbr = g.vNeighbor[v][k];
				for (int p = 0; p < uIndex[nbr].size(); p++)
					if (uIndex[nbr][p] > i && visited[uIndex[nbr][p]] != rec) {
						visited[uIndex[nbr][p]] = rec;
						int a = min(asg[i].alpha, asg[uIndex[nbr][p]].alpha);
						int b = min(asg[i].beta, asg[uIndex[nbr][p]].beta);
						if (a < 1 || b < 1)
							continue;
						if (block[a].size() <= b)
							block[a].resize(b + 1);
						block[a][b].push_back(make_pair(i, uIndex[nbr][p]));
						sum++;
					}
			}
		}
		rec++;
	}
	//cout << sum << endl;
}


