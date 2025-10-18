#pragma once
#include"utility.h"

class BiGraph {
public:
	int n1, n2, m, maxAlpha, maxBeta, maxK, maxUDeg, maxVDeg;
	vector<int>uDeg, vDeg;
	vector<vector<int>> uNeighbor, vNeighbor;

	vector<tuple<int, int, int>>cluster;
	map<tuple<int, int, int>, int>cMap;
	vector<vector<int>> uNumber, vNumber, uBeta, vBeta;

	unordered_map<int, int>umap, vmap;
	//vector<vector<set<pair<int, int>>>>block;
	set<pair<int, int>> block1, block2;
	unordered_set<int>gamSet;
	vector<tuple<int, bool, set<pair<bool, int>>>>tempSet;

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
	int decompose();
	int coreDecompose();
	int alphaDecompose(int alpha, vector<int>& leftDeg, vector<int>& rightDeg, vector<bool>& uDelete,
		vector<bool>& vDelete, vector<vector<int>>& uDAG,
		vector<vector<int>>& vDAG, vector<vector<pair<int, bool>>>& order);
	int betaDecompose(int beta, vector<int>& leftDeg, vector<int>& rightDeg, vector<bool>& uDelete,
		vector<bool>& vDelete, vector<vector<int>>& uDAG,
		vector<vector<int>>& vDAG, vector<vector<pair<int, bool>>>& order);
	int comBN(pair<int, bool>ab, vector<vector<int>>& uDAG, vector<vector<int>>& vDAG,
		vector<vector<pair<int, bool>>>& order);
	int updateNumber(int d, bool isAlpha, int p, bool isLeft, int x);

public:
	int verifyCom(BiGraph& g, vector<bool>& leftResult, vector<bool>& rightResult, int alpha, int beta, int q, bool isLeft);
	int verifyCore(BiGraph& g, vector<int>& leftResult, vector<int>& rightResult, int alpha, int beta);
};


BiGraph::BiGraph() {}
BiGraph::BiGraph(string str) {
	n1 = 0;
	n2 = 0;
	m = 0;
	loadGraph(str);
	cout << "n1: " << n1 << " n2: " << n2 << endl;
	cout << "m: " << m << endl;
}

void BiGraph::addEdge(int u, int v) {
	m++;
	uNeighbor[u].push_back(v);
	vNeighbor[v].push_back(u);
	uDeg[u]++;
	vDeg[v]++;
	if (maxUDeg < uDeg[u])
		maxUDeg = uDeg[u];
	if (maxVDeg < vDeg[v])
		maxVDeg = vDeg[v];
}

void BiGraph::deleteEdge(int u, int v) {
	m--;
	auto it = remove(uNeighbor[u].begin(), uNeighbor[u].end(), v);
	uNeighbor[u].erase(it, uNeighbor[u].end());
	uDeg[u]--;
	it = remove(vNeighbor[v].begin(), vNeighbor[v].end(), u);
	vNeighbor[v].erase(it, vNeighbor[v].end());
	vDeg[v]--;
}

void BiGraph::loadGraph(string str) {
	string gstr = str + "/graph.txt";
	string estr = str + "/edge.txt";
	FILE* g = fopen(gstr.c_str(), "r");
	FILE* e = fopen(estr.c_str(), "r");
	fscanf(g, "%lld %lld %lld", &n1, &n2, &m);
	init(n1, n2);
	int u, v;
	while ((fscanf(e, "%d %d", &u, &v)) != EOF)
		addEdge(u, v);
	fclose(g);
	fclose(e);
}

void BiGraph::init(int n1, int n2) {
	this->n1 = n1;
	this->n2 = n2;
	m = 0;

	maxUDeg = 0;
	maxVDeg = 0;
	uNeighbor.resize(n1 + 1);
	vNeighbor.resize(n2 + 1);
	uDeg.resize(n1 + 1, 0);
	vDeg.resize(n2 + 1, 0);
	uNumber.resize(n1 + 1);
	vNumber.resize(n2 + 1);
	uBeta.resize(n1 + 1);
	vBeta.resize(n2 + 1);
}

void BiGraph::print() {
	for (int i = 1; i <= n1; i++) {
		cout << i << " : ";
		for (int j = 0; j < uNeighbor[i].size(); j++)
			cout << uNeighbor[i][j] << " ";
		cout << endl;
	}
	for (int i = 1; i <= n2; i++) {
		cout << i << " : ";
		for (int j = 0; j < vNeighbor[i].size(); j++)
			cout << vNeighbor[i][j] << " ";
		cout << endl;
	}
}


int BiGraph::coreDecompose() {
	vector<int> leftQ;
	vector<int> rightQ;

	int num1 = n1 + 1;
	vector<int> leftR(num1);
	for (int i = 0; i < leftR.size(); i++)
		leftR[i] = i;

	int leftRTnum = 0;
	vector<int> leftRT(num1);

	int num2 = n2 + 1;//
	vector<int> rightR(num2);
	for (int i = 0; i < rightR.size(); i++)
		rightR[i] = i;

	int rightRTnum = 0;
	vector<int> rightRT(num2);

	vector<bool>leftDel(num1, false), rightDel(num2, false);
	vector<int>leftDeg = uDeg, rightDeg = vDeg;

	int kc = 1;
	for (kc = 1; kc <= maxVDeg + 1; kc++) {
		bool stop = true;
		leftRTnum = 0;
		for (int i = 0; i < num1; i++) {
			int u = leftR[i];
			if (!leftDel[u]) {
				stop = false;
				leftRT[leftRTnum] = u;
				leftRTnum++;
				if (leftDeg[u] < kc) {
					leftQ.push_back(u);
				}
			}
		}
		swap(leftR, leftRT);
		num1 = leftRTnum;
		if (stop)
			break;
		stop = true;
		rightRTnum = 0;
		for (int i = 0; i < num2; i++) {
			int v = rightR[i];
			if (!rightDel[v]) {
				stop = false;
				rightRT[rightRTnum] = v;
				rightRTnum++;
				if (rightDeg[v] < kc) {
					rightQ.push_back(v);
				}
			}
		}
		swap(rightR, rightRT);
		num2 = rightRTnum;
		if (stop)
			break;
		while (!leftQ.empty() || !rightQ.empty()) {

			for (auto j = leftQ.begin(); j != leftQ.end(); j++) {
				int u = *j;
				if (leftDel[u])
					continue;
				for (int k = 0; k < uNeighbor[u].size(); k++) {
					int v = uNeighbor[u][k];
					if (rightDel[v])
						continue;
					rightDeg[v]--;

					if (rightDeg[v] == 0) {
						rightDel[v] = true;
					}
					if (rightDeg[v] < kc) {
						rightQ.push_back(v);
					}
				}
				leftDeg[u] = 0;
				leftDel[u] = true;
			}
			leftQ.clear();

			for (auto j = rightQ.begin(); j != rightQ.end(); j++) {
				int v = *j;
				if (rightDel[v])
					continue;
				for (int k = 0; k < vNeighbor[v].size(); k++) {
					int u = vNeighbor[v][k];
					if (leftDel[u]) continue;

					leftDeg[u]--;
					if (leftDeg[u] == 0) {
						leftDel[u] = true;
					}
					if (leftDeg[u] < kc) {
						leftQ.push_back(u);
					}
				}
				rightDeg[v] = 0;
				rightDel[v] = true;
			}
			rightQ.clear();
		}
	}
	return kc - 2;
}
int BiGraph::decompose() {
	vector<bool>uDel(n1 + 1, false), vDel(n2 + 1, false);
	vector<vector<int>>uDAG, vDAG;
	vector<vector<pair<int, bool>>> order;


	maxK = coreDecompose();
	cout << "max k-core:" << maxK << endl;
	auto start = chrono::system_clock::now();
	auto end = chrono::system_clock::now();
	chrono::duration<double> time = end - start;

	for (int alpha = 1; alpha <= maxK; alpha++) {
		alphaDecompose(alpha, uDeg, vDeg, uDel, vDel, uDAG, vDAG, order);

		comBN({ alpha,true }, uDAG, vDAG, order);
		/*if (alpha % 10 == 0)
			cout << alpha << endl;
		else
			*/
		if (alpha == 1)
			maxBeta = order.size();
		uDAG.clear();
		vDAG.clear();
		order.clear();
	}
	for (int beta = 1; beta <= maxK; beta++) {
		betaDecompose(beta, uDeg, vDeg, uDel, vDel, uDAG, vDAG, order);

		comBN({ beta,false }, uDAG, vDAG, order);
		/*if (beta % 10 == 0)
			cout << beta << endl;
		else */
		if (beta == 1)
			maxAlpha = order.size() + maxK;
		uDAG.clear();
		vDAG.clear();
		order.clear();
	}

	for (int i = 1; i <= n1; i++)
		if (uBeta[i].size() > 0) {
			int na = uNumber[i][uNumber[i].size() - 1];
			int nb = uBeta[i][uBeta[i].size() - 1];
			if (get<1>(cluster[na]) == get<1>(cluster[nb])) {
				block2.insert({ nb,uNumber[i][uNumber[i].size() - 1] });
				uNumber[i][uNumber[i].size() - 1] = nb;
			}
			else
				uNumber[i].push_back(uBeta[i][uBeta[i].size() - 1]);
			for (int j = uBeta[i].size() - 2; j >= 0; j--)
				uNumber[i].push_back(uBeta[i][j]);
		}
	for (int i = 1; i <= n2; i++)
		if (vBeta[i].size() > 0) {
			int na = vNumber[i][vNumber[i].size() - 1];
			int nb = vBeta[i][vBeta[i].size() - 1];
			if (get<1>(cluster[na]) == get<1>(cluster[nb])) {
				block2.insert({ nb,vNumber[i][vNumber[i].size() - 1] });
				vNumber[i][vNumber[i].size() - 1] = nb;
			}
			else
				vNumber[i].push_back(vBeta[i][vBeta[i].size() - 1]);
			for (int j = vBeta[i].size() - 2; j >= 0; j--)
				vNumber[i].push_back(vBeta[i][j]);
		}
	return 0;
}
int BiGraph::alphaDecompose(int alpha, vector<int>& leftDeg, vector<int>& rightDeg, vector<bool>& uDelete,
	vector<bool>& vDelete, vector<vector<int>>& uDAG,
	vector<vector<int>>& vDAG, vector<vector<pair<int, bool>>>& order) {

	vector<int>leftQ, rightQ;

	for (int i = 1; i < uDelete.size(); i++)
		if (!uDelete[i] && leftDeg[i] < alpha) {
			leftQ.push_back(i);
			uDelete[i] = true;
		}
	while (!leftQ.empty() || !rightQ.empty()) {
		for (int i = 0; i < leftQ.size(); i++) {
			int u = leftQ[i];
			for (int j = 0; j < uNeighbor[u].size(); j++) {
				int v = uNeighbor[u][j];
				if (vDelete[v])
					continue;
				rightDeg[v]--;
				if (rightDeg[v] == 0) {
					rightQ.push_back(v);
					vDelete[v] = true;
				}
			}
		}
		leftQ.clear();
		for (int i = 0; i < rightQ.size(); i++) {
			int v = rightQ[i];
			for (int j = 0; j < vNeighbor[v].size(); j++) {
				int u = vNeighbor[v][j];
				if (uDelete[u])
					continue;
				leftDeg[u]--;
				if (leftDeg[u] < alpha) {
					leftQ.push_back(u);
					uDelete[u] = true;
				}
			}
		}
		rightQ.clear();
	}
	//初始化
	int num = n2;
	int nextNum = 0;
	vector<pair<int, bool>>bfsQ;
	vector<int>remain(num);
	vector<int>nextRemain(num);
	vector<bool>leftDel = uDelete;
	vector<bool>rightDel = vDelete;
	vector<int>degL = leftDeg;
	vector<int>degR = rightDeg;

	//bi-core number
	uDAG.resize(n1 + 1);
	vDAG.resize(n2 + 1);

	for (int i = 0; i < remain.size(); i++)
		remain[i] = i + 1;
	for (int beta = 1; beta <= maxVDeg + 1; beta++) {
		vector<pair<int, bool>> bh;
		nextNum = 0;
		for (int i = 0; i < num; i++) {
			int v = remain[i];
			if (!rightDel[v]) {
				if (degR[v] <= beta) {
					bfsQ.push_back({ v,false });
					rightDel[v] = true;
					bh.push_back({ v,false });
					for (int i = 0; i < bfsQ.size(); i++) {
						int p = bfsQ[i].first;
						if (bfsQ[i].second) {
							for (int j = 0; j < uNeighbor[p].size(); j++)
								if (!rightDel[uNeighbor[p][j]]) {
									degR[uNeighbor[p][j]]--;
									uDAG[p].push_back(uNeighbor[p][j]);
									if (degR[uNeighbor[p][j]] == beta) {
										bfsQ.push_back({ uNeighbor[p][j],false });
										rightDel[uNeighbor[p][j]] = true;
										bh.push_back({ uNeighbor[p][j],false });
									}
								}
						}
						else {
							for (int j = 0; j < vNeighbor[p].size(); j++)
								if (!leftDel[vNeighbor[p][j]]) {
									degL[vNeighbor[p][j]]--;
									vDAG[p].push_back(vNeighbor[p][j]);
									if (degL[vNeighbor[p][j]] < alpha) {
										bfsQ.push_back({ vNeighbor[p][j],true });
										leftDel[vNeighbor[p][j]] = true;
										bh.push_back({ vNeighbor[p][j],true });
									}
								}
						}
					}
					bfsQ.clear();
				}
				else {
					nextRemain[nextNum] = v;
					nextNum++;
				}
			}
		}
		order.push_back(bh);
		bh.clear();
		swap(remain, nextRemain);
		num = nextNum;
		if (nextNum == 0)
			break;
	}
	return 0;
}
int BiGraph::betaDecompose(int beta, vector<int>& leftDeg, vector<int>& rightDeg, vector<bool>& uDelete,
	vector<bool>& vDelete, vector<vector<int>>& uDAG, vector<vector<int>>& vDAG,
	vector<vector<pair<int, bool>>>& order) {
	//
	vector<int>leftQ, rightQ;
	//
	for (int i = 1; i < uDelete.size(); i++)
		if (!uDelete[i] && leftDeg[i] <= maxK) {
			leftQ.push_back(i);
			uDelete[i] = true;
		}
	for (int i = 1; i < vDelete.size(); i++)
		if (!vDelete[i] && rightDeg[i] < beta) {
			rightQ.push_back(i);
			vDelete[i] = true;
		}
	while (!leftQ.empty() || !rightQ.empty()) {
		for (int i = 0; i < leftQ.size(); i++) {
			int u = leftQ[i];
			for (int j = 0; j < uNeighbor[u].size(); j++) {
				int v = uNeighbor[u][j];
				if (vDelete[v])
					continue;
				rightDeg[v]--;
				if (rightDeg[v] < beta) {
					rightQ.push_back(v);
					vDelete[v] = true;
				}
			}
		}
		leftQ.clear();
		for (int i = 0; i < rightQ.size(); i++) {
			int v = rightQ[i];
			for (int j = 0; j < vNeighbor[v].size(); j++) {
				int u = vNeighbor[v][j];
				if (uDelete[u])
					continue;
				leftDeg[u]--;
				if (leftDeg[u] <= maxK) {
					leftQ.push_back(u);
					uDelete[u] = true;
				}
			}
		}
		rightQ.clear();
	}
	int num = n1;
	int nextNum = 0;
	vector<pair<int, bool>>bfsQ;
	vector<int>remain(num);
	vector<int>nextRemain(num);
	vector<bool>leftDel = uDelete;
	vector<bool>rightDel = vDelete;
	vector<int>degL = leftDeg;
	vector<int>degR = rightDeg;

	uDAG.resize(n1 + 1);
	vDAG.resize(n2 + 1);

	for (int i = 0; i < remain.size(); i++)
		remain[i] = i + 1;
	for (int alpha = maxK + 1; alpha <= maxUDeg + 1; alpha++) {
		vector<pair<int, bool>> bh;
		nextNum = 0;
		for (int i = 0; i < num; i++) {
			int u = remain[i];
			if (!leftDel[u]) {
				if (degL[u] <= alpha) {
					bfsQ.push_back({ u,true });
					leftDel[u] = true;
					bh.push_back({ u,true });
					for (int i = 0; i < bfsQ.size(); i++) {
						int p = bfsQ[i].first;
						if (bfsQ[i].second) {
							for (int j = 0; j < uNeighbor[p].size(); j++)
								if (!rightDel[uNeighbor[p][j]]) {
									degR[uNeighbor[p][j]]--;
									uDAG[p].push_back(uNeighbor[p][j]);
									if (degR[uNeighbor[p][j]] < beta) {
										bfsQ.push_back({ uNeighbor[p][j],false });
										rightDel[uNeighbor[p][j]] = true;
										bh.push_back({ uNeighbor[p][j],false });
									}
								}
						}
						else {
							for (int j = 0; j < vNeighbor[p].size(); j++)
								if (!leftDel[vNeighbor[p][j]]) {
									degL[vNeighbor[p][j]]--;
									vDAG[p].push_back(vNeighbor[p][j]);
									if (degL[vNeighbor[p][j]] == alpha) {
										bfsQ.push_back({ vNeighbor[p][j],true });
										leftDel[vNeighbor[p][j]] = true;
										bh.push_back({ vNeighbor[p][j],true });
									}
								}
						}
					}
					bfsQ.clear();
				}
				else {
					nextRemain[nextNum] = u;
					nextNum++;
				}
			}
		}
		order.push_back(bh);
		bh.clear();
		swap(remain, nextRemain);
		num = nextNum;
		if (nextNum == 0)
			break;
	}
	return 0;
}

int BiGraph::comBN(pair<int, bool>ab, vector<vector<int>>& uDAG, vector<vector<int>>& vDAG,
	vector<vector<pair<int, bool>>>& order) {

	UnionFind uf;
	umap.clear();
	vmap.clear();

	for (int i = order.size() - 1; i >= 0; i--) {
		tempSet.clear();
		for (int j = order[i].size() - 1; j >= 0; j--) {
			auto p = order[i][j];
			if (p.second) {
				set<int>roots;
				set<pair<bool, int>>mark;
				for (int k = 0; k < uDAG[p.first].size(); k++)
					roots.insert(uf.Find(vmap[uDAG[p.first][k]]));
				umap[p.first] = uf.BatchUnite(roots, p.second, p.first, mark);
				tempSet.push_back({ p.first,p.second,mark });
			}
			else {
				set<int>roots;
				set<pair<bool, int>>mark;
				for (int k = 0; k < vDAG[p.first].size(); k++)
					roots.insert(uf.Find(umap[vDAG[p.first][k]]));
				vmap[p.first] = uf.BatchUnite(roots, p.second, p.first, mark);
				tempSet.push_back({ p.first,p.second,mark });
			}
		}
		gamSet.clear();
		for (auto j : tempSet) {
			auto p = get<0>(j);
			auto is = get<1>(j);
			auto con = get<2>(j);
			int alpha, beta, gamma;
			if (ab.second) {
				alpha = ab.first;
				beta = i + 1;
			}
			else {
				alpha = i + maxK + 1;
				beta = ab.first;
			}
			if (is)
				gamma = uf.Find(umap[p]);
			else
				gamma = uf.Find(vmap[p]);
			if (gamSet.find(gamma) == gamSet.end()) {
				gamSet.insert(gamma);
				cMap[{ alpha, beta, gamma }] = cluster.size();
				cluster.push_back({ alpha,beta,gamma });
			}
			updateNumber(ab.first, ab.second, p, is, cMap[{ alpha, beta, gamma }]);
			for (auto k : con) {
				int nbr = 0;
				if (k.first) {
					if (uBeta[k.second].size() > 0)
						nbr = uBeta[k.second][uBeta[k.second].size() - 1];
					else
						nbr = uNumber[k.second][uNumber[k.second].size() - 1];
				}
				else {
					if (vBeta[k.second].size() > 0)
						nbr = vBeta[k.second][vBeta[k.second].size() - 1];
					else
						nbr = vNumber[k.second][vNumber[k.second].size() - 1];
				}
				if ((get<0>(cluster[nbr]) > alpha && get<1>(cluster[nbr]) == beta) ||
					(get<0>(cluster[nbr]) == alpha && get<1>(cluster[nbr]) > beta))
					block1.insert({ cMap[{ alpha, beta, gamma }],nbr });
			}
		}
	}
	return 0;
}

int BiGraph::updateNumber(int d, bool isAlpha, int p, bool isLeft, int x) {
	if (isAlpha && isLeft) {
		if (uNumber[p].size() > 0) {
			int nf = uNumber[p][uNumber[p].size() - 1];
			if (get<1>(cluster[nf]) == get<1>(cluster[x])) {
				uNumber[p][uNumber[p].size() - 1] = x;
				block2.insert({ x,nf });
			}
			else {
				uNumber[p].push_back(x);
			}
		}
		else
			uNumber[p].push_back(x);
	}
	else if (isAlpha && !isLeft) {
		if (vNumber[p].size() > 0) {
			int nf = vNumber[p][vNumber[p].size() - 1];
			if (get<1>(cluster[nf]) == get<1>(cluster[x])) {
				vNumber[p][vNumber[p].size() - 1] = x;
				block2.insert({ x,nf });
			}
			else {
				vNumber[p].push_back(x);
			}
		}
		else
			vNumber[p].push_back(x);
	}
	else if (!isAlpha && isLeft) {
		if (uBeta[p].size() > 0) {
			int nf = uBeta[p][uBeta[p].size() - 1];
			if (get<1>(cluster[nf]) == get<1>(cluster[x])) {
				uBeta[p][uBeta[p].size() - 1] = x;
				block2.insert({ x,nf });
			}
			else {
				uBeta[p].push_back(x);
			}
		}
		else
			uBeta[p].push_back(x);
	}
	else {
		if (vBeta[p].size() > 0) {
			int nf = vBeta[p][vBeta[p].size() - 1];
			if (get<1>(cluster[nf]) == get<1>(cluster[x])) {
				vBeta[p][vBeta[p].size() - 1] = x;
				block2.insert({ x,nf });
			}
			else {
				vBeta[p].push_back(x);
			}
		}
		else
			vBeta[p].push_back(x);
	}
	return 0;
}

int BiGraph::verifyCom(BiGraph& g, vector<bool>& leftResult, vector<bool>& rightResult, int alpha, int beta, int q, bool isLeft) {
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
		cout << "Verify: false";
		if (lower > 0) {
			cout << "  lower: " << lower;
			cout << "  Input:" << alpha << "  " << beta << "  " << q;
		}
		if (more > 0) {
			cout << "  more: " << more;
			cout << "  Input:" << alpha << "  " << beta << "  " << q;
		}
		cout << endl;
		return 0;
	}
}

int BiGraph::verifyCore(BiGraph& g, vector<int>& leftResult, vector<int>& rightResult, int alpha, int beta) {
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
	int low = 0;
	for (int i = 1; i <= g.n1; i++)
		if ((leftVerify[i] && leftResult[i] == 0) || (!leftVerify[i] && leftResult[i] != 0)) {
			low++;
			//cout << "U" << i << endl;
		}
	for (int i = 1; i <= g.n2; i++)
		if ((rightVerify[i] && rightResult[i]==0)||(!rightVerify[i] && rightResult[i] != 0)){
			low++;
			//cout << "V" << i << endl;
		}
	if (low > 0) {
		cout << "verify: false  low: "<<low << endl;
		return 0;
	}
	cout << "verify: true" << endl;
	return 0;

}
