#include"bigraph.h"

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

		//start = chrono::system_clock::now();
		dealDAG({ alpha,true }, uDAG, vDAG, order);
		//end = chrono::system_clock::now();
		//time += end - start;

		if (alpha % 20 == 0)
			cout << alpha << endl;
		else if (alpha == 1)
			maxBeta = order.size();

		uDAG.clear();
		vDAG.clear();
		order.clear();
	}
	for (int beta = 1; beta <= maxK; beta++) {
		betaDecompose(beta, uDeg, vDeg, uDel, vDel, uDAG, vDAG, order);

		//start = chrono::system_clock::now();
		dealDAG({ beta,false }, uDAG, vDAG, order);
		//end = chrono::system_clock::now();
		//time += end - start;


		if (beta % 20 == 0)
			cout << beta << endl;
		else if (beta == 1)
			maxAlpha = order.size() + maxK;

		uDAG.clear();
		vDAG.clear();
		order.clear();
	}

	//cout << "分解: " << time.count() << endl;

	for (int i = 1; i <= n1; i++) 
		if(uBeta[i].size()>0) {
			if (get<1>(uNumber[i][uNumber[i].size() - 1]) == get<1>(uBeta[i][uBeta[i].size() - 1])) {
				get<0>(uNumber[i][uNumber[i].size() - 1]) = get<0>(uBeta[i][uBeta[i].size() - 1]);
				get<2>(uNumber[i][uNumber[i].size() - 1]) = get<2>(uBeta[i][uBeta[i].size() - 1]);
			}
			else
				uNumber[i].push_back(uBeta[i][uBeta[i].size() - 1]);
			for (int j = uBeta[i].size() - 2; j >= 0; j--)
				uNumber[i].push_back(uBeta[i][j]);
		}
	for (int i = 1; i <= n2; i++) 
		if (vBeta[i].size() > 0) {
			if (get<1>(vNumber[i][vNumber[i].size() - 1]) == get<1>(vBeta[i][vBeta[i].size() - 1])) {
				get<0>(vNumber[i][vNumber[i].size() - 1]) = get<0>(vBeta[i][vBeta[i].size() - 1]);
				get<2>(vNumber[i][vNumber[i].size() - 1]) = get<2>(vBeta[i][vBeta[i].size() - 1]);
			}
			else
				vNumber[i].push_back(vBeta[i][vBeta[i].size() - 1]);
			for (int j = vBeta[i].size() - 2; j >= 0; j--)
				vNumber[i].push_back(vBeta[i][j]);
	}
	return 0;
}
int BiGraph::coreDecompose() {
	vector<int> leftQ;//左侧删除队列
	vector<int> rightQ;//右侧删除队列

	int num1 = n1 + 1;//左侧节点数量
	vector<int> leftR(num1);//用于存储当前迭代中尚未被删除的顶点
	for (int i = 0; i < leftR.size(); i++)
		leftR[i] = i;

	int leftRTnum = 0;//用于记录在当前迭代中尚未被删除的顶点数量。
	vector<int> leftRT(num1);//用于临时存储在当前迭代中尚未被删除的顶点

	int num2 = n2 + 1;//
	vector<int> rightR(num2);//用于存储当前迭代中尚未被删除的顶点
	for (int i = 0; i < rightR.size(); i++)
		rightR[i] = i;

	int rightRTnum = 0;//用于记录在当前迭代中尚未被删除的顶点数量。
	vector<int> rightRT(num2);//用于临时存储在当前迭代中尚未被删除的顶点

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
int BiGraph::alphaDecompose(int alpha, vector<int>& leftDeg, vector<int>& rightDeg,vector<bool>& uDelete, 
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
					bfsQ.push_back({v,false});
					rightDel[v] = true;
					bh.push_back({ v,false });
					for (int i = 0; i < bfsQ.size(); i++) {
						int p = bfsQ[i].first;
						if (bfsQ[i].second) {
							for(int j=0;j<uNeighbor[p].size();j++)
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
	vector<bool>& vDelete, vector<vector<int>>&uDAG,vector<vector<int>>&vDAG, 
	vector<vector<pair<int, bool>>>&order) {
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
	//初始化
	int num = n1;
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
int BiGraph::dealDAG(pair<int,bool>ab, vector<vector<int>>& uDAG, vector<vector<int>>& vDAG, 
	vector<vector<pair<int, bool>>>& order) {

	UnionFind uf;
	unordered_map<int,int>umap, vmap;

	for (int i = order.size() - 1; i >=0; i--) {
		for (int j = order[i].size() - 1; j >= 0; j--) {
			auto p = order[i][j];
			if (p.second) {
				set<int>roots;
				for (int k = 0; k < uDAG[p.first].size(); k++)
					roots.insert(uf.Find(vmap[uDAG[p.first][k]]));
				umap[p.first] =  uf.BatchUnite(roots);
			}
			else {
				set<int>roots;
				for (int k = 0; k < vDAG[p.first].size(); k++)
					roots.insert(uf.Find(umap[vDAG[p.first][k]]));
				vmap[p.first] =  uf.BatchUnite(roots);
			}
		}
		if (ab.second) {
			for (int j = order[i].size() - 1; j >= 0; j--) {
				if (order[i][j].second)
					updateNumber(ab.first, i + 1, uf.Find(umap[order[i][j].first]), order[i][j], true);
				else
					updateNumber(ab.first, i + 1, uf.Find(vmap[order[i][j].first]), order[i][j], true);
			}
		}
		else {
			for (int j = order[i].size() - 1; j >= 0; j--) {
				if(order[i][j].second)
					updateNumber(i + 1 + maxK, ab.first, uf.Find(umap[order[i][j].first]), order[i][j], false);
				else
					updateNumber(i + 1 + maxK, ab.first, uf.Find(vmap[order[i][j].first]), order[i][j], false);
			}
		}
	}
	return 0;
}

int BiGraph::updateNumber(int alpha, int beta, int x, pair<int, bool>p, bool isAlpha) {
	// uNumber, vNumber;
	if (isAlpha) {
		if (p.second) {
			if (uNumber[p.first].size() > 0 && get<1>(uNumber[p.first][uNumber[p.first].size() - 1]) == beta) {
				get<0>(uNumber[p.first][uNumber[p.first].size() - 1]) = alpha;
				get<2>(uNumber[p.first][uNumber[p.first].size() - 1]) = x;
			}
			else
				uNumber[p.first].emplace_back(alpha, beta, x);
		}
		else {
			if (vNumber[p.first].size() > 0 && get<1>(vNumber[p.first][vNumber[p.first].size() - 1]) == beta) {
				get<0>(vNumber[p.first][vNumber[p.first].size() - 1]) = alpha;
				get<2>(vNumber[p.first][vNumber[p.first].size() - 1]) = x;
			}
			else
				vNumber[p.first].emplace_back(alpha, beta, x);
		}
	}
	else {
		if (p.second) {
			if (uBeta[p.first].size() > 0 && get<0>(uBeta[p.first][uBeta[p.first].size() - 1]) == alpha) {
				get<1>(uBeta[p.first][uBeta[p.first].size() - 1]) = beta;
				get<2>(uBeta[p.first][uBeta[p.first].size() - 1]) = x;
			}
			else
				uBeta[p.first].emplace_back(alpha, beta, x);
		}
		else {
			if (vBeta[p.first].size() > 0 && get<0>(vBeta[p.first][vBeta[p.first].size() - 1]) == alpha) {
				get<1>(vBeta[p.first][vBeta[p.first].size() - 1]) = beta;
				get<2>(vBeta[p.first][vBeta[p.first].size() - 1]) = x;
			}
			else
				vBeta[p.first].emplace_back(alpha, beta, x);
		}
	}
	return 0;
}