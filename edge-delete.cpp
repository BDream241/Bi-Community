#include"maintenance.h"

void Maintenance::edgeDelete(int u, int v) {
	int falg = 0;
	for (int i = 0; i < uNeighbor[u].size(); i++)
		if (uNeighbor[u][i] == v)
			falg = 1;
	if (falg == 0)
		return;
	g.deleteEdge(u, v);
	uNeighbor[u].erase(std::remove(uNeighbor[u].begin(), uNeighbor[u].end(), v), uNeighbor[u].end());
	vNeighbor[v].erase(std::remove(vNeighbor[v].begin(), vNeighbor[v].end(), u), vNeighbor[v].end());

	//分析需要考虑的双核数
	vector<int>ubn = uNumber[u], vbn = vNumber[v];
	for (int i = 0; i < ubn.size(); i++) {
		int alpha = asg.asg[ubn[i]].alpha;
		int beta = asg.asg[ubn[i]].beta;
		int Bv = findOff(v, false, alpha, true);
		if ( Bv < beta)
			continue;
		vector<pair<bool, int>> R = neighborBaseDel(u, true, alpha, beta);
		alphaDecrease(u, v, alpha, beta, R);
		betaDecrease(u, v, alpha, beta, R);
	}
	for (int i = 0; i < vbn.size(); i++) {
		int alpha = asg.asg[vbn[i]].alpha;
		int beta = asg.asg[vbn[i]].beta;
		int Bu = findOff(u, true, beta, false);
		if (Bu < alpha)
			continue;
		vector<pair<bool, int>> R = neighborBaseDel(v, false, alpha, beta);
		alphaDecrease2(u, v, alpha, beta, R);
		betaDecrease2(u, v, alpha, beta, R);
	}
}
void Maintenance::alphaDecrease(int u, int v, int alpha, int beta, vector<pair<bool, int>>& R) {
	//top-beta
	vector<int>topX;
	for (int i = 0; i < uNeighbor[u].size(); i++)
		topX.push_back(findOff(uNeighbor[u][i], false, alpha -1, true));
	sort(topX.begin(), topX.end(), [](int a, int b) {return a > b; });
	int Lu;
	if (alpha-1 < topX.size())
		Lu = topX[alpha -1];
	else
		return ;
	//update
	if (R.size() == 0)
		return;
	dealC2(alpha - 1, Lu, R);
}
void Maintenance::betaDecrease(int u, int v, int alpha, int beta, vector<pair<bool, int>>& R) {
	//top-beta
	vector<int>topX;
	for (int i = 0; i < uNeighbor[u].size(); i++)
		topX.push_back(findOff(uNeighbor[u][i], false, alpha, true));
	sort(topX.begin(), topX.end(), [](int a, int b) {return a > b; });
	int Lu ;
	if (alpha < topX.size())
		Lu = topX[alpha ];
	else
		return;

	//update
	if (R.size() == 0)
		return;
	dealC2(alpha, beta - 1, R);
}
void Maintenance::alphaDecrease2(int u, int v, int alpha, int beta, vector<pair<bool, int>>& R) {
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
	dealC2(alpha-1, beta, R);
}
void Maintenance::betaDecrease2(int u, int v, int alpha, int beta, vector<pair<bool, int>>& R) {
	//top-beta
	vector<int>topX;
	for (int i = 0; i < vNeighbor[v].size(); i++)
		topX.push_back(findOff(vNeighbor[v][i], true, beta-1, false));
	sort(topX.begin(), topX.end(), [](int a, int b) {return a > b; });
	int Lu;
	if (beta-1 < topX.size())
		Lu = topX[beta-1];
	else
		return;

	//update
	if (R.size() == 0)
		return;
	dealC2(Lu, beta - 1, R);
}
vector<pair<bool, int>> Maintenance::neighborBaseDel(int u, bool isLeft, int a, int b) {
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
				for (int j = 0; j < vNumber[i].size(); j ++) 
					if (asg.asg[vNumber[i][j]].alpha >= a  && asg.asg[vNumber[i][j]].beta >= b) {
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
					if (asg.asg[uNumber[i][j]].alpha >= a && asg.asg[uNumber[i][j]].beta >= b) {
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
void Maintenance::addCand(vector<pair<bool, int>>& R, int w, bool isLeft, int a, int b, queue<pair<bool, int>>& S, vector<bool>& uS, vector<bool>& vS) {
	if (isLeft) {
		uC[w] = true;
		R.push_back({ true,w });
		for (auto i : uNeighbor[w]) {
			if (!vT[i] && findOff(i, false, a, true) == b && !vS[i]) {
				S.push({ false,i });
			}
			if (vT[i] && !vC[i]) {
				vDeg[i]--;
				if(vDeg[i]<b)
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
void Maintenance::dealC2(int a, int b, vector<pair<bool, int>>R) {
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
	int node = -1;
	if (C.first.size() > 0) {
		for (auto w : uNeighbor[C.first[0]])
			for (int i = 0; i < asg.vIndex[w].size(); i++)
				if (asg.asg[asg.vIndex[w][i]].alpha == a && asg.asg[asg.vIndex[w][i]].beta == b) {
					node = asg.vIndex[w][i];
					break;
				}
	}
	else if (C.second.size() > 0) {
		int node = -1;
		for (auto w : vNeighbor[C.second[0]])
			for (int i = 0; i < asg.uIndex[w].size(); i++)
				if (asg.asg[asg.uIndex[w][i]].alpha == a && asg.asg[asg.uIndex[w][i]].beta == b) {
					node = asg.uIndex[w][i];
					break;
				}
	}
	if (node == -1) {
		node = asg.asg.size();
		Node* newNode = new Node;
		newNode->alpha = a;
		newNode->beta = b;
		asg.asg.push_back(*newNode);
	}
	nodeQ.insert(node);
	asg.asg[node].neighbor.clear();
	for (int i = 0; i < C.first.size(); i++) {
		update(C.first[i], true, node);
		asg.uIndex[C.first[i]].push_back(node);
		asg.asg[node].left.push_back(C.first[i]);
		for (auto j : asg.uIndex[C.first[i]])
			if (asg.asg[j].alpha >= a && asg.asg[j].beta >= b) {
				asg.asg[j].left.erase(std::remove(asg.asg[j].left.begin(), asg.asg[j].left.end(), C.first[i]), asg.asg[j].left.end());
				nodeQ.insert(j);
				asg.asg[j].neighbor.clear();
			}
	}
	for (int i = 0; i < C.second.size(); i++) {
		update(C.second[i], false, node);
		asg.vIndex[C.second[i]].push_back(node);
		asg.asg[node].right.push_back(C.second[i]);
		for (auto j : asg.vIndex[C.second[i]])
			if (asg.asg[j].alpha >= a && asg.asg[j].beta >= b) {
				asg.asg[j].right.erase(std::remove(asg.asg[j].right.begin(), asg.asg[j].right.end(), C.second[i]), asg.asg[j].right.end());
				nodeQ.insert(j);
				asg.asg[j].neighbor.clear();
			}
	}
}