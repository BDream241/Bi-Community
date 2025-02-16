#include"maintenance.h"

Maintenance::Maintenance() {}
Maintenance::Maintenance(BiGraph& g, ASG& a) {
	this->asg = a;
	this->g = g;
	n1 = g.n1;
	n2 = g.n2;
	uNeighbor = g.uNeighbor;
	vNeighbor = g.vNeighbor;
	uNumber = asg.uIndex;
	vNumber = asg.vIndex;
}
void Maintenance::init(BiGraph& g, ASG& a) {
	this->asg = a;
	this->g = g;
	n1 = g.n1;
	n2 = g.n2;
	uNeighbor = g.uNeighbor;
	vNeighbor = g.vNeighbor;
	uNumber = asg.uIndex;
	vNumber = asg.vIndex;
}
Maintenance::Maintenance(string adds, ASG& myAsg, BiGraph& graph) {
	vector<pair<int, int>>edges;
	FILE* fp = fopen(adds.c_str(), "r");
	int u, v;
	while (fscanf(fp, "%ld %ld", &u, &v) != EOF)
		edges.push_back({ u,v });
	fclose(fp);
	auto start = chrono::system_clock::now();
	init(graph, myAsg);
	for (auto i : edges) {
		//cout << i.first << " " << i.second << endl;
		edgeInsert(i.first, i.second);
	}
	asg.creatSEviaNode(g, nodeQ);
	asg.addEdge();
	auto end = chrono::system_clock::now();
	chrono::duration<double> time = end - start;
	cout << "insert time: " << time.count() * 1000 << "ms" << endl;
}
void Maintenance::edgeInsert(int u, int v) {
	for (int i = 0; i < uNeighbor[u].size(); i++)
		if (uNeighbor[u][i] == v)
			return;
	g.addEdge(u, v);
	uNeighbor[u].push_back(v);
	vNeighbor[v].push_back(u);

	vector<pair<int, int>>H;
	//分析需要考虑的双核数
	vector<int> ubn = uNumber[u], vbn = vNumber[v];

	for (int i = 0; i < ubn.size(); i++) {
		int alpha = asg.asg[ubn[i]].alpha;
		int beta = asg.asg[ubn[i]].beta;
		int Bu = findOff(u, true, alpha + 1, true);
		int Bv = findOff(v, false, alpha + 1, true);
		if (Bv != 0 && Bv >= Bu)
			alphaIncrease(u, v, alpha, beta, H);
		Bv = findOff(v, false, alpha, true);
		if (Bv != 0 && Bv > beta)
			betaIncrease(u, v, alpha, beta, H);
	}
	for (int i = 0; i < vbn.size(); i++) {
		int alpha = asg.asg[vbn[i]].alpha;
		int beta = asg.asg[vbn[i]].beta;
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
void Maintenance::alphaIncrease(int u, int v, int alpha, int beta, vector<pair<int, int>>& H) {
	//top-beta
	vector<int>topX;
	for (int i = 0; i < uNeighbor[u].size(); i++)
		topX.push_back(findOff(uNeighbor[u][i], false, alpha + 1, true));
	sort(topX.begin(), topX.end(), [](int a, int b) {return a > b;});
	int Lu;
	if (alpha+1 < topX.size())
		Lu = topX[alpha + 1];
	else
		return ;
	//update
	bool had = false;
	for(int i=0;i<H.size();i++)
		if (H[i].first == alpha + 1 && H[i].second == Lu+1) {
			had = true;
			break;
		}
	if (!had) {
		//bool flag = update(u, true, alpha + 1, Lu+1);
		neighborBaseIns(u, true,alpha + 1, Lu + 1);
		H.push_back(make_pair(alpha + 1, Lu+1));
	}
}
void Maintenance::betaIncrease(int u, int v, int alpha, int beta, vector<pair<int, int>>& H) {
	//top-beta
	vector<int>topX;
	for (int i = 0; i < uNeighbor[u].size(); i++)
		topX.push_back(findOff(uNeighbor[u][i], false, alpha, true));
	sort(topX.begin(), topX.end(), [](int a, int b) {return a > b; });
	int Lu;
	if (alpha < topX.size())
		Lu = topX[alpha];
	else
		return ;
	//update
	bool had = false;
	for (int i = 0; i < H.size(); i++)
		if (H[i].first == alpha && H[i].second == Lu+1) {
			had = true;
			break;
		}
	if (!had) {
		//bool flag = update(u, true, alpha, Lu+1);
		neighborBaseIns(u, true, alpha, Lu+1);
		H.push_back(make_pair(alpha, Lu+1));
	}
}
void Maintenance::alphaIncrease2(int u, int v, int alpha, int beta, vector<pair<int, int>>& H) {
	//top-beta
	vector<int>topX;
	for (int i = 0; i < vNeighbor[v].size(); i++)
		topX.push_back(findOff(vNeighbor[v][i], true, beta, false));
	sort(topX.begin(), topX.end(), [](int a, int b) {return a > b; });
	int Lu;
	if (beta < topX.size())
		Lu = topX[beta];
	else
		return ;
	//update
	bool had = false;
	for (int i = 0; i < H.size(); i++)
		if (H[i].first == Lu+1 && H[i].second == beta) {
			had = true;
			break;
		}
	if (!had) {
		//bool flag = update(v, false, Lu+1, beta);
		neighborBaseIns2(v, false, Lu+1, beta);
		H.push_back(make_pair(Lu+1, beta));
	}
}
void Maintenance::betaIncrease2(int u, int v, int alpha, int beta, vector<pair<int, int>>& H) {
	//top-beta
	vector<int>topX;
	for (int i = 0; i < vNeighbor[v].size(); i++)
		topX.push_back(findOff(vNeighbor[v][i], true, beta+1, false));
	sort(topX.begin(), topX.end(), [](int a, int b) {return a > b; });
	int Lu;
	if (beta + 1 < topX.size())
		Lu = topX[beta + 1];
	else
		return;
	//update
	bool had = false;
	for (int i = 0; i < H.size(); i++)
		if (H[i].first == Lu+1 && H[i].second == beta+1) {
			had = true;
			break;
		}
	if (!had) {
		//bool flag = update(v, false, Lu+1, beta+1);
		neighborBaseIns2(v, false, Lu+1, beta+1);
		H.push_back(make_pair(Lu+1, beta+1));
	}
}
void Maintenance::neighborBaseIns(int u, bool isLeft, int a, int b) {
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
		R.push_back({is,w});
		S.pop();
		queue<pair<bool, int>>temp = S;
		if (is) {
			uC[w] = true;
			uT[w] = true;
			for (auto i : uNeighbor[w]) {
				for (auto j : vNumber[i]) {
					if (!vC[i] && asg.asg[j].alpha>=a && asg.asg[j].beta>=b-1) {
						uDeg[w]++;
						break;
					}
					if (!vT[i] && asg.asg[j].alpha == a && asg.asg[j].beta == b - 1) {
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
					if (!uC[i] && asg.asg[j].alpha >= a && asg.asg[j].beta >= b - 1) {
						vDeg[w]++;
						break;
					}
					if (!uT[i] && asg.asg[j].alpha == a && asg.asg[j].beta == b - 1) {
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
void Maintenance::neighborBaseIns2(int u, bool isLeft, int a, int b) {
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
					if (!vC[i] && asg.asg[j].alpha >= a-1 && asg.asg[j].beta >= b) {
						uDeg[w]++;
						break;
					}
					if (!vT[i] && asg.asg[j].alpha == a-1 && asg.asg[j].beta == b) {
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
					if (!uC[i] && asg.asg[j].alpha >= a && asg.asg[j].beta >= b - 1) {
						vDeg[w]++;
						break;
					}
					if (!uT[i] && asg.asg[j].alpha == a && asg.asg[j].beta == b - 1) {
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
void Maintenance::deleteCand(int w, bool isLeft, int a, int b) {
	if (isLeft) {
		uC[w] = false;
		for (auto i : uNeighbor[w]) 
			if(vC[i]){
				vDeg[i]--;
				if (vDeg[i] < b)
					deleteCand(i, false, a, b);
			}
	}
	else {
		vC[w] = false;
		for (auto i : vNeighbor[w]) 
			if(uC[i]) {
				uDeg[i]--;
				if (uDeg[i] < b)
					deleteCand(i, true, a, b);
			}
	}
}
bool Maintenance::update(int u, bool isLeft, int node) {
	if (isLeft) {
		for (int i = 0; i < uNumber[u].size(); i++)
			if (asg.asg[uNumber[u][i]].alpha >= asg.asg[node].alpha 
				&& asg.asg[uNumber[u][i]].beta >= asg.asg[node].beta)
				return false;
		uNumber[u].push_back(node);
		return true;
	}
	else {
		for (int i = 0; i < vNumber[u].size(); i++)
			if (asg.asg[vNumber[u][i]].alpha >= asg.asg[node].alpha
				&& asg.asg[vNumber[u][i]].beta >= asg.asg[node].beta)
				return false;
		vNumber[u].push_back(node);
		return true;
	}
	
}
int Maintenance::findOff(int u, bool isLeft, int value, bool isAlpha) {
	int off = 0;
	if (isLeft) {
		if (isAlpha) {
			for(int i=0;i<uNumber[u].size(); i += 2)
				if (asg.asg[uNumber[u][i]].alpha >= value && asg.asg[uNumber[u][i]].beta>off) 
					off = asg.asg[uNumber[u][i]].beta;
		}
		else {
			for (int i = 0; i < uNumber[u].size(); i += 2)
				if (asg.asg[uNumber[u][i]].alpha >= off && asg.asg[uNumber[u][i]].beta > value)
					off = asg.asg[uNumber[u][i]].alpha;
		}
	}
	else {
		if (isAlpha) {
			for (int i = 0; i < vNumber[u].size(); i += 2)
				if (asg.asg[vNumber[u][i]].alpha >= value && asg.asg[vNumber[u][i]].beta > off)
					off = asg.asg[vNumber[u][i]].beta;
		}
		else {
			for (int i = 0; i < vNumber[u].size(); i += 2)
				if (asg.asg[vNumber[u][i]].alpha >= off && asg.asg[vNumber[u][i]].beta > value)
					off = asg.asg[vNumber[u][i]].alpha;
		}
	}
	return off;
}
void Maintenance::dealC(int a,int b, vector<pair<bool, int>>R) {
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
			if (asg.asg[j].alpha <= a && asg.asg[j].beta <= b) {
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
			if (asg.asg[j].alpha <= a && asg.asg[j].beta <= b) {
				asg.asg[j].right.erase(std::remove(asg.asg[j].right.begin(), asg.asg[j].right.end(), C.second[i]), asg.asg[j].right.end());
				nodeQ.insert(j);
				asg.asg[j].neighbor.clear();
			}
	}
}

