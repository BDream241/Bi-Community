#include"bigraph.h"

BiGraph::BiGraph(){}
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

void BiGraph::printBiCore() {
	for (int i = 1; i <= n1; i++) {
		cout << i << " : ";
		for (auto j : uNumber[i])
			cout << " (" << get<0>(j) << "," << get<1>(j) << "," << get<2>(j) << ") ";
		cout << endl;
	}
	for (int i = 1; i <= n2; i++) {
		cout << i << " : ";
		for (auto j : vNumber[i])
			cout << " (" << get<0>(j) << "," << get<1>(j) << "," << get<2>(j) << ") ";
		cout << endl;
	}
}

