#pragma once
#include"utility.h"
#include"bigraph.h"
#include"dgraph.h"
#include"abIndex.h"

int main(int argc, char* argv[]) {
	string exec_type = argv[1];
	string dataset = argv[2];

	if (exec_type == "decompose") {
		BiGraph graph(dataset);
		graph.decompose();
	}
	else if (exec_type == "community-search") {
		int q = atoi(argv[3]);
		int alpha = atoi(argv[4]);
		int beta = atoi(argv[5]);

		BiGraph graph(dataset);
		graph.decompose();
		abIndex ab(graph);
		ab.creatSGOpt(graph);
		vector<bool>leftResult(graph.n1 + 1, false);
		vector<bool>rightResult(graph.n2 + 1, false);
		leftResult.resize(graph.n1 + 1, false);
		rightResult.resize(graph.n2 + 1, false);
		ab.queryOpt(q, alpha, beta, leftResult, rightResult, graph.maxK);
		cout << q << " " << alpha << " " << beta << ": ";
		graph.verifyCom(graph, leftResult, rightResult, q, alpha, beta, 1);
		leftResult.clear();
		rightResult.clear();
	}
	else if (exec_type == "bi-core-search") {
		int alpha = atoi(argv[3]);
		int beta = atoi(argv[4]);

		BiGraph graph(dataset);
		graph.decompose();
		abIndex ab(graph);
		ab.creatSGOpt(graph);
		ab.creatTree();
		vector<int>leftResult(graph.n1 + 1, false);
		vector<int>rightResult(graph.n2 + 1, false);
		leftResult.resize(graph.n1 + 1, 0);
		rightResult.resize(graph.n2 + 1, 0);
		leftResult.resize(graph.n1 + 1, 0);
		rightResult.resize(graph.n2 + 1, 0);
		ab.queryCore(alpha, beta, leftResult, rightResult, graph.maxK);
		cout  << alpha << " " << beta << ": ";
		graph.verifyCore(graph, leftResult, rightResult, alpha, beta);
		leftResult.clear();
		rightResult.clear();
	}
}

