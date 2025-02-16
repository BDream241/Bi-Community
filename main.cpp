#include"maintenance.h"
#include <random>

int main() {
	int arg = 2;
	string str = "../data/WC";
	//core分解
	if (arg == 0) {
		BiGraph graph(str);
		cout << "max k-core: " << graph.coreDecompose() << endl;
	}
	//alpha、beta分解
	else if (arg == 1) {
		BiGraph graph(str);
		graph.decompose();
	}
	//query
	else if (arg == 2) {
		BiGraph graph(str);
		auto start = chrono::system_clock::now();
		graph.decompose();
		auto end = chrono::system_clock::now();
		chrono::duration<double> time = end - start;
		cout << "decompose time: " << time.count() << endl;
		start = chrono::system_clock::now();
		ASG  myAsg(graph);
		MHG myMhG(myAsg.asg, graph.maxK, graph.maxAlpha, graph.maxBeta, graph.n1, graph.n2);
		end = chrono::system_clock::now();
		time = end - start;
		cout << "bulid time: " << time.count() << endl;
		string strq = str + "/query.txt";
		FILE* fq = fopen(strq.c_str(), "r");
		bool isLeft;
		int q, alpha, beta;
		vector<int>aList, bList, nList, isList;
		while (fscanf(fq, "%d%d%d%d", &q, &isLeft, &alpha, &beta) != EOF) {
			aList.push_back(alpha);
			bList.push_back(beta);
			nList.push_back(q);
			isList.push_back(isLeft);
		}
		fclose(fq);
		vector<bool>leftResult(graph.n1 + 1, false);
		vector<bool>rightResult(graph.n2 + 1, false);
		for (int i = 0; i < 5; i++) {
			start = chrono::system_clock::now();
			for (int j = 0; j < 100; j++) {
				leftResult.resize(graph.n1 + 1, false);
				rightResult.resize(graph.n2 + 1, false);
				myMhG.QueryViaMHG(aList[i * 100 + j], bList[i * 100 + j], nList[i * 100 + j], isList[i * 100 + j], leftResult, rightResult);
				//verifyCom(graph, leftResult, rightResult, aList[i * 100 + j], bList[i * 100 + j], nList[i * 100 + j], isList[i * 100 + j]);
				leftResult.clear();
				rightResult.clear();
			}
			end = chrono::system_clock::now();
			time = end - start;
			cout << i + 1 << "th time: " << time.count() * 1000 << endl;
		}
	}
	else if (arg == 3) {
		BiGraph graph(str);
		graph.decompose();
		ASG  myAsg(graph);
		string adds = str + "/insert.txt";
		Maintenance mt(adds, myAsg, graph);
	}
	return 0;
}