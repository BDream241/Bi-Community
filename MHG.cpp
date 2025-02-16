#include"MHG.h"

MHG::MHG(){}
MHG::MHG(vector<Node>&asg, int k, int a, int b, int n1, int n2) {
	qtree.resize(1);
	vector<vector<int>>grid(3);
	for (int i = 0; i < asg.size(); i++) {
		if (asg[i].alpha > k)
			grid[2].push_back(i);
		else if (asg[i].beta > k)
			grid[0].push_back(i);
		else
			grid[1].push_back(i);
	}
	vector<int>had(asg.size(), -1);
	for (int i = 0; i < 3; i++) {
		for (auto j : grid[i])
			had[j] = 0;
		for (auto j : grid[i]) {
			if (had[j] != 0)
				continue;
			vector<int>bfs;
			bfs.push_back(j);
			had[j] = 1;
			for(auto k:bfs)
				for(auto nbr:asg[k].neighbor)
					if (had[nbr] == 0) {
						had[nbr] = 1;
						bfs.push_back(nbr);
					}
			if (bfs.size() < 16)
				for (auto k : bfs)
					qtree[0].leaf.push_back(k);
			else {
				treeNode* tnode = new treeNode;
				tnode->leaf = bfs;
				tnode->father =0;
				qtree[0].children.push_back(qtree.size());
				qtree.push_back(*tnode);
			}
		}
	}
	had.clear();
	for (int i = 0; i < qtree.size(); i++)
		gridding(asg, i);
	vector<int>order(asg.size(), -1);
	vector<int>father(asg.size(), -1);
	num = 0;
	gf.resize(asg.size());
	tf.resize(qtree.size());
	pre(0, order, father);
	for (int i = 0; i < asg.size(); i++) {
		gf[order[i]] = asg[i];
		for (int j = 0; j < gf[order[i]].neighbor.size(); j++)
			gf[order[i]].neighbor[j] = order[gf[order[i]].neighbor[j]];
	}
	order.clear();
	uIndex.resize(n1 + 1);
	vIndex.resize(n2 + 1);
	for (int i = 0; i < gf.size(); i++) {
		for (auto u : gf[i].left)
			uIndex[u].push_back(i);
		for (auto v : gf[i].right)
			vIndex[v].push_back(i);
	}
	connection(father);
	grid.clear();
	order.clear();
	father.clear();
	qtree.clear();
}
void MHG::gridding(vector<Node>& asg, int c) {
	int minAlpha = asg[qtree[c].leaf[0]].alpha, maxAlpha = asg[qtree[c].leaf[0]].alpha;
	int minBeta = asg[qtree[c].leaf[0]].beta, maxBeta = asg[qtree[c].leaf[0]].beta;
	for (int i = 0; i < qtree[c].leaf.size(); i++) {
		if (minAlpha > asg[qtree[c].leaf[i]].alpha)
			minAlpha = asg[qtree[c].leaf[i]].alpha;
		if (maxAlpha < asg[qtree[c].leaf[i]].alpha)
			maxAlpha = asg[qtree[c].leaf[i]].alpha;
		if (minBeta > asg[qtree[c].leaf[i]].beta)
			minBeta = asg[qtree[c].leaf[i]].beta;
		if (maxBeta < asg[qtree[c].leaf[i]].beta)
			maxBeta = asg[qtree[c].leaf[i]].beta;
	}
	qtree[c].x1 = minAlpha;
	qtree[c].y1 = minBeta;
	qtree[c].x2 = maxAlpha;
	qtree[c].y2 = maxBeta;
	if (qtree[c].leaf.size() < 8)
		return ;
	if (maxAlpha - minAlpha <= 1 || maxBeta - minBeta <= 1)
		return ;
	vector<int>temp;
	int xmid = (minAlpha + maxAlpha) / 2;
	int ymid = (minBeta + maxBeta) / 2;
	vector<bool>inThis(asg.size(), false);
	vector<bool>visited(asg.size(), false);
	for (int i = 0; i < qtree[c].leaf.size(); i++)
		inThis[qtree[c].leaf[i]] = true;
	int tempX1, tempX2, tempY1, tempY2;
	for (int i = 0; i < qtree[c].leaf.size(); i++) {
		int s = qtree[c].leaf[i];
		if (visited[s])
			continue;
		visited[s] = true;
		if (asg[s].alpha >= xmid && asg[s].beta >= ymid) {
			tempX1 = xmid; tempX2 = maxAlpha;
			tempY1 = ymid; tempY2 = maxBeta;
		}
		else if (asg[s].alpha <= xmid && asg[s].beta >= ymid) {
			tempX1 = minAlpha; tempX2 = xmid;
			tempY1 = ymid; tempY2 = maxBeta;
		}
		else if (asg[s].alpha <= xmid && asg[s].beta <= ymid) {
			tempX1 = minAlpha; tempX2 = xmid;
			tempY1 = minBeta; tempY2 = ymid;
		}
		else {
			tempX1 = xmid; tempX2 = maxAlpha;
			tempY1 = minBeta; tempY2 = ymid;
		}
		treeNode* t = new treeNode;
		t->leaf.push_back(s);
		for (int j = 0; j < t->leaf.size(); j++) {
			for (int k = 0; k < asg[t->leaf[j]].neighbor.size(); k++) {
				int par = asg[t->leaf[j]].neighbor[k];
				if (inThis[par] && !visited[par])
					if (asg[par].alpha >= tempX1 && asg[par].alpha <= tempX2)
						if (asg[par].beta >= tempY1 && asg[par].beta <= tempY2) {
							visited[par] = true;
							t->leaf.push_back(par);
						}
			}
		}
		if (t->leaf.size() < 8)
			for (auto w : t->leaf)
				temp.push_back(w);
		else {
			t->father = c;
			qtree[c].children.push_back(qtree.size());
			qtree.push_back(*t);
		}
	}
	swap(qtree[c].leaf, temp);
}
void MHG::pre(int c, vector<int>& order, vector<int>& father) {
	tf[c].x1 = qtree[c].x1;
	tf[c].y1 = qtree[c].y1;
	tf[c].x2 = qtree[c].x2;
	tf[c].y2 = qtree[c].y2;
	tf[c].begin = num;
	for (auto i : qtree[c].leaf) {
		father[num] = c;
		order[i] = num++;
	}
	for (auto i : qtree[c].children)
		pre(i, order, father);
	tf[c].end = num;
}
void MHG::connection(vector<int>&father) {
	for (int i = 0; i < tf.size(); i++) {
		set<int>temp;
		for (int j = tf[i].begin; j < tf[i].end; j++)
			for (auto nbr : gf[j].neighbor)
				if (nbr < tf[i].begin || nbr >= tf[i].end)
					temp.insert(nbr);
		for (auto j : temp)
			tf[i].slide.push_back(j);
	}
	stairs.resize(gf.size());
	for (int i = 0; i < gf.size();i++) {
		vector<int>temp;
		for (auto j : gf[i].neighbor) {
			if (gf[j].alpha >= gf[i].alpha && gf[j].beta >= gf[i].beta && tf[father[j]].x1 >= gf[i].alpha && tf[father[j]].y1 >= gf[i].beta) {
				int p = father[j];
				while (tf[qtree[p].father].x1 >= gf[i].alpha && tf[qtree[p].father].y1 >= gf[i].beta && qtree[p].father != 0)
					p = qtree[p].father;
				stairs[i].push_back(p);
			}
			else
				temp.push_back(j);
		}
		swap(temp, gf[i].neighbor);
		temp.clear();
	}
	for (int i = 0; i < tf.size(); i++) {
		vector<int>temp;
		for (auto j : tf[i].slide) {
			if (gf[j].alpha >= tf[i].x1 && gf[j].beta >= tf[i].y1 && tf[father[j]].x1 >= tf[i].x1 && tf[father[j]].y1 >= tf[i].y1) {
				int p = father[j];
				while (tf[qtree[p].father].x1 >= tf[i].x1 && tf[qtree[p].father].y1 >= tf[i].y1&& qtree[p].father!=0)
					p = qtree[p].father;
				tf[i].neighbor.push_back(p);
			}
			else
				temp.push_back(j);
		}
		swap(temp, tf[i].slide);
		temp.clear();
	}
}
void MHG::Query(int alpha, int beta, int u, bool isU, vector<bool>& leftResult, vector<bool>& rightResult) {
	vector<bool>had(gf.size(), false);
	queue<int>Q;
	if (isU) {
		for (auto i : uIndex[u]) {
			had[i] = true;
			if (gf[i].alpha >= alpha && gf[i].beta >= beta) 
				Q.push(i);
		}
	}
	else {
		for (auto i : vIndex[u]) {
			had[i] = true;
			if (gf[i].alpha >= alpha && gf[i].beta >= beta) 
				Q.push(i);
		}
	}
	while (!Q.empty()) {
		int w = Q.front();
		Q.pop();
		for (auto i : gf[w].neighbor)
			if (!had[i]) {
				had[i] = true;
				if (gf[i].alpha >= alpha && gf[i].beta >= beta) 
					Q.push(i);
			}
		for (auto j : gf[w].left)
			leftResult[j] = true;
		for (auto j : gf[w].right)
			rightResult[j] = true;
	}
}
void MHG::QueryViaMHG(int alpha, int beta, int u, bool isU, vector<bool>& leftResult, vector<bool>& rightResult) {
	vector<bool>gHad(gf.size(), false), mHad(tf.size(), false);
	queue<int>gQ, mQ;
	if (isU) {
		for (auto i : uIndex[u]) {
			gHad[i] = true;
			if (gf[i].alpha >= alpha && gf[i].beta >= beta) 
				gQ.push(i);
		}
	}
	else {
		for (auto i : vIndex[u]) {
			gHad[i] = true;
			if (gf[i].alpha >= alpha && gf[i].beta >= beta) 
				gQ.push(i);
		}
	}
	while (!gQ.empty() || !mQ.empty()) {
		while (!gQ.empty()) {
			int w = gQ.front();
			gQ.pop();
			for (auto i : gf[w].left)
				leftResult[i] = true;
			for (auto i : gf[w].right)
				rightResult[i] = true;
			for (auto i : stairs[w]) 
				if (!mHad[i]) {
					mHad[i] = true;
					mQ.push(i);
				}
			for (auto i : gf[w].neighbor)
				if (!gHad[i]) {
					gHad[i] = true;
					if (gf[i].alpha >= alpha && gf[i].beta >= beta)
						gQ.push(i);
				}
		}
		while (!mQ.empty()) {
			int m = mQ.front();
			mQ.pop();
			for (int node = tf[m].begin; node < tf[m].end; node++)
				if (!gHad[node]) {
					gHad[node] = true;
					for (auto i : gf[node].left)
						leftResult[i] = true;
					for (auto i : gf[node].right)
						rightResult[i] = true;
				}
			for (auto w : tf[m].slide)
				if (!gHad[w]) {
					gHad[w] = true;
					if (gf[w].alpha >= alpha && gf[w].beta >= beta)
						gQ.push(w);
				}
			for (auto nbr : tf[m].neighbor)
				if (!mHad[nbr]) {
					mHad[nbr] = true;
					mQ.push(nbr);
				}
		}
	}
}