// R >= L, 左边完全匹配
// 最小权完备匹配
// 板子为最小权，改成最大全边权取负
/**
 * @param a 图，邻接矩阵，下标0开始
 * @return pair(v,ans). v为匹配的值, ans为左边的每个点匹配右边的哪个点
 * 
 * 用法：\Daimayuan\graphmid.cpp 3971行  二分图最大权匹配
*/
template<class T>
pair<T, vector<int>> hungarian(const vector<vector<T>> &a) {
	if (a.empty()) return {0, {}};
	int n = a.size() + 1, m = a[0].size() + 1;
	vector<T> u(n), v(m); // 顶标
	vector<int> p(m), ans(n - 1);
	for (int i = 1; i < n; i++) {
		p[0] = i;
		int j0 = 0;
		vector<T> dist(m, numeric_limits<T>::max());
		vector<int> pre(m, -1);
		vector<bool> done(m + 1);
		do { // dijkstra
			done[j0] = true;
			int i0 = p[j0], j1;
			T delta = numeric_limits<T>::max();
			for (int j = 1; j < m; j++) if (!done[j]) {
				auto cur = a[i0 - 1][j - 1] - u[i0] - v[j];
				if (cur < dist[j]) dist[j] = cur, pre[j] = j0;
				if (dist[j] < delta) delta = dist[j], j1 = j;
			}
			for (int j = 0; j < m; j++) {
				if (done[j]) u[p[j]] += delta, v[j] -= delta;
				else dist[j] -= delta;
			}
			j0 = j1;
		} while (p[j0]);
		while (j0) { // update alternating path
			int j1 = pre[j0];
			p[j0] = p[j1], j0 = j1;
		}
	}
	for (int j = 1; j < m; j++) {
		if (p[j]) ans[p[j] - 1] = j - 1;
	}
	return {-v[0], ans}; // min cost
}