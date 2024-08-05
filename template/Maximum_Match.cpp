/*
二分图最小路径覆盖 = 最大独立集 = 总节点数 - 最大匹配数，最小点覆盖 = 最大匹配数。
任意图中，最大独立集 + 最小点覆盖集 = V ， 最大团 = 补图的最大独立集。
独立集：点集，两两间没有边相连
二分图最大匹配：在图中选出一些边，使得这些边没有公共顶点，且边的数量最大
团：完全子图
最小点覆盖：选最少的点，满足每条边至少有一个端点被选。
最小路径覆盖：找出最少的路径，使得这些路径经过了所有的点。
最小路径覆盖分为最小不相交路径覆盖和最小可相交路径覆盖。
*/


namespace maximum_matching{
vector<int> G[N];
// col[i]i在左右  lef[i] 右边的节点i匹配左边节点lef[i]
int n, m, col[N], vis[N], lef[N];
void dfs(int u)
{
	for(auto v:G[u]){
		if(!col[v]){
			col[v] = 3-col[u];
			dfs(v);
		}
	}
}
bool find(int u)
{
	for(auto v:G[u]) if(!vis[v]){
		vis[v]=1;
		if(!lef[v] || find(lef[v])){
			lef[v]=u;
			return true;
		}
	}
	return false;
}
int match()	// 二分图最大匹配
{
	int ans = 0;
	for(int i = 1; i <= n; i++){	// 跑二分图一边
		if(col[i]==1){
			for(int j = 1; j <= n; j++)vis[j]=0;
			if(find(i))ans++;
		}
	}
	return ans;
}
};
