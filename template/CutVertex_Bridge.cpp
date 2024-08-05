namespace cut_vertex{
int dfn[N], low[N], idx, iscut[N];
void dfs(int u, int fa)
{
	dfn[u] = low[u] = ++idx;
	int child = 0;	// 孩子数
	for(int v : G[u]){
		if(!dfn[v]){
			dfs(v, u);
			child++;
			low[u] = min(low[u], low[v]);
			if(low[v] >= dfn[u])iscut[u] = 1;	// 子树v跳不出去 u是割点
		}
		else if(v != fa){
			low[u] = min(low[u], dfn[v]);	// 不能跳多步 只能dfn[v]
		}
	}
	if(u == 1 && child <= 1)iscut[u] = 0;
}
}
namespace bridge{
int dfn[N], low[N], idx;
vector<int> bridge;
void dfs(int u, int id)	// 通过id边下来
{
	dfn[u] = low[u] = ++idx;
	for(auto [v, idd] : G[u]){
		if(!dfn[v]) dfs(v, idd);	// 儿子
		if(idd != id) low[u] = min(low[u], low[v]);	// 不是id的边
	}
	if(dfn[u] == low[u] && id != 0){
		bridge.push_back(id);
		/* 可以和强连通分量一样把边双找出来 */
	}
}
}