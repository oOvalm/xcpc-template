/*
点权
*/
namespace treepoufen{
vector<int> G[N];
// hson重孩  top[u]点u所在重链的最上面点 dep深度 fa父亲 siz子树大小
int n, hson[N], top[N], dep[N], fa[N], siz[N];
// dfs2才做dfs序
int l[N], r[N], dfn, pos[N];
void dfs1(int u, int f=0)
{
	siz[u] = 1;
	hson[u] = -1;
	fa[u] = f;
	for(auto v : G[u]) if(v != f){
		dep[v] = dep[u] + 1;
		dfs1(v, u);
		siz[u] += siz[v];
		if(hson[u] == -1 || siz[hson[u]] < siz[v])
			hson[u] = v;
	}
}
void dfs2(int u, int t = 1)  // t 重链头
{
	top[u] = t;
	l[u] = ++dfn;   // 同一重链的dfs序是连续的
	pos[dfn] = u;
	if(hson[u] != -1)	//如果有重孩子
		dfs2(hson[u], t);	// 同一条重链
	
	for(auto v:G[u]){
		if(v != fa[u] && v != hson[u])
			dfs2(v, v);	// 通过轻边连接
	}
	r[u] = dfn;
}
int LCA(int u, int v)   // lca过程中可以整段重链信息一起查询
{
	// 跳到一条重链
	while(top[u] != top[v])
	{
		if(dep[top[u]] < dep[top[v]]){
            v = fa[top[v]];
        }
        else{
            u = fa[top[u]];
        }
	}
	if(dep[u] < dep[v])return u;
	return v;
}
Node query(int u, int v)
{
	Node ans = {-INF, 0};
	while(top[u] != top[v]){
		if(dep[top[u]] < dep[top[v]]){
			assert(l[top[v]] <= l[v]);
			ans = ans + query(1,1,n,l[top[v]],l[v]);
			v = fa[top[v]];
		}
		else{
			assert(l[top[u]] <= l[u]);
			ans = ans + query(1,1,n,l[top[u]],l[u]);
			u = fa[top[u]];
		}
	}
	ans = ans + query(1,1,n,min(l[u],l[v]),max(l[u],l[v]));
	return ans;
}
};