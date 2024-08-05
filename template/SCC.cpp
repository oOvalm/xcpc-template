namespace SCC{
int dfn[N], low[N], idx, ins[N], cnt, bel[N];
stack<int> st;
void dfs(int u)
{
	dfn[u] = low[u] = ++idx;
	ins[u] = 1;
	st.push(u);
	for(auto v:G[u]){
		if(!dfn[v])dfs(v);
		if(ins[v])low[u] = min(low[u], low[v]);
	}
	if(low[u] == dfn[u]){
		cnt++;
		while(true){
			int x = st.top(); st.pop();
			ins[x]=0;
			bel[x] = cnt;
			if(x==u)break;
		}
	}
}
}