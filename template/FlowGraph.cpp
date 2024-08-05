const int V = 2e4+100;
const int E = 5e4+100;
struct FlowGraph{
	int s, t, vtot, etot;
	int head[V], dis[V], cur[V];
	struct edge{
		int v, nxt;
		ll cap;
	}e[E<<1];
	void add(int u, int v, ll f, ll ff = 0){
		e[etot] = {v, head[u], f}; head[u] = etot++;
		e[etot] = {u, head[v], ff}; head[v] = etot++;
	}
	void init(int _s, int _t, int _vtot){
		s = _s, t = _t, vtot = _vtot;
		etot = 0;
		For(i,1,vtot)head[i] = -1;
	}
	bool bfs()
	{
		For(i,1,vtot){
			dis[i] = 0;
			cur[i] = head[i];
		}
		queue<int> q;
		q.push(s); dis[s] = 1;
		while(q.size()){
			int u = q.front(); q.pop();
			for(int i = head[u]; ~i; i = e[i].nxt){
				if(e[i].cap && !dis[e[i].v]){
					int v = e[i].v;
					dis[v] = dis[u] + 1;
					if(v == t)return true;
					q.push(v);
				}
			}
		}
		return false;
	}
	ll dfs(int u, ll mx)
	{
		if(u == t)return mx;
		ll flow = 0;
		for(int i = cur[u]; ~i; cur[u] = i = e[i].nxt){
			if(e[i].cap && dis[e[i].v] == dis[u] + 1){
				int f = dfs(e[i].v, min(mx, e[i].cap));
				e[i].cap -= f;
				e[i^1].cap += f;
				mx -= f;
				flow += f;
				if(!mx)break;
			}
		}
		if(!flow)dis[u] = -1;
		return flow;
	}
	ll dinic()
	{
		int flow = 0;
		while(bfs())
			flow += dfs(s, numeric_limits<ll>::max());
		return flow;
	}
};