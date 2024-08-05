const int V = 410;
const int E = 15010;
struct MinCostFlow{
	int s, t, vtot, etot;
	int head[V];
	int pre[V];		// 记路径
	ll dis[V], flow, cost;
	bool inq[V];
	struct edge{
		int v, nxt;
		ll f, c;
	}e[E<<1];
	void add(int u, int v, ll f, ll c, ll ff = 0){
		e[etot] = {v, head[u], f, c}; head[u] = etot++;
		e[etot] = {u, head[v], ff, -c}; head[v] = etot++;
	}
	void init(int _s, int _t, int _vtot){
		s = _s, t = _t, vtot = _vtot;
		etot = 0;
		For(i,0,vtot)head[i] = -1;
	}
	bool spfa(){
		ll inf = numeric_limits<ll>::max()/2;// ll最大值/2
		For(i,1,vtot){
			dis[i] = inf;
			inq[i] = false;
			pre[i] = -1;
		}
		dis[s] = 0;
		inq[s] = true;
		queue<int> q;
		q.push(s);
		while(!q.empty()){
			int u = q.front();q.pop();
			inq[u] = false;
			for(int i = head[u]; ~i; i = e[i].nxt){
				int v = e[i].v;
				if(e[i].f && dis[v] > dis[u]+e[i].c){
					dis[v] = dis[u]+e[i].c;
					pre[v] = i;
					if(!inq[v]){
						// if(++cnt[v] >= vtot)->有负环	// cnt[]为松弛次数
						inq[v] = true;
						q.push(v);
					}
				}
			}
		}
		return dis[t] != inf;
	}
	void augment(){
		int u = t;
		ll f = numeric_limits<ll>::max();
		while(~pre[u]){	// 找路径上最小值
			f = min(f, e[pre[u]].f);
			u = e[pre[u]^1].v;
		}
		flow += f;
		cost += f*dis[t];
		u = t;
		while(~pre[u]){		// 路劲上减掉流的(反向边加上)
			e[pre[u]].f -= f;
			e[pre[u]^1].f += f;
			u = e[pre[u]^1].v;
		}
	}
	pii solve(){
		flow = cost = 0;
		while(spfa())augment();
		return make_pair(flow, cost);
	}
};