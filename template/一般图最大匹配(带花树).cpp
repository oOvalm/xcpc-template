// 带花树算法，求一般图最大匹配 O(n m^2) 不满
// 洛谷 n=1000 m=5e4   < 50ms
int p[N];	// 并查集
int find(int x){return x==p[x]?x:p[x]=find(p[x]);}
int match[N];	// 和 i 匹配的点
int pre[N], col[N];	// i的前驱, i的颜色, 
int tim, dfn[N];	// tim要全局，或者清空dfn(只是记个时间, 不是真的dfn)
bool augment(int s)	// 找增广路
{
	For(i,1,n){
		pre[i] = col[i] = 0;
		p[i] = i;
	}
	queue<int> q;	// 只有白色点才会入队

	auto LCA = [&](int x, int y){
		++tim;
		x = find(x), y = find(y);	// 因为带了并查集，每个点均摊O(1)
		while(dfn[x] != tim){	// 两个点都经过的点就是lca
			dfn[x] = tim;
			x = find(pre[match[x]]);// lca一定是白的，可以隔着跳
			if(y)swap(x, y);	// 轮流跳
		}
		return x;
	};

	auto blossom = [&](int x, int y, int lc){
		while(find(x) != lc){
			pre[x] = y, y = match[x];
			if(col[y] == 2){
				col[y] = 1;
				q.push(y);	// 染黑入队
			}
			if(find(x) == x) p[x] = lc;
			if(find(y) == y) p[y] = lc;
			x = pre[y];
		}
	};

	q.push(s);
	col[s] = 1;	// 1为白色 2为黑色
	while(q.size()) {	// bfs
		int u = q.front(); q.pop();
		for(int v:G[u]){
			if(find(v) == find(u) || col[v]==2)continue;
			if(!col[v]){	// v未被染色
				pre[v] = u, col[v] = 2;
				if(!match[v]){// v没被匹配, 找到一条增广路(s->v)
					while(v){
						int ppre = match[pre[v]]; // 不能用pre[pre[v]]，在开花改乱了pre
						match[v] = pre[v];
						match[pre[v]] = v;
						v = ppre;
					}
					return true;
				}
				col[match[v]] = 1;
				q.push(match[v]);
			}
			else{	// col[v]==1 找到奇环了
				int lc = LCA(u, v);
				// 把这个环的所有点染成一样的颜色，当成一个点找增广路
				blossom(u, v, lc);	// 缩花
				blossom(v, u, lc);
			}
		}
	}
	return false;
}
signed main()
{
	// 随便匹配一些边
	for(auto [u, v]:e)if(!match[u] && !match[v]){
		match[u] = v, match[v] = u;
		ans++;
	}
	// 类似匈牙利, 每次从一个未匹配的点找增广路
	For(i,1,n)if(!match[i])
		ans += augment(i);
}