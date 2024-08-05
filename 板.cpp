#include<bits/stdc++.h>
#define ll long long
#define ull unsigned long long
#define int ll
#define fr first
#define se second
#define INF 0x3f3f3f3f    //1e9
#define For(i,a,b) for(int i = (a); (i) <= (b); ++i)
#define Rep(i,a,b) for(int i = (a); (i) >= (b); --i)
using namespace std;
typedef pair<int,int> pii;
#ifdef OVAL
const ll N = 2e3+10;
#else
const ll N = 2e5+10;
#endif


template <typename T> void inline read(T &x) {              //快读 用于读整数
    int f = 1; x = 0; char s = getchar();
    while (s < '0' || s > '9') { if (s == '-') f = -1; s = getchar();}
    while (s <= '9' && s >= '0') x = (x << 1) + (x << 3) + (s ^ 48), s = getchar();      //x*10+s-'0'
    x *= f;
}


mt19937_64 rnd(chrono::steady_clock::now().time_since_epoch().count());
// 树状数组
template<typename T>
struct BIT{
	T tre[N];
	BIT(){memset(tre,0,sizeof(tre));}
	int lowbit(int x){return x&-x;}
	void add(int x,T d)
	{
		while(x<N){
			tre[x] += d;
			x += lowbit(x);
		}
	}
	T query(int x)
	{
		T res = 0;
		while(x){
			res += tre[x];
			x -= lowbit(x);
		}
		return res;
	}
};

namespace Mat{
const int N=100;
struct Matrix{
	vector<vector<int>> a;
	int n,m;
	Matrix(){a.resize(1, vector<int>(1, 0)); n = m = 1;}
	Matrix(int _n, int _m){
		n = _n, m = _m;
		a.resize(n, vector<int>(m, 0));
	}
	Matrix(Matrix const& x){
		n = x.n, m = x.m;
		a.resize(n, vector<int>(m, 0));
		for(int i = 0; i < n; i++)
			for(int j = 0; j < m; j++)
				a[i][j] = x.a[i][j];
	}

	void unit(){
		for(int i = 0; i < n; i++)
			for(int j = 0; j < m; j++){
				if(i == j)a[i][j] = 1;
				else a[i][j] = 0;
			}
	}
	Matrix operator = (int x){
		Matrix res(n, m);
		for(int i = 0; i < n; i++)
			for(int j = 0; j < m; j++)
				if(i == j)res.a[i][j] = x;
				else res.a[i][j] = 0;
		return res;
	}
	Matrix operator = (Matrix x){
		this->a = x.a;
		this->n = x.n;
		this->m = x.m;
		return *this;
	}
	Matrix operator + (Matrix x)const{
		assert(n==x.n && m==x.m);
		Matrix res(n, m);
		for(int i = 0; i < n; i++)
			for(int j = 0; j < m; j++)
				res.a[i][j] = a[i][j] + x.a[i][j];
		return res;
	}
	Matrix operator - (Matrix x)const{
		assert(n==x.n && m==x.m);
		Matrix res(n, m);
		for(int i = 0; i < n; i++)
			for(int j = 0; j < m; j++)
				res.a[i][j] = a[i][j] - x.a[i][j];
		return res;
	}
	Matrix operator * (Matrix x)const{
		assert(m==x.n);
		Matrix res(n, x.m);
		for(int i = 0; i < n; i++)
			for(int k = 0; k < m; k++) if(a[i][k])
				for(int j = 0; j < x.m; j++) if(x.a[k][j])
					res.a[i][j] += a[i][k] * x.a[k][j];
		return res;
	}
	Matrix operator += (Matrix x){(*this) = (*this)+x; return *this;}
	Matrix operator -= (Matrix x){(*this) = (*this)-x; return *this;}
	Matrix operator *= (Matrix x){(*this) = (*this)*x; return *this;}
};
}

/* 线性基 */
namespace Linear_Basis{
const int B = 63;	// 最高位数
struct linear_basis
{
	ll num[B+1];
	bool insert(ll x)	// 将x插进线性基 true插入成功
	{
		for(int i = B; i >= 0; i--) if(x & (1ll<<i)) {
			if(num[i] == 0){
				num[i] = x;
				return true;
			}
			x ^= num[i];
		}
		return false;
	}
	ll querymin(ll x)	// 如果x为0 线性基可以表示x
	{
		for(int i = B; i >= 0; i--) {
			x = min(x, x ^ num[i]);
		}
		return x;
	}
};
}


/****马拉车 manacher*****/
int p[N];// 回文半径 [i-p[i]+1,i+p[i]-1]是回文串
int manacher(string& t)
{
	string s=" #";//每个字符间插个原串不存在的字符 头尾插一个
	for(int i = 0; i < t.size(); i++){
		s += t[i];
		s += "#";
	}
	int n=s.size();
	n--;
	For(i,1,n)p[i]=0;
	int m=0,r=0;
	for(int i = 1; i <= n; i++){
		if(i>r)p[i]=1;
		else p[i] = min(p[2*m-i], r-i+1); // i+p[i]-1 <= r
		while(i-p[i]>0 && i+p[i]<=n && s[i-p[i]]==s[i+p[i]]){//暴力扩
			p[i]++;
		}
		if(i+p[i]-1>r){//更新右端点和中点
			m=i;r=i+p[i]-1;
		}
	}
	int res = 0;//最大回文半径 对应插入后回文串长度为2*res-1  对应原串 res-1
	For(i,1,n){
		res = max(res,p[i]);
	}
	return res-1;
}

/***后缀数组***/
namespace SA{
/**
 * n 长度
 * m 字符集大小
 * s 字符串	0 base
 * sa 字典序第i小是哪个后缀
 * rk 后缀i的排名
 * ht lcp(sa[i], sa[i-1]);
*/
int sa[N], rk[N], ht[N];
void getSA(string& s, int m=128)
{
	static int xx[N], yy[N], c[N];
	int *x=xx, *y=yy;
	int n = s.size();
	for(int i = 0; i < m; i++) c[i] = 0;
	for(int i = 0; i < n; i++) c[x[i] = s[i]]++;
	for(int i = 1; i < m; i++) c[i] += c[i - 1];
	for(int i = n - 1; i >= 0; i--) sa[--c[x[i]]] = i;
	for(int k = 1; k < n; k <<= 1) {
		int p = 0;
		for(int i = n - 1; i >= n - k; i--) y[p++] = i;
		for(int i = 0; i < n; i++) if(sa[i] >= k) y[p++] = sa[i] - k;
		for(int i = 0; i < m; i++) c[i] = 0;
		for(int i = 0; i < n; i++) c[x[y[i]]]++;
		for(int i = 1; i < m; i++) c[i] += c[i - 1];
		for(int i = n - 1; i >= 0; i--) sa[--c[x[y[i]]]] = y[i];
		swap(x, y);
		p = 1; x[sa[0]] = 0; y[n] = -1;
		for(int i = 1; i < n; i++) {
			if(y[sa[i - 1]] == y[sa[i]] && y[sa[i - 1] + k] == y[sa[i] + k])
				x[sa[i]] = p - 1;
			else 
				x[sa[i]] = p++;
		}
		if(p == n)break;
		m = p;
	}
	// 计算ht
	for(int i = 0; i < n; i++)rk[sa[i]] = i;
	int k = 0;
	for(int i = 0; i < n; i++) {
		k = max(k-1, 0ll);
		if(rk[i] == 0)continue;
		int j = sa[rk[i] - 1];
		while(s[i + k] == s[j + k]) k++;
		ht[rk[i]] = k;
	}
}
}


namespace dsuOntree{
vector<int> G[N];
// siz 子树大小 lef,rig dfs序子树区间 pos[u]u在dfs序的位置， heavy重儿子 idx dfs序
int siz[N], lef[N], pos[N], rig[N], heavy[N], dfn;
int n;
// 预处理信息
void dfs1(int u, int fa)
{
	siz[u]=1;	//子树大小
	lef[u]=++dfn;pos[dfn]=u;// dfs序
	for(auto v : G[u]) if(v != fa){
		dfs1(v, u);
		siz[u] += siz[v];
		if(heavy[u] == -1 || siz[v] > siz[heavy[u]])	// 重儿子
			heavy[u] = v;
	}
	rig[u]=dfn;
}
void dfs2(int u, int fa, int keep)
{
	for(auto v : G[u])		// 先计算所有轻儿子 并且不保存信息
		if(v != fa && v!= heavy[u])
			dfs2(v, u, 0);

	if(heavy[u] != -1)// 不是叶子   
		dfs2(heavy[u], u, 1);	// 计算重儿子 保存重儿子信息
	
	auto add = [&](int u){
		/*添加函数*/
	};
	auto del = [&](int u){
		/*删除函数*/
	};

	// 轻儿子加入重儿子  暴力将轻儿子信息加入
	for(auto v : G[u])
		if(v != fa && v!= heavy[u]){
			For(i,lef[v],rig[v]){	// 用dfs序减少递归的时间消耗
				add(pos[i]);
			}
		}
	//加入自己
	add(u);
	//  如果不保存信息 就清空删除
	if(!keep){
		// maxcnt = sumcnt = 0;
		// For(i,lef[u],rig[u])
		// 	del(pos[i]);
	}
}
};



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
void dfs2(int u, int t)  // t 重链头
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
};


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

/***二分图最大匹配***/
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




/*
模拟退火
可调
初始温度 T 
结束温度 eps
降温系数 一般为接近1的数字
参考：
1e-1精度：下界1e-3或1e-4 降温系数0.8
1e-5精度：下界1e-7或1e-8 降温系数0.99
*/
namespace Simulated_Annealing
{
struct Point
{
	double x, y, z;
	Point(){}
	Point(double _x, double _y, double _z):x(_x),y(_y),z(_z){}
	bool operator == (Point b){return abs(x-b.x)<1e-6 && abs(y-b.y)<1e-6 && abs(z-b.z)<1e-6;}
};
double dis(Point a, Point b)
{
	double dx = a.x-b.x, dy = a.y-b.y, dz = a.z-b.z;
	return sqrt(dx*dx+dy*dy+dz*dz);
}
Point p[N];
int n;
double maxdis(Point a)
{
	double res = 0;
	for(int i = 1; i <= n; i++){
		res = max(res, dis(a, p[i]));
	}
	return res;
}
bool out(Point a)
{
	return abs(a.x) > 1e5 || abs(a.y) > 1e5 || abs(a.z) > 1e5;
}
Point ansp;
double ans;
void SA()	// 模拟退火
{
	double T = 550, eps = 1e-10;
	while(T>eps){	// 结束温度
		Point tmp(ansp);
		// 随机到的点在答案范围内 (也可以不加，温度高容易随机到边界外 可能会卡死)
		while(tmp == ansp || out(tmp)){		
			tmp.x = ansp.x + T*(2*rand()-RAND_MAX);		// 根据当前温度随机选
			tmp.y = ansp.y + T*(2*rand()-RAND_MAX);
			tmp.z = ansp.z + T*(2*rand()-RAND_MAX);
		}
		double tmpans = maxdis(tmp);	// 计算答案
		if(tmpans - ans < 1e-6){		// 比当前答案更优 更新答案
			ans = tmpans;
			ansp=tmp;
		}
		// 比当前答案不优 有一定概率选择不优的答案
		else if(exp(-abs(tmpans-ans)/T) * RAND_MAX > rand()){	
			ans = tmpans;
			ansp=tmp;
		}
		T *= 0.998;	//降温系数
	}
}
void tuihuo()
{
	ansp = Point(0,0,0);	// 随便指定一个答案
	ans = maxdis(ansp);
	// 做多几次退火  
	while(clock() < 1000)SA();		// 卡时间
	// for(int i = 0; i < 100; i++)SA();	// 固定次数
}
};



namespace Flow{

const int V = 1010;		// 点数
const int E = 100010;	// 边数
template<typename T>
class FlowGraph{
	int s, t, n, m;
	int head[V];
	int dis[V], cur[V];
	struct Edge{
		int v, nxt;
		T cap;
	}e[E<<1];		// 边编号从0开始
	bool bfs()
	{
		For(i,1,n){
			dis[i] = 0;
			cur[i] = head[i];
		}
		queue<int> q;
		q.push(s); dis[s] = 1;
		while(!q.empty()){
			int u = q.front(); q.pop();
			for(int i = head[u]; ~i; i = e[i].nxt){
				if(e[i].cap && !dis[e[i].v]){	// 还有流量才算
					int v = e[i].v;
					dis[v] = dis[u] + 1;
					if(v == t) return true;	// 存在最短路
					q.push(v);
				}
			}
		}
		return false;
	}
	T dfs(int u, T mx)
	{
		if(u == t)return mx;
		T flow = 0;
		// cur前面的边不用再跑了
		for(int i = cur[u]; ~i; cur[u] = i = e[i].nxt){
			if(e[i].cap && dis[e[i].v] == dis[u] + 1){
				T f = dfs(e[i].v, min(mx, e[i].cap));
				e[i].cap -= f;
				e[i^1].cap += f;
				mx -= f;
				flow += f;
				if(!mx)break;// 流完了
			}
		}
		if(!flow)dis[u] = -1;// 流不出去
		return flow;
	}
	public:
	T dinic(){
		T flow = 0;
		while(bfs())
			flow += dfs(s, numeric_limits<T>::max());	// ll 最大值
		return flow;
	}
	void init(int _s, int _t, int _vtot){
		s = _s, t = _t, n = _vtot;
		m = 0;
		For(i, 1, n) head[i] = -1;
	}
	void add(int u, int v, T f){
		e[m] = {v, head[u], f}; head[u] = m++;
		e[m] = {u, head[v], 0}; head[v] = m++;
	}

};

}