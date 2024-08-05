#include<bits/stdc++.h>
using namespace std;
typedef pair<double,int> pdi;
const int N = 1e5+10, mod = 39989, inf = 0x3f3f3f3f;
#define EPS 1e-9
struct Line{
	// 斜率为0时直接 k=0 b=max(y)
	double k, b;	// y = kx+b;
}line[N];
double f(int id, int x){return line[id].k*x + line[id].b;}
int sign(double a){return a < -EPS ? -1 : a > EPS;}
int cmp(double a, double b){return sign(a-b);}
struct LiChao{
	int tr[N<<2];	// tr[u]为u节点区间中电最优线段的下标
	#define lson (u<<1)
	#define rson (u<<1|1)
	// 插入下标为id的直线
	void insertLine(int u, int l, int r, int id){
		int &uid = tr[u], mid = (l+r)>>1;
		int midcmp = cmp(f(id, mid), f(uid, mid));
		// 如果新插入线段更优 || 相等&&编号更小，则u的最优线段为插入线段，原线段下传
		if(midcmp > 0 || (midcmp == 0 && id < uid))
			swap(id, uid);
		int lcmp = cmp(f(id, l), f(uid, l)), rcmp = cmp(f(id, r), f(uid, r));
		// 递归将id插入可能是最优线段的子节点中，最多只有1个if是成立的
		if(lcmp > 0 || (lcmp == 0 && id < uid)) insertLine(lson, l, mid, id);
		if(rcmp > 0 || (rcmp == 0 && id < uid)) insertLine(rson, mid+1, r, id);
	}
	// 插入下标为id 区间为[il, ir]的线段
	void insertSeg(int u, int l, int r, int il, int ir, int id){
		if(il <= l && r <= ir){	// 当前区间被线段完全覆盖，插
			insertLine(u, l, r, id);
			return;
		}
		int mid = (l+r)>>1;
		if(il <= mid) insertSeg(lson, l, mid, il, ir, id);
		if(ir > mid) insertSeg(rson, mid+1, r, il, ir, id);
	}
	// pair的比较函数
	pdi pmax(pdi x, pdi y){
		if(cmp(x.first, y.first) > 0)return x;
		else if(cmp(x.first, y.first) < 0)return y;
		else if(x.second < y.second)return x;
		else return y;
	}
	// 查询x=x的直线上最高的点在那个线段上，返回{y, id}
	pdi query(int u, int l, int r, int x){
		if(r < x || l > x) return {-inf, 0};	// 线段为空返回0
		if(l==r) return {f(tr[u], x), tr[u]};
		int mid = (l+r)>>1;
		double value = f(tr[u], x);
		return pmax({value, tr[u]}, 
			pmax(query(lson, l, mid, x), query(rson, mid+1, r, x)));
	}
}T;
signed main()
{
	ios::sync_with_stdio(false),cin.tie(0),cout.tie(0);
	int q, lastans = 0, cnt = 0;
	cin >> q;
	while(q--){
		int op, k, x[2], y[2];
		cin >> op;
		if(op==0){
			cin >> k;
			k = (k+lastans-1)%mod+1;
			lastans = T.query(1,1,mod,k).second;
			cout << lastans << '\n';
		}
		else{
			const int M = 1000000000;
			for(int i = 0; i <= 1; i++){
				cin >> x[i] >> y[i];
				x[i] = (x[i]+lastans-1)%mod+1;
				y[i] = (y[i]+lastans-1)%M+1;
			}
			if(x[0] > x[1]){
				swap(x[0], x[1]);
				swap(y[0], y[1]);
			}
			++cnt;
			if(x[1] == x[0]){	// 特判斜率不存在
				line[cnt] = {0, (double)max(y[1], y[0])};
			}
			else{
				double kk = 1.0*(y[1]-y[0])/(x[1]-x[0]);
				line[cnt] = {kk, y[0]-kk*x[0]};
			}
			T.insertSeg(1,1,mod,x[0],x[1],cnt);
		}
	}
	return 0;
}