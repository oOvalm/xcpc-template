//	P3380 【模板】二逼平衡树（树套树）
#include<bits/stdc++.h>		// 加上离散化
#define ll long long
#define For(i,a,b) for(int i = a; i <= b; ++i)
using namespace std;
const int N = 5e5;
int lowbit(int x){return x&-x;}
struct Node{
	int ls, rs, val;
}tr[N<<5];
vector<int> li;
int n, q, a[N], rt[N], m, cnt, b[N];
int nodeR[N], nodeL[N], cntR, cntL;
void pushup(int u)
{
	tr[u].val = tr[tr[u].ls].val + tr[tr[u].rs].val;
}
void _modify(int &u, int l, int r, int pos, int val)
{
	if(!u)u = ++cnt;	// 如果本来没点就开
	if(l==r){
		tr[u].val += val;
		return;
	}
	int mid = (l+r)>>1;
	if(pos <= mid) _modify(tr[u].ls, l, mid, pos, val);
	else _modify(tr[u].rs, mid+1, r, pos, val);
	pushup(u);
}
// op=3		a[pos] += val
void add(int pos, int val, int num)
{
	while(pos <= n) {	// 树状数组那些节点的线段树
		_modify(rt[pos], 0, m, num, val);
		pos += lowbit(pos);
	}
}
// op=1		k在[l,r]的排名
int _query_kRank(int l, int r, int k)
{
	if(l==r)return 0; 		// 同数字的不算数
	int mid = (l+r)>>1, lsum = 0;
	if(k<=mid){
		For(i,1,cntR) nodeR[i] = tr[nodeR[i]].ls;
		For(i,1,cntL) nodeL[i] = tr[nodeL[i]].ls;
		return _query_kRank(l, mid, k);
	}
	else{
		For(i,1,cntR)lsum += tr[tr[nodeR[i]].ls].val, nodeR[i] = tr[nodeR[i]].rs;
		For(i,1,cntL)lsum -= tr[tr[nodeL[i]].ls].val, nodeL[i] = tr[nodeL[i]].rs;
		return lsum + _query_kRank(mid+1, r, k);
	}
}
// op=2		第k名
int _query_kth(int l, int r, int k)
{
	if(l==r)return l;
	int mid = (l+r)>>1;
	int lsum = 0;
	For(i,1,cntR)lsum += tr[tr[nodeR[i]].ls].val;
	For(i,1,cntL)lsum -= tr[tr[nodeL[i]].ls].val;
	if(lsum >= k){
		For(i,1,cntL) nodeL[i] = tr[nodeL[i]].ls;
		For(i,1,cntR) nodeR[i] = tr[nodeR[i]].ls;
		return _query_kth(l, mid, k);
	}
	else{
		For(i,1,cntL) nodeL[i] = tr[nodeL[i]].rs;
		For(i,1,cntR) nodeR[i] = tr[nodeR[i]].rs;
		return _query_kth(mid+1, r, k-lsum);
	}
}
void get_Root(int l, int r)
{
	cntR = cntL = 0;
	for(int i = l-1; i; i -= lowbit(i))
		nodeL[++cntL] = rt[i];
	for(int i = r; i; i -= lowbit(i))
		nodeR[++cntR] = rt[i];
}
struct Que{
	int op, l, r, k;
}que[N];
signed main()
{
	ios::sync_with_stdio(false),cin.tie(0),cout.tie(0);
	cin >> n >> q;
	For(i,1,n){
		cin >> a[i];
		li.push_back(a[i]);
	}
	For(i,1,q){
		cin >> que[i].op;
		if(que[i].op == 3)cin >> que[i].l >> que[i].k;
		else cin >> que[i].l >> que[i].r >> que[i].k;
		li.push_back(que[i].k);
	}
	sort(li.begin(), li.end());
	li.erase(unique(li.begin(), li.end()), li.end());
	m = li.size();
	For(i,1,n){
		b[i] = lower_bound(li.begin(), li.end(), a[i])-li.begin()+1;
		add(i, 1, b[i]);
	}
	For(i,1,q){
		auto [op, l, r, k] = que[i];
		if(op == 3){	// 修改 a[l] = k
			add(l, -1, b[l]);		// 将原本删掉
			b[l] = lower_bound(li.begin(), li.end(), k)-li.begin()+1;
			add(l, 1, b[l]);
		}
		else {
			get_Root(l, r);
			if(op == 1) {		// 查区间[l, r] k的排名
				k = lower_bound(li.begin(), li.end(), k)-li.begin()+1;
				cout << _query_kRank(0, m, k)+1 << '\n';
			}
			else if(op == 2) 	// 查区间[l, r]第k名
				cout << li[_query_kth(0, m, k)-1] << '\n';
			else if(op == 4){	// 查区间[l, r] <k的最大数字
				k = lower_bound(li.begin(), li.end(), k)-li.begin()+1;
				int rk = _query_kRank(0, m, k);
				if(rk == 0)cout << "-2147483647\n";
				else {
					get_Root(l, r);
					cout << li[_query_kth(0, m, rk)-1] << '\n';
				}
			}
			else if(op == 5){		// 查区间[l, r] >k的最小数字
				k = lower_bound(li.begin(), li.end(), k)-li.begin()+1;
				int rk = _query_kRank(0, m, k+1);	// k+1能返回<=k的数量
				if(rk == r-l+1)cout << "2147483647\n";
				else {
					get_Root(l, r);
					cout << li[_query_kth(0, m, rk+1)-1] << '\n';
				}
			}
			else assert(0);
		}
	}
	return 0;
}