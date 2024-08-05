/*
P3834 【模板】可持久化线段树 2
主席树求 静态区间第k小
*/
#include<bits/stdc++.h>
#define For(i,a,b) for(int i = a; i <= b; i++)
using namespace std;
const int N = 2e7+10;
struct Node{
	int ls, rs, num;
}tr[N];
int rt[N], cnt;
void pushup(int u) {
	tr[u].num = tr[tr[u].ls].num+tr[tr[u].rs].num;
}
void build(int &u, int l, int r)
{
	u = ++cnt;
	if(l==r)return;
	int mid = (l+r)>>1;
	build(tr[u].ls, l, mid);
	build(tr[u].rs, mid+1, r);
	pushup(u);
}
void add(int &u, int pre, int l, int r, int pos, int x)
{
	u = ++cnt;		// 新建节点
	tr[u] = tr[pre];	// 继承原节点的信息
	if(l==r){
		assert(pos == l);
		tr[u].num+=x;
		return;
	}
	int mid = (l+r)>>1;
	// 递归下含有pos的那个子树
	if(pos <= mid)add(tr[u].ls, tr[pre].ls, l, mid, pos, x);
	else add(tr[u].rs, tr[pre].rs, mid+1, r, pos, x);
	pushup(u);
}
int n, a[N], q;
int query(int u, int v, int l, int r, int k)
{
	// cout << u << ' ' << v << ' ' << l << ' ' << r << ' ' << tr[u].num << ' ' << tr[v].num << '\n';
	if(l==r)return l;
	int mid = (l+r)>>1;
	int lsnum = tr[tr[u].ls].num-tr[tr[v].ls].num;
	if(lsnum >= k) return query(tr[u].ls, tr[v].ls, l, mid, k);
	else return query(tr[u].rs, tr[v].rs, mid+1, r, k-lsnum);
}
signed main()
{
	cin >> n >> q;
	vector<int> li;
	For(i,1,n){
		cin >> a[i];
		li.push_back(a[i]);
	}
	sort(li.begin(), li.end());	// 离散化
	li.erase(unique(li.begin(), li.end()), li.end());
	int m = li.size();
	build(rt[0], 1, m);		// 原数组的线段树
	For(i,1,n){
		int pos = lower_bound(li.begin(), li.end(), a[i])-li.begin()+1;
		add(rt[i], rt[i-1], 1, m, pos, 1);	// 插一个数一个版本
	}
	while(q--){
		int l, r, k;
		cin >> l >> r >> k;
		// 将r的树和l-1的树做差
		cout << li[query(rt[r], rt[l-1], 1, m, k)-1] << '\n';
	}
}