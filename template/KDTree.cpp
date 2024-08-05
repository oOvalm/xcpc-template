namespace KDTree1{	// P4148 简单题
struct P{
	int x, y, val;
}ps[N];
struct Node
{
	int sum, sz;	// 总和，点数
	int ls, rs, d;	// 左右孩，根据第d维划分(这里是x/y)
	int L, R, D, U;	// 这个节点管的边界
}tr[N];
int cnt, arr[N];	// 重构时用
#define lson (tr[u].ls)
#define rson (tr[u].rs)
int n, idx, root;
double alpha = 0.8;	// 有一个子树大小占>=alpha就重构
bool isbad(int u){
	return max(tr[lson].sz, tr[rson].sz) >= alpha*tr[u].sz;
}
void update(int u)	// 更新u结点的信息
{
	tr[u].sz = tr[lson].sz + tr[rson].sz + 1;	// 还有u自己
    // 开新点的时候u的编号就是这个点的下标
	tr[u].sum = tr[lson].sum + tr[rson].sum + ps[u].val;
	tr[u].L = tr[u].R = ps[u].x;
	tr[u].D = tr[u].U = ps[u].y;
	if(lson){	// 有左孩
		tr[u].L=min(tr[u].L,tr[lson].L),tr[u].R=max(tr[u].R,tr[lson].R);
		tr[u].D=min(tr[u].D,tr[lson].D),tr[u].U=max(tr[u].U,tr[lson].U);
	}
	if(rson){	// 有右孩
		tr[u].L=min(tr[u].L,tr[rson].L),tr[u].R=max(tr[u].R,tr[rson].R);
		tr[u].D=min(tr[u].D,tr[rson].D),tr[u].U=max(tr[u].U,tr[rson].U);
	}
}
bool cmpx(int a, int b){return ps[a].x < ps[b].x;}
bool cmpy(int a, int b){return ps[a].y < ps[b].y;}
void middfs(int u){
	if(!u)return;
	middfs(lson);
	arr[++cnt] = u;
	middfs(rson);
}
int build(int l, int r){
	if(l>r)return 0;
	int mid = (l+r)>>1;
	// 计算横纵坐标的方差，谁大按谁分
	double avg1 = 0, avg2 = 0, var1 = 0, var2 = 0;
	For(i,l,r){
		avg1 += ps[arr[i]].x;
		avg2 += ps[arr[i]].y;
	}
	avg1 /= (r-l+1), avg2 /= (r-l+1);
	For(i,l,r){
		var1 += (avg1 - ps[arr[i]].x)*(avg1 - ps[arr[i]].x);
		var2 += (avg2 - ps[arr[i]].y)*(avg2 - ps[arr[i]].y);
	}
	if(var1 > var2){
		nth_element(arr+l, arr+mid, arr+r+1, cmpx);
		tr[arr[mid]].d = 1;
	}
	else{
		nth_element(arr+l, arr+mid, arr+r+1, cmpy);
		tr[arr[mid]].d = 2;
	}
	int u = arr[mid];
	lson = build(l, mid-1);
	rson = build(mid+1, r);
	update(u);
	return u;
}
void rebuild(int& u){
	cnt = 0;
	middfs(u);
	u = build(1, cnt);
}
void insert(int& u, int id)
{
	if(!u){	// 当前节点位空
		u = id;
		update(u);
		return;
	}
	if(tr[u].d == 1){	// 以x为分界线的
		if(ps[id].x <= ps[u].x) insert(tr[u].ls, id);	// 在分界线右边
		else insert(tr[u].rs, id);
	}
	else{
		if(ps[id].y <= ps[u].y) insert(tr[u].ls, id);
		else insert(tr[u].rs, id);
	}
	update(u);
	if(isbad(u))rebuild(u);	// 如果当前结点过于不平衡，重构该子树
}
#define y1 alskdjalskdjaslkdj
int x1, x2, y1, y2;
int query(int u)	// 查询矩形内点的权值和
{
	if(!u || x2 < tr[u].L || x1 > tr[u].R || y1 > tr[u].U || y2 < tr[u].D)return 0;
	if(x1 <= tr[u].L && tr[u].R <= x2 && y1 <= tr[u].D && tr[u].U <= y2)return tr[u].sum;
	int res = 0;
	if(x1 <= ps[u].x && ps[u].x <= x2 && y1 <= ps[u].y && ps[u].y <= y2)
		res = ps[u].val;
	return query(lson) + query(rson) + res;
}
signed main()
{
	ios::sync_with_stdio(false),cin.tie(0),cout.tie(0);
	cin >> n;
	// 可以用类似rebuild的方法先build n个点
	int lastans = 0, op, x, y, val;
	while(cin >> op){
		if(op==1){
			cin >> x >> y >> val;
			x ^= lastans, y ^= lastans, val ^= lastans;
			ps[++idx] = {x, y, val};
			insert(root, idx);
		}
		else if(op==2){
			cin >> x1 >> y1 >> x2 >> y2;
			x1 ^= lastans, x2 ^= lastans, y1 ^= lastans, y2 ^= lastans;
			lastans = query(root);
			cout << lastans << '\n';
		}
		if(op==3)break;
	}
	return 0;
}
}


namespace KDTree2{

}