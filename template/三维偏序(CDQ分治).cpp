/*
$n$个点，每个点有一个坐标$(a_i, b_i, c_i)$ (可能会有重点), 
对于每个点$i$求有多少个$j$满足
$a_j\leq a_i$, $b_j \leq b_i$, $c_j \leq c_i$, $j \neq i$。
*/
#include<bits/stdc++.h>
#define ll long long
#define int ll
#define For(i,a,b) for(int i = a; i <= b; ++i)
using namespace std;
const int N = 2e5+10;
BIT<int> T;
struct Node{
	int x, y, z, cnt, ans; // cnt:重合点数, ans:这个点的答案
	bool operator < (const Node& b)const{
		return array<int, 3>{x, y, z} < array<int, 3>{b.x, b.y, b.z};
	}
	bool operator == (const Node& b){
		return !((*this)<b)&&!(b<(*this));
	}
}a[N], p[N];
bool cmpy(const Node& a, const Node& b){
	return array<int, 2>{a.y, a.z} < array<int, 2>{b.y, b.z};
}
int n, k, ans[N];
void cdq(int l, int r)
{
	if(l==r)return;
	int mid = (l+r)>>1, pos = l;
	cdq(l, mid);cdq(mid+1, r);
	// 合并两边(第二维)
	sort(p+l, p+mid+1, cmpy);
	sort(p+mid+1, p+r+1, cmpy);
	// 树状数组做第三维
	For(i, mid+1, r){
		// 把可以的第二维都进树状数组
		while(pos <= mid && p[pos].y <= p[i].y){
			T.add(p[pos].z, p[pos].cnt);
			pos++;
		}
		p[i].ans += T.query(p[i].z);
	}
	For(i,l,pos-1)T.add(p[i].z, -p[i].cnt);
}
signed main()
{
	cin >> n >> k;
	int idx = 0;
	For(i,1,n) cin >> a[i].x>>a[i].y>>a[i].z;
	sort(a+1,a+1+n);	// 第一维排序
	p[++idx] = a[1];
	p[idx].cnt = 1;
	For(i,2,n){
		if(p[idx] == a[i])p[idx].cnt++;
		else{
			p[++idx] = a[i];
			p[idx].cnt = 1;
		}
	}
	cdq(1, idx);
	For(i,1,idx){
		ans[p[i].ans+p[i].cnt-1] += p[i].cnt;
	}
	For(i,0,n-1)cout << ans[i] << '\n';
}