/**
 * f[i] = min(g[j]+w(i,j))
 * w(i,j)满足四边形不等式
*/
void solve(int l, int r, int optl, int optr){
	if(l>r)return;
	int mid = (l+r)>>1;	// 算mid的opt
	int optm = optl;	// 随便定个决策点
	int rr = min(optr, mid-1);	// 根据题意来，这里表示mid的决策点一定在mid前面
	f[mid] = g[mid];		// 这应该是f[mid] = g[mid]+w(optl+1,mid)
	For(i,optl,rr){	// for在最优决策的区间
		if(f[mid] > g[i] + w(i+1, mid)){	// 如果用i更新mid更优，决策点更新为i
			f[mid] = g[i]+w(i+1, mid);
			optm = i;
		}
	}
    // 递归解决左右的f
    // 因为决策单调，左区间的决策点在[optl, optm]，右区间同理
	solve(l,mid-1,optl,optm);
	solve(mid+1,r,optm,optr);
}