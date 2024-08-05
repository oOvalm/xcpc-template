namespace VirtualTree{
int k = 0, m = 0, tot = 0;
cin >> k;
For(i,1,k){	// k个关键点
	cin >> b[i];
	imp[c] = 1;
}
auto cmp = [&](int a, int b){return dfn[a]<dfn[b];}
m = tot = k;
sort(b+1, b+1+tot, cmp);		// dfs序排序
For(i,2,m) b[++tot] = LCA(b[i], b[i-1]);	// 相邻关键点的lca
b[++tot] = 1;	// 加根节点
sort(b+1, b+1+tot, cmp);
tot = unique(b+1, b+1+tot)-b-1;	// 去重
For(i,1,tot-1){
	int lc = LCA(b[i], b[i+1]);
	H[lc].push_back(b[i+1]);	// 建树 上到下的边
}
/***work***/
// 清空
For(i,1,tot){
	H[b[i]].clear();
	imp[b[i]] = 0;
}
}