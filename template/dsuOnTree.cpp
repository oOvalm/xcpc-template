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