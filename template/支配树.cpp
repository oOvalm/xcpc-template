/*
给定一个有向图，给定起点S，终点T，要求在所有从S到T的路径中，必须经过的点有哪些(即各条路径上点集的交集），称为支配点。
简而言之，如果删去某个点（即支配点），则不存在一条路径能够使S到达T。由支配点构成的树叫做支配树
支配树是一颗由起点S为根的树，根到树上任意一点路径经过的点都是它的支配点。
对于每颗子树，都符合上述性质（即子树内一点到该子树的根路径经过的点都是该子树的根的支配点）
使用add_edge(u, v)加有向边u->v
solve(root)，以root为根求支配树
idom[i] -- i 为支配树上的边
*/
class DominatorTree{
public:
    std::vector <std::vector <int>> e, _e, tmp;	// e原图 _e反图
    std::vector <int> dfn, inv;	// dfn[x] x的dfs序
    int dfncnt;
    std::vector <int> sdom, idom;	// sdom[x] x的半支配点， idom[x] x的最近支配点
    std::vector <int> fa, father, value;

    explicit DominatorTree(int n):dfncnt(0){
        int sz = n + 10;
        e.resize(sz);
        _e.resize(sz);
        tmp.resize(sz);
        dfn.resize(sz);
        inv.resize(sz);
        sdom.resize(sz);
        idom.resize(sz);
        fa.resize(sz);
        father.resize(sz);
        value.resize(sz);
    }

    void add_edge(int u, int v){
        e[u].emplace_back(v);
        _e[v].emplace_back(u);
    }

    int min(int u, int v){
        return dfn[u] < dfn[v] ? u : v;
    }

    int find(int u){
        if (fa[u] == u) return u;
        int f = fa[u];
        fa[u] = find(fa[u]);
        if (dfn[sdom[value[f]]] < dfn[sdom[value[u]]]) value[u] = value[f];
        return fa[u];
    }

    void dfs(int u){
        dfn[u] = ++ dfncnt;
        inv[dfncnt] = u;
        for (auto v : e[u]){
            if (dfn[v]){
                continue;
            }
            father[v] = u;
            dfs(v);
        }
    }

    void solve(int rt){
        dfs(rt);
        for (int i = 1; i <= dfncnt; ++ i){
            fa[inv[i]] = value[inv[i]] = sdom[inv[i]] = inv[i];
        }
        for (int i = dfncnt; i >= 2; -- i){
            int u = inv[i];
            for (auto v : _e[u]){
                if (!dfn[v]) continue;
                if (dfn[v] < i){
                    sdom[u] = min(sdom[u], v);
                }
                else{
                    find(v);
                    sdom[u] = min(sdom[u], sdom[value[v]]);
                }
            }
            fa[u] = father[u];
            tmp[sdom[u]].emplace_back(u);
            int pa = fa[u];
            for (auto v : tmp[pa]){
                find(v);
                idom[v] = value[v];
            }
            tmp[pa].clear();
        }
        for (int i = 2; i <= dfncnt; ++ i){
            int u = inv[i];
            idom[u] = sdom[u] == sdom[idom[u]] ? sdom[u] : idom[idom[u]];
        }
    }
};