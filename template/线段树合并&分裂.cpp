void merge(int &u, int v, int l, int r)
{
	if(!v)return;	// v树为空
	if(!u){	// u树为空，当前节点直接是v
		u = v;
		return;
	}
	if(l==r){	// 合并节点信息
		tr[u].cnt += tr[v].cnt;
		tr[u].val = l;
		return;
	}
    // 暴力合并左右子树
	int mid = (l+r)>>1;
	merge(ls(u), ls(v), l, mid);
	merge(rs(u), rs(v), mid+1, r);
	update(u);
}
void split(int& u, int& v, int l, int r, int ql, int qr)
{
	if(ql <= l && r <= qr){	// 当前区间为要分裂的区间
		u = v;
		v = 0;	// v树删掉这个节点
		return;
	}
	if(!v)return;
	if(!u)u = ++tot;
	int mid = (l+r)>>1;
	if(ql <= mid)split(ls(u), ls(v), l, mid, ql, qr);
	if(qr > mid)split(rs(u), rs(v), mid+1, r, ql, qr);
	update(u);
	update(v);
}