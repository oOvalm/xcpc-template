template<class T>
struct BIT{
	T tre[N];
	BIT(){memset(tre,0,sizeof(tre));}
	int lowbit(int x){return x&-x;}
	void add(int x,T d)
	{
		while(x<N){
			tre[x] += d;
			x += lowbit(x);
		}
	}
	T query(int x)
	{
		T res = 0;
		while(x){
			res += tre[x];
			x -= lowbit(x);
		}
		return res;
	}
	T query(int l, int r){
		if(l > r) return 0;
		if(l==1)return query(r);
		return query(r)-query(l-1);
	}
};