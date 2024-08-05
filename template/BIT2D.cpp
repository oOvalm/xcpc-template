
const int N = 3e3+10;
// op为更新树状数组的方法 可以是add, max, min...
template<class T, T(*op)(T, T)>
class BIT2D{
	T tr[N][N];
	int lowbit(int x){return x&-x;}
public:
	BIT2D(){memset(tr, 0, sizeof(tr));}
	BIT2D(T x){
		For(i,0,N-1)For(j,0,N-1)tr[i][j] = x;
	}
	void add(int x, int y, T d){	// (x, y)这点+=d
		while(x<N){
			for(int i = y; i < N; i += lowbit(i))
				tr[x][i] = op(tr[x][i], d);
			x += lowbit(x);
		}
	}
	T query(int x, int y, T ori){
		T res = ori;
		while(x){
			for(int i = y; i; i -= lowbit(i))
				res = op(res, tr[x][i]);
			x -= lowbit(x);
		}
		return res;
	}
};
ll maxx (ll x, ll y){return max(x, y);}
ll minn (ll x, ll y){return min(x, y);}
BIT2D<ll, maxx> t1;
BIT2D<ll, minn> t2(LINF);