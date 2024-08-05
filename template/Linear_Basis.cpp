namespace Linear_Basis{
const int B = 63;	// 最高位数
struct linear_basis
{
	ll num[B+1];
	bool insert(ll x)	// 将x插进线性基 true插入成功
	{
		for(int i = B; i >= 0; i--) if(x & (1ll<<i)) {
			if(num[i] == 0){
				num[i] = x;
				return true;
			}
			x ^= num[i];
		}
		return false;
	}
	ll querymin(ll x)	// 如果x为0 线性基可以表示x
	{
		for(int i = B; i >= 0; i--) {
			x = min(x, x ^ num[i]);
		}
		return x;
	}
};
}