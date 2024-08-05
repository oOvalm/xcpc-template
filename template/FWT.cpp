/**
 * op = 0 1 2 时分别为 异或 与或
 * 不需要取模就将mod去掉
*/
void FWT(int* a, int n, int op){
	for(int d = 1; d < n; d <<= 1){
		for(int m = (d<<1), i = 0; i < n; i += m){
			for(int j = 0; j < d; j++){
				int x = a[i+j], y = a[i+j+d];
				//xor: a[i+j]=x+y,a[i+j+d]=x−y;
				//and: a[i+j]=x+y;
				//or:  a[i+j+d]=x+y;
				if(op == 0) a[i+j] = (x+y)%mod, a[i+j+d] = (x-y+mod)%mod;
				else if(op == 1) a[i+j] = (x+y)%mod;
				else a[i+j+d] = (x+y)%mod;
			}
		}
	}
}
void UFWT(int* a, int n, int op){
	// 这的3个for和上面完全一样
	for(int d = 1; d < n; d <<= 1){
		for(int m = (d<<1), i = 0; i < n; i += m){
			for(int j = 0; j < d; j++){
				int x = a[i+j], y = a[i+j+d];
				//xor: a[i+j]=(x+y)/2,a[i+j+d]=(x−y)/2;
				//and: a[i+j]=x−y;
				//or:  a[i+j+d]=y−x;
				if(op == 0) a[i+j] = (x+y)%mod*inv2%mod, a[i+j+d] = (x-y+mod)%mod*inv2%mod;
				else if(op == 1) a[i+j] = (x-y+mod)%mod;
				else a[i+j+d] = (y-x+mod)%mod;
			}
		}
	}
}
// a · b -> c, a,b长度都为n 下标为[0, n-1]
void mul(int* a, int* b, int* c, int n, int op)
{
	FWT(a, n, op), FWT(b, n, op);
	For(i,0,n-1) c[i] = a[i]*b[i]%mod;
	UFWT(c, n, op);
}