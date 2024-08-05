
namespace ntt{
/*******NTT 板子******/
const int mod = 998244353;
int g=3;//mod的原根
int gi=ksm(g,mod-2); //g 的逆元
int rev[N];
void ntt(ll A[], int tot, int type)
{
	for(int i = 0; i < tot; i++)
		if(i<rev[i]) swap(A[i], A[rev[i]]);
	for(int mid = 1; mid < tot; mid <<= 1){
		ll gn = ksm(type>0 ? g : gi, (mod-1)/(mid<<1));//与下面ifelse等价 可整除
		// if(type==1)gn = ksm(g, (mod-1)/(mid<<1));
		// else gn = ksm(gi, (mod-1)/(mod<<1));
		int R = mid<<1;
		for(int i = 0; i < tot; i+=R){
			ll g0 = 1;
			for(int j = i; j < mid+i; j++, g0=g0*gn%mod){
				ll x = A[j], y=g0*A[j+mid]%mod;
				A[j] = (x+y)%mod;
				A[j+mid] = (x-y+mod)%mod;
			}
		}
	}
	ll invtot = ksm(tot, mod-2);
	if(type==-1){
		for(int i = 0; i < tot; i++)
			A[i] = A[i]*invtot%mod;
	}
}
void polymul_ntt(ll a[], ll b[],int n,int m)
{
	int tot=1;
	int l=0;	//求rev
	while(tot <= n+m){
		tot<<=1; l++;
	}
	for(int i = 0; i <= tot; i++)
		rev[i] = (rev[i>>1]>>1) | ((i&1)<<(l-1));
	// 清空后面的系数
	For(i,n+1,tot)a[i] = 0;
	For(i,m+1,tot)b[i] = 0;
	ntt(a,tot,1);
	ntt(b,tot,1);
	for(int i = 0; i <= tot; i++)
		a[i] = a[i]*b[i]%mod;
	ntt(a,tot,-1);
}
vector<int> mul_all(int l, int r, vector<vector<int>>& v)
{
	if(l==r){
		return v[l];
	}
	static int a[N], b[N];
	int mid = (l+r)>>1;
	auto lp = mul_all(l, mid, v), rp = mul_all(mid+1, r, v);
	int n = lp.size()-1, m = rp.size()-1;
	For(i,0,n)a[i] = lp[i];
	For(i,0,m)b[i] = rp[i];
	polymul_ntt(a, b, n, m);
	vector<int> res(n+m+1);
	For(i,0,n+m)res[i] = a[i];
	return res;
}
//*******多项式求逆*******
/**maybe问题一****/
// 设 g0 满足  g0*f ≡ 1 (mod x^(n/2))  n/2上取整
// 1移项 --> 平方 模数变为 (mod x^n) --> 乘个g --> 已知 f*g=1
// g-2g0+f*g0^2 ≡ 0 (mod x^n)
// g ≡ g0*(2-f*g0) (mod x^n)
// 最高次项为len-1 的多项式f  返回多项式的逆g  倍增法求解  需要2个辅助数组c d
ll c[N],d[N];
void polyinv(int len,ll f[], ll g[]) 
{
	if(len==1){
		g[0]=ksm(f[0],mod-2);
		return;
	}
	polyinv((len+1)>>1,f,g);
	int tot=1,l=0;
	while(tot<(len<<1)){
		tot<<=1;l++;
	}
	For(i,0,tot-1)rev[i] = (rev[i>>1]>>1) | ((i&1)<<(l-1));
	For(i,0,tot-1){
		c[i]=f[i];
		if(i<((len+1)>>1))d[i]=g[i];
		else d[i]=0;
	}
	ntt(c,1,tot);
	ntt(d,1,tot);
	For(i,0,tot-1){
		c[i] = (d[i]*((2-c[i]*d[i]%mod+mod)%mod)%mod)%mod;
	}
	ntt(c,-1,tot);	//ntt里已经除了tot
	for(int i = 0; i < len; i++)
		g[i] = c[i];
}
}
