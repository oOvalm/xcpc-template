#include<bits/stdc++.h>
#define ll long long
#define ull unsigned long long
#define fr first
#define se second
#define PI acos(-1)
#define INF 0x3f3f3f3f    //1e9
#define For(i,a,b) for(int i = (a); (i) <= (b); ++i)
#define Rep(i,a,b) for(int i = (a); (i) >= (b); --i)
using namespace std;
typedef pair<int,int> pii;
const ll N = 1e6+10;
const ll mod = 1e9+7;

// 欧几里得
ll gcd(ll a, ll b) {
	if (b == 0) return a;
	else return gcd(b, a % b);
}

// 快速幂
ll ksm(ll a, ll b, ll m = mod) {//记得定义mod
	ll res = 1;
	a %= m;
	while(b) {
		if (b & 1) res = res * a % m;
		b >>= 1;
		a = a * a % m;
	}
	return res;
}
template<typename T>
T ksm(T a, ll b) {
	T res(1);
	while(b) {
		if (b & 1) res = res * a;
		b >>= 1;
		a = a * a;
	}
	return res;
}

// 整数除法上下取整(含负数)
ll flor(ll x,ll y)// x/y下取整
{
	if(y==0)exit(123);
	if(x%y==0)return x/y;
	if(x>0)return x/y;
	return x/y-1;
}
ll cil(ll x,ll y)// x/y上取整
{
	if(y==0)exit(123);
	if(x%y==0)return x/y;
	if(x<0)return x/y;
	return x/y+1;
}

// 组合数
ll C(int n, int m)
{
	if(n<m)return 0;
	if (m < n - m) m = n - m;
	ll ans = 1;
	for (int i = m + 1; i <= n; i++) ans *= i;
	for (int i = 1; i <= n - m; i++) ans /= i;
	return ans;
}

// 扩展欧几里得
ll exgcd(ll a, ll b, ll &x, ll &y) {//返回gcd(a,b)
	if (b == 0) {
		x = 1;
		y = 0;
		return a;
	}
	ll d = exgcd(b, a % b, y, x);
	y -= (a / b) * x;
	return d;
}

/*
类欧几里得
求 \sum_{i=0}^n \lfloor \frac{ai+b}{c} \rfloor
*/
int floor_sum(int a, int b, int c, int n)
{
	int res = 0;
	if(a >= c)
	{
		res += n * (n + 1) * (a / c) / 2;
		a %= c;
	}
	if(b >= c)
	{
		res += (n + 1) * (b / c);
		b %= c;
	}
	int m = (a * n + b) / c;
	if(m == 0) return res;
	res += n * m - floor_sum(c, c - b - 1, a, m - 1);
	return res;
}


// 解方程 ax+by=c   最小非负 x 解
int solve_equation(ll a, ll b, ll c, ll& x, ll& y)
{
	// ax+by=gcd(a,b)
	int g = exgcd(a, b, x, y);
	if(c % g != 0){		// 无解
		x = y = 0;
		return 0;
	}
	int t = b / g, k = a / g, d = c / g;
	x *= d;		// 将方程转为 ax+by=c
	x = (x % t + t) % t;
	y = (c - a * x) / b;
	/*
	// 最小非负 y
	y *= d;
	y = (y % k + k) % k;
	x = (c - b * y) / a;
	*/
	return g;
}

// 求 a * x = b (mod m) 的解	
ll modequ(ll a, ll b, ll m) {	//等价于ax+my=b
	ll x, y;
	ll d = exgcd(a, m, x, y);	// ax+my=gcd(a,m);
	if (b % d != 0) return -1;
	m /= d; a /= d; b /= d;		//a'x+m'y=d'   x=x*d' y=y*d'
	x = x * b % m;
	if (x < 0) x += m;
	return x;
}


// 求解指数方程  a^x ≡ b (mod m) 最小x
// 复杂度为O(sqrt(m)) 用map多一个log
ll Coprime_BSGS(ll a, ll b, ll m){	// 要求: gcd(a, m)==1
	// a^(T*q) = b*a^r (mod m)
	assert(__gcd(a, m)==1);
	int T = sqrt(m)+2;
	ll v = ksm(a, T, m);
	unordered_map<ll,ll> mp;
	ll cur = 1;
	For(q,1,T){
		cur = cur*v%m;
		if(!mp.count(cur)) mp[cur] = q;
	}
	ll ans = m+1;
	cur = b;
	For(r,1,T){
		cur = cur*a%m;
		if(mp.count(cur)){
			ans = min(mp[cur]*T-r, ans);
		}
	}
	if(ans >= m)return -1;
	return ans;
}
ll BSGS(ll a, ll b, ll m){	// (a, m)可以不互质
	// a^x = b (mod m)
	// a^M * a^(x-M) = b (mod m)
	// 要让 m/gcd(a^M, m) 与 a 互质, 大概是log
	// gcd(a^M,m) = gcd(a^M%m, m)	模掉自己gcd不变
	if(b==1)return 0;
	ll cur = 1, M = 30;
	For(i,1,M){
		cur = cur*a%m;
		if(cur==b) return i;
	}
	int g = __gcd(cur, m);
	if(b % g != 0)return -1;// 无解
	cur /= g, m /= g, b /= g;
	a %= m;
	b = modequ(cur, b, m);
	ll x = Coprime_BSGS(a, b, m);
	if(x == -1)return -1;
	return x + M;
}


// 合并两个同余方程
// 初值  b=1, a=0
// x = a (mod b)
// x = c (mod d)
void merge(ll &a, ll &b, ll c, ll d) { // d <= 10^9
	/*
	x=a(mod b)  -->  x= bt+a
	-->  bt = c - a(mod d)
	-->  bt = dy+(c-a)  -->  bt-dy=c-a	(t,y为未知数)
	*/
	if (a == -1 && b == -1) return;
	ll x, y;
	ll g = exgcd(b, d, x, y);	//bx+dy=gcd(b,d)
	//bx = g(mod d)
	if ((c - a) % g != 0) { // 无解
		a = b = -1;
		return;
	}
	d /= g; // d'	//理论上b/=g 但是后面用不到b
	ll t0 = ((c - a) / g) % d * x % d;
	if (t0 < 0) t0 += d;
	// t = t0 (mod d')
	a = b * t0 + a;
	b = b * d;
}


namespace Linear_Shai{
// 线性筛
int p[N],prim[N],tot; // p[i]为i的最小质因子  prim[i]为第i个质数
void init(){
	p[1] = 1;
	for (int i = 2; i < N; i++) {
		if (!p[i]) p[i] = i, prim[++tot] = i;
		for (int j = 1; j <= tot && prim[j] * i < N; j++) {
			p[i * prim[j]] = prim[j];
			if (p[i] == prim[j]) break;
		}
	}
}

// 线性筛筛任意积性函数
// pe[i]表示i中最小值因子的指数次幂，即pe[i]=p1^a1
int p[N], prim[N], pe[N], tot, mu[N];
void init(){
	pe[1]=p[1]=1;
	For(i,2,N-1){
		if(!p[i])pe[i]=p[i]=prim[++tot]=i;
		for(int j = 1; j <= tot && prim[j]*i < N; j++){
			p[prim[j]*i] = prim[j];
			if(p[i]==prim[j]){
				pe[prim[j]*i] = pe[i]*prim[j];
				break;
			}
			pe[prim[j]*i] = prim[j];
			/*
			欧拉函数:
			if(prim[j]==p[i]){phi[i*prim[j]]=phi[i]*prim[j];break;}
			phi[i*prim[j]]=phi[i]*phi[prim[j]];
			莫比乌斯函数:
			if(prim[j]==p[i]){mu[i*prim[j]]=0;break;}
			mu[i*prim[j]]=-mu[i];
			*/
		}
	}
	// 算积性函数
	mu[1]=1;
	For(i,2,N-1)
		if(pe[i] == i) {
			if(p[i]==i) mu[i] = -1;	// 质数
			else mu[i] = 0;	// 素数幂
		}
		else mu[i] = mu[pe[i]] * mu[i/pe[i]];	// 其他
}
}

/** 取模意义下开根 **/
namespace sqrt_mod{
	/*
	满足式子 x^2 ≡ a mod p a的个数有O((p+1)/2)个

	当 p%4==3 时 x = a^{(p+1)/4}%p 为上式子的一个解
	*/
	ll W,W0;
	const ll mod = 1000000007;
	struct comp{
		ll x,y;
		comp(ll xx=0,ll yy=0){x=xx,y=yy;}
	};
	comp operator * (comp a,comp b){
		return comp((a.x*b.x%mod+a.y*b.y%mod*W%mod)%mod,(a.y*b.x%mod+a.x*b.y%mod)%mod);
	}
	inline ll ksm(ll a, ll b){
		ll t=1;a%=mod;
		while (b){
			if (b&1) t=t*a%mod;
			b>>=1;
			a=a*a%mod;
		}
		return t;
	}
	inline comp ksm(comp a, ll b){
		comp t=comp(1,0);
		while (b){
			if (b&1) t=t*a;
			b>>=1;a=a*a;
		}
		return t;
	}
	int legendre(int n) { // 无解 欧拉准则
		if(!(n % mod)) return 0;
		return ksm(n, (mod - 1) / 2);
	}
	// O(log(mod));
	inline pair<ll,ll> Cipolla(ll x){	// y^2 = x % mod
		// if (ksm(x,(mod-1)>>1)==mod-1) return make_pair(-1,-1);	// 无解 欧拉准则
		if(legendre(x)==-1)return make_pair(-1, -1);		// 与上一行等价
		srand(time(NULL));
		if(x==0)return make_pair(0, 0);
		while (1){
			W0=rand()%mod;W=(W0*W0%mod-x+mod)%mod;
			if (ksm(W,(mod-1)>>1)==mod-1) break;
		}
		comp a=ksm(comp(W0,1),(mod+1)>>1);
		return make_pair(min(a.x,mod-a.x),max(a.x,mod-a.x));	// 互为相反数
	}
}


namespace Gauss_Solve{
int a[N][N];	// 增广矩阵
int ans[N], n;
void gaosi()
{
	int r = 1;
	For(i,1,n){//解第i个未知数
		int p = r;
		For(j,r,n)	// 找到一个这一项系数不为0的式子
			if(fabs(a[j][i])>fabs(a[p][i])){
				p = j; break;
			}
		if(fabs(a[p][i])<1e-6)continue;
		For(j,1,n+1)	// 把式子换到第r行
			swap(a[p][j],a[r][j]);
		Rep(j,n+1,i)	// 将a[r][i]变成1
			a[r][j] /= a[r][i];
		// 把其他式子的第i项消掉
		For(j,i+1,n){
			double cha = a[j][i];
			For(k,i,n+1)
				a[j][k] -= cha*a[r][k];
		}
		r++;
	}
	ans[n] = a[n][n+1];
	Rep(i,n-1,1){
		double res = a[i][n+1];
		For(j,i+1,n)
			res -= ans[j]*a[i][j];
		ans[i] = res;
	}
	if(r==n+1)
		For(i,1,n){
			if(fabs(ans[i])<1e-7)ans[i]=0;//防-0.00
			cout << ans[i] << "\n";
		}
	else cout << "No Solution";
}
}


namespace Lagrange{
// n+1 个点可以唯一确定一个n次多项式
// 已知n个点(x[i], y[i])
// 求f(xx)   O(n^2)
int Lagrange(int n, int *x, int *y, int k)	
{
	int ans = k;
	For(i,1,n){
		int zi=y[i], mu=1;
		For(j,1,n) if(i!=j){
			zi = zi*(k-x[j]+mod)%mod;
			mu = mu*(x[i]-x[j]+mod)%mod;
		}
		ans = (ans + zi*ksm(mu, mod-2)%mod)%mod;
	}
	return ans;
}

// 当x取值为连续整数时， 可以优化到O(n)
int pre[N], inv[N];	// (k-x[i])的前缀积和 1/(k-x[i])
int jie[N], ijie[N];
int Lagrange_On(int n, int *x, int *y, int k)
{
	// 要保证 k > x[n] (x要满足 x[i+1] - x[i] = 1)
	int ans = 0;
	pre[0]=1;
	For(i,1,n)pre[i] = pre[i-1]*(k-x[i])%mod;
	inv[n] = ksm(pre[n], mod-2);
	Rep(i,n,1){
		inv[i-1] = inv[i]*(k-x[i])%mod;
		inv[i] = inv[i]*pre[i-1]%mod;
	}
	For(i,1,n){
		int num = ijie[i-1] * ijie[n-i]%mod * inv[i]%mod * pre[n]%mod;
		int f = ((n-i)&1?-1:1);
		ans = (ans + f*num*y[i]%mod)%mod;
	}
	return ans<0?ans+mod:ans;
}
}






































































































/***** FFT 板子 *******/
namespace fft{
// 精度保证：系数总和小于 1e15 开longduoble
// 先调用fft_init()
#define rep(i, a, b) for (int i = a; i < (int)b; i++)
// FFT_MAXN = 2^k
// fft_init() to precalc FFT_MAXN-th roots
typedef long double db;
const int FFT_MAXN=262144, N = 301000;
const db pi=acosl(-1.);
struct cp{
	db a,b;
	cp operator+(const cp&y)const{return (cp){a+y.a,b+y.b};}
	cp operator-(const cp&y)const{return (cp){a-y.a,b-y.b};}
	cp operator*(const cp&y)const{return (cp){a*y.a-b*y.b,a*y.b+b*y.a};}
	cp operator!()const{return (cp){a,-b};};
}nw[FFT_MAXN+1];int bitrev[FFT_MAXN];
void dft(cp*a,int n,int flag=1){
	int d=0;while((1<<d)*n!=FFT_MAXN)d++;
	rep(i,0,n)if(i<(bitrev[i]>>d))swap(a[i],a[bitrev[i]>>d]);
	for (int l=2;l<=n;l<<=1){
		int del=FFT_MAXN/l*flag;
		for (int i=0;i<n;i+=l){
			cp *le=a+i,*ri=a+i+(l>>1),*w=flag==1?nw:nw+FFT_MAXN;
			rep(k,0,l>>1){
				cp ne=*ri**w;
				*ri=*le-ne,*le=*le+ne;
				le++,ri++,w+=del;
			}
		}
	}
	if(flag!=1)rep(i,0,n)a[i].a/=n,a[i].b/=n;
}
void fft_init(){
	int L=0;while((1<<L)!=FFT_MAXN)L++;
	bitrev[0]=0;rep(i,1,FFT_MAXN)bitrev[i]=bitrev[i>>1]>>1|((i&1)<<(L-1));
	nw[0]=nw[FFT_MAXN]=(cp){1,0};
	rep(i,0,FFT_MAXN+1)nw[i]=(cp){cosl(2*pi/FFT_MAXN*i),sinl(2*pi/FFT_MAXN*i)};	//very slow
}

void convo(db*a,int n,db*b,int m,db*c){
	static cp f[FFT_MAXN>>1],g[FFT_MAXN>>1],t[FFT_MAXN>>1];
	int N=2;while(N<=n+m)N<<=1;
	rep(i,0,N)
		if(i&1){
			f[i>>1].b=(i<=n)?a[i]:0.0;
			g[i>>1].b=(i<=m)?b[i]:0.0;
		}else{
			f[i>>1].a=(i<=n)?a[i]:0.0;
			g[i>>1].a=(i<=m)?b[i]:0.0;
		}
	dft(f,N>>1);dft(g,N>>1);
	int del=FFT_MAXN/(N>>1);
	cp qua=(cp){0,0.25},one=(cp){1,0},four=(cp){4,0},*w=nw;
	rep(i,0,N>>1){
		int j=i?(N>>1)-i:0;
		t[i]=(four*!(f[j]*g[j])-(!f[j]-f[i])*(!g[j]-g[i])*(one+*w))*qua;
		w+=del;
	}
	dft(t,N>>1,-1);
	rep(i,0,n+m+1)c[i]=(i&1)?t[i>>1].a:t[i>>1].b;
}
}


namespace ntt{
/*******NTT 板子******/
int g=3;//mod的原根
int gi=ksm(g,mod-2); //g 的逆元
int rev[N];
ll a[2*N],b[2*N];
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
	ntt(a,tot,1);
	ntt(b,tot,1);
	for(int i = 0; i <= tot; i++)
		a[i] = a[i]*b[i]%mod;
	ntt(a,tot,-1);
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

/*
多项式	NTT
*/
namespace vector_ntt{
#define poly vector<ll>
int lim,rev[N],G=3,Ginv=ksm(G,mod-2);
void polyinit(int n)
{
    for (lim = 1; lim < n; lim <<= 1);
    for (int i = 0; i < lim; i ++) rev[i] = (rev[i >> 1] >> 1) | (i & 1 ? lim >> 1 : 0);
}
void NTT(poly &f, int op)
{
    for (int i = 0; i < lim; i ++)
    {
        if (i < rev[i]) swap(f[i], f[rev[i]]);
    }
    for (int mid = 1; mid < lim; mid <<= 1)
    {
        int Gn = ksm(op == 1 ? G : Ginv, (mod - 1) / (mid << 1));
        for (int i = 0; i < lim; i += mid * 2)
        {
            for (int j = 0, G0 = 1; j < mid; j ++, G0 = G0 * Gn % mod)
            {
                int x = f[i + j], y = G0 * f[i + j + mid] % mod;
                f[i + j] = (x + y) % mod, f[i + j + mid] = (x - y + mod) % mod;
            }
        }
    }
    if (op == -1)
    {
        int inv = ksm(lim, mod - 2);
        for (int i = 0; i < lim; i ++) f[i] = f[i] * inv % mod;
    }
}
poly operator + (poly f, int x)
{
    f[0] = (f[0] + x) % mod;
    return f;
}
poly operator - (poly f, int x)
{
    f[0] = (f[0] - x + mod) % mod;
    return f;
}
poly operator * (poly f, int x)
{
    for (int i = 0; i < (int)f.size(); i ++) f[i] = f[i] * x % mod;
    return f;
}
poly operator + (poly f, poly g)
{
    poly res = f;
    res.resize(max((int)f.size(), (int)g.size()));
    for (int i = 0; i < (int)g.size(); i ++) res[i] = (res[i] + g[i]) % mod;
    return res;
}
poly operator - (poly f, poly g)
{
    poly res = f;
    res.resize(max((int)f.size(), (int)g.size()));
    for (int i = 0; i < (int)g.size(); i ++) res[i] = (res[i] - g[i] + mod) % mod;
    return res;
}
poly operator * (poly f, poly g)
{
    int n = f.size() + g.size() - 1;
    polyinit(n), f.resize(lim), g.resize(lim);
    NTT(f, 1), NTT(g, 1);
    for (int i = 0; i < lim; i ++) f[i] = f[i] * g[i] % mod;
    NTT(f, -1), f.resize(n);
    return f;
}
}




