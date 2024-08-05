/***** FFT 板子 *******/
namespace fft{
// 精度保证：系数总和小于 1e15 开longduoble
// 精度够(maybe总和小于1e10)可以不开longdouble
// 先调用fft_init()
// FFT_MAXN = 2^k
// fft_init() to precalc FFT_MAXN-th roots
typedef long double db;
const int FFT_MX=262144, N = 301000;
const db pi=acosl(-1.);
struct cp{
	db a,b;
	cp operator+(const cp&y)const{return (cp){a+y.a,b+y.b};}
	cp operator-(const cp&y)const{return (cp){a-y.a,b-y.b};}
	cp operator*(const cp&y)const{return (cp){a*y.a-b*y.b,a*y.b+b*y.a};}
	cp operator!()const{return (cp){a,-b};};
}nw[FFT_MX+1];
int bitrev[FFT_MX];

void dft(cp*a,int n,int flag=1){
	int d=0; 
	while((1<<d)*n!=FFT_MX) d++;
	For(i,0,n-1) if(i < (bitrev[i]>>d))
		swap(a[i], a[bitrev[i]>>d]);
	for (int l=2;l<=n;l<<=1){
		int del=FFT_MX/l*flag;
		for (int i=0;i<n;i+=l){
			cp *le=a+i, *ri=a+i+(l>>1), *w=(flag==1) ? nw : nw+FFT_MX;
			For(k,0,(l>>1)-1){
				cp ne=*ri**w;
				*ri=*le-ne, *le=*le+ne;
				le++, ri++, w+=del;
			}
		}
	}
	if(flag!=1) For(i,0,n-1)
		a[i].a/=n, a[i].b/=n;
}
void fft_init(){
	int L=0;
	while((1<<L)!=FFT_MX) L++;
	bitrev[0]=0;
	For(i,1,FFT_MX-1)
		bitrev[i] = bitrev[i>>1]>>1 | ((i&1)<<(L-1));
	nw[0]=nw[FFT_MX]=(cp){1,0};
	For(i,0,FFT_MX)nw[i] = (cp){cosl(2*pi/FFT_MX*i), sinl(2*pi/FFT_MX*i)};	//very slow
}
// n, m 分别为a, b的最高次幂, 数组a的范围为[0, n], b为[0, m]
void polymul(db*a, int n, db*b, int m, db*c) {
	static cp f[FFT_MX>>1],g[FFT_MX>>1],t[FFT_MX>>1];
	int N=2;
	while(N<=n+m)N<<=1;
	For(i,0,N-1)	// 此N非全局的N
		if(i&1){
			f[i>>1].b=(i<=n)?a[i]:0.0;
			g[i>>1].b=(i<=m)?b[i]:0.0;
		}else{
			f[i>>1].a=(i<=n)?a[i]:0.0;
			g[i>>1].a=(i<=m)?b[i]:0.0;
		}
	dft(f,N>>1); dft(g,N>>1);
	int del=FFT_MX/(N>>1);
	cp qua=(cp){0,0.25}, one=(cp){1,0}, four=(cp){4,0}, *w=nw;
	For(i,0,(N>>1)-1){
		int j=i?(N>>1)-i:0;
		t[i] = (four*!(f[j]*g[j])-(!f[j]-f[i])*(!g[j]-g[i])*(one+*w))*qua;
		w+=del;
	}
	dft(t,N>>1,-1);
	For(i,0,n+m)c[i]=(i&1) ? t[i>>1].a : t[i>>1].b;
}
}