namespace Pollard_rho{
int maxfac;
ll ksm(ll a, ll b, ll m){
	ll res = 1;
	a %= m;
	while(b) {
		if (b & 1) res = (__int128)res * a % m;
		b >>= 1;
		a = (__int128)a * a % m;
	}
	return res;
}
bool miller_rabin(ll p)	// 判质数
{
	if (p < 2) return 0;
	if (p <= 3) return 1;
	ll d = p - 1, r = 0;
	while (!(d & 1)) ++r, d >>= 1;  // 将d处理为奇数
	for (ll k = 0; k < 10; ++k) {
		ll a = rand() % (p - 2) + 2;
		ll x = ksm(a, d, p);
		if (x == 1 || x == p - 1) continue;
		for (int i = 0; i < r - 1; ++i) {
			x = (__int128)x * x % p;
			if (x == p - 1) break;
		}
		if (x != p - 1) return 0;
	}
	return 1;
}
ll Pollard_Rho(ll x) {	// 找x的一个因子
	ll s = 0, t = 0;
	ll c = (ll)rand() % (x - 1) + 1;
	int step = 0, goal = 1;
	ll val = 1;
	for (goal = 1;; goal *= 2, s = t, val = 1){  // 倍增优化
		for (step = 1; step <= goal; ++step){
			t = ((__int128)t * t + c) % x;
			val = (__int128)val * abs(t - s) % x;
			if ((step % 127) == 0) {
				ll d = gcd(val, x);
				if (d > 1) return d;
			}
		}
		ll d = gcd(val, x);
		if (d > 1) return d;
	}
}
void findfac(ll n){
	if(n==1 || n<=maxfac)return;	// 剪枝
	if(miller_rabin(n)){	// 质数
		maxfac = max(maxfac, n);	// 找最大质因子
		return;
	}
	ll m=n;
	while(m==n)m=Pollard_Rho(n);
	while(n%m==0)n/=m;
	findfac(m), findfac(n);
}
}