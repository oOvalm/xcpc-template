namespace SA{
/**
 * n 长度
 * m 字符集大小
 * s 字符串	0 base
 * sa 字典序第i小是哪个后缀
 * rk 后缀i的排名
 * ht lcp(sa[i], sa[i-1]);
*/
int sa[N], rk[N], ht[N];
void getSA(string& s, int m=128)
{
	static int xx[N], yy[N], c[N];
	int *x=xx, *y=yy;
	int n = s.size();
	for(int i = 0; i < m; i++) c[i] = 0;
	for(int i = 0; i < n; i++) c[x[i] = s[i]]++;
	for(int i = 1; i < m; i++) c[i] += c[i - 1];
	for(int i = n - 1; i >= 0; i--) sa[--c[x[i]]] = i;
	for(int k = 1; k < n; k <<= 1) {
		int p = 0;
		for(int i = n - 1; i >= n - k; i--) y[p++] = i;
		for(int i = 0; i < n; i++) if(sa[i] >= k) y[p++] = sa[i] - k;
		for(int i = 0; i < m; i++) c[i] = 0;
		for(int i = 0; i < n; i++) c[x[y[i]]]++;
		for(int i = 1; i < m; i++) c[i] += c[i - 1];
		for(int i = n - 1; i >= 0; i--) sa[--c[x[y[i]]]] = y[i];
		swap(x, y);
		p = 1; x[sa[0]] = 0; y[n] = -1;
		for(int i = 1; i < n; i++) {
			if(y[sa[i - 1]] == y[sa[i]] && y[sa[i - 1] + k] == y[sa[i] + k])
				x[sa[i]] = p - 1;
			else  x[sa[i]] = p++;
		}
		if(p == n)break;
		m = p;
	}
	// 计算ht
	for(int i = 0; i < n; i++)rk[sa[i]] = i;
	int k = 0;
	for(int i = 0; i < n; i++) {
		k = max(k-1, 0ll);
		if(rk[i] == 0)continue;
		int j = sa[rk[i] - 1];
		while(s[i + k] == s[j + k]) k++;
		ht[rk[i]] = k;
	}
}
}