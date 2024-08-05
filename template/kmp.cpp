// 1-base
void get_next(const string& s)
{
    nxt[1] = 0;
    nxt[2] = 1;
    int j = 1, i = 2, n = s.size()-1;
    while(i<=n){
        if(j==0 || s[i]==s[j]){
            i++; j++;
            nxt[i] = j;
        }
        else j = nxt[j];
    }
}
void kmp()
{
	string s, t;
	int n = s.size()-1, m = t.size()-1;
	get_next(s);
	For(i,1,m){
		while(j && t[i]!=s[j])j = nxt[j];
		if(t[i] == s[j])j++;
		// 出现一次
		if(j == n+1)x++, j = nxt[j];
		if(j==0)j++;
	}
}