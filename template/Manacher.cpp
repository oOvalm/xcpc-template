int p[N];// 回文半径 [i-p[i]+1,i+p[i]-1]是回文串
int manacher(string& t)
{
	string s=" #";//每个字符间插个原串不存在的字符 头尾插一个
	for(int i = 0; i < t.size(); i++){
		s += t[i];
		s += "#";
	}
	int n=s.size();
	n--;
	For(i,1,n)p[i]=0;
	int m=0,r=0;
	for(int i = 1; i <= n; i++){
		if(i>r)p[i]=1;
		else p[i] = min(p[2*m-i], r-i+1); // i+p[i]-1 <= r
		while(i-p[i]>0 && i+p[i]<=n && s[i-p[i]]==s[i+p[i]]){//暴力扩
			p[i]++;
		}
		if(i+p[i]-1>r){//更新右端点和中点
			m=i;r=i+p[i]-1;
		}
	}
	int res = 0;//最大回文半径 对应插入后回文串长度为2*res-1  对应原串 res-1
	For(i,1,n){
		res = max(res,p[i]);
	}
	return res-1;
}
