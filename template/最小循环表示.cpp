string getmin(string a)
{
	int n = a.size();
	a = " " + a + a;
	// cout << a << "\n";
	int i=1,j=2;
	while(j<=n){
		int k = 0;
		while(k<n && a[i+k]==a[j+k])k++;
		if(a[i+k]>a[j+k])i=i+k+1;
		else j=j+k+1;
		if(i==j)j++;
		if(i>j)swap(i,j);
	}
	string res = "";
	For(x,0,n-1){
		res += a[x+i];
	}
	return res;
}