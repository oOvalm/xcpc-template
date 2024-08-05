
#include<iostream>
#include<cstring>
#include<cmath>
#include<algorithm>
#include<vector>
#include<iomanip>
#define ll long long
#define fr first
#define se second
#define INF 0x3f3f3f3f    //1e9
#define For(i,a,b) for(int i = (a); (i) <= (b); ++i)
#define Rep(i,a,b) for(int i = (a); (i) >= (b); --i)
using namespace std;
typedef pair<int,int> pii;
const ll N = 2e2+10;
int n;
double a[N][N];
double ans[N];
void gaosi(void)
{
	int r = 1;
	For(i,1,n){//解第几个未知数
		int p = r;
		For(j,r,n)
			if(fabs(a[j][i])>fabs(a[p][i]))
				{p = j;break;}
		if(fabs(a[p][i])<1e-6)continue;
		For(j,1,n+1)
			swap(a[p][j],a[r][j]);
		Rep(j,n+1,i)
			a[r][j] /= a[r][i];
		For(j,i+1,n){
			double cha = a[j][i];
			For(k,i,n+1)
				a[j][k] -= cha*a[r][k];
		}
		// display();
		// cout << "\n";
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
void display()
{
	For(i,1,n){
		For(j,1,n+1)printf("%10f ",a[i][j]);
		printf("\n");
	}
}
int main()
{
	// ios::sync_with_stdio(false),cin.tie(0),cout.tie(0);
	// freopen("in.txt", "r", stdin);
	cout << fixed << setprecision(2);
	cin >> n;
	For(i,1,n)
		For(j,1,n+1)
			cin >> a[i][j];
	gaosi();
	return 0;
}