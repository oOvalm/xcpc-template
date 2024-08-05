namespace Mat{
const int M = 4;
struct Matrix{
	int a[M][M];
	Matrix(){memset(a,0,sizeof(a));}
	Matrix(int x){
		memset(a,0,sizeof(a));
		For(i,0,M-1)a[i][i]=x;
	}
	void unit(){
		memset(a,0,sizeof(a));
		For(i,0,M-1)a[i][i]=1;
	}
	void A(){
		
	}
	void a0(){
		
	}
	Matrix operator * (Matrix& o){
		Matrix res;
		For(i,0,M-1)For(j,0,M-1)if(a[i][j])
			For(k,0,M-1){
				res.a[i][k] += a[i][j]*o.a[j][k];
			}
		return res;
	}
	void print(){
		For(i,0,M-1){
			For(j,0,M-1)cout << a[i][j] << ' ';
			cout << '\n';
		}
	}
};
}