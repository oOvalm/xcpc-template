namespace Mat{
const int N=100;
struct Matrix{
	vector<vector<int>> a;
	int n,m;
	Matrix(){a.resize(1, vector<int>(1, 0)); n = m = 1;}
	Matrix(int _n, int _m){
		n = _n, m = _m;
		a.resize(n, vector<int>(m, 0));
	}
	Matrix(Matrix const& x){
		n = x.n, m = x.m;
		a.resize(n, vector<int>(m, 0));
		for(int i = 0; i < n; i++)
			for(int j = 0; j < m; j++)
				a[i][j] = x.a[i][j];
	}

	void unit(){
		for(int i = 0; i < n; i++)
			for(int j = 0; j < m; j++){
				if(i == j)a[i][j] = 1;
				else a[i][j] = 0;
			}
	}
	Matrix operator = (int x){
		Matrix res(n, m);
		for(int i = 0; i < n; i++)
			for(int j = 0; j < m; j++)
				if(i == j)res.a[i][j] = x;
				else res.a[i][j] = 0;
		return res;
	}
	Matrix operator = (Matrix x){
		this->a = x.a;
		this->n = x.n;
		this->m = x.m;
		return *this;
	}
	Matrix operator + (Matrix x)const{
		assert(n==x.n && m==x.m);
		Matrix res(n, m);
		for(int i = 0; i < n; i++)
			for(int j = 0; j < m; j++)
				res.a[i][j] = a[i][j] + x.a[i][j];
		return res;
	}
	Matrix operator - (Matrix x)const{
		assert(n==x.n && m==x.m);
		Matrix res(n, m);
		for(int i = 0; i < n; i++)
			for(int j = 0; j < m; j++)
				res.a[i][j] = a[i][j] - x.a[i][j];
		return res;
	}
	Matrix operator * (Matrix x)const{
		assert(m==x.n);
		Matrix res(n, x.m);
		for(int i = 0; i < n; i++)
			for(int k = 0; k < m; k++) if(a[i][k])
				for(int j = 0; j < x.m; j++) if(x.a[k][j])
					res.a[i][j] += a[i][k] * x.a[k][j];
		return res;
	}
	Matrix operator += (Matrix x){(*this) = (*this)+x; return *this;}
	Matrix operator -= (Matrix x){(*this) = (*this)-x; return *this;}
	Matrix operator *= (Matrix x){(*this) = (*this)*x; return *this;}
};
}