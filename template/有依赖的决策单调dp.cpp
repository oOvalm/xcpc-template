int dp[N];
pii que[N];		// 维护i..n从哪转移最优
int head, tail;	// 双端队列头尾
int que_size(){return tail-head+1;}

int calc(int i, int j)	// 从dp[j]转移到dp[i]的代价
{
	db val = 1;
	For(_,1,p){	// 这循环变量别写成ij
		val *= abs(pre[i]-pre[j]-l);
	}
	return dp[j] + val;
}

tail = -1, head = 0;	// 闭区间
que[++tail] = {1, 0};	// 1... 从0转移过来
For(i,1,n){
	if(que_size() && que[head].fr < i)que[head].fr++;
	// 如果队头区间空了就pop
	if(que_size()>=2 && que[head].fr >= que[head+1].fr)head++;
	dp[i] = calc(i, que[head].se);
	// cout << i << ' ' << dp[i] << ' ' << que[head].se << '\n';
	// 如果用i更新最后一段第一个数更优，就pop调最后一段
	while(que_size() && calc(que[tail].fr, que[tail].se) >= calc(que[tail].fr, i))
		tail--;
	if(que_size()){
		// 在最后一段上二分，找出i更新更好的后缀
		int l = que[tail].fr, r = n, res = n+1;
		// 最后一段表示的区间是[l, r]
		while(l<=r){
			int mid = (l+r)>>1;
			// i更新mid更优
			if(calc(mid, i) < calc(mid, que[tail].se)){
				res = mid;
				r = mid-1;
			}
			else l = mid+1;
		}
		if(res <= n)que[++tail] = {res, i};
	}
	else que[++tail] = {i, i};
}