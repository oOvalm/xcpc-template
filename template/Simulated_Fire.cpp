
/*
模拟退火
可调
初始温度 T 
结束温度 eps
降温系数 一般为接近1的数字
参考：
1e-1精度：下界1e-3或1e-4 降温系数0.8
1e-5精度：下界1e-7或1e-8 降温系数0.99
*/
namespace Simulated_Annealing
{
struct Point
{
	double x, y, z;
	Point(){}
	Point(double _x, double _y, double _z):x(_x),y(_y),z(_z){}
	bool operator == (Point b){return abs(x-b.x)<1e-6 && abs(y-b.y)<1e-6 && abs(z-b.z)<1e-6;}
};
double dis(Point a, Point b)
{
	double dx = a.x-b.x, dy = a.y-b.y, dz = a.z-b.z;
	return sqrt(dx*dx+dy*dy+dz*dz);
}
Point p[N];
int n;
double maxdis(Point a)
{
	double res = 0;
	for(int i = 1; i <= n; i++){
		res = max(res, dis(a, p[i]));
	}
	return res;
}
bool out(Point a)
{
	return abs(a.x) > 1e5 || abs(a.y) > 1e5 || abs(a.z) > 1e5;
}
Point ansp;
double ans;
void SA()	// 模拟退火
{
	double T = 550, eps = 1e-10;
	while(T>eps){	// 结束温度
		Point tmp(ansp);
		// 随机到的点在答案范围内 (也可以不加，温度高容易随机到边界外 可能会卡死)
		while(tmp == ansp || out(tmp)){		
			tmp.x = ansp.x + T*(2*rand()-RAND_MAX);		// 根据当前温度随机选
			tmp.y = ansp.y + T*(2*rand()-RAND_MAX);
			tmp.z = ansp.z + T*(2*rand()-RAND_MAX);
		}
		double tmpans = maxdis(tmp);	// 计算答案
		if(tmpans - ans < 1e-6){		// 比当前答案更优 更新答案
			ans = tmpans;
			ansp=tmp;
		}
		// 比当前答案不优 有一定概率选择不优的答案
		else if(exp(-abs(tmpans-ans)/T) * RAND_MAX > rand()){	
			ans = tmpans;
			ansp=tmp;
		}
		T *= 0.998;	//降温系数
	}
}
void tuihuo()
{
	ansp = Point(0,0,0);	// 随便指定一个答案
	ans = maxdis(ansp);
	// 做多几次退火  
	while(clock() < 1000)SA();		// 卡时间
	// for(int i = 0; i < 100; i++)SA();	// 固定次数
}
};



