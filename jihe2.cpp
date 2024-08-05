#include <bits/stdc++.h>
#define ll long long
#define fr first
#define se second
#define For(i,a,b) for(int i = (a); i <= (b); ++i)
using namespace std;
typedef pair<int, int> pii;
const ll N = 5010;
/*
# 任意多边形面积
$$ S=\frac{|\sum_{i=1}^n x_iy_{i+1}-x_{i+1}y_{i}|}{2} $$
若多边形的点逆时针顺序为$p_0,p_1,p_2,\dots p_{n-1}$则有 ($\times$为叉积)
$$ S=\frac{1}{2}\sum_{i=0}^{n-1}p_i \times p_{(i+1) \mod n} $$
*/
namespace geometry2{
#define db double
#define EPS 1e-9
#define PI acos(-1)
// 判断a的符号 返回值为 -1,0,1
inline int sign(db a) { return a < -EPS ? -1 : a > EPS; }
// 比较 a 和 b 的大小(带eps) a<b -1; a==b 0; a>b 1
inline int cmp(db a, db b) { return sign(a - b); }
// 点/向量
struct P{
	db x, y;
	int id;
	P(db _x = 0, db _y = 0) : x(_x), y(_y) {}

	P operator+(P p) const{ return P(x + p.x, y + p.y); }
	P operator-(P p) const{ return P(x - p.x, y - p.y); }
	P operator*(db a) { return P(x * a, y * a); }
	P operator/(db a) { return P(x / a, y / a); }
	// 按照pair(x,y)排序
	bool operator<(P p) const{
		int op = cmp(x, p.x);
		if (op != 0)return op == -1;
		return cmp(y, p.y) == -1;
	}
	// 两点重叠(带eps)
	bool operator==(P p) const { return cmp(x, p.x) == 0 && cmp(y, p.y) == 0; }
	bool operator!=(P p) const { return !((*this) == p);}
	// 点乘 (*this) · p
	db dotmul(P p) { return x * p.x + y * p.y; }
	// 叉乘	(*this) × p
	db xmul(P p) { return x * p.y - y * p.x; }
	// 模长平方
	db abs2() { return x * x + y * y; }
	double abs() { return sqrt(abs2()); }
	// 点到点p的距离
	double disTo(P p) { return (*this - p).abs(); }
	double disTo2(P p) { return (*this - p).abs2(); }
	// 极角 [-PI,PI]
	double alpha() { return atan2(y, x); }
	// 逆时针旋转90度
	P rota90() { return P(-y, x); }
	// 单位向量
	P unit() { return (*this) / abs(); }
	// 点是否在x轴上方 or 在x正半轴 (极角排序) 在上方返回1
	int quad() const { return sign(y) == 1 || (sign(y) == 0 && sign(x) >= 0); }
	// 逆时针旋转 a (弧度)
	P rota(db a) { return P(x * cos(a) - y * sin(a), x * sin(a) + y * cos(a)); }
	void read(){cin >> x >> y;}
	void print(){cout << "(" << x << ", " << y << ")\n";}
};
// 向量(p1p2) · 向量(p1p3)
db dotmul(P p1, P p2, P p3) {
	return (p2.x - p1.x) * (p3.x - p1.x) + (p2.y - p1.y) * (p3.y - p1.y);
}
// 向量(p1p2) × 向量(p1p3)
db xmul(P p1, P p2, P p3) {
	return (p2.x - p1.x) * (p3.y - p1.y) - (p3.x - p1.x) * (p2.y - p1.y);
}
// 叉乘的符号
int xmulOP(P p1, P p2, P p3) { return sign(xmul(p1, p2, p3)); }


// 直线p1p2 与 直线q1q2 是否相交
bool LLIsCross(P p1, P p2, P q1, P q2){
	db a1 = xmul(q1, q2, p1), a2 = -xmul(q1, q2, p2);
	return sign(a1 + a2);
}
// 直线p1p2 与 直线q1q2 的交点
P LLCross(P p1, P p2, P q1, P q2){
	db a1 = xmul(q1, q2, p1), a2 = -xmul(q1, q2, p2);
	return (p1 * a2 + p2 * a1) / (a1 + a2);
}
// 区间[l1,r1], [l2,r2]是否相交 (相交返回1)
bool intersect(db l1, db r1, db l2, db r2){
	if (l1 > r1) swap(l1, r1);
	if (l2 > r2) swap(l2, r2);
	return !(cmp(r1, l2) == -1 || cmp(r2, l1) == -1);
}

// 线段 p1p2, q1q2 是否相交
bool isSegCross(P p1, P p2, P q1, P q2){
	if (intersect(p1.x, p2.x, q1.x, q2.x) && intersect(p1.y, p2.y, q1.y, q2.y) 
		&& xmulOP(p1, p2, q1) * xmulOP(p1, p2, q2) <= 0 
		&& xmulOP(q1, q2, p1) * xmulOP(q1, q2, p2) <= 0)return true;
	return false;
}
// 线段 p1p2, q1q2 是否严格相交 交点不在端点
bool isSegCross_strict(P p1, P p2, P q1, P q2){
	if (xmulOP(p1, p2, q1) * xmulOP(p1, p2, q2) < 0 
		&& xmulOP(q1, q2, p1) * xmulOP(q1, q2, p2) < 0)return true;
	return false;
}
// m 是否在[a,b]之间
bool isInside(db a, db b, db m){
	if (a > b) swap(a, b);
	return cmp(a, m) <= 0 && cmp(m, b) <= 0;
}
// m 是否在a b为对角线的矩形区域内
bool isInside(P a, P b, P m){
	return isInside(a.x, b.x, m.x) && isInside(a.y, b.y, m.y);
}
// m 是否在线段ab上
bool isonSeg(P a, P b, P m){
	return xmulOP(a, b, m) == 0 && isInside(a, b, m);
}
// m是否严格在线段ab上
bool isonSeg_strict(P a, P b, P m){
	return xmulOP(a, b, m) && sign((m - a).dotmul(a - b)) * sign((m - b).dotmul(a - b)) < 0;
}

// 点 q 到直线 p1p2 的投影 p1 != p2
P proj(P p1, P p2, P q){
	P dir = p2 - p1;
	return p1 + dir * (dir.dotmul(q - p1) / dir.abs2());
}
// 点 q 关于直线 p1p2 的轴对称点
P reflect(P p1, P p2, P q){
	return proj(p1, p2, q) * 2 - q;
}
// q 到线段 p1p2 的最小距离
db nearest(P p1, P p2, P q){
	if (p1 == p2) return p1.disTo(q);
	P h = proj(p1, p2, q);
	if (isInside(p1, p2, h)) return q.disTo(h);
	return min(p1.disTo(q), p2.disTo(q));
}
// 线段 p1p2 q1q2 的最小距离
db Segdis(P p1, P p2, P q1, P q2){
	if (isSegCross(p1, p2, q1, q2))return 0;
	return min(min(nearest(p1, p2, q1), nearest(p1, p2, q2)),
			   min(nearest(q1, q2, p1), nearest(q1, q2, p2)));
}
//极角排序
// 此cmp1对同极角的不作区分
// 要加const的话xmul也要加
bool cmp1(P a, P b){
	// 不在x轴同边 x轴上方和x正半轴排前
	if (a.quad() != b.quad())return a.quad() > b.quad(); 
	return sign(a.xmul(b)) > 0;		// b要在a的逆时针方向
};

/**********************多边形**********************/
/**多边形端点都需要先极角排序**/
// 求多边形面积(简单多边形)
db area(vector<P> &ps){
	db res = 0;
	for (int i = 0; i < ps.size(); i++)
		res += ps[i].xmul(ps[(i + 1) % ps.size()]);
	return res / 2;
}
// 点包含 判断点p是否在多边形内 (2:内 1:在边上 0:外)(简单多边形)
int contain(const vector<P> &ps, P p){
	int n = ps.size(), res = 0;
	for (int i = 0; i < n; i++){
		P u = ps[i], v = ps[(i + 1) % n];
		if (isonSeg(u, v, p))return 1;
		if (cmp(u.y, v.y) <= 0)swap(u, v);
		if (cmp(p.y, u.y) > 0 || cmp(p.y, v.y) <= 0)continue;
		res ^= (xmulOP(p, u, v) > 0); //在p左边的有多少 p向左射线与多边形多少条边相交
	}
	return res * 2;
}
// 求凸包 给定点击 求最小的多边形使得所有点在多边形上/内
vector<P> convexHull(vector<P> ps) { // 只要顶点

	int n = ps.size();
	if (n <= 1)return ps;
	sort(ps.begin(), ps.end()); // 按(x,y)字典序排 先排x再排y Point的小于号
	vector<P> res(n * 2);
	int k = 0;
	for (int i = 0; i < n; res[k++] = ps[i++]) //求下凸包
		while (k > 1 && xmulOP(res[k - 2], res[k - 1], ps[i]) <= 0)k--;
	for (int i = n - 2, t = k; i >= 0; res[k++] = ps[i--]) // 求上凸包
		while (k > t && xmulOP(res[k - 2], res[k - 1], ps[i]) <= 0)k--;
	res.resize(k - 1); // ps[0]被加入了2次(去掉一个)
	return res;
}
// 要顶点和边上的点
vector<P> convexHull_notStrict(vector<P> ps){
	// ps要去重
	int n = ps.size();
	if (n <= 1)return ps;
	sort(ps.begin(), ps.end()); // 按(x,y)字典序排 先排x再排y Point的小于号
	vector<P> res(n * 2);
	int k = 0;
	for (int i = 0; i < n; res[k++] = ps[i++]) //求下凸包
		while (k > 1 && xmulOP(res[k - 2], res[k - 1], ps[i]) < 0) k--;
	for (int i = n - 2, t = k; i >= 0; res[k++] = ps[i--]) // 求上凸包
		while (k > t && xmulOP(res[k - 2], res[k - 1], ps[i]) < 0) k--;
	res.resize(k - 1); // ps[0]被加入了2次(去掉一个)
	return res;
}
// 切凸多边形(极角序) 用有向直线q1q2切凸多边形 返回切完后q1q2左边(逆时针方向)的多边形
vector<P> convexCut(const vector<P> &ps, P q1, P q2){
	vector<P> res;
	int n = ps.size();
	for (int i = 0; i < n; i++){
		P p1 = ps[i], p2 = ps[(i + 1) % n];
		int d1 = xmulOP(q1, q2, p1), d2 = xmulOP(q1, q2, p2);
		if (d1 >= 0)res.push_back(p1); // p1在直线左边
		if (d1 * d2 < 0)
			res.push_back(LLCross(p1, p2, q1, q2)); // p1p2 穿过直线 交点为新多边形的一个顶点
	}
	return res;
}
// 旋转卡壳(极角序) 凸多边形最远两点距离 / 凸多边形卡在两平行直线间 转一圈后 平行直线最大的距离
db convexDiameter(vector<P> ps){
	int n = ps.size();
	if (n <= 1)return 0;
	int is = 0, js = 0;
	for (int k = 1; k < n; k++){
		is = (ps[k] < ps[is] ? k : is); //选最大的
		js = (ps[js] < ps[k] ? k : js); //选最小的
	}
	int i = is, j = js;
	db res = ps[i].disTo(ps[j]);
	do{
		if ((ps[(i + 1) % n] - ps[i]).xmul(ps[(j + 1) % n] - ps[j]) >= 0)
			j = (j + 1) % n;
		else i = (i + 1) % n;
		res = max(res, ps[i].disTo(ps[j]));
	} while (i != is || j != js);
	return res;
}
// 最近点对 ps先按x后按y排序
// min_dist(ps, 0, ps.size()-1)
db min_dist(vector<P>& ps, int l, int r)
{
	if(r-l<=5){
		db res = 1e100;
		For(i,l,r)For(j,i+1,r)
			res = min(res, ps[i].disTo(ps[j]));
		return res;
	}
	int m = (l+r)>>1;
	db res = min(min_dist(ps, l, m), min_dist(ps, m+1, r));
	vector<P> qs;
	For(i,l,r) if(abs(ps[i].x-ps[m].x) <= res)
		qs.push_back(ps[i]);
	sort(qs.begin(), qs.end(), [](P a, P b){return a.y<b.y;});
	For(i,1,(ll)qs.size()-1){
		for(int j = i-1; j>=0 && qs[j].y>=qs[i].y-res;j--)
			res = min(res, qs[i].disTo(qs[j]));
	}
	return res;
}


// 半平面交
struct Line{
	P a,b;
	Line(P _a=P(0,0), P _b=P(0,0)):a(_a),b(_b){}
	// 倾斜角(有向线段)
	double angle(){
		P v=b-a;
		return atan2(v.y,v.x);
	}
	//直线斜率排序(有向)
	bool operator < (const Line L)const{
		P va=b-a,vb=L.b-L.a;
		double A=atan2(va.y,va.x), B=atan2(vb.y,vb.x);
		if(cmp(A, B)==0)return va.xmul(L.b-a)>=0;
		return cmp(A,B)<0;
	}
};
// L1,L2 交点是否在 L0 右边
bool onRight(Line L0,P p){
	return xmulOP(L0.a,L0.b,p)<0;//三点呈顺时针
}
// 传入多边形的边 返回面积 面积为0点集不为空
double HalfPlaneInter(vector<Line> ls){
	sort(ls.begin(),ls.end());
	int cnt=0,n=ls.size();
	for(int i = 0; i < n-1; i++){// 斜率相同的直线要靠左的一条
		if(cmp(ls[i].angle(), ls[i+1].angle())==0)continue;
		ls[cnt++]=ls[i];
	}
	ls[cnt++]=ls[n-1];
	ls.resize(cnt);
	deque<Line> q(n);
	deque<P> p(n);
	int front=0,last=0;//[front,last]
	q[0]=ls[0];//第一个半平面
	for(int i = 1; i < cnt; i++){
		while(last-front>0 && onRight(ls[i], p[last-1]))last--;
		while(last-front>0 && onRight(ls[i], p[front]))front++;
		q[++last]=ls[i];
		if(last-front>0)p[last-1]=LLCross(q[last-1].a, q[last-1].b, q[last].a, q[last].b);
	}
	while(q.size()>1 && onRight(q[front], p[last-1]))last--;
	if(last-front<=1)return -1;
	p[last]=LLCross(q[last].a, q[last].b, q[front].a, q[front].b);
	vector<P> ps;// 核
	for(int i = front; i <= last; i++)ps.push_back(p[i]);
	return area(ps); // 面积为0多边形可能存在一条线或一个点的核
}

// 三角形ABC 内心 (角平分线交点 内切圆)
P inCenter(P A, P B, P C){
	double a = (B - C).abs(), b = (C - A).abs(), c = (A - B).abs();
	return (A * a + B * b + C * c) / (a + b + c);
}
// 三角形ABC 外心 (中垂线交点 外接圆)
P outCenter(P A, P B, P C){
	P b = B - A, c = C - A;
	double d1 = b.abs2(), d2 = c.abs2(), d3 = 2 * b.xmul(c);
	return A - P(b.y * d2 - c.y * d1, c.x * d1 - b.x * d2) / d3;
}
// 三角形ABC 垂心 (垂线交点)
P orthoCenter(P A, P B, P C){
	P ba = B - A, ca = C - A, bc = B - C;
	double Y = ba.y * ca.y * bc.y,
		   a = ca.x * ba.y - ba.x * ca.y,
		   x0 = (Y + ca.x * ba.y * B.x - ba.x * ca.y * C.x) / a,
		   y0 = -ba.x * (x0 - C.x) / ba.y + ca.y;
	return {x0, y0};
}

/**********************圆相关**********************/
// 圆o1 o2的关系 {0,1,2,3,4}分别为 {包含,内切,相交,外切,外离}
int type(P o1, db r1, P o2, db r2){
	db d = o1.disTo(o2);
	if (cmp(d, r1 + r2) == 1)return 4; // 外离
	if (cmp(d, r1 + r2) == 0)return 3; // 外切
	if (cmp(d, abs(r1 - r2)) == 1)return 2; // 相交
	if (cmp(d, abs(r1 - r2)) == 0)return 1; // 内切
	return 0;	  // 包含
}
//圆与直线p1p2的交点 相切也返回两个点
vector<P> CirLine(P o, db r, P p1, P p2){
	if (cmp(abs((o - p1).xmul(p2 - p1) / p1.disTo(p2)), r) > 0)return {};
	db x = (p1-o).dotmul(p2-p1), y = (p2-p1).abs2(), d = x*x - y*((p1-o).abs2() - r*r);
	d = max(d, (db)0.0);
	P m = p1 - (p2-p1) * (x/y), dr = (p2-p1) * (sqrt(d)/y);
	return {m - dr, m + dr}; //沿着p1p2
}
// 两圆交点 要先判断两圆是否一样
vector<P> CirCir(P o1, db r1, P o2, db r2){
	db d = o1.disTo(o2);
	if (cmp(d, r1 + r2) == 1)return {};
	if (cmp(d, abs(r1 - r2)) == -1)return {};
	d = min(d, r1 + r2);
	db y = (r1*r1 + d*d - r2*r2) / (2*d), x = sqrt(r1*r1 - y*y);
	P dr = (o2-o1).unit();
	P q1 = o1 + dr*y, q2 = dr.rota90() * x;
	return {q1 - q2, q1 + q2}; //沿着o1逆时针
}
// 求切线 r2:两圆外切线  -r2:两圆内切线   r2=0:点到圆的切线
vector<pair<P, P>> tanCC(P o1, db r1, P o2, db r2){
	P d = o2 - o1;
	db dr = r1 - r2, d2 = d.abs2(), h2 = d2 - dr * dr;
	if (sign(d2) == 0 || sign(h2) < 0)return {}; // o1==o2 || 内含 无切线
	h2 = max(0.0, h2);
	vector<pair<P, P>> res;
	for (db sign : {-1, 1}){
		P v = (d * dr + d.rota90() * sqrt(h2) * sign) / d2; // v:圆心垂直切线的单位方向向量
		res.push_back({o1 + v * r1, o2 + v * r2});
	}
	if (sign(h2) == 0)res.pop_back();
	return res;
}
//最小圆覆盖 找一个最小的圆可以覆盖所有点 O(n)
pair<P, db> min_circle(vector<P> ps){
	random_shuffle(ps.begin(), ps.end()); // 随机化 srand()一下
	int n = ps.size();
	P o = ps[0];
	db r = 0;
	for (int i = 1; i < n; i++)if (o.disTo(ps[i]) > r + EPS){
		o = ps[i], r = 0;
		for (int j = 0; j < i; j++)if (o.disTo(ps[j]) > r + EPS){
			o = (ps[i] + ps[j]) / 2;
			r = o.disTo(ps[i]);
			for (int k = 0; k < j; k++)if (o.disTo(ps[k]) > r + EPS){
				o = outCenter(ps[i], ps[j], ps[k]); //三角形外心
				r = o.disTo(ps[i]);
			}
		}
	}
	return {o, r};
}

// 向量p1 p2的夹角
db rad(P p1, P p2) { return atan2l(p1.xmul(p2), p1.dotmul(p2)); }
// 圆与三角形交的面积 圆心和三角形第三个顶点为原点  (unsure)
// 面积为有向面积 p1p2逆时针为正
db area_CirTri(db r, P p1, P p2){ 
	vector<P> is = CirLine(P(0, 0), r, p1, p2);
	if (is.empty())return r * r * rad(p1, p2) / 2;
	bool b1 = cmp(p1.abs2(), r * r) == 1, b2 = cmp(p2.abs2(), r * r) == 1;
	if (b1 && b2){ // 两点都在圆外
		if (sign((p1 - is[0]).dotmul(p2 - is[0])) <= 0 &&
			sign((p1 - is[0]).dotmul(p2 - is[0])) <= 0) // 在圆的不同边
			return r * r * (rad(p1, is[0]) + rad(is[1], p2)) / 2 + is[0].xmul(is[1]) / 2;
		else return r * r * rad(p1, p2) / 2; // 在圆的同一边
	}
	// p1在圆外 is[0]是靠近p1的交点
	if (b1) return (r * r * rad(p1, is[0]) + is[0].xmul(p2)) / 2;
	// p2在圆外 is[0]是靠近p1的交点
	if (b2) return (r * r * rad(is[1], p2) + p1.xmul(is[1])) / 2;
	return p1.xmul(p2) / 2;		// 两点都在圆内 三角形面积
}

// 圆面积并 (格林公式, 积分求面积并)
db intergal(db x, db y, db r, db L, db R){
	return r*r*(R-L) + x*r*(sinl(R)-sinl(L)) + y*r*(-cosl(R)+cosl(L));
}
db calc_area_cir(P c, db r, db L, db R){
	return intergal(c.x, c.y, r, L, R) / 2;
}
db norm(db x){ //将角度化到[0,2*PI]
	while (x < 0) x += 2 * PI;
	while (x > 2 * PI) x -= 2 * PI;
	return x;
}
P cs[N]; //圆心
db rs[N];	 //圆的半径
int m;		 //圆的数量
void work(){
	vector<int> cand;
	//去除包含的圆 重叠的保留1个
	for (int i = 0; i < m; i++){
		bool ok = 1;
		for (int j = 0; j < m; j++)if (i != j){
			if (rs[j] > rs[i] + EPS && rs[i] + cs[i].disTo(cs[j]) <= rs[j] + EPS){ //圆包含
				ok = 0;break;
			}
			if (cs[i] == cs[j] && cmp(rs[i], rs[j]) == 0 && j < i){ //圆重叠
				ok = 0;break;
			}
		}
		if (ok)cand.push_back(i);
	}
	m = cand.size();
	for (int i = 0; i < m; i++)cs[i] = cs[cand[i]], rs[i] = rs[cand[i]];

	db area = 0;

	for (int i = 0; i < m; i++){
		vector<pair<db, int>> ev = {{0, 0}, {2 * PI, 0}}; //扫描线
		int cur = 0;
		for (int j = 0; j < m; j++)if (j != i){
			vector<P> ret = CirCir(cs[i], rs[i], cs[j], rs[j]);
			if (!ret.empty()){
				db l = (ret[0] - cs[i]).alpha();
				db r = (ret[1] - cs[i]).alpha();
				l = norm(l), r = norm(r); //转化到[0,2*PI]
				ev.push_back({l, 1});
				ev.push_back({r, -1});
				if (l > r)++cur;
			}
		}

		sort(ev.begin(), ev.end());
		for (int j = 0; j < ev.size() - 1; j++){
			cur += ev[j].se;
			if (cur == 0){
				area += calc_area_cir(cs[i], rs[i], ev[j].fr, ev[j + 1].fr);
			}
		}
	}
	// return area;
}

// 平面最近点对
bool cmp_y(const P& a, const P& b) { return a.y < b.y; }
P a[N];
double mindist=1e99;
int ansa, ansb;// 最近点对的下标
void upd_ans(P pa,P pb){
	double dis = pa.disTo(pb);
	if(cmp(dis, mindist)<0)mindist=dis, ansa=pa.id, ansb=pb.id;
}
void NearestPoint(int l,int r){
	if(r-l<=3){
		for(int i = l; i <= r; i++)
			for(int j = i+1; j <= r; j++)
				upd_ans(a[i], a[j]);
		sort(a+l,a+r+1,cmp_y);return;
	}
	int mid = (l+r)>>1;
	db midx = a[mid].x;
	NearestPoint(l,mid);
	NearestPoint(mid+1,r);
	inplace_merge(a+l, a+mid+1, a+r+1, cmp_y);//归并两个数组

	static P t[N];
	int tsz=0;
	for (int i = l; i <= r; ++i){
		if (abs(a[i].x - midx) < mindist) {
			for (int j = tsz - 1; j >= 0 && a[i].y - t[j].y < mindist; --j)
				upd_ans(a[i], t[j]);
			t[tsz++] = a[i];
		}
	}
}
void work(int n){
	sort(a+1,a+1+n);mindist=1e99;
	NearestPoint(1,n);
}

// 最小矩形覆盖 返回{面积,矩形(点逆时针 且y最小的在第一个)} （凸包ps)
pair<db,vector<P>> min_rectangle(vector<P> ps)
{
	vector<P> pans;
	int n=ps.size();
	for(int i = 0; i < n; i++)
		ps.push_back(ps[i]);//复制一份
	int l=0,u=0,r=0;
	P d=ps[1]-ps[0];// 最小矩形必有一条边与多边形一边重合
	for(int i = 0; i < n; i++){//找一个最初解
		if(cmp(d.dotmul(ps[i]), d.dotmul(ps[l]))<0)l=i;
		if(cmp(d.dotmul(ps[i]), d.dotmul(ps[r]))>0)r=i;
		if(cmp(d.xmul(ps[i]), d.xmul(ps[u]))>0)u=i;
	}
	db ans=1e99;//矩形面积
	for(int i = 0; i < n; i++){
		d=ps[i+1]-ps[i];
		while(cmp(d.dotmul(ps[l+1]), d.dotmul(ps[l]))<0)l++;//叉积最小
		while(cmp(d.dotmul(ps[r+1]), d.dotmul(ps[r]))>0)r++;//叉积最大
		while(cmp(d.xmul(ps[u+1]), d.xmul(ps[u]))>0)u++;// 最高的点
		assert(l<2*n && r<2*n && u<2*n);
		db x1=d.dotmul(ps[r]-ps[l])/d.abs();
		db x2=d.xmul(ps[u]-ps[i])/d.abs();
		db tmp=x1*x2;
		if(cmp(tmp,ans)<0){// 设ps[i]ps[i+1]矩形的右下边
			ans=tmp;
			P p1=LLCross(ps[i], ps[i]+d, ps[r], ps[r]+d.rota90());// 右顶点 
			P p2=LLCross(ps[r], ps[r]+d.rota90(), ps[u], ps[u]-d);// 上顶点
			P p3=LLCross(ps[u], ps[u]-d, ps[l], ps[l]-d.rota90());// 左顶点
			P p4=LLCross(ps[l], ps[l]-d.rota90(), ps[i], ps[i]+d);// 下顶点
			pans = {p1,p2,p3,p4};// 已按逆时针
		}
	}
	auto cmp1 = [](const P a,const P b){
		if(cmp(a.y,b.y)==0)return cmp(a.x,b.x)<0;
		return cmp(a.y,b.y)<0;
	};
	//按y关键字排序最小的放第一个
	rotate(pans.begin(), min_element(pans.begin(),pans.end(),cmp1), pans.end());
	return {ans,pans};
}
}