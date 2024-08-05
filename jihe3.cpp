#define eps 1e-8
#define zero(x) (((x)>0?(x):-(x))<eps)
struct point3{double x,y,z;};
struct line3{point3 a,b;};
struct plane3{point3 a,b,c;};
point3 xmult(point3 u,point3 v){//叉乘
	point3 ret;
	ret.x=u.y*v.z-v.y*u.z;
	ret.y=u.z*v.x-u.x*v.z;
	ret.z=u.x*v.y-u.y*v.x;
	return ret;
}
//计算 dot product product product product U . V
double dmult(point3 u,point3 v){
	return u.x*v.x+u.y*v.y+u.z*v.z;
}
//矢量差 U - V
point3 subt(point3 u,point3 v){
	point3 ret;
	ret.x=u.x-v.x;
	ret.y=u.y-v.y;
	ret.z=u.z-v.z;
	return ret;
}
//取法向量
point3 pvec(plane3 s){
	return xmult(subt(s.a,s.b),subt(s.b,s.c));
}
point3 pvec(point3 s1,point3 s2,point3 s3){
	return xmult(subt(s1,s2),subt(s2,s3));
}
point3 intersection(point3 l1,point3 l2,point3 s1,point3 s2,point3 s3){
	point3 ret=pvec(s1,s2,s3);
	double t=(ret.x*(s1.x-l1.x)+ret.y*(s1.y-l1.y)+ret.z*(s1.z-l1.z))/
	(ret.x*(l2.x-l1.x)+ret.y*(l2.y-l1.y)+ret.z*(l2.z-l1.z));
	ret.x=l1.x+(l2.x-l1.x)*t;
	ret.y=l1.y+(l2.y-l1.y)*t;
	ret.z=l1.z+(l2.z-l1.z)*t;
	return ret;
}
double vlen(point3 p){
	return sqrt(p.x*p.x+p.y*p.y+p.z*p.z);
}
int dots_inline(point3 p1,point3 p2,point3 p3){
	return vlen(xmult(subt(p1,p2),subt(p2,p3)))<eps;
}
//直线平行
int parallel(line3 u,line3 v){
return vlen(xmult(subt(u.a,u.b),subt(v.a,v.b)))<eps;
}
int parallel(point3 u1,point3 u2,point3 v1,point3 v2){
return vlen(xmult(subt(u1,u2),subt(v1,v2)))<eps;
}
//平面平行
int parallel(plane3 u,plane3 v){
	return vlen(xmult(pvec(u),pvec(v)))<eps;
}
int parallel(point3 u1,point3 u2,point3 u3,point3 v1,point3 v2,point3 v3){
	return vlen(xmult(pvec(u1,u2,u3),pvec(v1,v2,v3)))<eps;
}
//线面平行
int parallel(line3 l,plane3 s){
	return zero(dmult(subt(l.a,l.b),pvec(s)));
}
int parallel(point3 l1,point3 l2,point3 s1,point3 s2,point3 s3){
	return zero(dmult(subt(l1,l2),pvec(s1,s2,s3)));
}

//面交线
line3 intersection(plane3 u,plane3 v){
	line3 ret;
	ret.a=parallel(v.a,v.b,u.a,u.b,u.c)?intersection(v.b,v.c,u.a,u.b,u.c):intersection(v.a,v.b,u.a,u.b,u.c);
	ret.b=parallel(v.c,v.a,u.a,u.b,u.c)?intersection(v.b,v.c,u.a,u.b,u.c):intersection(v.c,v.a,u.a,u.b,u.c);
	return ret;
}
line3 intersection(point3 u1,point3 u2,point3 u3,point3 v1,point3 v2,point3 v3){
	line3 ret;
	ret.a=parallel(v1,v2,u1,u2,u3)?intersection(v2,v3,u1,u2,u3):intersection(v1,v2,u1,u2,u3);
	ret.b=parallel(v3,v1,u1,u2,u3)?intersection(v2,v3,u1,u2,u3):intersection(v3,v1,u1,u2,u3);
	return ret;
}
double ptoplane(point3 p,plane3 s){
	return fabs(dmult(pvec(s),subt(p,s.a)))/vlen(pvec(s));
}
double ptoplane(point3 p,point3 s1,point3 s2,point3 s3){
	return fabs(dmult(pvec(s1,s2,s3),subt(p,s1)))/vlen(pvec(s1,s2,s3));
}

//计算直线与平面交点,注意事先判断是否平行,并保证三点不共线!
point3 intersection(line3 l, plane3 s)
{
	point3 ret = pvec(s);
	double t = (ret.x * (s.a.x - l.a.x) + ret.y * (s.a.y - l.a.y) + ret.z * (s.a.z - l.a.z)) /
			   (ret.x * (l.b.x - l.a.x) + ret.y * (l.b.y - l.a.y) + ret.z * (l.b.z - l.a.z));
	ret.x = l.a.x + (l.b.x - l.a.x) * t;
	ret.y = l.a.y + (l.b.y - l.a.y) * t;
	ret.z = l.a.z + (l.b.z - l.a.z) * t;
	return ret;
}
point3 intersection(point3 l1, point3 l2, point3 s1, point3 s2, point3 s3)
{
	point3 ret = pvec(s1, s2, s3);
	double t = (ret.x * (s1.x - l1.x) + ret.y * (s1.y - l1.y) + ret.z * (s1.z - l1.z)) /
			   (ret.x * (l2.x - l1.x) + ret.y * (l2.y - l1.y) + ret.z * (l2.z - l1.z));
	ret.x = l1.x + (l2.x - l1.x) * t;
	ret.y = l1.y + (l2.y - l1.y) * t;
	ret.z = l1.z + (l2.z - l1.z) * t;
	return ret;
}










