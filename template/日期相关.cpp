struct Date{
	int year, month, day;
};
int days[12]={31,28,31,30,31,30,31,31,30,31,30,31};

// 闰年
inline int leap(int year){
	return (year%4==0&&year%100!=0)||year%400==0;
}

//返回指定日期是星期几
int weekday(Date a){
	int tm = a.month>=3 ? (a.month-2) : (a.month+10);
	int ty = a.month>=3 ? a.year : (a.year-1);
	return (ty+ty/4-ty/100+ty/400+(int)(2.6*tm-0.2)+a.day)%7;
}

// 日期转天数 0/1/1 为第一天
int date2int(Date a){
	int ret = a.year*365+(a.year-1)/4-(a.year-1)/100+(a.year-1)/400;
	days[1] += leap(a.year);
	for(int i=0; i<a.month-1; ret+=days[i++]){}
	days[1] = 28;
	return ret+a.day;
}

//天数偏移转日期
Date int2date(int a){
	Date ret;
	ret.year = a/146097*400;
	for(a%=146097; a>=365+leap(ret.year); a-=365+leap(ret.year),ret.year++);
	days[1] += leap(ret.year);
	for(ret.month=1; a>=days[ret.month-1]; a-=days[ret.month-1],ret.month++);
	days[1] = 28;
	ret.day = a+1;
	return ret;
}