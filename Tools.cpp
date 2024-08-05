#include <bits/stdc++.h>
#define ll long long
#define ull unsigned long long
#define int ll
#define fr first
#define se second
#define INF 0x3f3f3f3f
#define LINF 0x3f3f3f3f3f3f3f3f
#define For(i, a, b) for (int i = (a); (i) <= (b); ++i)
#define Rep(i, a, b) for (int i = (a); (i) >= (b); --i)
using namespace std;
typedef pair<int, int> pii;
const ll N = 2e5 + 10;

double abs(vector<double> &v) // 绝对偏差
{
	double ave = 0;
	for (auto c : v)
		ave += c;
	ave /= v.size();
	double res = 0;
	for (auto c : v)
		res += (fabs(ave - c));
	return res / v.size();
}

double average(const vector<double> &v) // 平均值
{
	double ave = 0;
	for (auto c : v)
		ave += c;
	ave /= v.size();
	return ave;
}

double StandardDeviation(const vector<double>& v)// 标准差
{
	double mean = average(v), res = 0;
	for(auto c:v)res += (mean - c) * (mean - c);
	res /= (double)v.size();
	res = sqrt(res);
	return res;
}


namespace myIO{
// 输出
void debug(){cerr << "\n";}
template<typename T, typename... Args>
void debug(T head, Args... tail)
{
	cerr << head << " ";
	debug(tail...);
}
#define dbg(...) cerr << "[" << #__VA_ARGS__ << "]: ", debug(__VA_ARGS__)

// 交互输入
void interScan(string s){}
template<class T, class ...Args>
void interScan(string s, T& head, Args&... args)
{
	cout << s;
	cin >> head;
	interScan("", args...);
}
}


#include <io.h>
void GetFileName(string path, vector<string> &filesName) //返回指定文件夹下所有文件名
{
	long hFile = 0;				 //文件句柄
	struct _finddata_t fileinfo; //定义文件信息结构体
	string p;
	if ((hFile = _findfirst(p.assign(path).append("\\*").c_str(), &fileinfo)) != -1) //使用函数_findfirst()打开文件并获取第一个文件名
	{
		do
		{
			if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0) //"."表示当前目录，".."表示父目录
				filesName.push_back(fileinfo.name);
		} while (_findnext(hFile, &fileinfo) == 0); //使用函数_findnext()继续获取其他文件名
		_findclose(hFile);							//使用函数_findclose()关闭文件夹
	}
}

void solve()
{

}

signed main()
{
	// ios::sync_with_stdio(false),cin.tie(0),cout.tie(0);
	// freopen("in.txt", "r", stdin);
	// freopen("out.txt", "w", stdout);
	int tt = 1;
	// cin >> tt;
	while (tt--)
		solve();
	return 0;
}