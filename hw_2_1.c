#include<stdio.h>
#include<math.h>

double bot=-5,top=5;
int n=7;
struct point{
	double x;
	double y;
	double y_;
};

double func(double x)
{
	return 1/(1+x*x);
}
double func_(double x)
{
	return -2*x/((1+x*x)*(1+x*x));
}

double lagrange(struct point *pbr,double x)  //输入指针和x值，返回x处y的拉格朗日插值 
{
	double y=0,y_s;  //psr pbr 代表small round 和 big round ，y_s代表y_single ，埃米尔特函数中意义相同 
	struct point *psr,*p_static=pbr;
	for(;pbr<p_static+n;pbr++)  
	{
		y_s=pbr->y;
		
		if(pbr!=p_static+n-1) 
		for(psr=pbr+1;psr!=pbr;)
		{
			y_s=y_s*(x-psr->x)/(pbr->x-psr->x);
			if(psr==p_static+n-1)
			psr=p_static;
			else
			psr++;
		}
		else
		for(psr=p_static;psr!=pbr;psr++)
			y_s=y_s*(x-psr->x)/(pbr->x-psr->x);
		y+=y_s;
	}
	return y;
}

double hermite(struct point *pbr,double x)
{
	double y=0,y_s,y_s_;
	struct point *psr,*p_static=pbr;
	for(;pbr<p_static+n;pbr++)  
	{
		y_s=1;
		y_s_=0;
		
		
		if(pbr!=p_static+n-1) 			//这一部分给L(x)赋值 
		for(psr=pbr+1;psr!=pbr;)
		{
			y_s=y_s*(x-psr->x)/(pbr->x-psr->x);
			if(psr==p_static+n-1)
			psr=p_static;
			else
			psr++;
		}
		else
		for(psr=p_static;psr!=pbr;psr++)
		y_s=y_s*(x-psr->x)/(pbr->x-psr->x);
			
			
		if(pbr!=p_static+n-1)			//此处给L'(x)赋值 
		for(psr=pbr+1;psr!=pbr;)
		{
			y_s_+=y_s/(x-psr->x);
			if(psr==p_static+n-1)
			psr=p_static;
			else
			psr++;
		}
		else
		for(psr=p_static;psr!=pbr;psr++)
		y_s_+=y_s/(x-psr->x);
		
		
		y+=(1-2*y_s_*(x-pbr->x))*y_s*y_s*pbr->y+(x-pbr->x)*y_s*y_s*pbr->y_;
	}
	return y;
}

double spline(struct point *p,double x)
{
	int j;
	double h[20] = { 0 }, m[20] = { 0 }, a[20] = { 0 }, b[20] = { 0 }, c[20][20] = { 0 }, r[20] = { 0 };
	double d1, d2;
	
	for (j = 1; j <= n; j++)h[j] = p[j].x - p[j-1].x;  //求系数矩阵 
	for (j = 1; j < n; j++) {
		a[j] = h[j] / (h[j] + h[j + 1]);
		b[j] = 6 / (h[j] + h[j + 1])*(((func(p[j + 1].x) - func(p[j].x)) / h[j + 1]) - ((func(p[j].x) - func(p[j - 1].x)) / h[j]));
		c[j][j-1] = a[j];
		c[j][j] = 2;
		c[j][j+1] = 1 - a[j];
		c[j][n+1] = b[j];
	}
	c[0][0] = 2;
	c[0][1] = 1;
	c[0][n + 1] = 6 * ((func(p[1].x) - func(p[0].x)) / h[1] - (func_(p[0].x))) / h[1];
	c[n][n] = 2;
	c[n][n - 1] = 1;
	c[n][n + 1] = 6 * ((func(p[n].x) - func(p[n - 1].x)) / h[n-1] - (func_(p[n].x))) / h[n];


	c[0][1] = c[0][1] / c[0][0]; //解方程组
	c[0][n + 1] = c[0][n + 1] / c[0][0];
	c[0][0] = 1;
	for (j = 1; j <= n; j++) {
		c[j][n + 1] = (c[j][n + 1] - c[j - 1][n + 1] * c[j][j - 1]) / (c[j][j] - c[j-1][j]*c[j][j-1]);
		if (j != n)c[j][j + 1] = c[j][j + 1] / (c[j][j] - c[j - 1][j] * c[j][j - 1]);
		c[j][j] = 1;
		c[j][j - 1] = 0;
	}
	m[n] = c[n][n + 1];
	for (j = n - 1; j >= 0; j--)m[j] = c[j][n + 1] - c[j][j + 1] * m[j + 1];

	for (j = 0; p[j].x <= x;) {
		if (p[j++].x == x)return func(x);
	}
	d1 = p[j].x - x;
	d2 = x - p[j - 1].x;
	return (m[j - 1] * d1*d1*d1 + m[j] * d2*d2*d2) / 6 / h[j]\
		+ (func(p[j].x) / h[j] - h[j] * m[j] / 6)*d2\
		+ (func(p[j - 1].x) / h[j] - h[j] * m[j - 1] / 6)*d1;
}

int main()
{
	struct point eql[n],cheb[n],*p;
	double x;
	int i;
	FILE *fp;
	for(p=eql;p<eql+n;p++)
	{
		if(p==eql)
		p->x=bot;
		else
		p->x=(p-1)->x+(top-bot)/(n-1);
		p->y=func(p->x);
		p->y_=func_(p->x);
	}
	for(p=cheb,i=0;p<cheb+n;p++,i++)
	{
		p->x=(top+bot)/2.0+(top-bot)*cos((2*i+1.0)/(2.0*n)*3.14)/2.0;
		p->y=func(p->x);
		p->y_=func_(p->x);
	}
	/*printf("%lf\n",hermite(eql,test));
	printf("%lf\n",lagrange(eql,test));
	printf("%lf\n",lagrange(cheb,test));
	printf("%lf\n",func(test));*/
	fp=fopen("data.csv","w");
	fprintf(fp,"%s,%s,%s,%s,%s,%s\n","x","origin","lagrange_eql","lagrange_cheb","hermite","spline");
	for(x=-5;x<=5;x+=0.01)
	fprintf(fp,"%lf,%lf,%lf,%lf,%lf,%lf\n",x,func(x),lagrange(eql,x),lagrange(cheb,x),hermite(eql,x),spline(eql,x));
	return 0;
}
