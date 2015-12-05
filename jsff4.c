/*************************************************************************************
						计算方法第四次作业
				二分法、牛顿法、双点弦截法解非线性方程
**************************************************************************************/

#include<stdio.h>
#include<math.h>

#define A_DICHOTOMY -10000
#define B_DICHOTOMY 10000
#define ACCURACY 10e-6
#define INI_NEWTON -100000
#define A_SECANT -1000
#define B_SECANT 10

double func(double x)
{
	return x*x-exp(x);
}

double func_(double x)
{
	return 2*x-exp(x);
}

double dichotomy(double a,double b);
double newton(double x1);
double secant(double x1,double x2);

int main(void)
{
	double a=A_DICHOTOMY,b=B_DICHOTOMY,ini=INI_NEWTON,as=A_SECANT,bs=B_SECANT;
	printf("dichotomy : %lf\n",dichotomy(a,b));
	printf("newton : %lf\n",newton(ini));
	printf("secant : %lf\n",secant(as,bs));
	return 0;
}

double dichotomy(double a,double b)
{
	double c;
	int n=0;
	while(1)
	{	
		n++;
		c=(a+b)/2;
		if(fabs(func(c))==0)
			break;
		else if(func(c)*func(a)>0)
			a=c;
		else
			b=c;
		if(fabs(b-a)<ACCURACY)
		{
			c=(a+b)/2;
			break;
		}
	}
	printf("dichotomy times : %d\n",n);
	return c;
}

double newton(double x1)
{
	double x2;
	int n=0;
	while(1)
	{
		n++;
		x2=x1-func(x1)/func_(x1);
		if(fabs(func(x2))<ACCURACY)
			break;
		x1=x2;
	}
	printf("newton times : %d\n",n);
	return x2;
}

double secant(double x1,double x2)
{
	double x3;
	int n=0;
	while(1)
	{
		n++;
		x3=x2-(x2-x1)*func(x2)/(func(x2)-func(x1));
		if(fabs(func(x3))<ACCURACY)
			break;
		x1=x2,x2=x3;
	}
	printf("secant times : %d\n",n);
	return x3;
}
