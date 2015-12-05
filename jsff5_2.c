#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<windows.h>
#define A 0.0 //积分上下限 
#define B 100.0
#define N 100 //区间数目
#define epsilon 1e-6
 
typedef struct node{
double fx;
struct node *next;
}NODE;
double X;
double func(double t);
double f(double t);
double simpson(NODE *header,int n); //输入储存有节点处函数值的链表 ，返回辛普森积分值 header是链表头，n是区间数 
double step_simpson(NODE *header); //输入储存有节点处函数值的链表，对链表插入中间值，调用simpson（）计算积分，直到满足精度要求 
double romberg(void);  
double trapzoid(NODE *header,int n);//参数与功能同于simpson函数 
void double_list(NODE *header,int n); //分半一个链表 n为插值后的节点数 
double gauss_raguel(void);

int main()
{
	NODE *header1,*header2,*p;
	int n=0;
	FILE *fp;
	double x;
	time_t start,end;
	X=1;
	header1=(NODE*)malloc(sizeof(NODE));
	header2=(NODE*)malloc(sizeof(NODE));
	p=header1;
	for(x=A,n=0;n<N+1;n++,x+=(B-A)/N) //将每个节点处的函数值储存在链表中 ，节点数是网格数加一 
	{
		p->fx=func(x);
		p->next=(NODE*)malloc(sizeof(NODE));
		p=p->next;
	}
	p->next=NULL;
	/*p=header2;
	for(x=A,n=0;n<2*N+1;n++,x+=(B-A)/(2*N)) //辛普森法用的节点，节点数是网格数的二倍加一 
	{
		p->fx=func(x);
		p->next=(NODE*)malloc(sizeof(NODE));
		p=p->next;
	}
	p->next=NULL;*/
	fp=fopen("jsff5_2_plot.csv","w");
	fprintf(fp,"%s,%s,%s,%s\n","x","trapzoid","inside","gauss_raguel");
	for(X=1;X<10;X++)
	{
	fprintf(fp,"%lf,",X);
	for(x=A,n=0,p=header1;n<N+1;n++,x+=(B-A)/N) 
	{
		p->fx=func(x);
		p=p->next;
	}
	start=clock();
	fprintf(fp,"%lf,",trapzoid(header1,N));
	fprintf(fp,"%lf,",tgamma(X));
	fprintf(fp,"%lf\n",gauss_raguel());
	end=clock();
	printf("%lf\n",(double)start);
	}
	return 0;
}
double gauss_raguel(void)
{
	return f(0.3225476896)*0.60315410434+f(1.7457611012)*0.35741869244+f(4.5366202969)*3.8887908515e-2+f(9.3950709123)*5.3929470556e-4;
}
double f(double t)
{
	return pow(t,X-1);
}
double func(double t)
{
	return exp(-t)*pow(t,X-1);
}

double simpson(NODE *header,int n)  
{
	double y;
	NODE *p;
	int flag=1;
	y=(1.0/6)*header->fx;
	p=header->next;
	while(p->next!=NULL)
	{
		if(flag%2==1)
			{
			y+=p->fx*4.0/6;
			}
		else
			{y+=p->fx*2.0/6;
			}
		p=p->next;
		flag++;
	}
	y+=p->fx*1.0/6;
	return y*(B-A)/n;
}
double trapzoid(NODE *header,int n)
{
	double y;
	NODE *p;
	y=0.5*header->fx;
	p=header->next;
	while(p->next!=NULL)
	{
		y+=p->fx;
		p=p->next;
	}
	y+=p->fx*0.5;
	return y*(B-A)/n;
}

double step_simpson(NODE *header)
{
	int n=1,flag;
	double f1,f2,x;
	NODE *p,*temp;
	f1=simpson(header,N);
	f2=f1;
	do{
		f1=f2;
		/*for(flag=1,x=A,p=header;flag<=pow(2,n)*N+1;flag++,x+=(B-A)/(pow(2,n)*N))
		{
			if(flag%2) //在两个旧节点之间插入新节点，旧节点不必再计算 
			{
				temp=p->next;   
				p->next=(NODE*)malloc(sizeof(NODE));
				p=p->next;
			}
			else
			{
				p->next=temp;
				p->fx=func(x);
				p=p->next;
			}		
		}
		p->next=NULL;
		*/
		double_list(header,pow(2,n)*N);
		double_list(header,pow(2,n)*N*2);
		f2=simpson(header,pow(2,n)*N);
		n++;
		}while(fabs(f2-f1)/15>epsilon);
		printf("change step simpson time : %d\n",n);
		return f2;
}
void double_list(NODE *header,int n)
{
	int flag;
	double x;
	NODE *p,*temp;
	for(flag=1,x=A,p=header;flag<=n+1;flag++,x+=(B-A)/n)
		{
			if(flag%2) //在两个旧节点之间插入新节点，旧节点不必再计算 
			{
				temp=p->next;   
				p->next=(NODE*)malloc(sizeof(NODE));
				p=p->next;
			}
			else
			{
				p->next=temp;
				p->fx=func(x);
				p=p->next;
			}		
		}
		p->next=NULL;
	return;
}
	
double romberg(void) 
{
	double T[4][2]={0};
	int n=0,flag;
	NODE *header;
	header=(NODE*)malloc(sizeof(NODE));
	header->next=(NODE*)malloc(sizeof(NODE));
	header->fx=func(A);
	header->next->fx=func(B);
	header->next->next=NULL;
	T[0][0]=trapzoid(header,pow(2,n++));
	while(fabs(T[3][1]-T[2][1])>epsilon || T[3][1]==0) 
	{
		double_list(header,pow(2,n));
		T[0][1]=trapzoid(header,pow(2,n++));
		for(flag=1;flag<4;flag++)  
		{
			if(T[flag][0]==0 && T[flag][1]==0)
				{
				T[flag][0]=(pow(4,flag)*T[flag-1][1]-T[flag-1][0])/(pow(4,flag)-1);
				break;
				}
			else if(T[flag][0]!=0 && T[flag][1]==0)
				T[flag][1]=(pow(4,flag)*T[flag-1][1]-T[flag-1][0])/(pow(4,flag)-1);
			else
				{
				T[flag][0]=T[flag][1];
				T[flag][1]=(pow(4,flag)*T[flag-1][1]-T[flag-1][0])/(pow(4,flag)-1);
				}
		}
		T[0][0]=T[0][1];
	}
	printf("romberg time : %d\n",n);
	return T[3][1];
}
