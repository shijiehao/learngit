#include<stdio.h>
#include<math.h>

#define N 15
#define T 10
#define O 1.8
void gauss(double hilbert[N][N+1],double* x);
void jacobi(double hilbert[N][N+1],double* x);
void gauss_seidel(double hilbert[N][N+1],double* x);
void sor(double hilbert[N][N+1],double* x);
int main()
{
	double hilbert[N][N+1],x_g[N]={0},x_j[N],x_gs[N],x_sor[N];
	int row,col;
	for(row=0;row<N;row++)
	for(col=0;col<N;col++)
	{
		hilbert[row][col]=1.0/(row+col+1);  //生成希尔伯特矩阵 
		hilbert[row][N]+=hilbert[row][col];
	}
	gauss(hilbert,x_g);	
	jacobi(hilbert,x_j);
	gauss_seidel(hilbert,x_gs);
	sor(hilbert,x_sor);
	printf("*****************************************************\n");
	printf("%15s%15s%15s%15s\n","Gauss","Jacobi","Gauss_Seidel","SOR");
	for(col=0;col<N;col++)
	{
		printf("%d",col);
		printf("%15lf  ",x_g[col]);
		printf("%15lf",x_j[col]);
		printf("%15lf",x_gs[col]);
		printf("%15lf",x_sor[col]);
		printf("\n");
	}
	printf("*****************************************************");
	return 0;
}

void gauss(double hilbert[N][N+1],double* x)
{
	int r1,r2,r3;
	double h[N][N+1],temp=0,temp1;
	
	for(r1=0;r1<N;r1++)
	for(r2=0;r2<N+1;r2++)
	h[r1][r2]=hilbert[r1][r2];
	
	for(r1=0;r1<N-1;r1++) //求三角矩阵 
		for(r2=r1+1;r2<N;r2++)
		{
			temp1=h[r2][r1]/h[r1][r1];
			for(r3=r1;r3<N+1;r3++)
				h[r2][r3]-=temp1*h[r1][r3];
		}
	for(r1=N-1;r1>=0;r1--) //求解 
	{
		temp=0;
		for(r2=0;r2<N;r2++)
		temp+=h[r1][r2]*x[r2];
		x[r1]=(h[r1][N]-temp)/h[r1][r1];
	}
}

void jacobi(double hilbert[N][N+1],double* x)
{ 
	int r1,r2,r3;
	double xi[N]={0};
	for(r1=0;r1<T;r1++)
	{
		for(r2=0;r2<N;r2++)
		{
			x[r2]=hilbert[r2][N]/hilbert[r2][r2];
			for(r3=0;r3<N;r3++)
			if(r3!=r2)
			x[r2]-=hilbert[r2][r3]*xi[r3]/hilbert[r2][r2];
		}
		for(r3=0;r3<N;r3++)
		xi[r3]=x[r3];
	}
}

void gauss_seidel(double hilbert[N][N+1],double* x)
{ 
	int r1,r2,r3;
	for(r1=0;r1<N;r1++)
		x[r1]=0;
	for(r1=0;r1<T;r1++)
	{
		for(r2=0;r2<N;r2++)
		{
			x[r2]=hilbert[r2][N]/hilbert[r2][r2];
			for(r3=0;r3<N;r3++)
			if(r3!=r2)
			x[r2]-=hilbert[r2][r3]*x[r3]/hilbert[r2][r2];
		}
	}
}

void sor(double hilbert[N][N+1],double* x)
{
	int r1,r2,r3;
	double xi[N]={0};
	for(r1=0;r1<T;r1++)
	{
		for(r2=0;r2<N;r2++)
		{
			x[r2]=hilbert[r2][N]/hilbert[r2][r2];
			for(r3=0;r3<N;r3++)
			if(r3!=r2)
			x[r2]-=hilbert[r2][r3]*x[r3]/hilbert[r2][r2];
			x[r2]=x[r2]*O+xi[r2]*(1-O);
			//printf("%lf ",x[r2]);
		}
		for(r3=0;r3<N;r3++)
		xi[r3]=x[r3];
	}
}

	
