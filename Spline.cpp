#include<stdio.h>
#include<conio.h>
#include<math.h>

double f1(double x);
double pn(double xx, double x[], int n);
double f0(double x);
double f1(double x);

int main() {
	int i, n = 6;
	double xx, x[10];
	double xn[1001], yn[1001];

	for (i = 0; i <= n; i++) x[i] = (double)i * 10 / 6 - 5;
	scanf_s("%lf", &xx, 7);
	if (xx < -5 || xx>5)printf("error");
	else printf("%lf", pn(xx, x, n));
 
	for (i = 0; i <= 1000; i++) {
		xn[i] = i*0.01 - 5;
		yn[i] = pn(xn[i], x, n);
		printf("\n%lf,%lf", xn[i], yn[i]);
	}

	return 0;
}

double pn(double xx, double x[], int n) {
	int j;
	double h[20] = { 0 }, m[20] = { 0 }, a[20] = { 0 }, b[20] = { 0 }, c[20][20] = { 0 }, r[20] = { 0 };
	double d1, d2;

	//求出方程组矩阵，未知数为m
	for (j = 1; j <= n; j++)h[j] = x[j] - x[j - 1];
	for (j = 1; j < n; j++) {
		a[j] = h[j] / (h[j] + h[j + 1]);
		b[j] = 6 / (h[j] + h[j + 1])*(((f0(x[j + 1]) - f0(x[j])) / h[j + 1]) - ((f0(x[j]) - f0(x[j - 1])) / h[j]));
		c[j][j-1] = a[j];
		c[j][j] = 2;
		c[j][j+1] = 1 - a[j];
		c[j][n+1] = b[j];
	}
	c[0][0] = 2;
	c[0][1] = 1;
	c[0][n + 1] = 6 * ((f0(x[1]) - f0(x[0])) / h[1] - (f1(x[0]))) / h[1];
	c[n][n] = 2;
	c[n][n - 1] = 1;
	c[n][n + 1] = 6 * ((f0(x[n]) - f0(x[n - 1])) / h[n-1] - (f1(x[n]))) / h[n];

	//用追赶法求解方程组
	c[0][1] = c[0][1] / c[0][0];
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

	//至此求出m
	for (j = 0; x[j] <= xx;) {
		if (x[j++] == xx)return f0(xx);
	}
	d1 = x[j] - xx;
	d2 = xx - x[j - 1];
	return (m[j - 1] * d1*d1*d1 + m[j] * d2*d2*d2) / 6 / h[j]\
		+ (f0(x[j]) / h[j] - h[j] * m[j] / 6)*d2\
		+ (f0(x[j - 1]) / h[j] - h[j] * m[j - 1] / 6)*d1;
}

double f0(double x) {
	return 1 / (x*x + 1);
}

double f1(double x) {
	return -2 * x / (x*x + 1) / (x*x + 1);
}