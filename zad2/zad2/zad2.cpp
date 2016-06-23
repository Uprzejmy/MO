// zad2.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "iostream"

using namespace std;

double taylorValue(double x, int n)
{
	double sum = 1.0; //pierwszy (n=0) wyraz szeregu
	double addend = 1.0;

	if (x == 0)
		return 1;
	else if (x > 0)
	{
		for (int i = 1; i <= n; ++i)
		{
			addend = addend * x / i;
			sum = sum + addend;
		}
	}

	else
	{
		return 1 / taylorValue(-x, n);
	}
	

	return sum;
}

void print(double x, int n)
{
	double accurate = exp(x);
	double rounded = taylorValue(x, n);
	double nextRounded = taylorValue(x, n + 1);

	cout << "x: " << x << endl;
	cout << "n: " << n << endl;

	cout << "Taylor(x): ";
	printf("%.10e", rounded);
	cout << endl;

	cout << "exp(x):    ";
	printf("%.10e", accurate);
	cout << endl;

	cout << "blad wzgledny: ";
	printf("%.10e", (rounded - accurate)/accurate);
	cout << endl;

	cout << "blad wzgledny (taylor(x,n+1)) - blad wzgledny (taylor(x,n)): ";
	printf("%.10e", ((nextRounded - accurate) / accurate) - ((rounded - accurate) / accurate));
	cout << endl;

	cout << endl;
}

int _tmain()
{
	int n = 150;

	for (double x = -30.0; x < 31.0; x = x + 1.0)
	{
		print(x, n);
	}

	return 0;
}

