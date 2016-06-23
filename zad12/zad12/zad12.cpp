// zad12.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#define _USE_MATH_DEFINES

#include <math.h>
#include <iostream>
#include <iomanip>

using namespace std;

double fx(double x)
{
	double absX = abs(x);
	return x / (1 + 25 * pow(absX, 3));
}

double lagrangeInterpolation(double* x0, double* y0, double x, int n)
{
	double product;
	double y = 0.0;

	for (int i = 0; i < n; i++)
	{
		product = 1.0;
		for (int j = 0; j < n; j++)
		{
			if (i != j)
			{
				product *= (x - x0[j]) / (x0[i] - x0[j]);
			}
		}
		y += product * y0[i];
	}

	return y;
}

void equalDistantNodes(double* x, double xp, double xk, int n)
{
	double interval = (xk - xp) / n;

	for (int i = 0; i <n; i++)
	{
		x[i] = xp + i*interval;
	}
}

void chebyshevNodes(double* x, int n)
{
	for (int k = 1; k <= n; k++)
	{
		x[k - 1] = cos((2.0*k+1.0) / (2.0*n+2.0) * M_PI);
	}
}

int _tmain(int argc, _TCHAR* argv[])
{
	int n = 10;

	double* xEqual = new double[n];
	double* xChebyshev = new double[n];
	double* yEqual = new double[n];
	double* yChebyshev = new double[n];

	equalDistantNodes(xEqual, -1.0, 1.0, n);
	chebyshevNodes(xChebyshev, n); //tylko dla (-1,1)

	for (int i = 0; i <n; i++)
	{
		yEqual[i] = fx(xEqual[i]);
		yChebyshev[i] = fx(xChebyshev[i]);
	}

	double xp = -1.0;
	double xk = 1.0;
	double nValues = 50.0;
	double interval = (xk - xp) / nValues;

	cout << "Ilosc wezlow: " << n << endl;
	cout << setw(6) << "x" << setw(15) << "y" << setw(15) << "lg rownolegle" << setw(17) << "y-lgR" << setw(15) << "lg Czebyszewa" << setw(17) << "y-lgC" << endl;
	for (double x = xp+interval; x <=xk; x+=interval)
	{
		if (abs(x)< 1e-15)
			x = 0.0;
		double lgE = lagrangeInterpolation(xEqual, yEqual, x, n);
		double lgC = lagrangeInterpolation(xChebyshev, yChebyshev, x, n);
		double y = fx(x);
		cout << std::setprecision(6) << setw(6) << x << setw(15) << y << setw(15) << lgE << setw(17) << y - lgE << setw(15) << lgC << setw(17) << y - lgC << endl;
	}

	delete[] xEqual;
	delete[] xChebyshev;
	delete[] yEqual;
	delete[] yChebyshev;

	return 0;
}