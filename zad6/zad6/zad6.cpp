// zad6.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <cstdio>

using namespace std;

void printVector(double* vector, int size);

void transformDiagonals(double* diagU, double* diagL, double* diagD, double* vectorB, int n)
{
	for (int i = 1; i < n; i++)
	{
		diagD[i] = diagD[i] - diagL[i-1] * diagU[i - 1] / diagD[i - 1];
		vectorB[i] = vectorB[i] - diagL[i-1] * vectorB[i - 1] / diagD[i - 1];
	}
}

void solveSet(double* diagU, double* diagD, double* vectorB, double* x, int n)
{
	x[n - 1] = vectorB[n - 1] / diagD[n - 1];
	
	for (int i = n - 2; i >= 0; i--)
	{
		x[i] = (vectorB[i] - diagU[i] * x[i + 1]) / diagD[i];
	}
}

void loadData(double* diagU, double* diagL, double* diagD, double* vectorB)
{
	int i = 0;
	diagU[i++] = 1.0 / 2.0;
	diagU[i++] = 1.0 / 4.0;
	diagU[i++] = 1.0 / 6.0;
	diagU[i++] = 1.0 / 8.0;
	diagU[i++] = 1.0 / 10.0;

	i = 0;
	diagD[i++] = 10.0;
	diagD[i++] = 20.0;
	diagD[i++] = 30.0;
	diagD[i++] = 30.0;
	diagD[i++] = 20.0;
	diagD[i++] = 10.0;

	i = 0;
	diagL[i++] = 1.0 / 3.0;
	diagL[i++] = 1.0 / 5.0;
	diagL[i++] = 1.0 / 7.0;
	diagL[i++] = 1.0 / 9.0;
	diagL[i++] = 1.0 / 11.0;

	i = 0;
	vectorB[i++] = 31.0;
	vectorB[i++] = 165.0 / 4.0;
	vectorB[i++] = 917.0 / 30.0;
	vectorB[i++] = 851.0 / 28.0;
	vectorB[i++] = 3637.0 / 90.0;
	vectorB[i++] = 332.0 / 11.0;
}

int _tmain(int argc, _TCHAR* argv[])
{
	int n = 6;

	double* diagU = new double[n - 1];
	double* diagL = new double[n - 1];
	double* diagD = new double[n];
	double* vectorB = new double[n];
	double* x = new double[n];

	loadData(diagU, diagL, diagD, vectorB);

	transformDiagonals(diagU, diagL, diagD, vectorB, n);
	solveSet(diagU, diagD, vectorB, x, n);

	cout << "wektor x" << endl;
	printVector(x,n);
	delete[] diagU;
	delete[] diagL;
	delete[] diagD;
	delete[] vectorB;
	delete[] x;

	return 0;
}

void printVector(double* vector, int size)
{
	int i = 0;

	for (i = 0; i < size; i++)
	{
		printf("%.2f \n", vector[i]);
	}

	cout << endl << endl;
}
