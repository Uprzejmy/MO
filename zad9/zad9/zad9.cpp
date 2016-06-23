// zad9.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#define _USE_MATH_DEFINES

#include <math.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include <locale>
#include <iomanip>

using namespace std;

//warunki zadania, wyliczone na kartce
double p = 1.0;
double q = 0.0;
double r = 1.0;
double s(double x){ return 2.0*sin(x); }
double alfa = 0.0;
double beta = 1.0;
double gamma = 0.0;
double fi = 0.0;
double psi = 1.0;
double teta = 0.0;

double xp = 0.0;
double xk = M_PI_2;

double numerovValues[100];
double conventionalValues[100];
double analyticValues[100];
double arguments[100];

void resultsToFile(string filename="values.csv")
{
	ofstream file(filename, std::ofstream::trunc);
	std::locale mylocale("");
	file.imbue(mylocale);

	for (int i = 0; i < 100; i++)
	{
		if (file.is_open())
		{
			file << arguments[i] << ";"
				 << analyticValues[i] << ";"
				 << conventionalValues[i] << ";"
				 << numerovValues[i] << endl;
		}
		else
		{
			cout << "Blad otwarcia pliku" << endl;
		}
	}

	file.close();
}

void createNodes(double* nodes, double h, double n)
{
	for (int i = 0; i < n; i++)
	{
		nodes[i] = h*i;
	}
}

double normMax(double* vector, int n)
{
	double max = abs(vector[0]);

	for (int i = 1; i < n; i++)
	{
		if (abs(vector[i]) > max)
			max = abs(vector[i]);
	}

	return max;
}

double analyticSolution(double x)
{
	return x*cos(x);
}

void transformDiagonals(double* diagU, double* diagL, double* diagD, double* vectorB, int n)
{
	for (int i = 1; i < n; i++)
	{
		diagD[i] = diagD[i] - diagL[i - 1] * diagU[i - 1] / diagD[i - 1];
		vectorB[i] = vectorB[i] - diagL[i - 1] * vectorB[i - 1] / diagD[i - 1];
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

void solveThomas(double* diagL, double* diagD, double* diagU, double* vectorB, double* y, int n)
{
	transformDiagonals(diagU, diagL, diagD, vectorB, n);
	solveSet(diagU, diagD, vectorB, y, n);
}

double numerovAproach(double h, double* nodes, int n)
{
	double *l, *d, *u, *b, *x, *error;

	l = new double[n];
	d = new double[n];
	u = new double[n];
	b = new double[n];
	x = new double[n];
	error = new double[n];

	u[0] = alfa / h;
	d[0] = beta - alfa / h;
	b[0] = -gamma;

	for (int i = 1; i < n - 1; i++)
	{

		l[i - 1] = p / (h * h) + r / 12.0;
		d[i] = (-2.0 * p) / (h * h) + r * 10.0 / 12.0;
		u[i] = p / (h * h) + r / 12.0;
		b[i] = -(s(nodes[i-1]) + 10.0*s(nodes[i]) + s(nodes[i+1])) / 12.0;
	}

	l[n - 2] = -fi / h;
	d[n - 1] = -fi / h + psi;
	b[n - 1] = -teta;

	solveThomas(l, d, u, b, x, n);

	if (n == 100)
	{
		for (int i = 0; i < n; i++)
		{
			numerovValues[i] = x[i];
			arguments[i] = nodes[i];
			analyticValues[i] = analyticSolution(arguments[i]);
		}
	}

	for (int i = 0; i < n; i++)
	{
		error[i] = abs(x[i] - analyticSolution(nodes[i]));
	}

	delete[] l;
	delete[] d;
	delete[] u;
	delete[] x;
	delete[] b;

	return normMax(error, n);
}

double threepointConventionalAproach(double h,double* nodes,int n)
{
	double *l, *d, *u, *b, *x, *error;

	l = new double[n];
	d = new double[n];
	u = new double[n];
	b = new double[n];
	x = new double[n];
	error = new double[n];

	u[0] = alfa / h;
	d[0] = beta - alfa / h;
	b[0] = -gamma;

	for (int i = 1; i < n - 1; i++)
	{
		l[i - 1] = p / (h * h) - q / (2.0 * h);
		d[i] = (-2.0 * p) / (h * h) + r;
		u[i] = p / (h * h) - q / (2.0 * h);
		b[i] = -s(nodes[i]);
	}

	l[n - 2] = -fi / h;
	d[n - 1] = -fi / h + psi;
	b[n - 1] = -teta;

	solveThomas(l, d, u, b, x, n);

	if (n == 100)
	{
		for (int i = 0; i < n; i++)
		{
			conventionalValues[i] = x[i];
		}
	}

	for (int i = 0; i < n; i++)
	{
		error[i] = abs(x[i] - analyticSolution(nodes[i]));	
	}

	delete[] l;
	delete[] d;
	delete[] u;
	delete[] x;
	delete[] b;

	return normMax(error, n);
}

void calculations()
{
	double h = 0.1;
	double* nodes;

	ofstream file("results.csv", std::ofstream::trunc);
	std::locale mylocale("");
	file.imbue(mylocale);

	int n = 100;
	//calculations
	while (h>9e-6)
	{
		nodes = new double[n];
		//h = (xk - xp) / (n - 1);
		h = xk / (n-1); //bo xp = 0
		createNodes(nodes, h, n);
		
		if (file.is_open())
		{
			cout << "h: " << h << endl;
			cout << "n: " << n << endl;
			file << std::setprecision(8) << log10(h) << ";";
			file << std::setprecision(8) << log10(threepointConventionalAproach(h,nodes,n)) << ";";
			file << std::setprecision(8) << log10(numerovAproach(h,nodes,n)) << ";";
			file << endl;
		}
		else
		{
			cout << "Blad otwarcia pliku" << endl;
			file.close();
			return;
		}

		if (h > 9e-4)
		{
			n += 50;
		}
		else
		{
			n += 500;
		}
		delete[] nodes;
	}

	file.close();
}

int _tmain(int argc, _TCHAR* argv[])
{
	calculations();

	resultsToFile();
	
	return 0;
}

