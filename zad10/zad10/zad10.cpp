// zad10.cpp : Defines the entry point for the console application.
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

double analyticFormula(double t)
{
	return 1.0 - exp(-10.0 * (t + atan(t)));
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

void createNodes(double* nodes, double h, double n)
{
	for (int i = 0; i < n; i++)
	{
		nodes[i] = h*i;
	}
}

double errorDirectEulerMethod(double* nodes, double h, int n)
{
	double normMax = 0.0;
	double error;

	double y = 0.0;

	for (int i = 1; i < n; i++)
	{
		y = y - h*(10.0*nodes[i] * nodes[i] + 20.0) / (nodes[i] * nodes[i] + 1)*(y - 1);

		error = abs(analyticFormula(nodes[i]) - y);
		if (error > normMax)
			normMax = error;
	}

	return normMax;
}

double errorTrapeziumMethod(double* nodes, double h, int n)
{
	double normMax = 0.0;
	double error;

	double y = 0.0;

	double tmpK0;
	double tmpK1;

	for (int i = 1; i < n; i++)
	{
		tmpK0 = (10.0* nodes[i] * nodes[i] + 20.0) / (nodes[i] * nodes[i] + 1.0)*h;
		tmpK1 = (10.0* nodes[i + 1] * nodes[i + 1] + 20.0) / (nodes[i + 1] * nodes[i + 1] + 1.0)*h;
		y = (2.0* y - tmpK0*(y - 1.0) + tmpK1) / (2.0 + tmpK1);

		error = abs(analyticFormula(nodes[i]) - y);
		if (error > normMax)
			normMax = error;
	}

	return normMax;
}

double errorIndirectEulerMethod(double* nodes, double h, int n)
{
	double normMax = 0.0;
	double error;

	double y = 0.0;

	double tmp;

	for (int i = 1; i < n; i++)
	{
		tmp = (10.0* nodes[i] * nodes[i] + 20.0) / (nodes[i] * nodes[i] + 1.0)*h;
		y = (y + tmp) / (1.0 + tmp);

		error = abs(analyticFormula(nodes[i]) - y);
		if (error > normMax)
			normMax = error;
	}

	return normMax;
}

void directEulerMethod(double* y, double* nodes, double h, int n)
{
	y[0] = 0.0;

	for (int i = 1; i < n; i++)
	{
		y[i] = y[i - 1] - h*(10.0*nodes[i] * nodes[i] + 20.0) / (nodes[i]*nodes[i] + 1)*(y[i - 1] - 1);
	}
}

void indirectEulerMethod(double* y, double* nodes, double h, int n)
{
	y[0] = 0.0;

	double tmp;

	for (int i = 1; i < n; i++)
	{
		tmp = (10.0* nodes[i] * nodes[i] + 20.0) / (nodes[i] * nodes[i] + 1.0)*h;
		y[i] = (y[i - 1] + tmp) / (1.0 + tmp);
	}
}

void trapeziumMethod(double* y, double* nodes, double h, int n)
{
	y[0] = 0.0;

	double tmpK0;
	double tmpK1;

	for (int i = 1; i < n; i++)
	{
		tmpK0 = (10.0* nodes[i] * nodes[i] + 20.0) / (nodes[i] * nodes[i] + 1.0)*h;
		tmpK1 = (10.0* nodes[i + 1] * nodes[i + 1] + 20.0) / (nodes[i + 1] * nodes[i + 1] + 1.0)*h;
		y[i] = (2.0* y[i - 1] - tmpK0*(y[i-1]-1.0)+tmpK1) / (2.0 + tmpK1);
	}
}

void chartsIndirectAndTrapezium()
{
	ofstream file("stabilnePosrTrap.csv", std::ofstream::trunc);
	std::locale mylocale("");
	file.imbue(mylocale);

	double* yIndirectEuler;
	double* yTrapezium;

	double h = 0.01;
	double* nodes;
	double tp = 0.0, tk = 3.0;
	int n = (tk - tp) / h;

	yIndirectEuler = new double[n];
	yTrapezium = new double[n];
	nodes = new double[n];

	double t = tp;
	for (int i = 0; i < n; i++)
	{
		nodes[i] = tp + i*h;
	}

	indirectEulerMethod(yIndirectEuler, nodes, h, n);
	trapeziumMethod(yTrapezium, nodes, h, n);

	/*
	cout.width(4);
	cout << "i: " << "|";
	cout.width(6);
	cout << "t: " << "|";
	cout.width(25);
	cout << "rozwiazanie analityczne" << "|";
	cout.width(25);
	cout << "metoda posrednia eulera" << "|";
	cout.width(25);
	cout << "metoda trapezow" << "|" << endl;
	*/

	for (int i = 0; i < n; i++)
	{
		if (file.is_open())
		{
			file << nodes[i] << ";";
			file << analyticFormula(nodes[i]) << ";";
			file << yIndirectEuler[i] << ";";
			file << yTrapezium[i] << ";" << endl;
		}
		else
		{
			cout << "Blad otwarcia pliku" << endl;
		}
	}

	file.close();

	delete[] yIndirectEuler;
	delete[] yTrapezium;
	delete[] nodes;
}

void chartDirect()
{
	ofstream fileBS("bezposrStab.csv", std::ofstream::trunc);
	std::locale mylocale("");
	fileBS.imbue(mylocale);

	double* yDirectEuler;

	double h = 0.01;
	double* nodes;
	double tp = 0.0, tk = 3.0;
	int n = (tk - tp) / h;

	yDirectEuler = new double[n];
	nodes = new double[n];

	for (int i = 0; i < n; i++)
	{
		nodes[i] = tp + i*h;
	}

	directEulerMethod(yDirectEuler, nodes, h, n);

	/*
	cout.width(4);
	cout << "i: " << "|";
	cout.width(6);
	cout << "t: " << "|";
	cout.width(25);
	cout << "rozwiazanie analityczne" << "|";
	cout.width(25);
	cout << "metoda posrednia eulera" << "|";
	cout.width(25);
	cout << "metoda trapezow" << "|" << endl;
	*/

	for (int i = 0; i < n; i++)
	{
		if (fileBS.is_open())
		{
			fileBS << nodes[i] << ";";
			fileBS << analyticFormula(nodes[i]) << ";";
			fileBS << yDirectEuler[i] << ";" << endl;
		}
		else
		{
			cout << "Blad otwarcia pliku" << endl;
		}
	}

	fileBS.close();

	delete[] yDirectEuler;
	delete[] nodes;

	//niestabilnosc

	ofstream fileBNS("bezposrNStab.csv", std::ofstream::trunc);
	fileBNS.imbue(mylocale);


	h = 0.3;
	tp = 0.0, tk = 4.0;
	n = (tk - tp) / h;

	yDirectEuler = new double[n];
	nodes = new double[n];

	for (int i = 0; i < n; i++)
	{
		nodes[i] = tp + i*h;
	}

	directEulerMethod(yDirectEuler, nodes, h, n);

	/*
	cout.width(4);
	cout << "i: " << "|";
	cout.width(6);
	cout << "t: " << "|";
	cout.width(25);
	cout << "rozwiazanie analityczne" << "|";
	cout.width(25);
	cout << "metoda posrednia eulera" << "|";
	cout.width(25);
	cout << "metoda trapezow" << "|" << endl;
	*/

	for (int i = 0; i < n; i++)
	{
		if (fileBNS.is_open())
		{
			fileBNS << nodes[i] << ";";
			fileBNS << analyticFormula(nodes[i]) << ";";
			fileBNS << yDirectEuler[i] << ";" << endl;
		}
		else
		{
			cout << "Blad otwarcia pliku" << endl;
		}
	}

	fileBNS.close();

	delete[] yDirectEuler;
	delete[] nodes;
}

void errors()
{
	ofstream file("errors.csv", std::ofstream::trunc);
	std::locale mylocale("");
	file.imbue(mylocale);

	double h = 0.1;
	int n;
	double tp = 0.0, tk = 1.5;

	double* nodes;
	
	//calculations
	while (h>1e-18)
	{
		n = (tk - tp) / h;
		if (n > 10000 || n<0)
			n = 10000;

		nodes = new double[n];
		double tmp = tp;
		for (int i = 0; i < n; i++)
		{
			nodes[i] = tmp;
			tmp += h;
		}

		if (file.is_open())
		{
			cout << "h: " << h << endl;
			//cout << "n: " << n << endl;
			file << std::setprecision(8) << log10(h) << ";";
			file << std::setprecision(8) << log10(errorDirectEulerMethod(nodes,h, n)) << ";";
			file << std::setprecision(8) << log10(errorIndirectEulerMethod(nodes, h, n)) << ";";
			file << std::setprecision(8) << log10(errorTrapeziumMethod(nodes, h, n)) << ";";
			file << endl;
		}
		else
		{
			cout << "Blad otwarcia pliku" << endl;
			file.close();
			return;
		}

		h /= 2.0;
		delete[] nodes;
	}

	file.close();
}

int _tmain(int argc, _TCHAR* argv[])
{
	chartsIndirectAndTrapezium();
	chartDirect();
	errors();
	
	return 0;

}

