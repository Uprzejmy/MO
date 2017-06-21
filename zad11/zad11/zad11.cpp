// zad11.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <math.h>
#include <fstream>
#include <locale>
#include <iomanip>
#include "calerf.h"

using namespace std;

//macierz trójdiagonalna jako wektory
double* L;
double* U;
double* D;

double* prevResult;
double* result;
double* analyticResult;

double DFactor = 1.0; //wspó³czynnik dyfuzji
double b = 10;//6.0*sqrt(2.0); //koniec przedzia³u dla zmiennej przestrzennej

double getMatrix(int i, int j)
{
	if (i == j - 1) 
	{
		return U[i];
	}
	else if (i == j) 
	{
		return D[i];
	}
	else if (i == j + 1) 
	{
		return L[j];
	}
	else 
	{
		return 0.0;
	}
}

double analyticSolution(double x, double t)
{
	return erf(x / 2 / sqrt(t));//sqrt(t)*D ale D = 1.0
}

void solveThomas(double* prevResult, double* result, int n)
{
	double* vectU = new double[n]; // eta
	double* vectB = new double[n]; // theta
	vectU[0] = D[0];
	for (int i = 1; i < n; ++i)
		vectU[i] = D[i] - (L[i - 1] * U[i - 1] / vectU[i-1]);

	vectB[0] = prevResult[0];
	for (int i = 1; i < n; ++i)
		vectB[i] = prevResult[i] - (L[i-1] * vectB[i - 1]) / vectU[i-1]; // nowe b

	result[n - 1] = vectB[n - 1]/vectU[n-1];
	for (int i = n - 2; i >= 0; --i)
		result[i] = (vectB[i] - U[i] * result[i + 1]) / vectU[i]; // propagacja wsteczna - rozwiazania

	delete[] vectU;
	delete[] vectB;
}

void laasonenDiscretisation(double h, double dt)
{
	ofstream file("Dane.csv", std::ofstream::trunc);
	std::locale mylocale("");
	file.imbue(mylocale);

	double tMax = 2.0;
	double lambda = dt / h / h;//raczej zawsze 1, chyba ze beda baddzo male kroki, wtedy moze byc blad
	int n = b / h; // n - ilosc wêz³ów zmiennej przestrzennej             //moze h-1, zobaczyc
	int k = tMax / dt; // k - ilosc wezlow zmiennej czasowej

	prevResult = new double[n];
	result = new double[n];
	analyticResult = new double[n];

	L = new double[n - 1];
	U = new double[n - 1];
	D = new double[n];

	double* xNodes = new double[n+1];
	double* tNodes = new double[k+1];

	//wyznaczam siatke przestrzenna, zeby nie powtarzac obliczen
	for (int i = 0; i < n+1; i++)
	{
		xNodes[i] = h*i;
	}

	//wyznaczam siatke czasowa, zeby nie powtarzac obliczen
	for (int j = 0; j < k+1; j++)
	{
		tNodes[j] = dt*j;
	}




	//warunek poczatkowy
	for (int i = 0; i < n; i++)
	{
		prevResult[i] = -1.0;
	}

	//ustawienie macierzy i wyrazu wolnego
	prevResult[0] = 0.0;
	prevResult[n - 1] = 1.0;

	D[0] = 1.0;
	U[0] = 0.0;
	//uzupelnianie macierzy A   gosc ma inaczej
	for (int i = 1; i < n-1; i++) 
	{
		U[i] = lambda;
		L[i-1] = lambda;
		D[i] = -(1.0 + 2.0*lambda);
	}

	//i jeszcze ostatni wiersz
	L[n - 2] = 0.0;
	D[n - 1] = 1.0;

	// wrzucam pierwszy wezel czasowy do pliku
	if (file.is_open())
	{
		file << tNodes[0];
	}
	else
	{
		cout << "Blad otwarcia pliku" << endl;
		return;
	}
	
	//wrzucam siatke przestrzenna do pliku
	for (int i = 0; i < n; i++)
	{
		file << xNodes[i] << ";";
	}
	file << endl;

	//wrzucam rozwiazania w chwili t=0
	for (int i = 0; i < n; i++)
	{
		file << prevResult[i] << ";";
	}

	file << endl;

	for (int j = 1; j<k; j++)
	{
		result[0] = 0.0;
		result[n - 1] = 1.0; // zgodnie z warunkami brzegowymi
		//t = j*dtt;

		prevResult[0] = 0.0;
		prevResult[n - 1] = 1.0; // zgodnie z warunkami brzegowymi

		solveThomas(prevResult,result,n);

		result[0] = 0.0;
		result[n - 1] = 1.0; // zgodnie z warunkami brzegowymi

		file << tNodes[j];

		for (int i = 0; i < n; i++)
		{
			file << result[i] << ";";

			prevResult[i] = -result[i];
		}

		file << endl;


		//cout << "test j: " << j << " k: " << k << endl;
	}

	delete[] xNodes;
	delete[] tNodes;

	delete[] prevResult;
	delete[] result;
	delete[] analyticResult;

	delete[] L;
	delete[] U;
	delete[] D;

	file.close();
}




int _tmain(int argc, _TCHAR* argv[])
{
	double h = 3e-2;
	double dt = h*h;

	laasonenDiscretisation(h, dt);

	return 0;
}

