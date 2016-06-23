// lab1.cpp: Okreœla punkt wejœcia dla aplikacji konsoli.
//

#include "stdafx.h"
#include <iostream>

using namespace std;

void floatFunction()
{
	float x = 1.0f;
	int n = 0;

	while (1.0f + x > 1.0f)
	{
		x = x / 2.0f;
		n++;
	}

	x *= 2.0f;
	n--;

	cout << "Liczba bitow mantysy dla zmiennej typu float: " << n << endl;
	printf("%.10e\n", x);
}

void doubleFunction()
{
	double x = 1.0;
	int n = 0;

	while (1.0 + x > 1.0)
	{
		x = x / 2.0;
		n++;
	}

	x *= 2.0;
	n--;

	printf("%.10e\n", x);
	cout << "Liczba bitow mantysy dla zmiennej typu double: " << n << endl;
}

int main()
{
	floatFunction();
	doubleFunction();

	return 0;
}

