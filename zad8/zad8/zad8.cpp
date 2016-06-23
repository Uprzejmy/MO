// zad8.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#define _USE_MATH_DEFINES

#include <math.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include <locale>


using namespace std;

template <class type1> type1 fx(type1 x)
{
	return cos(x);
}

template <class type1> type1 derivative(type1 x)
{
	return -sin(x);
}

template <class type1> type1 forwardDifferenceBipoint(type1 x1, type1 x2, type1 h)
{
	return (fx(x2) - fx(x1)) / h;
}

template <class type1> type1 forwardDifferenceTripoint(type1 x0, type1 x1, type1 x2, type1 h)
{
	return (-fx(x2) + static_cast <type1>(4.0) * fx(x1) - static_cast <type1>(3.0) * fx(x0)) / static_cast<type1>(2.0) / h;
}

template <class type1> type1 backwardDifferenceBipoint(type1 x0, type1 x1, type1 h)
{
	return (fx(x1) - fx(x0)) / h;
}

template <class type1> type1 backwardDifferenceTripoint(type1 x0, type1 x1, type1 x2, type1 h)
{
	return (static_cast <type1>(3.0) * fx(x2) - static_cast <type1>(4.0) * fx(x1) + fx(x0)) / static_cast<type1>(2.0) / h;
}

template <class type1> type1 centralDifferenceBipoint(type1 x0, type1 x2, type1 h)
{
	return (fx(x2) - fx(x0)) / static_cast<type1>(2.0) / h;
}

template <class type1>type1 calculateEpsilon()
{
	type1 epsilon = static_cast <type1> (1.0);
	type1 one = static_cast <type1> (1.0);
	type1 two = static_cast <type1> (2.0);

	while (one + epsilon > one)
	{
		epsilon /= two;
	}

	return epsilon*two;
}

template <class type1>int countIterations()
{
	int i = 1;
	type1 epsilon = calculateEpsilon<type1>();
	type1 h = static_cast <type1> (0.1);
	type1 ten = static_cast <type1> (10.0);

	while (h > epsilon)
	{
		++i;
		h /= ten;
	}

	return i;
}

template <class type1> 
void matrixToCsv(type1** matrix, int size,string filename)
{
	ofstream file(filename, std::ofstream::trunc);
	std::locale mylocale("");
	file.imbue(mylocale);

	for (int i = 0; i < size; i++)
	{
		if (file.is_open())
		{
			for (int j = 0; j < 10; j++)
			{
				file << matrix[i][j] << ";";
				cout << matrix[i][j] << ";";
			}
		}

		file << endl;
		cout << endl;
		
	}

	file.close();
}

template <class type1> void calculations(string filename)
{
	type1 a = static_cast <type1> (1.0);
	type1 ab = static_cast <type1> (M_PI_4);
	type1 b = static_cast <type1> (M_PI_2);
	type1 h = static_cast <type1> (0.1);
	type1 ten = static_cast <type1> (10.0);
	int size = countIterations<type1>();
	
	//allocate space for a result
	type1** resultMatrix = new type1*[size];
	for (int i = 0; i < size; i++)
	{
		resultMatrix[i] = new type1[10];
	}

	//calculations
	for (int i = 0; i < size; i++)
	{
		//step
		resultMatrix[i][0] = log10(h);

		//left
		resultMatrix[i][1] = log10(abs(derivative(a) - forwardDifferenceBipoint(a, a+h,h)));
		resultMatrix[i][2] = log10(abs(derivative(a) - forwardDifferenceTripoint(a, a+h, a+h+h,h)));

		//center
		resultMatrix[i][3] = log10(abs(derivative(ab) - forwardDifferenceBipoint(ab, ab + h,h)));
		resultMatrix[i][4] = log10(abs(derivative(ab) - forwardDifferenceTripoint(ab - h, ab, ab + h,h)));
		resultMatrix[i][5] = log10(abs(derivative(ab) - centralDifferenceBipoint(ab-h, ab+h,h)));
		resultMatrix[i][6] = log10(abs(derivative(ab) - backwardDifferenceBipoint(ab - h, ab,h)));
		resultMatrix[i][7] = log10(abs(derivative(ab) - backwardDifferenceTripoint(ab-h, ab, ab+h,h)));

		//right
		resultMatrix[i][8] = log10(abs(derivative(b) - backwardDifferenceBipoint(b-h, b,h)));
		resultMatrix[i][9] = log10(abs(derivative(b) - backwardDifferenceTripoint(b-h-h, b-h, b,h)));

		h /= ten;
	}

	matrixToCsv(resultMatrix, size, filename);

	//deallocate space for a result
	for (int i = 0; i < size; i++)
	{
		delete[] resultMatrix[i];
	}

	delete[] resultMatrix;
}


int _tmain(int argc, _TCHAR* argv[])
{
	calculations<double>("double.csv");
	calculations<float>("float.csv");

	return 0;
}

