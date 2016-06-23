// zad5.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <cstdio>

using namespace std;

void printMatrix(double matrixA[4][4]);
void printVector(double vector[4]);
int findAbsMax(double matrix[4][4], int i, int j);
void swapRows(double matrix[4][4], int j1, int j2, double vector[4]);
void partialChoice(double matrix[4][4], int i, int j, double vector[4]);
void gaussianEliminate(double matrix[4][4], int i, int j);
void solveLowerTriangular(double matrix[4][4], double y[4], double b[4]);
void solveUpperTriangular(double matrix[4][4], double x[4], double b[4]);

void gaussianEliminate(double matrix[4][4], int ip, int jp)
{
	for (int i = ip + 1; i < 4; i++)
	{
		matrix[i][jp] = matrix[i][jp] / matrix[ip][jp];
		for (int j = jp + 1; j < 4; j++)
		{
			matrix[i][j] = matrix[i][j] - matrix[i][jp] * matrix[ip][j];
		}
	}
	
}

void partialChoice(double matrix[4][4], int i, int j, double vector[4])
{
	swapRows(matrix, findAbsMax(matrix, i, j), j, vector);
}

void swapRows(double matrix[4][4], int j1, int j2, double vector[4])
{
	double tmp;

	for (int i = 0; i < 4; i++)
	{
		tmp = matrix[j1][i];
		matrix[j1][i] = matrix[j2][i];
		matrix[j2][i] = tmp;
	}

	tmp = vector[j1];
	vector[j1] = vector[j2];
	vector[j2] = tmp;
}

int findAbsMax(double matrix[4][4], int i, int j)
{
	int k = 0;
	int indexMax = i;
	if (indexMax == 3)
	{
		return indexMax;
	}

	for (k = i+1; k < 4; k++)
	{
		if (abs(matrix[k][j]) > abs(matrix[indexMax][j]))
		{
			indexMax = k;
		}
	}



	return indexMax;
}

void solveLowerTriangular(double matrix[4][4], double y[4], double b[4])
{
	for (int i = 0; i < 4; i++)
	{
		double sum = 0.0;
		for (int j = 0; j <= i - 1; j++)
		{
			sum += matrix[i][j] * y[j];
		}
		y[i] = (b[i] - sum); // /1.0 (bo na przekatnej sa jedynki)
	}
}

void solveUpperTriangular(double matrix[4][4], double x[4], double b[4])
{
	for (int i = 3; i >= 0; i--)
	{
		double sum = 0.0;
		for (int j = i+1; j < 4; j++)
		{
			sum += matrix[i][j] * x[j];
		}
		x[i] = (b[i] - sum) / matrix[i][i];
	}
}

int _tmain(int argc, _TCHAR* argv[])
{
	double matrixA[4][4] = 
	{
		{ 1.0,   -20.0,   30.0,  -4.0},
		{ 2.0,   -40.0,  -6.0,    50.0},
		{ 9.0,   -180.0,  11.0,  -12.0},
		{-16.0,   15.0,  -140.0,  13.0}
	};

	double vectorB[4] = { 35, 104, -366, -354 };

	double y[4];
	double x[4];

	cout << "\nA Matrix\n\n";
	printMatrix(matrixA);

	for (int i = 0; i < 4; i++)
	{
		
		if (matrixA[i][i] == 0)
		{
			cout << "Partial Choice " << i << ", " << i << endl;
			partialChoice(matrixA, i, i, vectorB);
		}
		
		printMatrix(matrixA);

		cout << "Gausssian eliminate " << i << ", " << i << endl;
		gaussianEliminate(matrixA, i, i);
		printMatrix(matrixA);
	}

	cout << "wektor B po zamianach" << endl;
	printVector(vectorB);

	solveLowerTriangular(matrixA, y, vectorB);

	cout << "wektor y" << endl;
	printVector(y);

	solveUpperTriangular(matrixA, x, y);

	cout << "wektor x" << endl;
	printVector(x);


	return 0;
}

void printMatrix(double matrix[4][4])
{
	int i = 0, j = 0;

	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			printf("%f ", matrix[i][j]);
		}
		cout << "\n";
	}
	cout << endl << endl;
}

void printVector(double vector[4])
{
	int i = 0;

	for (i = 0; i < 4; i++)
	{
		printf("%f ", vector[i]);
	}

	cout << endl << endl;
}

