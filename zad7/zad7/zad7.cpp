// zad7.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <cstdio>

using namespace std;

void createMatrix(double*** matrix, int n);
void deleteMatrix(double** matrix, int n);
void printVector(double* vector, int size);
void printMatrix(double** matrix, int n);
void loadData(double ** A, double* B, double* x);
void jacobi(double** A, double* B, double* x0, int n);
void gauss(double** A, double* B, double* x0, int n);
void calculateNextJacobiX(double* D, double* xn, double* B, double** LU, double* x0, int n);
void calculateNextGaussX(double** LD, double* xn, double* B, double** U, double* x0, int n);
void calculateNextSorX(double** LD, double* xn, double* B, double** LU, double* x0, int n);
void solveDiagonalSet(double* D, double* xn, double* B, int n);
void solveLowerTriangularSet(double** LD, double* xn, double* B, int n);
double normMax(double* vector1, double* vector2, int n);
void residuum(double** A, double* xn, double* B, double* xRes, int n);
double normMaxResiduum(double** A, double* xn, double* B, int n);


void createMatrix(double*** matrix, int n)
{
	int j;

	*matrix = new double* [n];

	for (j = 0; j<n; j++)
	{
		(*matrix)[j] = new double [n];
	}

	return;
}

void deleteMatrix(double** matrix, int n)
{
	for (int j = 0; j<n; j++)
	{
		delete[](matrix[j]);
	}

	delete[](matrix);
}

void printVector(double* vector, int size)
{
	int i = 0;

	for (i = 0; i < size; i++)
	{
		printf("%.8f ", vector[i]);
	}
}

void printMatrix(double** matrix, int n)
{
	printf("\n");
	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < n; i++)
		{
			printf("%lf ", matrix[j][i]);
		}
		printf("\n");
	}
	printf("\n");
}

void loadData(double** A,double* B, double* x)
{
	A[0][0] = 100.0; A[0][1] = -1.0;  A[0][2] =  2.0;  A[0][3] = -3.0;
	A[1][0] =  1.0;  A[1][1] = 200.0; A[1][2] = -4.0;  A[1][3] =  5.0;
	A[2][0] = -2.0;  A[2][1] =  4.0;  A[2][2] = 300.0; A[2][3] = -6.0;
	A[3][0] =  3.0;  A[3][1] = -5.0;  A[3][2] =  6.0;  A[3][3] = 400.0;

	B[0] =  116;
	B[1] = -226;
	B[2] =  912;
	B[3] = -1174;

	x[0] = 2;
	x[1] = 2;
	x[2] = 2;
	x[3] = 2;
}

void jacobi(double** A, double* B, double* X, int n)
{
	int nJacob = 50;
	double tolArg = 1e-8;
	double tolRes = 1e-10;
	double** LU;
	double* D = new double[n];
	double* xn = new double[n];
	double* x0 = new double[n];

	createMatrix(&LU, n);

	//rozklad A na (L + U) + D
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i == j)
			{
				D[i] = A[i][j];
				LU[i][j] = 0.0;
			}
			else
			{
				LU[i][j] = A[i][j];
			}
		}
	}

	//inicjalizacja xn na potrzeby algorytmu
	for (int i = 0; i < n; i++)
	{
		xn[i] = X[i];
	}

	int i = 0;

	cout << endl << "Metoda Jacobiego: " << endl << endl;

	do
	{
		for (int j = 0; j < n; j++)
		{
			x0[j] = xn[j];
		}

		calculateNextJacobiX(D, xn, B, LU, x0, n);

		cout << "iteracja: " << i << " wektor x: ";
		printVector(xn, n);
		cout << " normMax(Xn - Xn-1): " << normMax(xn, x0, n) << " normMaxRes(Ax-B): " << normMaxResiduum(A, xn, B, n) << endl;

		i++;

	} while ((i < nJacob) && (normMax(xn, x0, n)>tolArg || normMaxResiduum(A, xn, B, n)>tolRes));

	deleteMatrix(LU, n);
	delete[] D;
	delete[] xn;
	delete[] x0;

}

void gauss(double** A, double* B, double* X, int n)
{
	int nGauss = 50;
	double tolArg = 1e-8;
	double tolRes = 1e-10;
	double** LD;
	double** U;
	double* xn = new double[n];
	double* x0 = new double[n];

	createMatrix(&LD, n);
	createMatrix(&U, n);

	//rozklad A na (L + D) + U
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i < j)
			{
				LD[i][j] = 0.0;
				U[i][j] = A[i][j];
			}
			else
			{
				LD[i][j] = A[i][j];
				U[i][j] = 0.0;
			}
		}
	}

	//inicjalizacja xn na potrzeby algorytmu
	for (int i = 0; i < n; i++)
	{
		xn[i] = X[i];
	}

	int i = 0;

	cout << endl << "Metoda Gaussa-Seidela: " << endl << endl;

	do
	{
		for (int j = 0; j < n; j++)
		{
			x0[j] = xn[j];
		}

		calculateNextGaussX(LD, xn, B, U, x0, n);

		cout << "iteracja: " << i << " wektor x: ";
		printVector(xn, n);
		cout << " normMax(Xn - Xn-1): " << normMax(xn, x0, n) << " normMaxRes(Ax-B): " << normMaxResiduum(A, xn, B, n) << endl;

		i++;

	} while ((i < nGauss) && (normMax(xn, x0, n)>tolArg || normMaxResiduum(A, xn, B, n)>tolRes));

	deleteMatrix(LD, n);
	deleteMatrix(U, n);
	delete[] xn;
	delete[] x0;

}

void sor(double** A, double* B, double* X, int n)
{
	double omega = 2.0;
	int nSor = 50;
	double tolArg = 1e-8;
	double tolRes = 1e-10;
	double** LD;
	double** DU;
	double* xn = new double[n];
	double* x0 = new double[n];

	createMatrix(&LD, n);
	createMatrix(&DU, n);

	//rozklad A  dopisac...
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i < j)
			{
				LD[i][j] = 0.0;
				DU[i][j] = A[i][j];
			}
			else if (i == j)
			{
				LD[i][j] = A[i][j] * omega;
				DU[i][j] = A[i][j] * (1.0 - omega);
			}
			else
			{
				LD[i][j] = A[i][j];
				DU[i][j] = 0.0;
			}
		}
	}

	//inicjalizacja xn na potrzeby algorytmu
	for (int i = 0; i < n; i++)
	{
		xn[i] = X[i];
	}

	int i = 0;

	cout << endl << "Metoda Sukcesywnej relaksacji: " << endl << endl;

	do
	{
		for (int j = 0; j < n; j++)
		{
			x0[j] = xn[j];
		}

		calculateNextSorX(LD, xn, B, DU, x0, n);

		cout << "iteracja: " << i << " wektor x: ";
		printVector(xn, n);
		cout << " normMax(Xn - Xn-1): " << normMax(xn, x0, n) << " normMaxRes(Ax-B): " << normMaxResiduum(A, xn, B, n) << endl;

		i++;

	} while ((i < nSor) && (normMax(xn, x0, n)>tolArg || normMaxResiduum(A,xn,B,n)>tolRes));

	deleteMatrix(LD, n);
	deleteMatrix(DU, n);
	delete[] xn;
	delete[] x0;
}

void calculateNextJacobiX(double* D, double* xn, double* B, double** LU, double* x0, int n)
{
	double* tmpB = new double[n];
	double sum;

	//prawa strona B-(L+U)*xn-1
	for (int i = 0; i < n; i++)
	{
		sum = 0.0;

		for (int j = 0; j < n; j++)
		{
			sum += LU[i][j] * x0[j];
		}

		tmpB[i] = B[i] - sum;
	}

	solveDiagonalSet(D, xn, tmpB, n);

	delete[] tmpB;
}

void calculateNextGaussX(double** LD, double* xn, double* B, double** U, double* x0, int n)
{
	double* tmpB = new double[n];
	double sum;

	//prawa strona B-U*xn-1
	for (int i = 0; i < n; i++)
	{
		sum = 0.0;
		for (int j = 0; j < n; j++)
		{
			sum += U[i][j] * x0[j];
		}

		tmpB[i] = B[i] - sum;
	}

	solveLowerTriangularSet(LD, xn, tmpB, n);

	delete[] tmpB;
}

void calculateNextSorX(double** LD, double* xn, double* B, double** DU, double* x0, int n)
{
	double* tmpB = new double[n];
	double sum;

	//prawa strona B-U*xn-1
	for (int i = 0; i < n; i++)
	{
		sum = 0.0;
		for (int j = 0; j < n; j++)
		{
			sum += DU[i][j] * x0[j];
		}

		tmpB[i] = B[i] - sum;
	}

	solveLowerTriangularSet(LD, xn, tmpB, n);

	delete[] tmpB;
}

void solveDiagonalSet(double* D, double* xn, double* B, int n)
{
	for (int i = 0; i < n; i++)
	{
		xn[i] = B[i] / D[i];
	}
}

void solveLowerTriangularSet(double** LD, double* xn, double* B, int n)
{
	double sum;
	xn[0] = B[0] / LD[0][0];
	for (int i = 1; i < n; i++)
	{
		sum = 0.0;
		for (int j = 0; j < i; j++)
		{
			sum += LD[i][j] * xn[j];
		}
		xn[i] = (B[i] - sum) / LD[i][i];
	}
}

double normMax(double* vector1, double* vector2,int n)
{
	double max = abs(vector1[0]-vector2[0]);

	for (int i = 1; i < n; i++)
	{
		if (abs(vector1[i] - vector2[i]) > max)
			max = abs(vector1[i] - vector2[i]);
	}

	return max;
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

void residuum(double** A, double* xn, double* B,double* xRes, int n)
{
	double sum;

	for (int i = 0; i < n; i++)
	{
		sum = 0.0;

		for (int j = 0; j < n; j++)
		{
			sum += A[i][j] * xn[j];
		}

		xRes[i] = sum - B[i];
	}
}

double normMaxResiduum(double** A, double* xn, double* B, int n)
{
	double* xRes = new double[n];

	residuum(A, xn, B, xRes, n);

	double norm = normMax(xRes,n);

	delete[] xRes;

	return norm;
}

int _tmain(int argc, _TCHAR* argv[])
{
	int n = 4;

	double** A;
	double* x0 = new double[n];
	double* B = new double[n];

	createMatrix(&A, n);

	loadData(A,B,x0);

	jacobi(A, B, x0,n);
	gauss(A, B, x0, n);
	sor(A, B, x0, n);

	deleteMatrix(A, 4);

	delete[] x0;
	delete[] B;

	return 0;
}

