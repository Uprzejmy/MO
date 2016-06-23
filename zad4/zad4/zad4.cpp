// zad4.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#define _USE_MATH_DEFINES
#include "math.h"

double Xn[3];
double Xn1[3];
double fx[3];
double j[3][3];
double fd[3][3]; 
double fdo[3][3];
double tmpVector[3];
double detJ;
int nNewton = 50;
double eps1 = 1e-8;
double eps2 = 1e-11;

void fxValues(double x, double y, double z);
void jValues(double x, double y, double z);
void detJValue(double x, double y, double z);
void fdValues(double x, double y, double z);
void fdoValues(double x, double y, double z);
void tmpVectorValues(double x, double y, double z);
void newton();
void calculateEpsilon();
void printVector(double vector[3]);
void printMatrix(double matrix[3][3]);
void printDouble(double x);
double normMax(double vector[3]);
double norm2(double vector[3]);


int _tmain(int argc, _TCHAR* argv[])
{
	newton();

	return 0;
}

void newton()
{
	int i = 0;

	double vectorMax[3];

	Xn1[0] = 1.0; //pierwsze punkty do przybli¿eñ
	Xn1[1] = 1.0;
	Xn1[2] = 1.0;

	do
	{
		Xn[0] = Xn1[0]; //pierwsze punkty do przybli¿eñ
		Xn[1] = Xn1[1];
		Xn[2] = Xn1[2];

		fxValues(Xn[0], Xn[1], Xn[2]);
		//printVector(fx);

		jValues(Xn[0], Xn[1], Xn[2]);
		//printMatrix(j);

		detJValue(Xn[0], Xn[1], Xn[2]);
		//printDouble(detJ);

		fdValues(Xn[0], Xn[1], Xn[2]);
		//printMatrix(fd);

		fdoValues(Xn[0], Xn[1], Xn[2]);

		tmpVectorValues(Xn[0], Xn[1], Xn[2]);

		Xn1[0] = Xn[0] - tmpVector[0];
		Xn1[1] = Xn[1] - tmpVector[1];
		Xn1[2] = Xn[2] - tmpVector[2];

		i++;

		//printf("i: %3d | xn-1: %.10e   | yn-1: %.10e   | zn-1: %.10e   | xn: %.10e   | yn: %.10e   | zn: %.10e    \n", i, Xn[0], Xn[1], Xn[2], Xn1[0], Xn1[1], Xn1[2]);
		printf("i: %3d | xn-1: %lf   | yn-1: %lf   | zn-1: %lf   | xn: %lf   | yn: %lf   | zn: %lf    \n", i, Xn[0], Xn[1], Xn[2], Xn1[0], Xn1[1], Xn1[2]);
		
		vectorMax[0] = Xn1[0] - Xn[0];
		vectorMax[1] = Xn1[1] - Xn[1];
		vectorMax[2] = Xn1[2] - Xn[2];

	} while ((i < nNewton) && (normMax(fx)>eps2) && (normMax(vectorMax)>eps1));

}

void fxValues(double x, double y, double z)
{
	fx[0] = x*y - 2.0;
	fx[1] = y / 2.0 - sin(M_PI / 4.0 - z);
	fx[2] = x*x + y*y + z*z - 4;
}

void jValues(double x, double y, double z)
{
	j[0][0] = y;
	j[0][1] = x;
	j[0][2] = 0;
	j[1][0] = 0;
	j[1][1] = 0.5;
	j[1][2] = cos(M_PI/4.0 - z);
	j[2][0] = 2.0*x;
	j[2][1] = 2.0*y;
	j[2][2] = 2.0*z;
}

void detJValue(double x, double y, double z)
{
	detJ = y*z + cos(M_PI / 4.0 - z)*(2.0 * x*x - 2.0 * y*y);
}

void fdValues(double x, double y, double z)
{
	fd[0][0] = j[1][1] * j[2][2] - j[1][2] * j[2][1];
	fd[0][1] = -(j[1][0]*j[2][2] - j[1][2]*j[2][0]);
	fd[0][2] = j[1][0] * j[2][1] - j[1][1] * j[2][0];
	fd[1][0] = -(j[0][1]*j[2][2] - j[2][1]*j[0][2]);
	fd[1][1] = j[0][0] * j[2][2] - j[2][0] * j[0][2];
	fd[1][2] = -(j[0][0]*j[2][1] - j[2][0]*j[0][1]);
	fd[2][0] = j[0][1] * j[1][2] - j[1][1] * j[0][2];
	fd[2][1] = -(j[0][0]*j[1][2] - j[1][0]*j[0][2]);
	fd[2][2] = j[0][0] * j[1][1] - j[1][0] * j[0][1];
}

void fdoValues(double x, double y, double z)
{
	fdo[0][0] = fd[0][0];
	fdo[0][1] = fd[1][0];
	fdo[0][2] = fd[2][0];
	fdo[1][0] = fd[0][1];
	fdo[1][1] = fd[1][1];
	fdo[1][2] = fd[2][1];
	fdo[2][0] = fd[0][2];
	fdo[2][1] = fd[1][2];
	fdo[2][2] = fd[2][2];
}

void tmpVectorValues(double x, double y, double z)
{
	tmpVector[0] = (fdo[0][0] * fx[0] + fdo[0][1] * fx[1] + fdo[0][2] * fx[2]) / detJ;
	tmpVector[1] = (fdo[1][0] * fx[0] + fdo[1][1] * fx[1] + fdo[1][2] * fx[2]) / detJ;
	tmpVector[2] = (fdo[2][0] * fx[0] + fdo[2][1] * fx[1] + fdo[2][2] * fx[2]) / detJ;
}

double normMax(double vector[3])
{
	double max = abs(vector[0]);

	if (abs(vector[1]) > max)
		max = abs(vector[1]);
	
	if (abs(vector[2]) > max)
		max = abs(vector[2]);

	return max;
}

double norm2(double vector[3])
{
	return sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
}

void printMatrix(double matrix[3][3])
{
	printf("\n");
	for (int j = 0; j < 3; j++)
	{
		for (int i = 0; i < 3; i++)
		{
			printf("%lf ", matrix[j][i]);
		}
		printf("\n");
	}
	printf("\n");
}

void printVector(double vector[3])
{
	printf("\n");
	for (int i = 0; i < 3; i++)
	{
		printf("%lf ", vector[i]);
	}
	printf("\n");
	printf("\n");
}

void printDouble(double x)
{
	printf("\n");
	printf("%lf \n", x);
	printf("\n");
}



