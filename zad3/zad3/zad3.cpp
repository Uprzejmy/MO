// zad3.cpp : Defines the entry point for the console application.
//




#include "stdafx.h"

#define _USE_MATH_DEFINES

#include "math.h"



int nPicard = 50;
int nBisection = 100;
int nNewton = 50;
int nSecant = 50;
double epsM;

void picard();
void bisection();
void calculateEpsilonM();
void newton();
void secant();



int _tmain(int argc, _TCHAR* argv[])
{
	calculateEpsilonM();
	picard();
	bisection();
	newton();
	secant();

	return 0;
}

double function1(double x)
{
	//return pow(cos(x / 4.0), 2) - x;
	return exp(x) - exp(-x) + x - 1;
}

double function1Fi(double x)
{
	//return pow(cos(x / 4.0), 2);
	return exp(-x) - exp(x) + 1;
}

double function1d(double x)
{
	//return cos(x / 2.0) / (-4.0) - 1;
	return exp(-x) + exp(x) + 1;
}

void picard()
{
	printf("\n\n-------------------- Picard --------------------\n");

	int i=0;
	double x0 = 0.5; //nieznacznie przed pi
	double fix0 = 0.5; //przypisanie na potrzeby algorytmu

	do
	{
		x0 = fix0;
		fix0 = function1Fi(x0);
		i++;

		printf("i: %3d | xn: %.13e   |   est(xn): %.13e   |   f(xn): %.13e \n", i, x0, abs(fix0-x0), fix0);

	} while ((i < nPicard) && (abs(fix0) > epsM) && (abs(fix0 - x0) > epsM));

}

void bisection()
{
	printf("\n\n-------------------- Bisection --------------------\n");

	int i = 0;
	double a = 0.2;
	double b = 4; 
	double c;
	double fc,fa,fb;

	if (function1(a)*function1(b) >= 0)
	{
		printf("zle dobrane konce przedzialow, funkcja musi mieæ w (a,b) przecinac ox nieparzyst¹ ilosc razy");
		return;
	}

	do
	{
		c = (a + b) / 2.0;
		fa = function1(a);
		fb = function1(b);
		fc = function1(c);
		if (((fa>0) && (fc>0)) || ((fa<0) && (fc<0)))
		{
			a = c;
		}
		else if (((fb>0) && (fc>0)) || ((fb<0) && (fc<0)))
		{
			b = c;
		}
		else
		{
			break;
		}
		i++;

		printf("i: %3d | xn: %.10e   |   est(xn): %.10e   |   f(xn): %.10e \n", i, c, abs((b-a)/2.0), fc);

	} while ((i < nBisection) && (abs(fc) > epsM) && (abs((b - a) / 2.0) > epsM));
}

void newton()
{
	printf("\n\n-------------------- Newton --------------------\n");

	int i = 0;
	double x0 = 0.5;
	double fx0 = 0.5;

	do
	{
		x0 = fx0;
		fx0 = x0 - (function1(x0) / function1d(x0));
		i++;

		printf("i: %3d | xn: %.10e   |   est(xn): %.10e   |   f(xn): %.10e \n", i, x0, abs(fx0-x0), fx0);

	} while ((i < nNewton) && (abs(fx0) > epsM) && (abs(fx0 - x0) > epsM));
}

void secant()
{
	printf("\n\n-------------------- Secant --------------------\n");

	int i = 0;
	double x0 = 0.4;
	double x1 = 0.5;
	double x2 = 0.6;

	do
	{
		x0 = x1;
		x1 = x2;
		x2 = x1 - (function1(x1) / (((function1(x1) - function1(x0))/(x1-x0))));
		i++;

		printf("i: %3d | xn: %.10e   |   est(x2): %.10e   |   f(xn): %.10e \n", i, x0, abs(x2-x1), x2);

	} while ((i < nSecant) && (abs(x2) > epsM) && (abs(x2 - x1) > epsM));
}

void calculateEpsilonM()
{
		double x = 1.0;

		while (1.0 + x > 1.0)
		{
			x = x / 2.0;
		}

		epsM = x * 2.0;
}