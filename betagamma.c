//Funkcje potrzebne do obliczenia p-value. Nie rozwodze siê zanadto z tym co tu sie dzieje, bo s¹ wykorzystywane tylko przez inne moje funkcje.


#include <math.h>
#include "betagamma.h"

#define ACCURACY 15 // dla x64: musi byæ wartoœæ 11, a funcke z math bez "l", np. sqrtl -> sqrt


double gamma(double val) // funkcja gamma liczona aproksymacj¹ Spouge'a
{

	long double z = (long double)val;
	long double sc = powl((z + ACCURACY), (z + 0.5));
	sc *= expl(-1.0 * (z + ACCURACY));
	sc /= z;

	long double fract = 1.0;
	long double ck;
	long double sum = atanl(1.0) * 8.0;


	for (int i = 1; i < ACCURACY; i++)
	{
		z++;
		ck = powl(ACCURACY - i, i - 0.5);
		ck *= expl(ACCURACY - i);
		ck /= fract;

		sum += (ck / z);

		fract *= (-1.0 * i);
	}

	return (double)(sum * sc);
}

double beta(double a, double b) // funkcja beta
{
	return (gamma(a) * gamma(b)) / gamma(a + b);
}

double reg_beta(double x, double a, double b)  // regularized beta function
{
	if (x > (a + 1.0) / (a + b + 2.0)) 
	{
		return (1.0 - reg_beta(1.0 - x, b, a));
	}

	const double lim = 1.0e-8;
	const double min = 1.0e-30;
	const double bta = beta(a, b);
	const double f_part = (pow(x, a) * pow((1 - x), b)) / (a*beta(a, b));
	
	double fract = 1.0, c = 1.0, d = 0.0;

	int m;
	for (int i = 0; i <= 200; ++i) {
		m = i / 2;

		double numerator;
		if (i == 0) {
			numerator = 1.0;
		}
		else if (i % 2 == 0) {
			numerator = (m*(b - m)*x) / ((a + 2.0*m - 1.0)*(a + 2.0*m)); 
		}
		else {
			numerator = -((a + m)*(a + b + m)*x) / ((a + 2.0*m)*(a + 2.0*m + 1));		}

		d = 1.0 + numerator * d;
		if (fabs(d) < min) d = min;
		d = 1.0 / d;

		c = 1.0 + numerator / c;
		if (fabs(c) < min) c = min;

		const double cd = c * d;
		fract *= cd;

		if (fabs(1.0 - cd) < lim) {
		return f_part * (fract - 1.0);
		}
	}

	return 1.0; // jeœli tu dotrzemy to mamy b³¹d - narazie zwracam 1, bo taka wartoœc nie mo¿e byæ rezultatem tej funcji, ale mo¿n to jakoœ inaczej obs³u¿yæ
}