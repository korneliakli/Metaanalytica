#include "betagamma.h"
#include <math.h>

/*
Test Pearsona sprawdza korelację pomiędzy dwiema zmiennymi.
Można robić zbiorcze korelacje pomiędzy większą ilością zmiennych, ale w dalszym ciągu jest to porównanie każdy z każdym.
Do obliczenia p-value korzysta z tej zamej funckji co t-test.
*/


double pearson(double x[], double y[], int size, double mean_x, double mean_y) //x[], y[] - tablica z wszystkimi wartościami zmiennych, size - rozmiar obu tablic (musi być taki sam, myślę, ze można to sprawdzać jeszcze przed wywołaniem funkcji), mean_x, mean_y - średnie wartości obu zmiennych
{
    double tempsum = 0;
    double sum_ysq = 0;
    double sum_xsq = 0;
    
    for(int i = 0; i < size; i++)
    {
        tempsum += (x[i] - mean_x) * (y[i] - mean_y);
        sum_xsq += pow((x[i] - mean_x), 2);
        sum_ysq += pow((y[i] - mean_y), 2);
    }
    
    return tempsum / (sqrt(sum_xsq) * sqrt(sum_ysq));
}


double calc_pval(double r, double df) // r - wartość obliczona z funkcji pearson(), df - stopnie swobody: parametr size z funkcji pearson - 2
{
	double t = r * sqrt(6 / (1 - pow(r, 2)));
	double x = df / (pow(t, 2) + df);

	return reg_beta(x, df / 2, 0.5);
}

//wynik prezentuje się jako wartość r zwrócona z funkcji pearson() (zawsze < 1) i p-value
