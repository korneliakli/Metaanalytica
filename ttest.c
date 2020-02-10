#include <math.h>
#include "betagamma.h"
#include "ttest.h"

/*
W t-teście chodzi o porównanie dwóch średnich.
Zaimplementowałam dwa rodzaje - dla jednej próby i dla prób niezależnych, jest jeszcze dla prób zależnych, zrobię ją jak zaakceptujesz te.
Skoro masz już funkcje liczące średnią i sd to wrzuciłam to od razu do parametrów.
Wyniki pezentuje się:  średnia1 +/- sd1 | średnia2 +/- sd2 (tudzież dla populacji w pojedynczym) | p-value
*/

//t-test dla pojedynczej próby
double single_ttest(double mean, double sd, double mean_pop) //mean - średnia z analizowanych danych, sd - odchylenie standardowe, mean_pop - średnia z porównywanej populacji (podaje użytkownik)
{
	return (mean - mean_pop) / sd;
}

//t-test dla prób niezależnych
double indep_ttest(double mean1, double mean2, double sd1, double sd2, int n1, int n2) //mean1, mean2 - średnie z porównywanych grup, sd1, sd2 - odchylenia standardowe, n1, n2 - liczebności grup
{
	double sd_temp1 = (n1 - 1) * pow(sd1, 2);
	double sd_temp2 = (n2 - 1) * pow(sd2, 2);
	double sd_temp3 = sd_temp1 + sd_temp2;
	double sd_temp4 = sd_temp3 / (n1 + n2 - 2);
	double sd_sqr = sqrt(sd_temp4 * ((1 / n1) + (1 / n2)));

	return (mean1 - mean2) / sd_sqr;
}

//obliczanie p-value
double calc_pval(double t, double df) // t - wartość obliczona z jednej z poprzednich funkcji, df - stopnie swobody: dla single_ttest = n - 1, dla indep_ttest = n1 + n2 - 2
{
	x = df / (pow(t, 2) + df);

	return reg_beta(x, df / 2, 0.5);
}
