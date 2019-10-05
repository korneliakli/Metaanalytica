#include <math.h>
#include "betagamma.h"

/*
W t-te�cie chodzi o por�wnanie dw�ch �rednich.
Zaimplementowa�am dwa rodzaje - dla jednej pr�by i dla pr�b niezale�nych, jest jeszcze dla pr�b zale�nych, zrobi� j� jak zaakceptujesz te.
Skoro masz ju� funkcje licz�ce �redni� i sd to wrzuci�am to od razu do parametr�w.
Wyniki pezentuje si�:  �rednia1 +/- sd1 | �rednia2 +/- sd2 (tudzie� dla populacji w pojedynczym) | p-value
*/

//t-test dla pojedynczej pr�by
double single_ttest(double mean, double sd, double mean_pop) //mean - �rednia z analizowanych danych, sd - odchylenie standardowe, mean_pop - �rednia z por�wnywanej populacji (podaje u�ytkownik)
{
	return (mean - mean_pop) / sd;
}

//t-test dla pr�b niezale�nych
double indep_ttest(double mean1, double mean2, double sd1, double sd2, int n1, int n2) //mean1, mean2 - �rednie z por�wnywanych grup, sd1, sd2 - odchylenia standardowe, n1, n2 - liczebno�ci grup
{
	double sd_temp1 = (n1 - 1) * pow(sd1, 2);
	double sd_temp2 = (n2 - 1) * pow(sd2, 2);
	double sd_temp3 = sd_temp1 + sd_temp2;
	double sd_temp4 = sd_temp3 / (n1 + n2 - 2);
	double sd_sqr = sqrt(sd_temp4 * ((1 / n1) + (1 / n2)));

	return (mean1 - mean2) / sd_sqr;
}

//obliczanie p-value
double calc_pval(double t, double df) // t - warto�� obliczona z jednej z poprzednich funkcji, df - stopnie swobody: dla single_ttest = n - 1, dla indep_ttest = n1 + n2 - 2
{
	x = df / (pow(t, 2) + df);

	return reg_beta(x, df / 2, 0.5);
}