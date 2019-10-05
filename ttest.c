#include <math.h>
#include "betagamma.h"

/*
W t-teœcie chodzi o porównanie dwóch œrednich.
Zaimplementowa³am dwa rodzaje - dla jednej próby i dla prób niezale¿nych, jest jeszcze dla prób zale¿nych, zrobiê j¹ jak zaakceptujesz te.
Skoro masz ju¿ funkcje licz¹ce œredni¹ i sd to wrzuci³am to od razu do parametrów.
Wyniki pezentuje siê:  œrednia1 +/- sd1 | œrednia2 +/- sd2 (tudzie¿ dla populacji w pojedynczym) | p-value
*/

//t-test dla pojedynczej próby
double single_ttest(double mean, double sd, double mean_pop) //mean - œrednia z analizowanych danych, sd - odchylenie standardowe, mean_pop - œrednia z porównywanej populacji (podaje u¿ytkownik)
{
	return (mean - mean_pop) / sd;
}

//t-test dla prób niezale¿nych
double indep_ttest(double mean1, double mean2, double sd1, double sd2, int n1, int n2) //mean1, mean2 - œrednie z porównywanych grup, sd1, sd2 - odchylenia standardowe, n1, n2 - liczebnoœci grup
{
	double sd_temp1 = (n1 - 1) * pow(sd1, 2);
	double sd_temp2 = (n2 - 1) * pow(sd2, 2);
	double sd_temp3 = sd_temp1 + sd_temp2;
	double sd_temp4 = sd_temp3 / (n1 + n2 - 2);
	double sd_sqr = sqrt(sd_temp4 * ((1 / n1) + (1 / n2)));

	return (mean1 - mean2) / sd_sqr;
}

//obliczanie p-value
double calc_pval(double t, double df) // t - wartoœæ obliczona z jednej z poprzednich funkcji, df - stopnie swobody: dla single_ttest = n - 1, dla indep_ttest = n1 + n2 - 2
{
	x = df / (pow(t, 2) + df);

	return reg_beta(x, df / 2, 0.5);
}