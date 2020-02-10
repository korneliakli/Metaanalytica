#pragma once

extern double single_ttest(double mean, double sd, double mean_pop);
extern double indep_ttest(double mean1, double mean2, double sd1, double sd2, int n1, int n2);
extern double calc_pval(double t, double df);