#pragma once

extern void linearRegression(int n, int m, double x[n][m], double y[n][1], double betas[m][1], double stdErrors[m][1]);
extern void calc_pval_linreg(int m, double stdErrors[m][1], double betas[m][1], double df, double pvals[m][1]);
