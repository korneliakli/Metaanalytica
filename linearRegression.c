#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "matrix.h"

double sse(int n, double y[n][1], double y_dash[n][1])
{
    double res = 0;
    
    for(int i = 0; i < n; i++)
    {
        res += pow((y[i][1] - y_dash[i][1]), 2);
    }
    
    return res;
}

void standardError(int n, int m, double y[n][1], double x[n][m], double x_inv[m][m], double betas[m][1], double stdErrors[m][1])
{
    double y_dash[n][1];
    
    multiplyMatrices(n, m, m, 1, x, betas, y_dash);
    double sumSqErr = sse(n, y, y_dash);
    
    double mse = sumSqErr/(n - m - 1);
    
    for(int i = 0; i < m; i++)
    {
        for(int j = 0; j < m; j++)
        {
            x_inv[i][j] = x_inv[i][j]*mse;
        }
    }
    
    
    diagonal(m, x_inv, stdErrors);
    
    for(int i = 0; i < m; i++)
    {
        stdErrors[i][1] = sqrt(stdErrors[i][1]);
    }
}

void linearRegression(int n, int m, double x[n][m], double y[n][1], double betas[m][1], double stdErrors[m][1])
{
    double transp_x[m][n];
    
    transpose(n, m, x, transp_x);
    
    double x_mult[m][m];
    
    multiplyMatrices(m, n, n, m, transp_x, x, x_mult);
    
    
    double x_inv[m][m];
    int P[m];
    
    LUPDecompose(m, x_mult, P);
    LUPInvert(m, x_mult, P, x_inv);
    
    
    double x_inv_transp[m][n];
    
    multiplyMatrices(m, m, m, n, x_inv, transp_x, x_inv_transp);
    multiplyMatrices(m, n, n, 1, x_inv_transp, y, betas);
    
    standardError(n, m, y, x, x_inv, betas, stdErrors);
}

double calc_pval(int m, double stdErrors[m][1], double betas[m][1], double df, double pvals[m][1]) // t - wartość obliczona z jednej z poprzednich funkcji, df - stopnie swobody: dla single_ttest = n - 1, dla indep_ttest = n1 + n2 - 2
{
	double x[m][1];
    double t[m][1];
    
    for(int i = 0; i < m; i++)
    {
        t[i][0] = betas[i][0] / stdErrors[i][0];
    }
    for(int i = 0; i < m; i++)
    {
         x[i][0] = df / (pow(t[i][0], 2) + df);
    }
    
    for(int i = 0; i < m; i++)
    {
	    pvals[i][0] = reg_beta(x[i][0], df / 2, 0.5);
    }
}
