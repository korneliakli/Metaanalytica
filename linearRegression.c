#include "matrix.h"


void linearRegression(int n, int m, double x[n][m], double y[n][1], double betas[m][1], double stdErrors[m][1], double p_vals[m][1])
{
    double transp_x[m][n];
    
    transpose(n, m, x, transp_x);
    
    double x_mult[m][m];
    
    multiplyMatrices(m, n, n, m, transp_x, x, x_mult);
    
    
    double x_inv[m][m];
    double P[m];
    
    LUPDecompose(m, x_mult, P);
    LUPInvert(m, x_mult, P, x_inv);
    
    
    double x_inv_transp[m][n];
    
    multiplyMatrices(m, m, m, n, x_inv, transp_x, x_inv_transp);
    multiplyMatrices(m, n, n, 1, x_inv_transp, y, betas);
}
