#pragma once

void display(int row, int col, double mat[row][col]);
void transpose(int n, int m, double matrix[n][m], double result[m][n]);
void generateMat(int n, int m, double mat[n][m]);
void switchRows(int n, double mat[n][n], int r1, int r2);
int LUPDecompose(int n, double mat[n][n], int P[n]);
void LUPInvert(int n, double mat[n][n], int P[n], double res[n][n]);
void multiplyMatrices(int n1, int m1, int n2, int m2, double mat1[n1][m1], double mat2[n2][m2], double result[n1][m2]);
double diagonal(int n, double matrix[n][n], double diag[n][1]);
