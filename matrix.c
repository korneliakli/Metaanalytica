#include<stdio.h>
#include<math.h>
#include<stdlib.h>

void display(int row, int col, double mat[row][col]) 
{ 
    for (int i = 0; i < row; i++) 
    { 
        for (int j = 0; j < col; j++) 
            printf("  %f", mat[i][j]); 
        printf("\n"); 
    } 
} 

void transpose(int n, int m, double matrix[n][m], double result[m][n])  
{
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < m; j++)
        {
            result[j][i] = matrix[i][j];
        }
    }
}

void generateMat(int n, int m, double mat[n][m])
{
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < m; j++)
        {
            mat[i][j] = (rand() % 100) + 1;
        }
    }
}

void switchRows(int n, double mat[n][n], int r1, int r2)
{
    double temp[n];
    
    for(int i = 0; i < n; i++)
    {
        temp[i] = mat[r1][i];
        mat[r1][i] = mat[r2][i];
        mat[r2][i] = temp[i];
    }
}

int LUPDecompose(int n, double mat[n][n], int P[n]) 
{
    int i, j, k, imax; 
    double maxA, absA;
  

    for (i = 0; i <= n; i++)
        P[i] = i; 

    for (i = 0; i < n; i++) 
    {
        maxA = 0.0;
        imax = i;

        for (k = i; k < n; k++) {
            absA = fabs(mat[k][i]);
            if (absA > maxA) { 
                maxA = absA;
                imax = k;
            }
        }

        if (maxA < 0) return 0; 

        if (imax != i) {
            
            j = P[i];
            P[i] = P[imax];
            P[imax] = j;
            switchRows(n, mat, i, imax);
            P[n]++;
        }

        for (j = i + 1; j < n; j++) {
            mat[j][i] /= mat[i][i];

            for (k = i + 1; k < n; k++)
                mat[j][k] -= mat[j][i] * mat[i][k];
        }
    }

    return 1;  
}



void LUPInvert(int n, double mat[n][n], int P[n], double res[n][n]) 
{
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            if (P[i] == j) 
                res[i][j] = 1.0;
            else
                res[i][j] = 0.0;

            for (int k = 0; k < i; k++)
                res[i][j] -= mat[i][k] * res[k][j];
        }

        for (int i = n - 1; i >= 0; i--) {
            for (int k = i + 1; k < n; k++)
                res[i][j] -= mat[i][k] * res[k][j];

            res[i][j] = res[i][j] / mat[i][i];
        }
    }
}


double LUPDeterminant(double **A, int *P, int N) {

    double det = A[0][0];

    for (int i = 1; i < N; i++)
        det *= A[i][i];

    if ((P[N] - N) % 2 == 0)
        return det; 
    else
        return -det;
}

void multiplyMatrices(int n1, int m1, int n2, int m2, double mat1[n1][m1], double mat2[n2][m2], double result[n1][m2])
{
    for(int i = 0; i < n1; i++)
	{
		for(int j = 0; j < m2; j++)
		{
		    result[i][j] = 0;
		}
	}
	
	for(int i = 0; i < n1; i++)
	{
		for(int j = 0; j < m1; j++)
		{
			for(int k = 0; k < m2; k++)
			{
				result[i][k] += mat1[i][j] * mat2[j][k];
			}
		}
	}

}

double diagonal(int n, double matrix[n][n], double diag[n][1])
{
    for(int i = 0; i < n; i++)
    {
        diag[i][0] = matrix[i][i];
    }
}
