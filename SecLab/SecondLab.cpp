#include <iostream>
#include <stdio.h>
#include <cmath>
#include "../Calc/headers/matrix.h"


#define matrixSize 8

using namespace std;

int main() {
    MATRIX(matrix_A);
    VECTOR(vect, matrixSize);
    float cond[matrixSize];
    float work[matrixSize];

    int ipvt[matrixSize];
    MATRIX(matrix_A_copy);
    VECTOR(copyVect, matrixSize);

    for (double p = 1; p > 0.0000001; p *= 0.1) {


        matrix_A[0][0] = p + 6;
        matrix_A[0][1] = 2;
        matrix_A[0][2] = 6;
        matrix_A[0][3] = 8;
        matrix_A[0][4] = -2;
        matrix_A[0][5] = 1;
        matrix_A[0][6] = 8;
        matrix_A[0][7] = -5;


        matrix_A[1][0] = 6;
        matrix_A[1][1] = -22;
        matrix_A[1][2] = -2;
        matrix_A[1][3] = -1;
        matrix_A[1][4] = 0;
        matrix_A[1][5] = 5;
        matrix_A[1][6] = -6;
        matrix_A[1][7] = 4;

        matrix_A[2][0] = -2;
        matrix_A[2][1] = -3;
        matrix_A[2][2] = -16;
        matrix_A[2][3] = 0;
        matrix_A[2][4] = 0;
        matrix_A[2][5] = -4;
        matrix_A[2][6] = 2;
        matrix_A[2][7] = -5;

        matrix_A[3][0] = 1;
        matrix_A[3][1] = 1;
        matrix_A[3][2] = 4;
        matrix_A[3][3] = 9;
        matrix_A[3][4] = 1;
        matrix_A[3][5] = 0;
        matrix_A[3][6] = 0;
        matrix_A[3][7] = -6;

        matrix_A[4][0] = 0;
        matrix_A[4][1] = 2;
        matrix_A[4][2] = 0;
        matrix_A[4][3] = 2;
        matrix_A[4][4] = -3;
        matrix_A[4][5] = -5;
        matrix_A[4][6] = 7;
        matrix_A[4][7] = 5;

        matrix_A[5][0] = 6;
        matrix_A[5][1] = -2;
        matrix_A[5][2] = -4;
        matrix_A[5][3] = 2;
        matrix_A[5][4] = -8;
        matrix_A[5][5] = -12;
        matrix_A[5][6] = 3;
        matrix_A[5][7] = -3;

        matrix_A[6][0] = -6;
        matrix_A[6][1] = -6;
        matrix_A[6][2] = 0;
        matrix_A[6][3] = -8;
        matrix_A[6][4] = 0;
        matrix_A[6][5] = 5;
        matrix_A[6][6] = -15;
        matrix_A[6][7] = 0;

        matrix_A[7][0] = 0;
        matrix_A[7][1] = 7;
        matrix_A[7][2] = 6;
        matrix_A[7][3] = 0;
        matrix_A[7][4] = -5;
        matrix_A[7][5] = -8;
        matrix_A[7][6] = -5;
        matrix_A[7][7] = -3;



        for (int z = 0; z <= matrixSize; z++)
            for (int u = 0; u <= matrixSize; u++)
                matrix_A_copy[z][u] = matrix_A[z][u];


        printf("\n*************************************************************************\n");
        printf("p = ");
        printf("%7f",p);
        printf("\r\nDefault matrix:\r\n");
        for (int counti = 0; counti < matrixSize; counti++) {
            for (int countj = 0; countj < matrixSize; countj++)
                if (matrix_A[counti][countj] >= 0)
                    printf(" %.5f  ", matrix_A[counti][countj]);
                else printf("%.5f  ", matrix_A[counti][countj]);
            printf("\r\n");
        }
        for (int countVec = 0; countVec < matrixSize; countVec++)
            copyVect[countVec] = vect[countVec];



        MATRIX(tempA);
        MATRIX(invertible_A);
        VECTOR(tempX, matrixSize);
        tempX[0] = 1;
        tempX[1] = 0;
        tempX[2] = 0;
        tempX[3] = 0;
        tempX[4] = 0;
        tempX[5] = 0;
        tempX[6] = 0;
        tempX[7] = 0;


        for (int counti = 0; counti < matrixSize; counti++)
            for (int countj = 0; countj < matrixSize; countj++)
                tempA[counti][countj] = matrix_A_copy[counti][countj];

        decomp(matrixSize, tempA, cond, ipvt, work);
        solve(matrixSize, tempA, tempX, ipvt);

        invertible_A[0][0] = tempX[0];
        invertible_A[0][1] = tempX[1];
        invertible_A[0][2] = tempX[2];
        invertible_A[0][3] = tempX[3];
        invertible_A[0][4] = tempX[4];
        invertible_A[0][5] = tempX[5];
        invertible_A[0][6] = tempX[6];
        invertible_A[0][7] = tempX[7];

        //*******
        tempX[0] = 0;
        tempX[1] = 1;
        tempX[2] = 0;
        tempX[3] = 0;
        tempX[4] = 0;
        tempX[5] = 0;
        tempX[6] = 0;
        tempX[7] = 0;

        for (int counti = 0; counti < matrixSize; counti++)
            for (int countj = 0; countj < matrixSize; countj++)
                tempA[counti][countj] = matrix_A_copy[counti][countj];

        decomp(matrixSize, tempA, cond, ipvt, work);
        solve(matrixSize, tempA, tempX, ipvt);


        invertible_A[1][0] = tempX[0];
        invertible_A[1][1] = tempX[1];
        invertible_A[1][2] = tempX[2];
        invertible_A[1][3] = tempX[3];
        invertible_A[1][4] = tempX[4];
        invertible_A[1][5] = tempX[5];
        invertible_A[1][6] = tempX[6];
        invertible_A[1][7] = tempX[7];


        for (int counti = 0; counti < matrixSize; counti++)
            for (int countj = 0; countj < matrixSize; countj++)
                tempA[counti][countj] = matrix_A_copy[counti][countj];

        //*******
        tempX[0] = 0;
        tempX[1] = 0;
        tempX[2] = 1;
        tempX[3] = 0;
        tempX[4] = 0;
        tempX[5] = 0;
        tempX[6] = 0;
        tempX[7] = 0;
        decomp(matrixSize, tempA, cond, ipvt, work);
        solve(matrixSize, tempA, tempX, ipvt);

        invertible_A[2][0] = tempX[0];
        invertible_A[2][1] = tempX[1];
        invertible_A[2][2] = tempX[2];
        invertible_A[2][3] = tempX[3];
        invertible_A[2][4] = tempX[4];
        invertible_A[2][5] = tempX[5];
        invertible_A[2][6] = tempX[6];
        invertible_A[2][7] = tempX[7];


        for (int counti = 0; counti < matrixSize; counti++)
            for (int countj = 0; countj < matrixSize; countj++)
                tempA[counti][countj] = matrix_A_copy[counti][countj];

        //*******
        tempX[0] = 0;
        tempX[1] = 0;
        tempX[2] = 0;
        tempX[3] = 1;
        tempX[4] = 0;
        tempX[5] = 0;
        tempX[6] = 0;
        tempX[7] = 0;

        decomp(matrixSize, tempA, cond, ipvt, work);
        solve(matrixSize, tempA, tempX, ipvt);

        invertible_A[3][0] = tempX[0];
        invertible_A[3][1] = tempX[1];
        invertible_A[3][2] = tempX[2];
        invertible_A[3][3] = tempX[3];
        invertible_A[3][4] = tempX[4];
        invertible_A[3][5] = tempX[5];
        invertible_A[3][6] = tempX[6];
        invertible_A[3][7] = tempX[7];



        for (int counti = 0; counti < matrixSize; counti++)
            for (int countj = 0; countj < matrixSize; countj++)
                tempA[counti][countj] = matrix_A_copy[counti][countj];

        //*******
        tempX[0] = 0;
        tempX[1] = 0;
        tempX[2] = 0;
        tempX[3] = 0;
        tempX[4] = 1;
        tempX[5] = 0;
        tempX[6] = 0;
        tempX[7] = 0;

        decomp(matrixSize, tempA, cond, ipvt, work);
        solve(matrixSize, tempA, tempX, ipvt);

        invertible_A[4][0] = tempX[0];
        invertible_A[4][1] = tempX[1];
        invertible_A[4][2] = tempX[2];
        invertible_A[4][3] = tempX[3];
        invertible_A[4][4] = tempX[4];
        invertible_A[4][5] = tempX[5];
        invertible_A[4][6] = tempX[6];
        invertible_A[4][7] = tempX[7];



        for (int counti = 0; counti < matrixSize; counti++)
            for (int countj = 0; countj < matrixSize; countj++)
                tempA[counti][countj] = matrix_A_copy[counti][countj];

        //*******
        tempX[0] = 0;
        tempX[1] = 0;
        tempX[2] = 0;
        tempX[3] = 0;
        tempX[4] = 0;
        tempX[5] = 1;
        tempX[6] = 0;
        tempX[7] = 0;

        decomp(matrixSize, tempA, cond, ipvt, work);
        solve(matrixSize, tempA, tempX, ipvt);

        invertible_A[5][0] = tempX[0];
        invertible_A[5][1] = tempX[1];
        invertible_A[5][2] = tempX[2];
        invertible_A[5][3] = tempX[3];
        invertible_A[5][4] = tempX[4];
        invertible_A[5][5] = tempX[5];
        invertible_A[5][6] = tempX[6];
        invertible_A[5][7] = tempX[7];



        for (int counti = 0; counti < matrixSize; counti++)
            for (int countj = 0; countj < matrixSize; countj++)
                tempA[counti][countj] = matrix_A_copy[counti][countj];

        //*******
        tempX[0] = 0;
        tempX[1] = 0;
        tempX[2] = 0;
        tempX[3] = 0;
        tempX[4] = 0;
        tempX[5] = 0;
        tempX[6] = 1;
        tempX[7] = 0;

        decomp(matrixSize, tempA, cond, ipvt, work);
        solve(matrixSize, tempA, tempX, ipvt);

        invertible_A[6][0] = tempX[0];
        invertible_A[6][1] = tempX[1];
        invertible_A[6][2] = tempX[2];
        invertible_A[6][3] = tempX[3];
        invertible_A[6][4] = tempX[4];
        invertible_A[6][5] = tempX[5];
        invertible_A[6][6] = tempX[6];
        invertible_A[6][7] = tempX[7];




        for (int counti = 0; counti < matrixSize; counti++)
            for (int countj = 0; countj < matrixSize; countj++)
                tempA[counti][countj] = matrix_A_copy[counti][countj];

        //*******
        tempX[0] = 0;
        tempX[1] = 0;
        tempX[2] = 0;
        tempX[3] = 0;
        tempX[4] = 0;
        tempX[5] = 0;
        tempX[6] = 0;
        tempX[7] = 1;

        decomp(matrixSize, tempA, cond, ipvt, work);
        solve(matrixSize, tempA, tempX, ipvt);

        invertible_A[7][0] = tempX[0];
        invertible_A[7][1] = tempX[1];
        invertible_A[7][2] = tempX[2];
        invertible_A[7][3] = tempX[3];
        invertible_A[7][4] = tempX[4];
        invertible_A[7][5] = tempX[5];
        invertible_A[7][6] = tempX[6];
        invertible_A[7][7] = tempX[7];

        //������������� ���������� �������
        double costul[matrixSize][matrixSize];
        for (int counti = 0; counti < matrixSize; counti++) {
            for (int countj = 0; countj < matrixSize; countj++) {
                costul[counti][countj] = invertible_A[countj][counti];
            }
        }

        for (int counti = 0; counti < matrixSize; counti++) {
            for (int countj = 0; countj < matrixSize; countj++) {
                invertible_A[counti][countj] = costul[counti][countj];
            }
        }



        printf("\r\nInvertible matrix:\r\n");
        for (int counti = 0; counti < matrixSize; counti++) {
            for (int countj = 0; countj < matrixSize; countj++)
                if (invertible_A[counti][countj] >= 0)
                    printf(" %.5f  ", invertible_A[counti][countj]);
                else printf("%.5f  ", invertible_A[counti][countj]);
            printf("\r\n");
        }




        printf("\r\n");
        printf("Residual matrix:\n");
        double result_of_multiplying_Matrix[matrixSize][matrixSize];
        for (int i = 0; i < matrixSize; ++i) {
            for (int j = 0; j < matrixSize; ++j) {
                if (i == j)
                    result_of_multiplying_Matrix[i][j] = -1;
                else
                    result_of_multiplying_Matrix[i][j] = 0;
                for (int k = 0; k < matrixSize; ++k)
                    result_of_multiplying_Matrix[i][j] += matrix_A_copy[i][k] * invertible_A[k][j];
                printf("%.5f  ", result_of_multiplying_Matrix[i][j]);

            }
            printf("\r\n");
        }



        double norma = 0;
        for (int i = 0; i < matrixSize; ++i)
            for (int j = 0; j < matrixSize; ++j)
                norma += result_of_multiplying_Matrix[i][j] * result_of_multiplying_Matrix[i][j];

        norma = sqrt(norma);

        printf("\nEvklids norm\r\n");
        printf("%.5f\n", norma);

        printf("Cond\r\n");
        printf("%.5f\n", cond[0]);


    }
    return 1;
}
