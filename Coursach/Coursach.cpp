#include <iostream>
#include <math.h>
#include "../Calc/headers/rkf45.h"
#include "../Calc/headers/quanc8.h"
#include "../Calc/headers/fmin.h"

using namespace ::std;

//�������� �������
int start = 36, xout = 46;
double x[] = {0.0, 0.303, -0.465, 0.592, -0.409, 0.164, 0.180};
double t[] = {0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
const int nodeDigit = 7;


//��������� �������
const double M = 1;//����� ��������
const double g = 9.81;
double L;//��������� ����� �������
const double divisorOfL = 0.90452424;
double K;//��������� �������
double y[4];// y[0] = y, y[1] = y', y[2] = z, y[3] = z'
double dy[4];//�������������


//��������� QUANC8
const int down = 0, up = 1;
const double abserr = 1.0e-12, relerr = 0;
double errest, flag;
int nofun;

//��������� RKF45
double h = 0.4, rBottom = 0, rUp = 2.4, tOut = 0, re = 1e-8, ae = 1e-8;
int iflag = 1, n = 4;
double work[27];
int iwork[30];

//��������� FMIN
double ERROR_FMIN = 0.1;
int flag_fmin = 1;
double K1;


char line[] = "--------------------------------------------------\n";


//��������������� ������� L
double funL(double x) {
    return cos(x * x);
}


//����������� ������� ��
void fun(double t, double *y, double *dy) {
    dy[0] = y[1];
    dy[1] = -K / M * y[0] - g * (1 - cos(y[2])) + (L + y[0]) * (y[3] * y[3]);
    dy[2] = y[3];
    dy[3] = -g / (L + y[0]) * sin(y[2]) - 2 / (L + y[0]) * y[1] * y[3];
}


//��������� ���������� RKF45
void setUpRKF45param() {
    iflag = 1;
    rBottom = 0;
    y[0] = 0;
    y[1] = 0;
    y[2] = 0;
    y[3] = 4;
}


//����� K ��� ������ fmin, ������ � ����� �����������
double calcFun(double k) {
    K = k;
    double sum = 0;
    printf("K = %.2f\n", K);
    for (int i = 0; i < nodeDigit; ++i) {
        setUpRKF45param();
        tOut = t[i];
        RKF45(fun, n, y, &rBottom, &tOut, &re, &ae, &iflag, work, iwork);
        sum += pow(y[0] - x[i], 2);
        printf("tout = %.1f \t ", t[i]);
        printf("x*  = %.6f \t", x[i]);
        printf("x[] = %.6f\n", y[0]);
    }
    printf("errors = %.16f\n", sum);
    printf(line);

    return sum;

}


int main() {
    //����������� � ������� �������� L
    quanc8(funL, down, up, abserr, relerr, &L, &errest, &nofun, &flag);
    L = L / divisorOfL;
    printf("L = ");
    printf("%.10f\n", L);
    printf(line);

    //������� K ��� ������ fmin
    K = fmin(start, xout, calcFun, ERROR_FMIN, K1, flag_fmin);
    printf("Answer:\nK = %.2f\n", K);


    return 0;

}
