#include <stdio.h>
#include "../Calc/headers/quanc8.h"
#include "../Calc/headers/lagrange.h"
#include "../Calc/headers/spline.h"
#include <math.h>
#include <cmath>

double current_x;


double Function(double t) {
    return exp(current_x * t) * sin(t);
}

using namespace std;

int main() {
    const double min_x = 0, max_x = 2, h = 0.2;
    const int firstNodeDigit = 1 + (max_x - min_x) / h;
    const int secondNodeDigit = 10;

    const int down = 0, up = 1;
    const double abserr = 1.0e-12, relerr = 0;
    double errest, flag;
    int nofun[firstNodeDigit];
    double x[firstNodeDigit], y[firstNodeDigit];//???? ? ???????? ??????? ? ?????
    for (int i = 0; i < firstNodeDigit; ++i) {
        x[i] = h * i;
        current_x = x[i];
        quanc8(Function, down, up, abserr, relerr, &y[i], &errest, nofun, &flag);
    }

    
    double splineFunction[secondNodeDigit], b[secondNodeDigit], c[secondNodeDigit], d[secondNodeDigit];
    double lagrangeFunction[secondNodeDigit];
    double defaultFunction[secondNodeDigit];

    
    spline(secondNodeDigit, x, y, b, c, d);

    
    for (int k = 0; k < secondNodeDigit; ++k) {
        double xk = ((k + 1) - 0.5) * h;
        current_x = xk;
        quanc8(Function, down, up, abserr, relerr, &defaultFunction[k], &errest, nofun, &flag);
        lagrangeFunction[k] = lagrange(secondNodeDigit, x, y, xk);
        splineFunction[k] = seval(secondNodeDigit, &xk, x, y, b, c, d);
    }

    
    printf("------------------------------------------------------------------------\n");
    printf("k |  xk  |    Values          |   Lagrange         |    Spline\n");
    printf("------------------------------------------------------------------------\n");
    for (int k = 0; k < secondNodeDigit; ++k) {
        printf("%-*.d| ", 2, k + 1);
        printf("%.2f | ", ((k + 1) - 0.5) * h);
        printf("%.16f | ", defaultFunction[k]);
        printf("%.16f | ", lagrangeFunction[k]);
        printf("%.16f", splineFunction[k]);
        printf("\n");
    }

    
    printf("---------------------------------------------------\n");
    printf("k |  xk  |   Lagrange          |   Spline             \n");
    printf("---------------------------------------------------\n");
    for (int k = 0; k < secondNodeDigit; ++k) {
        printf("%-*.d| ", 2, k + 1);
        printf("%.2f | ", ((k + 1) - 0.5) * h);
        printf("%.16f | ", abs(defaultFunction[k] - lagrangeFunction[k]));
        printf("%.16f | ", abs(defaultFunction[k] - splineFunction[k]));
        printf("\n");
    }

    
    
    printf("flag\n");
    printf("%.16f\n", 0.232323232);
    printf("errest\n");
    printf("%.16f\n", errest);
    printf("nofun\n");
    for (int j = 0; j < firstNodeDigit; ++j) {
        printf("%.16d\n", nofun[j]);

    }


}



