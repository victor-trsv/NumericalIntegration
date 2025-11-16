#include "Integrator.h"

double Integrator::trapezoid(double (*f)(double), double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.5 * (f(a) + f(b));
    for (int i = 1; i < n; ++i) sum += f(a + i * h);
    return sum * h;
}

double Integrator::simpson(double (*f)(double), double a, double b, int n) {
    if (n % 2 == 1) n++;
    double h = (b - a) / n;
    double sum = f(a) + f(b);
    for (int i = 1; i < n; ++i) {
        double x = a + i * h;
        sum += (i % 2 == 0 ? 2 : 4) * f(x);
    }
    return sum * h / 3.0;
}

int Integrator::order_trapezoid() { return 2; }
int Integrator::order_simpson() { return 4; }
