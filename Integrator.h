#pragma once
#include <cmath>

// Класс Integrator, реализующий методы трапеций и Симпсона
class Integrator {
public:
    // Метод интегрирования по формуле трапеций
    double trapezoid(double (*f)(double), double a, double b, int n);

    // Метод интегрирования по формуле Симпсона
    double simpson(double (*f)(double), double a, double b, int n);

    // Порядок метода
    int order_trapezoid();
    int order_simpson();
};
