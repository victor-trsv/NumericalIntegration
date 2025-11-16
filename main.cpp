#include <iostream>
#include <vector>
#include <iomanip>
#include <windows.h>
#include "Integrator.h"
#include "richardson.h"

int main() {
    SetConsoleOutputCP(1251);

    auto f = [](double x){ return std::exp(x); };
    double a = -0.07, b = 0.73;
    double I_star = std::exp(b) - std::exp(a);

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "I* (analytical) = " << I_star << "\n\n";

    std::vector<int> N = {4,8,16,32};
    std::vector<double> h_values;
    for (int n : N) h_values.push_back((b-a)/(double)n);

    Integrator integrator;

    std::vector<double> I_trap, I_simp;
    for (int n : N) {
        I_trap.push_back(integrator.trapezoid(f,a,b,n));
        I_simp.push_back(integrator.simpson(f,a,b,n));
    }

    auto printTable = [&](const std::vector<double>& I, int k, const std::string& name){
        std::cout << name << " (k=" << k << ")\n";
        std::cout << "| h       | I*-I^h      | (I*-I^h)/(I*-I^{h/2}) | (I^{h/2}-I^h)/(2^k-1) | I^R        | I*-I^R      |\n";
        std::cout << "|---------|-------------|------------------------|----------------------|------------|-------------|\n";
        for(size_t i=0;i<I.size()-1;i++){
            double h = h_values[i], Ih=I[i], Ih2=I[i+1];
            double error_h = I_star - Ih, error_h2 = I_star - Ih2;
            double ratio_errors = error_h / error_h2;
            double runge = (Ih2 - Ih)/(std::pow(2,k)-1);
            double IR = richardson(Ih,Ih2,k);
            double error_R = I_star - IR;
            std::cout << "| " << h << " | "
                      << error_h << " | "
                      << std::setw(22) << ratio_errors << " | "
                      << std::setw(20) << runge << " | "
                      << IR << " | "
                      << error_R << " |\n";
        }
        double p = std::log2(std::abs((I_star-I[1])/(I_star-I[2])));
        std::cout << "\nПриближённый порядок аппроксимации = " << p << "\n\n";
    };

    printTable(I_trap, integrator.order_trapezoid(), "TRAPEZOIDAL RULE");
    printTable(I_simp, integrator.order_simpson(), "SIMPSON RULE");

    std::cout << "Значения интеграла для сравнения:\n";
    std::cout << "Trapezoidal: ";
    for(double val : I_trap) std::cout << val << " ";
    std::cout << "\nSimpson:     ";
    for(double val : I_simp) std::cout << val << " ";
    std::cout << "\n";

    return 0;
}
