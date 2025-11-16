#pragma once

// Экстраполяция Ричардсона
inline double richardson(double Ih, double Ih2, int k) {
    return Ih2 + (Ih2 - Ih) / (std::pow(2, k) - 1);
}
