#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

double exactSolution(double x) {
    return std::tan(x) - x;
}

double f(double x, double u) {
    return std::pow(u + x, 2);
}

double newtonMethod(double x, double u_prev, double h) {
    double u = u_prev; 
    const double tol = 1e-8; 
    const int max_iter = 100; 
    for (int iter = 0; iter < max_iter; ++iter) {
        double fu = u - u_prev - h / 2 * (f(x, u_prev) + f(x + h, u));
        double dfu = 1 - h / 2 * 2 * (u + x + h); 
        double delta = fu / dfu;
        u -= delta;
        if (std::abs(delta) < tol) break;
    }
    return u;
}

std::vector<double> implicitTrapezoidal(double a, double b, double h, double u0) {   // Неявный метод трапеций
    int n = static_cast<int>((b - a) / h);
    std::vector<double> u(n + 1);
    u[0] = u0;
    for (int i = 0; i < n; ++i) {
        double x = a + i * h;
        u[i + 1] = newtonMethod(x, u[i], h);
    }
    return u;
}

std::vector<double> rungeKutta2(double a, double b, double h, double u0) {   // Метод Рунге-Кутты второго порядка
    int n = static_cast<int>((b - a) / h);
    std::vector<double> u(n + 1);
    u[0] = u0;
    for (int i = 0; i < n; ++i) {
        double x = a + i * h;
        double k1 = f(x, u[i]);
        double k2 = f(x + h / 2, u[i] + h / 2 * k1);
        u[i + 1] = u[i] + h * k2;
    }
    return u;
}

double computeMaxError(const std::vector<double>& x_values, const std::vector<double>& y_values) {
    double max_error = 0.0;
    for (size_t i = 0; i < x_values.size(); ++i) {
        double error = std::abs(exactSolution(x_values[i]) - y_values[i]);
        if (error > max_error) {
            max_error = error;
        }
    }
    return max_error;
}

int main() {
    setlocale(LC_ALL, "ru");
    double a = 0.0, b = 1.0;
    double h1 = 0.1;         
    double h2 = h1 / 2;     
    double u0 = 0.0;       

    std::vector<double> u_trap_h1 = implicitTrapezoidal(a, b, h1, u0);
    std::vector<double> u_trap_h2 = implicitTrapezoidal(a, b, h2, u0);

    std::vector<double> u_rk_h1 = rungeKutta2(a, b, h1, u0);
    std::vector<double> u_rk_h2 = rungeKutta2(a, b, h2, u0);

    std::vector<double> x_values_h1(static_cast<size_t>((b - a) / h1) + 1);
    std::vector<double> x_values_h2(static_cast<size_t>((b - a) / h2) + 1);
    for (size_t i = 0; i < x_values_h1.size(); ++i) {
        x_values_h1[i] = a + i * h1;
    }
    for (size_t i = 0; i < x_values_h2.size(); ++i) {
        x_values_h2[i] = a + i * h2;
    }

    double error_trap_h1 = computeMaxError(x_values_h1, u_trap_h1);
    double error_trap_h2 = computeMaxError(x_values_h2, u_trap_h2);
    double error_rk_h1 = computeMaxError(x_values_h1, u_rk_h1);
    double error_rk_h2 = computeMaxError(x_values_h2, u_rk_h2);

    // Правило Рунге
    double p = 2; 
    double R_trap = std::abs(error_trap_h1 - error_trap_h2) / (std::pow(2, p) - 1);
    double R_rk = std::abs(error_rk_h1 - error_rk_h2) / (std::pow(2, p) - 1);

    std::cout << "Результаты для неявного метода трапеций:" << std::endl;
    std::cout << "Шаг h1 = " << h1 << ", Максимальная ошибка = " << error_trap_h1 << std::endl;
    std::cout << "Шаг h2 = " << h2 << ", Максимальная ошибка = " << error_trap_h2 << std::endl;
    std::cout << "Оценка по правилу Рунге: " << R_trap << std::endl;
    std::cout << std::endl;

    std::cout << "Результаты для метода Рунге-Кутты (порядок 2):" << std::endl;
    std::cout << "Шаг h1 = " << h1 << ", Максимальная ошибка = " << error_rk_h1 << std::endl;
    std::cout << "Шаг h2 = " << h2 << ", Максимальная ошибка = " << error_rk_h2 << std::endl;
    std::cout << "Оценка по правилу Рунге: " << R_rk << std::endl;

    return 0;
}
