#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

double f(double x) {
    double numerator = pow(x - 1, 2);
    double denominator = pow(x, 2) + exp(x) + 1;
    return numerator / denominator;
}

double exact_integral() {
    return 2 + log(2.0 / (5 + exp(2)));
}

double trapezoidal(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.5 * (f(a) + f(b));

    for (int i = 1; i < n; i++) {
        sum += f(a + i * h);
    }

    return h * sum;
}

double simpson_modified(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.0;

    for (int i = 0; i < n; i++) {
        double x_i = a + i * h;
        double x_i_half = a + (i + 0.5) * h;
        double x_i_plus_1 = a + (i + 1) * h;

        sum += f(x_i) + 4 * f(x_i_half) + f(x_i_plus_1);
    }

    return h / 6 * sum;
}

void compute_integral(const string& method_name,
    double (*method)(double, double, int),
    double a, double b, double eps, int order) {
    std::cout << "\n" << method_name << "\n";
    std::cout << setw(15) << "N" << setw(15) << "h"
        << setw(25) << "Приближенное значение интеграла" << setw(25) << "Ошибка погрешности\n";

    int n = 2;
    double Q_prev = method(a, b, n);
    std::cout << setw(15) << n << setw(15) << (b - a) / n
        << setw(25) << Q_prev << setw(25) << "N/A\n";

    while (true) {
        n *= 2;
        double Q = method(a, b, n);
        double error = abs(Q_prev - Q) / (pow(2, order) - 1);

        std::cout << setw(15) << n << setw(15) << (b - a) / n
            << setw(25) << Q << setw(25) << error << std::endl;

        if (error < eps || abs(Q - Q_prev) < 1e-15) 
            break;
        Q_prev = Q;
    }
}

int main() {
    setlocale(LC_ALL, "ru");
    const double a = 0, b = 2;
    const double eps = 1e-7;

    // Метод трапеций (порядок точности 2)
    compute_integral("Метод трапеций", trapezoidal, a, b, eps, 2);

    // Модифицированный метод Симпсона (порядок точности 4)
    compute_integral("Метод Симпсона", simpson_modified, a, b, eps, 4);

    std::cout << "\nТочное значение: " << exact_integral() << std::endl;

    return 0;
}

// Задание 1