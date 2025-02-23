#include <iostream>
#include <cmath>
#include <iomanip>
#include <clocale>

double f(double x) {
	return pow(4, x) - 5 * x - 2;
}

double f_(double x) {
	return pow(4, x) * log(4) - 5;
}

double dih(double a, double b, double eps) {
	double x;
	int k = 0;
	while (fabs(b - a) > 2 * eps) {
		x = (a + b) / 2;
		std::cout << "k=" << k << ", a=" << a << ", b=" << b << ", f(a)=" << f(a) << ", f(b)=" << f(b) << ", b-a=" << b - a << '\n';
		if (f(x) == 0)
			break;
		else if (f(x) * f(a) < 0)
			b = x;
		else
			a = x;
		k++;
	}
	return x;
}

double Newton(double x_k, double eps) {
	double x_k1;
	int k = 0;
	for (;;) {
		x_k1 = x_k - f(x_k) / f_(x_k);
		std::cout << "k=" << k << " , x_k=" << x_k << ", |x_k - x_k-1|=" << fabs(x_k - x_k1) << std::endl;
		if (fabs(x_k1 - x_k) < eps)
			break;
		x_k = x_k1;
		k++;
	}
	return x_k1;
}

double Newton2(double x_k, double eps) {
	double x_k1;
	double der = f_(x_k);
	int k = 0;
	for (;;) {
		x_k1 = x_k - f(x_k) / der;
		std::cout << "k=" << k << " , x_k=" << x_k << ", |x_k - x_k-1|=" << fabs(x_k - x_k1) << std::endl;
		if (fabs(x_k1 - x_k) < eps)
			break;
		x_k = x_k1;
		k++;
	}
	return x_k1;
}

double Sek(double x_k0, double x_k, double eps) {
	double x_k1;
	int k = 0;
	for (;;) {
		x_k1 = x_k - f(x_k) * (x_k - x_k0) / (f(x_k) - f(x_k0));
		std::cout << "k=" << k << ", x_k=" << x_k << ", |x_k - x_k-1|=" << fabs(x_k - x_k0) << std::endl;
		if (fabs(x_k1 - x_k) < eps) 
			break;
		x_k0 = x_k;
		x_k = x_k1;
		k++;
	}
	return x_k;
}

int main() {
	setlocale(LC_ALL, "ru");

	double a = 0.0, b = 2.0;
	double eps_dih = 1e-2;
	double eps = 1e-7;

	std::cout << "Метод Дихотомии: " << '\n';
	double dikhotomiya_root = dih(a, b, eps_dih);
	std::cout << "Корень, найденный методом дихотомии: " << dikhotomiya_root << '\n' << '\n';

	std::cout << "Метод Ньютона: " << '\n';
	double newton_root = Newton(dikhotomiya_root, eps);
	std::cout << "Корень, найденный методом Ньютона: " << newton_root << '\n' << '\n';

	
	std::cout << "Метод Ньютона с постоянной производной: " << '\n'; 
	double newton2_root = Newton2(dikhotomiya_root, eps);
	std::cout << "Корень, найденный методом Ньютона с постоянной поизводной: " << newton2_root << '\n' << '\n';

	std::cout << "Метод Секущих: " << '\n';
	double sec_root = Sek(dikhotomiya_root, dikhotomiya_root + 0.1, eps);
	std::cout << "Корень, найденный методом секущих: " << sec_root << '\n' << '\n';
}
