#include <iostream>
#include <cmath> 
#include <vector>
#include <iomanip>

const double PI = 3.14159265358979323846;

double f1(double x) {
	return sin(2 * x) * log(x + 5);
}

double f2(double x) {
	return sqrt(2 * fabs(x) + x * x);
}

std::vector<std::vector<double>> dividedDifferences(const std::vector<double>& x, const std::vector<double>& y) {
	int n = x.size();
	std::vector<std::vector<double>> table(n, std::vector<double>(n, 0.0));

	for (int i = 0; i < n; i++) {
		table[i][0] = y[i];
	}

	for (int j = 1; j < n; j++) { 
		for (int i = 0; i < n - j; i++) { // i - строка, вычисляем f[x_i, ..., x_{i+j}]
			table[i][j] = (table[i + 1][j - 1] - table[i][j - 1]) / (x[i + j] - x[i]);
		}
	}

	return table;
}

double P(const std::vector<double>& x, const std::vector<std::vector<double>>& table, double xi) {
	int n = x.size();
	double result = table[0][0];
	double product = 1.0;

	for (int j = 1; j < n; j++) {
		product *= (xi - x[j - 1]);
		result += table[0][j] * product;
	}
	return result;
}

void printDividedDifferences(const std::vector<std::vector<double>>& table) {
	int n = table.size();
	std::cout << "Таблица разделённых разностей:\n";
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n - i; j++) {
			std::cout << std::setw(10) << table[i][j] << " ";
		}
		std::cout << "\n";
	}
}


std::vector<double> equalNodes(int n, double a, double b) {
	std::vector<double> nodes(n);
	for (int i = 0; i < n; i++) {
		nodes[i] = a + i * ((b - a) / (n - 1));
	}
	return nodes;
}

std::vector<double> chebyshevNodes(int n, double a, double b) {
	std::vector<double> nodes(n);
	for (int i = 0; i < n; i++) {
		nodes[i] = cos(PI * (2 * i + 1) / (2 * n));  // Узлы Чебышёва на [-1,1]
		nodes[i] = a + (nodes[i] + 1) * (b - a) / 2; // Преобразование на [a, b]
	}
	return nodes;
}


int main() {
	setlocale(LC_ALL, "ru");
	int a = -2, b = 2;
	int n_values[] = { 5, 10, 15, 20, 30 };

	for (int n : n_values) {
		// Равноотстоящие узлы
		std::vector<double> x_eq = equalNodes(n, a, b);
		std::vector<double> y1_eq(n), y2_eq(n);
		for (int i = 0; i < n; i++) {
			y1_eq[i] = f1(x_eq[i]);
			y2_eq[i] = f2(x_eq[i]);
		}

		std::vector<std::vector<double>> table_f1_eq = dividedDifferences(x_eq, y1_eq);
		std::vector<std::vector<double>> table_f2_eq = dividedDifferences(x_eq, y2_eq);

		double max_error_P1 = 0.0, max_error_P2 = 0.0;


		std::vector<double> P1(n + 1), P2(n + 1);
		std::vector<double> C1(n + 1), C2(n + 1);

		for (double xi = a; xi <= b; xi += 0.01) {
			double P1_n = P(x_eq, table_f1_eq, xi);
			double P2_n = P(x_eq, table_f2_eq, xi);
			max_error_P1 = std::max(max_error_P1, fabs(P1_n - f1(xi)));
			max_error_P2 = std::max(max_error_P2, fabs(P2_n - f2(xi)));
		}

		std::cout << "n = " << n << ":\n";
		std::cout << "max |P1,n(xi) - f1(xi)| = " << max_error_P1 << "\n";
		std::cout << "max |P2,n(xi) - f2(xi)| = " << max_error_P2 << "\n";

		// Чебышёвские узлы
		std::vector<double> x_ch = chebyshevNodes(n, a, b);
		std::vector<double> y1_ch(n), y2_ch(n);
		for (int i = 0; i < n; i++) {
			y1_ch[i] = f1(x_ch[i]);
			y2_ch[i] = f2(x_ch[i]);
		}

		std::vector<std::vector<double>> table_f1_ch = dividedDifferences(x_ch, y1_ch);
		std::vector<std::vector<double>> table_f2_ch = dividedDifferences(x_ch, y2_ch);

		double max_error_C1 = 0.0;
		double max_error_C2 = 0.0;

		for (double xi = -2; xi <= 2; xi += 0.01) {
			double C1_n = P(x_ch, table_f1_ch, xi);
			double C2_n = P(x_ch, table_f2_ch, xi);
			max_error_C1 = std::max(max_error_C1, fabs(C1_n - f1(xi)));
			max_error_C2 = std::max(max_error_C2, fabs(C2_n - f2(xi)));
		}

		std::cout << "max |C1,n(xi) - f1(xi)| = " << max_error_C1 << "\n";
		std::cout << "max |C2,n(xi) - f2(xi)| = " << max_error_C2 << "\n";	
	}

	return 0;
}