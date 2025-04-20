#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

const int n = 15;
const double a = -2.0;
const double b = 2.0;
const double h = (double)(b - a) / n;

double f(double x) {
	return sin(2 * x) * log(x + 5);
}

double df(double x) {
	return 2 * cos(2 * x) * log(x + 5) + sin(2 * x) / (x + 5);
}

double d2f(double x) {
	double term1 = -4 * sin(2 * x) * log(x + 5);
	double term2 = 4 * cos(2 * x) / (x + 5);
	double term3 = -sin(2 * x) / pow(x + 5, 2);
	return term1 + term2 + term3;
}

std::vector<double> runThrowMethod(
	const std::vector<double>& A,
	const std::vector<double>& B,
	const std::vector<double>& C,
	const std::vector<double>& F
) {
	int n = C.size() - 1;
	std::vector<double> alpha(n + 2), beta(n + 2), M(n + 1);

	alpha[1] = B[0] / C[0];
	beta[1] = F[0] / C[0];

	for (int i = 1; i < n; i++) {
		double denom = C[i] - A[i] * alpha[i];
		alpha[i + 1] = B[i] / denom;
		beta[i + 1] = (F[i] + A[i] * beta[i]) / denom;
	}

	beta[n + 1] = (F[n] + A[n] * beta[n]) / (C[n] - A[n] * alpha[n]);
	M[n] = beta[n + 1];

	for (int i = n - 1; i >= 0; i--) {
		M[i] = alpha[i + 1] * M[i + 1] + beta[i + 1];
	}

	return M;
}

double S3(const std::vector<double>& M, const std::vector<double>& x, const std::vector<double>& h, const std::vector<double>& y, int n, double X) {
	int k;
	for (int i = 1; i <= n; ++i) {
		if (X >= x[i - 1] && X <= x[i]) {
			k = i;
			break;
		}
	}
	return M[k - 1] * std::pow(x[k] - X, 3) / (6 * h[k]) +
		M[k] * std::pow(X - x[k - 1], 3) / (6 * h[k]) +
		(y[k - 1] - M[k - 1] * h[k] * h[k] / 6) * (x[k] - X) / h[k] +
		(y[k] - M[k] * h[k] * h[k] / 6) * (X - x[k - 1]) / h[k];
}

int main() {
	setlocale(LC_ALL, "ru");
	std::vector<double>x(n + 1), y(n + 1), h_arr(n + 1);
	std::vector<double> A(n + 1), B(n + 1), C(n + 1), F(n + 1);
	std::vector<double> alpha(n + 2), beta(n + 2), M(n + 1);

	for (int i = 0; i <= n; ++i) {
		x[i] = a + i * h;
		y[i] = f(x[i]);
	}

	for (int i = 1; i <= n; ++i) {
		h_arr[i] = x[i] - x[i - 1];
	}

	C[0] = 1;
	C[n] = h_arr[n] / 3;
	for (int i = 1; i < n; i++) {
		C[i] = (h_arr[i] + h_arr[i + 1]) / 3;
	}

	A[n] = h_arr[n] / 6;
	B[0] = 0;
	for (int i = 1; i < n; i++) {
		A[i] = h_arr[i] / 6;
		B[i] = h_arr[i + 1] / 6;
	}

	F[0] = d2f(x[0]);  // S''(x0) = f''(x0)
	F[n] = df(x[n]) - (y[n] - y[n-1]) / h_arr[n];  // S'(xn) = f'(xn)
	for (int i = 1; i < n; i++) {
		F[i] = (y[i + 1] - y[i]) / h_arr[i + 1] - (y[i] - y[i - 1]) / h_arr[i];
	}

	M = runThrowMethod(A, B, C, F);

	std::ofstream fout("output.txt");
	double max_diff = 0;
	for (int i = 0; i <= 100; i++) {
		double X = a + i * (b - a) / 100;
		double fx = f(X);
		double sx = S3(M, x, h_arr, y, n, X);
		fout << X << " " << fx << " " << sx << std::endl;
		max_diff = std::max(max_diff, fabs(fx - sx));
	}

	fout.close();
	std::cout << "Максимальная погрешность: " << max_diff << std::endl;

	return 0;
}
