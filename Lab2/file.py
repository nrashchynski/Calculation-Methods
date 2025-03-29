import numpy as np
import matplotlib.pyplot as plt
import re


# Первая функция f1
def f1(x):
    return np.sin(2 * x) * np.log(x + 5)


# Вторая функция f2
def f2(x):
    return np.sqrt(2 * np.abs(x) + x ** 2)


# Интерполяционный многочлен Ньютона
def newton_interpolation(x, y, xi):
    n = len(x)
    coeffs = np.copy(y)

    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            coeffs[i] = (coeffs[i] - coeffs[i - 1]) / (x[i] - x[i - j])

    result = coeffs[0]
    product = 1.0

    for i in range(1, n):
        product *= (xi - x[i - 1])
        result += coeffs[i] * product

    return result


# Равноотстоящие узлы
def equidistant_nodes(n):
    return np.linspace(-2, 2, n)


# Чебышёвские узлы
def chebyshev_nodes(n, a=-2, b=2):
    nodes = np.cos((2 * np.arange(n) + 1) * np.pi / (2 * n))
    return a + (nodes + 1) * (b - a) / 2


# Очистка имени файла
def sanitize_filename(filename):
    return re.sub(r'[^\w\-_\.]', '_', filename)


# Основная функция
def main():
    n_values = [3, 10, 20]
    x_plot = np.linspace(-2, 2, 400)

    functions = [(f1, "f1(x) = sin(2x) * ln(x+5)"),
                 (f2, "f2(x) = sqrt(2|x| + x^2)")]

    for func, label in functions:
        f_values = func(x_plot)

        for n in n_values:
            # Равноотстоящие узлы
            x_eq = equidistant_nodes(n)
            y_eq = func(x_eq)
            P_n_values = [newton_interpolation(x_eq, y_eq, xi) for xi in x_plot]

            # Чебышёвские узлы
            x_ch = chebyshev_nodes(n)
            y_ch = func(x_ch)
            C_n_values = [newton_interpolation(x_ch, y_ch, xi) for xi in x_plot]

            # Построение графиков
            plt.figure(figsize=(8, 6))
            plt.title(f"Интерполяция {label} для n = {n}")
            plt.xlabel("x")
            plt.ylabel("y")

            plt.plot(x_plot, f_values, label=label, color='red', linewidth=2)
            plt.plot(x_plot, P_n_values, label=f"P_n({n}), равноотстоящие", linestyle='--', color='blue')
            plt.plot(x_plot, C_n_values, label=f"C_n({n}), Чебышёвские", linestyle='--', color='green')

            plt.legend()
            plt.grid()
            filename = sanitize_filename(f"interpolation_{label}_n_{n}.png")
            plt.savefig(filename)
            plt.close()


if __name__ == "__main__":
    main()
