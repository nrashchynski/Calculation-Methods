import numpy as np
import matplotlib.pyplot as plt


def exact_solution(x):
    return np.tan(x) - x


def f(x, u):
    return (u + x) ** 2


def newton_method(x, u_prev, h, tol=1e-8, max_iter=100):
    u = u_prev
    for _ in range(max_iter):
        fu = u - u_prev - h / 2 * (f(x, u_prev) + f(x + h, u))
        dfu = 1 - h / 2 * 2 * (u + x + h)
        delta = fu / dfu
        u -= delta
        if abs(delta) < tol:
            break
    return u


def implicit_trapezoidal(a, b, h, u0):
    n = int((b - a) / h)
    x_values = np.linspace(a, b, n + 1)
    u_values = np.zeros(n + 1)
    u_values[0] = u0
    for i in range(n):
        x = x_values[i]
        u_values[i + 1] = newton_method(x, u_values[i], h)
    return x_values, u_values


def runge_kutta_2(a, b, h, u0):
    n = int((b - a) / h)
    x_values = np.linspace(a, b, n + 1)
    u_values = np.zeros(n + 1)
    u_values[0] = u0
    for i in range(n):
        x = x_values[i]
        k1 = f(x, u_values[i])
        k2 = f(x + h / 2, u_values[i] + h / 2 * k1)
        u_values[i + 1] = u_values[i] + h * k2
    return x_values, u_values


def compute_max_error(x_values, y_values):
    exact = exact_solution(x_values)
    errors = np.abs(exact - y_values)
    return np.max(errors)


if __name__ == "__main__":
    a, b = 0.0, 1.0
    h1, h2 = 0.1, 0.05
    u0 = 0.0

    x_trap_h1, u_trap_h1 = implicit_trapezoidal(a, b, h1, u0)
    x_trap_h2, u_trap_h2 = implicit_trapezoidal(a, b, h2, u0)

    x_rk_h1, u_rk_h1 = runge_kutta_2(a, b, h1, u0)
    x_rk_h2, u_rk_h2 = runge_kutta_2(a, b, h2, u0)

    error_trap_h1 = compute_max_error(x_trap_h1, u_trap_h1)
    error_trap_h2 = compute_max_error(x_trap_h2, u_trap_h2)
    error_rk_h1 = compute_max_error(x_rk_h1, u_rk_h1)
    error_rk_h2 = compute_max_error(x_rk_h2, u_rk_h2)

    p = 2
    R_trap = abs(error_trap_h1 - error_trap_h2) / (2 ** p - 1)
    R_rk = abs(error_rk_h1 - error_rk_h2) / (2 ** p - 1)

    print("Результаты для неявного метода трапеций:")
    print(f"Шаг h1 = {h1}, Максимальная ошибка = {error_trap_h1}")
    print(f"Шаг h2 = {h2}, Максимальная ошибка = {error_trap_h2}")
    print(f"Оценка по правилу Рунге: {R_trap}")
    print()

    print("Результаты для метода Рунге-Кутты (порядок 2):")
    print(f"Шаг h1 = {h1}, Максимальная ошибка = {error_rk_h1}")
    print(f"Шаг h2 = {h2}, Максимальная ошибка = {error_rk_h2}")
    print(f"Оценка по правилу Рунге: {R_rk}")

    # Построение графиков
    plt.figure(figsize=(12, 6))

    # График для шага h1
    plt.subplot(1, 2, 1)
    plt.plot(x_trap_h1, u_trap_h1, label="Неявный метод трапеций (h1)", marker='o', linestyle='--')
    plt.plot(x_rk_h1, u_rk_h1, label="Метод Рунге-Кутты (h1)", marker='x', linestyle='-.')
    plt.plot(x_trap_h1, exact_solution(x_trap_h1), label="Точное решение", linestyle='-', color='black')
    plt.title("Сравнение решений (шаг h1)")
    plt.xlabel("x")
    plt.ylabel("u(x)")
    plt.legend()
    plt.grid()

    # График для шага h2
    plt.subplot(1, 2, 2)
    plt.plot(x_trap_h2, u_trap_h2, label="Неявный метод трапеций (h2)", marker='o', linestyle='--')
    plt.plot(x_rk_h2, u_rk_h2, label="Метод Рунге-Кутты (h2)", marker='x', linestyle='-.')
    plt.plot(x_trap_h2, exact_solution(x_trap_h2), label="Точное решение", linestyle='-', color='black')
    plt.title("Сравнение решений (шаг h2)")
    plt.xlabel("x")
    plt.ylabel("u(x)")
    plt.legend()
    plt.grid()

    plt.tight_layout()
    plt.show()
    
