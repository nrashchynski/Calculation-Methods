import numpy as np
from scipy.integrate import quad
from sympy import symbols, integrate, legendre_poly, diff, solve


def gauss_quadrature(f, a, b, k=5):
    x = symbols('x')
    Pk = legendre_poly(k, x)
    nodes = solve(Pk, x)
    nodes = [float(n.evalf()) for n in nodes]

    weights = []
    for xi in nodes:
        dPk = diff(Pk, x).subs(x, xi)
        wi = 2 / ((1 - xi ** 2) * dPk ** 2)
        weights.append(float(wi.evalf()))

    scaled_nodes = [0.5 * (b - a) * xi + 0.5 * (a + b) for xi in nodes]
    scaled_weights = [0.5 * (b - a) * wi for wi in weights]

    integral = sum(w * f(xi) for xi, w in zip(scaled_nodes, scaled_weights))

    return integral, scaled_nodes, scaled_weights


def f(x):
    return (x - 1) ** 2 / (x ** 2 + np.exp(x) + 1)


def exact_integral():
    return 2 + np.log(2 / (5 + np.exp(2)))


a, b = 0, 2
k = 5

approx, nodes, weights = gauss_quadrature(f, a, b, k)
exact = exact_integral()
error = abs(approx - exact)

print(f"=== Метод Гаусса с {k} узлами ===")
print(f"Приближенное значение: {approx:.12f}")
print(f"Точное значение:      {exact:.12f}")
print(f"Абсолютная погрешность: {error:.2e}")
print("\nУзлы и веса:")
for xi, wi in zip(nodes, weights):
    print(f"Узел: {xi:.6f}, Вес: {wi:.6f}")

scipy_value, _ = quad(f, a, b)
print(f"\nПроверка (scipy.quad): {scipy_value:.12f}")
