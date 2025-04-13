import matplotlib.pyplot as plt
import numpy as np

x_vals = []
f_vals = []
s_vals = []

with open('C:/Users/Lenovo/source/CMLabs/CMLab3/output.txt', 'r') as file:
    for line in file:
        parts = line.strip().split()
        if len(parts) == 3:
            x, f, s = map(float, parts)
            x_vals.append(x)
            f_vals.append(f)
            s_vals.append(s)

x_vals = np.array(x_vals)
f_vals = np.array(f_vals)
s_vals = np.array(s_vals)

# Построение графика
plt.figure(figsize=(10, 6))
plt.plot(x_vals, f_vals, label='Исходная функция f(x)', color='blue')
plt.plot(x_vals, s_vals, label='Интерполяционный сплайн S(x)', color='red', linestyle='--')
plt.title('Интерполяция кубическим сплайном')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid(True)
plt.savefig("spline_plot.png")
plt.show()
