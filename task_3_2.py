import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

X = [1.0, 1.9, 2.8, 3.7, 4.6]
Y = [2.4142, 1.01818, 0.50953, 0.11836, -0.24008]

def cubic_spline(X, Y):
    n = len(X)
    
    # Система для решения коэффициентов
    h = np.diff(X)  # Разность между соседними точками X
    alpha = np.zeros(n - 1)
    
    # Строим правую часть системы (alpha)
    for i in range(1, n - 1):
        alpha[i] = (3 / h[i]) * (Y[i+1] - Y[i]) - (3 / h[i-1]) * (Y[i] - Y[i-1])
    
    # Решаем систему для коэффициентов S''(X)
    A = np.zeros((n, n))
    b = np.zeros(n)
    
    A[0, 0] = 1
    A[n-1, n-1] = 1
    for i in range(1, n - 1):
        A[i, i-1] = h[i-1]
        A[i, i] = 2 * (h[i-1] + h[i])
        A[i, i+1] = h[i]
        b[i] = alpha[i]
    
    # Решаем систему A * c = b для c (вторые производные на узловых точках)
    c = np.linalg.solve(A, b)
    
    # Находим коэффициенты для полиномов на интервалах
    a = Y[:-1]
    b = (np.array(Y[1:]) - np.array(Y[:-1])) / h - h * (2 * c[:-1] + c[1:]) / 3
    d = (c[1:] - c[:-1]) / (3 * h)
    
    # Полиномы для каждого интервала
    x = sp.symbols('x')
    polynomials = []
    for i in range(n - 1):
        polynomials.append(sp.expand((a[i] + b[i] * (x - X[i]) + c[i] * (x - X[i])**2 + d[i] * (x - X[i])**3)))
    
    return polynomials

# Вычисление кубического сплайна
polynomials = cubic_spline(X, Y)

# Печать полиномов для каждого интервала
i = 0
for poly in polynomials:
    print(f"\nX=[{X[i]}, {X[i+1]}]:\n")
    sp.pprint(poly.evalf(4))
    i += 1

X_star = 2.66666667

for i in range(len(X) - 1):
    if X_star >= X[i] and X_star <= X[i + 1]:
        print(f"\nНайденное значение: ", polynomials[i].subs(sp.symbols('x'), X_star).evalf(5))


# Построение графика
x = sp.symbols('x')
x_vals = np.linspace(min(X), max(X), 500)
y_vals = []

# Вычисление значений каждого полинома на интервале
for i in range(len(X) - 1):
    interval_x = np.linspace(X[i], X[i + 1], 100)
    interval_y = [polynomials[i].subs(x, xi).evalf() for xi in interval_x]
    y_vals.extend(interval_y)
    plt.plot(interval_x, interval_y, label=f"Полином {i+1}")

# Отображение исходных точек
plt.scatter(X, Y, color="red", label="Исходные точки")

# Настройки графика
plt.title("Кубический сплайн")
plt.xlabel("X")
plt.ylabel("Y")
plt.legend()
plt.grid(True)
plt.show()