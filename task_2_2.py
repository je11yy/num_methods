import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# x - cosy = 1
# y - lg(x + 1) = 1

def f(x):
    if x[0] <= -1:
        raise ValueError(f"Невозможно вычислить логарифм: x[0] + 1 должно быть > 0, текущее значение x[0]: {x[0]}")
    return [x[0] - np.cos(x[1]) - 1, x[1] - np.log10(x[0] + 1) - 1]

f1 = lambda x: np.arccos(x - 1) if -1 <= x - 1 <= 1 else np.nan
f2 = lambda x: np.log10(x + 1) + 1 if x > -1 else np.nan

def g(x):
    if x[0] <= -1:
        raise ValueError(f"Невозможно вычислить логарифм: x[0] + 1 должно быть > 0, текущее значение x[0]: {x[0]}")
    return [np.cos(x[1]) + 1, np.log10(x[0] + 1) + 1]

def dxfunction(x):
    return [1, -(1 / (x * np.log(10)))]

def dyfunction(y):
    return [np.sin(y), 1]

def iterations_method(x, epsilon):
    iteration = 0
    while True:
        try:
            x_new = g(x)
        except ValueError as e:
            print("Ошибка в методе итераций:", e)
            return []
        delta = np.array(x_new) - np.array(x)
        if np.linalg.norm(delta) <= epsilon:
            break
        x = x_new
        iteration += 1
    return x, iteration

def newton_method(x, epsilon):
    iteration = 0
    while True:
        try:
            dx = dxfunction(x[0])
            dy = dyfunction(x[1])
            J = np.array([[dx[0], dy[0]], [dx[1], dy[1]]])
            if np.linalg.det(J) == 0:
                raise ValueError(f"Матрица Якобиана не вырождена")
            F = f(x)
        except ValueError as e:
            print("Ошибка в методе Ньютона:", e)
            return []
        
        x_new = x - np.linalg.inv(J) @ F
        if np.linalg.norm(np.array(x_new - x)) < epsilon:
            break

        x = x_new
        iteration += 1
    return x, iteration

# Начальные значения
x = [1.25, 1.25]
epsilon = 1e-10

# Простая итерация
result_simple, iteration = iterations_method(x, epsilon)
if len(result_simple) != 0:
    print(f"Метод простой итерации: x = {result_simple[0]:.6f}, y = {result_simple[1]:.6f}")
    print(f"Количество итераций: {iteration}")

# Метод Ньютона
result_newton, iteration = newton_method(x, epsilon)
if len(result_newton) != 0:
    print(f"Метод Ньютона: x = {result_newton[0]:.6f}, y = {result_newton[1]:.6f}")
    print(f"Количество итераций: {iteration}")

solution = fsolve(f, x)
print(f"Решение с использованием fsolve: x = {solution[0]:.6f}, y = {solution[1]:.6f}")

# Графическое представление
x_values = np.linspace(-1, 3, 400)
y1_values = [f1(x) for x in x_values]
y2_values = [f2(x) for x in x_values]

plt.figure(figsize=(10, 6))
plt.plot(x_values, y1_values, label='x - cos(y) = 1', color='blue')
plt.plot(x_values, y2_values, label='y - log(x + 1) = 1', color='red')
plt.scatter(x[0], x[1], color='green', label='Начальное приближение', zorder=5)
plt.xlabel('x')
plt.ylabel('y')
plt.title('График системы уравнений')
plt.axhline(0, color='black', lw=2)
plt.axvline(0, color='black', lw=2)

plt.legend()
plt.grid()
plt.show()