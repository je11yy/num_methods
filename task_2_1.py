# x^3 + x^2 - x - 0.5 = 0

import math
import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, diff, solve, Abs

function = lambda x : x ** 3 + x ** 2 - x - 0.5
dfunction = lambda x: 3 * (x ** 2) + 2 * x - 1
d2function = lambda x: 6 * x + 2
g = lambda x: math.sqrt(x + 0.5 - x ** 3)
dg = lambda x: (1.0 - 3 * (x ** 2)) / math.sqrt(x + 0.5 - x ** 3)

# метод простых итераций
def iterations_method(x0, a, b, epsilon):
    def get_q(a, b, steps=1000):
        return max(dg(a), dg(b))
    q = get_q(a, b)
    if q >= 1:
        raise ValueError("Условие сходимости метода простых итераций не выполняется (q >= 1).")

    x = x0
    iteration = 0
    while True:
        x_new = g(x)
        iteration += 1
        if (q / (1.0 - q)) * abs(x_new - x) < epsilon:
            break
        x = x_new
    return x, iteration

# метод дихотомии
def dichotomy_method(a, b, epsilon):
    if function(a) * function(b) >= 0:
        raise ValueError("Метод дихотомии: на концах отрезка нет смены знака функции")
    iteration = 0
    while (b - a) / 2.0 > epsilon:
        midpoint = (a + b) / 2.0
        iteration += 1
        if function(midpoint) == 0:
            return midpoint, iteration
        elif function(a) * function(midpoint) < 0:
            b = midpoint
        else:
            a = midpoint
    return (a + b) / 2.0, iteration

# поиск подходящего приближения
def find_approx(a, b, epsilon):
    if (function(a) * function(b) >= 0):
        raise ValueError("Невозможно найти приближение на данном отрезке")
    x = a
    while (x <= b):
        if (function(x) * d2function(x) > 0):
            return x
        x += epsilon
    raise ValueError("Невозможно найти приближение на данном отрезке")

# метод Ньютона
def newton_method(a, b, epsilon):
    x = find_approx(a, b, epsilon)
    iteration = 0
    while True:
        x_next = x - function(x) / dfunction(x)
        iteration += 1
        if abs(x_next - x) < epsilon:
            break
        x = x_next
    return x, iteration

def check(a, b, epsilon):
    x = a
    while (x <= b):
        y = g(x)
        print(y)
        if not (y >= a and y <= b):
            raise ValueError("Метод секущих не сходится [1]")
        if not (np.abs(dg(x)) < 1):
            raise ValueError("Метод секущих не сходится [2]")
        x += epsilon

def secant_method(x0, x1, epsilon):
    # check(x0, x1, epsilon)
    iteration = 0
    while True:
        f_x0 = function(x0)
        f_x1 = function(x1)
        x2 = x1 - f_x1 * (x1 - x0) / (f_x1 - f_x0)
        iteration += 1
        if abs(x2 - x1) < epsilon:
            break
        x0, x1 = x1, x2
    return x2, iteration

epsilon = 1e-6
a, b = 0.5, 1.0
x0 = 1.0

formatter = lambda arg: format(arg, "0.6f")

root_simple, iter_simple = iterations_method(x0, a, b, epsilon)
print(f"Метод простой итерации: корень = ", formatter(root_simple), f", итерации = {iter_simple}")

root_bisection, iter_bisection = dichotomy_method(a, b, epsilon)
print(f"Метод дихотомии: корень = ", formatter(root_bisection), f", итерации = {iter_bisection}")

root_newton, iter_newton = newton_method(a, b, epsilon)
print(f"Метод Ньютона: корень = ", formatter(root_newton), f", итерации = {iter_newton}")

root_secant, iter_secant = secant_method(a, b, epsilon)
print(f"Метод секущих: корень = ", formatter(root_secant), f", итерации = {iter_secant}")

# Рисуем график
x = np.linspace(-2, 2, 400)
y = function(x)

plt.style.use('dark_background')

plt.figure(figsize=(10, 6), facecolor='black')
plt.plot(x, y, label='f(x) = x^3 + x^2 - x - 0.5', color='cyan')

plt.axhline(0, color='white', lw=2)  # Ось X
plt.axvline(0, color='white', lw=2)  # Ось Y

# Добавление стрелок на оси
arrowprops = dict(width=1, headwidth=5, headlength=8, color='white')
plt.annotate('', xy=(2, 0), xytext=(1.8, 0), arrowprops=arrowprops)  # Стрелка на оси X
plt.annotate('', xy=(0, 5), xytext=(0, 4), arrowprops=arrowprops)  # Стрелка на оси Y

plt.title('График функции', color='white')
plt.xlabel('x', color='white')
plt.ylabel('y', color='white')
# plt.grid(color='gray', linestyle='--', linewidth=0.5)
# plt.grid(True)

plt.xlim(-3, 3)
plt.ylim(-3, 10)

plt.legend(facecolor='black', edgecolor='white', fontsize=12, loc='upper left')

# Добавление текста с результатами
text_str = (f"             Корни методов\n"
             f"Метод простых итераций: {root_simple:.6f}\n"
             f"Метод дихотомии: {root_bisection:.6f}\n"
             f"Метод Ньютона: {root_newton:.6f}\n"
             f"Метод секущих: {root_secant:.6f}")

# Размещение текста на графике
plt.text(-2.8, 4, text_str, fontsize=12, color='white', bbox=dict(facecolor='black'))

plt.show()