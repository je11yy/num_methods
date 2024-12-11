import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Исходная задача Коши: y'' = f(x, y, y')
def f(x, y, dy):
    return (x * (x**2 - 1.0) * dy + (x**2 + 1.0) * y) / (x**2)

# разбиваем на систему двух диф уравнений:
# y' = z
# z' = x * (x^2 - 1) * z + (x^2 + 1) * y
# y(1) = 1 + exp(1/2)
# z = 2* exp(1 / 2) - 1
# [1, 2]
# h = 0.1

# Точное решение
def exact_solution(x):
    return (1.0 / x) * (1.0 + np.exp((x**2) / 2.0))

# Начальные условия
x0, y0, dy0 = 1.0, 1 + np.exp(1 / 2), 2 * np.exp(1 / 2) - 1
x_end = 2.0
h = 0.1

# Метод Эйлера (явный)
def euler_method(f, x0, y0, dy0, h, x_end):
    x_vals = np.arange(x0, x_end + h, h)
    y_vals = [y0]
    dy_vals = [dy0]
    
    for i in range(1, len(x_vals)):
        x = x_vals[i - 1]
        y = y_vals[-1]
        dy = dy_vals[-1]
        
        ddy = f(x, y, dy)
        y_vals.append(y + h * dy)
        dy_vals.append(dy + h * ddy)
    
    return x_vals, np.array(y_vals)

# Метод Рунге-Кутты 4-го порядка
def runge_kutta_method(f, x0, y0, dy0, h, x_end):
    x_vals = np.arange(x0, x_end + h, h)
    y_vals = [y0]
    dy_vals = [dy0]
    
    for i in range(1, len(x_vals)):
        x = x_vals[i - 1]
        y = y_vals[-1]
        dy = dy_vals[-1]
        
        k1 = h * dy
        l1 = h * f(x, y, dy)
        
        k2 = h * (dy + l1 / 2)
        l2 = h * f(x + h / 2, y + k1 / 2, dy + l1 / 2)
        
        k3 = h * (dy + l2 / 2)
        l3 = h * f(x + h / 2, y + k2 / 2, dy + l2 / 2)
        
        k4 = h * (dy + l3)
        l4 = h * f(x + h, y + k3, dy + l3)
        
        y_next = y + (k1 + 2 * k2 + 2 * k3 + k4) / 6
        dy_next = dy + (l1 + 2 * l2 + 2 * l3 + l4) / 6
        y_vals.append(y_next)
        dy_vals.append(dy_next)
    
    return x_vals, np.array(y_vals)

# Метод Адамса 4-го порядка
def adams_method(f, x0, y0, dy0, h, x_end):
    # Используем метод Рунге-Кутты для начальных 3-х точек
    x_vals, y_vals = runge_kutta_method(f, x0, y0, dy0, h, x0 + 3 * h)
    dy_vals = [dy0]
    x_vals = list(x_vals)
    y_vals = list(y_vals)

    for i in range(1, len(x_vals)):
        x = x_vals[i - 1]
        y = y_vals[i - 1]
        dy = dy_vals[-1]
        ddy = f(x, y, dy)
        dy_vals.append(dy + h * ddy)

    # Итерации метода Адамса
    while x_vals[-1] < x_end:
        y3, y2, y1, y0 = y_vals[-4], y_vals[-3], y_vals[-2], y_vals[-1]
        dy3, dy2, dy1, dy0 = dy_vals[-4], dy_vals[-3], dy_vals[-2], dy_vals[-1]
        
        x = x_vals[-1]
        ddy0 = f(x, y0, dy0)
        
        # Метод Адамса 4-го порядка
        y_next = y0 + h * (55 * dy0 - 59 * dy1 + 37 * dy2 - 9 * dy3) / 24
        dy_next = dy0 + h * (55 * ddy0 - 59 * f(x - h, y1, dy1) + 
                             37 * f(x - 2 * h, y2, dy2) - 9 * f(x - 3 * h, y3, dy3)) / 24
        
        x_vals.append(x + h)
        y_vals.append(y_next)
        dy_vals.append(dy_next)
    
    return np.array(x_vals), np.array(y_vals)


# Оценка погрешности методом Рунге-Ромберга
def runge_romberg(y_h, y_h2, h, p):
    return (y_h2 - y_h) / (2**p - 1)

# Решение задачи
x_euler, y_euler = euler_method(f, x0, y0, dy0, h, x_end)
x_rk, y_rk = runge_kutta_method(f, x0, y0, dy0, h, x_end)
x_adams, y_adams = adams_method(f, x0, y0, dy0, h, x_end)

# Точное решение
num = int((x_end - x0) / h)
x_exact = np.linspace(x0, x_end, num + 1)
y_exact = exact_solution(x_exact)

# Оценка погрешности методом Рунге-Ромберга
h2 = h / 2
_, y_rk_h2 = runge_kutta_method(f, x0, y0, dy0, h2, x_end)
error_runge_romberg = runge_romberg(y_rk, y_rk_h2[::2], h, 4)

# Подготовка данных для таблицы
data = {
    "x": x_euler,
    "Точное значение": exact_solution(x_euler),
    "Метод Эйлера": y_euler,
    "Метод Рунге-Кутты": y_rk,
    "Метод Адамса": y_adams,
    "Погрешность Эйлера": np.abs(y_euler - exact_solution(x_euler)),
    "Погрешность Рунге-Кутты": np.abs(y_rk - exact_solution(x_euler)),
    "Погрешность Адамса": np.abs(y_adams - exact_solution(x_euler))
}

# Создание таблицы
df = pd.DataFrame(data)

# Вывод таблицы
print("Таблица результатов:")
print(df.to_string(index=False))

# Построение графиков
plt.figure(figsize=(10, 6))
plt.plot(x_exact, y_exact, label="Точное решение", color="black", linestyle="--")
plt.plot(x_euler, y_euler, label="Метод Эйлера", marker="o")
plt.plot(x_rk, y_rk, label="Метод Рунге-Кутты 4-го порядка", marker="s")
plt.plot(x_adams, y_adams, label="Метод Адамса 4-го порядка", marker="^")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.title("Решение задачи Коши")
plt.grid()
plt.show()
