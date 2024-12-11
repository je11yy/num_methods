
# y = (3x + 4) / (2x + 7)
X0 = -2.0
Xk = 2.0
h1 = 1.0
h2 = 0.5
y = lambda x: (3.0 * x + 4.0) / (2.0 * x + 7.0)

def rectangle_method(start, end, step):
    n = int((end - start) / step)
    x = start
    result = 0
    for _ in range(n):
        prev_x = x
        x += step
        result += y((prev_x + x) / 2)
    return result * step

def trapezoid_method(start, end, step):
    n = int((end - start) / step)
    x = start
    result = 0
    for _ in range(n):
        prev_x = x
        x += step
        result += (y(x) + y(prev_x))
    return result * 0.5 * step

def sympson_method(start, end, step):
    n = int((end - start) / step)
    x = start
    result = y(x) + y(end)

    for i in range(1, n):
        x += step
        if i % 2 == 1:
            result += 4 * y(x)
        else:
            result += 2 * y(x)

    return result * (step / 3.0)

def runge_romberg(I_h1, I_h2, p, h1, h2):
    return (I_h2 - I_h1) / ((h1 / h2) ** p - 1.0)

# Вычисляем значения для h1 и h2 для каждого метода
I_rect_h1 = rectangle_method(X0, Xk, h1)
I_rect_h2 = rectangle_method(X0, Xk, h2)

I_trap_h1 = trapezoid_method(X0, Xk, h1)
I_trap_h2 = trapezoid_method(X0, Xk, h2)

I_simp_h1 = sympson_method(X0, Xk, h1)
I_simp_h2 = sympson_method(X0, Xk, h2)

# Оценка погрешности
error_rect = runge_romberg(I_rect_h1, I_rect_h2, 1, h1, h2)
error_trap = runge_romberg(I_trap_h1, I_trap_h2, 2, h1, h2)
error_simp = runge_romberg(I_simp_h1, I_simp_h2, 4, h1, h2)

print(f"\nМетод прямоугольников\n")
print(f"C шагом 1: {I_rect_h1}")
print(f"C шагом 0.5: {I_rect_h2}")
print(f"Погрешность: {error_rect}")
print(f"Приближенное значение: {I_rect_h2 + error_rect}\n")

print(f"Метод трапеций\n")
print(f"C шагом 1: {I_trap_h1}")
print(f"C шагом 0.5: {I_trap_h2}")
print(f"Погрешность: {error_trap}")
print(f"Приближенное значение: {I_trap_h2 + error_trap}\n")

print(f"Метод Симпсона\n")
print(f"C шагом 1: {I_simp_h1}")
print(f"C шагом 0.5: {I_simp_h2}")
print(f"Погрешность: {error_simp}")
print(f"Приближенное значение: {I_simp_h2 + error_simp}\n")

print(f"Точное значение: 1.777\n")