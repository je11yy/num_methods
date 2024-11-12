import numpy as np
import matplotlib.pyplot as plt

x = np.array([1.0, 1.9, 2.8, 3.7, 4.6, 5.5])
y = np.array([2.4142, 1.0818, 0.50953, 0.11836, -0.24008, -0.66818])

def build_normal_system(x, y, degree):
    # Инициализируем матрицу нормальной системы и вектор правой части
    A = np.zeros((degree + 1, degree + 1))
    B = np.zeros(degree + 1)
    
    # Заполняем матрицу A и вектор B по формуле нормальной системы МНК
    for k in range(degree + 1):
        for i in range(degree + 1):
            A[k, i] = np.sum(x ** (k + i))
        B[k] = np.sum(y * (x ** k))
    
    return A, B

def fit_polynomial(x, y, degree):
    A, B = build_normal_system(x, y, degree)
    # Решаем систему линейных уравнений для коэффициентов многочлена
    coeffs = np.linalg.solve(A, B)
    return coeffs

def get_polynomial(x, y, degree):
    coeffs = fit_polynomial(x, y, degree)
    poly = np.poly1d(coeffs[::-1]) # обратный порядок коэффициентов
    y_approx = poly(x)
    error = np.sum((y - y_approx)**2)
    return poly, error

# a) Многочлен первой степени
poly_1, error_1 = get_polynomial(x, y, 1)

# b) Многочлен второй степени
poly_2, error_2 = get_polynomial(x, y, 2)

print(f"Сумма квадратов ошибок для многочлена 1-й степени: {error_1:.5f}")
print(f"Сумма квадратов ошибок для многочлена 2-й степени: {error_2:.5f}")
print(f"Приближающий многочлен 1 степени:\n")
print(poly_1)
print(f"\nПриближающий многочлен 2 степени:\n")
print(poly_2)

# Построение графиков
x_plot = np.linspace(min(x), max(x), 100)
y_poly_1 = poly_1(x_plot)
y_poly_2 = poly_2(x_plot)

plt.figure(figsize=(10, 6))
plt.plot(x, y, 'o', label='Табличная функция', color='black')
plt.plot(x_plot, y_poly_1, label='Приближающий многочлен 1-й степени', color='blue')
plt.plot(x_plot, y_poly_2, label='Приближающий многочлен 2-й степени', color='red')

plt.xlabel('x')
plt.ylabel('y')
plt.title('Аппроксимация табличной функции')
plt.legend()
plt.grid()
plt.show()