import math
from scipy.interpolate import lagrange
import sympy as sp
import numpy as np

Xa = [math.pi / 8, 2 * math.pi / 8, 3 * math.pi / 8, 4 * math.pi / 8]
Xb = [math.pi / 8, 5 * math.pi / 16, 3 * math.pi / 8, math.pi / 2]
X_star = math.pi / 3

# Лагранж
# Ньютон
# Погрешность интерполяции

f = lambda x: 1.0 / math.tan(x)
x_sym = sp.symbols('x')
f_sym = 1 / sp.tan(x_sym)

def get_Y(X):
    return [f(x) for x in X]

def Lagrange_polynomial(X, Y):
    def composition(X, i, n):
        result = 1
        for j in range(n):
            if (i != j):
                result *= (x_sym - X[j]) / (X[i] - X[j])
        return result
    n = len(X)
    lagrange_poly = sum(Y[i] * composition(X, i, n) for i in range(n))
    return sp.simplify(lagrange_poly).evalf(4)

def Newton_polynomial(X, Y):
    # Функция для вычисления разделенных разностей
    def divided_differences(X, Y):
        n = len(X)
        # Создаем таблицу разделенных разностей и заполняем первую колонку значениями Y
        coef = np.zeros([n, n])
        coef[:,0] = Y

        # Вычисляем разделенные разности
        for j in range(1, n):
            for i in range(n - j):
                coef[i][j] = (coef[i+1][j-1] - coef[i][j-1]) / (X[i+j] - X[i])
        
        # Возвращаем только первую строку, содержащую коэффициенты для многочлена Ньютона
        return coef[0, :]

    n = len(X)
    coefficients = divided_differences(X, Y)
    newton_poly = coefficients[0]
    term = 1

    for i in range(1, n):
        term *= (x_sym - X[i-1])
        newton_poly += coefficients[i] * term

    return newton_poly.evalf(4)

lagrange_a = Lagrange_polynomial(Xa, get_Y(Xa))
lagrange_b = Lagrange_polynomial(Xb, get_Y(Xb))

print(f"\nLagrange Polynomial (a):\n")
sp.pprint(lagrange_a)
print(f"\n\nLagrange Polynomial (b):\n")
sp.pprint(lagrange_b)

newton_a = Newton_polynomial(Xa, get_Y(Xa))
newtone_b = Newton_polynomial(Xb, get_Y(Xb))

print(f"\n\nNewton Polynomial (a):\n")
sp.pprint(newton_a)
sp.pprint(sp.simplify(newton_a))
print(f"\n\nNewton Polynomial (b):\n")
sp.pprint(newtone_b)
sp.pprint(sp.simplify(newtone_b))

print(f"\n\nCheck (a):\n")
lag = lagrange(Xa, get_Y(Xa))
print(lag)
print(f"\nCheck (b):\n")
lag = lagrange(Xb, get_Y(Xb))
print(f"{lag}\n")

f_exact = f(X_star)
P_X_star = lagrange_a.subs(x_sym, X_star)
interpolation_error = abs(f_exact - P_X_star)
print(f"\nExact function value at X_star = {f_exact}")
print(f"Interpolated value at X_star = {P_X_star}")
P_X_star = newton_a.subs(x_sym, X_star)
print(f"Interpolated value at X_star = {P_X_star}")
print(f"Interpolation error = {interpolation_error}")