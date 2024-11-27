import numpy as np
x = [0.0, 0.1, 0.2, 0.3, 0.4]
y = [1.0, 1.1052, 1.2214, 1.3499, 1.4918]
x_star = 0.2

def find_i(x, x_star):
    for i in range(len(x)):
        if (x[i] == x_star):
            return i
    raise ValueError("Can't calculate")


# первый порядок точности
def find_left_diff(x, y, index):
    return (y[index] - y[index - 1]) / (x[index] - x[index - 1])

# первый порядок точности
def find_right_diff(x, y, index):
    return (y[index + 1] - y[index]) / (x[index + 1] - x[index])

# второй порядок точности (первая производная)
# если сетка равномерная, то можно найти, как полусумму
def find_first_diff(x, y, x_star, index):
    left_diff = find_left_diff(x, y, index)
    right_diff = find_right_diff(x, y, index)
    part = 2 * x_star - x[index - 1] - x[index]
    return left_diff + ((right_diff - left_diff) / (x[index + 1] - x[index - 1])) * part

def find_second_diff(x, y, index):
    left_diff = find_left_diff(x, y, index)
    right_diff = find_right_diff(x, y, index)
    return 2.0 * ((right_diff - left_diff) / (x[index + 1] - x[index - 1]))

index = find_i(x, x_star)
first_diff = find_first_diff(x, y, x_star, index)
second_diff = find_second_diff(x, y, index)

left_diff = find_left_diff(x, y, index)
right_diff = find_right_diff(x, y, index)

print(f"Первая производная: {first_diff:.6f}")
print(f"Проверка, что полусумма: {((left_diff + right_diff) / 2.0):.6f}")
print(f"Вторая производная: {second_diff:.6f}")

print(f"Левосторонняя производная: {left_diff:.6f}")
print(f"Правосторонняя производная: {right_diff:.6f}")

coefficients = np.polyfit(x, y, 2)  
polynomial = np.poly1d(coefficients)
polynomial_derivative = polynomial.deriv()
print(f"Проверка первой производной: {polynomial_derivative(x_star):.6f}")
polynomial_derivative = polynomial_derivative.deriv()
print(f"Проверка второй производной: {polynomial_derivative(x_star):.6f}")