import math

coefficients = [
    [0, -14, -6, -78],
    [-9, 15, -1, -73],
    [1, -11, 1, -38],
    [-7, 12, 3, 77],
    [6, -7, 0, 91]
]

# прямой ход
def first_half(coefficients, size):
    p_array = [0 for i in range(0, size)]
    q_array = [0 for i in range(0, size)]

    p_array[0] = -(coefficients[0][2] / coefficients[0][1])
    q_array[0] = coefficients[0][3] / coefficients[0][1]

    for i in range(1, size - 1):
        # проверка устойчивости
        a = coefficients[i][0]
        b = coefficients[i][1]
        c = coefficients[i][2]

        if (a != 0) and (c != 0) and (abs(b) < a + c):
            raise ValueError()

        # отдельно выносим вычисление знаменателя, чтобы не считать его два раза
        denominator = coefficients[i][1] + coefficients[i][0] * p_array[i - 1]
        p_array[i] = (-coefficients[i][2]) / denominator

        # проверка, что Pi <= 1
        if abs(p_array[i]) > 1:
            raise ValueError()

        q_array[i] = (coefficients[i][3] - coefficients[i][0] * q_array[i - 1]) / denominator

    p_array[size - 1] = 0
    q_array[size - 1] = (coefficients[size - 1][3] - coefficients[size - 1][0] * q_array[size - 2]) / (coefficients[size - 1][1] + coefficients[size - 1][0] * p_array[size - 2])

    return [p_array, q_array]

# обратный ход
def second_half(p_array, q_array, size):
    x_array = [0 for i in range(0, size)]
    x_array[size - 1] = q_array[size - 1]

    for i in range(size - 2, 0, -1):
        x_array[i] = p_array[i] * x_array[i + 1] + q_array[i]
    
    x_array[0] = p_array[0] * x_array[1] + q_array[0]

    return x_array

# полное решение
def solution(coefficients, size):
    p_array, q_array = first_half(coefficients, size)
    return [second_half(p_array, q_array, size), find_determinant(coefficients, p_array, size)]

# печать решения
def print_solution(x_array):
    counter = 1
    for element in x_array:
        print("x" + str(counter) + ": " + format(element, '.6f'))
        counter += 1

# поиск определителя
def find_determinant(coefficients, p_array, size):
    result = 1
    for i in range(0, size):
        result *= coefficients[i][1] + coefficients[i][0] * p_array[i - 1]
    return result

try:
    x_array, determinant = solution(coefficients, 5)
except ValueError:
    print("Stability check failed")
else:
    print_solution(x_array)
    print("determinant: " + str(determinant))