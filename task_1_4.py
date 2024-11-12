from math import sqrt
import numpy as np

matrix = np.array([
    [ 8, 2, -1 ],
    [ 2, -5, -8 ],
    [ -1, -8, -5 ],
], dtype=float)

epsilon = 1e-10

def rotation_method(matrix, epsilon):
    if not np.array_equal(matrix, matrix.T):
        raise ValueError("Метод вращения Якоби: матрица должна быть симметрической")
    
    n = matrix.shape[0]

    # собственные значение
    eigen_values = np.copy(matrix)
    # собственные векторы, единичная матрица
    eigen_vectors = np.eye(n)

    # поиск максимального элемента по модулю вне диагонали
    def max_off_diagonal(matrix):
        max_val = 0
        p, q = 0, 1
        for i in range(n):
            for j in range(i + 1, n):
                if abs(matrix[i, j]) > abs(max_val):
                    max_val = matrix[i, j]
                    p, q = i, j
        return p, q
    
    def off_diagonal_sum_squares(matrix):
        return sqrt(np.sum(matrix**2) - np.sum(np.diag(matrix)**2))

    iterations = 0

    while off_diagonal_sum_squares(eigen_values) >= epsilon:
        p, q = max_off_diagonal(eigen_values)
        
        if eigen_values[p, p] == eigen_values[q, q]:
            angle = np.pi / 4
        else:
            angle = 0.5 * np.arctan(2 * eigen_values[p, q] / (eigen_values[p, p] - eigen_values[q, q]))
        
        sin, cos = np.sin(angle), np.cos(angle)
        
        rotation = np.eye(n)
        rotation[p, p] = cos
        rotation[q, q] = cos
        rotation[p, q] = -sin
        rotation[q, p] = sin

        # @ - матричное умножение (обычное * просто перемножает элементы матриц)
        eigen_values = rotation.T @ eigen_values @ rotation
        eigen_vectors = eigen_vectors @ rotation

        iterations += 1

    eigenvalues_diag = np.diag(eigen_values)
    sorted_indices = np.argsort(eigenvalues_diag)[::-1]  # Сортировка по убыванию
    sorted_eigenvalues = eigenvalues_diag[sorted_indices]
    sorted_eigenvectors = eigen_vectors[:, sorted_indices]

    return sorted_eigenvalues, sorted_eigenvectors, iterations

def power_method(matrix, epsilon):
    y = np.ones_like(matrix[0])
    eigen_value = y[0]

    # Умножаем матрицу на текущий вектор
    iterations = 0

    while True:
        z = matrix @ y
        eigen_value_old = eigen_value
        eigen_value = z[0] / y[0]
        y = z / np.linalg.norm(z)

        iterations += 1

        # Проверка сходимости
        if abs(eigen_value_old - eigen_value) < epsilon:
            break

    return eigen_value, y, iterations

np.set_printoptions(precision=8, suppress=True, formatter={'all': lambda x: f'{x:0.3f}'})

print("Метод вращений якоби\n")
eigen_values, eigen_vectors, iterations = rotation_method(matrix, epsilon)
print(f"Собственные значения:\n{eigen_values}\nСобственные векторы:\n{eigen_vectors}\n")
print(f"Количество итераций: {iterations}\n")

print("Степенной метод\n")
max_eigen_value, y, iterations = power_method(matrix, epsilon)
print(f"Максимальное собственное значение по модулю (спектральный радиус):\n", format(max_eigen_value, "0.3f"))
print(f"Собственный вектор:\n{y}")
print(f"Количество итераций: {iterations}\n")

print("\nПРОВЕРКА\n")
real_eigen_values, real_eigen_vectors = np.linalg.eig(matrix)
print(f"Собственные значения:\n{real_eigen_values}\nСобственные векторы:\n{real_eigen_vectors}")