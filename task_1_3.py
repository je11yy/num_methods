import numpy as np

matrix = np.array([
    [26, -9, -8, 8],
    [9, -21, -2, 8],
    [-3, 2, -18, 8],
    [1, -6, -1, 11]
], dtype=float)

vector = np.array([ 20, -164, 140, -81 ], dtype=float)

epsilon = 1e-6

def norm(matrix):
    max_sum = 0.0
    for line in matrix:
        summary = sum(line)
        if max_sum < summary:
            max_sum = summary
    return max_sum

def solve(A, b, epsilon):
    n = len(A)
    assert A.shape == (n, n), "Матрица A должна быть квадратной"

    B = np.zeros_like(A)
    c = np.zeros_like(b)

    # Формирование матрицы B и вектора c
    for i in range(n):
        B[i] = -A[i] / A[i, i]
        B[i, i] = 0  # Элементы на главной диагонали зануляем
        c[i] = b[i] / A[i, i]
    
    print(f"Норма B: {norm(B)}")
    
    print("Метод простых итераций")
    x, iterations = solve_with_iterations(B, c, epsilon)
    print(f"Решение: {x}\nВыполнено итераций: {iterations}")
    print()

    print("Метод Зейделя")
    x, iterations = solve_with_seidel(B, c, epsilon)
    print(f"Решение: {x}\nВыполнено итераций: {iterations}")
    
def solve_with_iterations(C, F, epsilon):
    # начальное приближение - вектор правых частей
    x = np.copy(F)
    iteration_count = 0
    while True:
        x_prev = np.copy(x)
        # Вычисляем новое приближение
        x = np.dot(C, x_prev) + F
        # Проверяем условие остановки (норма разности решений)
        if np.linalg.norm(x - x_prev, ord=np.inf) < epsilon:
            break
        
        iteration_count += 1

    return x, iteration_count

def solve_with_seidel(C, F, epsilon):
    x = np.copy(F)
    n = len(C)
    # Итерационный процесс
    iteration_count = 0

    while True:
        x_prev = np.copy(x)
        # Вычисляем новое приближение
        for i in range(n):
            part = 0.0
            for j in range(n):
                part += (C[i][j] * x[j])
            part += F[i]
            x[i] = part
        
        # Проверяем условие остановки (норма разности решений)
        if np.linalg.norm(x - x_prev, ord=np.inf) < epsilon:
            break
        
        iteration_count += 1

    return x, iteration_count

solve(matrix, vector, epsilon)