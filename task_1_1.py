coefficients_a = [
    [-1.0, -7.0, -3.0, -2.0],
    [-8.0, 1.0, -9.0, 0.0],
    [8.0, 2.0, -5.0, -3.0],
    [-5.0, 3.0, 5.0, -9.0]
]

coefficients_b = [ -12.0, -60.0, -91.0, -43.0 ]

def swap_rows(matrix_a, matrix_b, identity_matrix, row_index_1, row_index_2):
    matrix_a[row_index_1], matrix_a[row_index_2] = matrix_a[row_index_2], matrix_a[row_index_1]
    identity_matrix[row_index_1], identity_matrix[row_index_2] = identity_matrix[row_index_2], identity_matrix[row_index_1]
    matrix_b[row_index_1], matrix_b[row_index_2] = matrix_b[row_index_2], matrix_b[row_index_1]

def divide_row(matrix_a, matrix_b, identity_matrix, row_index, divider):
    matrix_a[row_index] = [ element / divider for element in matrix_a[row_index] ]
    identity_matrix[row_index] = [ element / divider for element in identity_matrix[row_index] ]
    matrix_b[row_index] /= divider

# сложение строки со строкой, умноженной на weight
def combine_rows(matrix_a, matrix_b, identity_matrix, row_index_1, row_index_2, weight):
    matrix_a[row_index_1] = [ (element_1 + element_2 * weight) for element_1, element_2 in zip(matrix_a[row_index_1], matrix_a[row_index_2]) ]
    identity_matrix[row_index_1] = [ (element_1 + element_2 * weight) for element_1, element_2 in zip(identity_matrix[row_index_1], identity_matrix[row_index_2]) ]
    matrix_b[row_index_1] += matrix_b[row_index_2] * weight

def get_triangular_matrix(matrix_a, matrix_b, identity_matrix):
    length = len(matrix_b)
    for column in range(length):
        # поиск максимального по модулю элемента в столбце
        current_row = None
        for row in range(column, len(matrix_a)):
            if current_row is None or abs(matrix_a[row][column]) > abs(matrix_a[current_row][column]):
                current_row = row
        # нет решения
        if current_row is None:
            raise ValueError()
        
        if current_row != column:
            # переставляем строку с макс элементом выше
            swap_rows(matrix_a, matrix_b, identity_matrix, current_row, column)
        
        # нормализация строки
        divide_row(matrix_a, matrix_b, identity_matrix, column, matrix_a[column][column])

        for row in range(column + 1, len(matrix_a)):
            combine_rows(matrix_a, matrix_b, identity_matrix, row, column, -matrix_a[row][column])

def Gauss(matrix_a, matrix_b):
    matrix_x = [0 for elem in matrix_b]
    for i in range(len(matrix_b) - 1, -1, -1):
        matrix_x[i] = matrix_b[i] - sum(x * a for x, a in zip(matrix_x[(i + 1):], matrix_a[i][(i + 1):]))
    print("\n".join("X{0} = {1:.6f}".format(i + 1, x) for i, x in enumerate(matrix_x)))

def print_matrix(matrix):
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            print(format(matrix[i][j], ".6f"), end=" ")
        print("")

# в треугольной матрице определитель - произведение всех элементов главной диагонали
def find_determinant(matrix):
    res = 1
    for i in range(len(matrix)):
        res *= matrix[i][i]
    return res

def get_identity_matrix(n):
    identity_matrix = [[0.0] * n for _ in range(n)]
    for i in range(n):
        identity_matrix[i][i] = 1.0
    return identity_matrix

# теперь надо привести матрицу к единичному виду, тогда в E останется обратная матрица
# то есть нужно обнулить треугольник сверху
def reverse_matrix(matrix, identity_matrix):
    length = len(matrix)
    for column in range(length -1, 0, -1):
        for row in range(column - 1, -1, -1):
            # из combine_rows
            weight = -matrix[row][column]
            matrix[row] = [ (element_1 + element_2 * weight) for element_1, element_2 in zip(matrix[row], matrix[column]) ]
            identity_matrix[row] = [ (element_1 + element_2 * weight) for element_1, element_2 in zip(identity_matrix[row], identity_matrix[column]) ]

identity_matrix = get_identity_matrix(len(coefficients_a))
get_triangular_matrix(coefficients_a, coefficients_b, identity_matrix)
Gauss(coefficients_a, coefficients_b)
print("determinant: " + str(format(find_determinant(coefficients_a), ".6f")))
print("reverse matrix: ")
reverse_matrix(coefficients_a, identity_matrix)
print_matrix(identity_matrix)