import numpy as np
import pandas as pd
import config
import os
import graphs


def Cp_function(T, a, b, c, d, e):
    return a + b * (T ** (-1)) + c * T + d * (T ** 2) + e * (T ** 3)

def H_function(T, coeffs, H_0):
    # Численный метод интегрирования (метод трапеций или простое суммирование)
    temperatures = np.linspace(config.T_0, T, 1000)  # 1000 точек для интегрирования
    cp_values = Cp_function(temperatures, *coeffs)
    integral = np.trapz(cp_values, temperatures)
    H = integral + H_0
    return H / 1000

def S_function(T, coeffs, S_0):
    temperatures = np.linspace(config.T_0, T, 1000)
    cp_values = Cp_function(temperatures, *coeffs)
    integral = np.trapz(cp_values / temperatures, temperatures)
    S = integral + S_0
    return S

def least_squares_fit_cp(T, Cp):
    T = np.array(T)
    Cp = np.array(Cp)
    X = np.column_stack([
        np.ones_like(T),       # a
        T**-1,                 # b / T
        T,                     # c * T
        T**2,                  # d * T^2
        T**3                   # e * T^3
    ])

    # Решение системы линейных уравнений: (XᵀX)a = Xᵀy
    XtX = np.dot(X.T, X)
    Xty = np.dot(X.T, Cp)
    coefficients = np.linalg.solve(XtX, Xty)
    return coefficients

def get_real_values(df):
    temperatures = df['Temperature'].values
    cp_values = df['Cp'].values
    S_values = df['S'].values
    H_values = df['H'].values
    return temperatures, cp_values, S_values, H_values

def get_low_values(values, temperatures):
    return values[temperatures <= config.max_low_temperature]

def get_high_values(values, temperatures):
    return values[(temperatures <= config.max_high_temperature) & (temperatures >= config.max_low_temperature)]

def correct_first_coefficient(Cp_value, T, coeffs):
    return Cp_value - (
        coeffs[1] / T +
        coeffs[2] * T +
        coeffs[3] * T**2 +
        coeffs[4] * T**3
    )

def print_data(name, Cp_low, Cp_high):
    print(f"Данные для " + name)
    print(f"Cp при T={config.max_low_temperature} K (низкий диапазон): {Cp_low:.3f}")
    print(f"Cp при T={config.max_low_temperature} K (высокий диапазон, после стыковки): {Cp_high:.3f}")

def get_properties(T, coeffs_low, coeffs_high, H_0, S_0):
    Cp_fit_low = Cp_function(T, *coeffs_low)
    Cp_fit_high = Cp_function(T, *coeffs_high)
    H_fit_low = H_function(T, coeffs_low, H_0)
    H_fit_high = H_function(T, coeffs_high, H_0)
    S_fit_low = S_function(T, coeffs_low, S_0)
    S_fit_high = S_function(T, coeffs_high, S_0)
    return Cp_fit_low, Cp_fit_high, H_fit_low, H_fit_high, S_fit_low, S_fit_high

def get_data(filename, name):
    df = pd.read_csv(filename)
    df['Temperature'] = pd.to_numeric(df['Temperature'], errors='coerce')
    df['Cp'] = pd.to_numeric(df['Cp'], errors='coerce')
    df['S'] = pd.to_numeric(df['S'], errors='coerce')
    df['H'] = pd.to_numeric(df['H'], errors='coerce')
    df = df.dropna()

    temperatures, cp_values, S_values, H_values = get_real_values(df)

    T_low = get_low_values(temperatures, temperatures)
    Cp_low = get_low_values(cp_values, temperatures)

    T_high = get_high_values(temperatures, temperatures)
    Cp_high = get_high_values(cp_values, temperatures)

    coeff_low = least_squares_fit_cp(T_low, Cp_low)
    coeff_high = least_squares_fit_cp(T_high, Cp_high)

    T_check = config.max_low_temperature
    Cp_dock_low = Cp_function(T_check, *coeff_low)
    
    # Корректируем свободный член второго уравнения для стыковки
    a_high_corrected = correct_first_coefficient(Cp_dock_low, T_check, coeff_high)
    coeff_high[0] = a_high_corrected 
    Cp_dock_high = Cp_function(T_check, *coeff_high)
    
    print_data(name, Cp_dock_low, Cp_dock_high)

    return temperatures, cp_values, S_values, H_values, coeff_low, coeff_high

def process_substance(output_dir, name, filename, H_0, S_0):
    (temperatures, cp_values, S_values, H_values, coeff_low, coeff_high) = get_data(filename, name)
    T_full = np.linspace(min(temperatures), max(temperatures), config.max_low_temperature)

    (Cp_fit_low, Cp_fit_high, H_fit_low, H_fit_high, S_fit_low, 
    S_fit_high) = get_properties(T_full, coeff_low, coeff_high, H_0, S_0)

    if not os.path.exists(output_dir + '/' + name):
        os.makedirs(output_dir + '/' + name)

    graphs.generate_plot_with_original(output_dir + '/' + name, 'Cp', temperatures, cp_values, T_full, Cp_fit_low, Cp_fit_high)
    graphs.generate_plot_with_original(output_dir + '/' + name, 'S', temperatures, S_values, T_full, S_fit_low, S_fit_high)
    graphs.generate_plot_with_original(output_dir + '/' + name, 'H', temperatures, H_values, T_full, H_fit_low, H_fit_high)

    return T_full, coeff_low, coeff_high
