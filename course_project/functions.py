import numpy as np
import pandas as pd
import config
import os
import graphs
from scipy.linalg import lstsq

def Cp_function_low(T, a, b, c, d, e):
    return a + b * (T ** (-2)) + c * T + d * (T ** 2) + e * (T ** 3)

def Cp_function_high(T, a, b, c, d, e):
    T_check = config.max_low_temperature
    return a + b * ((1 / (T**2)) - (1 / (T_check ** 2))) + c * (T - T_check) + d * ((T**2) - (T_check**2)) + e * ((T**3) - (T_check**3))

def H_function_low(T, a, b, c, d , e, H_0):
    x = T
    integral_1 = a * x - b / x + c * (x**2 / 2) + d * (x**3 / 3) + e * (x**4 / 4)
    x = config.T_0
    integral_2 = a * x - b / x + c * (x**2 / 2) + d * (x**3 / 3) + e * (x**4 / 4)
    H = integral_1 - integral_2 + H_0
    return H / 1000

def H_function_high(T, a, b, c, d, e, H_0):
    x = T
    T_check = config.max_low_temperature
    integral_1 = a * x + b * ((-1 / x) - (x / (T_check**2))) + c * (((x**2) / 2) - (x * T_check)) + d * (((x**3) / 3) - (T_check**2) * x) + e * (((x**4) / 4) - ((T_check**3) * x))
    x = config.T_0
    integral_2 = a * x + b * ((-1 / x) - (x / (T_check**2))) + c * (((x**2) / 2) - (x * T_check)) + d * (((x**3) / 3) - (T_check**2) * x) + e * (((x**4) / 4) - ((T_check**3) * x))
    H = integral_1 - integral_2 + H_0
    return H / 1000

def S_function_low(T, a, b, c, d , e, S_0):
    x = T
    integral_1 = a * np.log(np.abs(x)) - b / (2 * (x)) + c * (x) + d * ((x**2) / 2) + e * ((x**3) / 3)
    x = config.T_0
    integral_2 = a * np.log(np.abs(x)) - b / (2 * (x)) + c * (x) + d * ((x**2) / 2) + e * ((x**3) / 3)
    S = integral_1 - integral_2 + S_0
    return S

def S_function_high(T, a, b, c, d, e, S_0):
    x = T
    T_check = config.max_low_temperature
    integral_1 = a * np.log(np.abs(x)) + b * (-1 / (2 * (x**2)) - np.log(np.abs(x)) / (T_check**2)) + c * (x - np.log(np.abs(x)) * T_check) + d * ((x**2 / 2) - (T_check**2) * np.log(np.abs(x))) + e * (x**3 / 3 - (T_check**3) * np.log(np.abs(x)))
    x = config.T_0
    integral_2 = a * np.log(np.abs(x)) + b * (-1 / (2 * (x**2)) - np.log(np.abs(x)) / (T_check**2)) + c * (x - np.log(np.abs(x)) * T_check) + d * ((x**2 / 2) - (T_check**2) * np.log(np.abs(x))) + e * (x**3 / 3 - (T_check**3) * np.log(np.abs(x)))
    S = integral_1 - integral_2 + S_0
    return S

def least_squares_fit_cp_2(T, Cp):
    T = np.array(T)
    Cp = np.array(Cp)
    X = np.column_stack([
        T**-2 - config.max_low_temperature**-2,             # b * (1 / T^2 - 1 / (T_docking)^2)
        T - config.max_low_temperature,                     # c * (T- T_docking)
        T**2 - config.max_low_temperature**2,               # d * (T^2- T_docking ^ 2)
        T**3 - config.max_low_temperature**3                # e * (T^3- T_docking ^ 3)
    ])

    # Решение системы линейных уравнений: (XᵀX)a = Xᵀy
    XtX = np.dot(X.T, X)
    Xty = np.dot(X.T, Cp)
    coefficients = np.linalg.solve(XtX, Xty)
    return coefficients


def least_squares_fit_cp(T, Cp):
    T = np.array(T)
    Cp = np.array(Cp)
    X = np.vstack([
        np.ones(len(T)),       # a
        T**-2,                 # b / T
        T,                     # c * T
        T**2,                  # d * T^2
        T**3                   # e * T^3
    ]).T

    # Решение системы линейных уравнений: (XᵀX)a = Xᵀy
    XtX = np.dot(X.T, X)
    Xty = np.dot(X.T, Cp)
    coefficients = np.linalg.solve(XtX, Xty)
    return coefficients

def get_real_values(df):
    temperatures = df['Temperature'].values
    cp_values = df['Cp'].values
    # cp_values = cp_values[temperatures >= 100]
    S_values = df['S'].values
    # S_values = S_values[temperatures >= 100]
    H_values = df['H'].values
    # H_values = H_values[temperatures >= 100]
    # temperatures = temperatures[temperatures >= 100]
    return temperatures, cp_values, S_values, H_values

def get_low_values(values, temperatures):
    return values[(temperatures <= config.max_low_temperature)]

def get_high_values(values, temperatures):
    return values[(temperatures <= config.max_high_temperature) & (temperatures >= config.max_low_temperature)]

def print_data(name, Cp_low, Cp_high):
    print(f"Данные для " + name)
    print(f"Cp при T={config.max_low_temperature} K (низкий диапазон): {Cp_low:.3f}")
    print(f"Cp при T={config.max_low_temperature} K (высокий диапазон, после стыковки): {Cp_high:.3f}")

def get_properties(T, coeffs_low, coeffs_high, H_0, S_0):
    T_low = T[(T <= config.max_low_temperature)]
    T_high = T[(T >= config.max_low_temperature) & (T <= config.max_high_temperature)]
    Cp_fit_low = Cp_function_low(T_low, *coeffs_low)
    Cp_fit_high = Cp_function_high(T_high, *coeffs_high)
    H_fit_low = H_function_low(T_low, *coeffs_low, H_0)
    H_fit_high = H_function_high(T_high, *coeffs_high, H_0)
    S_fit_low = S_function_low(T_low, *coeffs_low, S_0)
    S_fit_high = S_function_high(T_high, *coeffs_high, S_0)
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

    T_check = config.max_low_temperature

    a_high = Cp_function_low(T_check, *coeff_low)
    coeff_high_tmp = least_squares_fit_cp_2(T_high, Cp_high - a_high)
    coeff_high = np.insert(coeff_high_tmp, 0, a_high)

    Cp_dock_low = Cp_function_low(T_check, *coeff_low)
    Cp_dock_high = Cp_function_high(T_check, *coeff_high)
    
    print_data(name, Cp_dock_low, Cp_dock_high)
    print(coeff_low)
    print(coeff_high)

    return temperatures, cp_values, S_values, H_values, coeff_low, coeff_high

def process_substance(output_dir, name, filename, H_0, S_0):
    (temperatures, cp_values, S_values, H_values, coeff_low, coeff_high) = get_data(filename, name)
    T_full = np.linspace(100, config.max_high_temperature, 100)

    (Cp_fit_low, Cp_fit_high, H_fit_low, H_fit_high, S_fit_low, 
    S_fit_high) = get_properties(T_full, coeff_low, coeff_high, H_0, S_0)

    print(((np.sum(Cp_fit_low) - np.sum(cp_values[temperatures < config.max_low_temperature]))**2) / np.size(Cp_fit_low))

    if not os.path.exists(output_dir + '/' + name):
        os.makedirs(output_dir + '/' + name)

    graphs.generate_plot_with_original(output_dir + '/' + name, 'Cp', temperatures, cp_values, T_full, Cp_fit_low, Cp_fit_high)
    graphs.generate_plot_with_original(output_dir + '/' + name, 'S', temperatures, S_values, T_full, S_fit_low, S_fit_high)
    graphs.generate_plot_with_original(output_dir + '/' + name, 'H', temperatures, H_values, T_full, H_fit_low, H_fit_high)

    return T_full, coeff_low, coeff_high
