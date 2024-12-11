import functions
import os
import graphs
import config

# Рассчитаем для смеси H2 и O2
def calculate_mixture_properties(T_full, coeff_h2_low, coeff_h2_high, coeff_o2_low, coeff_o2_high, ratio_h2, ratio_o2):
    # Расчет Cp для смеси
    Cp_mixture_low = (ratio_h2 * functions.Cp_function(T_full, *coeff_h2_low) + ratio_o2 * functions.Cp_function(T_full, *coeff_o2_low)) / (ratio_h2 + ratio_o2)
    Cp_mixture_high = (ratio_h2 * functions.Cp_function(T_full, *coeff_h2_high) + ratio_o2 * functions.Cp_function(T_full, *coeff_o2_high)) / (ratio_h2 + ratio_o2)
    
    # Расчет H для смеси
    H_mixture_low = (ratio_h2 * functions.H_function(T_full, coeff_h2_low, config.H2_H_0) + ratio_o2 * functions.H_function(T_full, coeff_o2_low, config.O2_H_0)) / (ratio_h2 + ratio_o2)
    H_mixture_high = (ratio_h2 * functions.H_function(T_full, coeff_h2_high, config.H2_H_0) + ratio_o2 * functions.H_function(T_full, coeff_o2_high, config.O2_H_0)) / (ratio_h2 + ratio_o2)
    
    # Расчет S для смеси
    S_mixture_low = (ratio_h2 * functions.S_function(T_full, coeff_h2_low, config.H2_S_0) + ratio_o2 * functions.S_function(T_full, coeff_o2_low, config.O2_S_0)) / (ratio_h2 + ratio_o2)
    S_mixture_high = (ratio_h2 * functions.S_function(T_full, coeff_h2_high, config.H2_S_0) + ratio_o2 * functions.S_function(T_full, coeff_o2_high, config.O2_S_0)) / (ratio_h2 + ratio_o2)
    
    return Cp_mixture_low, Cp_mixture_high, H_mixture_low, H_mixture_high, S_mixture_low, S_mixture_high

# генерация графиков для смесей
def generate_mixture_plots(T_full, coeff_h2_low, coeff_h2_high, coeff_o2_low, coeff_o2_high, output_dir):
    output_dir += '/mixtures'

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    ratios = [(1, 1), (1, 2), (2, 1), (1, 3), (3, 1)]
    for ratio in ratios:
        ratio_h2, ratio_o2 = ratio
        Cp_mixture_low, Cp_mixture_high, H_mixture_low, H_mixture_high, S_mixture_low, S_mixture_high = calculate_mixture_properties(T_full, coeff_h2_low, coeff_h2_high, coeff_o2_low, coeff_o2_high, ratio_h2, ratio_o2)
        
        graphs.generate_plot_for_mixture(output_dir, 'Cp', T_full, ratio_h2, ratio_o2, Cp_mixture_low, Cp_mixture_high)
        graphs.generate_plot_for_mixture(output_dir, 'H', T_full, ratio_h2, ratio_o2, H_mixture_low, H_mixture_high)
        graphs.generate_plot_for_mixture(output_dir, 'S', T_full, ratio_h2, ratio_o2, S_mixture_low, S_mixture_high)


    