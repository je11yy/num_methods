import matplotlib.pyplot as plt
import os
import config

def generate_plot_with_original(output_dir, name, temperatures, values, T_full, fit_low, fit_high):
    mask_low = (T_full <= config.max_low_temperature)
    T_low = T_full[mask_low]

    # Фильтруем T_full и fit_high для высокого диапазона
    mask_high = (T_full >= config.max_low_temperature) & (T_full <= config.max_high_temperature)
    T_high = T_full[mask_high]

    values = values[(temperatures < config.max_high_temperature)]
    temperatures = temperatures[temperatures < config.max_high_temperature]
    
    plt.figure(figsize=(10, 6))
    plt.scatter(temperatures, values, color='black', label='Реальные данные', alpha=0.5, s=10)
    plt.plot(T_low, fit_low, color='blue', label='Низкий диапазон температур')
    plt.plot(T_high, fit_high, color='red', label='Высокий диапазон температур')
    plt.xlabel('Температура (K)')
    plt.ylabel(name)
    plt.title('Зависимость ' + name + ' от температуры')
    plt.legend()
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, name + '.png'))
    plt.close() 

def generate_plot_for_mixture(output_dir, name, T_full, ratio_h2, ratio_o2, fit_low, fit_high):
    plt.figure(figsize=(10, 6))
    plt.plot(T_full[(T_full > 100) & (T_full < config.max_low_temperature)], fit_low[(T_full < config.max_low_temperature)], color='blue', label='Низкий диапазон температур')
    plt.plot(T_full[(T_full > 100) & (T_full > config.max_low_temperature) & (T_full < config.max_high_temperature)], fit_high[(T_full > config.max_low_temperature) & (T_full < config.max_high_temperature)], color='red', label='Высокий диапазон температур')
    plt.xlabel('Температура (K)')
    plt.ylabel(name)
    plt.title('Зависимость ' + name + ' смеси H2 и O2 от температуры')
    plt.legend()
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, name + '(' + str(int(ratio_h2)) + '-' + str(int(ratio_o2)) + ').png'))
    plt.close()