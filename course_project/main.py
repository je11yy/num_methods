import functions
import os
import mixture_properties
import config

output_dir = 'graphs'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

T_full, h2_coeffs_low, h2_coeffs_high = functions.process_substance(output_dir, 'H2', config.H2_data_filename, config.H2_H_0, config.H2_S_0)
T_full, o2_coeffs_low, o2_coeffs_high = functions.process_substance(output_dir, 'O2', config.O2_data_filename, config.O2_H_0, config.O2_S_0)

mixture_properties.generate_mixture_plots(T_full, h2_coeffs_low, h2_coeffs_high, o2_coeffs_low, o2_coeffs_high, output_dir)