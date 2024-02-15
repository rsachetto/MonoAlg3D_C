def calc_sigma (intra_conductivity, extra_conductivity):
	return (intra_conductivity*extra_conductivity)/(intra_conductivity+extra_conductivity)

# N-benchmark conductivities
sigma_i_l = 0.17		# Intracellular longitudinal {S/m}
sigma_e_l = 0.62		# Extracellular longitudinal {S/m}
sigma_i_t = 0.019		# Intracellular transversal {S/m}
sigma_e_t = 0.24		# Extracellular transversal {S/m}

print("sigma_l = %g" % (calc_sigma(sigma_i_l, sigma_e_l)))
print("sigma_t = %g" % (calc_sigma(sigma_i_t, sigma_e_t)))
