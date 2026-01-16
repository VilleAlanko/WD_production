import numpy as np

PDF_set = 'NNPDF40_nlo_pch_as_01180'
which_cross_sections_included = 'D'

num_members = 100
N_data = 10

y_eta_lept = np.loadtxt('input/theory_values/MC/variation/eta_lept_' + which_cross_sections_included + '_' + PDF_set + '_vals.txt', delimiter=',')
y_best_eta_lept = np.loadtxt('input/theory_values/MC/best/eta_lept_' + which_cross_sections_included + '_' + PDF_set + '_best.txt')
y_exp_eta_lept = np.loadtxt('input/experimental_values/eta_lept_' + which_cross_sections_included + '.txt')
C_inverse_eta_lept = np.loadtxt('input/covariance_matrix/eta_lept_' + which_cross_sections_included + '.txt', delimiter=' ')

y_pTD = np.loadtxt('input/theory_values/MC/variation/pTD_' + which_cross_sections_included + '_' + PDF_set + '_vals.txt', delimiter=',')
y_best_pTD = np.loadtxt('input/theory_values/MC/best/pTD_' + which_cross_sections_included + '_' + PDF_set + '_best.txt')
y_exp_pTD = np.loadtxt('input/experimental_values/pTD_' + which_cross_sections_included + '.txt')
C_inverse_pTD = np.loadtxt('input/covariance_matrix/pTD_' + which_cross_sections_included + '.txt', delimiter=' ')

chi_k_squared = np.zeros(num_members)

for k in range(num_members):
    for i in range(int(N_data / 2) - 0):
        for j in range(int(N_data / 2) - 0):
            chi_k_squared[k] += 1.645**2 * (y_eta_lept[k, i] - y_exp_eta_lept[i]) * C_inverse_eta_lept[i, j] * (y_eta_lept[k, j] - y_exp_eta_lept[j])
            #print((y_eta_lept[k, i] - y_exp_eta_lept[i]) * C_inverse_eta_lept[i, j] * (y_eta_lept[k, j] - y_exp_eta_lept[j]), k, i, j)
    for i in range(int(N_data / 2) - 0):
        for j in range(int(N_data / 2) - 0):
            chi_k_squared[k] += 1.645**2 * (y_pTD[k, i] - y_exp_pTD[i]) * C_inverse_pTD[i, j] * (y_pTD[k, j] - y_exp_pTD[j])

omega_k_chi_squared = np.zeros(num_members)

sum_in_omega = 0.

for i in range(num_members):
    sum_in_omega += chi_k_squared[i]**((N_data * 1. - 1. - 0) / 2.) * np.exp(-chi_k_squared[i] / 2.)

for k in range(num_members):
    omega_k_chi_squared[k] = chi_k_squared[k]**((N_data * 1. - 1. - 0) / 2.) * np.exp(-chi_k_squared[k] / 2.) / (1. / (num_members * 1.) * sum_in_omega)

    if omega_k_chi_squared[k] > 1:
        print(omega_k_chi_squared[k], k)

np.savetxt('output/omega.txt', omega_k_chi_squared)
    
