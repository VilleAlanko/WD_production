import numpy as np
from scipy.stats import chi2

PDF_sets = ['CT18ANLO', 'MSHT20nlo_as118', 'NNPDF40_nlo_pch_as_01180']
which_cross_sections_included = 'D'

dof = 10

print()
print('Which cross sections included:', which_cross_sections_included)
print()

for PDF_set in PDF_sets:
    if (PDF_set != 'NNPDF40_nlo_pch_as_01180'):
        y_best_eta_lept = np.loadtxt('input/theory_values/HESSIAN/best/eta_lept_' + which_cross_sections_included + '_' + PDF_set + '_best.txt')
        y_best_pTD = np.loadtxt('input/theory_values/HESSIAN/best/pTD_' + which_cross_sections_included + '_' + PDF_set + '_best.txt')
    else:
        y_best_eta_lept = np.loadtxt('input/theory_values/MC/best/eta_lept_' + which_cross_sections_included + '_' + PDF_set + '_best.txt')
        y_best_pTD = np.loadtxt('input/theory_values/MC/best/pTD_' + which_cross_sections_included + '_' + PDF_set + '_best.txt')

    y_exp_eta_lept = np.loadtxt('input/experimental_values/eta_lept_' + which_cross_sections_included + '.txt')
    y_exp_pTD = np.loadtxt('input/experimental_values/pTD_' + which_cross_sections_included + '.txt')

    N = len(y_best_eta_lept)

    C_inverse_eta_lept = np.loadtxt('input/covariance_matrix/eta_lept_' + which_cross_sections_included + '.txt', delimiter=' ')
    C_inverse_pTD = np.loadtxt('input/covariance_matrix/pTD_' + which_cross_sections_included + '.txt', delimiter=' ')

    chi_squared = 0.

    for i in range(N):
        for j in range(N):
            chi_squared += (y_exp_eta_lept[i] - y_best_eta_lept[i]) * C_inverse_eta_lept[i, j] * (y_exp_eta_lept[j] - y_best_eta_lept[j])
            chi_squared += (y_exp_pTD[i] - y_best_pTD[i]) * C_inverse_pTD[i, j] * (y_exp_pTD[j] - y_best_pTD[j])

    print(PDF_set)
    print('Chi squared:', chi_squared)

    p_val = chi2.sf(chi_squared, dof)

    print('p value:', p_val)
    print()
