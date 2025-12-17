import numpy as np
from numpy.linalg import inv

PDF_set = 'CT18ANLO'
which_cross_sections_included = 'both'

y_plus_eta_lept = np.loadtxt('input/theory_values/HESSIAN/variation/eta_lept_' + which_cross_sections_included + '_' + PDF_set + '_plus.txt', delimiter=',')
y_minus_eta_lept = np.loadtxt('input/theory_values/HESSIAN/variation/eta_lept_' + which_cross_sections_included + '_' + PDF_set + '_minus.txt', delimiter=',')
y_plus_pTD = np.loadtxt('input/theory_values/HESSIAN/variation/pTD_' + which_cross_sections_included + '_' + PDF_set + '_plus.txt', delimiter=',')
y_minus_pTD = np.loadtxt('input/theory_values/HESSIAN/variation/pTD_' + which_cross_sections_included + '_' + PDF_set + '_minus.txt', delimiter=',')

y_best_eta_lept = np.loadtxt('input/theory_values/HESSIAN/best/eta_lept_' + which_cross_sections_included + '_' + PDF_set + '_best.txt')
y_best_pTD = np.loadtxt('input/theory_values/HESSIAN/best/pTD_' + which_cross_sections_included + '_' + PDF_set + '_best.txt')

y_exp_eta_lept = np.loadtxt('input/experimental_values/eta_lept_' + which_cross_sections_included + '.txt')
y_exp_pTD = np.loadtxt('input/experimental_values/pTD_' + which_cross_sections_included + '.txt')

D_eta_lept = (y_plus_eta_lept - y_minus_eta_lept) / 2.
D_pTD = (y_plus_pTD - y_minus_pTD) / 2.

if (PDF_set == 'MSHT20nlo_as118'):
    D_eta_lept = D_eta_lept * 1.645
    D_pTD = D_pTD * 1.645

C_inverse_eta_lept = np.loadtxt('input/covariance_matrix/eta_lept.txt', delimiter=' ')
C_inverse_pTD = np.loadtxt('input/covariance_matrix/pTD.txt', delimiter=' ')

t = np.sqrt(10.)

# Number of members
N = len(D_eta_lept[0, :])
print(N)
# Number of data points / 2
M = len(D_eta_lept[:, 0])

B = np.zeros((N, N))

for k in range(N):
    for n in range(N):
        if (k == n):
            B[k, n] += t**2
        for i in range(M):
            for j in range(M):
                B[k, n] += D_eta_lept[i, k] * C_inverse_eta_lept[i, j] * D_eta_lept[j, n]
                B[k, n] += D_pTD[i, k] * C_inverse_pTD[i, j] * D_pTD[j, n]

a = np.zeros(N)

for k in range(N):
    for i in range(M):
        for j in range(M):
            a[k] += D_eta_lept[i, k] * C_inverse_eta_lept[i, j] * (y_best_eta_lept[j] - y_exp_eta_lept[j])
            a[k] += D_pTD[i, k] * C_inverse_pTD[i, j] * (y_best_pTD[j] - y_exp_pTD[j])

wmin = -1. * inv(B) @ a

P = sum(wmin**2)
print('P / delta chi^2 =', P)

np.savetxt('output/wmin_' + PDF_set + '_' + which_cross_sections_included + '.txt', wmin)

eps, v = np.linalg.eig(B)

dw = np.zeros((N, N))

for i in range(N):
    for k in range(N):
        dw[i, k] = v[k][i] * np.sqrt(1. / eps[k]) * t

np.savetxt('output/dw_' + PDF_set + '_' + which_cross_sections_included + '.txt', dw)

delta_y = np.zeros(2 * M)

for i in range(M):
    delta_y[i] = y_best_eta_lept[i] - y_exp_eta_lept[i]
    delta_y[M + i] = y_best_pTD[i] - y_exp_pTD[i]






