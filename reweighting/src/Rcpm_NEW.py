import numpy as np
import matplotlib.pyplot as plt


def eta_lept(which_cross_sections_included, PDF_set, num_members):
    if (PDF_set == 'NNPDF40_nlo_pch_as_01180'):
        old_member_vals = np.loadtxt('input/theory_values/MC/variation/eta_lept_' + \
                            which_cross_sections_included + '_' + PDF_set + '_vals.txt', delimiter=',')
    else:
        old_member_vals = np.loadtxt('input/theory_values/HESSIAN/variation/eta_lept_' + \
                            which_cross_sections_included + '_' + PDF_set + '_vals.txt', delimiter=',')
    
    omega = np.loadtxt('output/omega.txt')
    
    Rcpm_eta_lept_OLD = np.zeros(5)
    Rcpm_eta_lept_NEW = np.zeros(5)

    for eta_lept_index in range(5):
        for member in range(num_members):
            Rcpm_eta_lept_NEW[eta_lept_index] += omega[member] * old_member_vals[member, eta_lept_index]
            Rcpm_eta_lept_OLD[eta_lept_index] += old_member_vals[member, eta_lept_index]
    
    Rcpm_eta_lept_OLD = Rcpm_eta_lept_OLD / (num_members * 1.)
    Rcpm_eta_lept_NEW = Rcpm_eta_lept_NEW / (num_members * 1.)

    for member in range(num_members):
        if (omega[member] > 1):
            plt.scatter([0, 1, 2, 3, 4], old_member_vals[member, :], color='green', zorder=2)
        else:
            plt.scatter([0, 1, 2, 3, 4], old_member_vals[member, :], color='black', zorder=1)
    
    plt.scatter([0, 1, 2, 3, 4], np.loadtxt('input/experimental_values/eta_lept_both.txt'), color='red', zorder=2)
    plt.show()

    np.savetxt('output/new_Rcpm_eta_lept/' + PDF_set + '/' + which_cross_sections_included + '/best.txt', Rcpm_eta_lept_NEW)

which_cross_sections_included = 'both'
PDF_set = 'NNPDF40_nlo_pch_as_01180'
#PDF_set = 'MSHT20nlo_as118'
num_members = 100

eta_lept(which_cross_sections_included, PDF_set, num_members)
    
