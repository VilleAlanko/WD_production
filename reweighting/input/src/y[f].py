import numpy as np

PDF_set = 'MSHT20nlo_as118'
n_members = 64
values_main_path = '../../WD_production/output/'

z_def = 'minus'
FF_set = 'opal'

#--------------------------------------------------------------------------------------------------------------------------------------#
#                                                        Rcpm as a function of eta_lept                                                #
#--------------------------------------------------------------------------------------------------------------------------------------#

Rcpm_plus_member = np.zeros(int(n_members / 2))
Rcpm_minus_member = np.zeros(int(n_members / 2))

for member_index in range(n_members):
    for eta_lept_index in range(5):
        WmDp_sigma = sum(sum(np.loadtxt(values_main_path + 'W-D+/NLO/' + z_def + '/' + FF_set + '/pdf_errs/' + PDF_set + \
                    '/' + str(member_index) + '_' + str(eta_lept_index) + '_vals.txt', delimiter=','))) - \
                    sum(sum(np.loadtxt(values_main_path + 'W-D+/subtraction/' + z_def + '/' + FF_set + '/pdf_errs/' + PDF_set + \
                    '/' + str(member_index) + '_' + str(eta_lept_index) + '_vals.txt', delimiter=',')))

        WpDm_sigma = sum(sum(np.loadtxt(values_main_path + 'W+D-/NLO/' + z_def + '/' + FF_set + '/pdf_errs/' + PDF_set + \
                    '/' + str(member_index) + '_' + str(eta_lept_index) + '_vals.txt'))) - \
                    sum(sum(np.loadtxt(values_main_path + 'W+D-/subtraction/' + z_def + '/' + FF_set + '/pdf_errs/' + PDF_set + \
                    '/' + str(member_index) + '_' + str(eta_lept_index) + '_vals.txt')))

        WmDstarp_sigma = sum(sum(np.loadtxt(values_main_path + 'W-Dstar+/NLO/' + z_def + '/' + FF_set + '/pdf_errs/' + PDF_set + \
                    '/' + str(member_index) + '_' + str(eta_lept_index) + '_vals.txt'))) - \
                    sum(sum(np.loadtxt(values_main_path + 'W-Dstar+/subtraction/' + z_def + '/' + FF_set + '/pdf_errs/' + PDF_set + \
                    '/' + str(member_index) + '_' + str(eta_lept_index) + '_vals.txt')))

        WpDstarm_sigma = sum(sum(np.loadtxt(values_main_path + 'W+Dstar-/NLO/' + z_def + '/' + FF_set + '/pdf_errs/' + PDF_set + \
                    '/' + str(member_index) + '_' + str(eta_lept_index) + '_vals.txt'))) - \
                    sum(sum(np.loadtxt(values_main_path + 'W+Dstar-/subtraction/' + z_def + '/' + FF_set + '/pdf_errs/' + PDF_set + \
                    '/' + str(member_index) + '_' + str(eta_lept_index) + '_vals.txt')))

        Rcpm_val = (WpDm_sigma + WpDstarm_sigma) / (WmDp_sigma + WmDstarp_sigma)

        if (member_index % 2 == 0):
            Rcpm_minus_member = np.append(Rcpm_minus_member, Rcpm_val)
        else:
            Rcpm_plus_member = np.append(Rcpm_plus_member, Rcpm_val)

Rcpm = np.zeros((2, int(n_members / 2)))

for member_index in range(int(n_members / 2)):
    Rcpm[0, member_index] = Rcpm_minus_member[member_index]
    Rcpm[1, member_index] = Rcpm_plus_member[member_index]

np.savetxt('theory_values/eta_lept.txt')

#--------------------------------------------------------------------------------------------------------------------------------------#
#                                                           Rcpm as a function of pTD                                                  #
#--------------------------------------------------------------------------------------------------------------------------------------#

Rcpm_plus_member = np.zeros(int(n_members / 2))
Rcpm_minus_member = np.zeros(int(n_members / 2))

for member_index in range(n_members):

    WmDp_sigma = 0
    WpDm_sigma = 0
    WmDstarp_sigma = 0
    WpDstarm_sigma = 0

    bin_index = 0
    for pTD_index in range(284):
        if (pTD_data_min + (pTD_index + 1 / 2) * pTD_data_bin_width > pTD_bins[bin_index + 1]):
            if (bin_index == 5):
                break
            else:
                Rcpm_val = (WpDm_sigma + WpDstarm_sigma) / (WmDp_sigma + WmDstarp_sigma)

                if (member_index % 2 == 0):
                    Rcpm_minus_member = np.append(Rcpm_minus_member, Rcpm_val)
                else:
                    Rcpm_plus_member = np.append(Rcpm_plus_member, Rcpm_val)

                bin_index += 1

                WmDp_sigma = 0
                WpDm_sigma = 0
                WmDstarp_sigma = 0
                WpDstarm_sigma = 0

        for eta_lept_index in range(5):
                WmDp_sigma += sum(np.loadtxt(values_main_path + 'W-D+/NLO/' + z_def + '/' + FF_set + '/pdf_errs/' + PDF_set + \
                            '/' + str(member_index) + '_' + str(eta_lept_index) + '_vals.txt', delimiter=',')[pTD_index, :]) - \
                            sum(np.loadtxt(values_main_path + 'W-D+/subtraction/' + z_def + '/' + FF_set + '/pdf_errs/' + PDF_set + \
                            '/' + str(member_index) + '_' + str(eta_lept_index) + '_vals.txt', delimiter=',')[pTD_index, :])

                WpDm_sigma += sum(np.loadtxt(values_main_path + 'W+D-/NLO/' + z_def + '/' + FF_set + '/pdf_errs/' + PDF_set + \
                            '/' + str(member_index) + '_' + str(eta_lept_index) + '_vals.txt', delimiter=',')[pTD_index, :]) - \
                            sum(np.loadtxt(values_main_path + 'W+D-/subtraction/' + z_def + '/' + FF_set + '/pdf_errs/' + PDF_set + \
                            '/' + str(member_index) + '_' + str(eta_lept_index) + '_vals.txt', delimiter=',')[pTD_index, :])

                WmDstarp_sigma += sum(np.loadtxt(values_main_path + 'W-Dstar+/NLO/' + z_def + '/' + FF_set + '/pdf_errs/' + PDF_set + \
                            '/' + str(member_index) + '_' + str(eta_lept_index) + '_vals.txt', delimiter=',')[pTD_index, :]) - \
                            sum(np.loadtxt(values_main_path + 'W-Dstar+/subtraction/' + z_def + '/' + FF_set + '/pdf_errs/' + PDF_set + \
                            '/' + str(member_index) + '_' + str(eta_lept_index) + '_vals.txt', delimiter=',')[pTD_index, :])

                WpDstarm_sigma += sum(np.loadtxt(values_main_path + 'W+Dstar-/NLO/' + z_def + '/' + FF_set + '/pdf_errs/' + PDF_set + \
                            '/' + str(member_index) + '_' + str(eta_lept_index) + '_vals.txt', delimiter=',')[pTD_index, :]) - \
                            sum(np.loadtxt(values_main_path + 'W+Dstar-/subtraction/' + z_def + '/' + FF_set + '/pdf_errs/' + PDF_set + \
                            '/' + str(member_index) + '_' + str(eta_lept_index) + '_vals.txt', delimiter=',')[pTD_index, :])

Rcpm = np.zeros((2, int(n_members / 2)))

for member_index in range(int(n_members / 2)):
    Rcpm[0, member_index] = Rcpm_minus_member[member_index]
    Rcpm[1, member_index] = Rcpm_plus_member[member_index]

np.savetxt('theory_values/pTD.txt')




        
