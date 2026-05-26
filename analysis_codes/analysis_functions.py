import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import copy
import matplotlib as mpl


mpl.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.weight": "bold"
})
mpl.rcParams['text.latex.preamble'] = r'''
\usepackage{amsmath}
\usepackage{xcolor}
'''

# opal or global
fragmentation_set = 'KKKS08_opal'
# minus or plus
z_def = 'minus'

if (fragmentation_set == 'KKKS08_opal'):
    frag_set_text = 'KKKS08 OPAL'
elif (fragmentation_set == 'SMSKA19'):
    frag_set_text = 'SMSKA19'
else:
    frag_set_text = ' KKKS08 GLOBAL'

scale_var_color = 'cornflowerblue'
pdf_err_color = 'salmon'
atlas_err_color = 'lightgray'

PDF_sets = ['CT18ANLO', 'MSHT20nlo_as118', 'NNPDF40_nlo_pch_as_01180']
#PDF_sets = ['CT18ANLO', 'NNPDF40_nlo_pch_as_01180', 'MSHT20nlo_as118', 'NNPDF40_nlo_pch_as_01180']
num_err_members_in_sets = [58, 64, 100]

num_etac_bins = 11

main_vals_directory = '/home/alankovh/Documents/WD_production/output/'
plots_directory = '/home/alankovh/Documents/WD_production/plots/13 TeV/'
reweighting_input_directory = '/home/alankovh/Documents/WD_production/reweighting/input/'

marker_color = 'black'
theory_edge_colors = ['red', 'blue']
markers = ['d', 'v', 'o']
theory_labels = ['CT18ANLO', 'MSHT20NLO', 'NNPDF4.0NLO (pch)']

pdf_centrals = [[[np.zeros((284, num_etac_bins)) for _ in range(5)] for _ in range(len(PDF_sets))] for _ in range(4)]


def compute_normalized_3D_values_for_a_pdf_member(process, PDF_index, member_index):
    process_index = -1
    if (process == 'W-D+'):
        process_index = 0
    elif (process == 'W-Dstar+'):
        process_index = 1
    elif (process == 'W+D-'):
        process_index = 2
    else:
        process_index = 3

    # Here "normalized" means that the difference between the member value and the central value is added to the central value of the scale variation run.
    member_vals_normalized = [np.zeros((284, num_etac_bins)) for _ in range(5)]

    PDF_set = PDF_sets[PDF_index]

    for eta_lept_index in range(5):
        scale_var_central = np.loadtxt(main_vals_directory + process + '/NLO/' + z_def + '/' + \
                                                fragmentation_set + '/frag_main_scale/scale_variation/' + PDF_set + '/central/' + \
                                                str(eta_lept_index) + '_vals.txt', delimiter=',') - \
                                np.loadtxt(main_vals_directory + process + '/subtraction/' + z_def + '/' + \
                                                fragmentation_set + '/frag_main_scale/scale_variation/' + PDF_set + '/central/' + \
                                                str(eta_lept_index) + '_vals.txt', delimiter=',')

        pdf_central = np.zeros((284, num_etac_bins))

        if (PDF_set != "NNPDF40_nlo_pch_as_01180"):
            pdf_central = np.loadtxt(main_vals_directory + process + '/NLO/' + z_def + '/' + \
                                                fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                                str(0) + '_' + str(eta_lept_index) + '_vals.txt', delimiter=',') - \
                                np.loadtxt(main_vals_directory + process + '/subtraction/' + z_def + '/' + \
                                                fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                                str(0) + '_' + str(eta_lept_index) + '_vals.txt', delimiter=',')
            
        elif (PDF_set == 'NNPDF40_nlo_pch_as_01180'):
            if (member_index < 2):
                member_sum = np.zeros((284, num_etac_bins))

                for member_index_here in range(1, num_err_members_in_sets[PDF_index] + 1):
                    member_sum += np.loadtxt(main_vals_directory + process + '/NLO/' + z_def + '/' + \
                                                    fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                                    str(member_index_here) + '_' + str(eta_lept_index) + '_vals.txt', delimiter=',') - \
                                    np.loadtxt(main_vals_directory + process + '/subtraction/' + z_def + '/' + \
                                                    fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                                    str(member_index_here) + '_' + str(eta_lept_index) + '_vals.txt', delimiter=',')

                pdf_central = member_sum / (num_err_members_in_sets[PDF_index] * 1.)
                pdf_centrals[process_index][PDF_index][eta_lept_index] = pdf_central

            else:
                pdf_central = pdf_centrals[process_index][PDF_index][eta_lept_index]

        else:
            print('NOT IMPLEMENTED YET.')
            exit()

        member_vals = np.loadtxt(main_vals_directory + process + '/NLO/' + z_def + '/' + \
                                        fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                        str(member_index) + '_' + str(eta_lept_index) + '_vals.txt', delimiter=',') - \
                        np.loadtxt(main_vals_directory + process + '/subtraction/' + z_def + '/' + \
                                        fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                        str(member_index) + '_' + str(eta_lept_index) + '_vals.txt', delimiter=',')

        member_vals_normalized[eta_lept_index] = scale_var_central + member_vals - pdf_central

    return member_vals_normalized


def compute_Rcpm_pdf_err_HESSIAN(PDF_index, Rcpm_central, which_cross_sections_included):
    Rcpm_err_up = 0.
    Rcpm_err_down = 0.

    for member_index in range(1, int(num_err_members_in_sets[PDF_index] / 2) + 1):
        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W+D-', PDF_index, 2 * member_index - 1)
        WpDm_pdf_plus = sum(sum(sum(member_vals_normalized)))
        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W+D-', PDF_index, 2 * member_index)
        WpDm_pdf_minus = sum(sum(sum(member_vals_normalized)))

        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W+Dstar-', PDF_index, 2 * member_index - 1)
        WpDstarm_pdf_plus = sum(sum(sum(member_vals_normalized)))
        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W+Dstar-', PDF_index, 2 * member_index)
        WpDstarm_pdf_minus = sum(sum(sum(member_vals_normalized)))

        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W-D+', PDF_index, 2 * member_index - 1)
        WmDp_pdf_plus = sum(sum(sum(member_vals_normalized)))
        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W-D+', PDF_index, 2 * member_index)
        WmDp_pdf_minus = sum(sum(sum(member_vals_normalized)))

        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W-Dstar+', PDF_index, 2 * member_index - 1)
        WmDstarp_pdf_plus = sum(sum(sum(member_vals_normalized)))
        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W-Dstar+', PDF_index, 2 * member_index)
        WmDstarp_pdf_minus = sum(sum(sum(member_vals_normalized)))

        if (which_cross_sections_included == 'both'):
            Rcpm_plus = (WpDm_pdf_plus + WpDstarm_pdf_plus) / (WmDp_pdf_plus + WmDstarp_pdf_plus)
            Rcpm_minus = (WpDm_pdf_minus + WpDstarm_pdf_minus) / (WmDp_pdf_minus + WmDstarp_pdf_minus)

            #if (Rcpm_plus < 0.94):
                #(member_index)
            plt.plot(Rcpm_plus, 3 - PDF_index, marker='o', zorder=6, color='purple', linewidth=0.1)
            plt.plot(Rcpm_minus, 3 - PDF_index, marker='o', zorder=6, color='orange', linewidth=0.1)

            #print("member:", member_index,"plus:", Rcpm_plus, "minus:", Rcpm_minus, "central:", Rcpm_central)

        elif (which_cross_sections_included == 'D'):
            Rcpm_plus = (WpDm_pdf_plus) / (WmDp_pdf_plus)
            Rcpm_minus = (WpDm_pdf_minus) / (WmDp_pdf_minus)
        else:
            Rcpm_plus = (WpDstarm_pdf_plus) / (WmDstarp_pdf_plus)
            Rcpm_minus = (WpDstarm_pdf_minus) / (WmDstarp_pdf_minus)

        Rcpm_err_up += max(Rcpm_plus - Rcpm_central, Rcpm_minus - Rcpm_central, 0.)**2
        Rcpm_err_down += max(Rcpm_central - Rcpm_plus, Rcpm_central - Rcpm_minus, 0.)**2

    Rcpm_err_up = np.sqrt(Rcpm_err_up)
    Rcpm_err_down = np.sqrt(Rcpm_err_down)

    if (PDF_sets[PDF_index] == 'CT18ANLO'):
        Rcpm_err_up = Rcpm_err_up / 1.645
        Rcpm_err_down = Rcpm_err_down / 1.645

    return Rcpm_err_up, Rcpm_err_down


def compute_Rcpm_pdf_err_MC(PDF_index, Rcpm_central, which_cross_sections_included):
    average = 0.

    Rcpm_err_vals = np.zeros(num_err_members_in_sets[PDF_index])

    for member_index in range(1, num_err_members_in_sets[PDF_index] + 1):
        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W+D-', PDF_index, member_index)
        WpDm = sum(sum(sum(member_vals_normalized)))

        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W+Dstar-', PDF_index, member_index)
        WpDstarm = sum(sum(sum(member_vals_normalized)))

        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W-D+', PDF_index, member_index)
        WmDp = sum(sum(sum(member_vals_normalized)))

        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W-Dstar+', PDF_index, member_index)
        WmDstarp = sum(sum(sum(member_vals_normalized)))

        if (which_cross_sections_included == 'both'):
            Rcpm = (WpDm + WpDstarm) / (WmDp + WmDstarp)
        elif (which_cross_sections_included == 'D'):
            Rcpm = (WpDm) / (WmDp)
        else:
            Rcpm = (WpDstarm) / (WmDstarp)

        plt.plot(Rcpm, 3 - PDF_index, marker='o', zorder=6, color='purple', linewidth=0.1)
        #print(Rcpm)

        Rcpm_err_vals[member_index - 1] = Rcpm

        average += Rcpm
    
    average = average / (num_err_members_in_sets[PDF_index] * 1.)

    Rcpm_central = average

    sum_in_error_formula = 0.

    for i in range(len(Rcpm_err_vals)):
        sum_in_error_formula += (Rcpm_err_vals[i] - Rcpm_central)**2

    Rcpm_err_up = np.sqrt(1. / (num_err_members_in_sets[PDF_index] * 1. - 1.) * sum_in_error_formula)
    Rcpm_err_down = Rcpm_err_up

    return Rcpm_central, Rcpm_err_up, Rcpm_err_down


def compute_Rcpm_pdf_err_eta_lept_HESSIAN(PDF_index, Rcpm_central, which_cross_sections_included):
    Rcpm_err_up = np.zeros(5)
    Rcpm_err_down = np.zeros(5)

    Rcpm_vals_plus_member = np.zeros((int(num_err_members_in_sets[PDF_index] / 2), 5))
    Rcpm_vals_minus_member = np.zeros((int(num_err_members_in_sets[PDF_index] / 2), 5))

    for member_index in range(1, int(num_err_members_in_sets[PDF_index] / 2) + 1):
        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W+D-', PDF_index, 2 * (member_index - 1) + 1)
        WpDm_pdf_plus_3D_vals = copy.deepcopy(member_vals_normalized)
        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W+D-', PDF_index, 2 * member_index)
        WpDm_pdf_minus_3D_vals = copy.deepcopy(member_vals_normalized)

        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W+Dstar-', PDF_index, 2 * (member_index - 1) + 1)
        WpDstarm_pdf_plus_3D_vals = copy.deepcopy(member_vals_normalized)
        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W+Dstar-', PDF_index, 2 * member_index)
        WpDstarm_pdf_minus_3D_vals = copy.deepcopy(member_vals_normalized)

        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W-D+', PDF_index, 2 * (member_index - 1) + 1)
        WmDp_pdf_plus_3D_vals = copy.deepcopy(member_vals_normalized)
        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W-D+', PDF_index, 2 * member_index)
        WmDp_pdf_minus_3D_vals = copy.deepcopy(member_vals_normalized)

        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W-Dstar+', PDF_index, 2 * (member_index - 1) + 1)
        WmDstarp_pdf_plus_3D_vals = copy.deepcopy(member_vals_normalized)
        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W-Dstar+', PDF_index, 2 * member_index)
        WmDstarp_pdf_minus_3D_vals = copy.deepcopy(member_vals_normalized)

        for eta_lept_index in range(5):
            WpDm_pdf_plus = sum(sum(WpDm_pdf_plus_3D_vals[eta_lept_index]))
            WpDm_pdf_minus = sum(sum(WpDm_pdf_minus_3D_vals[eta_lept_index]))

            WpDstarm_pdf_plus = sum(sum(WpDstarm_pdf_plus_3D_vals[eta_lept_index]))
            WpDstarm_pdf_minus = sum(sum(WpDstarm_pdf_minus_3D_vals[eta_lept_index]))

            WmDp_pdf_plus = sum(sum(WmDp_pdf_plus_3D_vals[eta_lept_index]))
            WmDp_pdf_minus = sum(sum(WmDp_pdf_minus_3D_vals[eta_lept_index]))

            WmDstarp_pdf_plus = sum(sum(WmDstarp_pdf_plus_3D_vals[eta_lept_index]))
            WmDstarp_pdf_minus = sum(sum(WmDstarp_pdf_minus_3D_vals[eta_lept_index]))

            if (which_cross_sections_included == 'both'):
                Rcpm_plus = (WpDm_pdf_plus + WpDstarm_pdf_plus) / (WmDp_pdf_plus + WmDstarp_pdf_plus)
                Rcpm_minus = (WpDm_pdf_minus + WpDstarm_pdf_minus) / (WmDp_pdf_minus + WmDstarp_pdf_minus)
            elif (which_cross_sections_included == 'D'):
                Rcpm_plus = (WpDm_pdf_plus) / (WmDp_pdf_plus)
                Rcpm_minus = (WpDm_pdf_minus) / (WmDp_pdf_minus)
            else:
                Rcpm_plus = WpDstarm_pdf_plus / WmDstarp_pdf_plus
                Rcpm_minus = WpDstarm_pdf_minus / WmDstarp_pdf_minus

            Rcpm_vals_plus_member[member_index - 1, eta_lept_index] = Rcpm_plus
            Rcpm_vals_minus_member[member_index - 1, eta_lept_index] = Rcpm_minus

            Rcpm_err_up[eta_lept_index] += max(Rcpm_plus - Rcpm_central[eta_lept_index],
                                Rcpm_minus - Rcpm_central[eta_lept_index], 0)**2

            Rcpm_err_down[eta_lept_index] += max(Rcpm_central[eta_lept_index] - Rcpm_plus,
                                Rcpm_central[eta_lept_index] - Rcpm_minus, 0)**2

    np.savetxt(reweighting_input_directory + 'theory_values/HESSIAN/variation/eta_lept_' + which_cross_sections_included + '_' + \
                PDF_sets[PDF_index] + '_plus.txt', Rcpm_vals_plus_member.T, delimiter=',')
    np.savetxt(reweighting_input_directory + 'theory_values/HESSIAN/variation/eta_lept_' + which_cross_sections_included + '_' + \
                PDF_sets[PDF_index] + '_minus.txt', Rcpm_vals_minus_member.T, delimiter=',')

    Rcpm_err_up = np.sqrt(Rcpm_err_up)
    Rcpm_err_down = np.sqrt(Rcpm_err_down)

    if (PDF_sets[PDF_index] == 'CT18ANLO'):
        Rcpm_err_up = Rcpm_err_up / 1.645
        Rcpm_err_down = Rcpm_err_down / 1.645
    
    np.savetxt(reweighting_input_directory + 'theory_values/HESSIAN/errors/eta_lept_' + which_cross_sections_included + '_' + \
                PDF_sets[PDF_index] + '_plus.txt', Rcpm_err_up)
    np.savetxt(reweighting_input_directory + 'theory_values/HESSIAN/errors/eta_lept_' + which_cross_sections_included + '_' + \
                PDF_sets[PDF_index] + '_minus.txt', Rcpm_err_down)

    return Rcpm_err_up, Rcpm_err_down


def compute_Rcpm_pdf_err_eta_lept_MC(PDF_index, Rcpm_central, which_cross_sections_included):
    Rcpm_vals = np.zeros((num_err_members_in_sets[PDF_index], 5))
    Rcpm_best = np.zeros(5)

    PDF_set = PDF_sets[PDF_index]

    sum_in_error_formula = np.zeros(5)

    for member_index in range(1, num_err_members_in_sets[PDF_index] + 1):
        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W+D-', PDF_index, member_index)
        WpDm_3D_vals = copy.deepcopy(member_vals_normalized)

        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W+Dstar-', PDF_index, member_index)
        WpDstarm_3D_vals = copy.deepcopy(member_vals_normalized)

        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W-D+', PDF_index, member_index)
        WmDp_3D_vals = copy.deepcopy(member_vals_normalized)

        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W-Dstar+', PDF_index, member_index)
        WmDstarp_3D_vals = copy.deepcopy(member_vals_normalized)

        for eta_lept_index in range(5):
            WpDm = sum(sum(WpDm_3D_vals[eta_lept_index]))
            WpDstarm = sum(sum(WpDstarm_3D_vals[eta_lept_index]))
            WmDp = sum(sum(WmDp_3D_vals[eta_lept_index]))
            WmDstarp = sum(sum(WmDstarp_3D_vals[eta_lept_index]))

            if (which_cross_sections_included == 'both'):
                Rcpm = (WpDm + WpDstarm) / (WmDp + WmDstarp)
            elif (which_cross_sections_included == 'D'):
                Rcpm = WpDm / WmDp
            else:
                Rcpm = WpDstarm / WmDstarp

            Rcpm_vals[member_index - 1, eta_lept_index] = Rcpm
            
            sum_in_error_formula[eta_lept_index] += (Rcpm - Rcpm_central[eta_lept_index])**2

    np.savetxt(reweighting_input_directory + 'theory_values/MC/variation/eta_lept_' + which_cross_sections_included + '_' + \
                PDF_sets[PDF_index] + '_vals.txt', Rcpm_vals, delimiter=',')

    Rcpm_err_up = np.sqrt(1. / (num_err_members_in_sets[PDF_index] * 1. - 1.) * sum_in_error_formula)
    
    Rcpm_err_down = Rcpm_err_up

    for eta_lept_index in range(5):
        Rcpm_best[eta_lept_index] = sum(Rcpm_vals[:, eta_lept_index]) / (num_err_members_in_sets[PDF_index] * 1.)

    np.savetxt(reweighting_input_directory + 'theory_values/MC/errors/eta_lept_' + which_cross_sections_included + '_' + \
                PDF_sets[PDF_index] + '_plus.txt', Rcpm_err_up)
    np.savetxt(reweighting_input_directory + 'theory_values/MC/errors/eta_lept_' + which_cross_sections_included + '_' + \
                PDF_sets[PDF_index] + '_minus.txt', Rcpm_err_down)

    return Rcpm_best, Rcpm_err_up, Rcpm_err_down


def compute_Rcpm_pdf_err_pTD_HESSIAN(PDF_index, Rcpm_central, which_cross_sections_included):
    Rcpm_err_up = np.zeros(5)
    Rcpm_err_down = np.zeros(5)

    Rcpm_vals_plus_member = np.zeros((int(num_err_members_in_sets[PDF_index] / 2), 5))
    Rcpm_vals_minus_member = np.zeros((int(num_err_members_in_sets[PDF_index] / 2), 5))

    pTD_min = 8.
    pTD_bin_width = 0.5
    pTD_bins = [8., 12., 20., 40., 80., 150.]

    for member_index in range(1, int(num_err_members_in_sets[PDF_index] / 2) + 1):
        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W+D-', PDF_index, 2 * (member_index - 1) + 1)
        WpDm_3D_vals_plus = np.array(member_vals_normalized)
        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W+D-', PDF_index, 2 * member_index)
        WpDm_3D_vals_minus = np.array(member_vals_normalized)

        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W+Dstar-', PDF_index, 2 * (member_index - 1) + 1)
        WpDstarm_3D_vals_plus = np.array(member_vals_normalized)
        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W+Dstar-', PDF_index, 2 * member_index)
        WpDstarm_3D_vals_minus = np.array(member_vals_normalized)

        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W-D+', PDF_index, 2 * (member_index - 1) + 1)
        WmDp_3D_vals_plus = np.array(member_vals_normalized)
        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W-D+', PDF_index, 2 * member_index)
        WmDp_3D_vals_minus = np.array(member_vals_normalized)

        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W-Dstar+', PDF_index, 2 * (member_index - 1) + 1)
        WmDstarp_3D_vals_plus = np.array(member_vals_normalized)
        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W-Dstar+', PDF_index, 2 * member_index)
        WmDstarp_3D_vals_minus = np.array(member_vals_normalized)

        WpDm_plus = 0.
        WpDstarm_plus = 0.
        WmDp_plus = 0.
        WmDstarp_plus = 0.
        WpDm_minus = 0.
        WpDstarm_minus = 0.
        WmDp_minus = 0.
        WmDstarp_minus = 0.

        bin_index = 0

        for pTD_index in range(285):
            if (pTD_min + (pTD_index + 1. / 2.) * pTD_bin_width > pTD_bins[bin_index + 1]):
                if (which_cross_sections_included == 'both'):
                    Rcpm_plus = (WpDm_plus + WpDstarm_plus) / (WmDp_plus + WmDstarp_plus)
                    Rcpm_minus = (WpDm_minus + WpDstarm_minus) / (WmDp_minus + WmDstarp_minus)
                elif (which_cross_sections_included == 'D'):
                    Rcpm_plus = (WpDm_plus) / (WmDp_plus)
                    Rcpm_minus = (WpDm_minus) / (WmDp_minus)
                else:
                    Rcpm_plus = WpDstarm_plus / WmDstarp_plus
                    Rcpm_minus = WpDstarm_minus / WmDstarp_minus
                
                if (PDF_index == 1 and bin_index == 4):
                    print(Rcpm_central[4], Rcpm_plus, Rcpm_minus)
                            
                Rcpm_vals_plus_member[member_index - 1, bin_index] = Rcpm_plus
                Rcpm_vals_minus_member[member_index - 1, bin_index] = Rcpm_minus

                Rcpm_err_up[bin_index] += max(Rcpm_plus - Rcpm_central[bin_index],
                                    Rcpm_minus - Rcpm_central[bin_index],
                                    0)**2

                Rcpm_err_down[bin_index] += max(Rcpm_central[bin_index] - Rcpm_plus,
                                    Rcpm_central[bin_index] - Rcpm_minus,
                                    0)**2
                
                WpDm_plus = 0.
                WpDstarm_plus = 0.
                WmDp_plus = 0.
                WmDstarp_plus = 0.
                WpDm_minus = 0.
                WpDstarm_minus = 0.
                WmDp_minus = 0.
                WmDstarp_minus = 0.
                
                bin_index += 1

                if (bin_index == 5):
                    break
                
            WpDm_plus += WpDm_3D_vals_plus[:, pTD_index, :].sum()
            WpDm_minus += WpDm_3D_vals_minus[:, pTD_index, :].sum()

            WmDp_plus += WmDp_3D_vals_plus[:, pTD_index, :].sum()
            WmDp_minus += WmDp_3D_vals_minus[:, pTD_index, :].sum()

            WpDstarm_plus += WpDstarm_3D_vals_plus[:, pTD_index, :].sum()
            WpDstarm_minus += WpDstarm_3D_vals_minus[:, pTD_index, :].sum()

            WmDstarp_plus += WmDstarp_3D_vals_plus[:, pTD_index, :].sum()
            WmDstarp_minus += WmDstarp_3D_vals_minus[:, pTD_index, :].sum()
    
    np.savetxt(reweighting_input_directory + 'theory_values/HESSIAN/variation/pTD_' + which_cross_sections_included + \
                '_' + PDF_sets[PDF_index] + '_plus.txt', Rcpm_vals_plus_member.T, delimiter=',')
    np.savetxt(reweighting_input_directory + 'theory_values/HESSIAN/variation/pTD_' + which_cross_sections_included + \
                '_' + PDF_sets[PDF_index] + '_minus.txt', Rcpm_vals_minus_member.T, delimiter=',')

    Rcpm_err_up = np.sqrt(Rcpm_err_up)
    Rcpm_err_down = np.sqrt(Rcpm_err_down)

    np.savetxt(reweighting_input_directory + 'theory_values/HESSIAN/errors/pTD_' + which_cross_sections_included + '_' + \
                PDF_sets[PDF_index] + '_plus.txt', Rcpm_err_up)
    np.savetxt(reweighting_input_directory + 'theory_values/HESSIAN/errors/pTD_' + which_cross_sections_included + '_' + \
                PDF_sets[PDF_index] + '_minus.txt', Rcpm_err_down)

    return Rcpm_err_up, Rcpm_err_down


def compute_Rcpm_pdf_err_pTD_MC(PDF_index, Rcpm_central, which_cross_sections_included):
    Rcpm_vals = np.zeros((num_err_members_in_sets[PDF_index], 5))
    Rcpm_best = np.zeros(5)

    sum_in_error_formula = np.zeros(5)

    pTD_data_min = 8.
    pTD_data_bin_width = 0.5
    pTD_bins = [8., 12., 20., 40., 80., 150.]

    for member_index in range(1, num_err_members_in_sets[PDF_index] + 1):
        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W+D-', PDF_index, member_index)
        WpDm_3D_vals = copy.deepcopy(member_vals_normalized)

        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W+Dstar-', PDF_index, member_index)
        WpDstarm_3D_vals = copy.deepcopy(member_vals_normalized)

        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W-D+', PDF_index, member_index)
        WmDp_3D_vals = copy.deepcopy(member_vals_normalized)

        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W-Dstar+', PDF_index, member_index)
        WmDstarp_3D_vals = copy.deepcopy(member_vals_normalized)

        WpDm = 0.
        WpDstarm = 0.
        WmDp = 0.
        WmDstarp = 0.

        bin_index = 0
        for pTD_index in range(284 + 1):
            if (pTD_data_min + (pTD_index + 1 / 2) * pTD_data_bin_width > pTD_bins[bin_index + 1]):
                if (which_cross_sections_included == 'both'):
                    Rcpm = (WpDm + WpDstarm) / (WmDp + WmDstarp)
                elif (which_cross_sections_included == 'D'):
                    Rcpm = WpDm / WmDp
                else:
                    Rcpm = WpDstarm / WmDstarp
                
                Rcpm_vals[member_index - 1, bin_index] = Rcpm

                sum_in_error_formula[bin_index] += (Rcpm - Rcpm_central[bin_index])**2

                WpDm = 0.
                WpDstarm = 0.
                WmDp = 0.
                WmDstarp = 0.

                bin_index += 1

                if (bin_index == 5):
                    break

            for eta_lept_index in range(5):
                WpDm += sum(WpDm_3D_vals[eta_lept_index][pTD_index, :])
                WpDstarm += sum(WpDstarm_3D_vals[eta_lept_index][pTD_index, :])
                WmDp += sum(WmDp_3D_vals[eta_lept_index][pTD_index, :])
                WmDstarp += sum(WmDstarp_3D_vals[eta_lept_index][pTD_index, :])

    np.savetxt(reweighting_input_directory + 'theory_values/MC/variation/pTD_' + which_cross_sections_included + '_' + \
                PDF_sets[PDF_index] + '_vals.txt', Rcpm_vals, delimiter=',')

    Rcpm_err_up = np.sqrt(1. / (num_err_members_in_sets[PDF_index] * 1. - 1.) * sum_in_error_formula)
    
    Rcpm_err_down = Rcpm_err_up

    for pTD_index in range(5):
        Rcpm_best[pTD_index] = sum(Rcpm_vals[:, pTD_index]) / (num_err_members_in_sets[PDF_index] * 1.)

    np.savetxt(reweighting_input_directory + 'theory_values/MC/errors/pTD_' + which_cross_sections_included + '_' + \
                PDF_sets[PDF_index] + '_plus.txt', Rcpm_err_up)
    np.savetxt(reweighting_input_directory + 'theory_values/MC/errors/pTD_' + which_cross_sections_included + '_' + \
                PDF_sets[PDF_index] + '_minus.txt', Rcpm_err_down)

    return Rcpm_best, Rcpm_err_up, Rcpm_err_down


def compute_general_3D_vals_NLO(PDF_set, num_err_members_in_set, process, load_scale_variation, load_pdf_errs,
                                subtraction_flag, FF_scale_choice, z_def, fragmentation_set):
    central_vals = [np.zeros((284, num_etac_bins)) for _ in range(5)]
    dd_vals = [np.zeros((284, num_etac_bins)) for _ in range(5)]
    uu_vals = [np.zeros((284, num_etac_bins)) for _ in range(5)]
    scales_vals =[central_vals, dd_vals, uu_vals]
    pdf_err_plus = [np.zeros((284, num_etac_bins)) for _ in range(5)]
    pdf_err_minus = [np.zeros((284, num_etac_bins)) for _ in range(5)]

    central_MCerrs = [np.zeros((284, num_etac_bins)) for _ in range(5)]
    dd_MCerrs = [np.zeros((284, num_etac_bins)) for _ in range(5)]
    uu_MCerrs = [np.zeros((284, num_etac_bins)) for _ in range(5)]
    scales_MCerrs =[central_MCerrs, dd_MCerrs, uu_MCerrs]

    for eta_lept_index in range(5):
        for pTD_index in range(284):
            for etaD_index in range(num_etac_bins):
                pdf_err_plus[eta_lept_index][pTD_index, etaD_index] = 0
                pdf_err_minus[eta_lept_index][pTD_index, etaD_index] = 0

    stop = False

    scale_names = ['central', 'dd', 'uu']
    for eta_lept_index in range(5):
        # Get central and scale_variation values and MC errors. The central value will be corrected below for MC PDF sets.
        for scale_index in range(3):
            if (load_scale_variation or scale_index == 0):
                try:
                    scales_vals[scale_index][eta_lept_index] = np.loadtxt(
                        main_vals_directory + process + '/NLO/' + z_def + '/' + \
                        fragmentation_set + '/' + FF_scale_choice + '/scale_variation/' + PDF_set + '/' + scale_names[scale_index] + '/' + \
                        str(eta_lept_index) + '_vals.txt', delimiter=',')
                except FileNotFoundError:
                    stop = True
                    break
                
                scales_MCerrs[scale_index][eta_lept_index] = np.loadtxt(
                    main_vals_directory + process + '/NLO/' + z_def + '/' + \
                    fragmentation_set + '/' + FF_scale_choice + '/scale_variation/' + PDF_set + '/' + scale_names[scale_index] + '/' + \
                    str(eta_lept_index) + '_errs.txt', delimiter=',')

                if (subtraction_flag is True):
                    scales_vals[scale_index][eta_lept_index] -= np.loadtxt(
                        main_vals_directory + process + '/subtraction/' + z_def + '/' + \
                        fragmentation_set + '/' + FF_scale_choice + '/scale_variation/' + PDF_set + '/' + scale_names[scale_index] + '/' + \
                        str(eta_lept_index) + '_vals.txt', delimiter=',')

                    scales_MCerrs[scale_index][eta_lept_index] -= np.loadtxt(
                        main_vals_directory + process + '/subtraction/' + z_def + '/' + \
                        fragmentation_set + '/' + FF_scale_choice + '/scale_variation/' + PDF_set + '/' + scale_names[scale_index] + '/' + \
                        str(eta_lept_index) + '_errs.txt', delimiter=',')

            if (stop):
                break

        if (PDF_set == "NNPDF40_nlo_pch_as_01180"):
            average = np.zeros((284, num_etac_bins))
                    
            for member in range(1, int(num_err_members_in_set + 1)):
                average += np.loadtxt(main_vals_directory + process + '/NLO/' + z_def + '/' + \
                                                        fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                                        str(member) + '_' + str(eta_lept_index) + '_vals.txt', delimiter=',')
                
                if (subtraction_flag is True):
                    average -= np.loadtxt(main_vals_directory + process + '/subtraction/' + z_def + '/' + \
                                            fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                            str(member) + '_' + str(eta_lept_index) + '_vals.txt', delimiter=',')

            average = average / (num_err_members_in_set * 1.)

            pdf_errs_central = np.loadtxt(main_vals_directory + process + '/NLO/' + z_def + '/' + \
                                            fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                            str(0) + '_' + str(eta_lept_index) + '_vals.txt', delimiter=',')

            if (subtraction_flag is True):
                pdf_errs_central -= np.loadtxt(main_vals_directory + process + '/subtraction/' + z_def + '/' + \
                                            fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                            str(0) + '_' + str(eta_lept_index) + '_vals.txt', delimiter=',')
        
            scales_vals[0][eta_lept_index] = average + scales_vals[0][eta_lept_index] - pdf_errs_central
    
        # Get pdf err values
        if (load_pdf_errs):
            pdf_errs_central = np.loadtxt(main_vals_directory + process + '/NLO/' + z_def + '/' + \
                                            fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                            str(0) + '_' + str(eta_lept_index) + '_vals.txt', delimiter=',')
            
            if (subtraction_flag is True):
                pdf_errs_central -= np.loadtxt(main_vals_directory + process + '/subtraction/' + z_def + '/' + \
                                            fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                            str(0) + '_' + str(eta_lept_index) + '_vals.txt', delimiter=',')
            
            if (PDF_set == 'CT18NLO' or PDF_set == 'CT18ANLO' or PDF_set == 'MSHT20nlo_as118'):
                # This loops through error members, but isn't directly the member id.
                for i in range(1, int(num_err_members_in_set / 2 + 1)):
                    pdf_err_member_plus = np.zeros((284, num_etac_bins))
                    pdf_err_member_minus = np.zeros((284, num_etac_bins))

                    pdf_err_member_plus = np.loadtxt(main_vals_directory + process + '/NLO/' + z_def + '/' + \
                                            fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                            str(2 * (i - 1) + 1) + '_' + str(eta_lept_index) + '_vals.txt', delimiter=',')

                    pdf_err_member_minus = np.loadtxt(main_vals_directory + process + '/NLO/' + z_def + '/' + \
                                            fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                            str(2 * i) + '_' + str(eta_lept_index) + '_vals.txt', delimiter=',')

                    if (subtraction_flag is True):
                        pdf_err_member_plus -= np.loadtxt(main_vals_directory + process + '/subtraction/' + z_def + '/' + \
                                                fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                                str(2 * (i - 1) + 1) + '_' + str(eta_lept_index) + '_vals.txt', delimiter=',')
                        pdf_err_member_minus -= np.loadtxt(main_vals_directory + process + '/subtraction/' + z_def + '/' + \
                                                fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                                str(2 * i) + '_' + str(eta_lept_index) + '_vals.txt', delimiter=',')

                    for pTD_index in range(284):
                        for etaD_index in range(num_etac_bins):
                            pdf_err_plus[eta_lept_index][pTD_index, etaD_index] += max(
                                pdf_err_member_plus[pTD_index, etaD_index] - \
                                pdf_errs_central[pTD_index, etaD_index],
                                pdf_err_member_minus[pTD_index, etaD_index] - \
                                pdf_errs_central[pTD_index, etaD_index], 0)**2
                            
                            pdf_err_minus[eta_lept_index][pTD_index, etaD_index] += max(
                                pdf_errs_central[pTD_index, etaD_index] - \
                                pdf_err_member_plus[pTD_index, etaD_index],
                                pdf_errs_central[pTD_index, etaD_index] - \
                                pdf_err_member_minus[pTD_index, etaD_index], 0)**2

                pdf_err_plus[eta_lept_index] = np.sqrt(pdf_err_plus[eta_lept_index])
                pdf_err_minus[eta_lept_index] = np.sqrt(pdf_err_minus[eta_lept_index])

                if (PDF_set == 'CT18ANLO'):
                    pdf_err_plus[eta_lept_index] = pdf_err_plus[eta_lept_index] / 1.645
                    pdf_err_minus[eta_lept_index] = pdf_err_minus[eta_lept_index] / 1.645

            if (PDF_set == 'NNPDF40_nlo_pch_as_01180'):
                sum_val = np.zeros((284, num_etac_bins))

                for member in range(1, int(num_err_members_in_set + 1)):
                    pdf_err_member = np.loadtxt(main_vals_directory + process + '/NLO/' + z_def + '/' + \
                                                            fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                                            str(member) + '_' + str(eta_lept_index) + '_vals.txt', delimiter=',')
                                                        
                    if (subtraction_flag is True):
                        pdf_err_member -= np.loadtxt(main_vals_directory + process + '/subtraction/' + z_def + '/' + \
                                                    fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                                    str(member) + '_' + str(eta_lept_index) + '_vals.txt', delimiter=',')

                    sum_val += (average - pdf_err_member)**2

                pdf_err_plus[eta_lept_index] = np.sqrt(1. / (num_err_members_in_set * 1. - 1.) * sum_val)
                pdf_err_minus[eta_lept_index] = pdf_err_plus[eta_lept_index]
    
    return scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus

