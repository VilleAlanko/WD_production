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

variation_flag = True
without_subtraction_flag = False
# opal or global
fragmentation_set = 'KKKS08_opal'
# minus or plus
z_def = 'minus'

process = "W-D+"

PDF_sets = ['CT18ANLO', 'MSHT20nlo_as118', 'NNPDF40_nlo_pch_as_01180']
#PDF_sets = ['CT18ANLO', 'NNPDF40_nlo_pch_as_01180', 'MSHT20nlo_as118', 'NNPDF40_nlo_pch_as_01180']
num_err_members_in_sets = [58, 64, 100]

num_etac_bins = 11
pdf_centrals = [[[np.zeros((284, num_etac_bins)) for _ in range(5)] for _ in range(len(PDF_sets))] for _ in range(4)]

scale_var_color = 'cornflowerblue'
pdf_err_color = 'salmon'
atlas_err_color = 'lightgray'

if (fragmentation_set == 'KKKS08_opal'):
    frag_set_text = 'KKKS08 OPAL'
elif (fragmentation_set == 'SMSKA19'):
    frag_set_text = 'SMSKA19'
else:
    frag_set_text = ' KKKS08 GLOBAL'

if (process == "W+D-"):
    process_text = "$W^+D^-$"
    atlas_index = 1
elif (process == "W-D+"):
    process_text = "$W^-D^+$"
    atlas_index = 0
elif (process == "W+Dstar-"):
    process_text = "$W^+D^{*-}$"
    atlas_index = 3
else:
    process_text = "$W^-D^{*+}$"
    atlas_index = 2

QCD_orders = ['LO', 'NLO']
line_colors = ['#D81B60', 'blue']
bar_colors = ['peachpuff', 'lightskyblue']

marker_color = 'black'
theory_edge_colors = ['red', 'blue']
markers = ['d', 'v', 'o']
theory_labels = ['CT18ANLO', 'MSHT20NLO', 'NNPDF4.0NLO (pch)']

#colors = np.array([, "yellow", "#FFC107", "#004D40", "#1E88E5"])

main_vals_directory = '/home/alankovh/Documents/WD_production/output/'
plots_directory = '/home/alankovh/Documents/WD_production/plots/13 TeV/'
reweighting_input_directory = '/home/alankovh/Documents/WD_production/reweighting/input/'

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
                    print("!")
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
                average = 0.
                n_sum = 0
                for member in range(1, int(num_err_members_in_set + 1)):
                    average += np.loadtxt(main_vals_directory + process + '/NLO/' + z_def + '/' + \
                                                            fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                                            str(member) + '_' + str(eta_lept_index) + '_vals.txt', delimiter=',')
                    
                    if (subtraction_flag is True):
                        average -= np.loadtxt(main_vals_directory + process + '/subtraction/' + z_def + '/' + \
                                                fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                                str(member) + '_' + str(eta_lept_index) + '_vals.txt', delimiter=',')

                average = average / num_err_members_in_set
                
                scales_vals[0][eta_lept_index] = average + scales_vals[0][eta_lept_index] - pdf_errs_central

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

        PDF_set = PDF_sets[PDF_index]

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

        member_vals_normalized[eta_lept_index] = member_vals - pdf_central + scale_var_central
    
    return member_vals_normalized


def compute_Rcpm_pdf_err_HESSIAN(PDF_index, Rcpm_central, which_cross_sections_included):
    Rcpm_err_plus = 0.
    Rcpm_err_minus = 0.

    for member_index in range(1, int(num_err_members_in_sets[PDF_index] / 2) + 1):
        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W+D-', PDF_index, 2 * (member_index - 1) + 1)
        WpDm_pdf_plus = sum(sum(sum(member_vals_normalized)))
        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W+D-', PDF_index, 2 * member_index)
        WpDm_pdf_minus = sum(sum(sum(member_vals_normalized)))

        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W+Dstar-', PDF_index, 2 * (member_index - 1) + 1)
        WpDstarm_pdf_plus = sum(sum(sum(member_vals_normalized)))
        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W+Dstar-', PDF_index, 2 * member_index)
        WpDstarm_pdf_minus = sum(sum(sum(member_vals_normalized)))

        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W-D+', PDF_index, 2 * (member_index - 1) + 1)
        WmDp_pdf_plus = sum(sum(sum(member_vals_normalized)))
        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W-D+', PDF_index, 2 * member_index)
        WmDp_pdf_minus = sum(sum(sum(member_vals_normalized)))

        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W-Dstar+', PDF_index, 2 * (member_index - 1) + 1)
        WmDstarp_pdf_plus = sum(sum(sum(member_vals_normalized)))
        member_vals_normalized = compute_normalized_3D_values_for_a_pdf_member('W-Dstar+', PDF_index, 2 * member_index)
        WmDstarp_pdf_minus = sum(sum(sum(member_vals_normalized)))

        if (which_cross_sections_included == 'both'):
            Rcpm_plus = (WpDm_pdf_plus + WpDstarm_pdf_plus) / (WmDp_pdf_plus + WmDstarp_pdf_plus)
            Rcpm_minus = (WpDm_pdf_minus + WpDstarm_pdf_minus) / (WmDp_pdf_minus + WmDstarp_pdf_minus)
        elif (which_cross_sections_included == 'D'):
            Rcpm_plus = (WpDm_pdf_plus) / (WmDp_pdf_plus)
            Rcpm_minus = (WpDm_pdf_minus) / (WmDp_pdf_minus)
        else:
            Rcpm_plus = (WpDstarm_pdf_plus) / (WmDstarp_pdf_plus)
            Rcpm_minus = (WpDstarm_pdf_minus) / (WmDstarp_pdf_minus)

        Rcpm_err_plus += max(Rcpm_plus - Rcpm_central, Rcpm_central - Rcpm_minus, 0)**2
        Rcpm_err_minus += max(Rcpm_central - Rcpm_plus, Rcpm_minus - Rcpm_central, 0)**2

    Rcpm_err_plus = np.sqrt(Rcpm_err_plus)
    Rcpm_err_minus = np.sqrt(Rcpm_err_minus)

    if (PDF_sets[PDF_index] == 'CT18ANLO'):
        Rcpm_err_plus = Rcpm_err_plus / 1.645
        Rcpm_err_minus = Rcpm_err_minus / 1.645

    return Rcpm_err_plus, Rcpm_err_minus


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

        Rcpm_err_vals[member_index - 1] = Rcpm

        average += Rcpm
    
    average = average / (num_err_members_in_sets[PDF_index] * 1.)

    Rcpm_central = average

    sum_in_error_formula = 0.

    for i in range(len(Rcpm_err_vals)):
        sum_in_error_formula += (Rcpm_err_vals[i] - Rcpm_central)**2

    Rcpm_err_plus = np.sqrt(1. / (num_err_members_in_sets[PDF_index] * 1. - 1.) * sum_in_error_formula)
    Rcpm_err_minus = Rcpm_err_plus

    return Rcpm_central, Rcpm_err_plus, Rcpm_err_minus


def compute_Rcpm_pdf_err_eta_lept_HESSIAN(PDF_index, Rcpm_central, which_cross_sections_included):
    Rcpm_err_plus = np.zeros(5)
    Rcpm_err_minus = np.zeros(5)

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

            Rcpm_err_plus[eta_lept_index] += max(Rcpm_plus - Rcpm_central[eta_lept_index],
                                Rcpm_central[eta_lept_index] - Rcpm_minus, 0)**2

            Rcpm_err_minus[eta_lept_index] += max(Rcpm_central[eta_lept_index] - Rcpm_plus,
                                Rcpm_minus - Rcpm_central[eta_lept_index], 0)**2

    np.savetxt(reweighting_input_directory + 'theory_values/HESSIAN/variation/eta_lept_' + which_cross_sections_included + '_' + \
                PDF_sets[PDF_index] + '_plus.txt', Rcpm_vals_plus_member.T, delimiter=',')
    np.savetxt(reweighting_input_directory + 'theory_values/HESSIAN/variation/eta_lept_' + which_cross_sections_included + '_' + \
                PDF_sets[PDF_index] + '_minus.txt', Rcpm_vals_minus_member.T, delimiter=',')

    Rcpm_err_plus = np.sqrt(Rcpm_err_plus)
    Rcpm_err_minus = np.sqrt(Rcpm_err_minus)

    if (PDF_sets[PDF_index] == 'CT18ANLO'):
        Rcpm_err_plus = Rcpm_err_plus / 1.645
        Rcpm_err_minus = Rcpm_err_minus / 1.645
    
    np.savetxt(reweighting_input_directory + 'theory_values/HESSIAN/errors/eta_lept_' + which_cross_sections_included + '_' + \
                PDF_sets[PDF_index] + '_plus.txt', Rcpm_err_plus)
    np.savetxt(reweighting_input_directory + 'theory_values/HESSIAN/errors/eta_lept_' + which_cross_sections_included + '_' + \
                PDF_sets[PDF_index] + '_minus.txt', Rcpm_err_minus)

    return Rcpm_err_plus, Rcpm_err_minus


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

    Rcpm_err_plus = np.sqrt(1. / (num_err_members_in_sets[PDF_index] * 1. - 1.) * sum_in_error_formula)
    
    Rcpm_err_minus = Rcpm_err_plus

    for eta_lept_index in range(5):
        Rcpm_best[eta_lept_index] = sum(Rcpm_vals[:, eta_lept_index]) / (num_err_members_in_sets[PDF_index] * 1.)

    np.savetxt(reweighting_input_directory + 'theory_values/MC/errors/eta_lept_' + which_cross_sections_included + '_' + \
                PDF_sets[PDF_index] + '_plus.txt', Rcpm_err_plus)
    np.savetxt(reweighting_input_directory + 'theory_values/MC/errors/eta_lept_' + which_cross_sections_included + '_' + \
                PDF_sets[PDF_index] + '_minus.txt', Rcpm_err_minus)

    return Rcpm_best, Rcpm_err_plus, Rcpm_err_minus


def compute_Rcpm_pdf_err_pTD_HESSIAN(PDF_index, Rcpm_central, which_cross_sections_included):
    Rcpm_err_plus = np.zeros(5)
    Rcpm_err_minus = np.zeros(5)

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
                
                Rcpm_vals_plus_member[member_index - 1, bin_index] = Rcpm_plus
                Rcpm_vals_minus_member[member_index - 1, bin_index] = Rcpm_minus

                Rcpm_err_plus[bin_index] += max(Rcpm_plus - Rcpm_central[bin_index],
                                    Rcpm_central[bin_index] - Rcpm_minus,
                                    0)**2

                Rcpm_err_minus[bin_index] += max(Rcpm_central[bin_index] - Rcpm_plus,
                                    Rcpm_minus - Rcpm_central[bin_index],
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

    Rcpm_err_plus = np.sqrt(Rcpm_err_plus)
    Rcpm_err_minus = np.sqrt(Rcpm_err_minus)

    np.savetxt(reweighting_input_directory + 'theory_values/HESSIAN/errors/pTD_' + which_cross_sections_included + '_' + \
                PDF_sets[PDF_index] + '_plus.txt', Rcpm_err_plus)
    np.savetxt(reweighting_input_directory + 'theory_values/HESSIAN/errors/pTD_' + which_cross_sections_included + '_' + \
                PDF_sets[PDF_index] + '_minus.txt', Rcpm_err_minus)

    return Rcpm_err_plus, Rcpm_err_minus


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

    Rcpm_err_plus = np.sqrt(1. / (num_err_members_in_sets[PDF_index] * 1. - 1.) * sum_in_error_formula)
    
    Rcpm_err_minus = Rcpm_err_plus

    for pTD_index in range(5):
        Rcpm_best[pTD_index] = sum(Rcpm_vals[:, pTD_index]) / (num_err_members_in_sets[PDF_index] * 1.)

    np.savetxt(reweighting_input_directory + 'theory_values/MC/errors/pTD_' + which_cross_sections_included + '_' + \
                PDF_sets[PDF_index] + '_plus.txt', Rcpm_err_plus)
    np.savetxt(reweighting_input_directory + 'theory_values/MC/errors/pTD_' + which_cross_sections_included + '_' + \
                PDF_sets[PDF_index] + '_minus.txt', Rcpm_err_minus)

    return Rcpm_best, Rcpm_err_plus, Rcpm_err_minus


def pTD_plot(PDF_sets, plot_errors_flag, theory_labels):
    # Atlas values. The rows from top to bottom are D+W-, D-W+, D*+W-, D*-W+.
    atlas_vals = np.array([[15.04, 15.34, 13.78, 5.13, 0.93],
                        [14.61, 15.12, 13.07, 4.84, 0.82],
                        [14.50, 15.88, 14.19, 5.42, 1.07],
                        [14.26, 15.60, 14.08, 5.11, 0.99]])

    atlas_vals_up_err = np.array([[0.19+0.76, 0.14+0.78, 0.12+0.92, 0.07+0.34, 0.04+0.09],
                                    [0.19+0.73, 0.15+0.75, 0.12+0.89, 0.07+0.31, 0.04+0.08],
                                    [0.26+0.85, 0.19+0.73, 0.16+0.68, 0.10+0.31, 0.05+0.10],
                                    [0.27+0.82, 0.20+0.74, 0.17+0.68, 0.10+0.30, 0.06+0.09]])

    atlas_vals_down_err = np.array([[0.19+0.72, 0.14+0.75, 0.12+0.85, 0.07+0.31, 0.04+0.08],
                                    [0.19+0.69, 0.15+0.72, 0.12+0.82, 0.07+0.29, 0.04+0.07],
                                    [0.26+0.79, 0.19+0.69, 0.16+0.64, 0.10+0.29, 0.05+0.09],
                                    [0.27+0.76, 0.20+0.70, 0.17+0.64, 0.10+0.28, 0.06+0.08]])

    font_size = 17
    axis_label_font_size = 21
    axis_font_size = 15
    legend_fontsize = 13

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1], 'hspace': 0}, figsize=(6, 6))

    pTD_bins = np.array([8., 12., 20., 40., 80., 150.])

    pTD_data_min = 8.
    pTD_data_bin_width = 0.5

    bin_widths = np.diff(pTD_bins)

    bin_midpoints = np.zeros(len(pTD_bins) - 1)
    bin_midpoints_linear_scale = np.zeros(len(pTD_bins) - 1)

    for i in range(len(pTD_bins) - 1):
        bin_midpoints[i] = np.sqrt(pTD_bins[i] * pTD_bins[i + 1])
        bin_midpoints_linear_scale[i] = (pTD_bins[i + 1] + pTD_bins[i]) / 2

    places_inside_bins_log_scale = np.zeros((len(PDF_sets), 5))
    places_inside_bins_right_log_scale = np.zeros((len(PDF_sets), 5))
    places_inside_bins_left_log_scale = np.zeros((len(PDF_sets), 5))
    places_inside_bins_linear_scale = np.zeros((len(PDF_sets), 5))

    for i in range(len(PDF_sets)):
        places_inside_bins_log_scale[i, :] = pTD_bins[:-1] * (pTD_bins[1:] / pTD_bins[:-1])**((i * 1. + 1) / (len(PDF_sets) * 1. + 1.))
        places_inside_bins_linear_scale[i, :] = pTD_bins[:-1] + (i * 1. + 1.) / (len(PDF_sets) * 1. + 1.) * (pTD_bins[1:] - pTD_bins[:-1])

    bar_widths = np.zeros((3, 5))

    width_parameter = np.array([1.023, 1.03, 1.04, 1.04, 1.035])
    bar_widths = places_inside_bins_log_scale * width_parameter - places_inside_bins_log_scale / width_parameter

    for i in range(len(PDF_sets)):
        places_inside_bins_left_log_scale[i, :] = places_inside_bins_log_scale[i, :] / width_parameter
        places_inside_bins_right_log_scale[i, :] = places_inside_bins_log_scale[i, :] * width_parameter

    for PDF_index in range(len(PDF_sets)):
        for QCD_order_index in range(1, 2):
            QCD_order = QCD_orders[QCD_order_index]
            PDF_set = PDF_sets[PDF_index]
            num_err_members_in_set = num_err_members_in_sets[PDF_index]

            if (plot_errors_flag):
                scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus = compute_general_3D_vals_NLO(PDF_set, num_err_members_in_set,
                                                                            process, True, True, True, 'frag_main_scale', z_def, fragmentation_set)
            else:
                scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus = compute_general_3D_vals_NLO(PDF_set, num_err_members_in_set,
                                                                            process, False, False, True, 'frag_main_scale', z_def, fragmentation_set)

            HISTO_central_sigma_vals = np.zeros(5)
            HISTO_scales_dd_sigma_vals = np.zeros(5)
            HISTO_scales_uu_sigma_vals = np.zeros(5)

            HISTO_central_sigma_MCerrs = np.zeros(5)
            HISTO_scales_dd_sigma_MCerrs = np.zeros(5)
            HISTO_scales_uu_sigma_MCerrs = np.zeros(5)

            HISTO_pdf_err_plus = np.zeros(5)
            HISTO_pdf_err_minus = np.zeros(5)

            scale_names = ['central', 'dd', 'uu']

            for eta_lept_index in range(5):
                bin_index = 0
                for pTD_index in range(284 + 1):
                    if (pTD_data_min + (pTD_index + 1 / 2) * pTD_data_bin_width > pTD_bins[bin_index + 1]):
                        if (bin_index < 4):
                            bin_index += 1
                        else:
                            break

                    HISTO_central_sigma_vals[bin_index] += sum(scales_vals[0][eta_lept_index][pTD_index, :])
                    HISTO_scales_dd_sigma_vals[bin_index] += sum(scales_vals[1][eta_lept_index][pTD_index, :])
                    HISTO_scales_uu_sigma_vals[bin_index] += sum(scales_vals[2][eta_lept_index][pTD_index, :])

                    HISTO_central_sigma_MCerrs[bin_index] += sum(scales_MCerrs[0][eta_lept_index][pTD_index, :])
                    HISTO_scales_dd_sigma_MCerrs[bin_index] += sum(scales_MCerrs[1][eta_lept_index][pTD_index, :])
                    HISTO_scales_uu_sigma_MCerrs[bin_index] += sum(scales_MCerrs[2][eta_lept_index][pTD_index, :])

                    HISTO_pdf_err_plus[bin_index] += sum(pdf_err_plus[eta_lept_index][pTD_index, :])
                    HISTO_pdf_err_minus[bin_index] += sum(pdf_err_minus[eta_lept_index][pTD_index, :])

            if (variation_flag):
                if (QCD_order == "NLO"):
                    print(PDF_set, HISTO_central_sigma_vals)
                    ax1.plot(places_inside_bins_log_scale[PDF_index, :], HISTO_central_sigma_vals, marker=markers[PDF_index],
                                color=marker_color, markersize=5, linestyle='none',
                                label=theory_labels[PDF_index], zorder=4)

                    if (plot_errors_flag):
                        ax1.bar(places_inside_bins_left_log_scale[PDF_index, :], HISTO_pdf_err_plus + HISTO_pdf_err_minus, width=bar_widths[PDF_index],
                                    bottom=HISTO_central_sigma_vals - HISTO_pdf_err_minus, color=pdf_err_color,
                                    zorder=3)

                        ax1.bar(places_inside_bins_right_log_scale[PDF_index, :], HISTO_scales_uu_sigma_vals - HISTO_scales_dd_sigma_vals, width=bar_widths[PDF_index],
                                    bottom=HISTO_central_sigma_vals + HISTO_scales_dd_sigma_vals, color=scale_var_color,
                                    zorder=2)

            else:
                #ax1.hlines(HISTO_central_sigma_vals, pTD_bins[0:-1], pTD_bins[1:], color=line_colors[QCD_order_index], zorder=5, label='LO')
                print('NOT IMPLEMENTED YET')
                exit()

            print(HISTO_central_sigma_vals)

        if (QCD_order == 'NLO'):
            ratios = np.zeros(5)
            ratios_ATLAS_var_down = np.zeros(5)
            ratios_ATLAS_var_up = np.zeros(5)
            ratios_theory_scale_var_var_down = np.zeros(5)
            ratios_theory_scale_var_var_up = np.zeros(5)
            ratios_theory_pdf_err_var_down = np.zeros(5)
            ratios_theory_pdf_err_var_up = np.zeros(5)

            for eta_lept_index in range(5):
                ratios[eta_lept_index] = HISTO_central_sigma_vals[eta_lept_index] / atlas_vals[atlas_index][eta_lept_index]

                ratios_ATLAS_var_down[eta_lept_index] = HISTO_central_sigma_vals[eta_lept_index] / \
                                                        (atlas_vals[atlas_index][eta_lept_index] + atlas_vals_up_err[atlas_index][eta_lept_index])
                                                    
                ratios_ATLAS_var_up[eta_lept_index] = HISTO_central_sigma_vals[eta_lept_index] / \
                                                        (atlas_vals[atlas_index][eta_lept_index] - atlas_vals_down_err[atlas_index][eta_lept_index])

                if (plot_errors_flag):
                    ratios_theory_scale_var_var_down[eta_lept_index] = (HISTO_central_sigma_vals[eta_lept_index] + HISTO_scales_dd_sigma_vals[eta_lept_index]) / \
                                                                        HISTO_central_sigma_vals[eta_lept_index]
                                                        
                    ratios_theory_scale_var_var_up[eta_lept_index] = (HISTO_central_sigma_vals[eta_lept_index] + HISTO_scales_uu_sigma_vals[eta_lept_index]) / \
                                                                        HISTO_central_sigma_vals[eta_lept_index]
                    
                    ratios_theory_pdf_err_var_down[eta_lept_index] = (HISTO_central_sigma_vals[eta_lept_index] - HISTO_pdf_err_minus[eta_lept_index]) / \
                                                                        HISTO_central_sigma_vals[eta_lept_index]
                                                        
                    ratios_theory_pdf_err_var_up[eta_lept_index] = (HISTO_central_sigma_vals[eta_lept_index] + HISTO_pdf_err_plus[eta_lept_index]) / \
                                                                        HISTO_central_sigma_vals[eta_lept_index]
                                                        
            ax2.plot(places_inside_bins_log_scale[PDF_index, :], ratios, zorder=4, marker=markers[PDF_index],
                        color=marker_color, markersize=5, linestyle='none')

            if (plot_errors_flag):
                scale_var_bar_plot = ax2.bar(places_inside_bins_right_log_scale[PDF_index, :], ratios_theory_scale_var_var_up - ratios_theory_scale_var_var_down,
                        bottom=ratios - 1. + ratios_theory_scale_var_var_down, width=bar_widths[PDF_index], color=scale_var_color,
                        linewidth=1, zorder=3)
                
                pdf_err_bar_plot = ax2.bar(places_inside_bins_left_log_scale[PDF_index, :], ratios_theory_pdf_err_var_up - ratios_theory_pdf_err_var_down,
                        bottom=ratios - 1. + ratios_theory_pdf_err_var_down, width=bar_widths[PDF_index], color=pdf_err_color,
                        linewidth=1, zorder=3)
    
    ax1.hlines(atlas_vals[atlas_index], pTD_bins[:-1], pTD_bins[1:], color='black', label='ATLAS', zorder=1)
    #ax1.hlines(np.array([12.37307495, 13.34427756, 11.75877281,  4.27762968,  0.82701554]), pTD_bins[:-1], pTD_bins[1:], color='orange', zorder=10)
    #ax2.hlines(np.array([12.37307495, 13.34427756, 11.75877281,  4.27762968,  0.82701554]) / atlas_vals[atlas_index], pTD_bins[:-1], pTD_bins[1:], color='orange', zorder=10)
    # Plot error bars for the Atlas values.
    ATLAS_uncertainty = ax1.bar(bin_midpoints_linear_scale, atlas_vals_down_err[atlas_index] + atlas_vals_up_err[atlas_index],
            bottom=atlas_vals[atlas_index] - atlas_vals_down_err[atlas_index], color=atlas_err_color,
            width=bin_widths, zorder=0, linewidth=1.5)
    ax2.bar(bin_midpoints_linear_scale, ratios_ATLAS_var_up - ratios_ATLAS_var_down,
                    bottom=1. - ratios + ratios_ATLAS_var_down, color=atlas_err_color,
                    width=bin_widths, zorder=0, linewidth=1.5)

    plt.xscale('log', base=10)
    ax1.set_ylim(0, 30)
    plt.xlim(8, 150)
    ax2.set_ylim(0.65, 1.15)

    ax2.set_xlabel(r'$p_T (D)$ [GeV]', fontsize=axis_label_font_size - 2)
    ax1.set_ylabel(r'Cross section [pb]', fontsize=axis_label_font_size)
    ax2.set_ylabel(r'$\frac{\mathrm{Theory}}{\mathrm{ATLAS}}$', 
               fontsize=axis_label_font_size * 1.3)

    ax1.set_yticks([5, 10, 15, 20, 25, 30])
    ax2.set_yticks([0.7, 0.8, 0.9, 1., 1.1])

    # Configure ticks to appear on all sides
    ax1.tick_params(direction='in', top=True, right=True)

    # Add minor ticks
    ax1.minorticks_on()
    ax1.tick_params(which='both', direction='in', top=True, right=True)

    # Add minor ticks
    ax2.minorticks_on()
    ax2.tick_params(which='both', direction='in', top=True, right=True)

    # Configure ticks to appear on all sides
    ax2.tick_params(direction='in', top=True, right=True)

    ax2.plot([8, 150], [1, 1], linewidth=1, color='black', zorder=1)

    legend1 = ax1.legend(loc='upper right', framealpha=1, fontsize=legend_fontsize + 1)
    if (plot_errors_flag):
        legend2 = ax1.legend([ATLAS_uncertainty, pdf_err_bar_plot, scale_var_bar_plot],
                            ["ATLAS error", "PDF error (68\% C.L.)", "Scale variation"], loc='lower left', framealpha=1, fontsize=legend_fontsize + 0.5)
        ax1.add_artist(legend1)

    info_y_vals_1 = 26.5
    info_y_vals_2 = 23.5
    info_y_vals_3 = 20.5

    info_x_vals_1 = 9
    info_x_vals_2 = 25

    ax1.text(info_x_vals_1, info_y_vals_1, process_text + '  OS-SS', fontsize=font_size)
    ax1.text(info_x_vals_1, info_y_vals_2, r'$\sqrt{s} = 13$ TeV', fontsize=font_size)
    ax1.text(info_x_vals_1, info_y_vals_3, frag_set_text, fontsize=font_size)

    ax1.tick_params(axis='both', which='major', labelsize=axis_font_size)
    ax2.tick_params(axis='both', which='major', labelsize=axis_font_size)

    plt.xticks(pTD_bins, [f'{tick:.0f}' for tick in pTD_bins])

    for i in range(1, len(pTD_bins) - 1):
        ax1.axvline(pTD_bins[i], color='gray', linewidth=0.5, ymax=0.57, zorder=0)

    for i in range(1, len(pTD_bins) - 1):
        ax2.axvline(pTD_bins[i], color='gray', linewidth=0.5, ymax=1, zorder=0)

    plt.tight_layout()
    if (len(PDF_sets) == 3):
        plt.savefig(plots_directory + process + '/' + fragmentation_set + "/" + process + '_pT.pdf')
    else:
        plt.savefig(plots_directory + process + '/' + fragmentation_set + "/" + process + '_pT_nlo_vs_nnlo_pdf.pdf')
    plt.show()


def pTD_effect_of_subtraction_plot():
    for PDF_index in range(1):
        PDF_set = PDF_sets[PDF_index]
        num_err_members_in_set = num_err_members_in_sets[PDF_index]

        # Atlas values. The rows from top to bottom are D+W-, D-W+, D*+W-, D*-W+.
        atlas_vals = np.array([[15.04, 15.34, 13.78, 5.13, 0.93],
                            [14.61, 15.12, 13.07, 4.84, 0.82],
                            [14.50, 15.88, 14.19, 5.42, 1.07],
                            [14.26, 15.60, 14.08, 5.11, 0.99]])

        atlas_vals_up_err = np.array([[0.19+0.76, 0.14+0.78, 0.12+0.92, 0.07+0.34, 0.04+0.09],
                                        [0.19+0.73, 0.15+0.75, 0.12+0.89, 0.07+0.31, 0.04+0.08],
                                        [0.26+0.85, 0.19+0.73, 0.16+0.68, 0.10+0.31, 0.05+0.10],
                                        [0.27+0.82, 0.20+0.74, 0.17+0.68, 0.10+0.30, 0.06+0.09]])

        atlas_vals_down_err = np.array([[0.19+0.72, 0.14+0.75, 0.12+0.85, 0.07+0.31, 0.04+0.08],
                                        [0.19+0.69, 0.15+0.72, 0.12+0.82, 0.07+0.29, 0.04+0.07],
                                        [0.26+0.79, 0.19+0.69, 0.16+0.64, 0.10+0.29, 0.05+0.09],
                                        [0.27+0.76, 0.20+0.70, 0.17+0.64, 0.10+0.28, 0.06+0.08]])

        font_size = 17
        axis_label_font_size = 19
        axis_font_size = 14
        legend_fontsize = 13

        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1], 'hspace': 0}, figsize=(6, 6))

        pTD_bins = np.array([8, 12, 20, 40, 80, 150])

        pTD_data_min = 8.
        pTD_data_bin_width = 0.5

        bin_widths = np.diff(pTD_bins)

        places_inside_bins = np.zeros((2, 5))
        places_inside_bins_right = np.zeros((2, 5))
        places_inside_bins_left = np.zeros((2, 5))

        places_inside_bins[0, :] = pTD_bins[:-1]**(2 / 3) * pTD_bins[1:]**(1 / 3)
        places_inside_bins[1, :] = pTD_bins[:-1]**(1 / 3) * pTD_bins[1:]**(2 / 3)

        bar_width_over_bin_width = 1. / 9.
        bar_width = bar_width_over_bin_width * bin_widths
        shifts = np.zeros(5)

        scalings = [0.93, 1.1]

        for i in range(5):
            shifts[i] = ((pTD_bins[i + 1] * 1.) / (pTD_bins[i] * 1.))**(bar_width_over_bin_width / 2.)

        for i in range(2):
            places_inside_bins_left[i, :] = places_inside_bins[i, :] / shifts
            places_inside_bins_right[i, :] = places_inside_bins[i, :] * shifts

        bin_midpoints = np.zeros(len(pTD_bins) - 1)
        bin_midpoints_linear_scale = np.zeros(len(pTD_bins) - 1)

        for i in range(len(pTD_bins) - 1):
            bin_midpoints[i] = np.sqrt(pTD_bins[i] * pTD_bins[i + 1])
            bin_midpoints_linear_scale[i] = (pTD_bins[i + 1] + pTD_bins[i]) / 2
        
        xmin = np.sqrt(pTD_bins[0:-1] * bin_midpoints)
        xmax = np.sqrt(pTD_bins[1:] * bin_midpoints)

        QCD_order = 'NLO'
        
        subtraction_flags = [True, False, True, False]
        FF_scale_choices = ['frag_main_scale', 'frag_main_scale', 'frag_initial_scale', 'frag_initial_scale']
        subtraction_flag_colors = ['red', 'blue', 'black', 'orange']
        subtraction_flag_labels = [r'With subtraction' + '\n' + r'($\mu_\text{frag} = M_W$)', r'Without subtraction' + '\n' + r'($\mu_\text{frag} = M_W$)',
                                    r'With subtraction' + '\n' + r'($\mu_\text{frag} = m_c$)', r'Without subtraction' + '\n' + r'($\mu_\text{frag} = m_c$)']

        without_subtraction_vals = np.zeros(5)
        with_subtraction_vals = np.zeros(5)
        initial_scale_with_subtraction = np.zeros(5)
        initial_scale_without_subtraction = np.zeros(5)


        for i in range(len(FF_scale_choices)):
            subtraction_flag = subtraction_flags[i]
            FF_scale_choice = FF_scale_choices[i]

            scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus = compute_general_3D_vals_NLO(PDF_set,
                    num_err_members_in_set, process, True, True, subtraction_flag, FF_scale_choice, z_def, fragmentation_set)

            HISTO_central_sigma_vals = np.zeros(5)

            scale_names = ['central', 'dd', 'uu']

            for eta_lept_index in range(5):
                bin_index = 0
                for pTD_index in range(284):
                    if (pTD_data_min + (pTD_index + 1 / 2) * pTD_data_bin_width > pTD_bins[bin_index + 1]):
                        if (bin_index == 5):
                            break
                        else:
                            bin_index += 1
                    
                    HISTO_central_sigma_vals[bin_index] += sum(scales_vals[0][eta_lept_index][pTD_index, :])

            if (i == 0):
                with_subtraction_vals = HISTO_central_sigma_vals
            elif (i == 1):
                without_subtraction_vals = HISTO_central_sigma_vals
            elif(i == 2):
                initial_scale_with_subtraction = HISTO_central_sigma_vals
            else:
                initial_scale_without_subtraction = HISTO_central_sigma_vals

            ax1.hlines(HISTO_central_sigma_vals, pTD_bins[:-1], pTD_bins[1:],
                        color=subtraction_flag_colors[i],
                        label=subtraction_flag_labels[i], zorder=2)

        ratios1 = initial_scale_with_subtraction / with_subtraction_vals
        ratios2 = without_subtraction_vals / with_subtraction_vals
        ratios3 = initial_scale_without_subtraction / with_subtraction_vals

        ax2.hlines(ratios1, pTD_bins[:-1], pTD_bins[1:], zorder=3, color=subtraction_flag_colors[2])
        ax2.hlines(ratios2, pTD_bins[:-1], pTD_bins[1:], zorder=3, color=subtraction_flag_colors[1])
        ax2.hlines(ratios3, pTD_bins[:-1], pTD_bins[1:], zorder=3, color=subtraction_flag_colors[3])

        ax2.hlines(np.zeros(5) + 1, pTD_bins[:-1], pTD_bins[1:], zorder=3, color=subtraction_flag_colors[2])

        #ax1.hlines(atlas_vals[atlas_index], np.sqrt(bin_midpoints * (bin_midpoints - bin_widths / 15.)),
        #            np.sqrt(bin_midpoints * (bin_midpoints + bin_widths / 15.)), zorder=3, color='black')

        #ax1.bar(bin_midpoints, atlas_vals_down_err[atlas_index] + atlas_vals_up_err[atlas_index],
        #        bottom=atlas_vals[atlas_index] - atlas_vals_down_err[atlas_index], width=bar_width, zorder=2, color=atlas_err_color, edgecolor='black',
        #        linewidth=1)

        #ax1.hlines(atlas_vals[atlas_index], pTD_bins[:-1], pTD_bins[1:], color='black', label='ATLAS', zorder=1)
        # Plot error bars for the Atlas values.
        #ax1.bar(bin_midpoints_linear_scale, atlas_vals_down_err[atlas_index] + atlas_vals_up_err[atlas_index],
        #        bottom=atlas_vals[atlas_index] - atlas_vals_down_err[atlas_index], color=atlas_err_color,
        #        width=bin_widths, zorder=0, linewidth=1.5)

        plt.xscale('log', base=10)
        ax1.set_ylim(0, 50)
        plt.xlim(8, 150)
        ax2.set_ylim(0.5, 2.4)

        ax2.set_xlabel(r'$p_T(D)$ [GeV]', fontsize=axis_label_font_size)
        ax1.set_ylabel('Cross section [pb]', fontsize=axis_label_font_size)
        ax2.set_ylabel('Ratio', fontsize=axis_label_font_size)

        ax1.set_yticks([10, 20, 30, 40, 50])
        #ax2.set_yticks([0.25, 0.5, 0.75, 1])

        # Configure ticks to appear on all sides
        ax1.tick_params(direction='in', top=True, right=True)

        # Add minor ticks
        ax1.minorticks_on()
        ax1.tick_params(which='both', direction='in', top=True, right=True)

        # Add minor ticks
        ax2.minorticks_on()
        ax2.tick_params(which='both', direction='in', top=True, right=True)

        # Configure ticks to appear on all sides
        ax2.tick_params(direction='in', top=True, right=True)

        ax2.plot([8, 150], [1, 1], color='red', zorder=5)

        legend1 = ax1.legend(loc='upper right', framealpha=1, fontsize=legend_fontsize)
        #legend2 = ax1.legend([pdf_err_bar_plot, scale_var_bar_plot], ["PDF uncertainty", "Scale variation"], loc='lower left', framealpha=1, fontsize=legend_fontsize)
        ax1.add_artist(legend1)

        info_y_vals_1 = 43.5
        info_y_vals_2 = 38.5
        info_y_vals_3 = 33.5
        info_y_vals_4 = 28.5

        info_x_vals_1 = 9
        info_x_vals_2 = 25

        ax1.text(info_x_vals_1, info_y_vals_1, process_text, fontsize=font_size)
        ax1.text(info_x_vals_1, info_y_vals_2, r'$\sqrt{s} = 13$ TeV', fontsize=font_size)
        ax1.text(info_x_vals_1, info_y_vals_3, theory_labels[PDF_index], fontsize=font_size)
        ax1.text(info_x_vals_1, info_y_vals_4, frag_set_text, fontsize=font_size)

        ax1.tick_params(axis='both', which='major', labelsize=axis_font_size)
        ax2.tick_params(axis='both', which='major', labelsize=axis_font_size)

        plt.xticks(pTD_bins, [f'{tick:.0f}' for tick in pTD_bins])

        for i in range(1, len(pTD_bins) - 1):
            ax1.axvline(pTD_bins[i], color='gray', linewidth=0.5, ymax=0.5, zorder=0)

        for i in range(1, len(pTD_bins) - 1):
            ax2.axvline(pTD_bins[i], color='gray', linewidth=0.5, ymax=1, zorder=0)

        plt.tight_layout()
        plt.savefig(plots_directory + process + '/' + fragmentation_set + "/" + process + '_pT_subtr_' + theory_labels[PDF_index] + '.pdf')
        plt.show()


def pTD_dynamic_FF_scale(process, PDF_set, num_err_members_in_set, fragmentation_set, z_def):
    # Atlas values. The rows from top to bottom are D+W-, D-W+, D*+W-, D*-W+.
    atlas_vals = np.array([[15.04, 15.34, 13.78, 5.13, 0.93],
                        [14.61, 15.12, 13.07, 4.84, 0.82],
                        [14.50, 15.88, 14.19, 5.42, 1.07],
                        [14.26, 15.60, 14.08, 5.11, 0.99]])

    atlas_vals_up_err = np.array([[0.19+0.76, 0.14+0.78, 0.12+0.92, 0.07+0.34, 0.04+0.09],
                                    [0.19+0.73, 0.15+0.75, 0.12+0.89, 0.07+0.31, 0.04+0.08],
                                    [0.26+0.85, 0.19+0.73, 0.16+0.68, 0.10+0.31, 0.05+0.10],
                                    [0.27+0.82, 0.20+0.74, 0.17+0.68, 0.10+0.30, 0.06+0.09]])

    atlas_vals_down_err = np.array([[0.19+0.72, 0.14+0.75, 0.12+0.85, 0.07+0.31, 0.04+0.08],
                                    [0.19+0.69, 0.15+0.72, 0.12+0.82, 0.07+0.29, 0.04+0.07],
                                    [0.26+0.79, 0.19+0.69, 0.16+0.64, 0.10+0.29, 0.05+0.09],
                                    [0.27+0.76, 0.20+0.70, 0.17+0.64, 0.10+0.28, 0.06+0.08]])

    font_size = 17
    axis_label_font_size = 19
    axis_font_size = 14
    legend_fontsize = 16

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1], 'hspace': 0}, figsize=(6, 6))

    pTD_bins = np.array([8, 12, 20, 40, 80, 150])

    pTD_data_min = 8.
    pTD_data_bin_width = 0.5

    bin_widths = np.diff(pTD_bins)

    places_inside_bins = np.zeros((2, 5))
    places_inside_bins_right = np.zeros((2, 5))
    places_inside_bins_left = np.zeros((2, 5))

    places_inside_bins[0, :] = pTD_bins[:-1]**(2 / 3) * pTD_bins[1:]**(1 / 3)
    places_inside_bins[1, :] = pTD_bins[:-1]**(1 / 3) * pTD_bins[1:]**(2 / 3)

    bar_width_over_bin_width = 1. / 9.
    bar_width = bar_width_over_bin_width * bin_widths
    shifts = np.zeros(5)

    scalings = [0.93, 1.1]

    for i in range(5):
        shifts[i] = ((pTD_bins[i + 1] * 1.) / (pTD_bins[i] * 1.))**(bar_width_over_bin_width / 2.)

    for i in range(2):
        places_inside_bins_left[i, :] = places_inside_bins[i, :] / shifts
        places_inside_bins_right[i, :] = places_inside_bins[i, :] * shifts

    bin_midpoints = np.zeros(len(pTD_bins) - 1)
    bin_midpoints_linear_scale = np.zeros(len(pTD_bins) - 1)

    for i in range(len(pTD_bins) - 1):
        bin_midpoints[i] = np.sqrt(pTD_bins[i] * pTD_bins[i + 1])
        bin_midpoints_linear_scale[i] = (pTD_bins[i + 1] + pTD_bins[i]) / 2
    
    xmin = np.sqrt(pTD_bins[0:-1] * bin_midpoints)
    xmax = np.sqrt(pTD_bins[1:] * bin_midpoints)

    QCD_order = 'NLO'
    
    FF_scale_choices = ['frag_main_scale', 'frag_meson_pT_scale']
    colors = ['red', 'blue']
    labels = [r'$\mu_\text{frag} = M_W$', r'$\mu_\text{frag} = p_T(D)$']

    scale_MW_vals = np.zeros(5)
    scale_pTD_vals = np.zeros(5)


    for i in range(len(FF_scale_choices)):
        FF_scale_choice = FF_scale_choices[i]

        scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus = compute_general_3D_vals_NLO(PDF_set,
                num_err_members_in_set, process, False, True, True, FF_scale_choice, z_def, fragmentation_set)

        HISTO_central_sigma_vals = np.zeros(5)

        scale_names = ['central', 'dd', 'uu']

        for eta_lept_index in range(5):
            bin_index = 0
            for pTD_index in range(284):
                if (pTD_data_min + (pTD_index + 1 / 2) * pTD_data_bin_width > pTD_bins[bin_index + 1]):
                    if (bin_index == 5):
                        break
                    else:
                        bin_index += 1
                
                HISTO_central_sigma_vals[bin_index] += sum(scales_vals[0][eta_lept_index][pTD_index, :])

        if (i == 0):
            scale_MW_vals = HISTO_central_sigma_vals
        elif (i == 1):
            scale_pTD_vals = HISTO_central_sigma_vals

        ax1.hlines(HISTO_central_sigma_vals, pTD_bins[:-1], pTD_bins[1:],
                    color=colors[i],
                    label=labels[i], zorder=2)

    ratio = scale_pTD_vals / scale_MW_vals

    ax2.hlines(ratio, pTD_bins[:-1], pTD_bins[1:], zorder=3, color=colors[1])

    plt.xscale('log', base=10)
    ax1.set_ylim(0, 26)
    plt.xlim(8, 150)
    ax2.set_ylim(0.9, 1.1)

    ax2.set_xlabel(r'$p_T(D)$ [GeV]', fontsize=axis_label_font_size)
    ax1.set_ylabel('Cross section [pb]', fontsize=axis_label_font_size)
    ax2.set_ylabel('Ratio', fontsize=axis_label_font_size)

    ax1.set_yticks([10, 20])
    #ax2.set_yticks([0.25, 0.5, 0.75, 1])

    # Configure ticks to appear on all sides
    ax1.tick_params(direction='in', top=True, right=True)

    # Add minor ticks
    ax1.minorticks_on()
    ax1.tick_params(which='both', direction='in', top=True, right=True)

    # Add minor ticks
    ax2.minorticks_on()
    ax2.tick_params(which='both', direction='in', top=True, right=True)

    # Configure ticks to appear on all sides
    ax2.tick_params(direction='in', top=True, right=True)

    ax2.plot([8, 150], [1, 1], color=colors[0], zorder=1)

    legend = ax1.legend(loc='upper right', framealpha=1, fontsize=legend_fontsize)

    info_y_vals_1 = 23
    info_y_vals_2 = 20.5
    info_y_vals_3 = 18
    info_y_vals_4 = 15.5

    info_x_vals_1 = 9
    info_x_vals_2 = 25

    ax1.text(info_x_vals_1, info_y_vals_1, process_text + r'$\quad$OS-SS', fontsize=font_size)
    ax1.text(info_x_vals_1, info_y_vals_2, r'$\sqrt{s} = 13$ TeV', fontsize=font_size)
    ax1.text(info_x_vals_1, info_y_vals_3, PDF_set, fontsize=font_size)
    ax1.text(info_x_vals_1, info_y_vals_4, frag_set_text, fontsize=font_size)

    ax1.tick_params(axis='both', which='major', labelsize=axis_font_size)
    ax2.tick_params(axis='both', which='major', labelsize=axis_font_size)

    plt.xticks(pTD_bins, [f'{tick:.0f}' for tick in pTD_bins])

    for i in range(1, len(pTD_bins) - 1):
        ax1.axvline(pTD_bins[i], color='gray', linewidth=0.5, ymax=0.55, zorder=0)

    for i in range(1, len(pTD_bins) - 1):
        ax2.axvline(pTD_bins[i], color='gray', linewidth=0.5, ymax=1, zorder=0)

    plt.tight_layout()
    plt.savefig(plots_directory + process + '/' + fragmentation_set + "/" + process + '_dymamic_FF_scale_' + PDF_set + '.pdf')
    plt.show()


def pTD_varying_FF_fit(PDF_set, num_err_members, process, plot_errors_flag, FF_fits, theory_labels_here):
    # Atlas values. The rows from top to bottom are D+W-, D-W+, D*+W-, D*-W+.
    atlas_vals = np.array([[15.04, 15.34, 13.78, 5.13, 0.93],
                        [14.61, 15.12, 13.07, 4.84, 0.82],
                        [14.50, 15.88, 14.19, 5.42, 1.07],
                        [14.26, 15.60, 14.08, 5.11, 0.99]])

    atlas_vals_up_err = np.array([[0.19+0.76, 0.14+0.78, 0.12+0.92, 0.07+0.34, 0.04+0.09],
                                    [0.19+0.73, 0.15+0.75, 0.12+0.89, 0.07+0.31, 0.04+0.08],
                                    [0.26+0.85, 0.19+0.73, 0.16+0.68, 0.10+0.31, 0.05+0.10],
                                    [0.27+0.82, 0.20+0.74, 0.17+0.68, 0.10+0.30, 0.06+0.09]])

    atlas_vals_down_err = np.array([[0.19+0.72, 0.14+0.75, 0.12+0.85, 0.07+0.31, 0.04+0.08],
                                    [0.19+0.69, 0.15+0.72, 0.12+0.82, 0.07+0.29, 0.04+0.07],
                                    [0.26+0.79, 0.19+0.69, 0.16+0.64, 0.10+0.29, 0.05+0.09],
                                    [0.27+0.76, 0.20+0.70, 0.17+0.64, 0.10+0.28, 0.06+0.08]])

    font_size = 16
    axis_label_font_size = 17
    axis_font_size = 13
    legend_fontsize = 13

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1], 'hspace': 0}, figsize=(6, 6))

    pTD_bins = np.array([8., 12., 20., 40., 80., 150.])

    pTD_data_min = 8.
    pTD_data_bin_width = 0.5

    bin_widths = np.diff(pTD_bins)

    bin_midpoints = np.zeros(len(pTD_bins) - 1)
    bin_midpoints_linear_scale = np.zeros(len(pTD_bins) - 1)

    for i in range(len(pTD_bins) - 1):
        bin_midpoints[i] = np.sqrt(pTD_bins[i] * pTD_bins[i + 1])
        bin_midpoints_linear_scale[i] = (pTD_bins[i + 1] + pTD_bins[i]) / 2
    
    xmin = np.sqrt(pTD_bins[0:-1] * bin_midpoints)
    xmax = np.sqrt(pTD_bins[1:] * bin_midpoints)

    places_inside_bins = np.zeros((len(FF_fits), 5))
    places_inside_bins_right = np.zeros((len(FF_fits), 5))
    places_inside_bins_left = np.zeros((len(FF_fits), 5))
    places_inside_bins_linear_scale = np.zeros((len(FF_fits), 5))

    for i in range(len(FF_fits)):
        places_inside_bins[i, :] = pTD_bins[:-1] * (pTD_bins[1:] / pTD_bins[:-1])**((i * 1. + 1) / (len(FF_fits) * 1. + 1.))
        places_inside_bins_linear_scale[i, :] = pTD_bins[:-1] + (i * 1. + 1.) / (len(FF_fits) * 1. + 1.) * (pTD_bins[1:] - pTD_bins[:-1])

    bar_width_over_bin_width = 1. / 9.
    bar_width = bar_width_over_bin_width * bin_widths
    shifts = np.zeros(5)

    for i in range(5):
        shifts[i] = ((pTD_bins[i + 1] * 1.) / (pTD_bins[i] * 1.))**(bar_width_over_bin_width / 2.)

    for i in range(len(FF_fits)):
        places_inside_bins_left[i, :] = places_inside_bins[i, :] / shifts
        places_inside_bins_right[i, :] = places_inside_bins[i, :] * shifts

    scalings = [0.9, 1., 1.18]

    for FF_index in range(len(FF_fits)):
        FF_fit = FF_fits[FF_index]

        for QCD_order_index in range(1, 2):
            QCD_order = QCD_orders[QCD_order_index]

            if (plot_errors_flag):
                scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus = compute_general_3D_vals_NLO(PDF_set, num_err_members,
                                                                            process, True, True, True, 'frag_main_scale', z_def, FF_fit)
            else:
                scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus = compute_general_3D_vals_NLO(PDF_set, num_err_members,
                                                                            process, False, False, True, 'frag_main_scale', z_def, FF_fit)

            HISTO_central_sigma_vals = np.zeros(5)
            HISTO_scales_dd_sigma_vals = np.zeros(5)
            HISTO_scales_uu_sigma_vals = np.zeros(5)

            HISTO_central_sigma_MCerrs = np.zeros(5)
            HISTO_scales_dd_sigma_MCerrs = np.zeros(5)
            HISTO_scales_uu_sigma_MCerrs = np.zeros(5)

            HISTO_pdf_err_plus = np.zeros(5)
            HISTO_pdf_err_minus = np.zeros(5)

            scale_names = ['central', 'dd', 'uu']

            for eta_lept_index in range(5):
                bin_index = 0
                for pTD_index in range(284 + 1):
                    if (pTD_data_min + (pTD_index + 1 / 2) * pTD_data_bin_width > pTD_bins[bin_index + 1]):
                        if (bin_index < 4):
                            bin_index += 1
                        else:
                            break

                    HISTO_central_sigma_vals[bin_index] += sum(scales_vals[0][eta_lept_index][pTD_index, :])
                    HISTO_scales_dd_sigma_vals[bin_index] += sum(scales_vals[1][eta_lept_index][pTD_index, :])
                    HISTO_scales_uu_sigma_vals[bin_index] += sum(scales_vals[2][eta_lept_index][pTD_index, :])

                    HISTO_central_sigma_MCerrs[bin_index] += sum(scales_MCerrs[0][eta_lept_index][pTD_index, :])
                    HISTO_scales_dd_sigma_MCerrs[bin_index] += sum(scales_MCerrs[1][eta_lept_index][pTD_index, :])
                    HISTO_scales_uu_sigma_MCerrs[bin_index] += sum(scales_MCerrs[2][eta_lept_index][pTD_index, :])

                    HISTO_pdf_err_plus[bin_index] += sum(pdf_err_plus[eta_lept_index][pTD_index, :])
                    HISTO_pdf_err_minus[bin_index] += sum(pdf_err_minus[eta_lept_index][pTD_index, :])

            if (variation_flag):
                if (QCD_order == "NLO"):
                    ax1.plot(places_inside_bins[FF_index, :], HISTO_central_sigma_vals, marker=markers[FF_index],
                                color=marker_color, markersize=5, linestyle='none',
                                label=theory_labels_here[FF_index], zorder=4)

                    if (plot_errors_flag):
                        ax1.bar(places_inside_bins_left[FF_index, :], HISTO_pdf_err_plus + HISTO_pdf_err_minus, width=bar_width * scalings[FF_index],
                                    bottom=HISTO_central_sigma_vals - HISTO_pdf_err_minus, color=pdf_err_color,
                                    zorder=3)

                        ax1.bar(places_inside_bins_right[FF_index, :], HISTO_scales_uu_sigma_vals - HISTO_scales_dd_sigma_vals, width=bar_width * scalings[FF_index],
                                    bottom=HISTO_central_sigma_vals + HISTO_scales_dd_sigma_vals, color=scale_var_color,
                                    zorder=2)

            else:
                #ax1.hlines(HISTO_central_sigma_vals, pTD_bins[0:-1], pTD_bins[1:], color=line_colors[QCD_order_index], zorder=5, label='LO')
                print('NOT IMPLEMENTED YET')
                exit()
            
            print(HISTO_central_sigma_vals)

        if (QCD_order == 'NLO'):
            ratios = np.zeros(5)
            ratios_ATLAS_var_down = np.zeros(5)
            ratios_ATLAS_var_up = np.zeros(5)
            ratios_theory_scale_var_var_down = np.zeros(5)
            ratios_theory_scale_var_var_up = np.zeros(5)
            ratios_theory_pdf_err_var_down = np.zeros(5)
            ratios_theory_pdf_err_var_up = np.zeros(5)

            for eta_lept_index in range(5):
                ratios[eta_lept_index] = HISTO_central_sigma_vals[eta_lept_index] / atlas_vals[atlas_index][eta_lept_index]

                ratios_ATLAS_var_down[eta_lept_index] = HISTO_central_sigma_vals[eta_lept_index] / \
                                                        (atlas_vals[atlas_index][eta_lept_index] + atlas_vals_up_err[atlas_index][eta_lept_index])
                                                    
                ratios_ATLAS_var_up[eta_lept_index] = HISTO_central_sigma_vals[eta_lept_index] / \
                                                        (atlas_vals[atlas_index][eta_lept_index] - atlas_vals_down_err[atlas_index][eta_lept_index])

                if (plot_errors_flag):
                    ratios_theory_scale_var_var_down[eta_lept_index] = (HISTO_central_sigma_vals[eta_lept_index] + HISTO_scales_dd_sigma_vals[eta_lept_index]) / \
                                                                        HISTO_central_sigma_vals[eta_lept_index]
                                                        
                    ratios_theory_scale_var_var_up[eta_lept_index] = (HISTO_central_sigma_vals[eta_lept_index] + HISTO_scales_uu_sigma_vals[eta_lept_index]) / \
                                                                        HISTO_central_sigma_vals[eta_lept_index]
                    
                    ratios_theory_pdf_err_var_down[eta_lept_index] = (HISTO_central_sigma_vals[eta_lept_index] - HISTO_pdf_err_minus[eta_lept_index]) / \
                                                                        HISTO_central_sigma_vals[eta_lept_index]
                                                        
                    ratios_theory_pdf_err_var_up[eta_lept_index] = (HISTO_central_sigma_vals[eta_lept_index] + HISTO_pdf_err_plus[eta_lept_index]) / \
                                                                        HISTO_central_sigma_vals[eta_lept_index]
                                                        
            ax2.plot(places_inside_bins[FF_index, :], ratios, zorder=4, marker=markers[FF_index],
                        color=marker_color, markersize=5, linestyle='none')

            if (plot_errors_flag):
                scale_var_bar_plot = ax2.bar(places_inside_bins_right[FF_index, :], ratios_theory_scale_var_var_up - ratios_theory_scale_var_var_down,
                        bottom=ratios - 1. + ratios_theory_scale_var_var_down, width=bar_width * scalings[FF_index], color=scale_var_color,
                        linewidth=1, zorder=3)
                
                pdf_err_bar_plot = ax2.bar(places_inside_bins_left[FF_index, :], ratios_theory_pdf_err_var_up - ratios_theory_pdf_err_var_down,
                        bottom=ratios - 1. + ratios_theory_pdf_err_var_down, width=bar_width * scalings[FF_index], color=pdf_err_color,
                        linewidth=1, zorder=3)
    
    ax1.hlines(atlas_vals[atlas_index], pTD_bins[:-1], pTD_bins[1:], color='black', label='ATLAS', zorder=1)
    # Plot error bars for the Atlas values.
    ATLAS_uncertainty = ax1.bar(bin_midpoints_linear_scale, atlas_vals_down_err[atlas_index] + atlas_vals_up_err[atlas_index],
            bottom=atlas_vals[atlas_index] - atlas_vals_down_err[atlas_index], color=atlas_err_color,
            width=bin_widths, zorder=0, linewidth=1.5)
    ax2.bar(bin_midpoints_linear_scale, ratios_ATLAS_var_up - ratios_ATLAS_var_down,
                    bottom=1. - ratios + ratios_ATLAS_var_down, color=atlas_err_color,
                    width=bin_widths, zorder=0, linewidth=1.5)

    plt.xscale('log', base=10)
    ax1.set_ylim(0, 34)
    plt.xlim(8, 150)
    ax2.set_ylim(0.75, 2.05)

    ax2.set_xlabel(r'$p_T(D)$ [GeV]', fontsize=axis_label_font_size)
    ax1.set_ylabel(r'Cross section [pb]', fontsize=axis_label_font_size)
    ax2.set_ylabel(r'$\frac{\mathrm{Theory}}{\mathrm{ATLAS}}$', 
               fontsize=axis_label_font_size * 1.3)

    ax1.set_yticks([5, 10, 15, 20, 25, 30])
    ax2.set_yticks([0.8, 1., 1.2, 1.4, 1.6, 1.8, 2.0])

    # Configure ticks to appear on all sides
    ax1.tick_params(direction='in', top=True, right=True)

    # Add minor ticks
    ax1.minorticks_on()
    ax1.tick_params(which='both', direction='in', top=True, right=True)

    # Add minor ticks
    ax2.minorticks_on()
    ax2.tick_params(which='both', direction='in', top=True, right=True)

    # Configure ticks to appear on all sides
    ax2.tick_params(direction='in', top=True, right=True)

    ax2.plot([8, 150], [1, 1], linewidth=1, color='black', zorder=1)

    legend1 = ax1.legend(loc='upper right', framealpha=1, fontsize=legend_fontsize + 1)
    if (plot_errors_flag):
        legend2 = ax1.legend([ATLAS_uncertainty, pdf_err_bar_plot, scale_var_bar_plot],
                            ["ATLAS error", "PDF error (68\% C.L.)", "Scale variation"], loc='lower left', framealpha=1, fontsize=legend_fontsize - 1)
        ax1.add_artist(legend1)

    info_y_vals_1 = 30.5
    info_y_vals_2 = 27.5
    info_y_vals_3 = 24.5

    info_x_vals_1 = 9
    info_x_vals_2 = 25

    ax1.text(info_x_vals_1, info_y_vals_1, process_text + '  OS-SS', fontsize=font_size)
    ax1.text(info_x_vals_1, info_y_vals_2, r'$\sqrt{s} = 13$ TeV', fontsize=font_size)
    ax1.text(info_x_vals_1, info_y_vals_3, 'CT18ANLO', fontsize=font_size)

    ax1.tick_params(axis='both', which='major', labelsize=axis_font_size)
    ax2.tick_params(axis='both', which='major', labelsize=axis_font_size)

    plt.xticks(pTD_bins, [f'{tick:.0f}' for tick in pTD_bins])

    for i in range(1, len(pTD_bins) - 1):
        ax1.axvline(pTD_bins[i], color='gray', linewidth=0.5, ymax=0.57, zorder=0)

    for i in range(1, len(pTD_bins) - 1):
        ax2.axvline(pTD_bins[i], color='gray', linewidth=0.5, ymax=1, zorder=0)

    plt.tight_layout()

    plt.savefig(plots_directory + process + '/' + PDF_set + '_pT_varying_FF_fit.pdf')
    plt.show()


def eta_lept_plot():
    # Atlas values. The rows from top to bottom are D+W-, D-W+, D*+W-, D*-W+. These values are the values used in the upper plot.
    atlas_vals = np.array([[12.27, 11.57, 10.41, 9.09, 6.85],
                            [11.87, 11.55, 10.09, 8.6, 6.25],
                            [12.18, 11.77, 10.61, 8.85, 7.22],
                            [12.52, 12.14, 10.29, 8.38, 6.55]])
    
    atlas_vals_up_err = np.array([[0.13+0.67, 0.12+0.63, 0.12+0.64, 0.11+0.45, 0.11+0.39],
                                    [0.13+0.65, 0.12+0.61, 0.12+0.61, 0.12+0.43, 0.11+0.37],
                                    [0.18+0.48, 0.17+0.53, 0.17+0.67, 0.16+0.42, 0.16+0.38],
                                    [0.18+0.50, 0.18+0.55, 0.18+0.64, 0.16+0.39, 0.16+0.37]])

    atlas_vals_down_err = np.array([[0.13+0.64, 0.12+0.61, 0.12+0.59, 0.11+0.43, 0.11+0.37],
                                    [0.13+0.62, 0.12+0.60, 0.12+0.57, 0.12+0.41, 0.11+0.35],
                                    [0.18+0.46, 0.17+0.50, 0.17+0.61, 0.16+0.40, 0.16+0.36],
                                    [0.18+0.48, 0.18+0.52, 0.18+0.58, 0.16+0.37, 0.16+0.34]])

    font_size = 17
    axis_label_font_size = 21
    axis_font_size = 15
    legend_fontsize = 13

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1], 'hspace': 0}, figsize=(6, 6))
    eta_lept_bins = np.array([0., 0.5, 1.0, 1.5, 2.0, 2.5])

    xmin = eta_lept_bins[0:-1] + 0.1
    xmax = eta_lept_bins[1:] - 0.1
    bin_midpoints = (xmin + xmax) / 2
    theory_val_places = [bin_midpoints - 0.5/4., bin_midpoints, bin_midpoints + 0.5/4.]

    for PDF_index in range(len(PDF_sets)):
        for QCD_order_index in range(1, 2):
            QCD_order = QCD_orders[QCD_order_index]
            PDF_set = PDF_sets[PDF_index]
            num_err_members_in_set = num_err_members_in_sets[PDF_index]

            scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus = compute_general_3D_vals_NLO(PDF_set, num_err_members_in_set,
                                                                        process, True, True, True, 'frag_main_scale', z_def, fragmentation_set)

            HISTO_central_sigma_vals = np.zeros(5)
            HISTO_scales_dd_sigma_vals = np.zeros(5)
            HISTO_scales_uu_sigma_vals = np.zeros(5)

            HISTO_central_sigma_MCerrs = np.zeros(5)
            HISTO_scales_dd_sigma_MCerrs = np.zeros(5)
            HISTO_scales_uu_sigma_MCerrs = np.zeros(5)

            HISTO_pdf_err_plus = np.zeros(5)
            HISTO_pdf_err_minus = np.zeros(5)

            for eta_lept_index in range(5):
                HISTO_central_sigma_vals[eta_lept_index] = sum(sum(scales_vals[0][eta_lept_index]))
                HISTO_scales_dd_sigma_vals[eta_lept_index] = sum(sum(scales_vals[1][eta_lept_index]))
                HISTO_scales_uu_sigma_vals[eta_lept_index] = sum(sum(scales_vals[2][eta_lept_index]))

                HISTO_central_sigma_MCerrs[eta_lept_index] = sum(sum(scales_MCerrs[0][eta_lept_index]))
                HISTO_scales_dd_sigma_MCerrs[eta_lept_index] = sum(sum(scales_MCerrs[1][eta_lept_index]))
                HISTO_scales_uu_sigma_MCerrs[eta_lept_index] = sum(sum(scales_MCerrs[2][eta_lept_index]))

                HISTO_pdf_err_plus[eta_lept_index] = sum(sum(pdf_err_plus[eta_lept_index]))
                HISTO_pdf_err_minus[eta_lept_index] = sum(sum(pdf_err_minus[eta_lept_index]))

            if (variation_flag and QCD_order == "NLO"):
                bar_width = 0.11
                pdf_err_bar_plot = ax1.bar(theory_val_places[PDF_index] - bar_width / 4., HISTO_pdf_err_plus + HISTO_pdf_err_minus, width=bar_width / 2.,
                            bottom=HISTO_central_sigma_vals - HISTO_pdf_err_minus, color=pdf_err_color, zorder=3)

                scale_var_bar_plot = ax1.bar(theory_val_places[PDF_index] + bar_width / 4., HISTO_scales_uu_sigma_vals - HISTO_scales_dd_sigma_vals,
                            width=bar_width / 2., bottom=HISTO_central_sigma_vals + HISTO_scales_dd_sigma_vals, color=scale_var_color, zorder=2)

                ax1.plot(theory_val_places[PDF_index], HISTO_central_sigma_vals, marker=markers[PDF_index],
                            color=marker_color, markersize=5, linestyle='none', label=theory_labels[PDF_index], zorder=4)
                
                #ax1.bar(theory_val_places[PDF_index], 2. * HISTO_pdf_err, width=0.5,
                #            bottom=HISTO_central_sigma_vals - HISTO_pdf_err, color='gray', zorder=2)
            else:
                print('NOT IMPLEMENTED YET.')
                exit()
                #ax1.hlines(HISTO_central_sigma_vals, theory_val_places[PDF_index] - 0.1, theory_val_places[PDF_index] + 0.1,
                #        color=theory_edge_colors[PDF_index], label=str(PDF_index), zorder=5)

        if (QCD_order == 'NLO'):
            ratios = np.zeros(5)
            ratios_ATLAS_var_down = np.zeros(5)
            ratios_ATLAS_var_up = np.zeros(5)
            ratios_theory_scale_var_var_down = np.zeros(5)
            ratios_theory_scale_var_var_up = np.zeros(5)
            ratios_theory_pdf_err_var_down = np.zeros(5)
            ratios_theory_pdf_err_var_up = np.zeros(5)


            for eta_lept_index in range(5):
                ratios[eta_lept_index] = HISTO_central_sigma_vals[eta_lept_index] / atlas_vals[atlas_index][eta_lept_index]

                ratios_ATLAS_var_down[eta_lept_index] = HISTO_central_sigma_vals[eta_lept_index] / \
                                                    (atlas_vals[atlas_index][eta_lept_index] + atlas_vals_up_err[atlas_index][eta_lept_index])
                ratios_ATLAS_var_up[eta_lept_index] = HISTO_central_sigma_vals[eta_lept_index] / \
                                                    (atlas_vals[atlas_index][eta_lept_index] - atlas_vals_down_err[atlas_index][eta_lept_index])

                ratios_theory_scale_var_var_down[eta_lept_index] = (HISTO_central_sigma_vals[eta_lept_index] + HISTO_scales_dd_sigma_vals[eta_lept_index]) / \
                                                                    HISTO_central_sigma_vals[eta_lept_index]
                ratios_theory_scale_var_var_up[eta_lept_index] = (HISTO_central_sigma_vals[eta_lept_index] + HISTO_scales_uu_sigma_vals[eta_lept_index]) / \
                                                                    HISTO_central_sigma_vals[eta_lept_index]
                
                ratios_theory_pdf_err_var_down[eta_lept_index] = (HISTO_central_sigma_vals[eta_lept_index] - HISTO_pdf_err_minus[eta_lept_index]) / \
                                                                    HISTO_central_sigma_vals[eta_lept_index]
                ratios_theory_pdf_err_var_up[eta_lept_index] = (HISTO_central_sigma_vals[eta_lept_index] + HISTO_pdf_err_plus[eta_lept_index]) / \
                                                                    HISTO_central_sigma_vals[eta_lept_index]
                                                    

            ax2.plot(theory_val_places[PDF_index], ratios, zorder=3, color=marker_color, marker=markers[PDF_index], linestyle='none', markersize=5)

            bar_width = 0.11
            ax2.bar(theory_val_places[PDF_index] + bar_width / 4., ratios_theory_scale_var_var_up - ratios_theory_scale_var_var_down,
                    bottom=ratios - 1. + ratios_theory_scale_var_var_down, width=bar_width / 2., color=scale_var_color,
                    linewidth=1, zorder=2)
            
            ax2.bar(theory_val_places[PDF_index] - bar_width / 4., ratios_theory_pdf_err_var_up - ratios_theory_pdf_err_var_down,
                    bottom=ratios - 1. + ratios_theory_pdf_err_var_down, width=bar_width / 2., color=pdf_err_color,
                    linewidth=1, zorder=2)

    ax1.hlines(atlas_vals[atlas_index], eta_lept_bins[:-1], eta_lept_bins[1:], color='black', label='ATLAS', zorder=1)
    ax2.bar(bin_midpoints, ratios_ATLAS_var_up - ratios_ATLAS_var_down,
                    bottom=1. - ratios + ratios_ATLAS_var_down, width=0.5, zorder=0, color=atlas_err_color)

    # Plot error bars for the Atlas values.
    ATLAS_uncertainty = ax1.bar(bin_midpoints, atlas_vals_down_err[atlas_index] + atlas_vals_up_err[0],
            bottom=atlas_vals[atlas_index] - atlas_vals_down_err[atlas_index], color=atlas_err_color, width=0.5, zorder=0)

    for i in range(1, len(eta_lept_bins) - 1):
        ax1.axvline(eta_lept_bins[i], color='gray', linewidth=0.5, ymax=0.6, zorder=0)

    for i in range(1, len(eta_lept_bins) - 1):
        ax2.axvline(eta_lept_bins[i], color='gray', linewidth=0.5, ymax=1, zorder=0)

    ax1.set_xscale('linear')
    ax1.set_ylim(2, 22)
    ax2.set_ylim(0.6, 1.1)
    plt.xlim(0, 2.5)

    ax2.set_xlabel(r'$|\eta_\mathrm{lepton}|$', fontsize=axis_label_font_size - 2)
    ax1.set_ylabel(r'Cross section [pb]', fontsize=axis_label_font_size)
    ax2.set_ylabel(r'$\frac{\mathrm{Theory}}{\mathrm{ATLAS}}$', 
            fontsize=axis_label_font_size * 1.3)

    ax1.tick_params(axis='both', which='major', labelsize=axis_font_size)
    ax2.tick_params(axis='both', which='major', labelsize=axis_font_size)

    #ax2.set_yticks([0.8, 1.])
    ax1.set_yticks([5, 10, 15, 20])

    info_xval_1 = 0.1
    info_yval_1 = 19.8
    info_yval_2 = 17.8
    info_yval_3 = 15.8

    ax1.text(info_xval_1, info_yval_1, process_text + '  OS-SS', fontsize=font_size)
    ax1.text(info_xval_1, info_yval_2, r'$\sqrt{s} = 13$ TeV', fontsize=font_size)
    ax1.text(info_xval_1, info_yval_3, frag_set_text, fontsize=font_size)

    # Configure ticks to appear on all sides
    ax1.tick_params(direction='in', top=True, right=True)

    # Add minor ticks
    ax1.minorticks_on()
    ax1.tick_params(which='both', direction='in', top=True, right=True)

    # Configure ticks to appear on all sides
    ax2.tick_params(direction='in', top=True, right=True)

    # Add minor ticks
    ax2.minorticks_on()
    ax2.tick_params(which='both', direction='in', top=True, right=True)

    ax2.plot([0, 2.5], [1, 1], linewidth=1, color='black', zorder=1)


    legend1 = ax1.legend(loc='upper right', framealpha=1, fontsize=legend_fontsize)
    legend2 = ax1.legend([ATLAS_uncertainty, pdf_err_bar_plot, scale_var_bar_plot],
                        ["ATLAS error", "PDF error (68\% C.L.)", "Scale variation"], loc='lower left', framealpha=1, fontsize=legend_fontsize)
    ax1.add_artist(legend1)

    plt.tight_layout()
    plt.savefig(plots_directory + process + '/' + fragmentation_set + "/" + process + '_eta_lept.pdf')
    plt.show()


def etaD_plot():
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1], 'hspace': 0}, figsize=(6, 6))

    font_size = 16
    axis_label_font_size = 17
    axis_font_size = 13
    legend_fontsize = 14

    etaD_bins = np.arange(0, 2.201, 0.2)
    print(etaD_bins)
    xmin = etaD_bins[0:-1] + 0.1
    xmax = etaD_bins[1:] - 0.1
    bin_midpoints = (xmin + xmax) / 2
    theory_val_places = [bin_midpoints - 0.2/4., bin_midpoints, bin_midpoints + 0.2/4.]

    for PDF_index in range(len(PDF_sets)):
        for QCD_order_index in range(1, 2):
            QCD_order = QCD_orders[QCD_order_index]
            PDF_set = PDF_sets[PDF_index]
            num_err_members_in_set = num_err_members_in_sets[PDF_index]

            scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus = compute_general_3D_vals_NLO(PDF_set, num_err_members_in_set,
                                                                        process, True, True, True, 'frag_main_scale', z_def, fragmentation_set)

            HISTO_central_sigma_vals = np.zeros(len(etaD_bins) - 1)
            HISTO_scales_dd_sigma_vals = np.zeros(len(etaD_bins) - 1)
            HISTO_scales_uu_sigma_vals = np.zeros(len(etaD_bins) - 1)

            HISTO_central_sigma_MCerrs = np.zeros(len(etaD_bins) - 1)
            HISTO_scales_dd_sigma_MCerrs = np.zeros(len(etaD_bins) - 1)
            HISTO_scales_uu_sigma_MCerrs = np.zeros(len(etaD_bins) - 1)

            HISTO_pdf_err_plus = np.zeros(len(etaD_bins) - 1)
            HISTO_pdf_err_minus = np.zeros(len(etaD_bins) - 1)

            HISTO_total_variation_scales_dd = np.zeros(len(etaD_bins) - 1)
            HISTO_total_variation_scales_uu = np.zeros(len(etaD_bins) - 1)

            HISTO_ratio_up = np.zeros(len(etaD_bins) - 1)
            HISTO_ratio_down = np.zeros(len(etaD_bins) - 1)

            for eta_lept_index in range(5):
                for etaD_index in range(len(etaD_bins) - 1):
                    HISTO_central_sigma_vals[etaD_index] += sum(scales_vals[0][eta_lept_index][:, etaD_index])
                    HISTO_scales_dd_sigma_vals[etaD_index] += sum(scales_vals[1][eta_lept_index][:, etaD_index])
                    HISTO_scales_uu_sigma_vals[etaD_index] += sum(scales_vals[2][eta_lept_index][:, etaD_index])

                    HISTO_central_sigma_MCerrs[etaD_index] += sum(scales_MCerrs[0][eta_lept_index][:, etaD_index])
                    HISTO_scales_dd_sigma_MCerrs[etaD_index] += sum(scales_MCerrs[1][eta_lept_index][:, etaD_index])
                    HISTO_scales_uu_sigma_MCerrs[etaD_index] += sum(scales_MCerrs[2][eta_lept_index][:, etaD_index])

                    HISTO_pdf_err_plus[etaD_index] += sum(pdf_err_plus[eta_lept_index][:, etaD_index])
                    HISTO_pdf_err_minus[etaD_index] += sum(pdf_err_minus[eta_lept_index][:, etaD_index])

            HISTO_ratio_up_pdf_err = (HISTO_central_sigma_vals + HISTO_pdf_err_plus) / HISTO_central_sigma_vals
            HISTO_ratio_down_pdf_err = (HISTO_central_sigma_vals - HISTO_pdf_err_minus) / HISTO_central_sigma_vals

            HISTO_ratio_up_scale_var = (HISTO_central_sigma_vals + HISTO_scales_dd_sigma_vals) / HISTO_central_sigma_vals
            HISTO_ratio_down_scale_var = (HISTO_central_sigma_vals + HISTO_scales_uu_sigma_vals) / HISTO_central_sigma_vals

            if (variation_flag and QCD_order == "NLO"):
                #print(HISTO_central_sigma_vals)
                bar_width = 0.04
                print(len(HISTO_pdf_err_plus))
                print(len(HISTO_pdf_err_minus))
                print(len(theory_val_places[PDF_index]))
                pdf_err_bar_plot = ax1.bar(theory_val_places[PDF_index] - bar_width / 4., HISTO_pdf_err_plus + HISTO_pdf_err_minus, width=bar_width / 2.,
                            bottom=HISTO_central_sigma_vals - HISTO_pdf_err_minus, color=pdf_err_color, zorder=2)
                scale_var_bar_plot = ax1.bar(theory_val_places[PDF_index] + bar_width / 4., HISTO_scales_uu_sigma_vals - HISTO_scales_dd_sigma_vals, width=bar_width / 2.,
                            bottom=HISTO_central_sigma_vals + HISTO_scales_dd_sigma_vals, color=scale_var_color, zorder=2)
                
                ax1.plot(theory_val_places[PDF_index], HISTO_central_sigma_vals, marker=markers[PDF_index],
                            color=marker_color, markersize=5, linestyle='none',
                            label=theory_labels[PDF_index], zorder=4)
                print(HISTO_ratio_up_scale_var)
                ax2.bar(theory_val_places[PDF_index] - bar_width / 4., HISTO_ratio_up_pdf_err - HISTO_ratio_down_pdf_err,
                            width=bar_width / 2., bottom=HISTO_ratio_down_pdf_err, color=pdf_err_color, zorder=2)
                ax2.bar(theory_val_places[PDF_index] + bar_width / 4., HISTO_ratio_up_scale_var - HISTO_ratio_down_scale_var,
                            width=bar_width / 2., bottom=HISTO_ratio_down_scale_var, color=scale_var_color, zorder=2)

    plt.xlabel(r'$|\eta_D|$', fontsize=axis_label_font_size)
    ax1.set_ylabel('Cross section [pb]', fontsize=axis_label_font_size)
    ax2.set_ylabel(r'$\frac{\text{Variation}}{\text{Central}}$', fontsize=axis_label_font_size * 1.3)

    info_xval_1 = 0.1
    info_xval_2 = 0.8
    info_yval_1 = 9.7
    info_yval_2 = 8.9
    info_yval_3 = 8.1

    ax1.tick_params(axis='both', which='major', labelsize=axis_font_size)
    ax1.text(info_xval_1, info_yval_1, process_text + '  OS-SS', fontsize=font_size)
    ax1.text(info_xval_1, info_yval_2, r'$\sqrt{s} = 13$ TeV', fontsize=font_size)
    ax1.text(info_xval_1, info_yval_3, frag_set_text, fontsize=font_size)

    plt.xlim(0, 2.2)
    ax1.set_ylim(2.3, 10.5)
    ax2.set_ylim(0.85, 1.15)

    ax2.plot([-1, 5], [1, 1], color='black', zorder=1)

    plt.xticks(etaD_bins)

    for i in range(1, len(etaD_bins) - 1):
        ax1.axvline(etaD_bins[i], color='gray', linewidth=0.5, ymax=0.6, zorder=0)
        ax2.axvline(etaD_bins[i], color='gray', linewidth=0.5, ymax=1, zorder=0)

    ax1.tick_params(direction='in', top=True, right=True)
    ax1.minorticks_on()
    ax1.tick_params(which='both', direction='in', top=True, right=True)
    ax2.minorticks_on()
    ax2.tick_params(which='both', direction='in', top=True, right=True)
    ax2.tick_params(direction='in', top=True, right=True)
    
    legend1 = ax1.legend(loc='center right', framealpha=1, fontsize=legend_fontsize, bbox_to_anchor=(0.99, 0.83))
    legend2 = ax1.legend([pdf_err_bar_plot, scale_var_bar_plot], ["PDF error (90% C.L.)", "Scale variation"], loc='lower left', framealpha=1, fontsize=legend_fontsize)
    ax1.add_artist(legend1)

    plt.savefig(plots_directory + process + '/' + fragmentation_set + "/" + process + '_etaD.pdf')
    plt.tight_layout()
    plt.show()


def z_variation(QCD_order, PDF_set):
    font_size = 17
    axis_label_font_size = 19
    axis_font_size = 14
    legend_fontsize = 13

    fig, ax = plt.subplots(figsize=(6, 5.5))

    vals = np.loadtxt(main_vals_directory + process + '/' + QCD_order + '/' + z_def + '/' + \
                                                    fragmentation_set + '/frag_main_scale/scale_variation/' + PDF_set + '/central/z_variation.txt')

    if (QCD_order == "NLO"):
        vals -= np.loadtxt(main_vals_directory + process + '/subtraction/' + z_def + '/' + \
                                                        fragmentation_set + '/frag_main_scale/scale_variation/' + PDF_set + '/central/z_variation.txt')

    N = len(vals)

    z_bins = np.linspace(0., 1., N + 1)

    z_point = 0.

    
    for i in range(N - 1):
        if (vals[i] / vals[0] > 0.999 and vals[i + 1] / vals[0] < 0.999):
            z_point = z_bins[i]
            print(z_point)
            #plt.plot([z_bins[i], z_bins[i]], [0., 100.], color='black', zorder=3, label='z = ' + str(round(z_point, 3)))
            break

    for i in range(N - 1):
        if (vals[i] / vals[0] > 0.99 and vals[i + 1] / vals[0] < 0.99):
            z_point = z_bins[i]
            print(z_point)
            #plt.plot([z_bins[i], z_bins[i]], [0., 100.], color='purple', zorder=3, label='z = ' + str(round(z_point, 3)))
            break
    

    plt.plot(z_bins[1:], vals, color='blue', zorder=2, label=r'$\sigma_\text{int}(z_\text{min}^\text{cut}$)')

    plt.plot([z_bins[1], 1.], [vals[0], vals[0]], color='red', zorder=1, label=r'$\sigma_\text{int}(z_\text{min}^\text{cut}$ = 0.05)')

    plt.xlabel(r'$z_\mathrm{min}^\text{cut}$', fontsize=axis_label_font_size)
    plt.ylabel('Cross section [pb]', fontsize=axis_label_font_size)

    plt.xlim(0.05, 1.)
    plt.ylim(0., 45.)

    plt.xticks([0.05, 0.2, 0.4, 0.6, 0.8, 1.])

    text_x = 0.1
    text_y1 = 24
    text_y2 = 20
    text_y3 = 16
    text_y4 = 12


    plt.text(text_x, text_y1, r'$W^-D^+$  OS-SS', fontsize=font_size)
    plt.text(text_x, text_y2, r'$\sqrt{s} = 13$ TeV', fontsize=font_size)
    plt.text(text_x, text_y3, 'CT18ANLO', fontsize=font_size)
    plt.text(text_x, text_y4, 'KKKS08 OPAL', fontsize=font_size)

    ax.tick_params(axis='both', which='major', labelsize=axis_font_size)

    plt.tick_params(direction='in', top=True, right=True)
    plt.minorticks_on()
    plt.tick_params(which='both', direction='in', top=True, right=True)

    ax.tick_params(axis='x', pad=8)
    #ax.tick_params(axis='y', pad=8)

    plt.legend(fontsize=15)

    fig.tight_layout()
    plt.savefig(plots_directory + process + '/' + fragmentation_set + "/" + process + '_z.pdf')

    plt.show()


def z_def_difference(process, PDF_set, num_err_members_in_set, PDF_errors_flag):
    fig, ax = plt.subplots(figsize=(6, 6))

    font_size = 16
    axis_label_font_size = 17
    axis_font_size = 13
    legend_fontsize = 14

    z_def_here = ['minus', 'plus']

    vals_minus = np.array([])
    vals_plus = np.array([])

    pTD_bins = [8., 12., 20., 40., 80., 150.]
    pTD_data_min = pTD_bins[0]
    pTD_data_bin_width = 0.5

    HISTO_central_sigma_vals = np.zeros(5)
    HISTO_scales_dd_sigma_vals = np.zeros(5)
    HISTO_scales_uu_sigma_vals = np.zeros(5)

    HISTO_central_sigma_MCerrs = np.zeros(5)
    HISTO_scales_dd_sigma_MCerrs = np.zeros(5)
    HISTO_scales_uu_sigma_MCerrs = np.zeros(5)

    HISTO_pdf_err_plus = np.zeros(5)
    HISTO_pdf_err_minus = np.zeros(5)

    for z_def_index in range(2):
        z_def = z_def_here[z_def_index]

        if (PDF_errors_flag):
            scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus = compute_general_3D_vals_NLO(PDF_set, num_err_members_in_set,
                                                                        process, True, True, True, 'frag_main_scale', z_def, fragmentation_set)
        else:
            scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus = compute_general_3D_vals_NLO(PDF_set, num_err_members_in_set,
                                                                        process, True, False, True, 'frag_main_scale', z_def, fragmentation_set)

        central_sigma_vals = np.zeros(5)
        scales_dd_sigma_vals = np.zeros(5)
        scales_uu_sigma_vals = np.zeros(5)

        pdf_err_plus = np.zeros(5)
        pdf_err_minus = np.zeros(5)

        scale_names = ['central', 'dd', 'uu']

        for eta_lept_index in range(5):
            bin_index = 0
            for pTD_index in range(284 + 1):
                if (pTD_data_min + (pTD_index + 1 / 2) * pTD_data_bin_width > pTD_bins[bin_index + 1]):
                    if (bin_index < 4):
                        bin_index += 1
                    else:
                        break

                central_sigma_vals[bin_index] += sum(scales_vals[0][eta_lept_index][pTD_index, :])
                scales_dd_sigma_vals[bin_index] += sum(scales_vals[1][eta_lept_index][pTD_index, :])
                scales_uu_sigma_vals[bin_index] += sum(scales_vals[2][eta_lept_index][pTD_index, :])
        
        if (z_def_index == 0):
            sign = 1.
        else:
            sign = -1.

        central_save = HISTO_central_sigma_vals.copy()
        dd_save = HISTO_scales_dd_sigma_vals.copy()
        uu_save = HISTO_scales_uu_sigma_vals.copy()

        HISTO_central_sigma_vals += central_sigma_vals * sign
        HISTO_scales_dd_sigma_vals += scales_dd_sigma_vals * sign
        HISTO_scales_uu_sigma_vals += scales_uu_sigma_vals * sign

        if (z_def_index == 1):
            HISTO_central_sigma_vals = HISTO_central_sigma_vals / central_save
            HISTO_scales_dd_sigma_vals = HISTO_scales_dd_sigma_vals / dd_save
            HISTO_scales_uu_sigma_vals = HISTO_scales_uu_sigma_vals / uu_save

        HISTO_pdf_err_plus += pdf_err_plus * sign
        HISTO_pdf_err_minus += pdf_err_minus * sign
    
    print(HISTO_central_sigma_vals)

    plt.hlines(HISTO_central_sigma_vals * 100., pTD_bins[0:-1], pTD_bins[1:], zorder=2, color='red', linewidth=2)

    for i in range(1, len(pTD_bins) - 2):
        plt.axvline(pTD_bins[i], color='gray', linewidth=0.5, zorder=0)
    for i in range(len(pTD_bins) - 2, len(pTD_bins) - 1):
        plt.axvline(pTD_bins[i], color='gray', linewidth=0.5, ymax=0.165, zorder=0)

    plt.plot([8., 150.], [0., 0.], color='black', zorder=1)

    text_x = 45
    text_y1 = 1.25
    text_y2 = 1.10
    text_y3 = 0.95
    text_y4 = 0.8

    plt.text(text_x, text_y1, process_text + '  OS-SS', fontsize=font_size)
    plt.text(text_x, text_y2, r'$\sqrt{s} = 13$ TeV', fontsize=font_size)
    plt.text(text_x, text_y3, 'CT18ANLO', fontsize=font_size)
    plt.text(text_x, text_y4, 'KKKS08 OPAL', fontsize=font_size)

    plt.xscale('log', base=10)
    plt.xticks(pTD_bins, [f'{tick:.0f}' for tick in pTD_bins])
    plt.xlim(8., 150.)
    plt.ylim(-0.3, 1.5)

    plt.tick_params(direction='in', top=True, right=True)
    plt.minorticks_on()
    plt.tick_params(which='both', direction='in', top=True, right=True)

    plt.tick_params(axis='both', which='major', labelsize=axis_font_size)

    ax.text(0.01, 1.01, r'$\times 10^{-2}$',
        transform=ax.transAxes,
        fontsize=13, va='bottom', ha='left')

    plt.ylabel(r'$\frac{\sigma(z_-) - \sigma(z_+)}{\sigma(z_-)}$', fontsize=axis_label_font_size * 1.3)
    plt.xlabel(r'$p_{T, D} \ \mathrm{[GeV]}$', fontsize=axis_label_font_size)
    plt.tight_layout()

    plt.savefig(plots_directory + process + '/' + fragmentation_set + "/" + process + '_z_def_difference.pdf')
    plt.show()


def Rcpm(which_cross_sections_included):
    font_size = 16
    axis_label_font_size = 18
    axis_font_size = 14
    legend_fontsize = 14

    fig, ax = plt.subplots(figsize=(6, 6))

    y_vals = [1, 2, 3]

    for PDF_index in range(len(PDF_sets)):
        PDF_set = PDF_sets[PDF_index]
        num_err_members_in_set = num_err_members_in_sets[PDF_index]

        process_here = "W-D+"
        scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus = compute_general_3D_vals_NLO(PDF_set, num_err_members_in_set,
                                                                    process_here, True, False, True, 'frag_main_scale', z_def, fragmentation_set)
        Wm_cross_section = sum(sum(sum(scales_vals[0])))
        Wm_scales_dd = sum(sum(sum(scales_vals[1])))
        Wm_scales_uu = sum(sum(sum(scales_vals[2])))
        Wm_MCerr = sum(sum(sum(scales_MCerrs[0])))

        process_here = "W-Dstar+"
        scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus = compute_general_3D_vals_NLO(PDF_set, num_err_members_in_set,
                                                                    process_here, True, False, True, 'frag_main_scale', z_def, fragmentation_set)
        Wm_star_cross_section = sum(sum(sum(scales_vals[0])))
        Wm_star_scales_dd = sum(sum(sum(scales_vals[1])))
        Wm_star_scales_uu = sum(sum(sum(scales_vals[2])))
        Wm_star_MCerr = sum(sum(sum(scales_MCerrs[0])))

        process_here = "W+D-"
        scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus = compute_general_3D_vals_NLO(PDF_set, num_err_members_in_set,
                                                                    process_here, True, False, True, 'frag_main_scale', z_def, fragmentation_set)
        Wp_cross_section = sum(sum(sum(scales_vals[0])))
        Wp_scales_dd = sum(sum(sum(scales_vals[1])))
        Wp_scales_uu = sum(sum(sum(scales_vals[2])))
        Wp_MCerr = sum(sum(sum(scales_MCerrs[0])))

        process_here = "W+Dstar-"
        scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus = compute_general_3D_vals_NLO(PDF_set, num_err_members_in_set,
                                                                    process_here, True, False, True, 'frag_main_scale', z_def, fragmentation_set)
        Wp_star_cross_section = sum(sum(sum(scales_vals[0])))
        Wp_star_scales_dd = sum(sum(sum(scales_vals[1])))
        Wp_star_scales_uu = sum(sum(sum(scales_vals[2])))
        Wp_star_MCerr = sum(sum(sum(scales_MCerrs[0])))

        if (which_cross_sections_included == 'both'):
            Rcpm = (Wp_cross_section + Wp_star_cross_section) / (Wm_cross_section + Wm_star_cross_section)

            Rcpm_scale_var_dd = (Wp_cross_section + Wp_scales_dd + Wp_star_cross_section + Wp_star_scales_dd) / \
                                (Wm_cross_section + Wm_scales_dd + Wm_star_cross_section + Wm_star_scales_dd)
            Rcpm_scale_var_uu = (Wp_cross_section + Wp_scales_uu + Wp_star_cross_section + Wp_star_scales_uu) / \
                                (Wm_cross_section + Wm_scales_uu + Wm_star_cross_section + Wm_star_scales_uu)

            Rcpm_scale_var_down = Rcpm - min(Rcpm_scale_var_dd, Rcpm_scale_var_uu)
            Rcpm_scale_var_up = max(Rcpm_scale_var_dd, Rcpm_scale_var_uu) - Rcpm

            Rcpm_MCerr_up = (Wp_cross_section + Wp_MCerr + Wp_star_cross_section + Wp_star_MCerr) / \
                            (Wm_cross_section - Wm_MCerr + Wm_star_cross_section - Wm_star_MCerr) - Rcpm
            Rcpm_MCerr_down = Rcpm - (Wp_cross_section - Wp_MCerr + Wp_star_cross_section - Wp_star_MCerr) / \
                                        (Wm_cross_section + Wm_MCerr + Wm_star_cross_section + Wm_star_MCerr)

            if (PDF_set == 'CT18NLO' or PDF_set == 'CT18ANLO' or PDF_set == 'MSHT20nlo_as118'):
                Rcpm_pdf_err_up, Rcpm_pdf_err_down = compute_Rcpm_pdf_err_HESSIAN(PDF_index, Rcpm, 'both')
            else:
                Rcpm, Rcpm_pdf_err_up, Rcpm_pdf_err_down = compute_Rcpm_pdf_err_MC(PDF_index, Rcpm, 'both')

        elif (which_cross_sections_included == 'D'):
            Rcpm = Wp_cross_section / Wm_cross_section

            Rcpm_scale_var_dd = (Wp_cross_section + Wp_scales_dd) / (Wm_cross_section + Wm_scales_dd)
            Rcpm_scale_var_uu = (Wp_cross_section + Wp_scales_uu) / (Wm_cross_section + Wm_scales_uu)

            Rcpm_scale_var_down = Rcpm - min(Rcpm_scale_var_dd, Rcpm_scale_var_uu)
            Rcpm_scale_var_up = max(Rcpm_scale_var_dd, Rcpm_scale_var_uu) - Rcpm

            print(PDF_set)
            print(Rcpm_scale_var_down)
            print(Rcpm_scale_var_up)

            Rcpm_MCerr_up = (Wp_cross_section + Wp_MCerr) / (Wm_cross_section - Wm_MCerr) - Rcpm
            Rcpm_MCerr_down = Rcpm - (Wp_cross_section - Wp_MCerr) / (Wm_cross_section + Wm_MCerr)

            if (PDF_set == 'CT18NLO' or PDF_set == 'CT18ANLO' or PDF_set == 'MSHT20nlo_as118'):
                Rcpm_pdf_err_up, Rcpm_pdf_err_down = compute_Rcpm_pdf_err_HESSIAN(PDF_index, Rcpm, 'D')
            else:
                Rcpm, Rcpm_pdf_err_up, Rcpm_pdf_err_down = compute_Rcpm_pdf_err_MC(PDF_index, Rcpm, 'D')

        else:
            Rcpm = Wp_star_cross_section / Wm_star_cross_section

            Rcpm_scale_var_dd = (Wp_star_cross_section + Wp_star_scales_dd) / (Wm_star_cross_section + Wm_star_scales_dd)
            Rcpm_scale_var_uu = (Wp_star_cross_section + Wp_star_scales_uu) / (Wm_star_cross_section + Wm_star_scales_uu)

            Rcpm_scale_var_down = Rcpm - min(Rcpm_scale_var_dd, Rcpm_scale_var_uu)
            Rcpm_scale_var_up = max(Rcpm_scale_var_dd, Rcpm_scale_var_uu) - Rcpm

            Rcpm_MCerr_up = (Wp_star_cross_section + Wp_star_MCerr) / (Wm_star_cross_section - Wm_star_MCerr) - Rcpm
            Rcpm_MCerr_down = Rcpm - (Wp_star_cross_section - Wp_star_MCerr) / (Wm_star_cross_section + Wm_star_MCerr)

            if (PDF_set == 'CT18NLO' or PDF_set == 'CT18ANLO' or PDF_set == 'MSHT20nlo_as118'):
                Rcpm_pdf_err_up, Rcpm_pdf_err_down = compute_Rcpm_pdf_err_HESSIAN(PDF_index, Rcpm, 'Dstar')
            else:
                Rcpm, Rcpm_pdf_err_up, Rcpm_pdf_err_down = compute_Rcpm_pdf_err_MC(PDF_index, Rcpm, 'Dstar')

        Rcpm_error_up = np.sqrt(Rcpm_scale_var_up**2 + Rcpm_MCerr_up**2 + Rcpm_pdf_err_up**2)
        Rcpm_error_down = np.sqrt(Rcpm_scale_var_down**2 + Rcpm_MCerr_down**2 + Rcpm_pdf_err_down**2)

        print('Rcpm (' + which_cross_sections_included + ') with ' + PDF_set + ': ' + str(round(Rcpm, 5)) + \
                '(+' + str(round(Rcpm_error_up, 5)) + '-' + str(round(Rcpm_error_down, 5)) + ').')
    
        ax.plot(Rcpm, 4 - y_vals[PDF_index], marker=markers[PDF_index], color=marker_color,
                linestyle='none', label=theory_labels[PDF_index], zorder=6)

        pdf_err = patches.Rectangle((Rcpm - Rcpm_pdf_err_down, 3 - PDF_index - 0.2),
                    Rcpm_pdf_err_down + Rcpm_pdf_err_up, 0.4, facecolor=pdf_err_color, zorder=5)
        ax.add_patch(pdf_err)

        #scale_var = patches.Rectangle((Rcpm - Rcpm_scale_var_down, 3 - PDF_index - 0.2),
        #            Rcpm_scale_var_down + Rcpm_scale_var_up, 0.4, facecolor='green', zorder=7)
        #ax.add_patch(scale_var)

        total_err = patches.Rectangle((Rcpm - Rcpm_error_down, 3 - PDF_index  - 0.2),
                                        Rcpm_error_down + Rcpm_error_up, 0.4, facecolor=scale_var_color, zorder=4)
        ax.add_patch(total_err)

    atlas_val = 0.971

    atlas_syst_up = 0.011
    atlas_syst_down = 0.011

    atlas_stat_up = 0.006
    atlas_stat_down = 0.006

    plt.plot([atlas_val, atlas_val], [0, 4], color='black', zorder=3)
    plt.plot([-100, 100], [4, 4], color='black', zorder=3)

    atlas_stat_err = patches.Rectangle((atlas_val - atlas_stat_down, -1), atlas_stat_up + atlas_stat_down, 5, facecolor="darkgray", alpha=1, zorder=2)
    ax.add_patch(atlas_stat_err)
    atlas_total_err = patches.Rectangle((atlas_val - np.sqrt(atlas_syst_down**2 + atlas_stat_down**2), -1),
                                        np.sqrt(atlas_syst_up**2 + atlas_stat_up**2) + np.sqrt(atlas_syst_down**2 + atlas_stat_down**2),
                                        5, facecolor="lightgray", alpha=1, zorder=1)
    ax.add_patch(atlas_total_err)

    plt.xlim(0.885, 1)
    plt.ylim(0, 6.7)

    if (which_cross_sections_included == 'both'):
        ax.set_xlabel(r'$R_c^\pm(D^\pm, D^{*\pm})$', fontsize=axis_label_font_size)
    elif (which_cross_sections_included == 'D'):
        ax.set_xlabel(r'$R_c^\pm(D^\pm)$', fontsize=axis_label_font_size)
    else:
        ax.set_xlabel(r'$R_c^\pm(D^{*\pm})$', fontsize=axis_label_font_size)

    info_xval_1 = 0.89
    info_yval_1 = 6.2
    info_yval_2 = 5.8
    info_yval_3 = 5.4

    """
    if (which_cross_sections_included == 'both'):
        ax.text(info_xval_1, info_yval_1, r'$D$, $D^*$', fontsize=font_size)
    elif (which_cross_sections_included == 'D'):
        ax.text(info_xval_1, info_yval_1, r'$D$', fontsize=font_size)
    else:
        ax.text(info_xval_1, info_yval_1, r'$D^*$', fontsize=font_size)
    """

    plt.tick_params(direction='in', top=True, right=True)
    plt.minorticks_on()
    plt.tick_params(which='both', direction='in', top=True, right=True)

    ax.text(info_xval_1, info_yval_1, r'$\sqrt{s} = 13$ TeV', fontsize=font_size)
    ax.text(info_xval_1, info_yval_2, frag_set_text, fontsize=font_size)
    ax.text(info_xval_1, info_yval_3, 'OS-SS', fontsize=font_size)

    legend1 = ax.legend(fontsize=legend_fontsize, bbox_to_anchor=(0.99, 0.88), loc='center right')
    legend2 = ax.legend([pdf_err, total_err], ["PDF error (68\% C.L.)", "Total theory error"], loc='center right',
                        bbox_to_anchor=(0.99, 0.69), framealpha=1, fontsize=legend_fontsize)
    legend3 = ax.legend([atlas_stat_err, atlas_total_err], ['ATLAS stat. error', 'ATLAS tot. error'],
                        loc='center left', bbox_to_anchor=(0.01, 0.69), framealpha=1, fontsize=legend_fontsize)
    ax.add_artist(legend1)
    ax.add_artist(legend2)

    ax.tick_params(axis='both', which='major', labelsize=axis_font_size)
    ax.set_yticklabels([])

    plt.tight_layout()

    plt.savefig(plots_directory + 'Rcpm/Rcpm_' + which_cross_sections_included + '.pdf')
    plt.show()


def Rcpm_pp_pPb():
    font_size = 16
    axis_label_font_size = 17
    axis_font_size = 13
    legend_fontsize = 14

    fig, ax = plt.subplots(figsize=(6, 6))

    PDF_sets_here = ['NNPDF30_nlo_as_01180_pp', 'NNPDF30_nlo_as_01180_pPb']

    Wm_cross_section = Wm_scales_dd = Wm_scales_uu = Wm_MCerr = Wm_star_cross_section = Wm_star_scales_dd = Wm_star_scales_uu = Wm_star_MCerr = \
         Wp_cross_section = Wp_scales_dd = Wp_scales_uu = Wp_MCerr = Wp_star_cross_section = Wp_star_scales_dd = Wp_star_scales_uu = Wp_star_MCerr = 0.0

    sign = 1.
    for PDF_index in range(len(PDF_sets_here)):
        sign -= PDF_index * 2.

        PDF_set = PDF_sets_here[PDF_index]
        num_err_members_in_set = num_err_members_in_sets[PDF_index]

        process_here = "W-D+"
        scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus = compute_general_3D_vals_NLO(PDF_set, num_err_members_in_set,
                                                                    process_here, False, False, True, 'frag_main_scale', z_def, fragmentation_set)
        Wm_cross_section += sum(sum(sum(scales_vals[0]))) / 2. * sign

        process_here = "W-Dstar+"
        scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus = compute_general_3D_vals_NLO(PDF_set, num_err_members_in_set,
                                                                    process_here, False, False, True, 'frag_main_scale', z_def, fragmentation_set)
        Wm_star_cross_section += sum(sum(sum(scales_vals[0]))) / 2. * sign

        process_here = "W+D-"
        scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus = compute_general_3D_vals_NLO(PDF_set, num_err_members_in_set,
                                                                    process_here, False, False, True, 'frag_main_scale', z_def, fragmentation_set)
        Wp_cross_section += sum(sum(sum(scales_vals[0]))) / 2. * sign

        process_here = "W+Dstar-"
        scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus = compute_general_3D_vals_NLO(PDF_set, num_err_members_in_set,
                                                                    process_here, False, False, True, 'frag_main_scale', z_def, fragmentation_set)
        Wp_star_cross_section += sum(sum(sum(scales_vals[0]))) / 2. * sign

        #print(Wm_cross_section)
        #print(Wm_star_cross_section)
        #print(Wp_cross_section)
        #print(Wp_star_cross_section)
        #print()

        Rcpm = (Wp_cross_section + Wp_star_cross_section) / (Wm_cross_section + Wm_star_cross_section)
        print(PDF_set, Rcpm)


def total_cross_section():
    font_size = 16
    axis_label_font_size = 17
    axis_font_size = 13
    legend_fontsize = 14.5

    if (process == 'W-D+'):
        atlas_val = 50.2

        atlas_stat_up = 0.2
        atlas_stat_down = 0.2

        atlas_syst_up = 2.4
        atlas_syst_down = 2.3
    elif (process == 'W+D-'):
        atlas_val = 48.5

        atlas_stat_up = 0.2
        atlas_stat_down = 0.2

        atlas_syst_up = 2.3
        atlas_syst_down = 2.2
    elif (process == 'W-Dstar+'):
        atlas_val = 51.1

        atlas_stat_up = 0.4
        atlas_stat_down = 0.4

        atlas_syst_up = 1.9
        atlas_syst_down = 1.8
    else:
        atlas_val = 50.0

        atlas_stat_up = 0.4
        atlas_stat_down = 0.4

        atlas_syst_up = 1.9
        atlas_syst_down = 1.8

    fig, ax = plt.subplots(figsize=(6, 6))
    y_vals = [1, 2, 3]
    for PDF_set_index in range(len(PDF_sets)):
        PDF_set = PDF_sets[PDF_set_index]
        num_err_members_in_set = num_err_members_in_sets[PDF_set_index]

        scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus = compute_general_3D_vals_NLO(PDF_set, num_err_members_in_set,
                                                                    process, True, True, True, 'frag_main_scale', z_def, fragmentation_set)

        val = sum(sum(sum(scales_vals[0])))

        print(process + " - " + PDF_set + ": " + str(val))

        scales_dd = sum(sum(sum(scales_vals[1])))
        scales_uu = sum(sum(sum(scales_vals[2])))

        scales_down = min(scales_dd, scales_uu)
        scales_up = max(scales_dd, scales_uu)

        pdf_err_up = sum(sum(sum(pdf_err_plus)))
        pdf_err_down = sum(sum(sum(pdf_err_minus)))

        plt.plot(val, 4 - y_vals[PDF_set_index], marker=markers[PDF_set_index], color=marker_color,
                linestyle='none', label=theory_labels[PDF_set_index], zorder=5)

        scale_var = patches.Rectangle((val + scales_uu, 3 - PDF_set_index - 0.25), scales_dd - scales_uu, 0.25, facecolor=scale_var_color, zorder=4)
        ax.add_patch(scale_var)

        pdf_err = patches.Rectangle((val - pdf_err_down, 3 - PDF_set_index), pdf_err_up + pdf_err_down, 0.25, facecolor=pdf_err_color, zorder=4)
        ax.add_patch(pdf_err)

    plt.plot([atlas_val, atlas_val], [0, 4], color='black', zorder=3)
    plt.plot([-100, 100], [4, 4], color='black', zorder=7)

    atlas_tot_err = patches.Rectangle((atlas_val - np.sqrt(atlas_syst_down**2 + atlas_stat_down**2), -1),
                                    np.sqrt(atlas_syst_down**2 + atlas_stat_down**2) + np.sqrt(atlas_syst_up**2 + atlas_stat_up**2),
                                    5, facecolor="lightgray", alpha=1, zorder=1)
    atlas_stat_err = patches.Rectangle((atlas_val - np.sqrt(atlas_stat_down**2), -1),
                                    np.sqrt(atlas_stat_down**2) + np.sqrt(atlas_stat_up**2),
                                    5, facecolor="darkgray", alpha=1, zorder=2)
    
    ax.add_patch(atlas_tot_err)
    ax.add_patch(atlas_stat_err)

    plt.xlim(35, 55)
    plt.ylim(0, 7.)

    ax.set_xlabel('Cross section [pb]', fontsize=axis_label_font_size)

    info_xval_1 = 36
    info_xval_2 = 1.
    info_yval_1 = 6.4
    info_yval_2 = 5.9
    info_yval_3 = 5.4

    ax.text(info_xval_1, info_yval_1, process_text + '  OS-SS', fontsize=font_size)
    ax.text(info_xval_1, info_yval_2, r'$\sqrt{s} = 13$ TeV', fontsize=font_size)
    ax.text(info_xval_1, info_yval_3, frag_set_text, fontsize=font_size)

    legend1 = ax.legend(fontsize=legend_fontsize, bbox_to_anchor=(1, 0.696), loc='center right')
    legend2 = ax.legend([pdf_err, scale_var], ["PDF error (68\% C.L.)", "Scale variation"], framealpha=1, fontsize=legend_fontsize, loc='upper right')
    legend3 = ax.legend([atlas_stat_err, atlas_tot_err], ["ATLAS stat. error", "ATLAS tot. error"], bbox_to_anchor=(0., 0.665), framealpha=1, fontsize=legend_fontsize, loc='center left')

    ax.add_artist(legend1)
    ax.add_artist(legend2)

    ax.tick_params(axis='both', which='major', labelsize=axis_font_size)
    ax.set_yticklabels([])

    ax.minorticks_on()
    ax.tick_params(which='both', direction='in', top=True, right=True)
    ax.tick_params(direction='in', top=True, right=True)

    plt.tight_layout()

    plt.savefig(plots_directory + process + '/' + fragmentation_set + "/" + process + '_integrated.pdf')
    plt.show()


def Rcpm_bin_integrated(kinematic_variable, which_cross_sections_included, plot_errors_flag, PDF_sets):
    font_size = 17
    axis_label_font_size = 20
    axis_font_size = 15
    legend_fontsize = 14
    
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1], 'hspace': 0}, figsize=(6, 6))

    if (kinematic_variable == 'pTD'):
        plt.xscale('log')
        bin_edges = np.array([8., 12., 20., 40., 80., 150.])
    else:
        bin_edges = np.array([0., 0.5, 1.0, 1.5, 2.0, 2.5])

    pTD_raw_bin_width = 0.5

    bin_widths = np.diff(bin_edges)

    bin_midpoints_log_scale = np.zeros(len(bin_edges) - 1)
    bin_midpoints_linear_scale = np.zeros(len(bin_edges) - 1)

    for i in range(len(bin_edges) - 1):
        bin_midpoints_log_scale[i] = np.sqrt(bin_edges[i] * bin_edges[i + 1])
        bin_midpoints_linear_scale[i] = (bin_edges[i + 1] + bin_edges[i]) / 2

    places_inside_bins_log_scale = np.zeros((3, 5))

    places_inside_bins_log_scale[0, :] = np.sqrt(bin_edges[:-1] * bin_midpoints_log_scale)
    places_inside_bins_log_scale[1, :] = bin_midpoints_log_scale
    places_inside_bins_log_scale[2, :] = np.sqrt(bin_edges[1:] * bin_midpoints_log_scale)

    bar_widths = np.zeros((3, 5))

    if (kinematic_variable == 'pTD'):
        width_parameter = 1.03
        bar_widths = places_inside_bins_log_scale * width_parameter - places_inside_bins_log_scale / width_parameter
    else:
        for i in range(3):
            for j in range(5):
                bar_widths[i, j] = 0.05

    places_inside_bins_linear_scale = np.zeros((3, 5))

    places_inside_bins_linear_scale[0, :] = (bin_midpoints_linear_scale + bin_edges[:-1]) / 2.
    places_inside_bins_linear_scale[1, :] = bin_midpoints_linear_scale
    places_inside_bins_linear_scale[2, :] = (bin_midpoints_linear_scale + bin_edges[1:]) / 2.

    #--------------------------------------------------------------------------------------------------------------------------------------#
    #                                                   ATLAS VALUES AND COVARIANCE MATRICES                                               #
    #--------------------------------------------------------------------------------------------------------------------------------------#
    #The rows from top to bottom are D+W-, D-W+, D*+W-, D*-W+.

    if (kinematic_variable == 'pTD'):
        atlas_vals = np.array([[15.04, 15.34, 13.78, 5.13, 0.93],
                            [14.61, 15.12, 13.07, 4.84, 0.82],
                            [14.50, 15.88, 14.19, 5.42, 1.07],
                            [14.26, 15.60, 14.08, 5.11, 0.99]])

        atlas_covariance_starless = np.array([[0.0169, 0.0360, 0.0476, 0.0203, 0.00417, 0.0196, 0.0335, 0.0466, 0.0180, 0.00782],
                                            [0.0889, 0.171, 0.215, 0.0951, 0.0190, 0.0916, 0.160, 0.210, 0.0989, 0.0180],
                                            [0.197, 0.597, 0.783, 0.232, 0.0486, 0.196, 0.566, 0.804, 0.210, 0.0466],
                                            [0.305, 0.565, 0.593, 0.181, 0.0341, 0.279, 0.619, 0.566, 0.160, 0.0335],
                                            [0.495, 0.298, 0.201, 0.110, 0.0164, 0.615, 0.279, 0.196, 0.0916, 0.0196],
                                            [0.0142, 0.0370, 0.0493, 0.0207, 0.00853, 0.0164, 0.0341, 0.0486, 0.0190, 0.00417],
                                            [0.105, 0.190, 0.236, 0.119, 0.0207, 0.110, 0.181, 0.232, 0.0951, 0.0203],
                                            [0.203, 0.622, 0.864, 0.236, 0.0493, 0.201, 0.593, 0.783, 0.215, 0.0476],
                                            [0.313, 0.664, 0.622, 0.190, 0.0370, 0.298, 0.565, 0.597, 0.171, 0.0360],
                                            [0.645, 0.313, 0.203, 0.105, 0.0142, 0.495, 0.305, 0.197, 0.0889, 0.0169]])

        atlas_covariance_star = np.array([[0.0262, 0.0289, 0.0265, 0.0126, 0.00452, 0.0262, 0.0302, 0.0253, 0.0112, 0.0115],
                                        [0.00164, 0.0229, 0.167, 0.0765, 0.0148, 0.000730, 0.0301, 0.158, 0.0999, 0.0112],
                                        [0.0606, 0.0991, 0.414, 0.172, 0.0277, 0.0486, 0.115, 0.495, 0.158, 0.0253],
                                        [0.507, 0.483, 0.108, 0.0346, 0.0357, 0.485, 0.612, 0.115, 0.0301, 0.0302],
                                        [0.561, 0.490, 0.0464, 0.00176, 0.0295, 0.727, 0.485, 0.0486, 0.000730, 0.0262],
                                        [0.0297, 0.0339, 0.0315, 0.0149, 0.0144, 0.0295, 0.0357, 0.0277, 0.0148, 0.00452],
                                        [0.00149, 0.0261, 0.177, 0.106, 0.0149, 0.00176, 0.0346, 0.172, 0.0765, 0.0126],
                                        [0.0470, 0.0975, 0.492, 0.177, 0.0315, 0.0464, 0.108, 0.414, 0.167, 0.0265],
                                        [0.508, 0.584, 0.0975, 0.0261, 0.0339, 0.490, 0.483, 0.0991, 0.0229, 0.0289],
                                        [0.818, 0.508, 0.0470, 0.00149, 0.0297, 0.561, 0.507, 0.0606, 0.00164, 0.0262]])

    else:
        atlas_vals = np.array([[12.27, 11.57, 10.41, 9.09, 6.85],
                                [11.87, 11.55, 10.09, 8.6, 6.25],
                                [12.18, 11.77, 10.61, 8.85, 7.22],
                                [12.52, 12.14, 10.29, 8.38, 6.55]])
        
        atlas_covariance_starless = np.array([[0.201, 0.183, 0.184, 0.126, 0.115, 0.193, 0.178, 0.179, 0.121, 0.154],
                                            [0.245, 0.229, 0.219, 0.165, 0.129, 0.237, 0.221, 0.212, 0.206, 0.121],
                                            [0.390, 0.349, 0.367, 0.223, 0.194, 0.376, 0.340, 0.402, 0.212, 0.179],
                                            [0.385, 0.396, 0.352, 0.236, 0.208, 0.375, 0.431, 0.340, 0.221, 0.178],
                                            [0.425, 0.392, 0.388, 0.249, 0.216, 0.463, 0.375, 0.376, 0.237, 0.193],
                                            [0.222, 0.215, 0.206, 0.142, 0.180, 0.216, 0.208, 0.194, 0.129, 0.115],
                                            [0.259, 0.243, 0.233, 0.221, 0.142, 0.249, 0.236, 0.223, 0.165, 0.126],
                                            [0.401, 0.363, 0.435, 0.233, 0.206, 0.388, 0.352, 0.367, 0.219, 0.184],
                                            [0.398, 0.457, 0.363, 0.243, 0.215, 0.392, 0.396, 0.349, 0.229, 0.183],
                                            [0.490, 0.398, 0.401, 0.259, 0.222, 0.425, 0.385, 0.390, 0.245, 0.201]])

        atlas_covariance_star = np.array([[0.0733, 0.129, 0.166, 0.0941, 0.0854, 0.0782, 0.138, 0.156, 0.0875, 0.168],
                                        [0.109, 0.156, 0.174, 0.129, 0.0951, 0.118, 0.164, 0.167, 0.192, 0.0875],
                                        [0.0898, 0.255, 0.359, 0.180, 0.169, 0.105, 0.267, 0.423, 0.167, 0.156],
                                        [0.154, 0.261, 0.281, 0.176, 0.146, 0.170, 0.352, 0.267, 0.164, 0.138],
                                        [0.197, 0.161, 0.107, 0.131, 0.0805, 0.296, 0.170, 0.105, 0.118, 0.0782],
                                        [0.0778, 0.139, 0.179, 0.100, 0.177, 0.0805, 0.146, 0.169, 0.0951, 0.0854],
                                        [0.116, 0.172, 0.187, 0.212, 0.100, 0.131, 0.176, 0.180, 0.129, 0.0941],
                                        [0.0992, 0.268, 0.458, 0.187, 0.179, 0.107, 0.281, 0.359, 0.174, 0.166],
                                        [0.145, 0.323, 0.268, 0.172, 0.139, 0.161, 0.261, 0.255, 0.156, 0.129],
                                        [0.273, 0.145, 0.0992, 0.116, 0.0778, 0.197, 0.154, 0.0898, 0.109, 0.0733]])
    
    #--------------------------------------------------------------------------------------------------------------------------------------#
    #                                                             FILLING HISTOGRAMS                                                       #
    #--------------------------------------------------------------------------------------------------------------------------------------#
    Rcpm_atlas = np.zeros(5)
    Rcpm_atlas_error = np.zeros(5)

    HISTO_Rcpm_central = [np.zeros(5) for _ in range(len(PDF_sets))]
    HISTO_Rcpm_scale_var_up = [np.zeros(5) for _ in range(len(PDF_sets))]
    HISTO_Rcpm_scale_var_down = [np.zeros(5) for _ in range(len(PDF_sets))]
    HISTO_Rcpm_MCerr_up = [np.zeros(5) for _ in range(len(PDF_sets))]
    HISTO_Rcpm_MCerr_down = [np.zeros(5) for _ in range(len(PDF_sets))]
    HISTO_Rcpm_pdf_err_up = [np.zeros(5) for _ in range(len(PDF_sets))]
    HISTO_Rcpm_pdf_err_down = [np.zeros(5) for _ in range(len(PDF_sets))]
    HISTO_Rcpm_error_up = [np.zeros(5) for _ in range(len(PDF_sets))]
    HISTO_Rcpm_error_down = [np.zeros(5) for _ in range(len(PDF_sets))]

    HISTO_ratios_central = [np.zeros(5) for _ in range(len(PDF_sets))]
    HISTO_ratios_error_up = [np.zeros(5) for _ in range(len(PDF_sets))]
    HISTO_ratios_error_down = [np.zeros(5) for _ in range(len(PDF_sets))]
    HISTO_ratios_pdf_error_up = [np.zeros(5) for _ in range(len(PDF_sets))]
    HISTO_ratios_pdf_error_down = [np.zeros(5) for _ in range(len(PDF_sets))]

    for pTD_index in range(5):
        Rcpm_atlas[pTD_index] = (atlas_vals[1][pTD_index] + atlas_vals[3][pTD_index]) / \
                                        (atlas_vals[0][pTD_index] + atlas_vals[2][pTD_index])

        A = atlas_vals[1][pTD_index] + atlas_vals[3][pTD_index]
        B = atlas_vals[0][pTD_index] + atlas_vals[2][pTD_index]

        var_sigma_p = atlas_covariance_starless[4 - pTD_index, 5 + pTD_index]
        var_sigma_m = atlas_covariance_starless[9 - pTD_index, pTD_index]
        var_sigma_sp = atlas_covariance_star[4 - pTD_index, 5 + pTD_index]
        var_sigma_sm = atlas_covariance_star[9 - pTD_index, pTD_index]

        corr_p_m = atlas_covariance_starless[4 - pTD_index, pTD_index]
        corr_sp_sm = atlas_covariance_star[4 - pTD_index, pTD_index]
        
        if (which_cross_sections_included == 'both'):
            Rcpm_atlas[pTD_index] = (atlas_vals[1][pTD_index] + atlas_vals[3][pTD_index]) / \
                                        (atlas_vals[0][pTD_index] + atlas_vals[2][pTD_index])

            Rcpm_atlas_error[pTD_index] = np.sqrt(1. / B**2 * (var_sigma_p + var_sigma_sp) + \
                                                         A**2 / B**4 * (var_sigma_m + var_sigma_sm) - \
                                                         2. * A / B**3 * (corr_p_m + corr_sp_sm))
        elif (which_cross_sections_included == 'D'):
            Rcpm_atlas[pTD_index] = atlas_vals[1][pTD_index] / atlas_vals[0][pTD_index]

            Rcpm_atlas_error[pTD_index] = np.sqrt(1. / B**2 * (var_sigma_p) + \
                                                         A**2 / B**4 * (var_sigma_m) - \
                                                         2. * A / B**3 * (corr_p_m))
        elif (which_cross_sections_included == 'Dstar'):
            Rcpm_atlas[pTD_index] = atlas_vals[3][pTD_index] / atlas_vals[2][pTD_index]

            Rcpm_atlas_error[pTD_index] = np.sqrt(1. / B**2 * (var_sigma_sp) + \
                                                         A**2 / B**4 * (var_sigma_sm) - \
                                                         2. * A / B**3 * (corr_sp_sm))
        else:
            print('ERROR: INVALID VALUE FOR "which_cross_sections_included".')
            exit(1)

    np.savetxt(reweighting_input_directory + 'experimental_values/' + kinematic_variable + '_' + which_cross_sections_included + '.txt', Rcpm_atlas, delimiter=',')
    
    for PDF_index in range(len(PDF_sets)):
        PDF_set = PDF_sets[PDF_index]
        num_err_members_in_set = num_err_members_in_sets[PDF_index]

        processes_here = np.array(['W+D-', 'W+Dstar-', 'W-D+', 'W-Dstar+'])

        HISTO_central_sigma_vals = [np.zeros(5) for _ in range(len(processes_here))]
        HISTO_scales_dd_sigma_vals = [np.zeros(5) for _ in range(len(processes_here))]
        HISTO_scales_uu_sigma_vals = [np.zeros(5) for _ in range(len(processes_here))]

        HISTO_central_sigma_MCerrs = [np.zeros(5) for _ in range(len(processes_here))]

        for process_index in range(len(processes_here)):
            scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus = compute_general_3D_vals_NLO(
                PDF_set, num_err_members_in_set, processes_here[process_index], True, False, True, 'frag_main_scale', z_def, fragmentation_set)

            
            if (kinematic_variable == 'pTD'):
                for eta_lept_index in range(5):
                    bin_index = 0
                    for pTD_index in range(284):
                        if (bin_edges[0] + (pTD_index + 1 / 2) * pTD_raw_bin_width > bin_edges[bin_index + 1]):
                            if (bin_index == 5):
                                break
                            else:
                                bin_index += 1
                        HISTO_central_sigma_vals[process_index][bin_index] += sum(scales_vals[0][eta_lept_index][pTD_index, :])
                        HISTO_scales_dd_sigma_vals[process_index][bin_index] += sum(scales_vals[1][eta_lept_index][pTD_index, :])
                        HISTO_scales_uu_sigma_vals[process_index][bin_index] += sum(scales_vals[2][eta_lept_index][pTD_index, :])
                        HISTO_central_sigma_MCerrs[process_index][bin_index] += sum(scales_MCerrs[0][eta_lept_index][pTD_index, :])
            else:
                for eta_lept_index in range(5):
                    HISTO_central_sigma_vals[process_index][eta_lept_index] = sum(sum(scales_vals[0][eta_lept_index]))
                    HISTO_scales_dd_sigma_vals[process_index][eta_lept_index] = sum(sum(scales_vals[1][eta_lept_index]))
                    HISTO_scales_uu_sigma_vals[process_index][eta_lept_index] = sum(sum(scales_vals[2][eta_lept_index]))
                    HISTO_central_sigma_MCerrs[process_index][eta_lept_index] = sum(sum(scales_MCerrs[0][eta_lept_index]))
        
        # ax1
        if (which_cross_sections_included == 'both'):
            for pTD_index in range(5):
                HISTO_Rcpm_central[PDF_index][pTD_index] = (HISTO_central_sigma_vals[0][pTD_index] + \
                                                        HISTO_central_sigma_vals[1][pTD_index]) / \
                                                        (HISTO_central_sigma_vals[2][pTD_index] + \
                                                        HISTO_central_sigma_vals[3][pTD_index])

                Rcpm_scales_dd = (HISTO_central_sigma_vals[0][pTD_index] + HISTO_scales_dd_sigma_vals[0][pTD_index] + \
                                            HISTO_central_sigma_vals[1][pTD_index] + HISTO_scales_dd_sigma_vals[1][pTD_index]) / \
                                            (HISTO_central_sigma_vals[2][pTD_index] + HISTO_scales_dd_sigma_vals[2][pTD_index] + \
                                            HISTO_central_sigma_vals[3][pTD_index] + HISTO_scales_dd_sigma_vals[3][pTD_index])
                
                Rcpm_scales_uu = (HISTO_central_sigma_vals[0][pTD_index] + HISTO_scales_uu_sigma_vals[0][pTD_index] + \
                                            HISTO_central_sigma_vals[1][pTD_index] + HISTO_scales_uu_sigma_vals[1][pTD_index]) / \
                                            (HISTO_central_sigma_vals[2][pTD_index] + HISTO_scales_uu_sigma_vals[2][pTD_index] + \
                                            HISTO_central_sigma_vals[3][pTD_index] + HISTO_scales_uu_sigma_vals[3][pTD_index])

                HISTO_Rcpm_scale_var_down[PDF_index][pTD_index] = HISTO_Rcpm_central[PDF_index][pTD_index] - \
                                                            min(Rcpm_scales_dd, Rcpm_scales_uu)
                HISTO_Rcpm_scale_var_up[PDF_index][pTD_index] = max(Rcpm_scales_dd, Rcpm_scales_uu) - \
                                                            HISTO_Rcpm_central[PDF_index][pTD_index]
                
                HISTO_Rcpm_MCerr_up[PDF_index][pTD_index] = (HISTO_central_sigma_vals[0][pTD_index] + \
                                                                    HISTO_central_sigma_MCerrs[0][pTD_index] + \
                                    HISTO_central_sigma_vals[1][pTD_index] + HISTO_central_sigma_MCerrs[1][pTD_index]) / \
                                    (HISTO_central_sigma_vals[2][pTD_index] - HISTO_central_sigma_MCerrs[2][pTD_index] + \
                                    HISTO_central_sigma_vals[3][pTD_index] - HISTO_central_sigma_MCerrs[3][pTD_index]) - \
                                    HISTO_Rcpm_central[PDF_index][pTD_index]
                
                HISTO_Rcpm_MCerr_down[PDF_index][pTD_index] = HISTO_Rcpm_central[PDF_index][pTD_index] - \
                                                                    (HISTO_central_sigma_vals[0][pTD_index] - \
                                                                    HISTO_central_sigma_MCerrs[0][pTD_index] + \
                                    HISTO_central_sigma_vals[1][pTD_index] - HISTO_central_sigma_MCerrs[1][pTD_index]) / \
                                    (HISTO_central_sigma_vals[2][pTD_index] + HISTO_central_sigma_MCerrs[2][pTD_index] + \
                                    HISTO_central_sigma_vals[3][pTD_index] + HISTO_central_sigma_MCerrs[3][pTD_index])

            if (plot_errors_flag):
                if (PDF_set == 'CT18NLO' or PDF_set == 'CT18ANLO' or PDF_set == 'MSHT20nlo_as118'):
                    if (kinematic_variable == 'pTD'):
                        HISTO_Rcpm_pdf_err_up[PDF_index], HISTO_Rcpm_pdf_err_down[PDF_index] = compute_Rcpm_pdf_err_pTD_HESSIAN(PDF_index,
                                                                                HISTO_Rcpm_central[PDF_index], 'both')
                    else:
                        HISTO_Rcpm_pdf_err_up[PDF_index], HISTO_Rcpm_pdf_err_down[PDF_index] = compute_Rcpm_pdf_err_eta_lept_HESSIAN(PDF_index,
                                                                                HISTO_Rcpm_central[PDF_index], 'both')
                else:
                    if (kinematic_variable == 'pTD'):
                        HISTO_Rcpm_central[PDF_index], HISTO_Rcpm_pdf_err_up[PDF_index], HISTO_Rcpm_pdf_err_down[PDF_index] = compute_Rcpm_pdf_err_pTD_MC(PDF_index,
                                                                                HISTO_Rcpm_central[PDF_index], 'both')
                    else:
                        HISTO_Rcpm_central[PDF_index], HISTO_Rcpm_pdf_err_up[PDF_index], HISTO_Rcpm_pdf_err_down[PDF_index] = compute_Rcpm_pdf_err_eta_lept_MC(PDF_index,
                                                                                HISTO_Rcpm_central[PDF_index], 'both')

        elif (which_cross_sections_included == 'D'):
            for pTD_index in range(5):
                HISTO_Rcpm_central[PDF_index][pTD_index] = HISTO_central_sigma_vals[0][pTD_index] / \
                                                        HISTO_central_sigma_vals[2][pTD_index]

                Rcpm_scales_dd = (HISTO_central_sigma_vals[0][pTD_index] + HISTO_scales_dd_sigma_vals[0][pTD_index]) / \
                                            (HISTO_central_sigma_vals[2][pTD_index] + HISTO_scales_dd_sigma_vals[2][pTD_index])
                
                Rcpm_scales_uu = (HISTO_central_sigma_vals[0][pTD_index] + HISTO_scales_uu_sigma_vals[0][pTD_index]) / \
                                            (HISTO_central_sigma_vals[2][pTD_index] + HISTO_scales_uu_sigma_vals[2][pTD_index])
                
                HISTO_Rcpm_scale_var_down[PDF_index][pTD_index] = HISTO_Rcpm_central[PDF_index][pTD_index] - \
                                                            min(Rcpm_scales_dd, Rcpm_scales_uu)
                HISTO_Rcpm_scale_var_up[PDF_index][pTD_index] = max(Rcpm_scales_dd, Rcpm_scales_uu) - \
                                                            HISTO_Rcpm_central[PDF_index][pTD_index]
                
                HISTO_Rcpm_MCerr_up[PDF_index][pTD_index] = (HISTO_central_sigma_vals[0][pTD_index] + \
                                                                    HISTO_central_sigma_MCerrs[0][pTD_index]) / \
                                    (HISTO_central_sigma_vals[2][pTD_index] - HISTO_central_sigma_MCerrs[2][pTD_index]) - \
                                    HISTO_Rcpm_central[PDF_index][pTD_index]
                
                HISTO_Rcpm_MCerr_down[PDF_index][pTD_index] = HISTO_Rcpm_central[PDF_index][pTD_index] - \
                                                                    (HISTO_central_sigma_vals[0][pTD_index] - \
                                                                    HISTO_central_sigma_MCerrs[0][pTD_index]) / \
                                    (HISTO_central_sigma_vals[2][pTD_index] + HISTO_central_sigma_MCerrs[2][pTD_index])

            if (plot_errors_flag):
                if (PDF_set == 'CT18NLO' or PDF_set == 'CT18ANLO' or PDF_set == 'MSHT20nlo_as118'):
                    if (kinematic_variable == 'pTD'):
                        HISTO_Rcpm_pdf_err_up[PDF_index], HISTO_Rcpm_pdf_err_down[PDF_index] = compute_Rcpm_pdf_err_pTD_HESSIAN(PDF_index,
                                                                                HISTO_Rcpm_central[PDF_index], 'D')
                    else:
                        HISTO_Rcpm_pdf_err_up[PDF_index], HISTO_Rcpm_pdf_err_down[PDF_index] = compute_Rcpm_pdf_err_eta_lept_HESSIAN(PDF_index,
                                                                                HISTO_Rcpm_central[PDF_index], 'D')
                else:
                    if (kinematic_variable == 'pTD'):
                        HISTO_Rcpm_central[PDF_index], HISTO_Rcpm_pdf_err_up[PDF_index], HISTO_Rcpm_pdf_err_down[PDF_index] = compute_Rcpm_pdf_err_pTD_MC(PDF_index,
                                                                                HISTO_Rcpm_central[PDF_index], 'D')
                    else:
                        HISTO_Rcpm_central[PDF_index], HISTO_Rcpm_pdf_err_up[PDF_index], HISTO_Rcpm_pdf_err_down[PDF_index] = compute_Rcpm_pdf_err_eta_lept_MC(PDF_index,
                                                                                HISTO_Rcpm_central[PDF_index], 'D')
        else:
            for pTD_index in range(5):
                HISTO_Rcpm_central[PDF_index][pTD_index] = HISTO_central_sigma_vals[1][pTD_index] / \
                                                        HISTO_central_sigma_vals[3][pTD_index]

                Rcpm_scales_dd = (HISTO_central_sigma_vals[1][pTD_index] + HISTO_scales_dd_sigma_vals[1][pTD_index]) / \
                                            (HISTO_central_sigma_vals[3][pTD_index] + HISTO_scales_dd_sigma_vals[3][pTD_index])
                
                Rcpm_scales_uu = (HISTO_central_sigma_vals[1][pTD_index] + HISTO_scales_uu_sigma_vals[1][pTD_index]) / \
                                            (HISTO_central_sigma_vals[3][pTD_index] + HISTO_scales_uu_sigma_vals[3][pTD_index])
                
                HISTO_Rcpm_scale_var_down[PDF_index][pTD_index] = HISTO_Rcpm_central[PDF_index][pTD_index] - \
                                                            min(Rcpm_scales_dd, Rcpm_scales_uu)
                HISTO_Rcpm_scale_var_up[PDF_index][pTD_index] = max(Rcpm_scales_dd, Rcpm_scales_uu) - \
                                                            HISTO_Rcpm_central[PDF_index][pTD_index]
                
                HISTO_Rcpm_MCerr_up[PDF_index][pTD_index] = (HISTO_central_sigma_vals[1][pTD_index] + \
                                                                HISTO_central_sigma_MCerrs[1][pTD_index]) / \
                                    (HISTO_central_sigma_vals[3][pTD_index] - HISTO_central_sigma_MCerrs[3][pTD_index]) - \
                                    HISTO_Rcpm_central[PDF_index][pTD_index]
                
                HISTO_Rcpm_MCerr_down[PDF_index][pTD_index] = HISTO_Rcpm_central[PDF_index][pTD_index] - \
                                                                    (HISTO_central_sigma_vals[1][pTD_index] - \
                                                                    HISTO_central_sigma_MCerrs[1][pTD_index]) / \
                                    (HISTO_central_sigma_vals[3][pTD_index] + HISTO_central_sigma_MCerrs[3][pTD_index])
            
            if (plot_errors_flag):
                if (PDF_set == 'CT18NLO' or PDF_set == 'CT18ANLO' or PDF_set == 'MSHT20nlo_as118'):
                    if (kinematic_variable == 'pTD'):
                        HISTO_Rcpm_pdf_err_up[PDF_index], HISTO_Rcpm_pdf_err_down[PDF_index] = compute_Rcpm_pdf_err_pTD_HESSIAN(PDF_index,
                                                                                HISTO_Rcpm_central[PDF_index], 'Dstar')
                    else:
                        HISTO_Rcpm_pdf_err_up[PDF_index], HISTO_Rcpm_pdf_err_down[PDF_index] = compute_Rcpm_pdf_err_eta_lept_HESSIAN(PDF_index,
                                                                                HISTO_Rcpm_central[PDF_index], 'Dstar')
                else:
                    if (kinematic_variable == 'pTD'):
                        HISTO_Rcpm_central[PDF_index], HISTO_Rcpm_pdf_err_up[PDF_index], HISTO_Rcpm_pdf_err_down[PDF_index] = compute_Rcpm_pdf_err_pTD_MC(PDF_index,
                                                                                HISTO_Rcpm_central[PDF_index], 'Dstar')
                    else:
                        HISTO_Rcpm_central[PDF_index], HISTO_Rcpm_pdf_err_up[PDF_index], HISTO_Rcpm_pdf_err_down[PDF_index] = compute_Rcpm_pdf_err_eta_lept_MC(PDF_index,
                                                                                HISTO_Rcpm_central[PDF_index], 'Dstar')
        
        if (plot_errors_flag):
            for pTD_index in range(5):
                HISTO_Rcpm_error_up[PDF_index][pTD_index] = np.sqrt(HISTO_Rcpm_scale_var_up[PDF_index][pTD_index]**2 + \
                                                                    HISTO_Rcpm_MCerr_up[PDF_index][pTD_index]**2 + \
                                                                    HISTO_Rcpm_pdf_err_up[PDF_index][pTD_index]**2)
                HISTO_Rcpm_error_down[PDF_index][pTD_index] = np.sqrt(HISTO_Rcpm_scale_var_down[PDF_index][pTD_index]**2 + \
                                                                        HISTO_Rcpm_MCerr_down[PDF_index][pTD_index]**2 + \
                                                                        HISTO_Rcpm_pdf_err_down[PDF_index][pTD_index]**2)

        # ax2
        for pTD_index in range(5):
            HISTO_ratios_central[PDF_index][pTD_index] = HISTO_Rcpm_central[PDF_index][pTD_index] / \
                                                                Rcpm_atlas[pTD_index]

            HISTO_ratios_error_up[PDF_index][pTD_index] = (HISTO_Rcpm_central[PDF_index][pTD_index] + \
                                                                HISTO_Rcpm_error_up[PDF_index][pTD_index]) / \
                                                                Rcpm_atlas[pTD_index] - \
                                                                HISTO_ratios_central[PDF_index][pTD_index]

            HISTO_ratios_error_down[PDF_index][pTD_index] = HISTO_ratios_central[PDF_index][pTD_index] - \
                                                                    (HISTO_Rcpm_central[PDF_index][pTD_index] - \
                                                                    HISTO_Rcpm_error_down[PDF_index][pTD_index]) / \
                                                                    Rcpm_atlas[pTD_index]

            HISTO_ratios_pdf_error_up[PDF_index][pTD_index] = (HISTO_Rcpm_central[PDF_index][pTD_index] + \
                                                                HISTO_Rcpm_pdf_err_up[PDF_index][pTD_index]) / \
                                                                Rcpm_atlas[pTD_index] - \
                                                                HISTO_ratios_central[PDF_index][pTD_index]

            HISTO_ratios_pdf_error_down[PDF_index][pTD_index] = HISTO_ratios_central[PDF_index][pTD_index] - \
                                                                    (HISTO_Rcpm_central[PDF_index][pTD_index] - \
                                                                    HISTO_Rcpm_pdf_err_down[PDF_index][pTD_index]) / \
                                                                    Rcpm_atlas[pTD_index]

        #--------------------------------------------------------------------------------------------------------------------------------------#
        #                                                    SAVING THE BEST VALUES FOR REWEIGHTING                                            #
        #--------------------------------------------------------------------------------------------------------------------------------------#
        if (PDF_set == 'CT18NLO' or PDF_set == 'CT18ANLO' or PDF_set == 'MSHT20nlo_as118'):
            np.savetxt(reweighting_input_directory + 'theory_values/HESSIAN/best/' + kinematic_variable + '_' + which_cross_sections_included + '_' + \
                   PDF_sets[PDF_index] + '_best.txt', HISTO_Rcpm_central[PDF_index], delimiter=',')
        else:
            np.savetxt(reweighting_input_directory + 'theory_values/MC/best/' + kinematic_variable + '_' + which_cross_sections_included + '_' + \
                   PDF_sets[PDF_index] + '_best.txt', HISTO_Rcpm_central[PDF_index], delimiter=',')
    #--------------------------------------------------------------------------------------------------------------------------------------#
    #                                                                   PLOTTING                                                           #
    #--------------------------------------------------------------------------------------------------------------------------------------#

    places_inside_bins_here = np.zeros((3, 5))
    if (kinematic_variable == 'pTD'):
        places_inside_bins_here = places_inside_bins_log_scale
    else:
        places_inside_bins_here = places_inside_bins_linear_scale
    
    # ATLAS
    ax1.hlines(Rcpm_atlas, bin_edges[:-1], bin_edges[1:], color='black', zorder=1, label='ATLAS')
    ATLAS_uncertainty = ax1.bar(bin_midpoints_linear_scale, 2 * Rcpm_atlas_error, bottom=Rcpm_atlas - Rcpm_atlas_error,
            width=bin_widths, color=atlas_err_color, zorder=0)
    
    # THEORY
    for PDF_index in range(len(PDF_sets)):
        ax1.plot(places_inside_bins_here[PDF_index, :], HISTO_Rcpm_central[PDF_index], marker=markers[PDF_index],
                    color=marker_color, markersize=5, linestyle='none',
                    label=theory_labels[PDF_index], zorder=4)

        if (plot_errors_flag):
            PDF_uncertainty = ax1.bar(places_inside_bins_here[PDF_index, :], HISTO_Rcpm_pdf_err_down[PDF_index] + HISTO_Rcpm_pdf_err_up[PDF_index],
                                        width=bar_widths[PDF_index], bottom=HISTO_Rcpm_central[PDF_index] - HISTO_Rcpm_pdf_err_down[PDF_index],
                                        color=pdf_err_color, zorder=3)
        
            theory_uncertainty = ax1.bar(places_inside_bins_here[PDF_index, :], HISTO_Rcpm_error_up[PDF_index] + HISTO_Rcpm_error_down[PDF_index],
                    bottom=HISTO_Rcpm_central[PDF_index] - HISTO_Rcpm_error_down[PDF_index],
                    width=bar_widths[PDF_index], color=scale_var_color, zorder=2)

    # THEORY / ATLAS
    for PDF_index in range(len(PDF_sets)):
        ax2.plot(places_inside_bins_here[PDF_index, :], HISTO_ratios_central[PDF_index], marker=markers[PDF_index],
                    color=marker_color, markersize=5, linestyle='none',
                    label=theory_labels[PDF_index], zorder=4)

        if (plot_errors_flag):
            ax2.bar(places_inside_bins_here[PDF_index, :], HISTO_ratios_error_up[PDF_index] + HISTO_ratios_error_down[PDF_index],
                    bottom=HISTO_ratios_central[PDF_index] - HISTO_ratios_error_down[PDF_index],
                    width=bar_widths[PDF_index], color=scale_var_color, zorder=2)
            
            ax2.bar(places_inside_bins_here[PDF_index, :], HISTO_ratios_pdf_error_up[PDF_index] + HISTO_ratios_pdf_error_down[PDF_index],
                    bottom=HISTO_ratios_central[PDF_index] - HISTO_ratios_pdf_error_down[PDF_index],
                    width=bar_widths[PDF_index], color=pdf_err_color, zorder=3)

    ax2.bar(bin_midpoints_linear_scale, 2 * Rcpm_atlas_error / Rcpm_atlas, bottom=1. - Rcpm_atlas_error / Rcpm_atlas,
        width=bin_widths, color=atlas_err_color, zorder=0)

    # DECORATIONS
    for i in range(1, len(bin_edges) - 3):
        if (kinematic_variable == 'pTD'):
            ax1.axvline(bin_edges[i], color='gray', linewidth=0.5, ymax=0.63, zorder=0)
        else:
            ax1.axvline(bin_edges[i], color='gray', linewidth=0.5, ymax=0.67, zorder=0)
    for i in range(len(bin_edges) - 3, len(bin_edges) - 1):
        ax1.axvline(bin_edges[i], color='gray', linewidth=0.5, ymax=0.60, zorder=0)
    for i in range(1, len(bin_edges) - 1):
        ax2.axvline(bin_edges[i], color='gray', linewidth=0.5, zorder=0)

    ax2.plot([-1, 151], [1, 1], color='black', zorder=1)

    #--------------------------------------------------------------------------------------------------------------------------------------#
    #                                                        MAKING THE PLOT LOOK PRETTY :)                                                #
    #--------------------------------------------------------------------------------------------------------------------------------------#
    ax2.set_yticks([0.85, 0.9, 0.95, 1., 1.05, 1.1])

    plt.xlim(bin_edges[0], bin_edges[-1])
    ax1.set_ylim(0.755, 1.2)
    ax2.set_ylim(0.83, 1.14)

    ax1.tick_params(axis='both', which='major', labelsize=axis_font_size)
    ax2.tick_params(axis='both', which='major', labelsize=axis_font_size)

    text_x = 0
    if (kinematic_variable == 'pTD'):
        plt.xlabel(r'$p_T (D)$ [GeV]', fontsize=axis_label_font_size)
        text_x = 9
        plt.xticks(bin_edges, [f'{tick:.0f}' for tick in bin_edges])
    else:
        plt.xlabel(r'$|\eta_\text{lepton}|$', fontsize=axis_label_font_size)
        text_x = 0.1
        plt.xticks(bin_edges, [f'{tick:.1f}' for tick in bin_edges])

    if (which_cross_sections_included == 'both'):
        ax1.set_ylabel(r'$R_c^\pm(D^\pm, D^{*\pm})$', fontsize=axis_label_font_size)
    elif (which_cross_sections_included == 'D'):
        ax1.set_ylabel(r'$R_c^\pm(D^\pm)$', fontsize=axis_label_font_size)
    else:
        ax1.set_ylabel(r'$R_c^\pm(D^{*\pm})$', fontsize=axis_label_font_size)

    ax2.set_ylabel(r'$\frac{\mathrm{Theory}}{\mathrm{ATLAS}}$', fontsize=axis_label_font_size + 6)

    legend1 = ax1.legend(loc='upper right', fontsize=legend_fontsize, framealpha=1)
    ax1.legend([ATLAS_uncertainty, PDF_uncertainty, theory_uncertainty],
                        ['ATLAS error', 'PDF error (68\% C.L.)', 'Total theory error'],
                        loc='lower left', fontsize=legend_fontsize, framealpha=1)

    ax1.add_artist(legend1)

    text_y1 = 1.15
    text_y2 = 1.11
    text_y3 = 1.07

    ax1.text(text_x, text_y1, r'$\sqrt{s} = 13$ TeV', fontsize=font_size)
    if (fragmentation_set == 'KKKS08_opal'):
        ax1.text(text_x, text_y2, 'KKKS08 OPAL', fontsize=font_size)
    elif (fragmentation_set == 'SMSKA19'):
        ax1.text(text_x, text_y2, 'SMSKA19', fontsize=font_size)
    ax1.text(text_x, text_y3, 'OS-SS', fontsize=font_size)
    
    ax1.tick_params(axis='both', which='major', labelsize=axis_font_size)
    ax2.tick_params(axis='both', which='major', labelsize=axis_font_size)
    ax1.tick_params(direction='in', top=True, right=True)
    ax1.minorticks_on()
    ax1.tick_params(which='both', direction='in', top=True, right=True)
    ax2.minorticks_on()
    ax2.tick_params(which='both', direction='in', top=True, right=True)
    ax2.tick_params(direction='in', top=True, right=True)

    plt.tight_layout()

    plt.savefig(plots_directory + 'Rcpm/Rcpm_' + kinematic_variable + '_' + which_cross_sections_included + '.pdf')
    plt.show()


def Rcpm_pTD_varying_FF_fit(PDF_set, which_cross_sections_included):
    font_size = 16
    axis_label_font_size = 17
    axis_font_size = 13
    legend_fontsize = 14

    scalings = [0.9, 1., 1.18]
    
    fig, ax = plt.subplots(figsize=(6, 6))

    pTD_bins = np.array([8, 12, 20, 40, 80, 150])

    pTD_data_min = 8.
    pTD_data_bin_width = 0.5

    bin_widths = np.diff(pTD_bins)

    FF_sets = ['opal', 'global']
    HISTO_Rcpm_central = [np.zeros(5) for _ in range(len(FF_sets))]

    for FF_index in range(len(FF_sets)):
        FF_set_here = FF_sets[FF_index]
        num_err_members_in_set = num_err_members_in_sets[FF_index]

        processes_here = np.array(['W+D-', 'W+Dstar-', 'W-D+', 'W-Dstar+'])

        HISTO_central_sigma_vals = [np.zeros(5) for _ in range(len(processes_here))]

        for process_index in range(len(processes_here)):
            scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus = compute_general_3D_vals_NLO(
                PDF_set, num_err_members_in_set, processes_here[process_index], False, False, True, 'frag_main_scale', z_def, FF_set_here)

            for eta_lept_index in range(5):
                bin_index = 0
                for pTD_index in range(284):
                    if (pTD_data_min + (pTD_index + 1 / 2) * pTD_data_bin_width > pTD_bins[bin_index + 1]):
                        if (bin_index == 5):
                            break
                        else:
                            bin_index += 1

                    HISTO_central_sigma_vals[process_index][bin_index] += sum(scales_vals[0][eta_lept_index][pTD_index, :])

        if (which_cross_sections_included == 'both'):
            for pTD_index in range(5):
                HISTO_Rcpm_central[FF_index][pTD_index] = (HISTO_central_sigma_vals[0][pTD_index] + \
                                                        HISTO_central_sigma_vals[1][pTD_index]) / \
                                                        (HISTO_central_sigma_vals[2][pTD_index] + \
                                                        HISTO_central_sigma_vals[3][pTD_index])
        elif (which_cross_sections_included == 'D'):
            for pTD_index in range(5):
                HISTO_Rcpm_central[FF_index][pTD_index] = HISTO_central_sigma_vals[0][pTD_index] / \
                                                        HISTO_central_sigma_vals[2][pTD_index]
        else:
            for pTD_index in range(5):
                HISTO_Rcpm_central[FF_index][pTD_index] = HISTO_central_sigma_vals[1][pTD_index] / \
                                                        HISTO_central_sigma_vals[3][pTD_index]


    varying_FF_fit_ratio = np.zeros(5)

    print(HISTO_Rcpm_central[0])
    print(HISTO_Rcpm_central[1])

    for pTD_index in range(len(pTD_bins) - 1):
        varying_FF_fit_ratio[pTD_index] = (HISTO_Rcpm_central[1][pTD_index] - HISTO_Rcpm_central[0][pTD_index]) / HISTO_Rcpm_central[1][pTD_index] * 1e3

    #--------------------------------------------------------------------------------------------------------------------------------------#
    #                                                                   PLOTTING                                                           #
    #--------------------------------------------------------------------------------------------------------------------------------------#

    print(varying_FF_fit_ratio)
    ax.hlines(varying_FF_fit_ratio, pTD_bins[:-1], pTD_bins[1:], color='red', zorder=1)

    # DECORATIONS
    for i in range(1, 3):
        ax.axvline(pTD_bins[i], color='gray', linewidth=0.5, ymax=0.5)
    for i in range(3, 5):
        ax.axvline(pTD_bins[i], color='gray', linewidth=0.5)

    ax.text(0.01, 1.01, r'$\times 10^{-3}$',
        transform=ax.transAxes,
        fontsize=14, va='bottom', ha='left')

    #--------------------------------------------------------------------------------------------------------------------------------------#
    #                                                        MAKING THE PLOT LOOK PRETTY :)                                                #
    #--------------------------------------------------------------------------------------------------------------------------------------#

    ax.set_xlim(8, 150)
    ax.set_ylim(0, 0.6)

    plt.xscale('log')

    plt.xlabel(r'$p_T (D)$', fontsize=axis_label_font_size)
    ax.set_ylabel(r'$\frac{R_c^\pm(\text{Global}) - R_c^\pm(\text{Opal})}{R_c^\pm(\text{Global})}$', fontsize=axis_label_font_size * 1.3)

    text_x = 11
    text_y1 = 3.5
    text_y2 = 3.1
    text_y3 = 2.7
    text_y4 = 2.3

    ax.text(text_x, text_y1, r'$R_c^\pm(D^\pm, D^{*\pm})$', fontsize=font_size)
    ax.text(text_x, text_y2, r'$\sqrt{s} = 13$ TeV', fontsize=font_size)
    ax.text(text_x, text_y3, 'CT18ANLO', fontsize=font_size)
    ax.text(text_x, text_y4, 'KKKS08', fontsize=font_size)

    ax.tick_params(direction='in', top=True, right=True)
    ax.minorticks_on()
    ax.tick_params(which='both', direction='in', top=True, right=True)

    ax.tick_params(axis='both', which='major', labelsize=axis_font_size)

    ax.set_yticks([1, 2, 3, 4])
    plt.xticks(pTD_bins, [f'{tick:.0f}' for tick in pTD_bins])

    plt.tight_layout()

    plt.savefig(plots_directory + 'Rcpm/Rcpm_pTD_' + which_cross_sections_included + '_varying_FF_fit.pdf')
    plt.show()


def compute_LO_integrated_cross_section(PDF_set, process, z_def, fragmentation_set):
    s = 0.

    for i in range(5):
        s += sum(sum(np.loadtxt(main_vals_directory + process + '/LO/' + z_def + '/' + \
                    fragmentation_set + '/scale_variation/' + PDF_set + '/central/' + str(i) + '_vals.txt', delimiter=',')))
    return s


def Rcpm_LO(which_cross_sections_included, PDF_set, z_def, fragmentation_set):
    Wp = 0.
    Wm = 0.

    if (which_cross_sections_included == 'D' or which_cross_sections_included == 'both'):
        Wm += compute_LO_integrated_cross_section(PDF_set, 'W-D+', z_def, fragmentation_set)
        Wp += compute_LO_integrated_cross_section(PDF_set, 'W+D-', z_def, fragmentation_set)
    if (which_cross_sections_included == 'Dstar' or which_cross_sections_included == 'both'):
        Wm += compute_LO_integrated_cross_section(PDF_set, 'W-Dstar+', z_def, fragmentation_set)
        Wp += compute_LO_integrated_cross_section(PDF_set, 'W+Dstar-', z_def, fragmentation_set)

    print(Wp / Wm)


def Rcpm_LO_bin_integrated(kinematic_quantity, which_cross_sections_included, PDF_set):
    pTD_bins = [8., 12., 20., 40., 80., 150.]
    pTD_data_bin_width = 0.5
    Wp = np.zeros(5)
    Wm = np.zeros(5)

    for eta_lept_index in range(5):
        bin_index = 0
        for pTD_index in range(284 + 1):
            if (pTD_bins[0] + (pTD_index + 1 / 2) * pTD_data_bin_width > pTD_bins[bin_index + 1]):
                if (bin_index < 4):
                    bin_index += 1
                else:
                    break

            Wm[bin_index] += sum(np.loadtxt(
                        main_vals_directory + 'W-D+/LO/' + z_def + '/' + \
                        fragmentation_set + '/scale_variation/' + PDF_set + '/central/' + \
                        str(eta_lept_index) + '_vals.txt', delimiter=',')[pTD_index, :])

            Wp[bin_index] += sum(np.loadtxt(
                        main_vals_directory + 'W+D-/LO/' + z_def + '/' + \
                        fragmentation_set + '/scale_variation/' + PDF_set + '/central/' + \
                        str(eta_lept_index) + '_vals.txt', delimiter=',')[pTD_index, :])
    
    Rcpm = Wp / Wm

    print(Rcpm)



#pTD_plot(['CT18ANLO', 'MSHT20nlo_as118', 'NNPDF40_nlo_pch_as_01180'], True, theory_labels)
#pTD_plot(['CT18ANLO', 'CT18ANNLO'], False, ['CT18ANLO', 'CT18ANNLO'])
#pTD_effect_of_subtraction_plot()
#pTD_dynamic_FF_scale('W-D+', 'CT18ANLO', 58, 'KKKS08_opal', 'minus')
#pTD_varying_FF_fit('CT18ANLO', 58, process, True, ['KKKS08_opal', 'KKKS08_global', 'SMSKA19'], ['KKKS08 OPAL', 'KKKS08 GLOBAL', 'SMSKA19'])
#eta_lept_plot()
#etaD_plot()
#z_variation("NLO", 'CT18ANLO')
#z_def_difference(process, 'CT18ANLO', 58, False)
#Rcpm('both')
#Rcpm('D')
#Rcpm('Dstar')
#Rcpm_pp_pPb()
#total_cross_section()
Rcpm_bin_integrated('pTD', 'both', True, ['CT18ANLO', 'MSHT20nlo_as118', 'NNPDF40_nlo_pch_as_01180'])
#Rcpm_pTD_varying_FF_fit('CT18ANLO', 'both')
#Rcpm_LO('D', 'MSHT20nlo_as118', 'minus', 'KKKS08_opal')
#Rcpm_LO_bin_integrated('pTD', 'both', 'MSHT20nlo_as118')
