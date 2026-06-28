import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import copy
import os
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
# Do subtraction
subtraction_flag = True
#
process = 'W-D+'

num_etac_bins = 22

if (fragmentation_set == 'KKKS08_opal'):
    frag_set_text = 'KKKS08 OPAL'
else:
    frag_set_text = ' KKKS08 GLOBAL'

PDF_sets = ['NNPDF30_nlo_as_01180_pp', 'NNPDF30_nlo_as_01180_pPb', 'EPPS21nlo_pp', 'EPPS21nlo_pPb']
#PDF_sets = ['EPPS21nlo_pp', 'EPPS21nlo_pPb']
#theory_labels = ['pp', 'pPb']
theory_labels = ['nNNPDF3.0NLO', 'EPPS21NLO']
num_err_members_in_sets = [200, 200, 106, 106]
#num_err_members_in_sets = [10, 10, 10, 10]

pdf_centrals = [[np.zeros((284, num_etac_bins)) for _ in range(len(PDF_sets))] for _ in range(4)]

markers = ['v', 's', 'o', 'd']
marker_color = 'black'

scale_var_color = 'cornflowerblue'
pdf_err_color = 'salmon'

if (process == "W+D-"):
    process_text = "$W^+D^-$"
elif (process == "W-D+"):
    process_text = "$W^-D^+$"
elif (process == "W+Dstar-"):
    process_text = "$W^+D^{*-}$"
else:
    process_text = "$W^-D^{*+}$"

main_vals_directory = '/home/alankovh/Documents/WD_production/output/'
plots_directory = '/home/alankovh/Documents/WD_production/plots/8,5 TeV/'


# W-D+  W+D-  W-Dstar+  W+Dstar-
atlas_vals = [50.2, 48.5, 51.1, 50.0]
atlas_stat_errs = [0.2, 0.2, 0.4, 0.4]
L_atlas = 140. * 10**3

efficiency = np.zeros(4)

L_new_experiment = 1.2 * 1e3

for i in range(4):
    efficiency[i] = atlas_vals[i] / (atlas_stat_errs[i]**2 * L_atlas)

print(efficiency)


def compute_general_2D_vals_NLO(PDF_set, num_err_members_in_set, process, load_scale_variation, load_pdf_errs,
                                subtraction_flag, FF_scale_choice, z_def, fragmentation_set):
    central_vals = np.zeros((284, num_etac_bins))
    dd_vals = np.zeros((284, num_etac_bins))
    uu_vals = np.zeros((284, num_etac_bins))
    scales_vals =[central_vals, dd_vals, uu_vals]
    pdf_err_plus = np.zeros((284, num_etac_bins))
    pdf_err_minus = np.zeros((284, num_etac_bins))

    central_MCerrs = np.zeros((284, num_etac_bins))
    dd_MCerrs = np.zeros((284, num_etac_bins))
    uu_MCerrs = np.zeros((284, num_etac_bins))
    scales_MCerrs =[central_MCerrs, dd_MCerrs, uu_MCerrs]

    for pTD_index in range(284):
        for etaD_index in range(num_etac_bins):
            pdf_err_plus[pTD_index, etaD_index] = 0
            pdf_err_minus[pTD_index, etaD_index] = 0

    stop = False

    scaling = 1. / 2.

    scale_names = ['central', 'dd', 'uu']
    # Get central and scale_variation values and MC errors. The central value will be corrected below for MC PDF sets.
    for scale_index in range(3):
        if (load_scale_variation or scale_index == 0):
            try:
                scales_vals[scale_index] = np.loadtxt(
                    main_vals_directory + process + '/NLO/' + z_def + '/' + \
                    fragmentation_set + '/' + FF_scale_choice + '/scale_variation/' + PDF_set + '/' + scale_names[scale_index] + '/0_m_charm_central_vals.txt', delimiter=',') * scaling
            except FileNotFoundError:
                stop = True
                print(main_vals_directory + process + '/NLO/' + z_def + '/' + \
                    fragmentation_set + '/' + FF_scale_choice + '/scale_variation/' + PDF_set + '/' + scale_names[scale_index] + '/0_m_charm_central_vals.txt')
                break
            
            scales_MCerrs[scale_index] = np.loadtxt(
                main_vals_directory + process + '/NLO/' + z_def + '/' + \
                fragmentation_set + '/' + FF_scale_choice + '/scale_variation/' + PDF_set + '/' + scale_names[scale_index] + '/0_m_charm_central_errs.txt', delimiter=',') * scaling
            
            if (subtraction_flag is True):
                scales_vals[scale_index] -= np.loadtxt(
                    main_vals_directory + process + '/subtraction/' + z_def + '/' + \
                    fragmentation_set + '/' + FF_scale_choice + '/scale_variation/' + PDF_set + '/' + scale_names[scale_index] + '/0_m_charm_central_vals.txt', delimiter=',') * scaling

                scales_MCerrs[scale_index] -= np.loadtxt(
                    main_vals_directory + process + '/subtraction/' + z_def + '/' + \
                    fragmentation_set + '/' + FF_scale_choice + '/scale_variation/' + PDF_set + '/' + scale_names[scale_index] + '/0_m_charm_central_errs.txt', delimiter=',') * scaling
            
        if (stop):
            break
        
    # Get pdf err values
    if (load_pdf_errs):
        pdf_errs_central = np.loadtxt(main_vals_directory + process + '/NLO/' + z_def + '/' + \
                                        fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                        str(0) + '_0_m_charm_central_vals.txt', delimiter=',') * scaling
        
        if (subtraction_flag is True):
            pdf_errs_central -= np.loadtxt(main_vals_directory + process + '/subtraction/' + z_def + '/' + \
                                        fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                        str(0) + '_0_m_charm_central_vals.txt', delimiter=',') * scaling
        
        if (PDF_set == 'EPPS21nlo_Ep' or PDF_set == 'EPPS21nlo_EPb'):
            # This loops through error members, but isn't directly the member id.
            for i in range(1, int(num_err_members_in_set / 2 + 1)):
                pdf_err_member_plus = np.zeros((284, num_etac_bins))
                pdf_err_member_minus = np.zeros((284, num_etac_bins))

                pdf_err_member_plus = np.loadtxt(main_vals_directory + process + '/NLO/' + z_def + '/' + \
                                        fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                        str(2 * (i - 1) + 1) + '_0_m_charm_central_vals.txt', delimiter=',') * scaling

                pdf_err_member_minus = np.loadtxt(main_vals_directory + process + '/NLO/' + z_def + '/' + \
                                        fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                        str(2 * i) + '_0_m_charm_central_vals.txt', delimiter=',') * scaling
                
                if (subtraction_flag is True):
                    pdf_err_member_plus -= np.loadtxt(main_vals_directory + process + '/subtraction/' + z_def + '/' + \
                                            fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                            str(2 * (i - 1) + 1) + '_0_m_charm_central_vals.txt', delimiter=',') * scaling
                    pdf_err_member_minus -= np.loadtxt(main_vals_directory + process + '/subtraction/' + z_def + '/' + \
                                            fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                            str(2 * i) + '_0_m_charm_central_vals.txt', delimiter=',') * scaling

                for pTD_index in range(284):
                    for etaD_index in range(num_etac_bins):
                        pdf_err_plus[pTD_index, etaD_index] += max(
                            pdf_err_member_plus[pTD_index, etaD_index] - \
                            pdf_errs_central[pTD_index, etaD_index],
                            pdf_err_member_minus[pTD_index, etaD_index] - \
                            pdf_errs_central[pTD_index, etaD_index], 0)**2
                        
                        pdf_err_minus[pTD_index, etaD_index] += max(
                            pdf_errs_central[pTD_index, etaD_index] - \
                            pdf_err_member_plus[pTD_index, etaD_index],
                            pdf_errs_central[pTD_index, etaD_index] - \
                            pdf_err_member_minus[pTD_index, etaD_index], 0)**2

            pdf_err_plus = np.sqrt(pdf_err_plus) / 1.645
            pdf_err_minus = np.sqrt(pdf_err_minus) / 1.645

        if (PDF_set == 'NNPDF30_nlo_as_01180_Np' or PDF_set == 'NNPDF30_nlo_as_01180_NPb'):
            average = 0.
            n_sum = 0
            for member in range(1, int(num_err_members_in_set + 1)):
                average += np.loadtxt(main_vals_directory + process + '/NLO/' + z_def + '/' + \
                                                        fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                                        str(member) + '_0_m_charm_central_vals.txt', delimiter=',') * scaling
                
                if (subtraction_flag is True):
                    average -= np.loadtxt(main_vals_directory + process + '/subtraction/' + z_def + '/' + \
                                            fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                            str(member) + '_0_m_charm_central_vals.txt', delimiter=',') * scaling

            average = average / num_err_members_in_set
            
            scales_vals[0] = average + scales_vals[0] - pdf_errs_central

            sum_val = np.zeros((284, num_etac_bins))

            for member in range(1, int(num_err_members_in_set + 1)):
                pdf_err_member = np.loadtxt(main_vals_directory + process + '/NLO/' + z_def + '/' + \
                                                        fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                                        str(member) + '_0_m_charm_central_vals.txt', delimiter=',') * scaling
                                                    
                if (subtraction_flag is True):
                    pdf_err_member -= np.loadtxt(main_vals_directory + process + '/subtraction/' + z_def + '/' + \
                                                fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                                str(member) + '_0_m_charm_central_vals.txt', delimiter=',') * scaling

                sum_val += (average - pdf_err_member)**2

            pdf_err_plus = np.sqrt(1. / (num_err_members_in_set * 1. - 1.) * sum_val)
            pdf_err_minus = pdf_err_plus
        


    return scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus


def compute_normalized_2D_values_for_a_pdf_member(process, PDF_index, PDF_set, member_index, FF_scale_choice):
    process_index = -1
    if (process == 'W-D+'):
        process_index = 0
    elif (process == 'W-Dstar+'):
        process_index = 1
    elif (process == 'W+D-'):
        process_index = 2
    else:
        process_index = 3

    num_err_members = num_err_members_in_sets[PDF_index]
    # Here "normalized" means that the difference between the member value and the central value is added to the central value of the scale variation run.
    member_vals_normalized = np.zeros((284, num_etac_bins))


    scale_var_central = np.loadtxt(main_vals_directory + process + '/NLO/' + z_def + '/' + \
                                            fragmentation_set + '/frag_main_scale/scale_variation/' + PDF_set + '/central/' + \
                                            '0_m_charm_central_vals.txt', delimiter=',') / 2. - \
                            np.loadtxt(main_vals_directory + process + '/subtraction/' + z_def + '/' + \
                                            fragmentation_set + '/' + FF_scale_choice + '/scale_variation/' + PDF_set + '/central/' + \
                                            '0_m_charm_central_vals.txt', delimiter=',') / 2.
    pdf_central = np.zeros((284, num_etac_bins))

    if (PDF_set == 'EPPS21nlo_Ep' or PDF_set == 'EPPS21nlo_EPb'):
        pdf_central = np.loadtxt(main_vals_directory + process + '/NLO/' + z_def + '/' + \
                                            fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                            str(0) + '_0_m_charm_central_vals.txt', delimiter=',') / 2. - \
                            np.loadtxt(main_vals_directory + process + '/subtraction/' + z_def + '/' + \
                                            fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                            str(0) + '_0_m_charm_central_vals.txt', delimiter=',') / 2.
    else:
        if (member_index < 2):
            member_sum = np.zeros((284, num_etac_bins))

            for member_index_here in range(1, num_err_members + 1):
                member_sum += np.loadtxt(main_vals_directory + process + '/NLO/' + z_def + '/' + \
                                                fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                                str(member_index_here) + '_0_m_charm_central_vals.txt', delimiter=',') / 2. - \
                                np.loadtxt(main_vals_directory + process + '/subtraction/' + z_def + '/' + \
                                                fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                                str(member_index_here) + '_0_m_charm_central_vals.txt', delimiter=',') / 2.

            pdf_central = member_sum / num_err_members

            pdf_centrals[process_index][PDF_index] = pdf_central
        else:
            pdf_central = pdf_centrals[process_index][PDF_index]


    member_vals = np.loadtxt(main_vals_directory + process + '/NLO/' + z_def + '/' + \
                                    fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                    str(member_index) + '_0_m_charm_central_vals.txt', delimiter=',') / 2. - \
                    np.loadtxt(main_vals_directory + process + '/subtraction/' + z_def + '/' + \
                                    fragmentation_set + '/pdf_errs/' + PDF_set + '/' + \
                                    str(member_index) + '_0_m_charm_central_vals.txt', delimiter=',') / 2.

    member_vals_normalized = member_vals - pdf_central + scale_var_central
    #member_vals_normalized = scale_var_central
    
    return member_vals_normalized


def compute_LO_integrated_cross_section(PDF_set, process, z_def, fragmentation_set):
    return sum(sum(np.loadtxt(main_vals_directory + process + '/LO/' + z_def + '/' + \
                    fragmentation_set + '/scale_variation/' + PDF_set + '/central/0_m_charm_central_vals.txt', delimiter=',')))


def pTD_plot(PDF_sets, plot_errors_flag, theory_labels):
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

    places_inside_bins = np.zeros((len(PDF_sets), 5))
    places_inside_bins_right = np.zeros((len(PDF_sets), 5))
    places_inside_bins_left = np.zeros((len(PDF_sets), 5))
    places_inside_bins_linear_scale = np.zeros((len(PDF_sets), 5))

    for i in range(len(PDF_sets)):
        places_inside_bins[i, :] = pTD_bins[:-1] * (pTD_bins[1:] / pTD_bins[:-1])**((i * 1. + 1) / (len(PDF_sets) * 1. + 1.))
        places_inside_bins_linear_scale[i, :] = pTD_bins[:-1] + (i * 1. + 1.) / (len(PDF_sets) * 1. + 1.) * (pTD_bins[1:] - pTD_bins[:-1])

    bar_width_over_bin_width = 1. / 4.
    bar_width = bar_width_over_bin_width * bin_widths
    shifts = np.zeros(5)

    for i in range(5):
        shifts[i] = ((pTD_bins[i + 1] * 1.) / (pTD_bins[i] * 1.))**(bar_width_over_bin_width / 2.)

    for i in range(len(PDF_sets)):
        places_inside_bins_left[i, :] = places_inside_bins[i, :] / shifts
        places_inside_bins_right[i, :] = places_inside_bins[i, :] * shifts

    scalings = [0.9, 1., 1.18]

    QCD_order = 'NLO'

    ratios = np.zeros(5)

    for PDF_index in range(len(PDF_sets)):
        scales_vals = [np.zeros((284, num_etac_bins)) for _ in range(3)]
        PDF_set = PDF_sets[PDF_index]
        num_err_members_in_set = num_err_members_in_sets[PDF_index]

        if (plot_errors_flag):
            scales_vals_raw, scales_MCerrs, pdf_err_plus, pdf_err_minus = compute_general_2D_vals_NLO(PDF_set, num_err_members_in_set,
                                                                        process, True, True, True, 'frag_main_scale', z_def, fragmentation_set)
        else:
            scales_vals_raw, scales_MCerrs, pdf_err_plus, pdf_err_minus = compute_general_2D_vals_NLO(PDF_set, num_err_members_in_set,
                                                                        process, False, False, True, 'frag_main_scale', z_def, fragmentation_set)

        for i in range(3):
            scales_vals[i] += scales_vals_raw[i]

        HISTO_central_sigma_vals = np.zeros(5)
        HISTO_scales_dd_sigma_vals = np.zeros(5)
        HISTO_scales_uu_sigma_vals = np.zeros(5)

        HISTO_central_sigma_MCerrs = np.zeros(5)
        HISTO_scales_dd_sigma_MCerrs = np.zeros(5)
        HISTO_scales_uu_sigma_MCerrs = np.zeros(5)

        HISTO_pdf_err_plus = np.zeros(5)
        HISTO_pdf_err_minus = np.zeros(5)

        scale_names = ['central', 'dd', 'uu']

        bin_index = 0
        for pTD_index in range(284 + 1):
            if (pTD_data_min + (pTD_index + 1 / 2) * pTD_data_bin_width > pTD_bins[bin_index + 1]):
                if (bin_index < 4):
                    bin_index += 1
                else:
                    break

            HISTO_central_sigma_vals[bin_index] += sum(scales_vals[0][pTD_index, :])
            HISTO_scales_dd_sigma_vals[bin_index] += sum(scales_vals[1][pTD_index, :])
            HISTO_scales_uu_sigma_vals[bin_index] += sum(scales_vals[2][pTD_index, :])

            HISTO_central_sigma_MCerrs[bin_index] += sum(scales_MCerrs[0][pTD_index, :])
            HISTO_scales_dd_sigma_MCerrs[bin_index] += sum(scales_MCerrs[1][pTD_index, :])
            HISTO_scales_uu_sigma_MCerrs[bin_index] += sum(scales_MCerrs[2][pTD_index, :])

            HISTO_pdf_err_plus[bin_index] += sum(pdf_err_plus[pTD_index, :])
            HISTO_pdf_err_minus[bin_index] += sum(pdf_err_minus[pTD_index, :])
        
        if (PDF_index == 0):
            ratios = HISTO_central_sigma_vals
        else:
            ratios = ratios / HISTO_central_sigma_vals

        ax1.plot(places_inside_bins[PDF_index, :], HISTO_central_sigma_vals, marker=markers[PDF_index],
                        color=marker_color, markersize=5, linestyle='none',
                        label=theory_labels[PDF_index], zorder=4)

        if (plot_errors_flag):
            ax1.bar(places_inside_bins[PDF_index, :], HISTO_pdf_err_plus + HISTO_pdf_err_minus, width=bar_width * scalings[PDF_index],
                        bottom=HISTO_central_sigma_vals - HISTO_pdf_err_minus, color=pdf_err_color,
                        zorder=3)

            ax1.bar(places_inside_bins[PDF_index, :], HISTO_scales_uu_sigma_vals - HISTO_scales_dd_sigma_vals, width=bar_width * scalings[PDF_index],
                        bottom=HISTO_central_sigma_vals + HISTO_scales_dd_sigma_vals, color=scale_var_color,
                        zorder=2)

    
    ratios = ratios**(-1)

    ax2.hlines(ratios, pTD_bins[:-1], pTD_bins[1:], color='black', zorder=1)


    """    
    
    
    ratios_ATLAS_var_down = np.zeros(5)
    ratios_ATLAS_var_up = np.zeros(5)
    ratios_theory_scale_var_var_down = np.zeros(5)
    ratios_theory_scale_var_var_up = np.zeros(5)
    ratios_theory_pdf_err_var_down = np.zeros(5)
    ratios_theory_pdf_err_var_up = np.zeros(5)

    for pTD_index in range(5):
        ratios[pTD_index] = HISTO_central_sigma_vals[eta_lept_index] / atlas_vals[atlas_index][eta_lept_index]

        if (plot_errors_flag):
            ratios_theory_scale_var_var_down[pTD_index] = (HISTO_central_sigma_vals[eta_lept_index] + HISTO_scales_dd_sigma_vals[eta_lept_index]) / \
                                                                HISTO_central_sigma_vals[eta_lept_index]
                                                
            ratios_theory_scale_var_var_up[pTD_index] = (HISTO_central_sigma_vals[eta_lept_index] + HISTO_scales_uu_sigma_vals[eta_lept_index]) / \
                                                                HISTO_central_sigma_vals[eta_lept_index]
            
            ratios_theory_pdf_err_var_down[pTD_index] = (HISTO_central_sigma_vals[eta_lept_index] - HISTO_pdf_err_minus[eta_lept_index]) / \
                                                                HISTO_central_sigma_vals[eta_lept_index]
                                                
            ratios_theory_pdf_err_var_up[pTD_index] = (HISTO_central_sigma_vals[eta_lept_index] + HISTO_pdf_err_plus[eta_lept_index]) / \
                                                                HISTO_central_sigma_vals[eta_lept_index]
                                                
    ax2.plot(places_inside_bins[PDF_index, :], ratios, zorder=4, marker=markers[PDF_index],
                color=marker_color, markersize=5, linestyle='none')

    if (plot_errors_flag):
        scale_var_bar_plot = ax2.bar(places_inside_bins_right[PDF_index, :], ratios_theory_scale_var_var_up - ratios_theory_scale_var_var_down,
                bottom=ratios - 1. + ratios_theory_scale_var_var_down, width=bar_width * scalings[PDF_index], color=scale_var_color,
                linewidth=1, zorder=3)
        
        pdf_err_bar_plot = ax2.bar(places_inside_bins_left[PDF_index, :], ratios_theory_pdf_err_var_up - ratios_theory_pdf_err_var_down,
                bottom=ratios - 1. + ratios_theory_pdf_err_var_down, width=bar_width * scalings[PDF_index], color=pdf_err_color,
                linewidth=1, zorder=3)
    
    ax1.hlines(atlas_vals[atlas_index], pTD_bins[:-1], pTD_bins[1:], color='black', label='ATLAS', zorder=1)
    # Plot error bars for the Atlas values.
    """

    plt.xscale('log', base=10)
    ax1.set_ylim(0, 20)
    plt.xlim(8, 150)
    ax2.set_ylim(0.95, 1.05)

    ax2.set_xlabel(r'$P_T(D)$ [GeV]', fontsize=axis_label_font_size)
    ax1.set_ylabel(r'$\mathrm{Cross}\ \mathrm{section}\ \mathrm{[pb]}$', fontsize=axis_label_font_size)
    ax2.set_ylabel(r'$R_A$', 
               fontsize=axis_label_font_size * 1.3)

    ax1.set_yticks([5, 10, 15, 20])
    ax2.set_yticks([0.95, 1., 1.05, 1.1])

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

    #legend1 = ax1.legend(loc='center right', framealpha=1, fontsize=legend_fontsize, bbox_to_anchor=(0.98, 0.8))
    #if (plot_errors_flag):
    #    legend2 = ax1.legend([ATLAS_uncertainty, pdf_err_bar_plot, scale_var_bar_plot],
    #                        ["ATLAS uncertainty", "PDF uncertainty", "Scale variation"], loc='lower left', framealpha=1, fontsize=legend_fontsize)
    #    ax1.add_artist(legend1)

    info_y_vals_1 = 18
    info_y_vals_2 = 16
    info_y_vals_3 = 14

    info_x_vals_1 = 9
    info_x_vals_2 = 25

    ax1.text(info_x_vals_1, info_y_vals_1, process_text + '  OS-SS', fontsize=font_size)
    ax1.text(info_x_vals_1, info_y_vals_2, r'$\sqrt{s} = 13$ TeV', fontsize=font_size)
    #ax1.text(info_x_vals_1, info_y_vals_3, frag_set_text, fontsize=font_size)

    ax1.tick_params(axis='both', which='major', labelsize=axis_font_size)
    ax2.tick_params(axis='both', which='major', labelsize=axis_font_size)

    plt.xticks(pTD_bins, [f'{tick:.0f}' for tick in pTD_bins])

    for i in range(1, len(pTD_bins) - 1):
        ax1.axvline(pTD_bins[i], linestyle='dashed', color='black', linewidth=0.5, ymax=0.57)

    for i in range(1, len(pTD_bins) - 1):
        ax2.axvline(pTD_bins[i], linestyle='dashed', color='black', linewidth=0.5, ymax=1)

    ax1.legend(loc='upper right', framealpha=1, fontsize=legend_fontsize)
    plt.tight_layout()

    filename = plots_directory + process + '/' + fragmentation_set + "/" + process + '_pT.pdf'
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    plt.savefig(filename)

    plt.show()


def etaD_plot(PDF_sets, plot_errors_flag, theory_labels):
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1], 'hspace': 0}, figsize=(7, 7))

    font_size = 16
    axis_label_font_size = 17
    axis_font_size = 13
    legend_fontsize = 14

    etaD_bins = np.arange(-2.2, 2.201, 0.2)
    xmin = etaD_bins[0:-1] + 0.1
    xmax = etaD_bins[1:] - 0.1
    bin_midpoints = (xmin + xmax) / 2
    theory_val_places = [bin_midpoints - 0.2/4., bin_midpoints, bin_midpoints + 0.2/4.]

    QCD_order = 'NLO'

    scales_vals = [np.zeros((284, num_etac_bins)) for _ in range(3)]

    HISTO_ratio = np.zeros(len(etaD_bins) - 1)

    for PDF_index in range(len(PDF_sets)):
        PDF_set = PDF_sets[PDF_index]
        num_err_members_in_set = num_err_members_in_sets[PDF_index]

        if (plot_errors_flag):
            scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus = compute_general_2D_vals_NLO(PDF_set, num_err_members_in_set,
                                                                        process, True, True, True, 'frag_main_scale', z_def, fragmentation_set)
        else:
            scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus = compute_general_2D_vals_NLO(PDF_set, num_err_members_in_set,
                                                                        process, False, False, True, 'frag_main_scale', z_def, fragmentation_set)

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

        for etaD_index in range(len(etaD_bins) - 1):
            HISTO_central_sigma_vals[etaD_index] += sum(scales_vals[0][:, etaD_index])
            HISTO_scales_dd_sigma_vals[etaD_index] += sum(scales_vals[1][:, etaD_index])
            HISTO_scales_uu_sigma_vals[etaD_index] += sum(scales_vals[2][:, etaD_index])

            HISTO_central_sigma_MCerrs[etaD_index] += sum(scales_MCerrs[0][:, etaD_index])
            HISTO_scales_dd_sigma_MCerrs[etaD_index] += sum(scales_MCerrs[1][:, etaD_index])
            HISTO_scales_uu_sigma_MCerrs[etaD_index] += sum(scales_MCerrs[2][:, etaD_index])

            HISTO_pdf_err_plus[etaD_index] += sum(pdf_err_plus[:, etaD_index])
            HISTO_pdf_err_minus[etaD_index] += sum(pdf_err_minus[:, etaD_index])

        HISTO_ratio_up_pdf_err = (HISTO_central_sigma_vals + HISTO_pdf_err_plus) / HISTO_central_sigma_vals
        HISTO_ratio_down_pdf_err = (HISTO_central_sigma_vals - HISTO_pdf_err_minus) / HISTO_central_sigma_vals

        HISTO_ratio_up_scale_var = (HISTO_central_sigma_vals + HISTO_scales_dd_sigma_vals) / HISTO_central_sigma_vals
        HISTO_ratio_down_scale_var = (HISTO_central_sigma_vals + HISTO_scales_uu_sigma_vals) / HISTO_central_sigma_vals

        if (PDF_index == 0):
            HISTO_ratio = HISTO_central_sigma_vals
        else:
            HISTO_ratio = HISTO_ratio / HISTO_central_sigma_vals

        colors_here = ['red', 'black', 'orange', 'purple']

        ax1.plot(bin_midpoints, HISTO_central_sigma_vals, marker=markers[PDF_index],
                        color=colors_here[PDF_index], markersize=5, linestyle='none',
                        label=theory_labels[PDF_index], zorder=4)

        if (plot_errors_flag):
            #print(HISTO_central_sigma_vals)
            bar_width = 0.04

            pdf_err_bar_plot = ax1.bar(theory_val_places[PDF_index] - bar_width / 4., HISTO_pdf_err_plus + HISTO_pdf_err_minus, width=bar_width / 2.,
                        bottom=HISTO_central_sigma_vals - HISTO_pdf_err_minus, color=pdf_err_color, zorder=2)
            scale_var_bar_plot = ax1.bar(theory_val_places[PDF_index] + bar_width / 4., HISTO_scales_uu_sigma_vals - HISTO_scales_dd_sigma_vals, width=bar_width / 2.,
                        bottom=HISTO_central_sigma_vals + HISTO_scales_dd_sigma_vals, color=scale_var_color, zorder=2)
            
            ax2.bar(theory_val_places[PDF_index] - bar_width / 4., HISTO_ratio_up_pdf_err - HISTO_ratio_down_pdf_err,
                        width=bar_width / 2., bottom=HISTO_ratio_down_pdf_err, color=pdf_err_color, zorder=2)
            ax2.bar(theory_val_places[PDF_index] + bar_width / 4., HISTO_ratio_up_scale_var - HISTO_ratio_down_scale_var,
                        width=bar_width / 2., bottom=HISTO_ratio_down_scale_var, color=scale_var_color, zorder=2)
    
    HISTO_ratio = 1. / HISTO_ratio

    ax2.plot(bin_midpoints, HISTO_ratio, color='black', marker=markers[1], linestyle='none')

    plt.xlabel(r'$|\eta_D|$', fontsize=axis_label_font_size)
    ax1.set_ylabel('Cross section [pb]', fontsize=axis_label_font_size)
    ax2.set_ylabel(r'$\frac{\text{pPb}}{\text{pp}}$', fontsize=axis_label_font_size * 1.3)

    info_xval_1 = 0.1
    info_xval_2 = 0.8
    info_yval_1 = 9.7
    info_yval_2 = 8.9
    info_yval_3 = 8.1

    ax1.tick_params(axis='both', which='major', labelsize=axis_font_size)
    ax1.text(info_xval_1, info_yval_1, process_text + '  OS-SS', fontsize=font_size)
    ax1.text(info_xval_1, info_yval_2, r'$\sqrt{s} = 13$ TeV', fontsize=font_size)
    #ax1.text(info_xval_1, info_yval_3, frag_set_text, fontsize=font_size)

    plt.xlim(0, 2.2)
    ax1.set_ylim(0.5, 3.)
    ax2.set_ylim(0.85, 1.15)

    ax2.plot([-3, 5], [1, 1], color='black', zorder=0)

    plt.xticks(etaD_bins, rotation=45, ha='right')

    for i in range(1, len(etaD_bins) - 1):
        ax1.axvline(etaD_bins[i], linestyle='dashed', color='black', linewidth=0.5, ymax=0.6)
        ax2.axvline(etaD_bins[i], linestyle='dashed', color='black', linewidth=0.5, ymax=1)

    ax1.tick_params(direction='in', top=True, right=True)
    ax1.minorticks_on()
    ax1.tick_params(which='both', direction='in', top=True, right=True)
    ax2.minorticks_on()
    ax2.tick_params(which='both', direction='in', top=True, right=True)
    ax2.tick_params(direction='in', top=True, right=True)

    plt.text(-2, 1.9, process, fontsize=font_size)
    
    legend1 = ax1.legend(loc='upper right', framealpha=1, fontsize=legend_fontsize)
    #legend2 = ax1.legend([pdf_err_bar_plot, scale_var_bar_plot], ["PDF uncertainty", "Scale variation"], loc='lower left', framealpha=1, fontsize=legend_fontsize)
    #ax1.add_artist(legend1)

    filename = plots_directory + process + '/' + fragmentation_set + "/" + process + '_' + PDF_sets[1] + '_etaD.pdf'
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    plt.savefig(filename)
    plt.tight_layout()
    plt.show()


def Rcpm(PDF_sets):
    font_size = 16
    axis_label_font_size = 17
    axis_font_size = 13
    legend_fontsize = 14

    fig, ax = plt.subplots(figsize=(6, 6))

    for PDF_index in range(len(PDF_sets)):
        PDF_set = PDF_sets[PDF_index]
        num_err_members_in_set = num_err_members_in_sets[PDF_index]

        process_here = "W-D+"
        scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus = compute_general_2D_vals_NLO(PDF_set, num_err_members_in_set,
                                                                    process_here, False, False, True, 'frag_main_scale', z_def, fragmentation_set)

        Wm_cross_section = sum(sum(scales_vals[0]))
        print(sum(sum(scales_vals[0])))

        process_here = "W-Dstar+"
        scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus = compute_general_2D_vals_NLO(PDF_set, num_err_members_in_set,
                                                                    process_here, False, False, True, 'frag_main_scale', z_def, fragmentation_set)
        Wm_star_cross_section = sum(sum(scales_vals[0]))
        print(sum(sum(scales_vals[0])))

        process_here = "W+D-"
        scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus = compute_general_2D_vals_NLO(PDF_set, num_err_members_in_set,
                                                                    process_here, False, False, True, 'frag_main_scale', z_def, fragmentation_set)
        Wp_cross_section = sum(sum(scales_vals[0]))
        print(sum(sum(scales_vals[0])))

        process_here = "W+Dstar-"
        scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus = compute_general_2D_vals_NLO(PDF_set, num_err_members_in_set,
                                                                    process_here, False, False, True, 'frag_main_scale', z_def, fragmentation_set)
        Wp_star_cross_section = sum(sum(scales_vals[0]))
        print(sum(sum(scales_vals[0])))

        #print(Wm_cross_section)
        #print(Wm_star_cross_section)
        #print(Wp_cross_section)
        #print(Wp_star_cross_section)
        #print()

        Rcpm = (Wp_cross_section + Wp_star_cross_section) / (Wm_cross_section + Wm_star_cross_section)
        print(PDF_set, Rcpm)


def Rcpm_LO(PDF_set, z_def, fragmentation_set):
    Wp = 0.
    Wm = 0.

    Wm = compute_LO_integrated_cross_section(PDF_set, 'W-D+', z_def, fragmentation_set)
    Wm += compute_LO_integrated_cross_section(PDF_set, 'W-Dstar+', z_def, fragmentation_set)
    Wp = compute_LO_integrated_cross_section(PDF_set, 'W+D-', z_def, fragmentation_set)
    Wp += compute_LO_integrated_cross_section(PDF_set, 'W+Dstar-', z_def, fragmentation_set)

    print(Wp / Wm)


def total_cross_section(W_sign):
    font_size = 19
    axis_label_font_size = 22
    axis_font_size = 17
    legend_fontsize = 17

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [4, 7], 'hspace': 0}, figsize=(6, 6))

    y_vals = [2.8, 2.2, 1.3, 0.7]

    if (W_sign == 'minus'):
        processes = ['W-D+', 'W-Dstar+']
    else:
        processes = ['W+D-', 'W+Dstar-']

    for PDF_set_index in range(2, len(PDF_sets)):
        PDF_set = PDF_sets[PDF_set_index]

        vals = np.zeros(num_err_members_in_sets[PDF_set_index] + 1)
        vals_star = np.zeros(num_err_members_in_sets[PDF_set_index] + 1)
        err_plus = 0.
        err_minus = 0.

        sum_quantity = 0.

        for member_index in range(num_err_members_in_sets[PDF_set_index] + 1):
            for process in processes:
                process_index = -1
                if (process == 'W-D+'):
                    process_index = 0
                elif (process == 'W-Dstar+'):
                    process_index = 1
                elif (process == 'W+D-'):
                    process_index = 2
                else:
                    process_index = 3

                num_err_members_in_set = num_err_members_in_sets[PDF_set_index]

                if (process_index == 0 or process_index == 2):
                    vals[member_index] += sum(sum(compute_normalized_2D_values_for_a_pdf_member(process, PDF_set_index, PDF_set, member_index, 'frag_main_scale'))) * 208.
                else:
                    vals_star[member_index] += sum(sum(compute_normalized_2D_values_for_a_pdf_member(process, PDF_set_index, PDF_set, member_index, 'frag_main_scale'))) * 208.
        
        sum_quantity = vals[0] + vals_star[0]
        
        if (PDF_set == 'EPPS21nlo_pp' or PDF_set == 'EPPS21nlo_pPb'):
            # This is not directly member_index.
            for member_index in range(1, int(num_err_members_in_sets[PDF_set_index] / 2) + 1):
                plus_val = vals[2 * (member_index - 1) + 1] + vals_star[2 * (member_index - 1) + 1]
                minus_val = vals[2 * member_index] + vals_star[2 * member_index]
                err_plus += max(plus_val - vals[0] - vals_star[0], vals[0] + vals_star[0] - minus_val, 0)**2
                err_minus += max(vals[0] + vals_star[0] - plus_val, minus_val - vals[0] - vals_star[0], 0)**2

                print(plus_val, member_index, PDF_set)
            err_plus = np.sqrt(err_plus) / 1.645
            err_minus = np.sqrt(err_minus) / 1.645

            sum_quantity = vals[0] + vals_star[0]

        else:
            sum_quantity = sum(vals[1:] + vals_star[1:]) / (len(vals) * 1. - 1.)

            sum_in_error_formula = np.sum((vals[1:] + vals_star[1:] - sum_quantity)**2)

            err_plus = np.sqrt(1. / (num_err_members_in_sets[PDF_set_index] * 1. - 1.) * sum_in_error_formula)
            err_minus = err_plus

            vals[0] = sum(vals[1:]) / (num_err_members_in_sets[PDF_set_index] * 1.)
            vals_star[0] = sum(vals_star[1:]) / (num_err_members_in_sets[PDF_set_index] * 1.)

        vals[0] = vals[0] * 1e-3
        vals_star[0] = vals_star[0] * 1e-3

        plt.plot([sum_quantity / 208., sum_quantity / 208.], [y_vals[PDF_set_index] - 0.25, y_vals[PDF_set_index] + 0.25], color='black', zorder=5, solid_capstyle='butt')

        err_expected = 0.

        if (PDF_set == 'NNPDF30_nlo_as_01180_pPb' or PDF_set == 'EPPS21nlo_pPb'):
            if (W_sign == 'minus'):
                efficiency_starless = efficiency[0]
                efficiency_star = efficiency[2]
            else:
                efficiency_starless = efficiency[1]
                efficiency_star = efficiency[3]

            N_expected = L_new_experiment * (efficiency_starless * vals[0] + efficiency_star * vals_star[0])
            err_expected = 1. / np.sqrt(N_expected * 1.) * sum_quantity

            expected_uncertainty = patches.Rectangle((sum_quantity / 208. - err_expected / 208., y_vals[PDF_set_index] - 0.25), 2. * err_expected / 208.,
                                                    0.25, facecolor='lightgray', zorder=3)

            ax2.add_patch(expected_uncertainty)
            
            print(err_minus / 208.)
            pdf_err = patches.Rectangle((sum_quantity / 208. - err_minus / 208., y_vals[PDF_set_index]), err_plus / 208. + err_minus / 208.,
                                        0.25, facecolor=pdf_err_color, zorder=4)

        else:
            pdf_err = patches.Rectangle((sum_quantity / 208. - err_minus / 208., y_vals[PDF_set_index] - 0.25), err_plus / 208. + err_minus / 208.,
                                        0.5, facecolor=pdf_err_color, zorder=4)

        ax2.add_patch(pdf_err)

        if (PDF_set == 'EPPS21nlo_pp' or PDF_set == 'NNPDF30_nlo_as_01180_pp'):
            ax2.text(sum_quantity / 208. - err_minus / 208. - 3.5, y_vals[PDF_set_index] - 0.1, r'$pp$', fontsize=font_size)
            if (PDF_set == 'NNPDF30_nlo_as_01180_pp'):
                ax2.text(21, y_vals[PDF_set_index] - 0.5, r'nNNPDF3.0NLO', fontsize=font_size)
            else:
                ax2.text(21, y_vals[PDF_set_index] - 0.5, r'EPPS21NLO', fontsize=font_size)
        else:
            ax2.text(sum_quantity / 208. - err_expected / 208. - 4, y_vals[PDF_set_index] - 0.1, r'$p$Pb', fontsize=font_size)

    plt.xlim(19, 58)
    ax2.set_ylim(0, 3.5)
    ax1.set_ylim(0.2, 4)

    ax2.set_xlabel('Cross section per nucleon [pb]', fontsize=axis_label_font_size)

    info_xval_1 = 21
    info_xval_2 = 1.
    info_yval_1 = 3
    info_yval_2 = 2
    info_yval_3 = 1

    if (W_sign == 'minus'):
        process_text = r'$W^-D^{(*)+}$'
    else:
        process_text = r'$W^+D^{(*)-}$'

    ax1.text(info_xval_1, info_yval_1, process_text + r'\quad OS-SS', fontsize=font_size)
    ax1.text(info_xval_1, info_yval_2, r'$\sqrt{s} = 8.5$ TeV', fontsize=font_size)
    ax1.text(info_xval_1, info_yval_3, frag_set_text, fontsize=font_size)

    # For legends
    expected_uncertainty = patches.Rectangle((0, 0), 0, 0, facecolor='lightgray', zorder=3)
    ax1.add_patch(expected_uncertainty)
    pdf_err = patches.Rectangle((0, 0), 0, 0, facecolor=pdf_err_color, zorder=4)
    ax1.add_patch(pdf_err)
    ax1.legend([pdf_err, expected_uncertainty], ["PDF error (68\% C.L.)", "Expected statistical\nmeasurement error\n(68\% C.L.)"],
                framealpha=1, fontsize=legend_fontsize, loc='lower right')

    ax2.tick_params(axis='both', which='major', labelsize=axis_font_size)
    ax1.set_yticklabels([])
    ax2.set_yticklabels([])

    ax2.minorticks_on()
    ax2.tick_params(
        which='both',    # major and minor ticks
        direction='in',
        left=False,
        right=False,
        bottom=True,
        top=True,
        labelsize=axis_font_size
    )
    ax1.tick_params(top=False, bottom=False, right=False, left=False)

    ax1.xaxis.set_zorder(100)
    ax1.yaxis.set_zorder(100)

    ax2.xaxis.set_zorder(100)
    ax2.yaxis.set_zorder(100)

    plt.tight_layout()

    filename = plots_directory + 'integrated/' + fragmentation_set + '/W_' + W_sign + '_integrated.pdf'
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    plt.savefig(filename)
    plt.show()


#pTD_plot(PDF_sets, True, ['pp', 'pPb'])
#etaD_plot(PDF_sets, False, ['pp', 'pPb'])
#Rcpm(PDF_sets)
total_cross_section('plus')

#Rcpm_LO('EPPS16nlo_pPb', 'minus', 'opal')