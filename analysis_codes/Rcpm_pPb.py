import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
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

main_vals_directory = '/home/alankovh/Documents/WD_production/output/'
plots_directory = '/home/alankovh/Documents/WD_production/plots/8,5 TeV/'

PDF_sets_raw = ['NNPDF30_nlo_as_01180', 'EPPS21nlo']
#PDF_sets_raw = ['EPPS21nlo']
num_err_members_in_sets = [200, 10]

num_etac_bins = 22
pdf_centrals = [[np.zeros((284, num_etac_bins)) for _ in range(len(PDF_sets_raw))] for _ in range(4)]

z_def = 'minus'
fragmentation_set = 'KKKS08_opal'

if (fragmentation_set == 'KKKS08_opal'):
    frag_set_text = 'KKKS08 OPAL'
else:
    frag_set_text = ' KKKS08 GLOBAL'

markers = [['v', '^'], ['s', 'o']]
marker_color = 'black'
theory_labels = [[r'NNPDF30_nlo ($pp$)', r'NNPDF30_nlo ($pPb$)'], [r'EPPS21nlo ($pp$)', r'EPPS21nlo ($pPb$)']]

scale_var_color = 'cornflowerblue'
pdf_err_color = 'salmon'

# W-D+  W+D-  W-Dstar+  W+Dstar-
atlas_vals = [50.2, 48.5, 51.1, 50.0]
atlas_stat_errs = [0.2, 0.2, 0.4, 0.4]
L_atlas = 140. * 10**3

efficiency = np.zeros(4)

L_new_experiment = 1.2

for i in range(4):
    efficiency[i] = atlas_vals[i] / (atlas_stat_errs[i]**2 * L_atlas)


def compute_normalized_2D_values_for_a_pdf_member(process, PDF_index, PDF_set, member_index):
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
                                            fragmentation_set + '/frag_main_scale/scale_variation/' + PDF_set + '/central/' + \
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

    return member_vals_normalized * 208.


def compute_Rcpm_pdf_err_MC(PDF_index, PDF_set, Rcpm_central, which_cross_sections_included):
    sum_in_error_formula = 0.
    average = 0.
    Rcpm_err_vals = np.zeros(num_err_members_in_sets[PDF_index])

    for member_index in range(1, num_err_members_in_sets[PDF_index] + 1):
        member_vals_normalized = compute_normalized_2D_values_for_a_pdf_member('W+D-', PDF_index, PDF_set, member_index)
        WpDm = sum(sum(member_vals_normalized))

        member_vals_normalized = compute_normalized_2D_values_for_a_pdf_member('W+Dstar-', PDF_index, PDF_set, member_index)
        WpDstarm = sum(sum(member_vals_normalized))

        member_vals_normalized = compute_normalized_2D_values_for_a_pdf_member('W-D+', PDF_index, PDF_set, member_index)
        WmDp = sum(sum(member_vals_normalized))

        member_vals_normalized = compute_normalized_2D_values_for_a_pdf_member('W-Dstar+', PDF_index, PDF_set, member_index)
        WmDstarp = sum(sum(member_vals_normalized))

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

    sum_in_error_formula = np.sum((Rcpm_err_vals - Rcpm_central)**2)

    Rcpm_err_plus = np.sqrt(1. / (num_err_members_in_sets[PDF_index] * 1. - 1.) * sum_in_error_formula) * 1.645
    Rcpm_err_minus = Rcpm_err_plus

    #if (PDF_set == 'NNPDF30_nlo_as_01180_Np'):
    #    plt.scatter(Rcpm_err_vals, np.zeros(len(Rcpm_err_vals)) + 2.2, color='black', zorder=10, s=0.1)
    #if (PDF_set == 'NNPDF30_nlo_as_01180_NPb'):
    #    plt.scatter(Rcpm_err_vals, np.zeros(len(Rcpm_err_vals)) + 1.8, color='black', zorder=10, s=0.1)

    return Rcpm_central, Rcpm_err_plus, Rcpm_err_minus


def compute_Rcpm_pdf_err_HESSIAN(PDF_index, PDF_set, Rcpm_central, which_cross_sections_included):
    Rcpm_err_plus = 0.
    Rcpm_err_minus = 0.

    for member_index in range(1, int(num_err_members_in_sets[PDF_index] / 2) + 1):
        member_vals_normalized = compute_normalized_2D_values_for_a_pdf_member('W+D-', PDF_index, PDF_set, 2 * (member_index - 1) + 1)
        WpDm_pdf_plus = sum(sum(member_vals_normalized))
        member_vals_normalized = compute_normalized_2D_values_for_a_pdf_member('W+D-', PDF_index, PDF_set, 2 * member_index)
        WpDm_pdf_minus = sum(sum(member_vals_normalized))

        member_vals_normalized = compute_normalized_2D_values_for_a_pdf_member('W+Dstar-', PDF_index, PDF_set, 2 * (member_index - 1) + 1)
        WpDstarm_pdf_plus = sum(sum(member_vals_normalized))
        member_vals_normalized = compute_normalized_2D_values_for_a_pdf_member('W+Dstar-', PDF_index, PDF_set, 2 * member_index)
        WpDstarm_pdf_minus = sum(sum(member_vals_normalized))

        member_vals_normalized = compute_normalized_2D_values_for_a_pdf_member('W-D+', PDF_index, PDF_set, 2 * (member_index - 1) + 1)
        WmDp_pdf_plus = sum(sum(member_vals_normalized))
        member_vals_normalized = compute_normalized_2D_values_for_a_pdf_member('W-D+', PDF_index, PDF_set, 2 * member_index)
        WmDp_pdf_minus = sum(sum(member_vals_normalized))

        member_vals_normalized = compute_normalized_2D_values_for_a_pdf_member('W-Dstar+', PDF_index, PDF_set, 2 * (member_index - 1) + 1)
        WmDstarp_pdf_plus = sum(sum(member_vals_normalized))
        member_vals_normalized = compute_normalized_2D_values_for_a_pdf_member('W-Dstar+', PDF_index, PDF_set, 2 * member_index)
        WmDstarp_pdf_minus = sum(sum(member_vals_normalized))

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

    return Rcpm_err_plus, Rcpm_err_minus


def Rcpm(which_cross_sections_included, PDF_sets_raw, PDF_prefix):
    font_size = 19
    axis_label_font_size = 22
    axis_font_size = 17
    legend_fontsize = 17

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [4, 7], 'hspace': 0}, figsize=(6, 6))

    y_vals = [2, 1]

    p_Pb_shift = [0.2, -0.2]
    proton_or_lead_array = ['pp', 'pPb']

    for PDF_index in range(1, len(PDF_sets_raw)):
        for proton_or_lead_index in range(2):
            proton_or_lead = proton_or_lead_array[proton_or_lead_index]
            PDF_set = PDF_sets_raw[PDF_index] + '_' + proton_or_lead
            num_err_members_in_set = num_err_members_in_sets[PDF_index]

            if (PDF_sets_raw[PDF_index] == 'EPPS21nlo'):
                process_here = "W-D+"
                Wm_cross_section = sum(sum(compute_normalized_2D_values_for_a_pdf_member(process_here, PDF_index, PDF_set, 0)))

                process_here = "W-Dstar+"
                Wm_star_cross_section = sum(sum(compute_normalized_2D_values_for_a_pdf_member(process_here, PDF_index, PDF_set, 0)))

                process_here = "W+D-"
                Wp_cross_section = sum(sum(compute_normalized_2D_values_for_a_pdf_member(process_here, PDF_index, PDF_set, 0)))

                process_here = "W+Dstar-"
                Wp_star_cross_section = sum(sum(compute_normalized_2D_values_for_a_pdf_member(process_here, PDF_index, PDF_set, 0)))
            else:
                Wm_cross_section = 0.
                Wm_star_cross_section = 0.
                Wp_cross_section = 0.
                Wp_star_cross_section = 0.

                for member_index in range(1, num_err_members_in_sets[PDF_index] + 1):
                    process_here = "W-D+"
                    Wm_cross_section += sum(sum(compute_normalized_2D_values_for_a_pdf_member(process_here, PDF_index, PDF_set, member_index)))

                    process_here = "W-Dstar+"
                    Wm_star_cross_section += sum(sum(compute_normalized_2D_values_for_a_pdf_member(process_here, PDF_index, PDF_set, member_index)))

                    process_here = "W+D-"
                    Wp_cross_section += sum(sum(compute_normalized_2D_values_for_a_pdf_member(process_here, PDF_index, PDF_set, member_index)))

                    process_here = "W+Dstar-"
                    Wp_star_cross_section += sum(sum(compute_normalized_2D_values_for_a_pdf_member(process_here, PDF_index, PDF_set, member_index)))
                
                Wm_cross_section = Wm_cross_section / (num_err_members_in_sets[PDF_index] * 1.)
                Wm_star_cross_section = Wm_star_cross_section / (num_err_members_in_sets[PDF_index] * 1.)
                Wp_cross_section = Wp_cross_section / (num_err_members_in_sets[PDF_index] * 1.)
                Wp_star_cross_section = Wp_star_cross_section / (num_err_members_in_sets[PDF_index] * 1.)

            if (proton_or_lead_index == 1):
                Wm_N_expected = efficiency[0] * L_new_experiment * Wm_cross_section + efficiency[2] * L_new_experiment * Wm_star_cross_section
                Wp_N_expected = efficiency[1] * L_new_experiment * Wp_cross_section + efficiency[3] * L_new_experiment * Wp_star_cross_section

                Wm_err_expected = 1. / np.sqrt(Wm_N_expected) * (Wm_cross_section + Wm_star_cross_section)
                Wp_err_expected = 1. / np.sqrt(Wp_N_expected) * (Wp_cross_section + Wp_star_cross_section)


                Rcpm_err_expected = np.sqrt(1. / (Wm_cross_section + Wm_star_cross_section)**2 * Wp_err_expected**2 + \
                                    (Wp_cross_section + Wp_star_cross_section)**2 / (Wm_cross_section + Wm_star_cross_section)**4 * Wm_err_expected**2)
                
                print('err expected', Rcpm_err_expected, PDF_set)

            PDF_set = PDF_sets_raw[PDF_index] + '_' + proton_or_lead_array[proton_or_lead_index]

            if (which_cross_sections_included == 'both'):
                Rcpm = (Wp_cross_section + Wp_star_cross_section) / (Wm_cross_section + Wm_star_cross_section)

                if (PDF_sets_raw[PDF_index] == 'EPPS21nlo'):
                    Rcpm_pdf_err_up, Rcpm_pdf_err_down = compute_Rcpm_pdf_err_HESSIAN(PDF_index, PDF_set, Rcpm, 'both')
                else:
                    Rcpm, Rcpm_pdf_err_up, Rcpm_pdf_err_down = compute_Rcpm_pdf_err_MC(PDF_index, PDF_set, Rcpm, 'both')
                
                print('pdf err', Rcpm_pdf_err_up, Rcpm_pdf_err_down, PDF_set)

            elif (which_cross_sections_included == 'D'):
                Rcpm = Wp_cross_section / Wm_cross_section

                if (PDF_sets_raw[PDF_index] == 'EPPS21nlo'):
                    Rcpm_pdf_err_up, Rcpm_pdf_err_down = compute_Rcpm_pdf_err_HESSIAN(PDF_index, PDF_set, Rcpm, 'D')
                else:
                    Rcpm, Rcpm_pdf_err_up, Rcpm_pdf_err_down = compute_Rcpm_pdf_err_MC(PDF_index, PDF_set, Rcpm, 'D')

            else:
                Rcpm = Wp_star_cross_section / Wm_star_cross_section

                if (PDF_sets_raw[PDF_index] == 'EPPS21nlo'):
                    Rcpm_pdf_err_up, Rcpm_pdf_err_down = compute_Rcpm_pdf_err_HESSIAN(PDF_index, PDF_set, Rcpm, 'Dstar')
                else:
                    Rcpm, Rcpm_pdf_err_up, Rcpm_pdf_err_down = compute_Rcpm_pdf_err_MC(PDF_index, PDF_set, Rcpm, 'Dstar')

            #print('Rcpm (' + which_cross_sections_included + ') with ' + PDF_set + ': ' + str(round(Rcpm, 5)) + \
            #        '(+' + str(round(Rcpm_pdf_err_up, 5)) + '-' + str(round(Rcpm_pdf_err_down, 5)) + ').')

            ax2.plot([Rcpm, Rcpm],
                    [y_vals[PDF_index] + p_Pb_shift[proton_or_lead_index] - 0.2, y_vals[PDF_index] + p_Pb_shift[proton_or_lead_index] + 0.2],
                    color='black', zorder=6, solid_capstyle='butt')

            pdf_err = patches.Rectangle((Rcpm - Rcpm_pdf_err_down, 2 - PDF_index - 0.2 + p_Pb_shift[proton_or_lead_index]),
                        Rcpm_pdf_err_down + Rcpm_pdf_err_up, 0.4, facecolor=pdf_err_color, zorder=5)
            ax2.add_patch(pdf_err)

            if (proton_or_lead_index == 1):
                expected_err = patches.Rectangle((Rcpm - Rcpm_err_expected, 2 - PDF_index - 0.2 + p_Pb_shift[proton_or_lead_index]),
                                    2 * Rcpm_err_expected, 0.4, facecolor='lightgray', zorder=3)
                ax2.add_patch(expected_err)

                ax2.text(Rcpm - Rcpm_err_expected - 0.07, y_vals[PDF_index] - 0.28, r'$p$Pb', fontsize=font_size)
            else:
                ax2.text(Rcpm - Rcpm_pdf_err_down - 0.05, y_vals[PDF_index] + 0.16, r'$pp$', fontsize=font_size)

            if (PDF_set == 'NNPDF30_nlo_as_01180_pp'):
                ax2.text(0.5, y_vals[0], r'nNNPDF3.0NLO', fontsize=font_size)
            elif (PDF_set == 'EPPS21nlo_pp'):
                ax2.text(0.5, y_vals[1], r'EPPS21NLO', fontsize=font_size)

    atlas_val = 0.971

    atlas_syst_up = 0.011
    atlas_syst_down = 0.011

    atlas_stat_up = 0.006
    atlas_stat_down = 0.006

    plt.xlim(0.47, 1.12)
    ax2.set_ylim(0, 3.5)
    ax1.set_ylim(0.2, 4)

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

    ax2.set_xlabel(r'$R_c^\pm (D^\pm, D^{*\pm})$', fontsize=axis_label_font_size)

    info_xval_1 = 0.5
    info_yval_1 = 3
    info_yval_2 = 2
    info_yval_3 = 1

    """
    if (which_cross_sections_included == 'both'):
        ax.text(info_xval_1, info_yval_1, r'$D$, $D^*$', fontsize=font_size)
    elif (which_cross_sections_included == 'D'):
        ax.text(info_xval_1, info_yval_1, r'$D$', fontsize=font_size)
    else:
        ax.text(info_xval_1, info_yval_1, r'$D^*$', fontsize=font_size)
    """

    ax1.text(info_xval_1, info_yval_1, r'$\sqrt{s} = 8.5$ TeV', fontsize=font_size)
    ax1.text(info_xval_1, info_yval_2, frag_set_text, fontsize=font_size)
    ax1.text(info_xval_1, info_yval_3, 'OS-SS', fontsize=font_size)

    pdf_err = patches.Rectangle((0, 0), 0, 0, facecolor=pdf_err_color, zorder=5)
    expected_err = patches.Rectangle((0, 0), 0, 0, facecolor='lightgray', zorder=3)
    ax1.add_patch(pdf_err)
    ax1.add_patch(expected_err)

    ax1.legend([pdf_err, expected_err], [r"PDF error (68\% C.L.)", "Expected measurement\nerror"], loc='upper right',
                framealpha=1, fontsize=legend_fontsize)

    plt.tight_layout()

    plt.subplots_adjust(left=0.05)

    plt.savefig(plots_directory + 'Rcpm/Rcpm_' + which_cross_sections_included + '_pPb.pdf')
    plt.show()


Rcpm('both', PDF_sets_raw, ['N', 'E'])