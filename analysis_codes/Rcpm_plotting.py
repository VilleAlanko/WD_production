import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl

import analysis_functions

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
num_err_members_in_sets = [58, 58, 100]

num_etac_bins = 11

main_vals_directory = '/home/alankovh/Documents/WD_production/output/'
plots_directory = '/home/alankovh/Documents/WD_production/plots/13 TeV/'
reweighting_input_directory = '/home/alankovh/Documents/WD_production/reweighting/input/'

marker_color = 'black'
theory_edge_colors = ['red', 'blue']
markers = ['d', 'v', 'o']
theory_labels = ['CT18ANLO', 'MSHT20NLO', 'NNPDF4.0NLO (pch)']

pdf_centrals = [[[np.zeros((284, num_etac_bins)) for _ in range(5)] for _ in range(len(PDF_sets))] for _ in range(4)]


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
        scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus = analysis_functions.compute_general_3D_vals_NLO(PDF_set, num_err_members_in_set,
                                                                    process_here, True, True, True, 'frag_main_scale', z_def, fragmentation_set)
        Wm_cross_section = sum(sum(sum(scales_vals[0])))
        Wm_scales_dd = sum(sum(sum(scales_vals[1])))
        Wm_scales_uu = sum(sum(sum(scales_vals[2])))
        Wm_MCerr = sum(sum(sum(scales_MCerrs[0])))

        process_here = "W-Dstar+"
        scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus = analysis_functions.compute_general_3D_vals_NLO(PDF_set, num_err_members_in_set,
                                                                    process_here, True, True, True, 'frag_main_scale', z_def, fragmentation_set)
        Wm_star_cross_section = sum(sum(sum(scales_vals[0])))
        Wm_star_scales_dd = sum(sum(sum(scales_vals[1])))
        Wm_star_scales_uu = sum(sum(sum(scales_vals[2])))
        Wm_star_MCerr = sum(sum(sum(scales_MCerrs[0])))

        process_here = "W+D-"
        scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus = analysis_functions.compute_general_3D_vals_NLO(PDF_set, num_err_members_in_set,
                                                                    process_here, True, True, True, 'frag_main_scale', z_def, fragmentation_set)
        Wp_cross_section = sum(sum(sum(scales_vals[0])))
        Wp_scales_dd = sum(sum(sum(scales_vals[1])))
        Wp_scales_uu = sum(sum(sum(scales_vals[2])))
        Wp_MCerr = sum(sum(sum(scales_MCerrs[0])))

        process_here = "W+Dstar-"
        scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus = analysis_functions.compute_general_3D_vals_NLO(PDF_set, num_err_members_in_set,
                                                                    process_here, True, True, True, 'frag_main_scale', z_def, fragmentation_set)
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
                Rcpm_pdf_err_up, Rcpm_pdf_err_down = analysis_functions.compute_Rcpm_pdf_err_HESSIAN(PDF_index, Rcpm, 'both')
            else:
                Rcpm, Rcpm_pdf_err_up, Rcpm_pdf_err_down = analysis_functions.compute_Rcpm_pdf_err_MC(PDF_index, Rcpm, 'both')

        elif (which_cross_sections_included == 'D'):
            Rcpm = Wp_cross_section / Wm_cross_section

            Rcpm_scale_var_dd = (Wp_cross_section + Wp_scales_dd) / (Wm_cross_section + Wm_scales_dd)
            Rcpm_scale_var_uu = (Wp_cross_section + Wp_scales_uu) / (Wm_cross_section + Wm_scales_uu)

            Rcpm_scale_var_down = Rcpm - min(Rcpm_scale_var_dd, Rcpm_scale_var_uu)
            Rcpm_scale_var_up = max(Rcpm_scale_var_dd, Rcpm_scale_var_uu) - Rcpm

            Rcpm_MCerr_up = (Wp_cross_section + Wp_MCerr) / (Wm_cross_section - Wm_MCerr) - Rcpm
            Rcpm_MCerr_down = Rcpm - (Wp_cross_section - Wp_MCerr) / (Wm_cross_section + Wm_MCerr)

            if (PDF_set == 'CT18NLO' or PDF_set == 'CT18ANLO' or PDF_set == 'MSHT20nlo_as118'):
                Rcpm_pdf_err_up, Rcpm_pdf_err_down = analysis_functions.compute_Rcpm_pdf_err_HESSIAN(PDF_index, Rcpm, 'D')
            else:
                Rcpm, Rcpm_pdf_err_up, Rcpm_pdf_err_down = analysis_functions.compute_Rcpm_pdf_err_MC(PDF_index, Rcpm, 'D')

        else:
            Rcpm = Wp_star_cross_section / Wm_star_cross_section

            Rcpm_scale_var_dd = (Wp_star_cross_section + Wp_star_scales_dd) / (Wm_star_cross_section + Wm_star_scales_dd)
            Rcpm_scale_var_uu = (Wp_star_cross_section + Wp_star_scales_uu) / (Wm_star_cross_section + Wm_star_scales_uu)

            Rcpm_scale_var_down = Rcpm - min(Rcpm_scale_var_dd, Rcpm_scale_var_uu)
            Rcpm_scale_var_up = max(Rcpm_scale_var_dd, Rcpm_scale_var_uu) - Rcpm

            Rcpm_MCerr_up = (Wp_star_cross_section + Wp_star_MCerr) / (Wm_star_cross_section - Wm_star_MCerr) - Rcpm
            Rcpm_MCerr_down = Rcpm - (Wp_star_cross_section - Wp_star_MCerr) / (Wm_star_cross_section + Wm_star_MCerr)

            if (PDF_set == 'CT18NLO' or PDF_set == 'CT18ANLO' or PDF_set == 'MSHT20nlo_as118'):
                Rcpm_pdf_err_up, Rcpm_pdf_err_down = analysis_functions.compute_Rcpm_pdf_err_HESSIAN(PDF_index, Rcpm, 'Dstar')
            else:
                Rcpm, Rcpm_pdf_err_up, Rcpm_pdf_err_down = analysis_functions.compute_Rcpm_pdf_err_MC(PDF_index, Rcpm, 'Dstar')

        Rcpm_error_up = np.sqrt(Rcpm_scale_var_up**2 + Rcpm_MCerr_up**2 + Rcpm_pdf_err_up**2)
        Rcpm_error_down = np.sqrt(Rcpm_scale_var_down**2 + Rcpm_MCerr_down**2 + Rcpm_pdf_err_down**2)

        print('Rcpm (' + which_cross_sections_included + ') with ' + PDF_set + ': ' + str(round(Rcpm, 5)) + \
                '(+' + str(round(Rcpm_error_up, 5)) + '-' + str(round(Rcpm_error_down, 5)) + ').')
    
        ax.plot(Rcpm, 4 - y_vals[PDF_index], marker=markers[PDF_index], color=marker_color,
                linestyle='none', label=theory_labels[PDF_index], zorder=6)

        pdf_err = patches.Rectangle((Rcpm - Rcpm_pdf_err_down, 3 - PDF_index - 0.2),
                    Rcpm_pdf_err_down + Rcpm_pdf_err_up, 0.4, facecolor=pdf_err_color, zorder=5)
        ax.add_patch(pdf_err)

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

    plt.tick_params(direction='in', top=True, right=True)
    plt.minorticks_on()
    plt.tick_params(which='both', direction='in', top=True, right=True)

    ax.text(info_xval_1, info_yval_1, r'$\sqrt{s} = 13$ TeV', fontsize=font_size)
    ax.text(info_xval_1, info_yval_2, frag_set_text, fontsize=font_size)
    ax.text(info_xval_1, info_yval_3, 'OS-SS', fontsize=font_size)

    legend1 = ax.legend(fontsize=legend_fontsize, bbox_to_anchor=(0.99, 0.88), loc='center right')
    legend2 = ax.legend([pdf_err, total_err], ["PDF error (68\% C.L.)", "Total theory error"], loc='center right',
                        bbox_to_anchor=(0.99, 0.69), framealpha=1, fontsize=legend_fontsize)
    ax.legend([atlas_stat_err, atlas_total_err], ['ATLAS stat. error', 'ATLAS tot. error'],
                        loc='center left', bbox_to_anchor=(0.01, 0.69), framealpha=1, fontsize=legend_fontsize)
    ax.add_artist(legend1)
    ax.add_artist(legend2)

    ax.tick_params(axis='both', which='major', labelsize=axis_font_size)
    ax.set_yticklabels([])

    plt.tight_layout()

    plt.savefig(plots_directory + 'Rcpm/Rcpm_' + which_cross_sections_included + '.pdf')
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

    print(places_inside_bins_log_scale[1, :])

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
            scales_vals, scales_MCerrs, pdf_err_plus, pdf_err_minus = analysis_functions.compute_general_3D_vals_NLO(
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
                        HISTO_Rcpm_pdf_err_up[PDF_index], HISTO_Rcpm_pdf_err_down[PDF_index] = analysis_functions.compute_Rcpm_pdf_err_pTD_HESSIAN(PDF_index,
                                                                                HISTO_Rcpm_central[PDF_index], 'both')
                    else:
                        HISTO_Rcpm_pdf_err_up[PDF_index], HISTO_Rcpm_pdf_err_down[PDF_index] = analysis_functions.compute_Rcpm_pdf_err_eta_lept_HESSIAN(PDF_index,
                                                                                HISTO_Rcpm_central[PDF_index], 'both')
                        
                    if (PDF_set == 'MSHT20nlo_as118'):
                        print('PDF error up for ' + PDF_set + ': ' + str(HISTO_Rcpm_pdf_err_up[PDF_index]))
                        print('PDF error down for ' + PDF_set + ': ' + str(HISTO_Rcpm_pdf_err_down[PDF_index]))
                else:
                    if (kinematic_variable == 'pTD'):
                        HISTO_Rcpm_central[PDF_index], HISTO_Rcpm_pdf_err_up[PDF_index], HISTO_Rcpm_pdf_err_down[PDF_index] = analysis_functions.compute_Rcpm_pdf_err_pTD_MC(PDF_index,
                                                                                HISTO_Rcpm_central[PDF_index], 'both')
                    else:
                        HISTO_Rcpm_central[PDF_index], HISTO_Rcpm_pdf_err_up[PDF_index], HISTO_Rcpm_pdf_err_down[PDF_index] = analysis_functions.compute_Rcpm_pdf_err_eta_lept_MC(PDF_index,
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
                        HISTO_Rcpm_pdf_err_up[PDF_index], HISTO_Rcpm_pdf_err_down[PDF_index] = analysis_functions.compute_Rcpm_pdf_err_pTD_HESSIAN(PDF_index,
                                                                                HISTO_Rcpm_central[PDF_index], 'D')
                    else:
                        HISTO_Rcpm_pdf_err_up[PDF_index], HISTO_Rcpm_pdf_err_down[PDF_index] = analysis_functions.compute_Rcpm_pdf_err_eta_lept_HESSIAN(PDF_index,
                                                                                HISTO_Rcpm_central[PDF_index], 'D')
                else:
                    if (kinematic_variable == 'pTD'):
                        HISTO_Rcpm_central[PDF_index], HISTO_Rcpm_pdf_err_up[PDF_index], HISTO_Rcpm_pdf_err_down[PDF_index] = analysis_functions.compute_Rcpm_pdf_err_pTD_MC(PDF_index,
                                                                                HISTO_Rcpm_central[PDF_index], 'D')
                    else:
                        HISTO_Rcpm_central[PDF_index], HISTO_Rcpm_pdf_err_up[PDF_index], HISTO_Rcpm_pdf_err_down[PDF_index] = analysis_functions.analysis_functions.compute_Rcpm_pdf_err_eta_lept_MC(PDF_index,
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
                        HISTO_Rcpm_pdf_err_up[PDF_index], HISTO_Rcpm_pdf_err_down[PDF_index] = analysis_functions.compute_Rcpm_pdf_err_pTD_HESSIAN(PDF_index,
                                                                                HISTO_Rcpm_central[PDF_index], 'Dstar')
                    else:
                        HISTO_Rcpm_pdf_err_up[PDF_index], HISTO_Rcpm_pdf_err_down[PDF_index] = analysis_functions.compute_Rcpm_pdf_err_eta_lept_HESSIAN(PDF_index,
                                                                                HISTO_Rcpm_central[PDF_index], 'Dstar')
                else:
                    if (kinematic_variable == 'pTD'):
                        HISTO_Rcpm_central[PDF_index], HISTO_Rcpm_pdf_err_up[PDF_index], HISTO_Rcpm_pdf_err_down[PDF_index] = analysis_functions.compute_Rcpm_pdf_err_pTD_MC(PDF_index,
                                                                                HISTO_Rcpm_central[PDF_index], 'Dstar')
                    else:
                        HISTO_Rcpm_central[PDF_index], HISTO_Rcpm_pdf_err_up[PDF_index], HISTO_Rcpm_pdf_err_down[PDF_index] = analysis_functions.compute_Rcpm_pdf_err_eta_lept_MC(PDF_index,
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


#Rcpm('both')
Rcpm_bin_integrated('pTD', 'both', True, ['CT18ANLO', 'MSHT20nlo_as118', 'NNPDF40_nlo_pch_as_01180'])