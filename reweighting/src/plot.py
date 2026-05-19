import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.integrate import quad
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, VPacker
import matplotlib as mpl
import matplotlib.ticker as ticker
from termcolor import colored

mpl.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.weight": "bold"
})
mpl.rcParams['text.latex.preamble'] = r'''
\usepackage{amsmath}
\usepackage{xcolor}
'''

log10x_min = -4
log10x_max = 0.
log10x_step = 0.001
Nx = int((log10x_max - log10x_min) / log10x_step + 1)
x = np.zeros(Nx)

for i in range(Nx):
    x[i] = 10**(log10x_min + i * log10x_step)

font_size = 17
axis_label_font_size = 18
axis_font_size = 14
legend_fontsize = 14

PDF_set_labels = ['CT18ANLO', 'MSHT20NLO', 'NNPDF4.0NLO (pch)']


def ratio_HESSIAN(PDF_set, which_cross_sections_included, flavor_and_valence_list, particle_names, tolerances, labels, colors, lower_bounds, upper_bounds):
    for i in range(len(flavor_and_valence_list)):
        plt.figure(figsize=(6, 4))

        for tolerance_index in range(len(tolerances)):
            tolerance = tolerances[tolerance_index]

            NEW_PDF = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + \
                                '/flavor_' + str(flavor_and_valence_list[i][0]) + '_t_' + tolerance + '_best.txt', delimiter=',', skiprows=1, max_rows=1)
            OLD_PDF = np.loadtxt('/home/alankovh/Documents/PDFs/Values/' + PDF_set + '/flavor_id_' + str(flavor_and_valence_list[i][0]) + 
                                    '/best_vals_valence_' + str(flavor_and_valence_list[i][1]) + '.txt', delimiter=',')
            
            if (flavor_and_valence_list[i][1] == 1):
                NEW_PDF -= np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + \
                                '/flavor_' + str(-flavor_and_valence_list[i][0]) + '_t_' + tolerance + '_best.txt', delimiter=',', skiprows=1, max_rows=1)

            OLD_ERR_MINUS = np.loadtxt('/home/alankovh/Documents/PDFs/Values/' + PDF_set + '/flavor_id_' + str(flavor_and_valence_list[i][0]) + 
                                    '/err_minus_valence_' + str(flavor_and_valence_list[i][1]) + '.txt', delimiter=',')
            OLD_ERR_PLUS = np.loadtxt('/home/alankovh/Documents/PDFs/Values/' + PDF_set + '/flavor_id_' + str(flavor_and_valence_list[i][0]) + 
                                    '/err_plus_valence_' + str(flavor_and_valence_list[i][1]) + '.txt', delimiter=',')

            print(tolerance_index)
            print(flavor_and_valence_list[i])
            print(sum(NEW_PDF))
            print(sum(OLD_PDF))

            plt.plot(x, NEW_PDF / OLD_PDF, zorder=5, color=colors[tolerance_index], label=labels[tolerance_index])

        plt.fill_between(x, (OLD_PDF + OLD_ERR_PLUS) / OLD_PDF, (OLD_PDF - OLD_ERR_MINUS) / OLD_PDF, zorder=3, color='lightgray', alpha=0.5)

        plt.plot([0, 2], [1, 1], color='black', zorder=4)

        plt.xlim([10**(-3), 10**(-1)])
        plt.ylim([lower_bounds[i], upper_bounds[i]])
        plt.xscale('log')

        plt.xlabel(r'$x$')
        plt.ylabel('Reweighted / Original', fontsize=axis_label_font_size)

        plt.legend(fontsize=15, framealpha=1, loc='lower left')
        
        plt.title(particle_names[i], fontsize=axis_label_font_size)
        plt.tight_layout()
        plt.savefig('plots/ratio/' + PDF_set + '_' + which_cross_sections_included + '_id_' + str(flavor_and_valence_list[i][0]) + 
                    '_valence_' + str(flavor_and_valence_list[i][1]) + '.pdf')

        plt.close()


def ratio_MC(PDF_set, which_cross_sections_included, flavor_and_valence_list, particle_names, lower_bounds, upper_bounds):
    for i in range(len(flavor_and_valence_list)):
        plt.figure(figsize=(6, 4))

        x = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + \
                    '/flavor_' + str(flavors[i]) + '_best.txt', delimiter=',', max_rows=1)
        NEW_PDF = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + \
                            '/flavor_' + str(flavor_and_valence_list[i][0]) + '_best.txt', delimiter=',', skiprows=1, max_rows=1)

        OLD_PDF = np.loadtxt('/home/alankovh/Documents/PDFs/Values/' + PDF_set + '/flavor_id_' + str(flavor_and_valence_list[i][0]) + 
                                '/best_vals_valence_' + str(flavor_and_valence_list[i][1]) + '.txt', delimiter=',')
        
        if (flavor_and_valence_list[i][1] == 1):
            NEW_PDF -= np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + \
                            '/flavor_' + str(-flavor_and_valence_list[i][0]) + '_best.txt', delimiter=',', skiprows=1, max_rows=1)
        
        #print(flavor_and_valence_list[i])
        #print("NEW:", NEW_PDF[:5])
        #print("OLD:", OLD_PDF[:5])

        print(flavor_and_valence_list[i])
        print(sum(NEW_PDF))
        print(sum(OLD_PDF))

        OLD_ERR_MINUS = np.loadtxt('/home/alankovh/Documents/PDFs/Values/' + PDF_set + '/flavor_id_' + str(flavor_and_valence_list[i][0]) + 
                                '/err_minus_valence_' + str(flavor_and_valence_list[i][1]) + '.txt', delimiter=',')
        OLD_ERR_PLUS = np.loadtxt('/home/alankovh/Documents/PDFs/Values/' + PDF_set + '/flavor_id_' + str(flavor_and_valence_list[i][0]) + 
                                '/err_plus_valence_' + str(flavor_and_valence_list[i][1]) + '.txt', delimiter=',')

        log10x_max = 0.
        log10x_step = 0.001
        N = int((log10x_max - log10x_min) / log10x_step) + 1

        #print(sum(NEW_PDF / OLD_PDF) / len(NEW_PDF))
        plt.plot(x, NEW_PDF / OLD_PDF, zorder=5, color='red')

        plt.fill_between(x, (OLD_PDF + OLD_ERR_PLUS) / OLD_PDF, (OLD_PDF - OLD_ERR_MINUS) / OLD_PDF, zorder=3, color='lightgray', alpha=0.5)

        plt.plot([0, 2], [1, 1], color='black', zorder=4)

        plt.xlim([10**(-3), 10**(-1)])
        plt.ylim([lower_bounds[i], upper_bounds[i]])
        plt.xscale('log')

        plt.xlabel(r'$x$')
        plt.ylabel('Reweighted / Original', fontsize=axis_label_font_size)

        #plt.legend(fontsize=15, framealpha=1, loc='lower left')
        
        plt.title(particle_names[i], fontsize=axis_label_font_size)
        plt.tight_layout()
        plt.savefig('plots/ratio/' + PDF_set + '_' + which_cross_sections_included + '_id_' + str(flavor_and_valence_list[i][0]) + 
                    '_valence_' + str(flavor_and_valence_list[i][1]) + '.pdf')

        plt.close()


def ratio_to_other_PDF(PDF_set, comparison_set, flavors, which_cross_sections_included):
    for flavor in flavors:
        plt.figure(figsize=(6, 4))

        x = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + \
                    '/flavor_' + str(flavor) + '_best.txt', delimiter=',', max_rows=1)
        NEW_PDF = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + \
                            '/flavor_' + str(flavor) + '_best.txt', delimiter=',', skiprows=1, max_rows=1)
        OLD_PDF = np.loadtxt('output/old_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + \
                            '/flavor_' + str(flavor) + '.txt', delimiter=',', skiprows=1, max_rows=1)
        comparison_PDF_vals = np.loadtxt('/home/alankovh/Documents/PDFs/Values/' + comparison_set + '/flavor_id_' + str(flavor) + '/best_vals.txt', delimiter=',')
        lul = np.loadtxt('/home/alankovh/Documents/PDFs/Values/MSHT20nlo_as118/flavor_id_' + str(flavor) + '/best_vals.txt', delimiter=',')
        
        NEW_PDF_interpolation = CubicSpline(x, NEW_PDF)
        OLD_PDF_interpolation = CubicSpline(x, OLD_PDF)

        # These are used for the values of the comparison set
        log10xmin = -4.
        log10xmax = 0.
        log10xstep = 0.001

        N = int((log10xmax - log10xmin) / log10xstep) + 1

        log10X_vals = np.linspace(log10xmin, log10xmax, N)

        X = 10**log10X_vals
        print(X)

        ratio_OLD = OLD_PDF[1000:] / comparison_PDF_vals
        ratio_NEW = NEW_PDF[1000:] / comparison_PDF_vals

        plt.plot(X, ratio_OLD, label='OLD')
        plt.plot(X, ratio_NEW, label='NEW')
        plt.plot([0, 1], [1, 1], color='black', zorder=0)
        #plt.plot(x_fine, lul / comparison_PDF_vals, color='red')

        plt.xlabel(r'$x$', fontsize=axis_label_font_size)
        plt.ylabel(PDF_set + ' / ' + comparison_set, fontsize=axis_label_font_size)

        plt.xlim([10**(log10xmin), 10**(log10xmax)])
        plt.ylim(0.9, 1.1)
        plt.xscale('log', base=10)

        plt.legend(fontsize=legend_fontsize, loc='upper left')

        plt.title(str(flavor), fontsize=axis_label_font_size)
        plt.tight_layout()
        plt.savefig('plots/ratio_to_other_PDF/' + PDF_set + '_' + comparison_set + '_' + which_cross_sections_included + '_' + str(flavor) + '.pdf')

        plt.close()



def ratio_of_ratio(flavor):
    plt.figure(figsize=(6, 4))

    x = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/flavor_' + str(flavor) + '_best.txt', delimiter=',', max_rows=1)
    NEW_PDF = np.loadtxt('../output/new_PDF_vals/' + PDF_set + '/flavor_' + str(flavor) + '_best.txt', delimiter=',', skiprows=1, max_rows=1)
    NEW_PDF_bar = np.loadtxt('../output/new_PDF_vals/' + PDF_set + '/flavor_' + str(-flavor) + '_best.txt', delimiter=',', skiprows=1, max_rows=1)
    OLD_PDF = np.loadtxt('output/old_PDF_vals/' + PDF_set + '/flavor_' + str(flavor) + '.txt', delimiter=',', skiprows=1, max_rows=1)
    OLD_PDF_bar = np.loadtxt('output/old_PDF_vals/' + PDF_set + '/flavor_' + str(-flavor) + '.txt', delimiter=',', skiprows=1, max_rows=1)

    log10x_max = 0
    log10x_step = 0.001
    N = int((log10x_max - log10x_min) / log10x_step) + 1

    plt.plot(x, (NEW_PDF / NEW_PDF_bar) / (OLD_PDF / OLD_PDF_bar), zorder=5, color='red')

    plt.plot([0, 2], [1, 1], color='black', zorder=0)

    plt.xlim([10**(-3), 10**(0)])
    plt.ylim([0.9, 1.1])
    plt.xscale('log')

    plt.ylabel(r'$(\text{new} / \overline{\text{new}}) / (\text{old} / \overline{\text{old}})$', fontsize=axis_label_font_size)
    
    plt.title(str(flavor), fontsize=font_size)

    plt.tight_layout()

    plt.savefig('plots/ratio_of_ratio/' + PDF_set + '_' + str(flavor) + '.pdf')

    plt.show()
    plt.close()


def absolute(tolerance, which_cross_sections_included, flavors, valences):
    for i in range(len(flavors)):
        flavor = flavors[i]
        valence = valences[i]
        plt.figure(figsize=(6, 4))

        x = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + '/flavor_' + str(flavor) + '_t_' + tolerance + '_best.txt', delimiter=',', max_rows=1)
        NEW_PDF = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + '/flavor_' + str(flavor) + '_t_' + tolerance + '_best.txt', delimiter=',', skiprows=1, max_rows=1)
        OLD_PDF = np.loadtxt('/home/alankovh/Documents/PDFs/Values/' + PDF_set + '/flavor_id_' + str(flavor) + '/best_vals_valence_' + str(valence) + '.txt', delimiter=',')

        if (valence == 1):
            NEW_PDF -= np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + '/flavor_' + str(-flavor) + '_t_' + tolerance + '_best.txt', delimiter=',', skiprows=1, max_rows=1)

        log10x_min = -4
        log10x_max = 0
        log10x_step = 0.001
        N = int((log10x_max - log10x_min) / log10x_step) + 1

        plt.plot(x, NEW_PDF, zorder=1, color='red', label='new')
        plt.plot(x, OLD_PDF, zorder=0, color='blue', label='old')

        plt.plot([0, 2], [0, 0], color='black', zorder=0)

        plt.xlim([10**(-3), 10**(-1)])
        plt.xscale('log')

        plt.ylabel(r'$x\cdot $ PDF', fontsize=axis_label_font_size)
        plt.legend()
        
        plt.title(str(flavor), fontsize=font_size)

        plt.tight_layout()

        plt.savefig('plots/absolute/' + PDF_set + '_flavor_' + str(flavor) + '_valence_' + str(valence) + '.pdf')

        plt.close()


def Rcpm_OLD_and_NEW_and_DATA(kinematic_quantity, which_cross_sections_included, tolerance):
    markers = ['s', 'd', 'v']
    marker_color = 'black'
    theory_labels = ['CT18ANLO', 'MSHT20NLO', 'NNPDF40NLO4.0 (pch)']
    PDF_sets = ['CT18ANLO', 'MSHT20nlo_as118', 'NNPDF40_nlo_pch_as_01180']
    old_color = 'cornflowerblue'
    new_color = 'salmon'

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1], 'hspace': 0}, figsize=(6, 6))

    #--------------------------------------------------------------------------------------------------------------------------------------#
    #                                                                IMPORTING VALUES                                                      #
    #--------------------------------------------------------------------------------------------------------------------------------------#

    DATA = np.loadtxt('input/experimental_values/' + kinematic_quantity + '_' + which_cross_sections_included + '.txt')
    DATA_errors = np.loadtxt('input/exp_errors/' + kinematic_quantity + '_' + which_cross_sections_included + '.txt')

    for PDF_index in range(len(PDF_sets)):
        PDF_set = PDF_sets[PDF_index]
        if (PDF_set == 'CT18ANLO' or PDF_set == 'MSHT20nlo_as118'):
            OLD = np.loadtxt('input/theory_values/HESSIAN/best/' + kinematic_quantity + '_' + which_cross_sections_included + '_' + PDF_set + '_best.txt')
            OLD_err_plus = np.loadtxt('input/theory_values/HESSIAN/errors/' + kinematic_quantity + '_' + \
                            which_cross_sections_included + '_' + PDF_set + '_plus.txt')
            OLD_err_minus = np.loadtxt('input/theory_values/HESSIAN/errors/' + kinematic_quantity + '_' + \
                            which_cross_sections_included + '_' + PDF_set + '_minus.txt')
        else:
            OLD = np.loadtxt('input/theory_values/MC/best/' + kinematic_quantity + '_' + which_cross_sections_included + '_' + PDF_set + '_best.txt')
            OLD_err_plus = np.loadtxt('input/theory_values/MC/errors/' + kinematic_quantity + '_' + \
                            which_cross_sections_included + '_' + PDF_set + '_plus.txt')
            OLD_err_minus = np.loadtxt('input/theory_values/MC/errors/' + kinematic_quantity + '_' + \
                            which_cross_sections_included + '_' + PDF_set + '_minus.txt')

        if (PDF_set != "NNPDF40_nlo_pch_as_01180"):
            NEW = np.loadtxt('output/new_Rcpm_' + kinematic_quantity + '/' + PDF_set + '/' + which_cross_sections_included + '/t_' + tolerance + '_best.txt')
            NEW_err_up = np.loadtxt('output/new_Rcpm_' + kinematic_quantity + '/' + PDF_set + '/' + which_cross_sections_included + '/t_' + tolerance + '_error_up.txt')
            NEW_err_down = np.loadtxt('output/new_Rcpm_' + kinematic_quantity + '/' + PDF_set + '/' + which_cross_sections_included + '/t_' + tolerance + '_error_down.txt')
        else:
            NEW = np.loadtxt('output/new_Rcpm_' + kinematic_quantity + '/' + PDF_set + '/' + which_cross_sections_included + '/best.txt')
            NEW_err_up = np.loadtxt('output/new_Rcpm_' + kinematic_quantity + '/' + PDF_set + '/' + which_cross_sections_included + '/error_up.txt')
            NEW_err_down = np.loadtxt('output/new_Rcpm_' + kinematic_quantity + '/' + PDF_set + '/' + which_cross_sections_included + '/error_down.txt')

        #--------------------------------------------------------------------------------------------------------------------------------------#
        #                                                           AX1 PLOTTING (THEORY)                                                      #
        #--------------------------------------------------------------------------------------------------------------------------------------#

        if (kinematic_quantity == 'eta_lept'):
            bin_edges = np.array([0., 0.5, 1., 1.5, 2., 2.5])
            bar_width = 0.11
            plt.xlabel(r'$|\eta(l)|$', fontsize=axis_label_font_size)
            xmin = bin_edges[0:-1]
            xmax = bin_edges[1:]
            bin_midpoints = (xmin + xmax) / 2
            theory_val_places = [bin_midpoints - 1. / 4. * 0.5, bin_midpoints, bin_midpoints + 1. / 4. * 0.5]

            
        elif (kinematic_quantity == 'pTD'):
            bin_edges = np.array([8., 12., 20., 40., 80., 150.])
            plt.xlabel(r'$p_T (D)$ [GeV]', fontsize=axis_label_font_size)

            plt.xscale('log', base=10)

            bin_widths = np.diff(bin_edges)
            bar_width_over_bin_width = 1. / 5.
            bar_width = bar_width_over_bin_width * bin_widths

            theory_val_places = [bin_edges[:-1] * (bin_edges[1:] / bin_edges[:-1])**(1. / 4.), 
                                bin_edges[:-1] * (bin_edges[1:] / bin_edges[:-1])**(2. / 4.),
                                bin_edges[:-1] * (bin_edges[1:] / bin_edges[:-1])**(3. / 4.)]
        
        ax1.plot(theory_val_places[PDF_index], (NEW + OLD) / 2., marker=markers[PDF_index],
                                color=marker_color, markersize=5, linestyle='none',
                                label=theory_labels[PDF_index], zorder=5)
        
        before_reweighting = ax1.bar(theory_val_places[PDF_index] - bar_width / 4., OLD_err_minus + OLD_err_plus, width=bar_width / 2.,
                                bottom=OLD - OLD_err_minus, color=old_color, zorder=3)
        
        after_reweighting = ax1.bar(theory_val_places[PDF_index] + bar_width / 4., NEW_err_up + NEW_err_down, width=bar_width / 2.,
                                bottom=NEW - NEW_err_down, color=new_color, zorder=3)

        ax1.hlines(NEW, theory_val_places[PDF_index], theory_val_places[PDF_index] + bar_width / 2., color='darkred', zorder=4)
        ax1.hlines(OLD, theory_val_places[PDF_index] - bar_width / 2., theory_val_places[PDF_index], color='darkblue', zorder=4)

        #--------------------------------------------------------------------------------------------------------------------------------------#
        #                                                            AX2 PLOTTING (THEORY)                                                     #
        #--------------------------------------------------------------------------------------------------------------------------------------#

        ratio_up_err_OLD = (OLD + OLD_err_plus) / DATA - OLD / DATA
        ratio_down_err_OLD = OLD / DATA - (OLD - OLD_err_minus) / DATA

        
        ratio_up_err_NEW = (NEW + NEW_err_up) / DATA - NEW / DATA
        ratio_down_err_NEW = NEW / DATA - (NEW - NEW_err_down) / DATA

        ax2.plot(theory_val_places[PDF_index], ((NEW + OLD) / 2.) / DATA, marker=markers[PDF_index], color=marker_color,
                    markersize=5, linestyle='none', zorder=5)

        ax2.hlines(NEW / DATA, theory_val_places[PDF_index], theory_val_places[PDF_index] + bar_width / 2., color='darkred', zorder=4)
        ax2.hlines(OLD / DATA, theory_val_places[PDF_index], theory_val_places[PDF_index] - bar_width / 2., color='darkblue', zorder=4)

        ax2.bar(theory_val_places[PDF_index] - bar_width / 4., ratio_up_err_OLD + ratio_down_err_OLD, width=bar_width / 2.,
                            bottom=OLD / DATA - ratio_down_err_OLD, color=old_color, zorder=3)
        ax2.bar(theory_val_places[PDF_index] + bar_width / 4., ratio_up_err_NEW + ratio_down_err_NEW, width=bar_width / 2.,
                            bottom=NEW / DATA - ratio_down_err_NEW, color=new_color, zorder=3)

    #--------------------------------------------------------------------------------------------------------------------------------------#
    #                                                              AX1 PLOTTING (DATA)                                                     #
    #--------------------------------------------------------------------------------------------------------------------------------------#

    ax1.hlines(DATA, bin_edges[:-1], bin_edges[1:], color='black', label='ATLAS')
    ATLAS_uncertainty = ax1.bar((bin_edges[1:] + bin_edges[:-1]) / 2., 2. * DATA_errors, bottom=DATA - DATA_errors, width=bin_edges[1:] - bin_edges[:-1],
            color='lightgray', zorder=0)

    #--------------------------------------------------------------------------------------------------------------------------------------#
    #                                                              AX2 PLOTTING (DATA)                                                     #
    #--------------------------------------------------------------------------------------------------------------------------------------#

    ax2.bar((bin_edges[1:] + bin_edges[:-1]) / 2., (DATA + DATA_errors) / DATA - (DATA - DATA_errors) / DATA,
            bottom=(DATA - DATA_errors) / DATA, color='lightgray', width=bin_edges[1:] - bin_edges[:-1], zorder=0)
    ax2.plot([-1, 200], [1, 1], color='black', zorder=1)

    #--------------------------------------------------------------------------------------------------------------------------------------#
    #                                                        MAKING THE PLOT LOOK PRETTY :)                                                #
    #--------------------------------------------------------------------------------------------------------------------------------------#

    plt.xlim(bin_edges[0], bin_edges[-1])
    ax1.set_ylim(0.755, 1.15)
    ax2.set_ylim(0.83, 1.14)
    if (which_cross_sections_included == 'both'):
        ax1.set_ylabel(r'$R_c^\pm(D^\pm, D^{*\pm})$', fontsize=axis_label_font_size)
    elif (which_cross_sections_included == 'D'):
        ax1.set_ylabel(r'$R_c^\pm$($D^\pm)$', fontsize=axis_label_font_size)

    ax2.set_ylabel(r'$\frac{\text{Theory}}{\text{ATLAS}}$', fontsize=axis_label_font_size * 1.3)

    ax2.set_yticks([0.85, 0.9, 0.95, 1., 1.05, 1.1])

    for i in range(1, len(bin_edges) - 3):
        ax1.axvline(bin_edges[i], color='gray', linewidth=0.5, ymax=0.75, zorder=0)
    for i in range(len(bin_edges) - 3, len(bin_edges) - 1):
        ax1.axvline(bin_edges[i], color='gray', linewidth=0.5, ymax=0.62, zorder=0)
    for i in range(1, len(bin_edges) - 1):
        ax2.axvline(bin_edges[i], color='gray', linewidth=0.5, zorder=0)

    text_y1 = 1.11
    text_y2 = 1.075
    if (kinematic_quantity == 'pTD'):
        text_x = 9
    else:
        text_x = 0.1

    ax1.text(text_x, text_y1, r'$\sqrt{s} = 13$ TeV', fontsize=font_size)
    ax1.text(text_x, text_y2, 'KKKS08 OPAL', fontsize=font_size)

    if (kinematic_quantity == 'eta_lept'):
        plt.xticks(bin_edges, [f'{tick:.1f}' for tick in bin_edges])
    elif (kinematic_quantity == 'pTD'):
        plt.xticks(bin_edges, [f'{tick:.0f}' for tick in bin_edges])

    ax1.tick_params(axis='both', which='major', labelsize=axis_font_size)
    ax2.tick_params(axis='both', which='major', labelsize=axis_font_size)

    ax1.tick_params(direction='in', top=True, right=True)
    ax1.minorticks_on()
    ax1.tick_params(which='both', direction='in', top=True, right=True)
    ax2.minorticks_on()
    ax2.tick_params(which='both', direction='in', top=True, right=True)
    ax2.tick_params(direction='in', top=True, right=True)

    legend1 = ax1.legend(loc='upper right', fontsize=legend_fontsize, framealpha=1)
    ax1.legend([ATLAS_uncertainty, before_reweighting, after_reweighting],
                        ['ATLAS uncertainty', 'Before reweighting', 'After reweighting'],
                        loc='lower left', fontsize=legend_fontsize, framealpha=1)
    ax1.add_artist(legend1)

    plt.tight_layout()

    plt.savefig('plots/Rcpm/Rcpm_' + kinematic_quantity + '_' + which_cross_sections_included + '_t_' + tolerance + '.pdf')

    plt.show()
    plt.close()


def strangeness_asymmetry(PDF_set, which_cross_sections_included, plot_errors_flag, plot_mem_vals, tolerance):
    if (PDF_set == 'NNPDF40_nlo_pch_as_01180'):
        PDF_index = 2
        num_members = 100
    if (PDF_set == 'MSHT20nlo_as118'):
        PDF_index = 1
    
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1], 'hspace': 0}, figsize=(6, 6))

    x = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + '/flavor_3_t_' + tolerance + '_best.txt', delimiter=',', max_rows=1)
    NEW_s = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + '/flavor_3_t_' + tolerance + '_best.txt', delimiter=',', skiprows=1, max_rows=1)
    NEW_sbar = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + '/flavor_-3_t_' + tolerance + '_best.txt', delimiter=',', skiprows=1, max_rows=1)
    OLD_s = np.loadtxt('output/old_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + '/flavor_3.txt', delimiter=',', skiprows=1, max_rows=1)
    OLD_sbar = np.loadtxt('output/old_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + '/flavor_-3.txt', delimiter=',', skiprows=1, max_rows=1)
    if (plot_errors_flag):
        OLD_asymmetry_error_plus = np.loadtxt('output/strangeness_asymmetry_errors/' + PDF_set + '_OLD_plus.txt', delimiter=',')
        OLD_asymmetry_error_minus = np.loadtxt('output/strangeness_asymmetry_errors/' + PDF_set + '_OLD_minus.txt', delimiter=',')
        NEW_asymmetry_error_plus = np.loadtxt('output/strangeness_asymmetry_errors/' + PDF_set + '_' + \
                                                which_cross_sections_included + '_NEW_t_' + tolerance + '_plus.txt', delimiter=',')
        NEW_asymmetry_error_minus = np.loadtxt('output/strangeness_asymmetry_errors/' + PDF_set + '_' + \
                                                which_cross_sections_included + '_NEW_t_' + tolerance + '_minus.txt', delimiter=',')

    asymmetry_OLD = OLD_s - OLD_sbar
    asymmetry_NEW = NEW_s - NEW_sbar

    for x_index in range(len(asymmetry_OLD)):
        if (asymmetry_NEW[x_index] != asymmetry_NEW[x_index]):
            asymmetry_NEW[x_index] = 0.
        if (asymmetry_OLD[x_index] != asymmetry_OLD[x_index]):
            asymmetry_OLD[x_index] = 0.

    ax1.plot(x, asymmetry_OLD * 10**2, color='black', zorder=2, label=PDF_set_labels[PDF_index])
    ax1.plot(x, asymmetry_NEW * 10**2, color='red', zorder=2, label='reweighted')

    if (plot_errors_flag):
        PDF_uncertainty = ax1.fill_between(x, (asymmetry_OLD - OLD_asymmetry_error_minus) * 10**2, (asymmetry_OLD + OLD_asymmetry_error_plus) * 10**2,
                        color='lightgray', zorder=0, rasterized=True)
        reweighted_PDF_uncertainty = ax1.fill_between(x, (asymmetry_NEW - NEW_asymmetry_error_minus) * 10**2, (asymmetry_NEW + NEW_asymmetry_error_plus) * 10**2,
                        color='salmon', zorder=1, alpha=0.6, rasterized=True)

    # ax2.plot(x, (asymmetry_NEW - asymmetry_OLD) * 10**2, color='red', zorder=2)
    ax2.plot(x, asymmetry_NEW / asymmetry_OLD, color='red', zorder=2)

    if (plot_errors_flag):
        ax2.fill_between(x, OLD_asymmetry_error_plus * 10**2,
                        -OLD_asymmetry_error_minus * 10**2, color='lightgray', zorder=0)
        ax2.fill_between(x, (asymmetry_NEW + NEW_asymmetry_error_plus - asymmetry_OLD) * 10**2,
                        (asymmetry_NEW - NEW_asymmetry_error_minus - asymmetry_OLD) * 10**2, color='salmon', zorder=1, alpha=0.6)

    print(1. / (asymmetry_OLD + OLD_asymmetry_error_plus))
    #if (plot_errors_flag):
    #    ax2.fill_between(x, (asymmetry_OLD - OLD_asymmetry_error_minus) / asymmetry_OLD,
    #                    (asymmetry_OLD + OLD_asymmetry_error_plus) / asymmetry_OLD, color='lightgray', zorder=0, rasterized=True)
    #    ax2.fill_between(x, (asymmetry_NEW + NEW_asymmetry_error_plus) / asymmetry_OLD,
    #                    (asymmetry_NEW - NEW_asymmetry_error_minus) / asymmetry_OLD, color='salmon', zorder=1, alpha=0.6, rasterized=True)

    if (plot_mem_vals):
        mem_vals = np.loadtxt('output/strangeness_asymmetry_errors/' + PDF_set + '_OLD_mem_vals.txt', delimiter=',')

        for x_index in range(len(asymmetry_OLD)):
            if x_index % 100 == 0:
                for mem in range(num_members):
                    ax1.scatter(x[x_index], mem_vals[mem, x_index] * 10**2, color='black', s=0.1)

    plt.xscale('log')

    plt.xlim(10**(-5), 1)
    if (PDF_set == 'MSHT20nlo_as118'):
        ax1.set_ylim(-0.7, 2.3)
        #ax2.set_ylim(-1e-1 * 6., 1e-0 * 0.9)
        ax2.set_ylim(-0.5, 3.5)
    elif (PDF_set == 'NNPDF40_nlo_pch_as_01180'):
        ax1.set_ylim(-6, 15)
        ax2.set_ylim(-3, 7)
    

    ax1.plot([10**(-6), 2], [0, 0], color='black', zorder=1)
    ax2.plot([10**(-6), 2], [1, 1], color='black', zorder=1)

    plt.xlabel(r'$x$', fontsize=axis_label_font_size)
    ax1.set_ylabel(r'$x(s - \overline{s})$', fontsize=axis_label_font_size)
    #ax2.set_ylabel('Reweighted\n - Original', fontsize=axis_label_font_size)
    ax2.set_ylabel(r'$\frac{\text{Reweighted}}{\text{Original}}$', fontsize=axis_label_font_size * 1.3)

    ax1.tick_params(axis='both', which='major', labelsize=axis_font_size)
    ax2.tick_params(axis='both', which='major', labelsize=axis_font_size)

    text_x = 10**(-5) * 1.8
    text_y1 = 1.6

    ax1.text(text_x, text_y1, r'$\mu_\text{fact} = M_W$', fontsize=font_size)

    if (which_cross_sections_included == 'D'):
        ax1.text(text_x, 1.9, r'$R_c^\pm (D^\pm)$', fontsize=font_size)
    else:
        ax1.text(text_x, 1.9, r'$R_c^\pm (D^\pm, D^{*\pm})$', fontsize=font_size)

    ax1.text(0.01, 1.01, r'$\times 10^{-2}$',
        transform=ax1.transAxes,
        fontsize=13, va='bottom', ha='left')

    legend1 = ax1.legend(loc='lower left', fontsize=legend_fontsize, framealpha=1, bbox_to_anchor=(0., 0.5))
    legend2 = ax1.legend([PDF_uncertainty, reweighted_PDF_uncertainty],
                        ['PDF error (90% C.L.)', 'Reweighted'],
                        loc='lower left', fontsize=legend_fontsize, framealpha=1, bbox_to_anchor=(0., 0.3))
    ax1.add_artist(legend1)

    ax2.set_yticks([0, 1, 2, 3])

    # Configure ticks to appear on all sides
    ax1.tick_params(direction='in', top=True, right=True)
    ax2.tick_params(direction='in', top=True, right=True)

    # Add minor ticks
    ax1.minorticks_on()
    ax1.tick_params(which='both', direction='in', top=True, right=True)
    ax2.minorticks_on()
    ax2.tick_params(which='both', direction='in', top=True, right=True)

    plt.tight_layout()

    plt.savefig('plots/strange_asymmetry/strangeness_asymmetry_' + PDF_set + '_' + which_cross_sections_included + '.pdf', dpi=300)
    plt.show()


def strangeness_asymmetry_v2(PDF_set, which_cross_sections_included, plot_errors_flag, plot_mem_vals):
    if (PDF_set == 'NNPDF40_nlo_pch_as_01180'):
        PDF_index = 2
        num_members = 100
    if (PDF_set == 'MSHT20nlo_as118'):
        PDF_index = 1
    
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1], 'hspace': 0}, figsize=(6, 6))

    x = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + '/flavor_3_best.txt', delimiter=',', max_rows=1)
    NEW_s = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + '/flavor_3_best.txt', delimiter=',', skiprows=1, max_rows=1)
    NEW_sbar = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + '/flavor_-3_best.txt', delimiter=',', skiprows=1, max_rows=1)
    OLD_s = np.loadtxt('output/old_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + '/flavor_3.txt', delimiter=',', skiprows=1, max_rows=1)
    OLD_sbar = np.loadtxt('output/old_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + '/flavor_-3.txt', delimiter=',', skiprows=1, max_rows=1)
    if (plot_errors_flag):
        OLD_asymmetry_error_plus = np.loadtxt('output/strangeness_asymmetry_errors/' + PDF_set + '_OLD_plus.txt', delimiter=',')
        OLD_asymmetry_error_minus = np.loadtxt('output/strangeness_asymmetry_errors/' + PDF_set + '_OLD_minus.txt', delimiter=',')
        NEW_asymmetry_error_plus = np.loadtxt('output/strangeness_asymmetry_errors/' + PDF_set + '_' + \
                                                which_cross_sections_included + '_NEW_plus.txt', delimiter=',')
        NEW_asymmetry_error_minus = np.loadtxt('output/strangeness_asymmetry_errors/' + PDF_set + '_' + \
                                                which_cross_sections_included + '_NEW_minus.txt', delimiter=',')

    asymmetry_OLD = OLD_s - OLD_sbar
    asymmetry_NEW = NEW_s - NEW_sbar

    for x_index in range(len(asymmetry_OLD)):
        if (asymmetry_NEW[x_index] != asymmetry_NEW[x_index]):
            asymmetry_NEW[x_index] = 0.
        if (asymmetry_OLD[x_index] != asymmetry_OLD[x_index]):
            asymmetry_OLD[x_index] = 0.

    ax1.plot(x, asymmetry_OLD * 10**2, color='black', zorder=2, label=PDF_set_labels[PDF_index])
    ax1.plot(x, asymmetry_NEW * 10**2, color='red', zorder=2, label='reweighted')

    if (plot_errors_flag):
        PDF_uncertainty = ax1.fill_between(x, (asymmetry_OLD - OLD_asymmetry_error_minus) * 10**2, (asymmetry_OLD + OLD_asymmetry_error_plus) * 10**2,
                        color='lightgray', zorder=0, rasterized=True)
        reweighted_PDF_uncertainty = ax1.fill_between(x, (asymmetry_NEW - NEW_asymmetry_error_minus) * 10**2, (asymmetry_NEW + NEW_asymmetry_error_plus) * 10**2,
                        color='salmon', zorder=1, alpha=0.6, rasterized=True)

    # ax2.plot(x, (asymmetry_NEW - asymmetry_OLD) * 10**2, color='red', zorder=2)
    ax2.plot(x, (asymmetry_NEW - asymmetry_OLD) * 10**2, color='red', zorder=2)

    if (plot_errors_flag):
        ax2.fill_between(x, OLD_asymmetry_error_plus * 10**2,
                        -OLD_asymmetry_error_minus * 10**2, color='lightgray', zorder=0)
        ax2.fill_between(x, (asymmetry_NEW + NEW_asymmetry_error_plus - asymmetry_OLD) * 10**2,
                        (asymmetry_NEW - NEW_asymmetry_error_minus - asymmetry_OLD) * 10**2, color='salmon', zorder=1, alpha=0.6)

    plt.xscale('log')

    plt.xlim(10**(-5), 1)
    if (PDF_set == 'MSHT20nlo_as118'):
        ax1.set_ylim(-0.7, 2.3)
        ax2.set_ylim(-0.7, 0.9)
        #ax2.set_ylim(-0.5, 3.5)
    elif (PDF_set == 'NNPDF40_nlo_pch_as_01180'):
        ax1.set_ylim(-6, 15)
        ax2.set_ylim(-3, 7)

    ax1.plot([10**(-6), 2], [0, 0], color='black', zorder=1)
    ax2.plot([10**(-6), 2], [0, 0], color='black', zorder=1)

    plt.xlabel(r'$x$', fontsize=axis_label_font_size)
    ax1.set_ylabel(r'$x(s - \overline{s})$', fontsize=axis_label_font_size)
    ax2.set_ylabel('Reweighted\n - Original', fontsize=axis_label_font_size)
    #ax2.set_ylabel(r'$\frac{\text{Reweighted}}{\text{Original}}$', fontsize=axis_label_font_size * 1.3)

    ax1.tick_params(axis='both', which='major', labelsize=axis_font_size)
    ax2.tick_params(axis='both', which='major', labelsize=axis_font_size)

    text_x = 10**(-5) * 1.8
    text_y1 = 1.6

    ax1.text(text_x, text_y1, r'$\mu_\text{fact} = M_W$', fontsize=font_size)

    if (which_cross_sections_included == 'D'):
        ax1.text(text_x, 1.9, r'$R_c^\pm (D^\pm)$', fontsize=font_size)
    else:
        ax1.text(text_x, 1.9, r'$R_c^\pm (D^\pm, D^{*\pm})$', fontsize=font_size)

    ax1.text(0.01, 1.01, r'$\times 10^{-2}$',
        transform=ax1.transAxes,
        fontsize=13, va='bottom', ha='left')

    legend1 = ax1.legend(loc='lower left', fontsize=legend_fontsize, framealpha=1, bbox_to_anchor=(0., 0.5))
    legend2 = ax1.legend([PDF_uncertainty, reweighted_PDF_uncertainty],
                        ['PDF error (68\% C.L.)', 'Reweighted'],
                        loc='lower left', fontsize=legend_fontsize, framealpha=1, bbox_to_anchor=(0., 0.3))
    ax1.add_artist(legend1)

    #ax2.set_yticks([0, 1, 2, 3])

    # Configure ticks to appear on all sides
    ax1.tick_params(direction='in', top=True, right=True)
    ax2.tick_params(direction='in', top=True, right=True)

    # Add minor ticks
    ax1.minorticks_on()
    ax1.tick_params(which='both', direction='in', top=True, right=True)
    ax2.minorticks_on()
    ax2.tick_params(which='both', direction='in', top=True, right=True)

    plt.tight_layout()

    plt.savefig('plots/strange_asymmetry/strangeness_asymmetry_' + PDF_set + '_' + which_cross_sections_included + '_v2.pdf', dpi=300)
    plt.show()


def weights(PDF_set, which_cross_sections_included, tolerance):
    plt.figure(figsize=(6, 4))
        
    weights = np.loadtxt('output/wmin_' + PDF_set + '_' + which_cross_sections_included + '_t_' + tolerance + '.txt')

    for i in range(len(weights)):
        if (abs(weights[i]) > 0.2):
            print(colored(i + 1, 'green'), colored(weights[i], 'green'))
        else:
            print(i + 1, weights[i])

    X = np.arange(1, 32.01, 1)

    plt.scatter(X, weights, color='black')

    plt.plot([-1, 60], [0, 0], color='black', zorder=0)

    plt.xlim(0, 33)

    plt.xlabel(r'PDF member')
    plt.ylabel('Weight', fontsize=axis_label_font_size)
    
    plt.tight_layout()
    plt.savefig('plots/weights/' + PDF_set + '_' + which_cross_sections_included + '.pdf')

    plt.close()


def Rcpm_approximation(PDF_set, which_cross_sections_included, tolerances, colors, labels):
    fig, ax = plt.subplots(figsize=(6, 6))
    
    s_old = np.loadtxt('/home/alankovh/Documents/PDFs/Values/' + PDF_set + '/flavor_id_3/best_vals_valence_0.txt', delimiter=',')
    sv_old = np.loadtxt('/home/alankovh/Documents/PDFs/Values/' + PDF_set + '/flavor_id_3/best_vals_valence_1.txt', delimiter=',')
    dv_old = np.loadtxt('/home/alankovh/Documents/PDFs/Values/' + PDF_set + '/flavor_id_1/best_vals_valence_1.txt', delimiter=',')
    
    error_plus = np.loadtxt('output/Rcpm_approximation_errors/' + PDF_set + '_OLD_plus.txt', delimiter=',')
    error_minus = np.loadtxt('output/Rcpm_approximation_errors/' + PDF_set + '_OLD_minus.txt', delimiter=',')

    epsilon = 0.222**2 / 0.975**2

    if (PDF_set != 'NNPDF40_nlo_pch_as_01180'):
        for tolerance_index in range(len(tolerances)):
            tolerance = tolerances[tolerance_index]

            x = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + \
                '/flavor_3_t_' + tolerances[0] + '_best.txt', delimiter=',', max_rows=1)
            s_new = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + \
                            '/flavor_3_t_' + tolerance + '_best.txt', delimiter=',', skiprows=1, max_rows=1)
            d_new = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + \
                                '/flavor_1_t_' + tolerance + '_best.txt', delimiter=',', skiprows=1, max_rows=1)
            sbar_new = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + \
                                '/flavor_-3_t_' + tolerance + '_best.txt', delimiter=',', skiprows=1, max_rows=1)
            dbar_new = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + \
                                '/flavor_-1_t_' + tolerance + '_best.txt', delimiter=',', skiprows=1, max_rows=1)

            plt.plot(x, 1. - (epsilon * (d_new - dbar_new) + (s_new - sbar_new)) / s_new, zorder=2,
                    label=labels[tolerance_index], color=colors[tolerance_index])
            
            print(tolerance)
            print(sum(s_new))
            print(sum(s_new - sbar_new))
            print(sum(d_new - dbar_new))
            print()
            print(sum(s_old))
            print(sum(sv_old))
            print(sum(dv_old))
            print()
            print()
    
    else:
        x = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + \
                '/flavor_3_best.txt', delimiter=',', max_rows=1)
        s_new = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + \
                        '/flavor_3_best.txt', delimiter=',', skiprows=1, max_rows=1)
        d_new = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + \
                            '/flavor_1_best.txt', delimiter=',', skiprows=1, max_rows=1)
        sbar_new = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + \
                            '/flavor_-3_best.txt', delimiter=',', skiprows=1, max_rows=1)
        dbar_new = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + \
                            '/flavor_-1_best.txt', delimiter=',', skiprows=1, max_rows=1)
        
        print(sum(s_new))
        print(sum(s_new - sbar_new))
        print(sum(d_new - dbar_new))
        print()
        print(sum(s_old))
        print(sum(sv_old))
        print(sum(dv_old))

        plt.plot(x, 1. - (epsilon * (d_new - dbar_new) + (s_new - sbar_new)) / s_new, zorder=2,
                label='Reweighted', color='red')

    plt.plot(x, 1. - (epsilon * dv_old + sv_old) / s_old, color='black', label='original', zorder=1)

    plt.fill_between(x, 1. - (epsilon * dv_old + sv_old) / s_old + error_plus,
                        1. - (epsilon * dv_old + sv_old) / s_old - error_minus,
                        color='lightgray', alpha=0.5, zorder=0)
    
    ax.xaxis.set_major_locator(ticker.FixedLocator([1e-2, 1e-1]))
    ax.set_xticklabels([r'$10^{-2}$', r'$10^{-1}$'])
    ax.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*0.1))
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())

    ax.set_xlim(1e-2, 1e-1)
    ax.set_xscale('log')
    ax.xaxis.set_major_locator(ticker.FixedLocator([1e-2, 1e-1]))
    ax.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*0.1))
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())

    ax.minorticks_on()
    ax.tick_params(which='both', direction='in', top=True, right=True, zorder=100, labelsize=16)

    plt.ylim([0.45, 1.1])

    plt.xlabel(r'$x$', fontsize=axis_label_font_size)
    plt.ylabel(r'$R_c^\pm$', fontsize=axis_label_font_size)
    
    plt.legend(loc='lower left', fontsize=15)
    plt.tight_layout()
    plt.savefig('plots/Rcpm_approximation/' + PDF_set + '_' + which_cross_sections_included + '.pdf')
    plt.show()


def Rcpm_approximation_v2(PDF_set, which_cross_sections_included, tolerances, colors, labels):
    fig, ax = plt.subplots(figsize=(6, 6))
    
    s_old = np.loadtxt('/home/alankovh/Documents/PDFs/Values/' + PDF_set + '/flavor_id_3/best_vals_valence_0.txt', delimiter=',')
    d_old = np.loadtxt('/home/alankovh/Documents/PDFs/Values/' + PDF_set + '/flavor_id_1/best_vals_valence_0.txt', delimiter=',')
    sv_old = np.loadtxt('/home/alankovh/Documents/PDFs/Values/' + PDF_set + '/flavor_id_3/best_vals_valence_1.txt', delimiter=',')
    dv_old = np.loadtxt('/home/alankovh/Documents/PDFs/Values/' + PDF_set + '/flavor_id_1/best_vals_valence_1.txt', delimiter=',')

    epsilon = 0.222**2 / 0.975**2

    if (PDF_set != 'NNPDF40_nlo_pch_as_01180'):
        for tolerance_index in range(len(tolerances)):
            tolerance = tolerances[tolerance_index]

            x = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + \
                '/flavor_3_t_' + tolerances[0] + '_best.txt', delimiter=',', max_rows=1)
            s_new = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + \
                            '/flavor_3_t_' + tolerance + '_best.txt', delimiter=',', skiprows=1, max_rows=1)
            d_new = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + \
                                '/flavor_1_t_' + tolerance + '_best.txt', delimiter=',', skiprows=1, max_rows=1)
            sbar_new = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + \
                                '/flavor_-3_t_' + tolerance + '_best.txt', delimiter=',', skiprows=1, max_rows=1)
            dbar_new = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + \
                                '/flavor_-1_t_' + tolerance + '_best.txt', delimiter=',', skiprows=1, max_rows=1)

            plt.plot(x, (epsilon * dbar_new + sbar_new) / (epsilon * d_new + s_new), zorder=2,
                    label=labels[tolerance_index], color=colors[tolerance_index])
            
            print(tolerance)
            print(sum(s_new))
            print(sum(s_new - sbar_new))
            print(sum(d_new - dbar_new))
            print()
            print(sum(s_old))
            print(sum(sv_old))
            print(sum(dv_old))
            print()
            print()
    
    else:
        x = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + \
                '/flavor_3_best.txt', delimiter=',', max_rows=1)
        s_new = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + \
                        '/flavor_3_best.txt', delimiter=',', skiprows=1, max_rows=1)
        d_new = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + \
                            '/flavor_1_best.txt', delimiter=',', skiprows=1, max_rows=1)
        sbar_new = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + \
                            '/flavor_-3_best.txt', delimiter=',', skiprows=1, max_rows=1)
        dbar_new = np.loadtxt('output/new_PDF_vals/' + PDF_set + '/' + which_cross_sections_included + \
                            '/flavor_-1_best.txt', delimiter=',', skiprows=1, max_rows=1)
        
        print(sum(s_new))
        print(sum(s_new - sbar_new))
        print(sum(d_new - dbar_new))
        print()
        print(sum(s_old))
        print(sum(sv_old))
        print(sum(dv_old))

        #plt.plot(x, , zorder=2,
        #        label='Reweighted', color='red')

    plt.plot(x, (epsilon * (d_old - dv_old) + (s_old - sv_old)) / (epsilon * d_old + s_old), color='black', label='original', zorder=1)

    #plt.fill_between(x, 1. - (epsilon * dv_old + sv_old) / s_old + error_plus,
    #                    1. - (epsilon * dv_old + sv_old) / s_old - error_minus,
    #                    color='lightgray', alpha=0.5, zorder=0)
    
    ax.xaxis.set_major_locator(ticker.FixedLocator([1e-2, 1e-1]))
    ax.set_xticklabels([r'$10^{-2}$', r'$10^{-1}$'])
    ax.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*0.1))
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())

    ax.set_xlim(1e-2, 1e-1)
    ax.set_xscale('log')
    ax.xaxis.set_major_locator(ticker.FixedLocator([1e-2, 1e-1]))
    ax.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*0.1))
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())

    ax.minorticks_on()
    ax.tick_params(which='both', direction='in', top=True, right=True, zorder=100, labelsize=16)

    plt.ylim([0.45, 1.1])

    plt.xlabel(r'$x$', fontsize=axis_label_font_size)
    plt.ylabel(r'$R_c^\pm$', fontsize=axis_label_font_size)
    
    plt.legend(loc='lower left', fontsize=15)
    plt.tight_layout()
    plt.savefig('plots/Rcpm_approximation/' + PDF_set + '_' + which_cross_sections_included + '.pdf')
    plt.show()



PDF_set = 'MSHT20nlo_as118'
#PDF_set = 'NNPDF40_nlo_pch_as_01180'
#PDF_set = 'CT18ANLO'
flavors = [1, 3, 21]
which_cross_sections_included = 'both'

#ratio_HESSIAN(PDF_set, which_cross_sections_included, [[1, 1], [3, 0], [21, 0], [1, 0]],
#              [r'$d_\text{valence}$', r'$s$', r'$g$', r'$d$'], ['3.16', '1', '0.1', '0.01', '0.001'],
#              [r'$t = \sqrt{10}$', r'$t = 1$', r'$t = 0.1$', r'$t = 0.01$', r'$0.001$'], ['red', 'green', 'blue', 'purple', 'orange'],
#              [0.9, 0.6, 0.9, 0.9], [1.1, 1.15, 1.05, 1.1])
#ratio_HESSIAN(PDF_set, which_cross_sections_included, [[1, 1], [3, 0], [3, 1], [21, 0], [1, 0]],
#              [r'$d_\text{valence}$', r'$s$', r'$s_\text{valence}$', r'$g$', r'$d$'], ['3.16', '1', '0.1', '0.01'],
#              [r'$t = \sqrt{10}$', r'$t = 1$', r'$t = 0.1$', r'$t = 0.01$'], ['red', 'green', 'blue', 'purple', 'orange'],
#              [0.9, 0.9, -1, 0.94, 0.9], [1.1, 1.1, 3, 1.03, 1.1])
#ratio_MC(PDF_set, which_cross_sections_included, [[1, 1], [3, 1], [3, 0], [21, 0]], [r'$d_\text{valence}$', r'$s_\text{valence}$', r'$s$', r'$g$'],
#                   [0.4, -3, 0.6, 0.94], [1.2, 7, 1.25, 1.03])
#ratio_to_other_PDF('MSHT20nlo_as118', 'CT18ANLO', flavors, which_cross_sections_included)
#ratio_of_ratio(1)
#ratio_of_ratio(3)
#absolute('0.01', 'both', [1, 1, 3, 3, 21], [1, 0, 1, 0, 0])
Rcpm_approximation(PDF_set, which_cross_sections_included, ['3.16', '1', '0.3', '0.01'], ['red', 'green', 'blue', 'orange', 'purple', 'pink'],
                   [r'$t = \sqrt{10}$', r'$t = 1$', r'$t = 0.3$', r'$t = 0.01$', r'$t = 0.001$', r'$t = 0.0001$'])
#Rcpm_approximation_v2(PDF_set, which_cross_sections_included, ['3.16', '1', '0.3'], ['red', 'green', 'blue', 'orange', 'purple', 'pink'],
#                   [r'$t = \sqrt{10}$', r'$t = 1$', r'$t = 0.3$', r'$t = 0.01$', r'$t = 0.001$', r'$t = 0.0001$'])
#Rcpm_OLD_and_NEW_and_DATA('eta_lept', which_cross_sections_included, '0.01')
#Rcpm_OLD_and_NEW_and_DATA('pTD', which_cross_sections_included, '0.1')
#strangeness_asymmetry_v2(PDF_set, which_cross_sections_included, True, False)
#weights(PDF_set, which_cross_sections_included, '1')

