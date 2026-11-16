import numpy as np
import matplotlib.pyplot as plt
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

fontsize = 17
axis_label_font_size = 19
axis_font_size = 14
legend_fontsize = 14

pathNLO = "../output/W-D+/NLO/minus/KKKS08_opal/frag_main_scale/scale_variation/"
pathLO = "../output/W-D+/subtraction/minus/KKKS08_opal/frag_main_scale/scale_variation/"

NLOdata = np.loadtxt(pathNLO + "CT18ANLOprecise/central/0_m_charm_central_vals.txt", delimiter=',')

normal = np.zeros(5)

pTD_data_min = 8
pTD_data_bin_width = 0.5
pTD_bins = np.array([8, 12, 20, 40, 80, 150])

for eta_lept_index in range(5):
    LOdata = np.loadtxt(pathLO + "CT18ANLO/central/" + str(eta_lept_index) + "_m_charm_central_vals.txt", delimiter=',')

    bin_index = 0
    for pTD_index in range(284 + 1):
        if (pTD_data_min + (pTD_index + 1 / 2) * pTD_data_bin_width > pTD_bins[bin_index + 1]):
            if (bin_index < 4):
                bin_index += 1
            else:
                break

        if (eta_lept_index == 0):
            print(sum(NLOdata[pTD_index, :]))
            normal[bin_index] += sum(NLOdata[pTD_index, :])

        #normal -=sum(LOdata[pTD_index, :]) * 2.

print(normal)

        


smallmc = np.array([11.65841617, 12.44684877, 11.04695346, 4.16169237, 0.75769082])
normal = np.array([11.64387262, 12.45569112, 11.03653085, 4.15223348, 0.75066562])

print(normal)

smallmc_no_large_log = np.array([10.24254336, 10.37854752, 8.38538374, 2.80562687, 0.47807623])

normal_no_large_log = np.array([10.44019221, 10.69029741, 8.76106419, 2.99114801, 0.51074882])

ratio1 = smallmc / normal
ratio2 = smallmc_no_large_log / normal_no_large_log




fig, ax = plt.subplots(figsize=(6, 4.5))

ax.hlines(ratio1, pTD_bins[:-1], pTD_bins[1:], zorder=3, color='red', label=r'$\mu_\text{frac} = M_W$')
ax.hlines(ratio2, pTD_bins[:-1], pTD_bins[1:], zorder=3, color='blue', label=r'$\mu_\text{frac} = M_W$' + '\n' + r'$(\text{large log double-counted})$')

plt.plot([8, 150], [1, 1], color='black', zorder=0)

ax.set_xlabel(r'$p_{T, D}$ [GeV]', fontsize=axis_label_font_size)
ax.set_ylabel(r'$\frac{\sigma(m_c = 0.65\text{ GeV})}{\sigma(m_c = 1.3\text{ GeV})}$', fontsize=axis_label_font_size * 1.2)

plt.xscale('log')

ax.minorticks_on()
ax.tick_params(
    which='both',    # major and minor ticks
    direction='in',
    left=True,
    right=True,
    bottom=False,
    top=False,
    labelsize=axis_font_size
)

for i in range(1, len(pTD_bins) - 1):
    ax.axvline(pTD_bins[i], color='gray', linewidth=0.5, ymax=0.7, zorder=0)

plt.xticks(pTD_bins, [f'{tick:.0f}' for tick in pTD_bins])

plt.xlim(8, 150)
ax.set_ylim(0.93, 1.05)

text_x = 9
text_y1 = 1.04
text_y2 = 1.03

ax.text(text_x, text_y1, 'CT18ANLO', fontsize=fontsize)

ax.legend(loc='upper right', fontsize=legend_fontsize)
plt.tight_layout()

plt.savefig('../plots/13 TeV/W-D+/KKKS08_opal/W-D+_KKKS08_opal_mc_variation.pdf')

plt.show()