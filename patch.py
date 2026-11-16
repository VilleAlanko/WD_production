# I lost the original script, this is a new one that does not have all the things the old one had.

import numpy as np
from pathlib import Path

QCD_order = "NLO"
process = "W+c-"
task = "scale_variation"
PDF_set = "NNPDF40_nnlo_as_01180_charmasy"
num_err_members = 100
# 0 (off), 0.1 -1 1
mcvar = 0

main_directory = "/home/alankovh/Documents/PhD_projects/WD_production/"

if (task == "scale_variation"):
    filename_raw = "W_cjet_nlo_" + PDF_set + "_80___80___cW_xi"

    if (mcvar == 0):
        raw_vals_path = main_directory + "raw_input/" + process + "/" + QCD_order + "/" + task + "/" + PDF_set + "/"
    elif (mcvar == 0.1):
        raw_vals_path = main_directory + "raw_input/" + process + "/" + QCD_order + "/" + task + "/" + PDF_set + "mcvar/"
    elif (mcvar == -1):
        raw_vals_path = main_directory + "raw_input/" + process + "/" + QCD_order + "/" + task + "/" + PDF_set + "mcvar/mc_variation_down/"
    else:
       exit(1)

    if (mcvar == 0):
        vals_path = main_directory + "input/" + process + "/" + QCD_order + "/" + task + "/" + PDF_set + "/central/"
    elif (mcvar == 0.1):
        vals_path = main_directory + "input/" + process + "/" + QCD_order + "/" + task + "/" + PDF_set + "mcvar/central/"
    elif (mcvar == -1):
        vals_path = main_directory + "input/" + process + "/" + QCD_order + "/" + task + "/" + PDF_set + "mcvar/central/"
    else:
        exit(1)
        

    for eta_lept_index in range(5):
        vals_list = []

        for i in range(3):
            raw_vals = np.loadtxt(raw_vals_path + str(i) + "/" + filename_raw + str(eta_lept_index) + ".txt", skiprows=5)
            
            if (i != 1):
                vals_list.append(raw_vals)
            else:
                raw_vals[:, 2] *= 1.
                vals_list.append(raw_vals)

        vals = np.concatenate(vals_list, axis=0)

        if (mcvar == 0.1 or mcvar == 0):
            filename = vals_path + str(eta_lept_index) + "_m_charm_central.txt"
        elif (mcvar == -1):
            filename = vals_path + str(eta_lept_index) + "_m_charm_down.txt"
        else:
            exit(1)

        file_path = Path(filename)
        file_path.parent.mkdir(parents=True, exist_ok=True)
        
        np.savetxt(file_path, vals, delimiter=" ")

if (task == "pdf_errs"):
    filename_raw = "W_cjet_nlo_80___80___a74__cW_xi"
    raw_vals_path = main_directory + "raw_input/" + process + "/" + QCD_order + "/" + task + "/pp/" + PDF_set + "/"
    vals_path = main_directory + "input/" + process + "/" + QCD_order + "/" + task + "/" + PDF_set + "/"

    
    for eta_lept_index in range(5):
        for member in range(num_err_members + 1):
            vals_list = []

            for i in range(4):
                if (member == 0):
                    raw_vals_filename = raw_vals_path + str(i) + "/" + filename_raw + str(eta_lept_index) + ".txt"
                else:
                    raw_vals_filename = raw_vals_path + str(i) + "/" + filename_raw + str(eta_lept_index) + "_Nca" + str(member) + ".txt"

                raw_vals = np.loadtxt(raw_vals_filename, skiprows=5)
                vals_list.append(raw_vals)

            vals = np.concatenate(vals_list, axis=0)

            if (member == 0):
                member0 = vals
            else:
                vals[:, 2] = member0[:, 2] - vals[:, 2]

            filename = vals_path + "Nca" + str(member) + "_" + str(eta_lept_index) + "_m_charm_central.txt"

            file_path = Path(filename)
            file_path.parent.mkdir(parents=True, exist_ok=True)
            
            np.savetxt(file_path, vals, delimiter=" ")

    
