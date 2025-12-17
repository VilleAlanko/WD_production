#include "INIReader.h"
#include "LHAPDF/LHAPDF.h"
#include <iostream>
#include <fstream>
#include <cstdlib>  // for setenv

using namespace LHAPDF;
using namespace std;

const double Q = 80.385;

const double log10x_min = -5;
const double log10x_max = 0;
const double log10x_step = 0.001;
const int N = int((log10x_max - log10x_min) / log10x_step) + 1;

vector <int> split_to_ints(const string& s, char delimiter) {
    vector <int> result;
    stringstream ss(s);
    string item;
    while (getline(ss, item, delimiter)) {
        if (!item.empty()) {
            result.push_back(stoi(item));
        }
    }
    return result;
}

vector <double> load_1D_double_array_from_txt(const string& filename)
{
    ifstream file(filename);

    if (!file.is_open())
    {
        cout << "Error in the function 'load_1D_double_array_from_txt': could not open the following file: " << filename << endl;
        exit(1);
    }

    vector<double> values;
    string line;

    while (getline(file, line))
    {
        stringstream ss(line);
        string value;

        while (getline(ss, value, ','))
        {
            values.push_back(stod(value));
        }
    }

    file.close();
    return values;
}

vector <vector <double>> load_2D_double_array_from_txt(string filename, char delimiter)
{
    ifstream file(filename);

    if (!file.is_open())
    {
        cout << "Error in the function 'load_2D_double_array_from_txt': could not open the following file: " << filename << endl;
        exit(1);
    }

    vector <vector <double>> values;
    string line;

    while (getline(file, line))
    {
        vector <double> row;
        stringstream ss(line);
        string value;

        while (getline(ss, value, delimiter))
        {
            row.push_back(stod(value));
        }

        values.push_back(row);
    }

    file.close();

    return values;
}

void HESSIAN(int num_err_members, vector <int> flavors, string PDF_set,
            string which_cross_sections_included, string kinematic_quantity_for_new_Rcpm, int N_bins_in_kinematic_quantity)
{
    vector <double> wmin = load_1D_double_array_from_txt("output/wmin_" + PDF_set + "_" + which_cross_sections_included + ".txt");
    vector <vector <double>> dw = load_2D_double_array_from_txt("output/dw_" + PDF_set + "_" + which_cross_sections_included + ".txt", ' ');

    vector <double> x_vals;

    for (int flavor_index = 0; flavor_index < int(flavors.size()); flavor_index++)
    {
        int flavor = flavors[flavor_index];

        array <double, N> PDF_best_OLD = {};
        vector <array<double, N>> PDF_plus_members_OLD(num_err_members / 2);
        vector <array<double, N>> PDF_minus_members_OLD(num_err_members / 2);

        for (int member_id = 0; member_id < num_err_members + 1; member_id++)
        {
            const PDF* pdf = mkPDF(PDF_set, member_id);
        
            for (int x_index = 0; x_index < N; x_index++)
            {
                double x = pow(10, log10x_min + log10x_step * x_index);

                if (flavor_index == 0 && member_id == 0)
                {
                    x_vals.push_back(x);
                }

                double xPDF_val = pdf->xfxQ(flavor, x, Q);
                
                if (member_id == 0)
                {
                    PDF_best_OLD[x_index] = xPDF_val;
                }
                else if (member_id % 2 == 0)
                {
                    PDF_minus_members_OLD[int(member_id / 2) - 1][x_index] = xPDF_val;
                }
                else
                {
                    PDF_plus_members_OLD[int((member_id - 1) / 2)][x_index] = xPDF_val;
                }
            }
        }

        array <double, N> best_member_NEW = {};
        vector <array <double, N>> PDF_plus_members_NEW(num_err_members / 2);
        vector <array <double, N>> PDF_minus_members_NEW(num_err_members / 2);

        for (int x_index = 0; x_index < N; x_index++)
        {
            best_member_NEW[x_index] = PDF_best_OLD[x_index];

            for (int k = 0; k < int(num_err_members / 2); k++)
            {
                best_member_NEW[x_index] += (PDF_plus_members_OLD[k][x_index] - PDF_minus_members_OLD[k][x_index]) / 2. * wmin[k];
            }
        }
        for (int x_index = 0; x_index < N; x_index++)
        {
            for (int k = 0; k < int(num_err_members / 2); k++)
            {
                PDF_plus_members_NEW[k][x_index] = best_member_NEW[x_index];
                PDF_minus_members_NEW[k][x_index] = best_member_NEW[x_index];
                
                for (int i = 0; i < int(num_err_members / 2); i++)
                {
                    PDF_plus_members_NEW[k][x_index] +=  (PDF_plus_members_OLD[i][x_index] -
                        PDF_minus_members_OLD[i][x_index]) / 2. * dw[i][k];
                    PDF_minus_members_NEW[k][x_index] -= (PDF_plus_members_OLD[i][x_index] -
                        PDF_minus_members_OLD[i][x_index]) / 2. * dw[i][k];
                }
            }
        }

        string filename = "output/old_PDF_vals/" + PDF_set + "/" + which_cross_sections_included + "/flavor_" + to_string(flavor) + ".txt";
        ofstream outfile1(filename, ios::out);

        for (int x_index = 0; x_index < N; x_index++)
        {
            outfile1 << x_vals[x_index];
            if (x_index != N - 1)
            {
                outfile1 << ",";
            }
        }
        outfile1 << endl;
        for (int x_index = 0; x_index < N; x_index++)
        {
            outfile1 << PDF_best_OLD[x_index];
            if (x_index != N - 1)
            {
                outfile1 << ",";
            }
        }
        outfile1.close();

        filename = "output/new_PDF_vals/" + PDF_set + "/" + which_cross_sections_included + "/flavor_" + to_string(flavor) + "_best.txt";
        ofstream outfile2(filename, ios::out);

        for (int x_index = 0; x_index < N; x_index++)
        {
            outfile2 << x_vals[x_index];
            if (x_index != N - 1)
            {
                outfile2 << ",";
            }
        }
        outfile2 << endl;
        for (int x_index = 0; x_index < N; x_index++)
        {
            outfile2 << best_member_NEW[x_index];
            if (x_index != N - 1)
            {
                outfile2 << ",";
            }
        }
        
        outfile2.close();

        filename = "output/new_PDF_vals/" + PDF_set + "/" + which_cross_sections_included + "/flavor_" + to_string(flavor) + "_plus.txt";
        ofstream outfile2plus(filename, ios::out);

        for (int member = 0; member < num_err_members / 2; member++)
        {
            for (int x_index = 0; x_index < N; x_index++)
            {
                outfile2plus << PDF_plus_members_NEW[member][x_index];
                if (x_index != N - 1)
                {
                    outfile2plus << ",";
                }
            }

            if (member != num_err_members - 1)
            {
                outfile2plus << endl;
            }
        }
        outfile2plus.close();

        filename = "output/new_PDF_vals/" + PDF_set + "/" + which_cross_sections_included + "/flavor_" + to_string(flavor) + "_minus.txt";
        ofstream outfile2minus(filename, ios::out);

        for (int member = 0; member < num_err_members / 2; member++)
        {
            for (int x_index = 0; x_index < N; x_index++)
            {
                outfile2minus << PDF_minus_members_NEW[member][x_index];
                if (x_index != N - 1)
                {
                    outfile2minus << ",";
                }
            }

            if (member != num_err_members - 1)
            {
                outfile2minus << endl;
            }
        }
        outfile2minus.close();
    }

    vector <double> Rcpm_best_OLD = load_1D_double_array_from_txt("input/theory_values/HESSIAN/best/" + 
                                                            kinematic_quantity_for_new_Rcpm + "_" + 
                                                            which_cross_sections_included + "_" + PDF_set + "_best.txt");
    vector <vector <double>> Rcpm_plus_members_OLD = load_2D_double_array_from_txt("input/theory_values/HESSIAN/variation/" + 
                                                        kinematic_quantity_for_new_Rcpm + "_" + 
                                                        which_cross_sections_included + "_" + PDF_set + "_plus.txt", ',');
    vector <vector <double>> Rcpm_minus_members_OLD = load_2D_double_array_from_txt("input/theory_values/HESSIAN/variation/" + 
                                                        kinematic_quantity_for_new_Rcpm + "_" + 
                                                        which_cross_sections_included + "_" + PDF_set + "_minus.txt", ',');

    vector <double> Rcpm_best_NEW(N_bins_in_kinematic_quantity);
    vector <vector <double>> Rcpm_plus_members_NEW(num_err_members / 2,
                                                       vector <double>(N_bins_in_kinematic_quantity));

    vector <vector <double>> Rcpm_minus_members_NEW(num_err_members / 2,
                                                       vector <double>(N_bins_in_kinematic_quantity));
    
    for (int bin_index = 0; bin_index < N_bins_in_kinematic_quantity; bin_index++)
    {
        Rcpm_best_NEW[bin_index] = Rcpm_best_OLD[bin_index];

        for (int k = 0; k < int(num_err_members / 2); k++)
        {
            Rcpm_best_NEW[bin_index] += (Rcpm_plus_members_OLD[bin_index][k] - Rcpm_minus_members_OLD[bin_index][k]) / 2. * wmin[k];
        }
    }

    for (int bin_index = 0; bin_index < N_bins_in_kinematic_quantity; bin_index++)
    {
        for (int k = 0; k < int(num_err_members / 2); k++)
        {
            Rcpm_plus_members_NEW[k][bin_index] = Rcpm_best_NEW[bin_index];
            Rcpm_minus_members_NEW[k][bin_index] = Rcpm_best_NEW[bin_index];

            for (int i = 0; i < int(num_err_members / 2); i++)
            {
                Rcpm_plus_members_NEW[k][bin_index] +=  (Rcpm_plus_members_OLD[bin_index][i] -
                    Rcpm_minus_members_OLD[bin_index][i]) / 2. * dw[i][k];
                Rcpm_minus_members_NEW[k][bin_index] -= (Rcpm_plus_members_OLD[bin_index][i] -
                    Rcpm_minus_members_OLD[bin_index][i]) / 2. * dw[i][k];
            }
        }
    }

    vector <double> Rcpm_err_up_NEW(N_bins_in_kinematic_quantity);
    vector <double> Rcpm_err_down_NEW(N_bins_in_kinematic_quantity);

    for (int bin_index = 0; bin_index < N_bins_in_kinematic_quantity; bin_index++)
    {
        for (int k = 0; k < int(num_err_members / 2); k++)
        {
            Rcpm_err_up_NEW[bin_index] += pow(max(max(Rcpm_plus_members_NEW[k][bin_index] - Rcpm_best_NEW[bin_index],
                                            Rcpm_minus_members_NEW[k][bin_index] - Rcpm_best_NEW[bin_index]), 0.), 2);
            Rcpm_err_down_NEW[bin_index] += pow(max(max(Rcpm_best_NEW[bin_index] - Rcpm_minus_members_NEW[k][bin_index],
                                            Rcpm_best_NEW[bin_index] - Rcpm_plus_members_NEW[k][bin_index]), 0.), 2);
        }
    }
    
    for (int bin_index = 0; bin_index < N_bins_in_kinematic_quantity; bin_index++)
    {
        Rcpm_err_up_NEW[bin_index] = sqrt(Rcpm_err_up_NEW[bin_index]);
        Rcpm_err_down_NEW[bin_index] = sqrt(Rcpm_err_down_NEW[bin_index]);
    }

    string filename = "output/new_Rcpm_" + kinematic_quantity_for_new_Rcpm + "/" + PDF_set + "/" + which_cross_sections_included + "/best.txt";
    ofstream outfile3(filename, ios::out);

    for (int bin_index = 0; bin_index < N_bins_in_kinematic_quantity; bin_index++)
    {
        outfile3 << Rcpm_best_NEW[bin_index];
        if (bin_index != N_bins_in_kinematic_quantity - 1)
        {
            outfile3 << endl;
        }
    }
    outfile3.close();

    filename = "output/new_Rcpm_" + kinematic_quantity_for_new_Rcpm + "/" + PDF_set + "/" + which_cross_sections_included + "/error_up.txt";
    ofstream outfile4(filename, ios::out);

    for (int bin_index = 0; bin_index < N_bins_in_kinematic_quantity; bin_index++)
    {
        outfile4 << Rcpm_err_up_NEW[bin_index];
        if (bin_index != N_bins_in_kinematic_quantity - 1)
        {
            outfile4 << endl;
        }
    }
    outfile4.close();

    filename = "output/new_Rcpm_" + kinematic_quantity_for_new_Rcpm + "/" + PDF_set + "/" + which_cross_sections_included + "/error_down.txt";
    ofstream outfile5(filename, ios::out);

    for (int bin_index = 0; bin_index < N_bins_in_kinematic_quantity; bin_index++)
    {
        outfile5 << Rcpm_err_down_NEW[bin_index];
        if (bin_index != N_bins_in_kinematic_quantity - 1)
        {
            outfile5 << endl;
        }
    }
    outfile5.close();
}

void MC(int num_err_members, vector <int> flavors, string PDF_set, string which_cross_sections_included,
        string kinematic_quantity_for_new_Rcpm, int N_bins_in_kinematic_quantity)
{
    vector <double> omega(num_err_members);

    ifstream file("output/omega.txt");
    if (!file) {
        cerr << "Error opening file\n";
        exit(1);
    }

    for (int i = 0; i < num_err_members; ++i) {
        file >> omega[i];
    }

    for (int flavor_index = 0; flavor_index < int(flavors.size()); flavor_index++)
    {
        int flavor = flavors[flavor_index];

        array <double, N> PDF_best_OLD = {};
        array <double, N> best_member_NEW = {};
        vector <vector <double>> err_members_NEW(num_err_members, vector <double>(N));
        vector <double> x_vals;
        x_vals.reserve(N);

        for (int x_index = 0; x_index < N; x_index++)
        {
            double x = pow(10, log10x_min + log10x_step * x_index);
            x_vals.push_back(x);
        }

        for (int member_id = 0; member_id < num_err_members + 1; member_id++)
        {
            const PDF* pdf = mkPDF(PDF_set, member_id);
        
            for (int x_index = 0; x_index < N; x_index++)
            {
                double x = pow(10, log10x_min + log10x_step * x_index);

                double xPDF_val = pdf->xfxQ(flavor, x, Q);
                
                if (member_id == 0)
                {
                    PDF_best_OLD[x_index] = xPDF_val;
                }
                else
                {
                    best_member_NEW[x_index] += omega[member_id - 1] * xPDF_val;
                    err_members_NEW[member_id - 1][x_index] = omega[member_id - 1] * xPDF_val;
                }
            }
        }

        for (int x_index = 0; x_index < N; x_index++)
        {
            best_member_NEW[x_index] = 1. / (num_err_members * 1.) * best_member_NEW[x_index];
        }

        vector <double> xPDF_sum_in_error_formula(N);

        for (int x_index = 0; x_index < N; x_index++)
        {
            for (int member = 0; member < num_err_members; member++)
            {
                xPDF_sum_in_error_formula[x_index] += omega[member] *
                    pow(err_members_NEW[member][x_index] / omega[member] - best_member_NEW[x_index], 2);
            }
        }

        for (int x_index = 0; x_index < N; x_index++)
        {
            xPDF_sum_in_error_formula[x_index] = sqrt(xPDF_sum_in_error_formula[x_index]);
        }

        vector <double> xPDF_error_NEW(N);

        for (int x_index = 0; x_index < N; x_index++)
        {
            xPDF_error_NEW[x_index] = sqrt(1. / (num_err_members * 1.) * xPDF_sum_in_error_formula[x_index]);
        }

        string filename = "output/old_PDF_vals/" + PDF_set + "/" + which_cross_sections_included + "/flavor_" + to_string(flavor) + ".txt";
        ofstream outfile1(filename, ios::out);

        for (int x_index = 0; x_index < N; x_index++)
        {
            outfile1 << x_vals[x_index];
            if (x_index != N - 1)
            {
                outfile1 << ",";
            }
        }
        outfile1 << endl;
        for (int x_index = 0; x_index < N; x_index++)
        {
            outfile1 << PDF_best_OLD[x_index];
            if (x_index != N - 1)
            {
                outfile1 << ",";
            }
        }
        outfile1.close();

        filename = "output/new_PDF_vals/" + PDF_set + "/" + which_cross_sections_included + "/flavor_" + to_string(flavor) + "_best.txt";
        ofstream outfile2(filename, ios::out);

        for (int x_index = 0; x_index < N; x_index++)
        {
            outfile2 << x_vals[x_index];
            if (x_index != N - 1)
            {
                outfile2 << ",";
            }
        }
        outfile2 << endl;
        for (int x_index = 0; x_index < N; x_index++)
        {
            outfile2 << best_member_NEW[x_index];
            if (x_index != N - 1)
            {
                outfile2 << ",";
            }
        }
        outfile2.close();
        
        filename = "output/new_PDF_vals/" + PDF_set + "/" + which_cross_sections_included + "/flavor_" + to_string(flavor) + "_errmem.txt";
        ofstream outfile2err(filename, ios::out);

        for (int member = 0; member < num_err_members; member++)
        {
            for (int x_index = 0; x_index < N; x_index++)
            {
                outfile2err << err_members_NEW[member][x_index];
                if (x_index != N - 1)
                {
                    outfile2err << ",";
                }
            }

            if (member != num_err_members - 1)
            {
                outfile2err << endl;
            }
        }
        outfile2err.close();
    }

    vector <double> Rcpm_best_OLD = load_1D_double_array_from_txt("input/theory_values/MC/best/" + 
                                        kinematic_quantity_for_new_Rcpm + "_" + 
                                        which_cross_sections_included + "_" + PDF_set + "_best.txt");
    vector <vector <double>> Rcpm_members_OLD = load_2D_double_array_from_txt("input/theory_values/MC/variation/" + 
                                                    kinematic_quantity_for_new_Rcpm + "_" + 
                                                    which_cross_sections_included + "_" + PDF_set + "_vals.txt", ',');


    vector <double> Rcpm_best_NEW(N_bins_in_kinematic_quantity);
    vector <double> Rcpm_sum_in_error_formula(N_bins_in_kinematic_quantity);

    
    for (int eta_lept_index = 0; eta_lept_index < N_bins_in_kinematic_quantity; eta_lept_index++)
    {
        for (int member = 0; member < num_err_members; member++)
        {
            Rcpm_best_NEW[eta_lept_index] += 1. / num_err_members * omega[member] * Rcpm_members_OLD[member][eta_lept_index];
            
        }
    }
    for (int eta_lept_index = 0; eta_lept_index < N_bins_in_kinematic_quantity; eta_lept_index++)
    {
        for (int member = 0; member < num_err_members; member++)
        {
            Rcpm_sum_in_error_formula[eta_lept_index] += omega[member] *
                pow(Rcpm_members_OLD[member][eta_lept_index] - Rcpm_best_NEW[eta_lept_index], 2);
        }
    }

    vector <double> Rcpm_error_NEW(N_bins_in_kinematic_quantity);

    for (int eta_lept_index = 0; eta_lept_index < N_bins_in_kinematic_quantity; eta_lept_index++)
    {
        Rcpm_error_NEW[eta_lept_index] = sqrt(1. / (num_err_members * 1.) * Rcpm_sum_in_error_formula[eta_lept_index]);
    }

    string filename = "output/new_Rcpm_" + kinematic_quantity_for_new_Rcpm + "/" + PDF_set + "/" + which_cross_sections_included + "/best.txt";
    ofstream outfile3(filename, ios::out);

    for (int bin_index = 0; bin_index < N_bins_in_kinematic_quantity; bin_index++)
    {
        outfile3 << Rcpm_best_NEW[bin_index];
        if (bin_index != N_bins_in_kinematic_quantity - 1)
        {
            outfile3 << endl;
        }
    }
    outfile3.close();

    filename = "output/new_Rcpm_" + kinematic_quantity_for_new_Rcpm + "/" + PDF_set + "/" + which_cross_sections_included + "/error_up.txt";
    ofstream outfile4(filename, ios::out);

    for (int bin_index = 0; bin_index < N_bins_in_kinematic_quantity; bin_index++)
    {
        outfile4 << Rcpm_error_NEW[bin_index];
        if (bin_index != N_bins_in_kinematic_quantity - 1)
        {
            outfile4 << endl;
        }
    }
    outfile4.close();

    filename = "output/new_Rcpm_" + kinematic_quantity_for_new_Rcpm + "/" + PDF_set + "/" + which_cross_sections_included + "/error_down.txt";
    ofstream outfile5(filename, ios::out);

    for (int bin_index = 0; bin_index < N_bins_in_kinematic_quantity; bin_index++)
    {
        outfile5 << Rcpm_error_NEW[bin_index];
        if (bin_index != N_bins_in_kinematic_quantity - 1)
        {
            outfile5 << endl;
        }
    }
    outfile5.close();
}

int main()
{
    LHAPDF::setVerbosity(0);

    INIReader reader("config.ini");

    if (reader.ParseError() < 0) {
        cout << "ERROR IN INITIALIZATION: CANNOT FIND 'config.ini'." << endl;
        exit(1);
    }

    string PDF_set = reader.Get("settings", "PDF_set", "none");
    if (PDF_set != "CT18ANLO" && PDF_set != "MSHT20nlo_as118" && PDF_set != "NNPDF40_nlo_pch_as_01180")
    {
        cout << "ERROR IN INITIALIZATION: PDF SET '" << PDF_set << "' IS NOT SUPPORTED." << endl;
        exit(1);
    }

    string flavors_raw = reader.Get("settings", "flavors", "");
    vector <int> flavors = split_to_ints(flavors_raw, ' ');

    string which_cross_sections_included = reader.Get("settings", "which_cross_sections_included", "none");
    if (which_cross_sections_included != "both" && which_cross_sections_included != "D" && which_cross_sections_included != "Dstar")
    {
        cout << "ERROR IN INITIALIZATION: 'which_cross_sections_included' cannot have the following value: " << which_cross_sections_included << endl;
        exit(1);
    }

    string kinematic_quantity_for_new_Rcpm = reader.Get("settings", "kinematic_quantity_for_new_Rcpm", "none");
    if (kinematic_quantity_for_new_Rcpm != "pTD" && kinematic_quantity_for_new_Rcpm != "eta_lept")
    {
        cout << "ERROR IN INITIALIZATION: 'kinematic_quantity_for_new_Rcpm' cannot have the following value: " << kinematic_quantity_for_new_Rcpm << endl;
        exit(1);
    }

    int num_err_members;

    string PDF_type;

    if (PDF_set == "CT18ANLO")
    {
        num_err_members = 58;
        PDF_type = "HESSIAN";
    }
    else if (PDF_set == "NNPDF40_nlo_pch_as_01180")
    {
        num_err_members = 100;
        PDF_type = "MC";
    }
    else if (PDF_set == "MSHT20nlo_as118")
    {
        num_err_members = 64;
        PDF_type = "HESSIAN";
    }
    else
    {
        cout << "ERROR: the PDF set " << PDF_set << " is not supported." << endl;
        exit(1);
    }

    if (PDF_type == "HESSIAN")
    {
        HESSIAN(num_err_members, flavors, PDF_set, which_cross_sections_included, kinematic_quantity_for_new_Rcpm, 5);
    }
    if (PDF_type == "MC")
    {
        MC(num_err_members, flavors, PDF_set, which_cross_sections_included, kinematic_quantity_for_new_Rcpm, 5);
    }

    return 0;
}