#include "LHAPDF/LHAPDF.h"
#include <iostream>
#include <fstream>
#include <cstdlib>  // for setenv
using namespace LHAPDF;
using namespace std;

const double mu_fact2 = pow(80.385, 2);

const double log10x_min = -5;
const double log10x_max = 0;
const double log10x_step = 0.001;
const int N = int((log10x_max - log10x_min) / log10x_step) + 1;

const int flavor = 3;

vector <double> load_1D_double_array_from_txt(string filename, int skiprows)
{
    ifstream file(filename);

    if (!file.is_open())
    {
        cout << "Error in the function 'load_1D_double_array_from_txt': could not open the following file: " << filename << endl;
        exit(1);
    }

    vector<double> values;
    string line;

    int row = 0;

    while (getline(file, line))
    {
        if (row > skiprows - 1)
        {
            stringstream ss(line);
            string value;

            while (getline(ss, value, ','))
            {
                values.push_back(stod(value));
            }
        }
        row += 1;
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

int OLD(string PDF_set, int num_err_members)
{
    vector <double> log10x_vals;

    array <double, N> best_vals;
    vector <vector <double>> err_member_vals;
    err_member_vals.resize(num_err_members, vector <double> (N));

    for (int member_id = 0; member_id < num_err_members + 1; member_id++)
    {
        const PDF* pdf = mkPDF(PDF_set, member_id);
    
        for (int x_index = 0; x_index < N; x_index++)
        {
            double x = pow(10, log10x_min + log10x_step * x_index);

            double xs_val = pdf->xfxQ2(flavor, x, mu_fact2);
            double xsbar_val = pdf->xfxQ2(-flavor, x, mu_fact2);

            if (member_id != 0)
            {
                err_member_vals[member_id - 1][x_index] = xs_val - xsbar_val;
            }
            else
            {
                best_vals[x_index] = xs_val - xsbar_val;
            }
        }
    }

    array <double, N> err_plus_vals;
    array <double, N> err_minus_vals;

    if (PDF_set == "CT18ANLO" || PDF_set == "MSHT20nlo_as118")
    {
        for (int x_index = 0; x_index < N; ++x_index)
        {
            double sum_plus_sq = 0.0;
            double sum_minus_sq = 0.0;

            for (int pair = 1; pair <= num_err_members / 2; ++pair)
            {
                double e1 = err_member_vals[2 * (pair - 1)][x_index];
                double e2 = err_member_vals[2 * (pair - 1) + 1][x_index];
                double b  = best_vals[x_index];

                double d1 = e1 - b;
                double d2 = e2 - b;

                double plus = max(max(d1, d2), 0.0);
                double minus = max(max(-d1, -d2), 0.0); // negative deviations -> positive magnitude for minus

                sum_plus_sq  += plus * plus;
                sum_minus_sq += minus * minus;
            }

            err_plus_vals[x_index]  = sqrt(sum_plus_sq);
            err_minus_vals[x_index] = sqrt(sum_minus_sq);
        }
    }
    else
    {
        array <double, N> sum_in_error_formula = {};

        for (int x_index = 0; x_index < N; x_index++)
        {
            for (int member = 0; member < num_err_members; member++)
            {
                sum_in_error_formula[x_index] += pow(err_member_vals[member][x_index] - best_vals[x_index], 2);
            }
        }

        for (int x_index = 0; x_index < N; x_index++)
        {
            err_plus_vals[x_index] = sqrt(1. / (num_err_members * 1.) * sum_in_error_formula[x_index]);

            if (PDF_set == "NNPDF40_nlo_pch_as_01180")
            {
                err_plus_vals[x_index] = err_plus_vals[x_index] * 1.645;
            }

            err_minus_vals[x_index] = err_plus_vals[x_index];
        }
    }
    

    string filename = "output/strangeness_asymmetry_errors/" + PDF_set + "_OLD_plus.txt";

    cout << filename << endl;
    ofstream outfile1(filename, ios::out);

    for (int x_index = 0; x_index < N; x_index++)
    {
        outfile1 << err_plus_vals[x_index];
        if (x_index != N - 1)
        {
            outfile1 << ",";
        }
    }
    outfile1.close();

    filename = "output/strangeness_asymmetry_errors/" + PDF_set + "_OLD_minus.txt";
    ofstream outfile2(filename, ios::out);

    for (int x_index = 0; x_index < N; x_index++)
    {
        outfile2 << err_minus_vals[x_index];
        if (x_index != N - 1)
        {
            outfile2 << ",";
        }
    }
    outfile2.close();

    filename = "output/strangeness_asymmetry_errors/" + PDF_set + "_OLD_mem_vals.txt";
    ofstream outfile3(filename, ios::out);

    for (int mem = 0; mem < num_err_members; mem++)
    {
        for (int x_index = 0; x_index < N; x_index++)
        {
            outfile3 << err_member_vals[mem][x_index];
            if (x_index != N - 1)
            {
                outfile3 << ",";
            }
        }
        if (mem != num_err_members - 1)
        {
            outfile3 << endl;
        }
    }
    outfile3.close();

    return 0;
}

int NEW(string PDF_set, int num_err_members, string which_cross_sections_included)
{
    vector <double> xs_best = load_1D_double_array_from_txt("output/new_PDF_vals/" + PDF_set + "/" +
                which_cross_sections_included + "/flavor_" + to_string(flavor) + "_best.txt", 1);
    vector <double> xsbar_best = load_1D_double_array_from_txt("output/new_PDF_vals/" + PDF_set + "/" +
                which_cross_sections_included + "/flavor_" + to_string(-flavor) + "_best.txt", 1);

    array <double, N> err_plus_vals = {};
    array <double, N> err_minus_vals = {};

    if (PDF_set == "CT18ANLO" || PDF_set == "MSHT20nlo_as118")
    {
        vector <vector <double>> xs_plus = load_2D_double_array_from_txt("output/new_PDF_vals/" + PDF_set + "/" + 
                which_cross_sections_included + "/flavor_" + to_string(flavor) + "_plus.txt", ',');
        vector <vector <double>> xs_minus = load_2D_double_array_from_txt("output/new_PDF_vals/" + PDF_set + "/" +
                    which_cross_sections_included + "/flavor_" + to_string(flavor) + "_minus.txt", ',');
        vector <vector <double>> xsbar_plus = load_2D_double_array_from_txt("output/new_PDF_vals/" + PDF_set + "/" +
                    which_cross_sections_included + "/flavor_" + to_string(-flavor) + "_plus.txt", ',');
        vector <vector <double>> xsbar_minus = load_2D_double_array_from_txt("output/new_PDF_vals/" + PDF_set + "/" +
                    which_cross_sections_included + "/flavor_" + to_string(-flavor) + "_minus.txt", ',');

        for (int x_index = 0; x_index < N; ++x_index)
        {
            double sum_plus_sq = 0.0;
            double sum_minus_sq = 0.0;

            for (int pair = 0; pair < num_err_members / 2; pair++)
            {
                double e1 = xs_plus[pair][x_index] - xsbar_plus[pair][x_index];
                double e2 = xs_minus[pair][x_index] - xsbar_minus[pair][x_index];
                double b  = xs_best[x_index] - xsbar_best[x_index];

                double d1 = e1 - b;
                double d2 = e2 - b;

                double plus = max(max(d1, d2), 0.0);
                double minus = max(max(-d1, -d2), 0.0);

                sum_plus_sq  += plus * plus;
                sum_minus_sq += minus * minus;
            }

            err_plus_vals[x_index]  = sqrt(sum_plus_sq);
            err_minus_vals[x_index] = sqrt(sum_minus_sq);
        }
    }
    else
    {
        vector <vector <double>> xs_errmem = load_2D_double_array_from_txt("output/new_PDF_vals/" + PDF_set + "/" + 
                which_cross_sections_included + "/flavor_" + to_string(flavor) + "_errmem.txt", ',');
        vector <vector <double>> xsbar_errmem = load_2D_double_array_from_txt("output/new_PDF_vals/" + PDF_set + "/" + 
                which_cross_sections_included + "/flavor_" + to_string(-flavor) + "_errmem.txt", ',');

        vector <double> omega(num_err_members);

        ifstream file("output/omega.txt");
        if (!file)
        {
            cerr << "Error opening file" << endl;
            exit(1);
        }

        for (int i = 0; i < num_err_members; ++i)
        {
            file >> omega[i];
        }

        array <double, N> sum_in_error_formula = {};

        for (int x_index = 0; x_index < N; x_index++)
        {
            for (int member = 0; member < num_err_members; member++)
            {
                sum_in_error_formula[x_index] += omega[member] *
                    pow((xs_errmem[member][x_index] - xsbar_errmem[member][x_index]) / omega[member] -
                    (xs_best[x_index] - xsbar_best[x_index]), 2);
            }
        }

        for (int x_index = 0; x_index < N; x_index++)
        {
            if (sum_in_error_formula[x_index] < 0)
            {
                cout << sum_in_error_formula[x_index] << endl;
            }
            err_plus_vals[x_index] = sqrt(1. / (num_err_members * 1.) * sum_in_error_formula[x_index]);

            if (PDF_set == "NNPDF40_nlo_pch_as_01180")
            {
                err_plus_vals[x_index] = err_plus_vals[x_index] * 1.645;
            }

            err_minus_vals[x_index] = err_plus_vals[x_index];
        }
        
    }

    string filename = "output/strangeness_asymmetry_errors/" + PDF_set + "_" + which_cross_sections_included + "_NEW_plus.txt";

    ofstream outfile1(filename, ios::out);

    for (int x_index = 0; x_index < N; x_index++)
    {
        outfile1 << err_plus_vals[x_index];
        if (x_index != N - 1)
        {
            outfile1 << ",";
        }
    }
    outfile1.close();

    filename = "output/strangeness_asymmetry_errors/" + PDF_set + "_" + which_cross_sections_included + "_NEW_minus.txt";
    ofstream outfile2(filename, ios::out);

    for (int x_index = 0; x_index < N; x_index++)
    {
        outfile2 << err_minus_vals[x_index];
        if (x_index != N - 1)
        {
            outfile2 << ",";
        }
    }
    outfile2.close();

    
    return 0;
}

int main()
{
    string PDF_set;
    cout << "Enter PDF set: ";
    cin >> PDF_set;

    string which_cross_sections_included = "both";
    //cout << "Which cross sections should be included? ";
    //cin >> which_cross_sections_included;

    int num_err_members;

    if (PDF_set == "CT18ANLO")
    {
        num_err_members = 58;
    }
    else if (PDF_set == "MSHT20nlo_as118")
    {
        num_err_members = 64;
    }
    else if (PDF_set == "NNPDF40_nlo_pch_as_01180")
    {
        num_err_members = 100;
    }
    else
    {
        cout << "ERROR: the PDF set " << PDF_set << " is not supported." << endl;
        exit(1);
    }

    OLD(PDF_set, num_err_members);
    NEW(PDF_set, num_err_members, which_cross_sections_included);

    return 0;
}