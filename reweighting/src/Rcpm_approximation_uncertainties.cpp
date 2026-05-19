#include "LHAPDF/LHAPDF.h"
#include <iostream>
#include <fstream>
#include <cstdlib>  // for setenv
using namespace LHAPDF;
using namespace std;

const double mu_fact2 = pow(80.385, 2);
const double epsilon = 0.222 * 0.222 / (0.975 * 0.975);

const double log10x_min = -4;
const double log10x_max = 0;
const double log10x_step = 0.001;
const int N = int((log10x_max - log10x_min) / log10x_step) + 1;

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

            double s = pdf->xfxQ2(3, x, mu_fact2);
            double sbar = pdf->xfxQ2(-3, x, mu_fact2);
            double d = pdf->xfxQ2(1, x, mu_fact2);
            double dbar = pdf->xfxQ2(-1, x, mu_fact2);

            if (member_id != 0)
            {
                err_member_vals[member_id - 1][x_index] = 1. - (epsilon * (d - dbar) + s - sbar) / s;
            }
            else
            {
                best_vals[x_index] = 1. - (epsilon * (d - dbar) + s - sbar) / s;
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

            if (PDF_set == "CT18ANLO")
            {
                err_plus_vals[x_index] /= 1.645;
                err_minus_vals[x_index] /= 1.645;
            }
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
                err_plus_vals[x_index] = err_plus_vals[x_index];
            }

            err_minus_vals[x_index] = err_plus_vals[x_index];
        }
    }
    

    string filename = "output/Rcpm_approximation_errors/" + PDF_set + "_OLD_plus.txt";

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

    filename = "output/Rcpm_approximation_errors/" + PDF_set + "_OLD_minus.txt";
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

    filename = "output/Rcpm_approximation_errors/" + PDF_set + "_OLD_mem_vals.txt";
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

int main()
{
    string PDF_set;
    cout << "Enter PDF set: ";
    cin >> PDF_set;

    string which_cross_sections_included;
    cout << "Which cross sections should be included? ";
    cin >> which_cross_sections_included;

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

    return 0;
}
