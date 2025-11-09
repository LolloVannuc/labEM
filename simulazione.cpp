#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TRandom.h>
#include <iostream>
#include <vector>
#include <cmath>
#include "TApplication.h"
#include <TRandom3.h>

int main(int argc, char **argv)
{

    TApplication app("app", &argc, argv);

    // VARIABILI PRINCIPALI

    int Nevents = 1000000;
    int Nbins = 100;
    int NRipetizioni = 200;

    double xmin = 0.;
    double xmax = 2.;

    double k_medio = 5.2;
    double phi_medio = 1.8;
    double b_medio = 0.2;

    double k_perc = 0.02;
    double phi_perc = 0.05;
    double b_perc = 0.01;

    double k_inc = k_medio * k_perc;
    double phi_inc = phi_medio * phi_perc;
    double b_inc = b_medio * b_perc;

    // DEFINIZIONE DELLE VARIABILI PRINCIPALI

    TH1F *h1 = new TH1F("h1", "h1", Nbins, xmin, xmax);
    TF1 *f1 = new TF1("simulazione", "cos([0] * x + [1])^2 + [2]", xmin, xmax);
    TH1F *h_ideale = new TH1F("h_ideale", "Valore teorico", Nbins, xmin, xmax);

    f1->SetNpx(500);

    // Vectori per somma e somma quadratica per bin 3.2

    std::vector<double> somma_per_bin(Nbins, 0.0);
    std::vector<double> somma_quadratica_per_bin(Nbins, 0.0);

    // 1. Nuovi accumulatori per il metodo 3.3
    std::vector<double> somma_smeering(Nbins, 0.0);
    std::vector<double> somma_q_smeering(Nbins, 0.0);

    // SCALATURA DELLA FUNZIONE (ora è sbagliata perchè i parametri variano)

    double integral_f = f1->Integral(xmin, xmax);
    double integral_h = h1->Integral();

    double scale_factor_area = integral_h / integral_f;

    double bin_width = h1->GetBinWidth(1);

    // Il fattore di scala finale deve includere la larghezza del bin per matchare il contenuto del bin
    double final_scale_factor = scale_factor_area * bin_width;

    // per definire il fattore di scala uso tstring perche così è come se fosse costante
    TString formula_base = "cos([0] * x + [1])^2 + [2]";
    TString scaled_formula_str = TString::Format("%f * (%s)", final_scale_factor, formula_base.Data());

    TF1 *f_scaled = new TF1("f_scaled", scaled_formula_str.Data(), 0., 2.);

    f_scaled->SetParameters(5.2, 1.8, 0.2);

    // CICLO DI RIPETIZIONI

    for (int rip = 0; rip < NRipetizioni; ++rip)
    {

        // GENERAZIONE PARAMETRI

        double k = gRandom->Gaus(k_medio, k_inc);
        double phi = gRandom->Gaus(phi_medio, phi_inc);
        double b = gRandom->Gaus(b_medio, b_inc);

        f1->SetParameters(k, phi, b);

        double integral_fluct = f1->Integral(xmin, xmax);

        // Riempimento dell'istogramma ideale con i valori teorici

        for (int bin = 1; bin <= Nbins; ++bin)
        {
            double x_center = h_ideale->GetBinCenter(bin);
            double f_value = f1->Eval(x_center);
            double bin_width = h_ideale->GetBinWidth(bin);
            double expected_count = f_value * Nevents * bin_width / integral_fluct; // Normalizzazione
            h_ideale->SetBinContent(bin, expected_count);
        }

        // RIEMPIMENTO ISTOGRAMMA

        for (int i = 0; i < Nevents; ++i)
        {
            double x = f1->GetRandom();
            h1->Fill(x);
        }

        // AGGIORNAMENTO SOMME PER BIN

        double SQ = 0.; // serve per 3.1

        for (int bin = 1; bin <= Nbins; ++bin)
        {

            // AGGIORNAMENTO PER IL 3.2

            double content = h1->GetBinContent(bin);
            int vector_index = bin - 1;
            somma_per_bin[vector_index] += content;
            somma_quadratica_per_bin[vector_index] += content * content;

            // AGGIORNAMENTO PER IL 3.3

            double expected_count = h_ideale->GetBinContent(bin);

            double sigma = (expected_count > 0) ? TMath::Sqrt(expected_count) : 0.0;

            double fluctuated_count = gRandom->Gaus(expected_count, sigma);

            somma_smeering[vector_index] += fluctuated_count;
            somma_q_smeering[vector_index] += fluctuated_count * fluctuated_count;

            // CALCOLO SCARTO QUADRATICO PER IL 3.1
            double diff{0.};

            diff = content - expected_count;
            SQ += diff * diff;
        }

        double RMS = TMath::Sqrt(SQ / Nbins); // RMS per il 3.1

        // CLEANING
        h1->Reset();
    }

    // ELABORAZIONE RISULTATI E OUTPUT

    for (int bin = 1; bin <= Nbins; ++bin)
    {
        int i = bin - 1; // Indice vettore

        // Calcolo finale 3.2
        double media_3_2 = somma_per_bin[i] / NRipetizioni;
        double media_q_3_2 = somma_quadratica_per_bin[i] / NRipetizioni;
        double varianza_3_2 = media_q_3_2 - (media_3_2 * media_3_2);
        double std_dev_3_2 = TMath::Sqrt(varianza_3_2);

        // Calcolo finale 3.3
        double media_3_3 = somma_smeering[i] / NRipetizioni;
        double media_q_3_3 = somma_q_smeering[i] / NRipetizioni;
        double varianza_3_3 = media_q_3_3 - (media_3_3 * media_3_3);
        double std_dev_3_3 = TMath::Sqrt(varianza_3_3);

        // Stampa

        int bin_limit = 10; // Numero di bin da stampare

        if (bin <= bin_limit)
        {
            std::cout << "--- Bin: " << bin << " ---" << std::endl;
            std::cout << "  3.2 (Rigeneraz): Media = " << media_3_2
                      << ", Dev. Std = " << std_dev_3_2 << std::endl;
            std::cout << "  3.3 (Smeering):  Media = " << media_3_3
                      << ", Dev. Std = " << std_dev_3_3 << std::endl;
        }
    }

    // FIT

    TH1F *h_fit = new TH1F("h_fit", "h_fit", Nbins, xmin, xmax);

    f1->SetParameters(k_medio, phi_medio, b_medio);

    TF1 *f_fit = new TF1("f_fit", "[3] * (cos([0] * x + [1])^2 + [2])", xmin, xmax);
    f_fit->SetNpx(500);
    double integral_pdf = f1->Integral(xmin, xmax);
    double norm_factor_guess = (Nevents * bin_width) / integral_pdf;


    f_fit->SetParameters(k_medio, phi_medio, b_medio, norm_factor_guess);

    for (int i = 0; i < Nevents; ++i)
    {
        double x = f1->GetRandom();
        h_fit->Fill(x);
    }
    // NORMALIZZAZIONE PER IL FIT

    
    h_fit->Fit(f_fit, "R");

    // CALCOLO SCARTI

    std::cout << "\n--- Calcolo Scarti (Residui - Punto 4) ---" << std::endl;
    TH1F *h_residui = new TH1F("h_residui", "Residui (Dati - Fit)", Nbins, xmin, xmax);

    for (int bin = 1; bin <= Nbins; ++bin)
    {
        double x_center = h_fit->GetBinCenter(bin);
        double dati_content = h_fit->GetBinContent(bin);
        double fit_value = f_fit->Eval(x_center); 

        double residuo = dati_content - fit_value;

        h_residui->SetBinContent(bin, residuo);

        // Stampa i primi 5 punti
        if (bin <= 5)
        {
            std::cout << "  Bin " << bin << ": Dati=" << dati_content
                      << ", Fit=" << fit_value
                      << ", Scarto=" << residuo << std::endl;
        }
    }

    
    // DRAWING HISTO E FUNCTION
    /*TCanvas *c1 = new TCanvas ("simulazione", "simulazione", 200, 200);

   f_scaled->SetLineColor(kRed);
   f_scaled->SetTitle("Istogramma e Funzione Scalata per Confronto");

   c1 -> Divide(2,1);
   c1 -> cd(1);
   f1 -> Draw();

   c1 -> cd(2);
   h1-> Draw();
   f_scaled -> Draw("SAME");



   app.Run();
   */

    return 0;
}
