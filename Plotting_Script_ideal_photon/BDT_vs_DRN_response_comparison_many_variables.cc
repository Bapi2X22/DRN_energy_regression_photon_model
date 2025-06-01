#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <vector>
#include <iostream>
#include <iomanip>


using namespace std;

double Cruijff(double *x, double *par) {
    double mean   = par[0];
    double sigmaL = par[1];
    double sigmaR = par[2];
    double alphaL = par[3];
    double alphaR = par[4];
    double norm   = par[5];

    double dx = x[0] - mean;
    double sigma2;

    if (dx < 0)
        sigma2 = 2.0 * sigmaL * sigmaL + alphaL * dx * dx;
    else
        sigma2 = 2.0 * sigmaR * sigmaR + alphaR * dx * dx;

    if (sigma2 <= 0 || !std::isfinite(sigma2))
        return 0.0; // Protect against division by zero or NaNs

    double exponent = -dx * dx / sigma2;

    if (!std::isfinite(exponent)) 
        return 0.0; // Again, catch bad exponent values

    return norm * std::exp(exponent);
}

const double xmin = 0.95;
const double xmax = 1.05;


std::tuple<TGraphErrors*, TF1*, std::vector<double>, double> FitAndGraphHistogram(TH1F* hist, const std::string& name, int color, double xmin, double xmax) {
    // if (!hist) return std::make_tuple(nullptr, nullptr, std::vector<double>());

    hist->Scale(1.0 / hist->Integral());
    hist->Rebin(200);

    int nBins = hist->GetNbinsX();
    TGraphErrors* graph = new TGraphErrors(nBins);

    for (int j = 1; j <= nBins; j++) {
        double x = hist->GetBinCenter(j);
        double y = hist->GetBinContent(j);
        double xErr = hist->GetBinWidth(j) / 2.0;
        double yErr = hist->GetBinError(j);
        graph->SetPoint(j - 1, x, y);
        graph->SetPointError(j - 1, xErr, yErr);
    }

    graph->SetMarkerStyle(20); 
    graph->SetMarkerColor(color);
    graph->SetLineColor(color);


    TF1* fitFunc = new TF1(("fit_" + name).c_str(), Cruijff, xmin, xmax, 6);
    fitFunc->SetParLimits(1, 0.00001, 0.1);
    fitFunc->SetParLimits(2, 0.00001, 0.1);
    fitFunc->SetParameters(1.0, 0.01, 0.01, 0.1, 0.1, hist->GetMaximum());

    double old_mean = hist->GetMean();
    double old_sigma = 0;
    for (int iter = 0; iter < 400; ++iter) {
        hist->Fit(fitFunc, "RE");
        double mean = fitFunc->GetParameter(0);
        double sigmaL = fitFunc->GetParameter(1);
        double sigmaR = fitFunc->GetParameter(2);
        double sigma = (sigmaL + sigmaR) / 2.0;

        bool mean_converged = std::abs(old_mean - mean) / mean < 0.00002;
        bool sigma_converged = std::abs(old_sigma - sigma) / sigma < 0.00002;

        if (mean_converged && sigma_converged) break;

        old_mean = mean;
        old_sigma = sigma;
    }

    hist->GetXaxis()->SetRangeUser(xmin, xmax);

    std::vector<double> params;
    for (int i = 0; i < 6; ++i)
        params.push_back(fitFunc->GetParameter(i));

    double chi2 = fitFunc->GetChisquare();
    int ndf = fitFunc->GetNDF();
    double chi2_ndf = (ndf != 0) ? chi2 / ndf : 0;

    fitFunc->SetLineColor(color);
    fitFunc->SetLineStyle(0);

    // ---------------------- DRAWING SECTION ----------------------
    TCanvas* c1 = new TCanvas(("c_" + name).c_str(), ("Fit: " + name).c_str(), 800, 600);
    c1->cd();

    double yMax = hist->GetMaximum();

    graph->GetYaxis()->SetRangeUser(0, yMax * 1.25);

    graph->Draw("AP");
    fitFunc->Draw("same");
    graph->GetXaxis()->SetLimits(xmin, xmax);
    graph->GetXaxis()->SetTitle("E_{corr}/E_{gen}");
    graph->GetYaxis()->SetTitle("Normalized");

    TLatex latex;
    latex.SetTextSize(0.025);
    latex.SetTextColor(color);

    double xText = xmin + (xmax - xmin) * 0.05;
    // double yMax = hist->GetMaximum();
    double yBase = yMax * 1.2;
    double yStep = yMax * 0.07;

    latex.DrawLatex(xText, yBase - 0 * yStep, Form("#mu = %.4f, #sigma_{L} = %.4f, #sigma_{R} = %.4f",
                                                    fitFunc->GetParameter(0),
                                                    fitFunc->GetParameter(1),
                                                    fitFunc->GetParameter(2)));

    latex.DrawLatex(xText, yBase - 1 * yStep, Form("#chi^{2}/NDF = %.4f", chi2_ndf));

    c1->SaveAs(Form("Fit_%s.png", name.c_str()));
        // ------------------------------------------------------------

    return std::make_tuple(graph, fitFunc, params, chi2_ndf);


}


// void GenerateGraphsBDT_DRN(const std::string& file,
//     const std::vector<std::pair<int, int>>& bins,
//     std::vector<TGraphErrors*>& graphs_mu_BDT,
//     std::vector<TGraphErrors*>& graphs_sigma_BDT,
//     std::vector<TGraphErrors*>& graphs_mu_DRN,
//     std::vector<TGraphErrors*>& graphs_sigma_DRN) {

// auto ProcessFile = [&](const std::string& fileName,
// std::vector<TGraphErrors*>& graphs_mu_out,
// std::vector<TGraphErrors*>& graphs_sigma_out) {


// TFile* file = TFile::Open(fileName.c_str(), "READ");
// if (!file || file->IsZombie()) {
// std::cerr << "Error: Cannot open " << fileName << std::endl;
// return;
// }


// int num_bins = bins.size();
// TGraphErrors* g_mu = new TGraphErrors(num_bins);
// TGraphErrors* g_sigma = new TGraphErrors(num_bins);

// int idx = 0;
// for (const auto& bin : bins) {
//         // std::string histName = "eta_bin_B_" + std::to_string(bin.first) + "_" + std::to_string(bin.second);
//         std::string histName = "pT_bin_B_" + std::to_string(bin.first) + "_" + std::to_string(bin.second);
//         TH1F *hist = (TH1F*) file->Get(histName.c_str());
//         if (!hist) continue;

//         auto [graph, fit, params, chi2ndf] = FitAndGraphHistogram(hist, histName, kRed, xmin, xmax);
//         double mu_B = params[0];
//         double sigma_B = (params[1] + params[2]) / 2.0;
//         double x_B = (bin.first + bin.second) / 2.0;
//         double width_B = (bin.second - bin.first) / 2.0;

//         g_mu->SetPoint(idx, x, mu);
//         g_mu->SetPointError(idx, width, 0);

//         g_sigma->SetPoint(idx, x, sigma / mu);
//         g_sigma->SetPointError(idx, width, 0);
//         ++idx;

//         std::string histName = "pT_bin_D_" + std::to_string(bin.first) + "_" + std::to_string(bin.second);
//         TH1F *hist = (TH1F*) file->Get(histName.c_str());
//         if (!hist) continue;

//         auto [graph, fit, params, chi2ndf] = FitAndGraphHistogram(hist, histName, isFile2 ? kMagenta : kBlue, xmin, xmax);
//         double mu = params[0];
//         double sigma = (params[1] + params[2]) / 2.0;
//         double x = (bin.first + bin.second) / 2.0;
//         double width = (bin.second - bin.first) / 2.0;

//         g_mu->SetPoint(idx, x, mu);
//         g_mu->SetPointError(idx, width, 0);

//         g_sigma->SetPoint(idx, x, sigma / mu);
//         g_sigma->SetPointError(idx, width, 0);
//         ++idx;
//     }
// }

// graphs_mu_out.push_back(g_mu);
// graphs_sigma_out.push_back(g_sigma);

// file->Close();
// };

// // Process first file for both BDT and DRN
// ProcessFile(file, true, false, graphs_mu_BDT, graphs_sigma_BDT);
// ProcessFile(file, false, true, graphs_mu_DRN, graphs_sigma_DRN);
// }

void ProcessGraphs(const std::string& fileName,
                   const std::string& tag, // "B" or "D"
                   const std::vector<std::pair<int, int>>& bins,
                   int color,
                   std::vector<TGraphErrors*>& graphs_mu_out,
                   std::vector<TGraphErrors*>& graphs_sigma_out) {

    TFile* file = TFile::Open(fileName.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file " << fileName << std::endl;
        return;
    }

    int numBins = bins.size();
    TGraphErrors* g_mu = new TGraphErrors(numBins);
    TGraphErrors* g_sigma = new TGraphErrors(numBins);

    for (size_t i = 0; i < bins.size(); ++i) {
        const auto& bin = bins[i];
        int binLow = bin.first;
        int binHigh = bin.second;
        double x = (binLow + binHigh) / 2.0;
        double width = (binHigh - binLow) / 2.0;

        std::string histName = "pT_bin_" + tag + "_" + std::to_string(binLow) + "_" + std::to_string(binHigh);
        TH1F* hist = (TH1F*) file->Get(histName.c_str());
        if (!hist) continue;

        auto [graph, fit, params, chi2ndf] = FitAndGraphHistogram(hist, histName, color, xmin, xmax);
        double mu = params[0];
        double sigma = (params[1] + params[2]) / 2.0;

        g_mu->SetPoint(i, x, mu);
        g_mu->SetPointError(i, width, 0);

        g_sigma->SetPoint(i, x, sigma / mu);
        g_sigma->SetPointError(i, width, 0);
    }

    graphs_mu_out.push_back(g_mu);
    graphs_sigma_out.push_back(g_sigma);
    file->Close();
}

// ProcessGraphs("your_file.root", "B", bins, kRed, graphs_mu_BDT, graphs_sigma_BDT);
// ProcessGraphs("your_file.root", "D", bins, kBlue, graphs_mu_DRN, graphs_sigma_DRN);

void PlotMuAndSigmaFromFiles(
    const std::string& fileName_p, const std::string& fileName_e,
    const std::vector<std::pair<int, int>>& bins, int i, const std::string& label, const std::string& quantity, double muMin, double muMax, double muPMin, double muPMax, double sigMin, double sigMax, double sigPMin, double sigPMax, const std::string& xAxisLabel, const std::string& detector
) {
    std::vector<TGraphErrors*> graphs_mu_BDT, graphs_sigma_BDT;
    std::vector<TGraphErrors*> graphs_mu_DRN, graphs_sigma_DRN;

    // GenerateGraphsFromTwoFiles(
    //     fileName_p, fileName_e,
    //     bins,
    //     graphs_mu_BDT, graphs_sigma_BDT,
    //     graphs_mu_DRN_1, graphs_sigma_DRN_1,
    //     graphs_mu_DRN_2, graphs_sigma_DRN_2
    // );

    ProcessGraphs("your_file.root", "B", bins, kRed, graphs_mu_BDT, graphs_sigma_BDT);
    ProcessGraphs("your_file.root", "D", bins, kBlue, graphs_mu_DRN, graphs_sigma_DRN);
    // === Create canvas for mu with ratio ===
    TString canvasNameMu = "c_mu";
    TCanvas* c_mu = new TCanvas(canvasNameMu, canvasNameMu, 800, 800);

    // Create pads
    TPad* pad1_mu = new TPad("pad1_mu", "pad1_mu", 0, 0.3, 1, 1.0);
    TPad* pad2_mu = new TPad("pad2_mu", "pad2_mu", 0, 0.0, 1, 0.3);
    pad1_mu->SetBottomMargin(0.02);
    pad2_mu->SetTopMargin(0.02);
    pad2_mu->SetBottomMargin(0.3);
    pad1_mu->SetLeftMargin(0.15); // Default is ~0.1
    pad2_mu->SetLeftMargin(0.15); // Default is ~0.1
    pad1_mu->Draw();
    pad2_mu->Draw();

    // === Upper pad for mu ===
    pad1_mu->cd();
    pad1_mu->SetGrid();

    graphs_mu_BDT[0]->SetLineColor(kRed);
    graphs_mu_BDT[0]->SetMarkerStyle(20);
    graphs_mu_BDT[0]->SetTitle("Electron_pT_Gun_Sample");
    graphs_mu_BDT[0]->SetMarkerColor(kRed);
    graphs_mu_BDT[0]->Draw("APL");
    graphs_mu_BDT[0]->GetYaxis()->SetRangeUser(muMin, muMax);
    graphs_mu_BDT[0]->GetYaxis()->SetTitle("Response(#mu)");
    graphs_mu_BDT[0]->GetXaxis()->SetLabelSize(0);

    graphs_mu_DRN[0]->SetMarkerColor(kBlue);
    graphs_mu_DRN[0]->SetLineColor(kBlue);
    graphs_mu_DRN[0]->SetMarkerStyle(21);
    graphs_mu_DRN[0]->Draw("PL SAME");

    TLatex* text_run3 = new TLatex();
    text_run3->SetNDC();                 // Use normalized coordinates
    text_run3->SetTextFont(62);          // Standard font
    text_run3->SetTextSize(0.035);       // Adjust size as needed
    text_run3->SetTextAlign(21);         // Horizontal center alignment
    text_run3->DrawLatex(0.5, 0.82, label.c_str());  // (x, y, text)


    // TLegend* leg_mu = new TLegend(0.6, 0.7, 0.88, 0.88);
    // leg_mu->AddEntry(graphs_mu_BDT[0], "BDT", "lp");
    // leg_mu->AddEntry(graphs_mu_DRN_1[0], "photon_DRN", "lp");
    // leg_mu->AddEntry(graphs_mu_DRN_2[0], "electron_DRN", "lp");
    // leg_mu->Draw();

    TLegend* leg_mu = new TLegend(0.65, 0.75, 0.85, 0.87); // Smaller box
    leg_mu->SetTextSize(0.025); // Smaller font size
    leg_mu->AddEntry(graphs_mu_BDT[0], "BDT", "lp");
    leg_mu->AddEntry(graphs_mu_DRN_1[0], "photon_DRN", "lp");
    leg_mu->AddEntry(graphs_mu_DRN_2[0], "electron_DRN", "lp");
    leg_mu->Draw();

    // === Lower pad for mu ratio ===
    pad2_mu->cd();
    int nPoints = graphs_mu_BDT[0]->GetN();
    TGraphErrors* ratio_BDT_mu = new TGraphErrors(nPoints);
    TGraphErrors* ratio_DRN1_mu = new TGraphErrors(nPoints);

    for (int k = 0; k < nPoints; ++k) {
        double x, y_BDT, y_DRN1, y_DRN2;
        double err_BDT, err_DRN1, err_DRN2;

        graphs_mu_BDT[0]->GetPoint(k, x, y_BDT);
        graphs_mu_DRN_1[0]->GetPoint(k, x, y_DRN1);
        graphs_mu_DRN_2[0]->GetPoint(k, x, y_DRN2);

        err_BDT = graphs_mu_BDT[0]->GetErrorY(k);
        err_DRN1 = graphs_mu_DRN_1[0]->GetErrorY(k);
        err_DRN2 = graphs_mu_DRN_2[0]->GetErrorY(k);

        ratio_BDT_mu->SetPoint(k, x, y_BDT / y_DRN2);
        ratio_BDT_mu->SetPointError(k, 0, (y_BDT / y_DRN2) * sqrt(pow(err_BDT / y_BDT, 2) + pow(err_DRN2 / y_DRN2, 2)));

        ratio_DRN1_mu->SetPoint(k, x, y_DRN1 / y_DRN2);
        ratio_DRN1_mu->SetPointError(k, 0, (y_DRN1 / y_DRN2) * sqrt(pow(err_DRN1 / y_DRN1, 2) + pow(err_DRN2 / y_DRN2, 2)));
    }

    ratio_BDT_mu->SetMarkerColor(kRed);
    ratio_BDT_mu->SetLineColor(kRed);
    ratio_BDT_mu->SetMarkerStyle(20);
    // ratio_BDT_mu->SetTitle(";#eta;Ratio to DRN2");
    ratio_BDT_mu->SetTitle(( ";" + xAxisLabel + ";Ratio to DRN2" ).c_str());


    ratio_DRN1_mu->SetMarkerColor(kBlue);
    ratio_DRN1_mu->SetLineColor(kBlue);
    ratio_DRN1_mu->SetMarkerStyle(21);

    ratio_BDT_mu->Draw("APL");
    ratio_BDT_mu->GetYaxis()->SetTitleSize(0.07);
    ratio_BDT_mu->GetYaxis()->SetTitleOffset(0.8);
    ratio_BDT_mu->GetYaxis()->SetLabelSize(0.07);
    ratio_BDT_mu->GetXaxis()->SetLabelSize(0.08);
    ratio_BDT_mu->GetXaxis()->SetTitleSize(0.1);
    ratio_BDT_mu->GetXaxis()->SetTitleOffset(1.0);
    ratio_BDT_mu->SetMinimum(muPMin);
    ratio_BDT_mu->SetMaximum(muPMax);

    ratio_DRN1_mu->Draw("PL SAME");

    // Save the mu canvas
    c_mu->SaveAs(Form("%s_vs_mu_comparison_set_%s_%d.png", quantity.c_str(), detector.c_str(), i));


    // === Repeat everything similarly for sigma ===

    TString canvasNameSigma = "c_sigma";
    TCanvas* c_sigma = new TCanvas(canvasNameSigma, canvasNameSigma, 800, 800);
    TPad* pad1_sigma = new TPad("pad1_sigma", "pad1_sigma", 0, 0.3, 1, 1.0);
    TPad* pad2_sigma = new TPad("pad2_sigma", "pad2_sigma", 0, 0.0, 1, 0.3);
    pad1_sigma->SetBottomMargin(0.02);
    pad2_sigma->SetTopMargin(0.02);
    pad2_sigma->SetBottomMargin(0.3);
    pad1_sigma->SetLeftMargin(0.15); // Default is ~0.1
    pad2_sigma->SetLeftMargin(0.15); // Default is ~0.1
    pad1_sigma->Draw();
    pad2_sigma->Draw();

    pad1_sigma->cd();
    pad1_sigma->SetGrid();

    graphs_sigma_BDT[0]->SetMarkerColor(kRed);
    graphs_sigma_BDT[0]->SetLineColor(kRed);
    graphs_sigma_BDT[0]->SetMarkerStyle(20);
    graphs_sigma_BDT[0]->SetTitle("Electron_pT_Gun_Sample");
    graphs_sigma_BDT[0]->Draw("APL");
    graphs_sigma_BDT[0]->GetYaxis()->SetRangeUser(sigMin, sigMax);
    graphs_sigma_BDT[0]->GetYaxis()->SetTitle("#sigma/#mu");
    graphs_sigma_BDT[0]->GetXaxis()->SetLabelSize(0);

    graphs_sigma_DRN[0]->SetMarkerColor(kBlue);
    graphs_sigma_DRN[0]->SetLineColor(kBlue);
    graphs_sigma_DRN[0]->SetMarkerStyle(21);
    graphs_sigma_DRN[0]->Draw("PL SAME");

    TLatex* text_run2 = new TLatex();
    text_run2->SetNDC();                 // Use normalized coordinates
    text_run2->SetTextFont(62);          // Standard font
    text_run2->SetTextSize(0.035);       // Adjust size as needed
    text_run2->SetTextAlign(21);         // Horizontal center alignment
    text_run2->DrawLatex(0.5, 0.82, label.c_str());  // (x, y, text)


    // TLegend* leg_sigma = new TLegend(0.6, 0.7, 0.88, 0.88);
    // leg_sigma->AddEntry(graphs_sigma_BDT[0], "BDT", "lp");
    // leg_sigma->AddEntry(graphs_sigma_DRN_1[0], "photon_DRN", "lp");
    // leg_sigma->AddEntry(graphs_sigma_DRN_2[0], "electron_DRN", "lp");
    // leg_sigma->Draw();

    // TLegend* leg_sigma = new TLegend(0.15, 0.75, 0.35, 0.87); // Left side, smaller box
    TLegend* leg_sigma = new TLegend(0.65, 0.75, 0.85, 0.87);
    leg_sigma->SetTextSize(0.025); // Smaller font size
    leg_sigma->AddEntry(graphs_sigma_BDT[0], "BDT", "lp");
    leg_sigma->AddEntry(graphs_sigma_DRN[0], "DRN", "lp");
    leg_sigma->Draw();


    pad2_sigma->cd();
    TGraphErrors* ratio_BDT_sigma = new TGraphErrors(nPoints);
    TGraphErrors* ratio_DRN1_sigma = new TGraphErrors(nPoints);

    for (int m = 0; m < nPoints; ++m) {
        double X, Y_BDT, Y_DRN1, Y_DRN2;
        double Err_BDT, Err_DRN1, Err_DRN2;

        graphs_sigma_BDT[0]->GetPoint(m, X, Y_BDT);
        graphs_sigma_DRN_1[0]->GetPoint(m, X, Y_DRN1);
        graphs_sigma_DRN_2[0]->GetPoint(m, X, Y_DRN2);

        Err_BDT = graphs_sigma_BDT[0]->GetErrorY(m);
        Err_DRN1 = graphs_sigma_DRN_1[0]->GetErrorY(m);
        Err_DRN2 = graphs_sigma_DRN_2[0]->GetErrorY(m);

        ratio_BDT_sigma->SetPoint(m, X, Y_BDT / Y_DRN2);
        ratio_BDT_sigma->SetPointError(m, 0, (Y_BDT / Y_DRN2) * sqrt(pow(Err_BDT / Y_BDT, 2) + pow(Err_DRN2 / Y_DRN2, 2)));

        ratio_DRN1_sigma->SetPoint(m, X, Y_DRN1 / Y_DRN2);
        ratio_DRN1_sigma->SetPointError(m, 0, (Y_DRN1 / Y_DRN2) * sqrt(pow(Err_DRN1 / Y_DRN1, 2) + pow(Err_DRN2 / Y_DRN2, 2)));
    }

    ratio_BDT_sigma->SetMarkerColor(kRed);
    ratio_BDT_sigma->SetLineColor(kRed);
    ratio_BDT_sigma->SetMarkerStyle(20);
    ratio_BDT_sigma->SetTitle(( ";" + xAxisLabel + ";Ratio to DRN2" ).c_str());

    ratio_DRN1_sigma->SetMarkerColor(kBlue);
    ratio_DRN1_sigma->SetLineColor(kBlue);
    ratio_DRN1_sigma->SetMarkerStyle(21);

    ratio_BDT_sigma->Draw("APL");
    ratio_BDT_sigma->GetYaxis()->SetTitleSize(0.07);
    ratio_BDT_sigma->GetYaxis()->SetTitleOffset(0.6);
    ratio_BDT_sigma->GetYaxis()->SetLabelSize(0.07);
    ratio_BDT_sigma->GetXaxis()->SetLabelSize(0.08);
    ratio_BDT_sigma->GetXaxis()->SetTitleSize(0.1);
    ratio_BDT_sigma->GetXaxis()->SetTitleOffset(1.0);
    ratio_BDT_sigma->SetMinimum(sigPMin);
    ratio_BDT_sigma->SetMaximum(sigPMax);

    ratio_DRN1_sigma->Draw("PL SAME");

    // Save the sigma canvas
    c_sigma->SaveAs(Form("%s_vs_sigma_comparison_set_%s_%d.png", quantity.c_str(), detector.c_str(), i));
}
