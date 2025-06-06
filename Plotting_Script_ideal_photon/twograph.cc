#include <TFile.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TF1.h>
#include <TMath.h>
#include <vector>

// Cruijff function definition
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
void SaveGraphWithFitData(TGraphErrors* g, TF1* fitFunc,
                          double xMin, double xMax,
                          const char* outFile = "fit_data.txt") {
    std::ofstream file(outFile);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open output file.\n";
        return;
    }

    file << "# x\t y_data\t y_fit\t sigma_y\t chi2_i\n";

    int N = g->GetN();
    int nPar = fitFunc->GetNpar();

    double chi2_total = 0.0;
    int validPoints = 0;

    for (int i = 0; i < N; ++i) {
        double x, y, sigma_y;
        g->GetPoint(i, x, y);

        // Skip points outside the specified x-range
        if (x < xMin || x > xMax)
            continue;

        sigma_y = g->GetErrorY(i);
        if (sigma_y == 0) {
            std::cerr << "Warning: zero error at point " << i << ", skipping in chi2 calculation.\n";
            continue;
        }

        double y_fit = fitFunc->Eval(x);
        double chi2_i = ((y - y_fit) * (y - y_fit)) / (sigma_y * sigma_y);

        file << x << "\t" << y << "\t" << y_fit << "\t" << sigma_y << "\t" << chi2_i << "\n";

        chi2_total += chi2_i;
        ++validPoints;
    }

    double chi2_ndf = (validPoints > nPar) ? chi2_total / (validPoints - nPar) : 0.0;

    std::cout << "Total chi2 (in range): " << chi2_total << "\n";
    std::cout << "NDF (valid points - nPar): " << (validPoints - nPar) << "\n";
    std::cout << "Chi2 / NDF = " << chi2_ndf << std::endl;

    file.close();
    std::cout << "Saved data to " << outFile << std::endl;
}

    double xmin = 0.98;
    double xmax = 1.02;
void superimposePlotsWithFits() {
    // std::vector<std::string> fileNames = {"Plot_UND4.root", "Plot_UND2.root", "Plot_UND3.root", "Plot_UND1.root"};
    std::vector<std::string> fileNames = {"Plot_ideal_photon_all.root"};
    std::vector<Color_t> colors = {kRed, kBlue, kGreen+2, kBlack};
    // std::vector<std::string> labels = {"None Varied", "Rechit threshold varied", "Noise threshold varied", "Both varied"};
    std::vector<std::string> labels = {"BDT", "DRN"};
    // std::vector<int> Entries = {802294, 803486, 803240, 803498};
    
    
    TMultiGraph *mg = new TMultiGraph();
    // TMultiGraph *mg_DRN = new TMultiGraph();
    // TLegend *legend_BDT = new TLegend(0.7, 0.7, 0.9, 0.9);
    // TLegend *legend_DRN = new TLegend(0.7, 0.7, 0.9, 0.9);
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);

    std::vector<TF1*> fitFuncs_BDT;
    std::vector<TF1*> fitFuncs_DRN;
    std::vector<std::vector<double>> fitParams_BDT;
    std::vector<std::vector<double>> fitParams_DRN;
    double chi2_ndf_BDT;
    double chi2_ndf_DRN;

    for (size_t i = 0; i < fileNames.size(); i++) {
        TFile *file = TFile::Open(fileNames[i].c_str(), "READ");
        if (!file || file->IsZombie()) {
            std::cerr << "Error: Cannot open " << fileNames[i] << "!" << std::endl;
            continue;
        }

        TH1F *BDT_ratio = (TH1F*) file->Get("BDT_corrbyGen");
        TH1F *DRN_ratio = (TH1F*) file->Get("DRN_corrbyGen");

        if (!BDT_ratio || !DRN_ratio) {
            std::cerr << "Error: Missing histograms in " << fileNames[i] << "!" << std::endl;
            file->Close();
            continue;
        }

        BDT_ratio->Scale(1.0 / BDT_ratio->Integral());
        DRN_ratio->Scale(1.0 / DRN_ratio->Integral());
        
        BDT_ratio->Rebin(10);
        DRN_ratio->Rebin(10);


        
        // DRN_ratio->GetYaxis()->SetRangeUser(0, 1.5 *MaxY);
        // BDT_ratio->GetXaxis()->SetRangeUser(xmin, xmax);
        // DRN_ratio->GetXaxis()->SetRangeUser(xmin, xmax);


        // //////////// Iterative fitting///////////////////////////////
        // int peakBin_BDT = BDT_ratio->GetMaximumBin();
        // double peakValue_BDT = BDT_ratio->GetXaxis()->GetBinCenter(peakBin_BDT);


        // int peakBin_DRN = DRN_ratio->GetMaximumBin();
        // double peakValue_DRN = DRN_ratio->GetXaxis()->GetBinCenter(peakBin_DRN);


        int nBins = BDT_ratio->GetNbinsX();
        TGraphErrors *graph_BDT = new TGraphErrors(nBins);
        TGraphErrors *graph_DRN = new TGraphErrors(nBins);

        for (int j = 1; j <= nBins; j++) {
            double x = BDT_ratio->GetBinCenter(j);
            double y_BDT = BDT_ratio->GetBinContent(j);
            double y_DRN = DRN_ratio->GetBinContent(j);
            double bin_err_BDT = BDT_ratio->GetBinError(j);
            double bin_err_DRN = DRN_ratio->GetBinError(j);
            // double err_BDT = (y_BDT > 0) ? 1.0 / sqrt(Entries[i] * y_BDT) : 0.0;
            // double err_DRN = (y_DRN > 0) ? 1.0 / sqrt(Entries[i] * y_DRN) : 0.0;
            double err_BDT = (y_BDT > 0) ?  bin_err_BDT/y_BDT: 0.0;
            double err_DRN = (y_DRN > 0) ?  bin_err_DRN/y_DRN: 0.0;

            double width  = BDT_ratio->GetBinWidth(j);       // bin width
            double xErr   = width / 2.0;                    // symmetric error

            double weight = exp( -pow((x - 1.0) / 1.0, 2) / 2.0 );

            graph_BDT->SetPoint(j - 1, x, y_BDT);
            graph_BDT->SetPointError(j - 1, xErr, bin_err_BDT);
            graph_BDT->SetMarkerSize(0.2);
    
            graph_DRN->SetPoint(j - 1, x, y_DRN);
            graph_DRN->SetPointError(j - 1, xErr, bin_err_DRN);
            graph_DRN->SetMarkerSize(0.2);

            // BDT_ratio->SetBinError(j, err_BDT);
            // DRN_ratio->SetBinError(j, err_DRN);

            // BDT_ratio->SetBinError(j, err_BDT/weight);
            // DRN_ratio->SetBinError(j, err_DRN/weight);
        }

        graph_BDT->SetMarkerStyle(20 + i);
        graph_BDT->SetMarkerColor(colors[i]);
        graph_BDT->SetLineColor(colors[i]);

        graph_DRN->SetMarkerStyle(21 + i);
        graph_DRN->SetMarkerColor(colors[i+1]);
        graph_DRN->SetLineColor(colors[i+1]);

        mg->Add(graph_BDT);
        mg->Add(graph_DRN);
        // legend_BDT->AddEntry(graph_BDT, labels[i].c_str(), "lp");
        // legend_DRN->AddEntry(graph_DRN, labels[i+1].c_str(), "lp");
        legend->AddEntry(graph_BDT, "BDT");
        legend->AddEntry(graph_DRN, "DRN");
        legend->SetTextSize(0.03); // or a smaller value


                        // Initialize fit function
        TF1 *fitFunc_BDT = new TF1(Form("fit_BDT_%zu", i), Cruijff, xmin, xmax, 6);
        TF1 *fitFunc_DRN = new TF1(Form("fit_DRN_%zu", i), Cruijff, xmin , xmax, 6);
        fitFunc_BDT->SetParLimits(1, 0.00001, 0.1);
        fitFunc_BDT->SetParLimits(2, 0.00001, 0.1);
        fitFunc_DRN->SetParLimits(1, 0.00001, 0.1);
        fitFunc_DRN->SetParLimits(2, 0.00001, 0.1);
        
        fitFunc_BDT->SetParameters(1.0, 0.01, 0.01, 0.1, 0.1, BDT_ratio->GetMaximum());
        fitFunc_DRN->SetParameters(1.0, 0.01, 0.01, 0.1, 0.1, DRN_ratio->GetMaximum());

        double old_mean_BDT = BDT_ratio->GetMean();
        double old_sigma_BDT = 0;
        
        for (int iter = 0; iter < 400; ++iter) {
            BDT_ratio->Fit(fitFunc_BDT, "RE");
        
            double mean_BDT = fitFunc_BDT->GetParameter(0);
            double sigmaL_BDT = fitFunc_BDT->GetParameter(1);
            double sigmaR_BDT = fitFunc_BDT->GetParameter(2);
            double sigma_BDT = (sigmaL_BDT+sigmaR_BDT)/2.0;
        
            bool mean_converged   = std::abs(old_mean_BDT - mean_BDT) / mean_BDT < 0.00002;
            bool sigma_converged = std::abs(old_sigma_BDT - sigma_BDT) / sigma_BDT < 0.00002;
        
            if (mean_converged && sigma_converged)
                break;
        
            old_mean_BDT = mean_BDT;
            old_sigma_BDT = sigma_BDT;
        }

        double old_mean_DRN = DRN_ratio->GetMean();
        double old_sigma_DRN = 0;
        
        for (int iter = 0; iter < 400; ++iter) {
            DRN_ratio->Fit(fitFunc_DRN, "RE");
        
            double mean_DRN = fitFunc_DRN->GetParameter(0);
            double sigmaL_DRN = fitFunc_DRN->GetParameter(1);
            double sigmaR_DRN = fitFunc_DRN->GetParameter(2);
            double sigma_DRN = (sigmaL_DRN+sigmaR_DRN)/2.0;
        
            bool mean_converged   = std::abs(old_mean_DRN - mean_DRN) / mean_DRN < 0.00002;
            bool sigma_converged = std::abs(old_sigma_DRN - sigma_DRN) / sigma_DRN < 0.00002;
        
            if (mean_converged && sigma_converged)
                break;
        
            old_mean_DRN = mean_DRN;
            old_sigma_DRN = sigma_DRN;
        }

        BDT_ratio->GetXaxis()->SetRangeUser(xmin, xmax);
        DRN_ratio->GetXaxis()->SetRangeUser(xmin, xmax);


        
        // Iterative fitting process
        // double previousParams[6] = {1.0, 0.01, 0.01, 0.1, 0.1, BDT_ratio->GetMaximum()};
        // double accuracy = 1e-6;
        // int maxIterations = 100, iteration = 0;
        // double paramDifference = accuracy + 1.0;

        // while (paramDifference > accuracy && iteration < maxIterations) {
        //     iteration++;
        //     BDT_ratio->Fit(fitFunc_BDT, "RQ");
        //     DRN_ratio->Fit(fitFunc_DRN, "RQ");

        //     paramDifference = 0;
        //     for (int j = 0; j < 6; j++) {
        //         double diff_BDT = fabs(fitFunc_BDT->GetParameter(j) - previousParams[j]);
        //         double diff_DRN = fabs(fitFunc_DRN->GetParameter(j) - previousParams[j]);
        //         paramDifference += diff_BDT + diff_DRN;
        //         previousParams[j] = fitFunc_BDT->GetParameter(j);
        //     }
        // }

        // BDT_ratio->Fit(fitFunc_BDT, "RQ");
        // DRN_ratio->Fit(fitFunc_DRN, "RQ");
        std::vector<double> params_BDT;
        std::vector<double> params_DRN;

        for (int j = 0; j < 6; j++) {
            params_BDT.push_back(fitFunc_BDT->GetParameter(j));
            params_DRN.push_back(fitFunc_DRN->GetParameter(j));
        }

        double chi2_BDT = fitFunc_BDT->GetChisquare();
        int ndf_BDT = fitFunc_BDT->GetNDF();
        chi2_ndf_BDT = (ndf_BDT != 0) ? chi2_BDT / ndf_BDT : 0;  // Avoid division by zero

        double chi2_DRN = fitFunc_DRN->GetChisquare();
        int ndf_DRN = fitFunc_DRN->GetNDF();
        chi2_ndf_DRN = (ndf_DRN != 0) ? chi2_DRN / ndf_DRN : 0;  // Avoid division by zero


        fitParams_BDT.push_back(params_BDT);
        fitParams_DRN.push_back(params_DRN);
        
        fitFunc_BDT->SetLineColor(colors[i]);
        fitFunc_BDT->SetLineStyle(i);
        fitFunc_DRN->SetLineColor(colors[i+1]);
        fitFunc_DRN->SetLineStyle(i);
        fitFuncs_BDT.push_back(fitFunc_BDT);
        fitFuncs_DRN.push_back(fitFunc_DRN);

        file->Close();
    }

    TGraphErrors* graph_BDT = (TGraphErrors*) mg->GetListOfGraphs()->At(0);
    TGraphErrors* graph_DRN = (TGraphErrors*) mg->GetListOfGraphs()->At(1);

    TF1* fitFunc_BDT = fitFuncs_BDT[0];
    TF1* fitFunc_DRN = fitFuncs_DRN[0];

    SaveGraphWithFitData(graph_BDT, fitFunc_BDT, xmin, xmax, "my_fit_results_BDT.txt");
    SaveGraphWithFitData(graph_DRN, fitFunc_DRN, xmin, xmax, "my_fit_results_DRN.txt");

    TCanvas *c_twoGraph = new TCanvas("c_twoGraph", "BDT vs DRN", 800, 600);
    c_twoGraph->cd();

        // Draw the BDT multigraph first
    double yMax = 0;
    TIter next(mg->GetListOfGraphs());
    TGraph* g;
    while ((g = (TGraph*)next())) {
        for (int i = 0; i < g->GetN(); ++i) {
            double x, y;
            g->GetPoint(i, x, y);
            if (y > yMax) yMax = y;
        }
    }
    mg->GetYaxis()->SetRangeUser(0, yMax * 1.25); // Add some margin (25%) on top

    
    // Draw the BDT multigraph first
    mg->Draw("AP");
    mg->GetXaxis()->SetLimits(xmin, xmax);
    mg->GetXaxis()->SetTitle("E_{corr}/E_{gen}");
    mg->GetYaxis()->SetTitle("Normalized");
    // mg->GetYaxis()->SetRangeUser(0, mg->GetYaxis()->GetXmax() * 0.25);
    mg->GetXaxis()->SetLimits(xmin, xmax);
    // mg_DRN->Draw("AP same");
    // mg_DRN->GetXaxis()->SetLimits(xmin, xmax);
    // mg_DRN->GetYaxis()->SetRangeUser(0, mg_DRN->GetYaxis()->GetXmax() * 0.14);

    // Draw the DRN multigraph on the same canvas
    
    // Add labels and legends
    TLatex *lateX = new TLatex();
    lateX->DrawTextNDC(0.10, 0.95, "BDT vs DRN response comparison_Photon_Ideal_IC");
    // legend_BDT->Draw();
    // legend_DRN->Draw("same");
    // TLegend *legend = new TLegend(0.65, 0.75, 0.85, 0.9);
    // legend->AddEntry(mg, "BDT", "P");
    // legend->AddEntry(mg_DRN, "DRN", "P");
    // legend->Draw();

    lateX->SetTextSize(0.025);
    // lateX->DrawTextNDC(0.17, 0.55, chi2_ndf_DRN);
    lateX->SetTextColor(colors[1]);
    lateX->DrawLatexNDC(0.17, 0.65, Form("#chi^{2}/NDF_DRN = %.4f", chi2_ndf_DRN));

    lateX->SetTextSize(0.025);
    // lateX->DrawTextNDC(0.17, 0.55, chi2_ndf_DRN);
    lateX->SetTextColor(colors[0]);
    lateX->DrawLatexNDC(0.17, 0.70, Form("#chi^{2}/NDF_BDT = %.4f", chi2_ndf_BDT));


    
    // Draw BDT fit functions
    int k = 0;
    for (auto fit : fitFuncs_BDT){
        fit->Draw("same");
        
        TLatex latex;
        latex.SetTextSize(0.025);
        latex.SetTextColor(colors[k]);
        
        double xTextPos = xmin + (xmax - xmin) * 0.05;  // 5% from left
        double yBase    = yMax * 1.2;                  // 95% of Y max (near top)
        double yStep    = yMax * 0.07;                  // spacing between lines
        
        // Display fit parameters next to the fit line
        latex.DrawLatex(xTextPos, yBase - 0 * yStep, Form("#mu = %.4f, #sigmaL = %.4f, #sigmaR = %.4f",
                                              fit->GetParameter(0),
                                              fit->GetParameter(1),
                                              fit->GetParameter(2)));
        k++;
    }
    
    // mg_DRN->Draw("P SAME");
    // mg_DRN->GetXaxis()->SetLimits(xmin, xmax);
    // mg_DRN->GetYaxis()->SetRangeUser(0, mg_DRN->GetYaxis()->GetXmax() * 0.14);
    // legend_DRN->Draw();
    // Draw DRN fit functions
    int l = 0; 
    for (auto fit : fitFuncs_DRN){ 
        fit->Draw("same");
        
        TLatex latex;
        latex.SetTextSize(0.025);
        latex.SetTextColor(colors[l+1]);
        
        double xTextPos = xmin + (xmax - xmin) * 0.05;  // 5% from left
        double yBase    = yMax * 1.2;                  // 95% of Y max (near top)
        double yStep    = yMax * 0.07;                  // spacing between lines
        
        // Display fit parameters next to the fit line
        latex.DrawLatex(xTextPos, yBase - 1 * yStep, Form("#mu = %.4f, #sigmaL = %.4f, #sigmaR = %.4f",
                                              fit->GetParameter(0),
                                              fit->GetParameter(1),
                                              fit->GetParameter(2)));
        l++;
    }
    legend->Draw();
    
    // Save the canvas
    c_twoGraph->SaveAs("BDT_vs_DRN_comparison_DY_photon.pdf");
    c_twoGraph->SaveAs("BDT_vs_DRN_comparison_DY_photon.png");

}
