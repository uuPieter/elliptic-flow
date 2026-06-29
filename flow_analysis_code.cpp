
#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include "TFile.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TLegend.h"
#include <numeric>
#include <iomanip>  
#include "TLine.h"
#include "TF1.h"


int main() {
    // Open the ROOT file //

    TFile *f = TFile::Open("/home/pieter/Downloads/Lood_result_train_691371.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return 1;
    }

    // create new Root file to save figures and th1d //
    TFile* file = TFile::Open("pieter_train_691371.root", "RECREATE");
    file->cd();

    // save between loop
    std::vector<int> histTypeIndex; // 0 = S+B, 1 = B, 2 = S
    std::vector<std::string> histLabels;
    std::vector<TH1D*> hFlowVec;
    std::vector<std::string> cenLabels = {"base0", "base_10","base_20","base_30","base_50","occuML_0", "occuML_10", "occuML_20", "occuML_30", "occuML_50", "occuMS_0", "occuMS_10", "occuMS_20", "occuMS_30", "occuMS_50"};
    std::vector<int> colors = {kRed+1, kBlue, kBlack, kMagenta+1, kCyan+1, kGreen+2, kOrange+1, kViolet+1, kAzure+1, kPink+1, kYellow+1, kGray+1, kRed-3, kBlue-3, kBlack-3};

    // start the loop over the different centralities //
    for (int ip = 0; ip < 1; ++ip) {
    
    // enter the correct subfolder of root file //
    TDirectory* dir = file->mkdir(Form("centrality%s", cenLabels[ip].c_str()));
    dir->cd(); 

        TString dirName = Form("dielectron_%s", cenLabels[ip].c_str());
    
        // read the data form origional root file //
    TDirectory *d_dielectron = (TDirectory*) f->Get(dirName);
    
    if (!d_dielectron) { std::cerr << "dielectron not found!\n"; return 1; }

    TDirectory *d_Pair = (TDirectory*) d_dielectron->Get("Pair");
    if (!d_Pair) { std::cerr << "Pair not found!\n"; return 1; }

    TDirectory *d_same = (TDirectory*) d_Pair->Get("same");
    if (!d_same) { std::cerr << "same not found!\n"; return 1; }

    //  get the uls
    TDirectory *d_uls = (TDirectory*) d_same->Get("uls");
    THnSparseD* hs_uls = (THnSparseD*) d_uls->Get("hs");
    // get the plus plus 
    TDirectory *d_lspp = (TDirectory*) d_same->Get("lspp");
    THnSparseD* hs_lspp = (THnSparseD*) d_lspp->Get("hs");
    // get the min min
    TDirectory *d_lsmm = (TDirectory*) d_same->Get("lsmm");
    THnSparseD* hs_lsmm = (THnSparseD*) d_lsmm->Get("hs");
    // get the mix event data
    TDirectory *d_mix = (TDirectory*) d_Pair->Get("mix");
    THnSparseD* hs_uls_mix = (THnSparseD*) ((TDirectory*)d_mix->Get("uls"))->Get("hs");
    THnSparseD* hs_lspp_mix = (THnSparseD*) ((TDirectory*)d_mix->Get("lspp"))->Get("hs");
    THnSparseD* hs_lsmm_mix = (THnSparseD*) ((TDirectory*)d_mix->Get("lsmm"))->Get("hs");

     // get the n events
     if (!d_dielectron) { std::cerr << "dielectron not found!\n"; return 1; }
     TDirectory *d_Event = (TDirectory*) d_dielectron->Get("Event");
     if (!d_Event) { std::cerr << "Pair not found!\n"; return 1; }
     TDirectory *d_after = (TDirectory*) d_Event->Get("after");
     if (!d_after) { std::cerr << "Pair not found!\n"; return 1; }
     TH2D* hZvtx = (TH2D*) d_after->Get("hZvtx"); 
     if (!hZvtx) { std::cerr << "hZvtx histogram not found!\n"; return 1; }
     double multiplicity = hZvtx->GetEntries();
     std::cout << "nEvents" << multiplicity << std::endl;



    // project on to 2D
    TH2D* h_proj_lsmm_2D = (TH2D*) hs_lsmm->Projection(1, 0);
    TH2D* h_proj_lspp_2D = (TH2D*) hs_lspp->Projection(1, 0);
    TH2D* h_proj_uls_2D = (TH2D*) hs_uls->Projection(1, 0);


    TH2D* h_proj_uls_mix = (TH2D*) hs_uls_mix->Projection(1, 0);
    TH2D* h_proj_lspp_mix = (TH2D*) hs_lspp_mix->Projection(1, 0);
    TH2D* h_proj_lsmm_mix = (TH2D*) hs_lsmm_mix->Projection(1, 0);

    // ======================= 3D projections =======================//
    //                projection of Mee Pt and DCA                   //
    // ======================= 3D projections =======================//
    TH3D* h_proj_uls_mix_3D = (TH3D*) hs_uls_mix->Projection(0, 1, 2);
    TH3D* h_proj_lspp_mix_3D = (TH3D*) hs_lspp_mix->Projection(0, 1, 2);
    TH3D* h_proj_lsmm_mix_3D = (TH3D*) hs_lsmm_mix->Projection(0, 1, 2);
    TH3D* h_proj_uls_3D = (TH3D*) hs_uls->Projection(0 , 1 , 2);  
    TH3D* h_proj_lspp_3D = (TH3D*) hs_lspp->Projection(0 , 1 , 2);  
    TH3D* h_proj_lsmm_3D = (TH3D*) hs_lsmm->Projection(0 , 1 , 2);  

    // ======================= 3D projections =======================//
    //                projection of Mee Pt and Qvec                   //
    // ======================= 3D projections =======================//

    double delta = 1e-4;
    


    // Make a clone to not touch original
THnSparseD* hs_uls_slice  = (THnSparseD*) hs_uls->Clone("hs_uls_slice");
THnSparseD* hs_lspp_slice = (THnSparseD*) hs_lspp->Clone("hs_lspp_slice");
THnSparseD* hs_lsmm_slice = (THnSparseD*) hs_lsmm->Clone("hs_lsmm_slice");


TH3D* h_proj_uls3D = (TH3D*) hs_uls_slice->Projection(0 , 1 , 4);  
TH3D* h_proj_lspp3D = (TH3D*) hs_lspp_slice->Projection(0 , 1 , 4);  
TH3D* h_proj_lsmm3D = (TH3D*) hs_lsmm_slice->Projection(0 , 1 , 4); 


    TAxis* xax = h_proj_uls_3D->GetXaxis();
    TAxis* yax = h_proj_uls_3D->GetYaxis();
    TAxis* zax = h_proj_uls_3D->GetZaxis();
    int nX = h_proj_uls_3D->GetNbinsX();
    int nY = h_proj_uls_3D->GetNbinsY();
    int nZ = h_proj_uls_3D->GetNbinsZ();
    
    // ======================= rebinning =======================//
    //        rebinning to get rid of under and oferflow        //
    // ======================= rebinning =======================//
    struct BinRange {
        double start;
        double end;
        double step;
    };
    
    std::vector<BinRange> xRanges = {
        {0.0, 0.45, 0.01},   
        {0.45, 1.2, 0.01},   
        {1.2, 2.8, 0.01},    
        {2.8, 3.2, 0.01},    
        {3.2, 4.0, 0.01}      
    };
    std::vector<BinRange> yRanges = {
        {0.0, 5.0, 0.1},   
        {5.0, 10.0, 0.1}  
    };
    std::vector<BinRange> zRanges = {
        {0.0, 5.0, 0.2}, 
        {5.0, 10.0, 0.5},  
        {10.0, 40.0, 1.0}
    };
    

    // --- X-axis bins ---
    std::vector<double> xbins;
    for (auto &[start, end, step] : xRanges) {
        if (xbins.empty()) xbins.push_back(start); 
        for (double x = start + step; x < end - 1e-4; x += step) {
            xbins.push_back(x);
        }
        xbins.push_back(end); 
    }
    
    // --- Y-axis bins ---
    std::vector<double> ybins;
    for (auto &[start, end, step] : yRanges) {
        if (ybins.empty()) ybins.push_back(start);
        for (double y = start + step; y < end - 1e-4; y += step) {
            ybins.push_back(y);
        }
        ybins.push_back(end);
    }
    
    // --- Z-axis bins ---
    std::vector<double> zbins;
    for (auto &[start, end, step] : zRanges) {
        if (zbins.empty()) zbins.push_back(start);
        for (double z = start + step; z < end - 1e-4; z += step) {
            zbins.push_back(z);
        }
        zbins.push_back(end);
    }
    
    auto FindBinLow2 = [&](TAxis* ax, double edge) {
        return ax->FindBin(edge + delta);
    };
    auto FindBinHigh2 = [&](TAxis* ax, double edge) {
        return ax->FindBin(edge - delta);
    };
    
    
    std::vector<std::pair<int,int>> meeSlices2;
    for (size_t i = 0; i < xbins.size() - 1; ++i) {
    
        int lo = FindBinLow2 (xax, xbins[i]);
        int hi = FindBinHigh2(xax, xbins[i+1]);
    
        lo = std::max(lo, 1);
        hi = std::min(hi, xax->GetNbins());
    
        meeSlices2.emplace_back(lo, hi);
    }
    
    std::vector<std::pair<int,int>> ptSlices2;
    for (size_t i = 0; i < ybins.size() - 1; ++i) {
    
        int lo = FindBinLow2 (yax, ybins[i]);
        int hi = FindBinHigh2(yax, ybins[i+1]);
    
        lo = std::max(lo, 1);
        hi = std::min(hi, yax->GetNbins());
    
        ptSlices2.emplace_back(lo, hi);
    }
    
    std::vector<std::pair<int,int>> DCASlices2;
    for (size_t i = 0; i < zbins.size() -1; ++i) {
    
        int lo = FindBinLow2 (zax, zbins[i]);
        int hi = FindBinHigh2(zax, zbins[i+1]);
    
        lo = std::max(lo, 1);
        hi = std::min(hi, zax->GetNbins());
    
        DCASlices2.emplace_back(lo, hi);
    }
    
    TH3D* h_rebinned_uls_3D  = new TH3D("h_rebinned_uls",  "ULS;X;Y;Z",xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data(),zbins.size()-1, zbins.data());
    TH3D* h_rebinned_lspp_3D = new TH3D("h_rebinned_lspp", "LS++;X;Y;Z",xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data(),zbins.size()-1, zbins.data());
    TH3D* h_rebinned_lsmm_3D = new TH3D("h_rebinned_lsmm", "LS--;X;Y;Z",xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data(),zbins.size()-1, zbins.data());
    TH3D* h_rebinned_uls_mix_3D = new TH3D("h_rebinned_uls_mix", "ULS;X;Y;Z",xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data(),zbins.size()-1, zbins.data());
    TH3D* h_rebinned_lspp_mix_3D = new TH3D("h_rebinned_lspp_mix", "LS++;X;Y;Z",xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data(),zbins.size()-1, zbins.data());
    TH3D* h_rebinned_lsmm_mix_3D = new TH3D("h_rebinned_lsmm_mix", "LS--;X;Y;Z",xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data(),zbins.size()-1, zbins.data());
    
    for (size_t i = 0; i < meeSlices2.size(); ++i) {
        for (size_t j = 0; j < ptSlices2.size(); ++j) {
            for (size_t k = 0; k < DCASlices2.size(); ++k) {
    
                double sum_uls  = 0.0;
                double sum_lspp = 0.0;
                double sum_lsmm = 0.0;
                double sum_uls_mix  = 0.0;
                double sum_lspp_mix = 0.0;
                double sum_lsmm_mix = 0.0;
    
                int xlo = meeSlices2[i].first;
                int xhi = meeSlices2[i].second;
    
                int ylo = ptSlices2[j].first;
                int yhi = ptSlices2[j].second;
    
                int zlo = DCASlices2[k].first;
                int zhi = DCASlices2[k].second;
    
                for (int ix = xlo; ix <= xhi; ++ix) {
                    for (int iy = ylo; iy <= yhi; ++iy) {
                        for (int iz = zlo; iz <= zhi; ++iz) {
    
                            sum_uls  += h_proj_uls_3D ->GetBinContent(ix,iy,iz);
                            sum_lspp += h_proj_lspp_3D->GetBinContent(ix,iy,iz);
                            sum_lsmm += h_proj_lsmm_3D->GetBinContent(ix,iy,iz);
                            sum_uls_mix  += h_proj_uls_mix_3D ->GetBinContent(ix,iy,iz);
                            sum_lspp_mix += h_proj_lspp_mix_3D->GetBinContent(ix,iy,iz);
                            sum_lsmm_mix += h_proj_lsmm_mix_3D->GetBinContent(ix,iy,iz);
    
                        }
                    }
                }
                h_rebinned_uls_3D ->SetBinContent(i+1, j+1, k+1, sum_uls);
                h_rebinned_lspp_3D->SetBinContent(i+1, j+1, k+1, sum_lspp);
                h_rebinned_lsmm_3D->SetBinContent(i+1, j+1, k+1, sum_lsmm);
                h_rebinned_uls_mix_3D->SetBinContent(i+1, j+1, k+1, sum_uls_mix);
                h_rebinned_lspp_mix_3D->SetBinContent(i+1, j+1, k+1, sum_lspp_mix);
                h_rebinned_lsmm_mix_3D->SetBinContent(i+1, j+1, k+1, sum_lsmm_mix);
            }
        }
    }

// ======================== calculate variables ========================//
//         calculate uls, signal, background, S/B and R in 2D           //
// ======================== calculate variables ========================//

TH2D* h_rebinned_uls_2D  = (TH2D*)h_rebinned_uls_3D->Project3D("yx");
TH2D* h_rebinned_lspp_2D = (TH2D*)h_rebinned_lspp_3D->Project3D("yx");
TH2D* h_rebinned_lsmm_2D = (TH2D*)h_rebinned_lsmm_3D->Project3D("yx");

std::cout << "Original integral 2D: "
          << h_rebinned_uls_2D->Integral() << std::endl;

std::cout << "Rebinned integral 3D: "
          << h_rebinned_uls_3D->Integral() << std::endl;

          std::cout << "Original integral 2D: "
          << h_rebinned_lspp_2D->Integral() << std::endl;

std::cout << "Rebinned integral 3D: "
          << h_rebinned_lspp_3D->Integral() << std::endl;

          std::cout << "Original integral 2D: "
          << h_rebinned_lsmm_2D->Integral() << std::endl;

std::cout << "Rebinned integral 3D: "
          << h_rebinned_lsmm_3D->Integral() << std::endl;

TH2D* h_rebinned_uls_mix_2D  = (TH2D*)h_rebinned_uls_mix_3D->Project3D("yx");
TH2D* h_rebinned_lspp_mix_2D = (TH2D*)h_rebinned_lspp_mix_3D->Project3D("yx");
TH2D* h_rebinned_lsmm_mix_2D = (TH2D*)h_rebinned_lsmm_mix_3D->Project3D("yx");

int nBinsXrebin = h_rebinned_uls_2D->GetNbinsX();
int nBinsYrebin = h_rebinned_uls_2D->GetNbinsY();

TH2D *hZ_2D_rebinned = new TH2D("hZ_2D_rebinned", "Z;X;Y",xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data());
TH2D *hB_2D_rebinned = new TH2D("hB_2D_rebinned", "B;X;Y",xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data());
TH2D *hS_over_B_2D_rebinned = new TH2D("hSoverB_2D_rebinned", "S/B;X;Y",xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data());
TH2D *hR_2D_rebinned = new TH2D("hR_2D_rebinned", "R;X;Y", xbins.size()-1, xbins.data(), ybins.size()-1, ybins.data());
TH2D *hN_2D_rebinned = new TH2D("hN_2D_rebinned", "N;X;Y",xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data());
h_rebinned_uls_mix_2D->Write();
h_rebinned_lspp_mix_2D->Write();
h_rebinned_lsmm_mix_2D->Write();
// get the entries from each individual bin and calculate the variables we are interessted in 
for (int ix = 1; ix <= nBinsXrebin; ix++) {
    for (int iy = 1; iy <= nBinsYrebin; iy++) {
        double N_pm_same = h_rebinned_uls_2D->GetBinContent(ix, iy);
        double N_pp_same = h_rebinned_lspp_2D->GetBinContent(ix, iy);
        double N_mm_same = h_rebinned_lsmm_2D->GetBinContent(ix, iy);

        double N_pm_mix = h_rebinned_uls_mix_2D->GetBinContent(ix, iy);
        double N_pp_mix = h_rebinned_lspp_mix_2D->GetBinContent(ix, iy);
        double N_mm_mix = h_rebinned_lsmm_mix_2D->GetBinContent(ix, iy);
        double N = N_pm_same;
        hN_2D_rebinned     ->SetBinContent(ix, iy, N );

        if (N_pp_mix > 0 && N_mm_mix > 0 && N_pm_mix > 0) {
            double R = N_pm_mix / (2.0 * std::sqrt(N_pp_mix * N_mm_mix));
            double B = 2.0  * std::sqrt(R*N_pp_same * R*N_mm_same);
            double Z = N_pm_same - B;
            double S_over_B = (B > 0) ? Z / B : 0.0;

            hR_2D_rebinned     ->SetBinContent(ix, iy, R );
            hZ_2D_rebinned     ->SetBinContent(ix, iy, Z );
            hB_2D_rebinned     ->SetBinContent(ix, iy, B );
            hS_over_B_2D_rebinned->SetBinContent(ix, iy, S_over_B);
            double A = N_pm_mix;
            double D = 2.0 * std::sqrt(N_pp_mix * N_mm_mix);
            // --- uncertainties ---
            double sigmaA = std::sqrt(N_pm_mix);
            
            // propagate denominator 
            double sigmaD = D * 0.5 * std::sqrt(
                1.0 / N_pp_mix + 1.0 / N_mm_mix
            );
            
            // now propagate ratio
            double sigmaR = std::abs(R) * std::sqrt(
                (sigmaA / A)*(sigmaA / A) +
                (sigmaD / D)*(sigmaD / D)
            );
            
            hR_2D_rebinned->SetBinError(ix, iy, sigmaR);
    
            // Background uncertainty
            // --- propagate C first ---
            double C = std::sqrt(N_pp_same * N_mm_same);
            double sigmaC = C * 0.5 * std::sqrt(
            1.0 / N_pp_same + 1.0 / N_mm_same
            );
    
            // --- now combine with R uncertainty ---
            double sigmaB = 0.0;
    
            if (R > 0 && C > 0 && sigmaR > 0 && sigmaC > 0) {
                sigmaB = B * std::sqrt(
                    (sigmaR / R)*(sigmaR / R) +
                    (sigmaC / C)*(sigmaC / C)
                );
            }
            if (sigmaB < 0 || std::isnan(sigmaB) || std::isinf(sigmaB)) sigmaB = 0;
            hB_2D_rebinned->SetBinError(ix, iy, sigmaB);
    
            // Z uncertainty
            double varZ = N_pm_same + sigmaB*sigmaB;
            if (varZ < 0 || std::isnan(varZ) || std::isinf(varZ)) varZ = 0;
            hZ_2D_rebinned->SetBinError(ix, iy, std::sqrt(varZ));
    
            // N uncertainty
            hN_2D_rebinned->SetBinError(ix, iy,  std::sqrt(N_pm_same));
                            // Signal over background uncertainty
            double sigmaSB = 0.0;

            if (Z > 0 && B > 0) {
                double sigmaZ = std::sqrt(varZ);

            
                double termZ = sigmaZ / Z;
                double termB = sigmaB / B;
                sigmaSB = (Z/B) * std::sqrt(termZ*termZ + termB*termB);
            if (std::isnan(sigmaSB) || std::isinf(sigmaSB)) sigmaSB = 0.0;
            }
            hS_over_B_2D_rebinned->SetBinError(ix, iy, sigmaSB);
        }
    }
} 
std::cout << "Original integral: "
          << hN_2D_rebinned->Integral() << std::endl;

std::cout << "Rebinned integral: "
          << h_rebinned_uls_3D->Integral() << std::endl;

hR_2D_rebinned->Write();
TCanvas *c2 = new TCanvas("c","2D plot",800,700);

c2->SetTicks(1,1);
c2->SetRightMargin(0.15);  // space for color bar
c2->SetLeftMargin(0.12);
c2->SetBottomMargin(0.12);

hR_2D_rebinned->SetStats(0);

// set Z axis range
hR_2D_rebinned->GetZaxis()->SetRangeUser(0.96, 1.04);
hR_2D_rebinned->GetYaxis()->SetTitle("p_{T} (GeV/c)");
hR_2D_rebinned->GetXaxis()->SetTitle("m_{ee} (GeV/c^{2})");

// axis styling (optional but nicer)
hR_2D_rebinned->GetXaxis()->SetTitleSize(0.045);
hR_2D_rebinned->GetYaxis()->SetTitleSize(0.045);
hR_2D_rebinned->GetXaxis()->SetLabelSize(0.04);
hR_2D_rebinned->GetYaxis()->SetLabelSize(0.04);

// draw
hR_2D_rebinned->Draw("COLZ");
TLatex latex;




c2->SaveAs("hR_2D_rebinned_pbpb.pdf");


double xMin_rebin = h_rebinned_uls_3D->GetXaxis()->GetXmin();
double xMax_rebin = h_rebinned_uls_3D->GetXaxis()->GetXmax();
double yMin_rebin = h_rebinned_uls_3D->GetYaxis()->GetXmin();
double yMax_rebin = h_rebinned_uls_3D->GetYaxis()->GetXmax();

TH1D* hR_projX = hR_2D_rebinned->ProjectionX("hR_projX");

TH1D* hR_avgX = (TH1D*)hR_projX->Clone("hR_avgX");
hR_avgX->Reset();

for (int ix = 1; ix <= nBinsXrebin; ix++) {
    double sumR = 0;
    double sumW = 0;

    for (int iy = 1; iy <= nBinsYrebin; iy++) {
        double R = hR_2D_rebinned->GetBinContent(ix, iy);
        double w = hN_2D_rebinned->GetBinContent(ix, iy);

        sumR += w * R;
        sumW += w;
    }

    double avg = (sumW > 0) ? sumR / sumW : 0;
    hR_avgX->SetBinContent(ix, avg);
}
hR_avgX->Write("hR_avgX");
TH1D* hZ_projX = hZ_2D_rebinned->ProjectionX("hZ_projX");
TH1D* hB_projX = hB_2D_rebinned->ProjectionX("hB_projX");
TH1D* hN_projX = hN_2D_rebinned->ProjectionX("hN_projX");
TH1D* hSoverB_projX = hS_over_B_2D_rebinned->ProjectionX("hSoverB_projX");
for (int ix = 1; ix <= nBinsXrebin; ix++) {
    double contentN = hN_projX->GetBinContent(ix);
    double contentZ = hZ_projX->GetBinContent(ix);
    double contentB = hB_projX->GetBinContent(ix);
    double contentSB = hSoverB_projX->GetBinContent(ix);
}
TH1D* hZ_projY = hZ_2D_rebinned->ProjectionY("hZ_projY");
TH1D* hB_projY = hB_2D_rebinned->ProjectionY("hB_projY");
TH1D* hR_projY = hR_2D_rebinned->ProjectionY("hR_projY");
TH1D* hN_projY = hN_2D_rebinned->ProjectionY("hN_projY");
TH1D* hSoverB_projY = hS_over_B_2D_rebinned->ProjectionY("hSoverB_projY");
double norm = (multiplicity > 0) ? 1.0 / multiplicity : 0.0;

hN_projX->Scale(norm);
hZ_projX->Scale(norm);
hB_projX->Scale(norm);

hN_projY->Scale(norm);
hZ_projY->Scale(norm);
hB_projY->Scale(norm);

hN_projX->Write("hN_mass_using2D");
hZ_projX->Write("hZ_mass_using2D");
hB_projX->Write("hB_mass_using2D");
hR_projX->Write("hR_mass_using2D");

hN_projY->Write("hN_pT_using2D");
hZ_projY->Write("hZ_pT_using2D");
hB_projY->Write("hB_pT_using2D");
hR_projY->Write("hR_pT_using2D");

 // ======================== calculate variables 1D Mee ========================//
 //       calculate uls, signal, background, S/B and R in 1D                    //
 // ======================== calculate variables 1D Mee ========================//

TH1D* h_rebinned_uls_1D_x  = h_rebinned_uls_3D->ProjectionX();  
TH1D* h_rebinned_lspp_1D_x = h_rebinned_lspp_3D->ProjectionX();
TH1D* h_rebinned_lsmm_1D_x = h_rebinned_lsmm_3D->ProjectionX();
TH1D* h_rebinned_uls_mix_1D_x  = h_rebinned_uls_mix_3D->ProjectionX();
TH1D* h_rebinned_lspp_mix_1D_x = h_rebinned_lspp_mix_3D->ProjectionX();
TH1D* h_rebinned_lsmm_mix_1D_x = h_rebinned_lsmm_mix_3D->ProjectionX();

h_rebinned_lspp_1D_x->Write("h_rebinned_lspp_1D_x");
h_rebinned_lsmm_1D_x->Write("h_rebinned_lsmm_1D_x");

TH1D *hZ_1D_mee_rebin = new TH1D("hZ_1D_mee_rebin", "Z;m_{ee};Z",xbins.size()-1, xbins.data());
TH1D *hB_1D_mee_rebin = new TH1D("hB_1D_mee_rebin", "B;m_{ee};B",xbins.size()-1, xbins.data());
TH1D *hSoverB_1D_mee_rebin = new TH1D("hSoverB_1D_mee_rebin", "S/B;m_{ee};S/B",xbins.size()-1, xbins.data());
TH1D *hR_1D_mee_rebin = new TH1D("hR_1D_mee_rebin", "R;m_{ee};R",xbins.size()-1, xbins.data());
TH1D *hN_1D_mee_rebin = new TH1D("hN_1D_mee_rebin", "N;m_{ee};N",xbins.size()-1, xbins.data());

    // Loop over bins
    for (int ix = 1; ix <= nBinsXrebin; ++ix) {
        double N_pm_same = h_rebinned_uls_1D_x->GetBinContent(ix);
        double N_pp_same = h_rebinned_lspp_1D_x->GetBinContent(ix);
        double N_mm_same = h_rebinned_lsmm_1D_x->GetBinContent(ix);
    
        double N_pm_mix = h_rebinned_uls_mix_1D_x->GetBinContent(ix);
        double N_pp_mix = h_rebinned_lspp_mix_1D_x->GetBinContent(ix);
        double N_mm_mix = h_rebinned_lsmm_mix_1D_x->GetBinContent(ix);
        double N = N_pm_same;
        hN_1D_mee_rebin     ->SetBinContent(ix, N );

            if (N_pp_mix > 0 && N_mm_mix > 0) {
                double R = N_pm_mix / (2.0 * std::sqrt(N_pp_mix * N_mm_mix));
                double B = 2.0 * R * std::sqrt(N_pp_same * N_mm_same);
                double Z = N_pm_same - B;
                double S_over_B = (B > 0) ? Z / B : 0.0;
                hR_1D_mee_rebin     ->SetBinContent(ix, R);
                hZ_1D_mee_rebin     ->SetBinContent(ix, Z);
                hB_1D_mee_rebin     ->SetBinContent(ix, B);
                hSoverB_1D_mee_rebin->SetBinContent(ix, S_over_B);

                double A = N_pm_mix;
                double D = 2.0 * std::sqrt(N_pp_mix * N_mm_mix);
                // --- uncertainties ---
                double sigmaA = std::sqrt(N_pm_mix);
                
                // propagate denominator 
                double sigmaD = D * 0.5 * std::sqrt(
                    1.0 / N_pp_mix + 1.0 / N_mm_mix
                );
                
                // now propagate ratio
                double sigmaR = std::abs(R) * std::sqrt(
                    (sigmaA / A)*(sigmaA / A) +
                    (sigmaD / D)*(sigmaD / D)
                );
                
                hR_1D_mee_rebin->SetBinError(ix, sigmaR);
        
                // Background uncertainty
                // --- propagate C first ---
                double C = std::sqrt(N_pp_same * N_mm_same);
                double sigmaC = C * 0.5 * std::sqrt(
                1.0 / N_pp_same + 1.0 / N_mm_same
                );
        
                // --- now combine with R uncertainty ---
                double sigmaB = 0.0;
        
                if (R > 0 && C > 0 && sigmaR > 0 && sigmaC > 0) {
                    sigmaB = B * std::sqrt(
                        (sigmaR / R)*(sigmaR / R) +
                        (sigmaC / C)*(sigmaC / C)
                    );
                }
                if (sigmaB < 0 || std::isnan(sigmaB) || std::isinf(sigmaB)) sigmaB = 0;
                hB_1D_mee_rebin->SetBinError(ix, sigmaB);
        
                // Z uncertainty
                double varZ = N_pm_same + sigmaB*sigmaB;
                if (varZ < 0 || std::isnan(varZ) || std::isinf(varZ)) varZ = 0;
                hZ_1D_mee_rebin->SetBinError(ix, std::sqrt(varZ));
        
                // N uncertainty
                hN_1D_mee_rebin->SetBinError(ix,  std::sqrt(N_pm_same));
                
                // Signal over background uncertainty
                double sigmaSB = 0.0;

                if (Z > 0 && B > 0) {
                    double sigmaZ = std::sqrt(varZ);
            
                    double termZ = sigmaZ / Z;
                    double termB = sigmaB / B;
                    sigmaSB = (Z/B) * std::sqrt(termZ*termZ + termB*termB);
                if (std::isnan(sigmaSB) || std::isinf(sigmaSB)) sigmaSB = 0.0;
                }
                hSoverB_1D_mee_rebin->SetBinError(ix, sigmaSB);
            
        }
    }
hZ_1D_mee_rebin->Scale(1.0 / multiplicity);
hB_1D_mee_rebin->Scale(1.0 / multiplicity);
hN_1D_mee_rebin->Scale(1.0 / multiplicity);


hR_1D_mee_rebin->Write();
hZ_1D_mee_rebin->Write();
hN_1D_mee_rebin->Write();
hB_1D_mee_rebin->Write();


 // ======================== calculate variables 1D pT ========================//
 //       calculate uls, signal, background, S/B and R in 1D                   //
 // ======================== calculate variables 1D pT ========================//

TH1D* h_rebinned_uls_1D_pt  = h_rebinned_uls_3D->ProjectionY();
TH1D* h_rebinned_lspp_1D_pt = h_rebinned_lspp_3D->ProjectionY();
TH1D* h_rebinned_lsmm_1D_pt = h_rebinned_lsmm_3D->ProjectionY();

TH1D* h_rebinned_uls_mix_1D_pt  = h_rebinned_uls_mix_3D->ProjectionY();
TH1D* h_rebinned_lspp_mix_1D_pt = h_rebinned_lspp_mix_3D->ProjectionY();
TH1D* h_rebinned_lsmm_mix_1D_pt = h_rebinned_lsmm_mix_3D->ProjectionY();
h_rebinned_lspp_1D_pt->Write();
h_rebinned_lsmm_1D_pt->Write();

TH1D *hZ_1D_pt_rebin_initial = new TH1D("hZ_1D_pt_rebin", "Z;p_{T};Z",nBinsYrebin, ybins.data());
TH1D *hB_1D_pt_rebin_initial = new TH1D("hB_1D_pt_rebin", "B;p_{T};B",nBinsYrebin, ybins.data());
TH1D *hSoverB_1D_pt_rebin_initial = new TH1D("hSoverB_1D_pt_rebin", "S/B;p_{T};S/B", nBinsYrebin, ybins.data());
TH1D *hR_1D_pt_rebin_initial = new TH1D( "hR_1D_pt_rebin", "R;p_{T};R", nBinsYrebin, ybins.data());
TH1D *hN_1D_pt_rebin_initial = new TH1D( "hN_1D_pt_rebin", "N;p_{T};N", nBinsYrebin, ybins.data());



for (int iy = 1; iy <= nBinsYrebin; ++iy) {

    double N_pm_same = h_rebinned_uls_1D_pt->GetBinContent(iy);
    double N_pp_same = h_rebinned_lspp_1D_pt->GetBinContent(iy);
    double N_mm_same = h_rebinned_lsmm_1D_pt->GetBinContent(iy);

    double N_pm_mix = h_rebinned_uls_mix_1D_pt->GetBinContent(iy);
    double N_pp_mix = h_rebinned_lspp_mix_1D_pt->GetBinContent(iy);
    double N_mm_mix = h_rebinned_lsmm_mix_1D_pt->GetBinContent(iy);
    double N = N_pm_same;
    hN_1D_pt_rebin_initial     ->SetBinContent(iy, N);

    if (N_pp_mix > 0 && N_mm_mix > 0) {
        double R = N_pm_mix / (2.0 * std::sqrt(N_pp_mix * N_mm_mix));
        double B = 2.0  * R * std::sqrt(N_pp_same * N_mm_same);
        double Z = N_pm_same - B;
        double S_over_B = (B > 0.0) ? Z / B : 0.0;
        hR_1D_pt_rebin_initial     ->SetBinContent(iy, R);
        hZ_1D_pt_rebin_initial     ->SetBinContent(iy, Z);
        hSoverB_1D_pt_rebin_initial->SetBinContent(iy, S_over_B);
        hB_1D_pt_rebin_initial     ->SetBinContent(iy, B);

        double A = N_pm_mix;
        double D = 2.0 * std::sqrt(N_pp_mix * N_mm_mix);
        // --- uncertainties ---
        double sigmaA = std::sqrt(N_pm_mix);
        
        // propagate denominator 
        double sigmaD = D * 0.5 * std::sqrt(
            1.0 / N_pp_mix + 1.0 / N_mm_mix
        );
        
        // now propagate ratio
        double sigmaR = std::abs(R) * std::sqrt(
            (sigmaA / A)*(sigmaA / A) +
            (sigmaD / D)*(sigmaD / D)
        );
        
        hR_1D_pt_rebin_initial->SetBinError(iy, sigmaR);

        // Background uncertainty
        // --- propagate C first ---
        double C = std::sqrt(N_pp_same * N_mm_same);
        double sigmaC = C * 0.5 * std::sqrt(
        1.0 / N_pp_same + 1.0 / N_mm_same
        );

        // --- now combine with R uncertainty ---
        double sigmaB = 0.0;

        if (R > 0 && C > 0 && sigmaR > 0 && sigmaC > 0) {
            sigmaB = B * std::sqrt(
                (sigmaR / R)*(sigmaR / R) +
                (sigmaC / C)*(sigmaC / C)
            );
        }
        if (sigmaB < 0 || std::isnan(sigmaB) || std::isinf(sigmaB)) sigmaB = 0;
        hB_1D_pt_rebin_initial->SetBinError(iy, sigmaB);

        // Z uncertainty
        double varZ = N_pm_same + sigmaB*sigmaB;
        if (varZ < 0 || std::isnan(varZ) || std::isinf(varZ)) varZ = 0;
        hZ_1D_pt_rebin_initial->SetBinError(iy, std::sqrt(varZ));

        // N uncertainty
        hN_1D_pt_rebin_initial->SetBinError(iy,  std::sqrt(N_pm_same));
        
        // Signal over background uncertainty
        double sigmaSB = 0.0;

        if (Z > 0 && B > 0) {
            double sigmaZ = std::sqrt(varZ);
    
            double termZ = sigmaZ / Z;
            double termB = sigmaB / B;
            sigmaSB = (Z/B) * std::sqrt(termZ*termZ + termB*termB);
        if (std::isnan(sigmaSB) || std::isinf(sigmaSB)) sigmaSB = 0.0;
        }
        hSoverB_1D_pt_rebin_initial->SetBinError(iy, sigmaSB);
    }
}
hZ_1D_pt_rebin_initial->Scale(1.0 / multiplicity);
hB_1D_pt_rebin_initial->Scale(1.0 / multiplicity);
hN_1D_pt_rebin_initial->Scale(1.0 / multiplicity);
hR_1D_pt_rebin_initial->Write("rfactor_pt");
hR_1D_mee_rebin->Write("rfactor_mee");
hZ_1D_pt_rebin_initial->Write();
hN_1D_pt_rebin_initial->Write();
hB_1D_pt_rebin_initial->Write();

// ======================================== Rebinning again ==================================== //
// here we rebin the already projected 1D histograms for better visualization and then plot them //
// ======================================= Rebinning again ===================================== // 





std::vector<BinRange> xRanges_visual = {
    {0.0, 0.45, 0.05},   
    {0.45, 1.2, 0.05},   
    {1.2, 2.8, 0.4},    
    {2.8, 3.2, 0.05},    
    {3.2, 4.0, 0.4}      
};
std::vector<BinRange> yRanges_visual = {
    {0.0, 5.0, 0.5},   
    {5.0, 10.0, 1.0}  
};

// --- X-axis bins ---
std::vector<double> xbins_visual;
for (auto &[start, end, step] : xRanges_visual) {
    if (xbins_visual.empty()) xbins_visual.push_back(start); 
    for (double x = start + step; x < end - 1e-4; x += step) {
        xbins_visual.push_back(x);
    }
    xbins_visual.push_back(end); 
}
std::cout << std::endl;
// --- Y-axis bins ---
std::vector<double> ybins_visual;
for (auto &[start, end, step] : yRanges_visual) {
    if (ybins_visual.empty()) ybins_visual.push_back(start);
    for (double y = start + step; y < end - 1e-4; y += step) {
        ybins_visual.push_back(y);
    }
    ybins_visual.push_back(end);
}

TAxis* xax_rebin = hZ_2D_rebinned->GetXaxis();
TAxis* yax_rebin = hZ_2D_rebinned->GetYaxis();

std::vector<std::pair<int,int>> meeSlices2_visual;

for (size_t i = 0; i < xbins_visual.size() - 1; ++i) {

    int lo = FindBinLow2(xax_rebin, xbins_visual[i]);
    int hi = FindBinHigh2(xax_rebin, xbins_visual[i+1]);

    lo = std::max(lo, 1);
    hi = std::min(hi, xax_rebin->GetNbins());

    meeSlices2_visual.emplace_back(lo, hi);
}


std::vector<std::pair<int,int>> ptSlices2_visual;
for (size_t i = 0; i < ybins_visual.size() - 1; ++i) {

    int lo = FindBinLow2(yax_rebin, ybins_visual[i]);
    int hi = FindBinHigh2(yax_rebin, ybins_visual[i+1]);

    lo = std::max(lo, 1);
    hi = std::min(hi, yax_rebin->GetNbins());

    ptSlices2_visual.emplace_back(lo, hi);
}
// =============================== correct 1D rebinning =============================== //
//         rebinning the TH1D origionating form a 2D R factor into 2D Z, B and N        
// =============================== correct 1D rebinning =============================== //
TH1D *hZ_1D_mee_visual = new TH1D("hZ_1D_mee_visual", "Z;M_{ee};Z",
    xbins_visual.size()-1, xbins_visual.data());

TH1D *hB_1D_mee_visual = new TH1D("hB_1D_mee_visual", "B;M_{ee};B",
    xbins_visual.size()-1, xbins_visual.data());

TH1D *hN_1D_mee_visual = new TH1D("hN_1D_mee_visual", "N;M_{ee};N",
    xbins_visual.size()-1, xbins_visual.data());

    TH1D *hZ_1D_pT_visual = new TH1D("hZ_1D_pT_visual", "Z;p_{T};Z",
        ybins_visual.size()-1, ybins_visual.data());
    
    TH1D *hB_1D_pT_visual = new TH1D("hB_1D_pT_visual", "B;p_{T};B",
        ybins_visual.size()-1, ybins_visual.data());
    
    TH1D *hN_1D_pT_visual = new TH1D("hN_1D_pT_visual", "N;p_{T};N",
        ybins_visual.size()-1, ybins_visual.data());

        for (size_t i = 0; i < meeSlices2_visual.size(); ++i) {

            double sumZ = 0, sumB = 0, sumN = 0;
            double errZ2 = 0, errB2 = 0, errN2 = 0;
        
            int xlo = meeSlices2_visual[i].first;
            int xhi = meeSlices2_visual[i].second;
        
            for (int ix = xlo; ix <= xhi; ++ix) {
        
                double Z = hZ_1D_mee_rebin->GetBinContent(ix);
                double B = hB_1D_mee_rebin->GetBinContent(ix);
                double N = hN_1D_mee_rebin->GetBinContent(ix);
        
                double eZ = hZ_1D_mee_rebin->GetBinError(ix);
                double eB = hB_1D_mee_rebin->GetBinError(ix);
                double eN = hN_1D_mee_rebin->GetBinError(ix);
        
                sumZ += Z;
                sumB += B;
                sumN += N;
        
                errZ2 += eZ * eZ;
                errB2 += eB * eB;
                errN2 += eN * eN;
            }
        
            hZ_1D_mee_visual->SetBinContent(i+1, sumZ);
            hB_1D_mee_visual->SetBinContent(i+1, sumB);
            hN_1D_mee_visual->SetBinContent(i+1, sumN);
        
            hZ_1D_mee_visual->SetBinError(i+1, sqrt(errZ2));
            hB_1D_mee_visual->SetBinError(i+1, sqrt(errB2));
            hN_1D_mee_visual->SetBinError(i+1, sqrt(errN2));
        }
        for (size_t i = 0; i < ptSlices2_visual.size(); ++i) {

            double sumZ = 0, sumB = 0, sumN = 0;
            double errZ2 = 0, errB2 = 0, errN2 = 0;
        
            int ylo = ptSlices2_visual[i].first;
            int yhi = ptSlices2_visual[i].second;
        
            for (int iy = ylo; iy <= yhi; ++iy) {
        
                double Z = hZ_1D_pt_rebin_initial->GetBinContent(iy);
                double B = hB_1D_pt_rebin_initial->GetBinContent(iy);
                double N = hN_1D_pt_rebin_initial->GetBinContent(iy);
        
                double eZ = hZ_1D_pt_rebin_initial->GetBinError(iy);
                double eB = hB_1D_pt_rebin_initial->GetBinError(iy);
                double eN = hN_1D_pt_rebin_initial->GetBinError(iy);
        
                sumZ += Z;
                sumB += B;
                sumN += N;
        
                errZ2 += eZ * eZ;
                errB2 += eB * eB;
                errN2 += eN * eN;
            }
        
            hZ_1D_pT_visual->SetBinContent(i+1, sumZ);
            hB_1D_pT_visual->SetBinContent(i+1, sumB);
            hN_1D_pT_visual->SetBinContent(i+1, sumN);
        
            hZ_1D_pT_visual->SetBinError(i+1, sqrt(errZ2));
            hB_1D_pT_visual->SetBinError(i+1, sqrt(errB2));
            hN_1D_pT_visual->SetBinError(i+1, sqrt(errN2));
        }
        for (int i = 1; i <= hZ_1D_mee_visual->GetNbinsX(); ++i) {

            double w = hZ_1D_mee_visual->GetBinWidth(i);
            if (w <= 0) continue;
        
            hZ_1D_mee_visual->SetBinContent(i, hZ_1D_mee_visual->GetBinContent(i)/w);
            hB_1D_mee_visual->SetBinContent(i, hB_1D_mee_visual->GetBinContent(i)/w);
            hN_1D_mee_visual->SetBinContent(i, hN_1D_mee_visual->GetBinContent(i)/w);
        
            hZ_1D_mee_visual->SetBinError(i, hZ_1D_mee_visual->GetBinError(i)/w);
            hB_1D_mee_visual->SetBinError(i, hB_1D_mee_visual->GetBinError(i)/w);
            hN_1D_mee_visual->SetBinError(i, hN_1D_mee_visual->GetBinError(i)/w);
        }
        for (int i = 1; i <= hZ_1D_pT_visual->GetNbinsX(); ++i) {

            double w = hZ_1D_pT_visual->GetBinWidth(i);
            if (w <= 0) continue;
        
            hZ_1D_pT_visual->SetBinContent(i, hZ_1D_pT_visual->GetBinContent(i)/w);
            hB_1D_pT_visual->SetBinContent(i, hB_1D_pT_visual->GetBinContent(i)/w);
            hN_1D_pT_visual->SetBinContent(i, hN_1D_pT_visual->GetBinContent(i)/w);
        
            hZ_1D_pT_visual->SetBinError(i, hZ_1D_pT_visual->GetBinError(i)/w);
            hB_1D_pT_visual->SetBinError(i, hB_1D_pT_visual->GetBinError(i)/w);
            hN_1D_pT_visual->SetBinError(i, hN_1D_pT_visual->GetBinError(i)/w);
        }
        hN_1D_mee_visual->Write("hN_1D_mee_visual");
hZ_1D_mee_visual->Write("hZ_1D_mee_visual");
hB_1D_mee_visual->Write("hB_1D_mee_visual");
hN_1D_pT_visual->Write("hN_1D_pT_visual");
hZ_1D_pT_visual->Write("hZ_1D_pT_visual");
hB_1D_pT_visual->Write("hB_1D_pT_visual");

TCanvas *RfactorpT = new TCanvas("c","1D plot",800,700);

RfactorpT->SetTicks(1,1);
RfactorpT->SetRightMargin(0.15);  // space for color bar
RfactorpT->SetLeftMargin(0.12);
RfactorpT->SetBottomMargin(0.12);

hR_1D_pt_rebin_initial->SetStats(0);
hR_1D_pt_rebin_initial->SetMarkerStyle(24);

// set Z axis range

hR_1D_pt_rebin_initial->GetXaxis()->SetTitle("p_{T} (GeV/c)");
hR_1D_pt_rebin_initial->GetYaxis()->SetTitle("R");
hR_1D_pt_rebin_initial->SetTitle("");

// axis styling (optional but nicer)
hR_1D_pt_rebin_initial->GetXaxis()->SetTitleSize(0.045);
hR_1D_pt_rebin_initial->GetYaxis()->SetTitleSize(0.045);
hR_1D_pt_rebin_initial->GetXaxis()->SetLabelSize(0.04);
hR_1D_pt_rebin_initial->GetYaxis()->SetLabelSize(0.04);
hR_1D_pt_rebin_initial->GetYaxis()->SetRangeUser(0.85, 1.22);

// draw
hR_1D_pt_rebin_initial->Draw("COLZ");
latex.SetNDC();
latex.SetTextSize(0.036);
latex.SetTextAlign(13);

latex.DrawLatex(0.22, 0.88, "O-O, #sqrt{s_{NN}} = 5.36 TeV");
latex.DrawLatex(0.22, 0.846, "0-90%");
latex.DrawLatex(0.22, 0.82, "p_{T} > 0.2 GeV/c, |#eta| < 0.8");

RfactorpT->SaveAs("hR_1D_rebinned_pT_OO.pdf");

TCanvas *Rfactormee = new TCanvas("c","1D plot mee",800,700);
Rfactormee->SetTicks(1,1);
Rfactormee->SetRightMargin(0.15);  // space for color bar
Rfactormee->SetLeftMargin(0.12);
Rfactormee->SetBottomMargin(0.12);


hR_1D_mee_rebin->SetStats(0);
hR_1D_mee_rebin->SetMarkerStyle(24);
hR_1D_mee_rebin->GetXaxis()->SetRangeUser(0.0, 1.0);

// set Z axis range

hR_1D_mee_rebin->GetXaxis()->SetTitle("m_{ee} (GeV/c^{2})");
hR_1D_mee_rebin->GetYaxis()->SetTitle("R");
hR_1D_mee_rebin->SetTitle("");

// axis styling (optional but nicer)
hR_1D_mee_rebin->GetXaxis()->SetTitleSize(0.045);
hR_1D_mee_rebin->GetYaxis()->SetTitleSize(0.045);
hR_1D_mee_rebin->GetXaxis()->SetLabelSize(0.04);
hR_1D_mee_rebin->GetYaxis()->SetLabelSize(0.04);
hR_1D_mee_rebin->GetYaxis()->SetRangeUser(0.94, 1.01);

// draw
hR_1D_mee_rebin->Draw("COLZ");
latex.SetNDC();
latex.SetTextSize(0.036);
latex.SetTextAlign(13);

latex.DrawLatex(0.5, 0.38, "O-O, #sqrt{s_{NN}} = 5.36 TeV");
latex.DrawLatex(0.5, 0.346, "0-90%");
latex.DrawLatex(0.5, 0.32, "p_{T} > 0.2 GeV/c, |#eta| < 0.8");

Rfactormee->SaveAs("hR_1D_rebinned_mee_OO.pdf");

TH1D *hZ_projX_rebinned = new TH1D("hZ_projX_rebinned","Z;X;Z",xbins_visual.size()-1,xbins_visual.data());
TH1D *hB_projX_rebinned = new TH1D("hB_projX_rebinned","B;X;B",xbins_visual.size()-1,xbins_visual.data());
TH1D *hN_projX_rebinned = new TH1D("hN_projX_rebinned","N;X;N",xbins_visual.size()-1,xbins_visual.data());
TH1D *hZ_projY_rebinned = new TH1D("hZ_projY_rebinned","Z;Y;Z",ybins_visual.size()-1,ybins_visual.data());
TH1D *hB_projY_rebinned = new TH1D("hB_projY_rebinned","B;Y;B",ybins_visual.size()-1,ybins_visual.data());
TH1D *hN_projY_rebinned = new TH1D("hN_projY_rebinned","N;Y;N",ybins_visual.size()-1,ybins_visual.data());


for (size_t i = 0; i < meeSlices2_visual.size(); ++i) {

    double sumZ = 0, sumB = 0, sumN = 0;
    double errZ2 = 0, errB2 = 0, errN2 = 0;

    int xlo = meeSlices2_visual[i].first;
    int xhi = meeSlices2_visual[i].second;

    for (int ix = xlo; ix <= xhi; ++ix) {


            double Z = hZ_projX->GetBinContent(ix);
            double B = hB_projX->GetBinContent(ix);
            double N = hN_projX->GetBinContent(ix);

            double eZ = hZ_projX->GetBinError(ix);
            double eB = hB_projX->GetBinError(ix);
            double eN = hN_projX->GetBinError(ix);

            sumZ += Z;  sumB += B;  sumN += N;
            errZ2 += eZ * eZ;
            errB2 += eB * eB;
            errN2 += eN * eN;
        
    }

    hZ_projX_rebinned->SetBinContent(i + 1, sumZ);
    hB_projX_rebinned->SetBinContent(i + 1, sumB);
    hN_projX_rebinned->SetBinContent(i + 1, sumN);

    hZ_projX_rebinned->SetBinError(i + 1, std::sqrt(errZ2));
    hB_projX_rebinned->SetBinError(i + 1, std::sqrt(errB2));
    hN_projX_rebinned->SetBinError(i + 1, std::sqrt(errN2));
}


for (size_t i = 0; i < ptSlices2_visual.size(); ++i) {

    double sumZ = 0, sumB = 0, sumN = 0;
    double errZ2 = 0, errB2 = 0, errN2 = 0;

    int ylo = ptSlices2_visual[i].first;
    int yhi = ptSlices2_visual[i].second;


    for (int iy = ylo; iy <= yhi; ++iy) {


            double Z = hZ_projY->GetBinContent(iy);
            double B = hB_projY->GetBinContent(iy);
            double N = hN_projY->GetBinContent(iy);

            double eZ = hZ_projY->GetBinError(iy);
            double eB = hB_projY->GetBinError(iy);
            double eN = hN_projY->GetBinError(iy);

            sumZ += Z;  sumB += B;  sumN += N;
            errZ2 += eZ * eZ;
            errB2 += eB * eB;
            errN2 += eN * eN;
        
    }

    hZ_projY_rebinned->SetBinContent(i + 1, sumZ);
    hB_projY_rebinned->SetBinContent(i + 1, sumB);
    hN_projY_rebinned->SetBinContent(i + 1, sumN);

    hZ_projY_rebinned->SetBinError(i + 1, std::sqrt(errZ2));
    hB_projY_rebinned->SetBinError(i + 1, std::sqrt(errB2));
    hN_projY_rebinned->SetBinError(i + 1, std::sqrt(errN2));
}

for (int i = 1; i <= hZ_projX_rebinned->GetNbinsX(); ++i) {
    double width = hZ_projX_rebinned->GetBinWidth(i);

    hZ_projX_rebinned->SetBinContent(i, hZ_projX_rebinned->GetBinContent(i) / width);
    hB_projX_rebinned->SetBinContent(i, hB_projX_rebinned->GetBinContent(i) / width);
    hN_projX_rebinned->SetBinContent(i, hN_projX_rebinned->GetBinContent(i) / width);

    hZ_projX_rebinned->SetBinError(i, hZ_projX_rebinned->GetBinError(i) / width);
    hB_projX_rebinned->SetBinError(i, hB_projX_rebinned->GetBinError(i) / width);
    hN_projX_rebinned->SetBinError(i, hN_projX_rebinned->GetBinError(i) / width);
}
for (int i = 1; i <= hZ_projY_rebinned->GetNbinsX(); ++i) {
    double width = hZ_projY_rebinned->GetBinWidth(i);

    hZ_projY_rebinned->SetBinContent(i, hZ_projY_rebinned->GetBinContent(i) / width);
    hB_projY_rebinned->SetBinContent(i, hB_projY_rebinned->GetBinContent(i) / width);
    hN_projY_rebinned->SetBinContent(i, hN_projY_rebinned->GetBinContent(i) / width);

    hZ_projY_rebinned->SetBinError(i, hZ_projY_rebinned->GetBinError(i) / width);
    hB_projY_rebinned->SetBinError(i, hB_projY_rebinned->GetBinError(i) / width);
    hN_projY_rebinned->SetBinError(i, hN_projY_rebinned->GetBinError(i) / width);
}
hZ_projX_rebinned->Write("hZ_mass_using2D_rebinned");
hN_projX_rebinned->Write("hN_mass_using2D_rebinned");
hB_projX_rebinned->Write("hB_mass_using2D_rebinned");
hZ_projY_rebinned->Write("hZ_pT_using2D_rebinned");
hN_projY_rebinned->Write("hN_pT_using2D_rebinned");
hB_projY_rebinned->Write("hB_pT_using2D_rebinned");

// ========================================= applying 2D R factor onto 3D LS ========================================= //
TH3D* hlspp_corr_3D = (TH3D*)h_rebinned_lspp_3D->Clone("hlspp_corr_3D");
hlspp_corr_3D->Reset();
TH3D* hlsmm_corr_3D = (TH3D*)h_rebinned_lsmm_3D->Clone("hlsmm_corr_3D");
hlsmm_corr_3D->Reset();

for (int ix = 1; ix <= h_rebinned_lspp_3D->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= h_rebinned_lspp_3D->GetNbinsY(); ++iy) {

        double R  = hR_2D_rebinned->GetBinContent(ix, iy);
        double eR = hR_2D_rebinned->GetBinError(ix, iy);

        for (int iz = 1; iz <= h_rebinned_lspp_3D->GetNbinsZ(); ++iz) {

            double Npp  = h_rebinned_lspp_3D->GetBinContent(ix, iy, iz);
            double eNpp = h_rebinned_lspp_3D->GetBinError(ix, iy, iz);

            double Nmm  = h_rebinned_lsmm_3D->GetBinContent(ix, iy, iz);
            double eNmm = h_rebinned_lsmm_3D->GetBinError(ix, iy, iz);

            double Bpp  = Npp;   // default: keep original value
            double Bmm  = Nmm;

            double eBpp = eNpp;  // default: keep original error
            double eBmm = eNmm;

            // Only apply correction if R is nonzero
            if (R != 0.0) {

                Bpp = R * Npp;
                Bmm = R * Nmm;

                eBpp = std::sqrt(
                    (Npp * eR) * (Npp * eR) +
                    (R * eNpp) * (R * eNpp)
                );

                eBmm = std::sqrt(
                    (Nmm * eR) * (Nmm * eR) +
                    (R * eNmm) * (R * eNmm)
                );
            }

            hlspp_corr_3D->SetBinContent(ix, iy, iz, Bpp);
            hlspp_corr_3D->SetBinError(ix, iy, iz, eBpp);

            hlsmm_corr_3D->SetBinContent(ix, iy, iz, Bmm);
            hlsmm_corr_3D->SetBinError(ix, iy, iz, eBmm);
        }
    }
}
hlsmm_corr_3D->Write("hlsmm_corr_3D");
hlspp_corr_3D->Write("hlspp_corr_3D");
TH2D* h_lspp_corr_2D_M_pT = (TH2D*)hlspp_corr_3D->Project3D("yx");
TH2D* h_lsmm_corr_2D_M_pT = (TH2D*)hlsmm_corr_3D->Project3D("yx");
TH2D* h_lspp_corr_2D_M_DCA = (TH2D*)hlspp_corr_3D->Project3D("zx");
TH2D* h_lsmm_corr_2D_M_DCA = (TH2D*)hlsmm_corr_3D->Project3D("zx");
TH2D* h_lspp_corr_2D_pT_DCA = (TH2D*)hlspp_corr_3D->Project3D("zy");
TH2D* h_lsmm_corr_2D_pT_DCA = (TH2D*)hlsmm_corr_3D->Project3D("zy");
h_lspp_corr_2D_pT_DCA->Write("h_lspp_corr_2D_pT_DCA");
h_lsmm_corr_2D_pT_DCA->Write("h_lsmm_corr_2D_pT_DCA");

TH2D* h_ULS_2D_M_pT = (TH2D*)h_rebinned_uls_3D->Project3D("yx");
TH2D* h_ULS_2D_M_DCA = (TH2D*)h_rebinned_uls_3D->Project3D("zx");
TH2D* h_ULS_2D_pT_DCA = (TH2D*)h_rebinned_uls_3D->Project3D("zy");

TH2D* hZ_2D_rebinned_M_pT = new TH2D("hZ_2D_rebinned_M_pT","Z;M;p_{T}",h_ULS_2D_M_pT->GetNbinsX(),h_ULS_2D_M_pT->GetXaxis()->GetXbins()->GetArray(),h_ULS_2D_M_pT->GetNbinsY(),h_ULS_2D_M_pT->GetYaxis()->GetXbins()->GetArray());
TH2D* hB_2D_rebinned_M_pT = new TH2D("hB_2D_rebinned_M_pT","Z;M;p_{T}",h_ULS_2D_M_pT->GetNbinsX(),h_ULS_2D_M_pT->GetXaxis()->GetXbins()->GetArray(),h_ULS_2D_M_pT->GetNbinsY(),h_ULS_2D_M_pT->GetYaxis()->GetXbins()->GetArray());

TH2D* hZ_2D_rebinned_M_DCA = new TH2D("hZ_2D_rebinned_M_DCA","Z;M;DCA",h_ULS_2D_M_DCA->GetNbinsX(),h_ULS_2D_M_DCA->GetXaxis()->GetXbins()->GetArray(),h_ULS_2D_M_DCA->GetNbinsY(),h_ULS_2D_M_DCA->GetYaxis()->GetXbins()->GetArray());
TH2D* hB_2D_rebinned_M_DCA = new TH2D("hB_2D_rebinned_M_DCA","Z;M;DCA",h_ULS_2D_M_DCA->GetNbinsX(),h_ULS_2D_M_DCA->GetXaxis()->GetXbins()->GetArray(),h_ULS_2D_M_DCA->GetNbinsY(),h_ULS_2D_M_DCA->GetYaxis()->GetXbins()->GetArray());

TH2D* hZ_2D_rebinned_pT_DCA = new TH2D("hZ_2D_rebinned_pT_DCA","Z;M;DCA",h_ULS_2D_pT_DCA->GetNbinsX(),h_ULS_2D_pT_DCA->GetXaxis()->GetXbins()->GetArray(),h_ULS_2D_pT_DCA->GetNbinsY(),h_ULS_2D_pT_DCA->GetYaxis()->GetXbins()->GetArray());
TH2D* hB_2D_rebinned_pT_DCA = new TH2D("hB_2D_rebinned_pT_DCA","Z;M;DCA",h_ULS_2D_pT_DCA->GetNbinsX(),h_ULS_2D_pT_DCA->GetXaxis()->GetXbins()->GetArray(),h_ULS_2D_pT_DCA->GetNbinsY(),h_ULS_2D_pT_DCA->GetYaxis()->GetXbins()->GetArray());

for (int ix = 1; ix <= h_lspp_corr_2D_M_pT->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= h_lspp_corr_2D_M_pT->GetNbinsY(); ++iy) {

        double Bpp = h_lspp_corr_2D_M_pT->GetBinContent(ix, iy);
        double Bmm = h_lsmm_corr_2D_M_pT->GetBinContent(ix, iy);

        double ULS = h_ULS_2D_M_pT->GetBinContent(ix, iy);
        double B = 2 * sqrt(Bpp * Bmm);
        double Z = ULS - B;
        hZ_2D_rebinned_M_pT->SetBinContent(ix, iy, Z);
        hB_2D_rebinned_M_pT->SetBinContent(ix, iy, B);
    }
}

for (int ix = 1; ix <= h_lspp_corr_2D_M_DCA->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= h_lspp_corr_2D_M_DCA->GetNbinsY(); ++iy) {

        double Bpp = h_lspp_corr_2D_M_DCA->GetBinContent(ix, iy);
        double Bmm = h_lsmm_corr_2D_M_DCA->GetBinContent(ix, iy);

        double ULS = h_ULS_2D_M_DCA->GetBinContent(ix, iy);
        double B = 2 * sqrt(Bpp * Bmm);
        double Z = ULS - B;
        hZ_2D_rebinned_M_DCA->SetBinContent(ix, iy, Z);
        hB_2D_rebinned_M_DCA->SetBinContent(ix, iy, B);
    }
}

for (int ix = 1; ix <= h_lspp_corr_2D_pT_DCA->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= h_lspp_corr_2D_pT_DCA->GetNbinsY(); ++iy) {

        double Bpp = h_lspp_corr_2D_pT_DCA->GetBinContent(ix, iy);
        double Bmm = h_lsmm_corr_2D_pT_DCA->GetBinContent(ix, iy);

        double ULS = h_ULS_2D_pT_DCA->GetBinContent(ix, iy);
        double B = 2 * sqrt(Bpp * Bmm);
        double Z = ULS - B;
        hZ_2D_rebinned_pT_DCA->SetBinContent(ix, iy, Z);
        hB_2D_rebinned_pT_DCA->SetBinContent(ix, iy, B);
    }
}
hZ_2D_rebinned_pT_DCA->Write("hZ_2D_rebinned_pT_DCA");
hB_2D_rebinned_pT_DCA->Write("hB_2D_rebinned_pT_DCA");

TH1D* hZ_1D_M_pT = (TH1D*)hZ_2D_rebinned_M_pT->ProjectionX();
TH1D* hZ_1D_M_DCA = (TH1D*)hZ_2D_rebinned_M_DCA->ProjectionX();
TH1D* hB_1D_M_pT = (TH1D*)hB_2D_rebinned_M_pT->ProjectionX();
TH1D* hB_1D_M_DCA = (TH1D*)hB_2D_rebinned_M_DCA->ProjectionX();
TH1D* hN_1D_M_pT = (TH1D*)h_ULS_2D_M_pT->ProjectionX();
TH1D* hN_1D_M_DCA = (TH1D*)h_ULS_2D_M_DCA->ProjectionX();
TH1D* hZ_1D_pT_DCA = (TH1D*)hZ_2D_rebinned_pT_DCA->ProjectionX();
TH1D* hB_1D_pT_DCA = (TH1D*)hB_2D_rebinned_pT_DCA->ProjectionX();
TH1D* hN_1D_pT_DCA = (TH1D*)h_ULS_2D_pT_DCA->ProjectionX();
TH1D* hZ_1D_pT_M = (TH1D*)hZ_2D_rebinned_M_pT->ProjectionY();
TH1D* hB_1D_pT_M = (TH1D*)hB_2D_rebinned_M_pT->ProjectionY();
TH1D* hN_1D_pT_M = (TH1D*)h_ULS_2D_M_pT->ProjectionY();

hZ_1D_M_pT->Scale(norm);
hZ_1D_M_DCA->Scale(norm);
hB_1D_M_pT->Scale(norm);
hB_1D_M_DCA->Scale(norm);
hN_1D_M_pT->Scale(norm);
hN_1D_M_DCA->Scale(norm);
hZ_1D_pT_DCA->Scale(norm);
hB_1D_pT_DCA->Scale(norm);
hN_1D_pT_DCA->Scale(norm);
hZ_1D_pT_M->Scale(norm);
hB_1D_pT_M->Scale(norm);
hN_1D_pT_M->Scale(norm);

hZ_2D_rebinned_M_pT->Write("hZ_2D_M_pT");
hZ_2D_rebinned_M_DCA->Write("hZ_2D_M_DCA");
hB_2D_rebinned_M_pT->Write("hB_2D_M_pT");
hB_2D_rebinned_M_DCA->Write("hB_2D_M_DCA");
hZ_2D_rebinned_pT_DCA->Write("hZ_2D_pT_DCA");
hB_2D_rebinned_pT_DCA->Write("hB_2D_pT_DCA");
h_ULS_2D_M_pT->Write("hN_2D_M_pT");
h_ULS_2D_M_DCA->Write("hN_2D_M_DCA");

hZ_1D_M_pT->Write("hZ_1D_M_pT");
hZ_1D_M_DCA->Write("hZ_1D_M_DCA");
hB_1D_M_pT->Write("hB_1D_M_pT");
hB_1D_M_DCA->Write("hB_1D_M_DCA");
hN_1D_M_pT->Write("hN_1D_M_pT");
hN_1D_M_DCA->Write("hN_1D_M_DCA");

hZ_1D_pT_DCA->Write("hZ_1D_pT_DCA");
hB_1D_pT_DCA->Write("hB_1D_pT_DCA");
hN_1D_pT_DCA->Write("hN_1D_pT_DCA");
hZ_1D_pT_M->Write("hZ_1D_pT_M");
hB_1D_pT_M->Write("hB_1D_pT_M");
hN_1D_pT_M->Write("hN_1D_pT_M");


TH1D *hZ_1D_rebinned_M_pT  = new TH1D("hZ_1D_rebinned_M_pT",  "Z;Y;Z", xbins_visual.size()-1, xbins_visual.data());
TH1D *hZ_1D_rebinned_M_DCA = new TH1D("hZ_1D_rebinned_M_DCA", "Z;X;Z", xbins_visual.size()-1, xbins_visual.data());
TH1D *hB_1D_rebinned_M_pT  = new TH1D("hB_1D_rebinned_M_pT",  "B;Y;B", xbins_visual.size()-1, xbins_visual.data());
TH1D *hB_1D_rebinned_M_DCA = new TH1D("hB_1D_rebinned_M_DCA", "B;X;B", xbins_visual.size()-1, xbins_visual.data());
TH1D *hN_1D_rebinned_M_pT  = new TH1D("hN_1D_rebinned_M_pT",  "N;Y;N", xbins_visual.size()-1, xbins_visual.data());
TH1D *hN_1D_rebinned_M_DCA = new TH1D("hN_1D_rebinned_M_DCA", "N;X;N", xbins_visual.size()-1, xbins_visual.data());

TH1D *hZ_1D_rebinned_pT_DCA = new TH1D("hZ_1D_rebinned_pT_DCA", "Z;Y;Z", ybins_visual.size()-1, ybins_visual.data());
TH1D *hB_1D_rebinned_pT_DCA = new TH1D("hB_1D_rebinned_pT_DCA", "B;Y;B", ybins_visual.size()-1, ybins_visual.data());
TH1D *hN_1D_rebinned_pT_DCA = new TH1D("hN_1D_rebinned_pT_DCA", "N;Y;N", ybins_visual.size()-1, ybins_visual.data());
TH1D *hZ_1D_rebinned_pT_M = new TH1D("hZ_1D_rebinned_pT_M", "Z;X;Z", ybins_visual.size()-1, ybins_visual.data());
TH1D *hB_1D_rebinned_pT_M = new TH1D("hB_1D_rebinned_pT_M", "B;X;B", ybins_visual.size()-1, ybins_visual.data());
TH1D *hN_1D_rebinned_pT_M = new TH1D("hN_1D_rebinned_pT_M", "N;X;N", ybins_visual.size()-1, ybins_visual.data());

for (size_t i = 0; i < meeSlices2_visual.size(); ++i) {

    double sumZ_pT = 0, sumB_pT = 0, sumN_pT = 0;
    double errZ2_pT = 0, errB2_pT = 0, errN2_pT = 0;
    double sumZ_dca = 0, sumB_dca = 0, sumN_dca = 0;
    double errZ2_dca = 0, errB2_dca = 0, errN2_dca = 0;

    int xlo = meeSlices2_visual[i].first;
    int xhi = meeSlices2_visual[i].second;

    for (int ix = xlo; ix <= xhi; ++ix) {

        double Z_pT = hZ_1D_M_pT->GetBinContent(ix);
        double B_pT = hB_1D_M_pT->GetBinContent(ix);
        double N_pT = hN_1D_M_pT->GetBinContent(ix);
        double Z_dca = hZ_1D_M_DCA->GetBinContent(ix);
        double B_dca = hB_1D_M_DCA->GetBinContent(ix);
        double N_dca = hN_1D_M_DCA->GetBinContent(ix);

        double eZ_dca = hZ_1D_M_DCA->GetBinError(ix);
        double eB_dca = hB_1D_M_DCA->GetBinError(ix);
        double eN_dca = hN_1D_M_DCA->GetBinError(ix);
        double eZ_pT = hZ_1D_M_pT->GetBinError(ix);
        double eB_pT = hB_1D_M_pT->GetBinError(ix);
        double eN_pT = hN_1D_M_pT->GetBinError(ix);

        sumZ_pT += Z_pT;
        sumB_pT += B_pT;
        sumN_pT += N_pT;
        sumZ_dca += Z_dca;
        sumB_dca += B_dca;
        sumN_dca += N_dca;

        errZ2_dca += eZ_dca * eZ_dca;
        errB2_dca += eB_dca * eB_dca;
        errN2_dca += eN_dca * eN_dca;
        errZ2_pT += eZ_pT * eZ_pT;
        errB2_pT += eB_pT * eB_pT;
        errN2_pT += eN_pT * eN_pT;
    }

    hZ_1D_rebinned_M_pT->SetBinContent(i + 1, sumZ_pT);
    hB_1D_rebinned_M_pT->SetBinContent(i + 1, sumB_pT);
    hN_1D_rebinned_M_pT->SetBinContent(i + 1, sumN_pT);

    hZ_1D_rebinned_M_pT->SetBinError(i + 1, std::sqrt(errZ2_pT));
    hB_1D_rebinned_M_pT->SetBinError(i + 1, std::sqrt(errB2_pT));
    hN_1D_rebinned_M_pT->SetBinError(i + 1, std::sqrt(errN2_pT));

    hZ_1D_rebinned_M_DCA->SetBinContent(i + 1, sumZ_dca);
    hB_1D_rebinned_M_DCA->SetBinContent(i + 1, sumB_dca);
    hN_1D_rebinned_M_DCA->SetBinContent(i + 1, sumN_dca);

    hZ_1D_rebinned_M_DCA->SetBinError(i + 1, std::sqrt(errZ2_dca));
    hB_1D_rebinned_M_DCA->SetBinError(i + 1, std::sqrt(errB2_dca));
    hN_1D_rebinned_M_DCA->SetBinError(i + 1, std::sqrt(errN2_dca));
}

for (size_t i = 0; i < ptSlices2_visual.size(); ++i) {

    double sumZ_M = 0, sumB_M = 0, sumN_M = 0;
    double errZ2_M = 0, errB2_M = 0, errN2_M = 0;

    double sumZ_DCA = 0, sumB_DCA = 0, sumN_DCA = 0;
    double errZ2_DCA = 0, errB2_DCA = 0, errN2_DCA = 0;

    int ylo = ptSlices2_visual[i].first;
    int yhi = ptSlices2_visual[i].second;

    for (int iy = ylo; iy <= yhi; ++iy) {

        // --- M case ---
        double ZM = hZ_1D_pT_M->GetBinContent(iy);
        double BM = hB_1D_pT_M->GetBinContent(iy);
        double NM = hN_1D_pT_M->GetBinContent(iy);

        double eZM = hZ_1D_pT_M->GetBinError(iy);
        double eBM = hB_1D_pT_M->GetBinError(iy);
        double eNM = hN_1D_pT_M->GetBinError(iy);

        sumZ_M += ZM;
        sumB_M += BM;
        sumN_M += NM;

        errZ2_M += eZM * eZM;
        errB2_M += eBM * eBM;
        errN2_M += eNM * eNM;

        // --- DCA case ---
        double ZD = hZ_1D_pT_DCA->GetBinContent(iy);
        double BD = hB_1D_pT_DCA->GetBinContent(iy);
        double ND = hN_1D_pT_DCA->GetBinContent(iy);

        double eZD = hZ_1D_pT_DCA->GetBinError(iy);
        double eBD = hB_1D_pT_DCA->GetBinError(iy);
        double eND = hN_1D_pT_DCA->GetBinError(iy);

        sumZ_DCA += ZD;
        sumB_DCA += BD;
        sumN_DCA += ND;

        errZ2_DCA += eZD * eZD;
        errB2_DCA += eBD * eBD;
        errN2_DCA += eND * eND;
    }

    // --- fill M ---
    hZ_1D_rebinned_pT_M->SetBinContent(i + 1, sumZ_M);
    hB_1D_rebinned_pT_M->SetBinContent(i + 1, sumB_M);
    hN_1D_rebinned_pT_M->SetBinContent(i + 1, sumN_M);

    hZ_1D_rebinned_pT_M->SetBinError(i + 1, std::sqrt(errZ2_M));
    hB_1D_rebinned_pT_M->SetBinError(i + 1, std::sqrt(errB2_M));
    hN_1D_rebinned_pT_M->SetBinError(i + 1, std::sqrt(errN2_M));

    // --- fill DCA ---
    hZ_1D_rebinned_pT_DCA->SetBinContent(i + 1, sumZ_DCA);
    hB_1D_rebinned_pT_DCA->SetBinContent(i + 1, sumB_DCA);
    hN_1D_rebinned_pT_DCA->SetBinContent(i + 1, sumN_DCA);

    hZ_1D_rebinned_pT_DCA->SetBinError(i + 1, std::sqrt(errZ2_DCA));
    hB_1D_rebinned_pT_DCA->SetBinError(i + 1, std::sqrt(errB2_DCA));
    hN_1D_rebinned_pT_DCA->SetBinError(i + 1, std::sqrt(errN2_DCA));
}

for (int i = 1; i <= hZ_1D_rebinned_M_pT->GetNbinsX(); ++i) {

    double width = hZ_1D_rebinned_M_pT->GetBinWidth(i);
    if (width <= 0) continue;

    hZ_1D_rebinned_M_pT->SetBinContent(i, hZ_1D_rebinned_M_pT->GetBinContent(i) / width);
    hB_1D_rebinned_M_pT->SetBinContent(i, hB_1D_rebinned_M_pT->GetBinContent(i) / width);
    hN_1D_rebinned_M_pT->SetBinContent(i, hN_1D_rebinned_M_pT->GetBinContent(i) / width);

    hZ_1D_rebinned_M_pT->SetBinError(i, hZ_1D_rebinned_M_pT->GetBinError(i) / width);
    hB_1D_rebinned_M_pT->SetBinError(i, hB_1D_rebinned_M_pT->GetBinError(i) / width);
    hN_1D_rebinned_M_pT->SetBinError(i, hN_1D_rebinned_M_pT->GetBinError(i) / width);
}

for (int i = 1; i <= hZ_1D_rebinned_M_DCA->GetNbinsX(); ++i) {

    double width = hZ_1D_rebinned_M_DCA->GetBinWidth(i);
    if (width <= 0) continue;

    hZ_1D_rebinned_M_DCA->SetBinContent(i, hZ_1D_rebinned_M_DCA->GetBinContent(i) / width);
    hB_1D_rebinned_M_DCA->SetBinContent(i, hB_1D_rebinned_M_DCA->GetBinContent(i) / width);
    hN_1D_rebinned_M_DCA->SetBinContent(i, hN_1D_rebinned_M_DCA->GetBinContent(i) / width);

    hZ_1D_rebinned_M_DCA->SetBinError(i, hZ_1D_rebinned_M_DCA->GetBinError(i) / width);
    hB_1D_rebinned_M_DCA->SetBinError(i, hB_1D_rebinned_M_DCA->GetBinError(i) / width);
    hN_1D_rebinned_M_DCA->SetBinError(i, hN_1D_rebinned_M_DCA->GetBinError(i) / width);
}
for (int i = 1; i <= hZ_1D_rebinned_pT_M->GetNbinsX(); ++i) {

    double width = hZ_1D_rebinned_pT_M->GetBinWidth(i);
    if (width <= 0) continue;

    hZ_1D_rebinned_pT_M->SetBinContent(i, hZ_1D_rebinned_pT_M->GetBinContent(i) / width);
    hB_1D_rebinned_pT_M->SetBinContent(i, hB_1D_rebinned_pT_M->GetBinContent(i) / width);
    hN_1D_rebinned_pT_M->SetBinContent(i, hN_1D_rebinned_pT_M->GetBinContent(i) / width);

    hZ_1D_rebinned_pT_M->SetBinError(i, hZ_1D_rebinned_pT_M->GetBinError(i) / width);
    hB_1D_rebinned_pT_M->SetBinError(i, hB_1D_rebinned_pT_M->GetBinError(i) / width);
    hN_1D_rebinned_pT_M->SetBinError(i, hN_1D_rebinned_pT_M->GetBinError(i) / width);
}
for (int i = 1; i <= hZ_1D_rebinned_pT_DCA->GetNbinsX(); ++i) {

    double width = hZ_1D_rebinned_pT_DCA->GetBinWidth(i);
    if (width <= 0) continue;

    hZ_1D_rebinned_pT_DCA->SetBinContent(i, hZ_1D_rebinned_pT_DCA->GetBinContent(i) / width);
    hB_1D_rebinned_pT_DCA->SetBinContent(i, hB_1D_rebinned_pT_DCA->GetBinContent(i) / width);
    hN_1D_rebinned_pT_DCA->SetBinContent(i, hN_1D_rebinned_pT_DCA->GetBinContent(i) / width);

    hZ_1D_rebinned_pT_DCA->SetBinError(i, hZ_1D_rebinned_pT_DCA->GetBinError(i) / width);
    hB_1D_rebinned_pT_DCA->SetBinError(i, hB_1D_rebinned_pT_DCA->GetBinError(i) / width);
    hN_1D_rebinned_pT_DCA->SetBinError(i, hN_1D_rebinned_pT_DCA->GetBinError(i) / width);
}

hZ_1D_rebinned_M_DCA->Write("hZ_1D_rebinned_M_DCA");
hZ_1D_rebinned_M_pT->Write("hZ_1D_rebinned_M_pT");
hB_1D_rebinned_M_DCA->Write("hB_1D_rebinned_M_DCA");
hB_1D_rebinned_M_pT->Write("hB_1D_rebinned_M_pT");
hN_1D_rebinned_M_DCA->Write("hN_1D_rebinned_M_DCA");
hN_1D_rebinned_M_pT->Write("hN_1D_rebinned_M_pT");

hZ_1D_rebinned_pT_DCA->Write("hZ_1D_rebinned_pT_DCA");
hB_1D_rebinned_pT_DCA->Write("hB_1D_rebinned_pT_DCA");
hN_1D_rebinned_pT_DCA->Write("hN_1D_rebinned_pT_DCA");
hZ_1D_rebinned_pT_M->Write("hZ_1D_rebinned_pT_M");
hB_1D_rebinned_pT_M->Write("hB_1D_rebinned_pT_M");
hN_1D_rebinned_pT_M->Write("hN_1D_rebinned_pT_M");

// ================================ final method ================================
TH1D* h_lspp_corr_1D_M = (TH1D*)hlspp_corr_3D->Project3D("x");
TH1D* h_lsmm_corr_1D_M = (TH1D*)hlsmm_corr_3D->Project3D("x");
TH1D* h_lspp_corr_1D_pT = (TH1D*)hlspp_corr_3D->Project3D("y");
TH1D* h_lsmm_corr_1D_pT = (TH1D*)hlsmm_corr_3D->Project3D("y");

TH1D* h_ULS_1D_M = (TH1D*)h_rebinned_uls_3D->Project3D("x");
TH1D* h_ULS_1D_pT = (TH1D*)h_rebinned_uls_3D->Project3D("y");

h_lspp_corr_1D_M->Write("h_lspp_corr_1D_M_compare_olaf");
h_lsmm_corr_1D_M->Write("h_lsmm_corr_1D_M_compare_olaf");
h_lspp_corr_1D_pT->Write("h_lspp_corr_1D_pT_compare_olaf");
h_lsmm_corr_1D_pT->Write("h_lsmm_corr_1D_pT_compare_olaf");



TH1D* hZ_1D_rebinned_M = new TH1D("hZ_1D_rebinned_M","Z;M;counts",h_ULS_1D_M->GetNbinsX(),h_ULS_1D_M->GetXaxis()->GetXbins()->GetArray());
TH1D* hB_1D_rebinned_M = new TH1D("hB_1D_rebinned_M","B;M;counts",h_ULS_1D_M->GetNbinsX(),h_ULS_1D_M->GetXaxis()->GetXbins()->GetArray());
TH1D* hN_1D_rebinned_M = new TH1D("hN_1D_rebinned_M","N;M;counts",h_ULS_1D_M->GetNbinsX(),h_ULS_1D_M->GetXaxis()->GetXbins()->GetArray());
TH1D* hZ_1D_rebinned_pT = new TH1D("hZ_1D_rebinned_pT","Z;M;counts",h_ULS_1D_pT->GetNbinsX(),h_ULS_1D_pT->GetXaxis()->GetXbins()->GetArray());
TH1D* hB_1D_rebinned_pT = new TH1D("hB_1D_rebinned_pT","B;M;counts",h_ULS_1D_pT->GetNbinsX(),h_ULS_1D_pT->GetXaxis()->GetXbins()->GetArray());
TH1D* hN_1D_rebinned_pT = new TH1D("hN_1D_rebinned_pT","N;M;counts",h_ULS_1D_pT->GetNbinsX(),h_ULS_1D_pT->GetXaxis()->GetXbins()->GetArray());


for (int ix = 1; ix <= h_lspp_corr_1D_M->GetNbinsX(); ++ix) {


        double Bpp = h_lspp_corr_1D_M->GetBinContent(ix);
        double eBpp = h_lspp_corr_1D_M->GetBinError(ix); 

        double Bmm = h_lsmm_corr_1D_M->GetBinContent(ix);
        double eBmm = h_lsmm_corr_1D_M->GetBinError(ix);

        double ULS = h_ULS_1D_M->GetBinContent(ix);
        double eULS = std::sqrt(ULS);

        double B = 0.0;
        double eB = 0.0; 



        B = 2.0 * std::sqrt(Bpp * Bmm);
    
        eB = std::sqrt((Bmm / Bpp) * eBpp * eBpp + (Bpp / Bmm) * eBmm * eBmm);
        
        double Z  = ULS - B;
        double eZ = std::sqrt(eULS * eULS + eB * eB);
    
        hB_1D_rebinned_M->SetBinContent(ix, B);
        hB_1D_rebinned_M->SetBinError(ix, eB);
    
        hZ_1D_rebinned_M->SetBinContent(ix, Z);
        hZ_1D_rebinned_M->SetBinError(ix, eZ);

        hN_1D_rebinned_M->SetBinContent(ix, ULS);
        hN_1D_rebinned_M->SetBinError(ix, eULS);
    }

for (int ix = 1; ix <= h_lspp_corr_1D_pT->GetNbinsX(); ++ix) {

    double Bpp  = h_lspp_corr_1D_pT->GetBinContent(ix);
    double eBpp = h_lspp_corr_1D_pT->GetBinError(ix);

    double Bmm  = h_lsmm_corr_1D_pT->GetBinContent(ix);
    double eBmm = h_lsmm_corr_1D_pT->GetBinError(ix);

    double ULS  = h_ULS_1D_pT->GetBinContent(ix);
    double eULS = std::sqrt(ULS);

    double B = 0.0;
    double eB = 0.0;

    B = 2.0 * std::sqrt(Bpp * Bmm);

    eB = std::sqrt((Bmm / Bpp) * eBpp * eBpp + (Bpp / Bmm) * eBmm * eBmm);

    double Z  = ULS - B;
    double eZ = std::sqrt(eULS * eULS + eB * eB);

    hB_1D_rebinned_pT->SetBinContent(ix, B);
    hB_1D_rebinned_pT->SetBinError(ix, eB);

    hZ_1D_rebinned_pT->SetBinContent(ix, Z);
    hZ_1D_rebinned_pT->SetBinError(ix, eZ);

    hN_1D_rebinned_pT->SetBinContent(ix, ULS);
    hN_1D_rebinned_pT->SetBinError(ix, eULS);
}
hZ_1D_rebinned_M->Scale(norm);
hB_1D_rebinned_M->Scale(norm);
hN_1D_rebinned_M->Scale(norm);
hZ_1D_rebinned_pT->Scale(norm);
hB_1D_rebinned_pT->Scale(norm);
hN_1D_rebinned_pT->Scale(norm);

hZ_1D_rebinned_M->Write("hZ_1D_M_compare_olaf");
hB_1D_rebinned_M->Write("hB_1D_M_compare_olaf");
hN_1D_rebinned_M->Write("hN_1D_M_compare_olaf");

hZ_1D_rebinned_pT->Write("hZ_1D_pT_compare_olaf");
hB_1D_rebinned_pT->Write("hB_1D_pT_compare_olaf");  
hN_1D_rebinned_pT->Write("hN_1D_pT_compare_olaf");



TH1D *hZ_1D_rebinned_visual_M = new TH1D("hZ_1D_rebinned_visual_M",  "Z;Y;Z", xbins_visual.size()-1, xbins_visual.data());
TH1D *hB_1D_rebinned_visual_M = new TH1D("hB_1D_rebinned_visual_M",  "B;Y;B", xbins_visual.size()-1, xbins_visual.data());
TH1D *hN_1D_rebinned_visual_M = new TH1D("hN_1D_rebinned_visual_M",  "N;Y;N", xbins_visual.size()-1, xbins_visual.data());

TH1D *hZ_1D_rebinned_visual_pT = new TH1D("hZ_1D_rebinned_visual_pT", "Z;Y;Z", ybins_visual.size()-1, ybins_visual.data());
TH1D *hB_1D_rebinned_visual_pT = new TH1D("hB_1D_rebinned_visual_pT", "B;Y;B", ybins_visual.size()-1, ybins_visual.data());
TH1D *hN_1D_rebinned_visual_pT = new TH1D("hN_1D_rebinned_visual_pT", "N;Y;N", ybins_visual.size()-1, ybins_visual.data());


for (size_t i = 0; i < meeSlices2_visual.size(); ++i) {

    double sumZ = 0, sumB = 0, sumN = 0;
    double errZ2 = 0, errB2 = 0, errN2 = 0;


    int xlo = meeSlices2_visual[i].first;
    int xhi = meeSlices2_visual[i].second;

    for (int ix = xlo; ix <= xhi; ++ix) {

        double Z = hZ_1D_rebinned_M->GetBinContent(ix);
        double B = hB_1D_rebinned_M->GetBinContent(ix);
        double N = hN_1D_rebinned_M->GetBinContent(ix);

        double eZ = hZ_1D_rebinned_M->GetBinError(ix);
        double eB = hB_1D_rebinned_M->GetBinError(ix);
        double eN = hN_1D_rebinned_M->GetBinError(ix);

        sumZ += Z;
        sumB += B;
        sumN += N;


        errZ2 += eZ * eZ;
        errB2 += eB * eB;
        errN2 += eN * eN;
    }

    hZ_1D_rebinned_visual_M->SetBinContent(i + 1, sumZ);
    hB_1D_rebinned_visual_M->SetBinContent(i + 1, sumB);
    hN_1D_rebinned_visual_M->SetBinContent(i + 1, sumN);

    hZ_1D_rebinned_visual_M->SetBinError(i + 1, std::sqrt(errZ2));
    hB_1D_rebinned_visual_M->SetBinError(i + 1, std::sqrt(errB2));
    hN_1D_rebinned_visual_M->SetBinError(i + 1, std::sqrt(errN2));
}

for (size_t i = 0; i < ptSlices2_visual.size(); ++i) {

    double sumZ = 0, sumB = 0, sumN = 0;
    double errZ2 = 0, errB2 = 0, errN2 = 0;


    int xlo = ptSlices2_visual[i].first;
    int xhi = ptSlices2_visual[i].second;

    for (int ix = xlo; ix <= xhi; ++ix) {

        double Z_pT = hZ_1D_rebinned_pT->GetBinContent(ix);
        double B_pT = hB_1D_rebinned_pT->GetBinContent(ix);
        double N_pT = hN_1D_rebinned_pT->GetBinContent(ix);

        double eZ = hZ_1D_rebinned_pT->GetBinError(ix);
        double eB = hB_1D_rebinned_pT->GetBinError(ix);
        double eN = hN_1D_rebinned_pT->GetBinError(ix);

        sumZ += Z_pT;
        sumB += B_pT;
        sumN += N_pT;


        errZ2 += eZ * eZ;
        errB2 += eB * eB;
        errN2 += eN * eN;
    }

    hZ_1D_rebinned_visual_pT->SetBinContent(i + 1, sumZ);
    hB_1D_rebinned_visual_pT->SetBinContent(i + 1, sumB);
    hN_1D_rebinned_visual_pT->SetBinContent(i + 1, sumN);

    hZ_1D_rebinned_visual_pT->SetBinError(i + 1, std::sqrt(errZ2));
    hB_1D_rebinned_visual_pT->SetBinError(i + 1, std::sqrt(errB2));
    hN_1D_rebinned_visual_pT->SetBinError(i + 1, std::sqrt(errN2));
}

for (int i = 1; i <= hZ_1D_rebinned_visual_M->GetNbinsX(); ++i) {

    double width = hZ_1D_rebinned_visual_M->GetBinWidth(i);
    if (width <= 0) continue;

    hZ_1D_rebinned_visual_M->SetBinContent(i, hZ_1D_rebinned_visual_M->GetBinContent(i) / width);
    hB_1D_rebinned_visual_M->SetBinContent(i, hB_1D_rebinned_visual_M->GetBinContent(i) / width);
    hN_1D_rebinned_visual_M->SetBinContent(i, hN_1D_rebinned_visual_M->GetBinContent(i) / width);

    hZ_1D_rebinned_visual_M->SetBinError(i, hZ_1D_rebinned_visual_M->GetBinError(i) / width);
    hB_1D_rebinned_visual_M->SetBinError(i, hB_1D_rebinned_visual_M->GetBinError(i) / width);
    hN_1D_rebinned_visual_M->SetBinError(i, hN_1D_rebinned_visual_M->GetBinError(i) / width);
}

for (int i = 1; i <= hZ_1D_rebinned_visual_pT->GetNbinsX(); ++i) {

    double width = hZ_1D_rebinned_visual_pT->GetBinWidth(i);
    if (width <= 0) continue;

    hZ_1D_rebinned_visual_pT->SetBinContent(i, hZ_1D_rebinned_visual_pT->GetBinContent(i) / width);
    hB_1D_rebinned_visual_pT->SetBinContent(i, hB_1D_rebinned_visual_pT->GetBinContent(i) / width);
    hN_1D_rebinned_visual_pT->SetBinContent(i, hN_1D_rebinned_visual_pT->GetBinContent(i) / width);

    hZ_1D_rebinned_visual_pT->SetBinError(i, hZ_1D_rebinned_visual_pT->GetBinError(i) / width);
    hB_1D_rebinned_visual_pT->SetBinError(i, hB_1D_rebinned_visual_pT->GetBinError(i) / width);
    hN_1D_rebinned_visual_pT->SetBinError(i, hN_1D_rebinned_visual_pT->GetBinError(i) / width);
}

hZ_1D_rebinned_visual_M->Write("hZ_1D_rebinned_M_compare_olaf");
hB_1D_rebinned_visual_M->Write("hB_1D_rebinned_M_compare_olaf");
hN_1D_rebinned_visual_M->Write("hN_1D_rebinned_M_compare_olaf");
hZ_1D_rebinned_visual_pT->Write("hZ_1D_rebinned_pT_compare_olaf");
hB_1D_rebinned_visual_pT->Write("hB_1D_rebinned_pT_compare_olaf");
hN_1D_rebinned_visual_pT->Write("hN_1D_rebinned_pT_compare_olaf");

TCanvas* cMcompare = new TCanvas("cMcompare", "", 800, 800);

cMcompare->SetBottomMargin(0.1); 
cMcompare->SetLeftMargin(0.15);
cMcompare->SetTicks(1,1);
cMcompare->SetLogy();
cMcompare->Draw();

// Optional styling (line colors / markers)
hZ_1D_rebinned_visual_M->SetMarkerStyle(24);
hZ_1D_rebinned_visual_M->SetMarkerColor(kBlue);
hZ_1D_rebinned_visual_M->SetLineColor(kBlack);
hB_1D_rebinned_visual_M->SetLineColor(kBlack);
hN_1D_rebinned_visual_M->SetLineColor(kBlack);
hZ_1D_rebinned_visual_M->Draw("E1P");
hZ_1D_rebinned_visual_M->SetStats(0);
hZ_1D_rebinned_visual_M->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
hZ_1D_rebinned_visual_M->GetYaxis()->SetTitle("#frac{1}{N_{ev}} #frac{dN}{dm_{ee}} (GeV/c^{2})^{-1}");
hZ_1D_rebinned_visual_M->SetTitle("");
hZ_1D_rebinned_visual_M->GetYaxis()->SetTitleOffset(1.9);
hZ_1D_rebinned_visual_M->GetXaxis()->SetRangeUser(0, 4);
hZ_1D_rebinned_visual_M->GetYaxis()->SetRangeUser(1e-6, 10);


hN_1D_rebinned_visual_M->SetMarkerStyle(25);
hB_1D_rebinned_visual_M->SetMarkerStyle(26);
hN_1D_rebinned_visual_M->SetMarkerColor(kGreen);
hB_1D_rebinned_visual_M->SetMarkerColor(kRed);

// Overlay others
hN_1D_rebinned_visual_M->Draw("E1P SAME");
hB_1D_rebinned_visual_M->Draw("E1P SAME");

// Legend
TLegend* legendM = new TLegend(0.6, 0.7, 0.79, 0.89);
legendM->AddEntry(hN_1D_rebinned_visual_M, "signal+background", "lp");
legendM->AddEntry(hB_1D_rebinned_visual_M, "background", "lp");
legendM->AddEntry(hZ_1D_rebinned_visual_M, "signal", "lp");
legendM->SetBorderSize(0);
legendM->SetTextSize(0.03);
legendM->Draw();

// Text labels

latex.SetNDC();
latex.SetTextSize(0.026);
latex.SetTextAlign(13);

latex.DrawLatex(0.32, 0.88, "Pb-Pb, #sqrt{s_{NN}} = 5.36 TeV");
latex.DrawLatex(0.32, 0.846, "Centrality: 0-90%");
latex.DrawLatex(0.32, 0.82, "p_{T} > 0.2 GeV/c,|#eta| < 0.8");

// Save
cMcompare->SaveAs(("mee_spectra_olaf_method_0_6_" + cenLabels[ip] + ".png").c_str());

// ======================================= projection of 2D histograms and rebinning =======================================
TH1D* hZ_1D_2DR_pT = hZ_2D_rebinned->ProjectionY("hZ_1D_pT");
TH1D* hB_1D_2DR_pT = hB_2D_rebinned->ProjectionY("hB_1D_pT");
TH1D* hN_1D_2DR_pT = hN_2D_rebinned->ProjectionY("hN_1D_pT");

TH1D* hZ_1D_2DR_mee = hZ_2D_rebinned->ProjectionX("hZ_1D");
TH1D* hB_1D_2DR_mee = hB_2D_rebinned->ProjectionX("hB_1D");
TH1D* hN_1D_2DR_mee = hN_2D_rebinned->ProjectionX("hN_1D");

hZ_1D_2DR_mee->Scale(norm);
hB_1D_2DR_mee->Scale(norm);
hN_1D_2DR_mee->Scale(norm);
hZ_1D_2DR_pT->Scale(norm);
hB_1D_2DR_pT->Scale(norm);
hN_1D_2DR_pT->Scale(norm);

hZ_1D_2DR_mee->Write("hZ_1D_2DR_M_compare_olaf");
hB_1D_2DR_mee->Write("hB_1D_2DR_M_compare_olaf");
hN_1D_2DR_mee->Write("hN_1D_2DR_M_compare_olaf");
hZ_1D_2DR_pT->Write("hZ_1D_2DR_pT_compare_olaf");
hB_1D_2DR_pT->Write("hB_1D_2DR_pT_compare_olaf");
hN_1D_2DR_pT->Write("hN_1D_2DR_pT_compare_olaf");

TH1D *hZ_1D_rebinned_visual_pT_2DR = new TH1D("hZ_1D_rebinned_visual_pT_2DR", "Z;p_{T};Z",ybins_visual.size()-1, ybins_visual.data());
TH1D *hB_1D_rebinned_visual_pT_2DR = new TH1D("hB_1D_rebinned_visual_pT_2DR", "B;p_{T};B",ybins_visual.size()-1, ybins_visual.data());
TH1D *hN_1D_rebinned_visual_pT_2DR = new TH1D("hN_1D_rebinned_visual_pT_2DR", "N;p_{T};N",ybins_visual.size()-1, ybins_visual.data());

TH1D *hZ_1D_rebinned_visual_M_2DR = new TH1D("hZ_1D_rebinned_visual_M_2DR", "Z;M_{ee};Z",xbins_visual.size()-1, xbins_visual.data());
TH1D *hB_1D_rebinned_visual_M_2DR = new TH1D("hB_1D_rebinned_visual_M_2DR", "B;M_{ee};B",xbins_visual.size()-1, xbins_visual.data());
TH1D *hN_1D_rebinned_visual_M_2DR = new TH1D("hN_1D_rebinned_visual_M_2DR", "N;M_{ee};N",xbins_visual.size()-1, xbins_visual.data());

for (size_t i = 0; i < ptSlices2_visual.size(); ++i) {

    double sumZ = 0, sumB = 0, sumN = 0;
    double errZ2 = 0, errB2 = 0, errN2 = 0;

    int ylo = ptSlices2_visual[i].first;
    int yhi = ptSlices2_visual[i].second;

    for (int iy = ylo; iy <= yhi; ++iy) {

        double Z = hZ_1D_2DR_pT->GetBinContent(iy);
        double B = hB_1D_2DR_pT->GetBinContent(iy);
        double N = hN_1D_2DR_pT->GetBinContent(iy);

        double eZ = hZ_1D_2DR_pT->GetBinError(iy);
        double eB = hB_1D_2DR_pT->GetBinError(iy);
        double eN = hN_1D_2DR_pT->GetBinError(iy);

        sumZ += Z;
        sumB += B;
        sumN += N;

        errZ2 += eZ * eZ;
        errB2 += eB * eB;
        errN2 += eN * eN;
    }

    hZ_1D_rebinned_visual_pT_2DR->SetBinContent(i+1, sumZ);
    hB_1D_rebinned_visual_pT_2DR->SetBinContent(i+1, sumB);
    hN_1D_rebinned_visual_pT_2DR->SetBinContent(i+1, sumN);

    hZ_1D_rebinned_visual_pT_2DR->SetBinError(i+1, sqrt(errZ2));
    hB_1D_rebinned_visual_pT_2DR->SetBinError(i+1, sqrt(errB2));
    hN_1D_rebinned_visual_pT_2DR->SetBinError(i+1, sqrt(errN2));
}
for (size_t i = 0; i < meeSlices2_visual.size(); ++i) {

    double sumZ = 0, sumB = 0, sumN = 0;
    double errZ2 = 0, errB2 = 0, errN2 = 0;

    int xlo = meeSlices2_visual[i].first;
    int xhi = meeSlices2_visual[i].second;

    for (int ix = xlo; ix <= xhi; ++ix) {

        double Z = hZ_1D_2DR_mee->GetBinContent(ix);
        double B = hB_1D_2DR_mee->GetBinContent(ix);
        double N = hN_1D_2DR_mee->GetBinContent(ix);

        double eZ = hZ_1D_2DR_mee->GetBinError(ix);
        double eB = hB_1D_2DR_mee->GetBinError(ix);
        double eN = hN_1D_2DR_mee->GetBinError(ix);

        sumZ += Z;
        sumB += B;
        sumN += N;

        errZ2 += eZ * eZ;
        errB2 += eB * eB;
        errN2 += eN * eN;
    }

    hZ_1D_rebinned_visual_M_2DR->SetBinContent(i+1, sumZ);
    hB_1D_rebinned_visual_M_2DR->SetBinContent(i+1, sumB);
    hN_1D_rebinned_visual_M_2DR->SetBinContent(i+1, sumN);

    hZ_1D_rebinned_visual_M_2DR->SetBinError(i+1, sqrt(errZ2));
    hB_1D_rebinned_visual_M_2DR->SetBinError(i+1, sqrt(errB2));
    hN_1D_rebinned_visual_M_2DR->SetBinError(i+1, sqrt(errN2));
}
for (int i = 1; i <= hZ_1D_rebinned_visual_pT_2DR->GetNbinsX(); ++i) {

    double width = hZ_1D_rebinned_visual_pT_2DR->GetBinWidth(i);
    if (width <= 0) continue;

    hZ_1D_rebinned_visual_pT_2DR->SetBinContent(i, hZ_1D_rebinned_visual_pT_2DR->GetBinContent(i) / width);
    hB_1D_rebinned_visual_pT_2DR->SetBinContent(i, hB_1D_rebinned_visual_pT_2DR->GetBinContent(i) / width);
    hN_1D_rebinned_visual_pT_2DR->SetBinContent(i, hN_1D_rebinned_visual_pT_2DR->GetBinContent(i) / width);

    hZ_1D_rebinned_visual_pT_2DR->SetBinError(i, hZ_1D_rebinned_visual_pT_2DR->GetBinError(i) / width);
    hB_1D_rebinned_visual_pT_2DR->SetBinError(i, hB_1D_rebinned_visual_pT_2DR->GetBinError(i) / width);
    hN_1D_rebinned_visual_pT_2DR->SetBinError(i, hN_1D_rebinned_visual_pT_2DR->GetBinError(i) / width);
}
for (int i = 1; i <= hZ_1D_rebinned_visual_M_2DR->GetNbinsX(); ++i) {

    double width = hZ_1D_rebinned_visual_M_2DR->GetBinWidth(i);
    if (width <= 0) continue;

    hZ_1D_rebinned_visual_M_2DR->SetBinContent(i, hZ_1D_rebinned_visual_M_2DR->GetBinContent(i) / width);
    hB_1D_rebinned_visual_M_2DR->SetBinContent(i, hB_1D_rebinned_visual_M_2DR->GetBinContent(i) / width);
    hN_1D_rebinned_visual_M_2DR->SetBinContent(i, hN_1D_rebinned_visual_M_2DR->GetBinContent(i) / width);

    hZ_1D_rebinned_visual_M_2DR->SetBinError(i, hZ_1D_rebinned_visual_M_2DR->GetBinError(i) / width);
    hB_1D_rebinned_visual_M_2DR->SetBinError(i, hB_1D_rebinned_visual_M_2DR->GetBinError(i) / width);
    hN_1D_rebinned_visual_M_2DR->SetBinError(i, hN_1D_rebinned_visual_M_2DR->GetBinError(i) / width);
}
hZ_1D_rebinned_visual_M_2DR->Write("hZ_1D_rebinned_visual_M_2DR");
hB_1D_rebinned_visual_M_2DR->Write("hB_1D_rebinned_visual_M_2DR");
hN_1D_rebinned_visual_M_2DR->Write("hN_1D_rebinned_visual_M_2DR");
hZ_1D_rebinned_visual_pT_2DR->Write("hZ_1D_rebinned_visual_pT_2DR");
hB_1D_rebinned_visual_pT_2DR->Write("hB_1D_rebinned_visual_pT_2DR");
hN_1D_rebinned_visual_pT_2DR->Write("hN_1D_rebinned_visual_pT_2DR");

TH2D *hZ_2D_rebinned_vis  = new TH2D("hZ_2D_rebinned_vis","Z;X;Y",xbins_visual.size()-1,xbins_visual.data(),ybins_visual.size()-1,ybins_visual.data());
TH2D *hB_2D_rebinned_vis  = new TH2D("hB_2D_rebinned_vis","B;X;Y",xbins_visual.size()-1,xbins_visual.data(),ybins_visual.size()-1,ybins_visual.data());
TH2D *hN_2D_rebinned_vis  = new TH2D("hN_2D_rebinned_vis","N;X;Y",xbins_visual.size()-1,xbins_visual.data(),ybins_visual.size()-1,ybins_visual.data());
TH2D *hSoverB_2D_rebinned_vis  = new TH2D("hSoverB_2D_rebinned_vis","N;X;Y",xbins_visual.size()-1,xbins_visual.data(),ybins_visual.size()-1,ybins_visual.data());
TH2D *hSigma_Z_2D_rebinned_vis = new TH2D("hSigma_Z_2D_rebinned_vis","σZ;X;Y",xbins_visual.size()-1,xbins_visual.data(),ybins_visual.size()-1,ybins_visual.data());
TH2D *hSigma_B_2D_rebinned_vis = new TH2D("hSigma_B_2D_rebinned_vis","σB;X;Y",xbins_visual.size()-1,xbins_visual.data(),ybins_visual.size()-1,ybins_visual.data());
TH2D *hSigma_N_2D_rebinned_vis = new TH2D("hSigma_N_2D_rebinned_vis","σN;X;Y",xbins_visual.size()-1,xbins_visual.data(),ybins_visual.size()-1,ybins_visual.data());
TH2D *hSigma_SoverB_2D_rebinned_vis = new TH2D("hSigma_N_2D_rebinned_vis","σN;X;Y",xbins_visual.size()-1,xbins_visual.data(),ybins_visual.size()-1,ybins_visual.data());

for(size_t i=0;i<meeSlices2_visual.size();++i){
    for(size_t j=0;j<ptSlices2_visual.size();++j){

        double sumZ=0;
        double sumB=0;
        double sumN=0;

        double errZ2=0;
        double errB2=0;
        double errN2=0;

        int xlo = meeSlices2_visual[i].first;
        int xhi = meeSlices2_visual[i].second;

        int ylo = ptSlices2_visual[j].first;
        int yhi = ptSlices2_visual[j].second;

        for(int ix=xlo; ix<=xhi; ++ix){
            for(int iy=ylo; iy<=yhi; ++iy){

                double Z = hZ_2D_rebinned->GetBinContent(ix,iy);
                double B = hB_2D_rebinned->GetBinContent(ix,iy);
                double N = hN_2D_rebinned->GetBinContent(ix,iy);

                double eZ = hZ_2D_rebinned->GetBinError(ix,iy);
                double eB = hB_2D_rebinned->GetBinError(ix,iy);
                double eN = hN_2D_rebinned->GetBinError(ix,iy);

                sumZ += Z;
                sumB += B;
                sumN += N;

                errZ2 += eZ*eZ;
                errB2 += eB*eB;
                errN2 += eN*eN;

            }
        }
        hZ_2D_rebinned_vis->SetBinContent(i+1,j+1,sumZ);
        hB_2D_rebinned_vis->SetBinContent(i+1,j+1,sumB);
        hN_2D_rebinned_vis->SetBinContent(i+1,j+1,sumN);

        hSigma_Z_2D_rebinned_vis->SetBinContent(i+1,j+1,std::sqrt(errZ2));
        hSigma_B_2D_rebinned_vis->SetBinContent(i+1,j+1,std::sqrt(errB2));
        hSigma_N_2D_rebinned_vis->SetBinContent(i+1,j+1,std::sqrt(errN2));

    }
}
double I1 = hN_2D_rebinned->Integral();
double I2 = hN_2D_rebinned_vis->Integral();

std::cout << "Integral original (rebinned): " << I1 << std::endl;
std::cout << "Integral visual rebin:        " << I2 << std::endl;
std::cout << "Difference:                   " << (I1 - I2) << std::endl;
hZ_2D_rebinned_vis->Write("hZ_2D_rebinned_vis");

for(int ix=1; ix<=hZ_2D_rebinned_vis->GetNbinsX(); ++ix){
    for(int iy=1; iy<=hZ_2D_rebinned_vis->GetNbinsY(); ++iy){

        double Z = hZ_2D_rebinned_vis->GetBinContent(ix,iy);
        double B = hB_2D_rebinned_vis->GetBinContent(ix,iy);

        double eZ = hSigma_Z_2D_rebinned_vis->GetBinContent(ix,iy);
        double eB = hSigma_B_2D_rebinned_vis->GetBinContent(ix,iy);

        if (B > 0 && Z > 0){
            double SB = Z/B;

            double eSB = SB * std::sqrt(
                (eZ*eZ)/(Z*Z) +
                (eB*eB)/(B*B)
            );

            hSoverB_2D_rebinned_vis->SetBinContent(ix,iy,SB);
            hSigma_SoverB_2D_rebinned_vis->SetBinContent(ix,iy,eSB);
        }

    }
}




 // ========================== projection from 2D to 1D ==========================//
 //                   here we project the 2D histograms to 1D                     //
 // ========================== projection from 2D to 1D ==========================//

int nBins_rebin_X = hZ_2D_rebinned_vis->GetNbinsX();
int nBins_rebin_Y = hZ_2D_rebinned_vis->GetNbinsY();

std::vector<double> x_rebinned(nBins_rebin_X); 
std::vector<double> ex_rebinned(nBins_rebin_X);  
std::vector<double> Z_rebinned(nBins_rebin_X, 0.0);
std::vector<double> eyZ_rebinned(nBins_rebin_X, 0.0);  
std::vector<double> N_rebinned(nBins_rebin_X, 0.0);
std::vector<double> eyN_rebinned(nBins_rebin_X, 0.0);
std::vector<double> B_rebinned(nBins_rebin_X, 0.0);
std::vector<double> eyB_rebinned(nBins_rebin_X, 0.0);
std::vector<double> S_over_B_rebinned(nBins_rebin_X, 0.0);
std::vector<double> eyS_over_B_rebinned(nBins_rebin_X, 0.0);
std::vector<double> Significance_mee_spectra(nBins_rebin_X, 0.0);
std::vector<double> eSignificance_mee_spectra(nBins_rebin_X, 0.0);



for (int ix = 0; ix < nBins_rebin_X; ++ix) {
x_rebinned[ix] = hZ_2D_rebinned_vis->GetXaxis()->GetBinCenter(ix + 1); 
}

for (int ix = 0; ix < nBins_rebin_X; ++ix) {
double low  = hZ_2D_rebinned_vis->GetXaxis()->GetBinLowEdge(ix + 1);
double high = hZ_2D_rebinned_vis->GetXaxis()->GetBinUpEdge(ix + 1);

x_rebinned[ix] = 0.5 * (low + high);   
ex_rebinned[ix]     = 0.5 * (high - low);   
}

for (int ix = 0; ix < nBins_rebin_X; ix++) {
double Z_sum = 0.0, B_sum = 0.0, N_sum = 0.0;
double sigmaZ2 = 0.0, sigmaB2 = 0.0, sigmaN2 = 0.0;

for (int iz = 0; iz < nBins_rebin_Y; iz++) {
    Z_sum += hZ_2D_rebinned_vis->GetBinContent(ix+1,iz+1);


    sigmaZ2 += hSigma_Z_2D_rebinned_vis->GetBinContent(ix+1,iz+1)*hSigma_Z_2D_rebinned_vis->GetBinContent(ix+1,iz+1);

}
double binWidth = hZ_2D_rebinned_vis->GetXaxis()->GetBinWidth(ix+1);
Z_rebinned[ix] = Z_sum/(multiplicity*binWidth);
eyZ_rebinned[ix] = std::sqrt(sigmaZ2)/(multiplicity*binWidth);

}

for (int ix = 0; ix < nBins_rebin_X; ix++) {
double Z_sum = 0.0, B_sum = 0.0, N_sum = 0.0;
double sigmaZ2 = 0.0, sigmaB2 = 0.0, sigmaN2 = 0.0;

for (int iz = 0; iz < nBins_rebin_Y; iz++) {
    Z_sum += hN_2D_rebinned_vis->GetBinContent(ix+1,iz+1);


    sigmaZ2 += hSigma_N_2D_rebinned_vis->GetBinContent(ix+1,iz+1)*hSigma_N_2D_rebinned_vis->GetBinContent(ix+1,iz+1);

}
double binWidth =  hZ_2D_rebinned_vis->GetXaxis()->GetBinWidth(ix+1);

N_rebinned[ix] = Z_sum/(multiplicity*binWidth);
eyN_rebinned[ix] = std::sqrt(sigmaZ2)/(multiplicity*binWidth);
}

for (int ix = 0; ix < nBins_rebin_X; ix++) {
double Z_sum = 0.0;
double sigmaZ2 = 0.0;

for (int iz = 0; iz < nBins_rebin_Y; iz++) {
    Z_sum += hB_2D_rebinned_vis->GetBinContent(ix+1,iz+1);


    sigmaZ2 += hSigma_B_2D_rebinned_vis->GetBinContent(ix+1,iz+1)*hSigma_B_2D_rebinned_vis->GetBinContent(ix+1,iz+1);

}
double binWidth  = hZ_2D_rebinned_vis->GetXaxis()->GetBinWidth(ix+1);

B_rebinned[ix] = Z_sum/(multiplicity*binWidth);
eyB_rebinned[ix] = std::sqrt(sigmaZ2)/(multiplicity*binWidth);
}

for (int ix = 0; ix < nBins_rebin_X; ++ix) {

    double Z_sum = 0.0;
    double B_sum = 0.0;
    double sigmaZ2 = 0.0;
    double sigmaB2 = 0.0;

    for (int iy = 0; iy < nBins_rebin_Y; ++iy) {

        double Z = hZ_2D_rebinned_vis->GetBinContent(ix+1, iy+1);
        double B = hB_2D_rebinned_vis->GetBinContent(ix+1, iy+1);

        double sigmaZ = hSigma_Z_2D_rebinned_vis->GetBinContent(ix+1, iy+1);
        double sigmaB = hSigma_B_2D_rebinned_vis->GetBinContent(ix+1, iy+1);

        Z_sum += Z;
        B_sum += B;

        sigmaZ2 += sigmaZ * sigmaZ;
        sigmaB2 += sigmaB * sigmaB;
    }

    if (B_sum > 0 && Z_sum > 0) {
        S_over_B_rebinned[ix] = Z_sum / B_sum;

        eyS_over_B_rebinned[ix] =
            (Z_sum / B_sum) *
            std::sqrt(
                (sigmaZ2 / (Z_sum * Z_sum)) +
                (sigmaB2 / (B_sum * B_sum))
            );
    }
}


for (int iz = 0; iz < nBins_rebin_X; ++iz) {

    double Z_sum = 0.0;
    double N_sum = 0.0;
    for (int ix = 0; ix < nBins_rebin_Y; ++ix) {

        double Z = hZ_2D_rebinned_vis->GetBinContent(iz+1, ix+1);
        double N = hN_2D_rebinned_vis->GetBinContent(iz+1, ix+1);


        Z_sum += Z;
        N_sum += N;
    }

    if (N_sum > 0 && Z_sum > 0) {
        double S = Z_sum / std::sqrt(N_sum);
        double eS = std::sqrt(
            (Z_sum / N_sum) +
            (Z_sum * Z_sum) / (4.0 * N_sum * N_sum)
        );
        double binWidth  = hZ_2D_rebinned_vis->GetXaxis()->GetBinWidth(iz+1);
        Significance_mee_spectra[iz]  = S/binWidth;
        eSignificance_mee_spectra[iz] = eS/binWidth;

    }
}

 // ========================== projection from 2D to 1D for pT ==========================//
 //                        here we project the 2D histograms to 1D                       //
 // ========================== projection from 2D to 1D for pT ==========================//
 
 std::vector<double> N_rebinned_pt(nBins_rebin_Y, 0.0);
 std::vector<double> eyN_rebinned_pt(nBins_rebin_Y, 0.0);
 std::vector<double> B_rebinned_pt(nBins_rebin_Y, 0.0);
 std::vector<double> eyB_rebinned_pt(nBins_rebin_Y, 0.0);
 std::vector<double> Z_rebinned_pt(nBins_rebin_Y, 0.0);
 std::vector<double> eyZ_rebinned_pt(nBins_rebin_Y, 0.0);
 std::vector<double> y_rebinned(nBins_rebin_X, 0.0);
 std::vector<double> ey_rebinned(nBins_rebin_X, 0.0);
 std::vector<double> S_over_B_rebinned_pt(nBins_rebin_Y, 0.0);
 std::vector<double> eyS_over_B_rebinned_pt(nBins_rebin_Y, 0.0);
 std::vector<double> Significance_pt_spectra(nBins_rebin_Y, 0.0);
 std::vector<double> eSignificance_pt_spectra(nBins_rebin_Y, 0.0);
 
 for (int ix = 0; ix < nBins_rebin_Y; ix++) {
 double Z_sum = 0.0, B_sum = 0.0, N_sum = 0.0;
 double sigmaZ2 = 0.0, sigmaB2 = 0.0, sigmaN2 = 0.0;
 
 for (int iz = 0; iz < nBins_rebin_X; iz++) {
     Z_sum += hN_2D_rebinned_vis->GetBinContent(iz+1,ix+1);
 
 
     sigmaZ2 += hSigma_N_2D_rebinned_vis->GetBinContent(iz+1,ix+1)*hSigma_N_2D_rebinned_vis->GetBinContent(iz+1,ix+1);
 
 }
 double binWidth  = hZ_2D_rebinned_vis->GetYaxis()->GetBinWidth(ix+1);
 N_rebinned_pt[ix] = Z_sum/(multiplicity*binWidth);
 eyN_rebinned_pt[ix] = std::sqrt(sigmaZ2)/(multiplicity*binWidth);
 }
 
 for (int ix = 0; ix < nBins_rebin_Y; ix++) {
 double Z_sum = 0.0;
 double sigmaZ2 = 0.0;
 
 for (int iz = 0; iz < nBins_rebin_X; iz++) {
     Z_sum += hB_2D_rebinned_vis->GetBinContent(iz+1,ix+1);
 
 
     sigmaZ2 += hSigma_B_2D_rebinned_vis->GetBinContent(iz+1,ix+1)*hSigma_B_2D_rebinned_vis->GetBinContent(iz+1,ix+1);
 
 }
 double binWidth  = hZ_2D_rebinned_vis->GetYaxis()->GetBinWidth(ix+1);
 B_rebinned_pt[ix] = Z_sum/(multiplicity*binWidth);
 eyB_rebinned_pt[ix] = std::sqrt(sigmaZ2)/(multiplicity*binWidth);
 }
 
 for (int iy = 0; iy < nBins_rebin_Y; ++iy) {
     double Z_sum = 0.0;
     double sigmaZ2 = 0.0;
 
     for (int ix = 0; ix < nBins_rebin_X; ++ix) {
         Z_sum += hZ_2D_rebinned_vis->GetBinContent(ix+1, iy+1);
         double sigma = hSigma_Z_2D_rebinned_vis->GetBinContent(ix+1, iy+1);
         sigmaZ2 += sigma*sigma;
     }
     double binWidth  = hZ_2D_rebinned_vis->GetYaxis()->GetBinWidth(iy+1);
     Z_rebinned_pt[iy] = Z_sum/(multiplicity*binWidth);
     eyZ_rebinned_pt[iy] = std::sqrt(sigmaZ2)/(multiplicity*binWidth);
 }
 
 
 for (int ix = 0; ix < nBins_rebin_Y; ++ix) {
     double low  = hZ_2D_rebinned_vis->GetYaxis()->GetBinLowEdge(ix + 1);
     double high = hZ_2D_rebinned_vis->GetYaxis()->GetBinUpEdge(ix + 1);
 
     y_rebinned[ix]  = 0.5 * (low + high);
     ey_rebinned[ix] = 0.5 * (high - low); 
 }
 
 
 
 for (int iy = 0; iy < nBins_rebin_Y; ++iy) {
 
     double Z_sum = 0.0;
     double B_sum = 0.0;
     double sigmaZ2 = 0.0;
     double sigmaB2 = 0.0;
 
     for (int ix = 0; ix < nBins_rebin_X; ++ix) {
 
         double Z = hZ_2D_rebinned_vis->GetBinContent(ix+1, iy+1);
         double B = hB_2D_rebinned_vis->GetBinContent(ix+1, iy+1);
 
         double sigmaZ = hSigma_Z_2D_rebinned_vis->GetBinContent(ix+1, iy+1);
         double sigmaB = hSigma_B_2D_rebinned_vis->GetBinContent(ix+1, iy+1);
 
         Z_sum += Z;
         B_sum += B;
 
         sigmaZ2 += sigmaZ * sigmaZ;
         sigmaB2 += sigmaB * sigmaB;
     }
 
     if (B_sum > 0 && Z_sum > 0) {
         S_over_B_rebinned_pt[iy] = Z_sum / B_sum;
 
         eyS_over_B_rebinned_pt[iy] =
             (Z_sum / B_sum) *
             std::sqrt(
                 (sigmaZ2 / (Z_sum * Z_sum)) +
                 (sigmaB2 / (B_sum * B_sum))
             );
     }
 }
 
 for (int ix = 0; ix < nBins_rebin_Y; ++ix) {
 
     double Z_sum = 0.0;
     double N_sum = 0.0;
     for (int iz = 0; iz < nBins_rebin_X; ++iz) {
 
         double Z = hZ_2D_rebinned_vis->GetBinContent(iz+1, ix+1);
         double N = hN_2D_rebinned_vis->GetBinContent(iz+1, ix+1);
 
 
         Z_sum += Z;
         N_sum += N;
     }
 
     if (N_sum > 0 && Z_sum > 0) {
         double S = Z_sum / std::sqrt(N_sum);
         double eS = std::sqrt(
             (Z_sum / N_sum) +
             (Z_sum * Z_sum) / (4.0 * N_sum * N_sum)
         );
         double binWidth  = hZ_2D_rebinned_vis->GetYaxis()->GetBinWidth(ix+1);
         Significance_pt_spectra[ix]  = S/(binWidth);
         eSignificance_pt_spectra[ix] = eS/binWidth;
 
     }
 }

// ================================== Create TH1Ds 1D spectra mee and than pT ===============================//
//           create and fill the newly calculated variables into TH1Ds for spectra mee and then pT           //
// ================================== Create TH1Ds 1D spectra mee and than pT ===============================//
std::vector<std::string> labels = {"Signal+Background", "Background_", "Signal_", "signal_over_background","significance"};

TH1D* hN = new TH1D(Form("hN_%s", cenLabels[ip].c_str()),Form("%s; m_{ee} (GeV/c^{2}); 1/N_{ev} dN/dm_{ee}", labels[0].c_str()),xbins_visual.size()-1, xbins_visual.data());
TH1D* hB = new TH1D(Form("hB_%s", cenLabels[ip].c_str()),Form("%s; m_{ee} (GeV/c^{2}); 1/N_{ev} dN/dm_{ee}", labels[1].c_str()),xbins_visual.size()-1, xbins_visual.data());
TH1D* hZ = new TH1D(Form("hZ_%s", cenLabels[ip].c_str()),Form("%s; m_{ee} (GeV/c^{2}); 1/N_{ev} dN/dm_{ee}", labels[2].c_str()),xbins_visual.size()-1, xbins_visual.data());
TH1D* hS_over_B = new TH1D(Form("hS_over_B_%s", cenLabels[ip].c_str()),Form("%s; m_{ee} (GeV/c^{2}); 1/N_{ev} dN/dm_{ee}", labels[3].c_str()),xbins_visual.size()-1, xbins_visual.data());
TH1D* hsignificance = new TH1D(Form("hSignificance_%s", cenLabels[ip].c_str()),Form("%s; m_{ee} (GeV/c^{2}); 1/N_{ev} dN/dm_{ee}", labels[4].c_str()),xbins_visual.size()-1, xbins_visual.data());

// --- Fill TH1Ds from your re-binned arrays ---
for (int ix = 0; ix < nBins_rebin_X; ++ix) {
    hN->SetBinContent(ix+1, N_rebinned[ix]);
    hN->SetBinError(ix+1, eyN_rebinned[ix]);
    hB->SetBinContent(ix+1, B_rebinned[ix]);
    hB->SetBinError(ix+1, eyB_rebinned[ix]);
    hZ->SetBinContent(ix+1, Z_rebinned[ix]);
    hZ->SetBinError(ix+1, eyZ_rebinned[ix]);
    hS_over_B->SetBinContent(ix+1, S_over_B_rebinned[ix]);
    hS_over_B->SetBinError(ix+1, eyS_over_B_rebinned[ix]);
    hsignificance->SetBinContent(ix+1, Significance_mee_spectra[ix]);
    hsignificance->SetBinError(ix+1, eSignificance_mee_spectra[ix]);
}

for (int ix = 1; ix <= hN->GetNbinsX(); ++ix) {
    if (!std::isfinite(hN->GetBinContent(ix))) hN->SetBinContent(ix, 1e-12);
    if (!std::isfinite(hN->GetBinError(ix))) hN->SetBinError(ix, 1e-12);
    if (!std::isfinite(hB->GetBinContent(ix))) hB->SetBinContent(ix, 1e-12);
    if (!std::isfinite(hB->GetBinError(ix))) hB->SetBinError(ix, 1e-12);
    if (!std::isfinite(hZ->GetBinContent(ix))) hZ->SetBinContent(ix, 1e-12);
    if (!std::isfinite(hZ->GetBinError(ix))) hZ->SetBinError(ix, 1e-12);
}

// --- Write to ROOT file ---
hN->Write(Form("spectra_mee_SplusB_%s", cenLabels[ip].c_str()), TObject::kOverwrite);
hB->Write(Form("spectra_mee_B_%s", cenLabels[ip].c_str()), TObject::kOverwrite);
hZ->Write(Form("spectra_mee_S_%s", cenLabels[ip].c_str()), TObject::kOverwrite);
hS_over_B->Write(Form("spectra_mee_S_over_B_%s", cenLabels[ip].c_str()), TObject::kOverwrite);
hsignificance->Write(Form("spectra_mee_significance_%s", cenLabels[ip].c_str()), TObject::kOverwrite);

// --- Create TH1Ds for pT ---
TH1D* hN_pt = new TH1D(Form("hN_pt_%s", cenLabels[ip].c_str()),Form("%s; p_{T} (GeV/c); 1/N_{ev} dN/dp_{T}", labels[0].c_str()),ybins_visual.size()-1, ybins_visual.data());
TH1D* hB_pt = new TH1D(Form("hB_pt_%s", cenLabels[ip].c_str()),Form("%s; p_{T} (GeV/c); 1/N_{ev} dN/dp_{T}", labels[1].c_str()),ybins_visual.size()-1, ybins_visual.data());
TH1D* hZ_pt = new TH1D(Form("hZ_pt_%s", cenLabels[ip].c_str()),Form("%s; p_{T} (GeV/c); 1/N_{ev} dN/dp_{T}", labels[2].c_str()),ybins_visual.size()-1, ybins_visual.data());
TH1D* hS_over_B_pt = new TH1D(Form("hS_over_B_pt_%s", cenLabels[ip].c_str()),Form("%s; p_{T} (GeV/c); 1/N_{ev} dN/dp_{T}", labels[3].c_str()),ybins_visual.size()-1, ybins_visual.data());
TH1D* hSignificance_pt = new TH1D(Form("%s", cenLabels[ip].c_str()),Form("%s; p_{T} (GeV/c); 1/N_{ev} dN/dp_{T}", labels[4].c_str()),ybins_visual.size()-1, ybins_visual.data());

// --- Fill TH1Ds from your re-binned pT arrays ---
for (int ix = 0; ix < nBins_rebin_Y; ++ix) {
    hN_pt->SetBinContent(ix+1, N_rebinned_pt[ix]);
    hN_pt->SetBinError(ix+1, eyN_rebinned_pt[ix]);
    hB_pt->SetBinContent(ix+1, B_rebinned_pt[ix]);
    hB_pt->SetBinError(ix+1, eyB_rebinned_pt[ix]);
    hZ_pt->SetBinContent(ix+1, Z_rebinned_pt[ix]);
    hZ_pt->SetBinError(ix+1, eyZ_rebinned_pt[ix]);
    hS_over_B_pt->SetBinContent(ix+1, S_over_B_rebinned_pt[ix]);
    hS_over_B_pt->SetBinError(ix+1, eyS_over_B_rebinned_pt[ix]);
    hSignificance_pt->SetBinContent(ix+1, Significance_pt_spectra[ix]);
    hSignificance_pt->SetBinError(ix+1, eSignificance_pt_spectra[ix]);
}

// --- Protect against NaN / Inf ---
for (int ix = 1; ix <= hN_pt->GetNbinsX(); ++ix) {
    if (!std::isfinite(hN_pt->GetBinContent(ix))) hN_pt->SetBinContent(ix, 1e-12);
    if (!std::isfinite(hN_pt->GetBinError(ix)))   hN_pt->SetBinError(ix, 1e-12);
    if (!std::isfinite(hB_pt->GetBinContent(ix))) hB_pt->SetBinContent(ix, 1e-12);
    if (!std::isfinite(hB_pt->GetBinError(ix)))   hB_pt->SetBinError(ix, 1e-12);
    if (!std::isfinite(hZ_pt->GetBinContent(ix))) hZ_pt->SetBinContent(ix, 1e-12);
    if (!std::isfinite(hZ_pt->GetBinError(ix)))   hZ_pt->SetBinError(ix, 1e-12);
}
// --- Write to ROOT file ---
hN_pt->Write(Form("spectra_pt_SplusB_%s", cenLabels[ip].c_str()), TObject::kOverwrite);
hB_pt->Write(Form("spectra_pt_B_%s", cenLabels[ip].c_str()), TObject::kOverwrite);
hZ_pt->Write(Form("spectra_pt_S_%s", cenLabels[ip].c_str()), TObject::kOverwrite);
hS_over_B_pt->Write(Form("spectra_pt_S_over_B_%s", cenLabels[ip].c_str()), TObject::kOverwrite);
hSignificance_pt->Write(Form("spectra_pt_Significance_%s", cenLabels[ip].c_str()), TObject::kOverwrite);

// =================================== create TGraphErrors ==================================//
//      here we create TGraphErrors for the spectra with error bars for mee and then pT      //
// =================================== create TGraphErrors ==================================//
TGraphErrors* gNerr_tot = new TGraphErrors(
    nBins_rebin_X,
    x_rebinned.data(),
    N_rebinned.data(),   // total values for chosen slice
    ex_rebinned.data(),
    eyN_rebinned.data() 
);
gNerr_tot->SetLineColor(kBlack);
gNerr_tot->SetMarkerColor(kGreen);
gNerr_tot->SetMarkerStyle(24);
gNerr_tot->SetTitle("");

TGraphErrors* gBerr_tot = new TGraphErrors(
    nBins_rebin_X,
    x_rebinned.data(),
    B_rebinned.data(),   // B values for chosen slice
    ex_rebinned.data(),
    eyB_rebinned.data() 
);
gBerr_tot->SetLineColor(kBlack);
gBerr_tot->SetMarkerColor(kRed);
gBerr_tot->SetMarkerStyle(25);
gBerr_tot->SetTitle("");

TGraphErrors* gSoverBerr_tot = new TGraphErrors(
    nBins_rebin_X,
    x_rebinned.data(),
    S_over_B_rebinned.data(),   // S/B values
    ex_rebinned.data(),
    eyS_over_B_rebinned.data()
);

gSoverBerr_tot->SetLineColor(kBlack);
gSoverBerr_tot->SetMarkerColor(kBlue);
gSoverBerr_tot->SetMarkerStyle(25);
gSoverBerr_tot->SetMarkerSize(1.5);  
gSoverBerr_tot->SetTitle("");

TGraphErrors* gsignificance_err_mee = new TGraphErrors(
    nBins_rebin_X,
    x_rebinned.data(),
    Significance_mee_spectra.data(),  
    ex_rebinned.data(),
    eSignificance_mee_spectra.data()
);

gsignificance_err_mee->SetLineColor(kBlack);
gsignificance_err_mee->SetMarkerColor(kBlue);
gsignificance_err_mee->SetMarkerStyle(25);
gsignificance_err_mee->SetMarkerSize(1.5);  
gsignificance_err_mee->SetTitle("");

    TGraphErrors* gZerr_rebinned = new TGraphErrors(
        nBins_rebin_X,
        x_rebinned.data(),
        Z_rebinned.data(),
        ex_rebinned.data(),
        eyZ_rebinned.data()          
    );
    gZerr_rebinned->SetLineColor(kBlack);
    gZerr_rebinned->SetMarkerColor(kBlue);
    gZerr_rebinned->SetMarkerStyle(26);
    gZerr_rebinned->SetTitle(""); 


    TGraphErrors* gerrS_over_B_rebinned_pt = new TGraphErrors(
        nBins_rebin_Y,
        y_rebinned.data(),
        S_over_B_rebinned_pt.data(),
        ey_rebinned.data(),
        eyS_over_B_rebinned_pt.data()
    );
    
    gerrS_over_B_rebinned_pt->SetLineColor(kBlack);
    gerrS_over_B_rebinned_pt->SetMarkerColor(kBlue);
    gerrS_over_B_rebinned_pt->SetMarkerStyle(26);
    gerrS_over_B_rebinned_pt->SetTitle("");
    
    
    TGraphErrors* gZerr_rebinned_pt = new TGraphErrors(
        nBins_rebin_Y,
        y_rebinned.data(),
        Z_rebinned_pt.data(),
        ey_rebinned.data(),
        eyZ_rebinned_pt.data()
    );
    
    gZerr_rebinned_pt->SetLineColor(kBlack);
    gZerr_rebinned_pt->SetMarkerColor(kBlue);
    gZerr_rebinned_pt->SetMarkerStyle(26);
    gZerr_rebinned_pt->SetTitle("");
    
    TGraphErrors* gNerr_rebin_pt = new TGraphErrors(
        nBins_rebin_Y,
        y_rebinned.data(),
        N_rebinned_pt.data(),  
        ey_rebinned.data(),
        eyN_rebinned_pt.data() 
    );
    gNerr_rebin_pt->SetLineColor(kBlack);
    gNerr_rebin_pt->SetMarkerColor(kGreen);
    gNerr_rebin_pt->SetMarkerStyle(24);
    gNerr_rebin_pt->SetTitle("");
    
    TGraphErrors* gsignificance_err_pt = new TGraphErrors(
        nBins_rebin_Y,
        y_rebinned.data(),
        Significance_pt_spectra.data(),   
        ey_rebinned.data(),
        eSignificance_pt_spectra.data()
    );
    
    gsignificance_err_pt->SetLineColor(kBlack);
    gsignificance_err_pt->SetMarkerColor(kBlue);
    gsignificance_err_pt->SetMarkerStyle(25);
    gsignificance_err_pt->SetMarkerSize(1.5);  
    gsignificance_err_pt->SetTitle("");
    
    TGraphErrors* gBerr_rebin_pt = new TGraphErrors(
        nBins_rebin_Y,
        y_rebinned.data(),
        B_rebinned_pt.data(), 
        ey_rebinned.data(),
        eyB_rebinned_pt.data() 
    );
    gBerr_rebin_pt->SetLineColor(kBlack);
    gBerr_rebin_pt->SetMarkerColor(kRed);
    gBerr_rebin_pt->SetMarkerStyle(25);
    gBerr_rebin_pt->SetTitle(""); 

// ================================== Plotting the spectra with TGraphErrors ==================================//
//            here we plot the spectra with error bars for mee and then pT using TGraphErrors                  //
// ================================== Plotting the spectra with TGraphErrors ==================================//

TCanvas* cNZrebinned = new TCanvas("cNZsmall", "", 800, 800);

cNZrebinned->SetBottomMargin(0.1); 
cNZrebinned->SetLeftMargin(0.13);
cNZrebinned->SetTicks(1,1);
cNZrebinned->SetLogy();
cNZrebinned->Draw();

gZerr_rebinned->Draw("AP");
gZerr_rebinned->SetStats(0);
gZerr_rebinned->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})"); 
gZerr_rebinned->GetYaxis()->SetTitle("#frac{1}{N_{ev}} #frac{dN}{dm_{ee}} (GeV/c^{2})^{-1}");
gZerr_rebinned->SetTitle("");
gZerr_rebinned->GetYaxis()->SetRangeUser(1e-6, 1e-1);  
gZerr_rebinned->GetXaxis()->SetRangeUser(0,4);

gNerr_tot->Draw("P SAME");
gBerr_tot->Draw("P SAME");


TLegend* legendRZ_rebinned = new TLegend(0.6, 0.7, 0.79, 0.89);
legendRZ_rebinned->AddEntry(gZerr_rebinned, "signal", "lp");
legendRZ_rebinned->AddEntry(gNerr_tot, "signal+background", "lp");
legendRZ_rebinned->AddEntry(gBerr_tot, "Background", "lp");
legendRZ_rebinned->SetBorderSize(0);
legendRZ_rebinned->SetTextSize(0.03);
legendRZ_rebinned->Draw();


latex.SetNDC();
latex.SetTextSize(0.03);
latex.SetTextAlign(13);
latex.DrawLatex(0.15, 0.88, "Pb-Pb #sqrt{s_{NN}} = 5.36 TeV");
latex.DrawLatex(0.15, 0.85, "0_90");

cNZrebinned->SaveAs(("mee_spectra_"+ cenLabels[ip]  + ".png").c_str());

TCanvas* cNZrebinned_S_over_B = new TCanvas("cNZsmall_S_over_B", "", 800, 800);
cNZrebinned_S_over_B->SetBottomMargin(0.1);
cNZrebinned_S_over_B->SetLeftMargin(0.11);
cNZrebinned_S_over_B->SetTicks(1,1);
cNZrebinned_S_over_B->SetLogy();
cNZrebinned_S_over_B->Draw();



gSoverBerr_tot->Draw("AP");
gSoverBerr_tot->SetStats(0);
gSoverBerr_tot->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
gSoverBerr_tot->GetYaxis()->SetTitle("S/B");
gSoverBerr_tot->SetTitle("");
gSoverBerr_tot->GetXaxis()->SetRangeUser(0,4);
gSoverBerr_tot->GetYaxis()->SetRangeUser(1e-3, 10);

latex.SetNDC();
latex.SetTextSize(0.03);
latex.SetTextAlign(13);
latex.DrawLatex(0.15, 0.88, "Pb-Pb #sqrt{s_{NN}} = 5.36 TeV");
latex.DrawLatex(0.15, 0.85, "0_90");

cNZrebinned_S_over_B->SaveAs(("mee_spectra_S_over_B_"+ cenLabels[ip]  + ".png").c_str());

TCanvas* csignificance_mee_spectra = new TCanvas("csignificance_mee_spectra", "", 800, 800);
csignificance_mee_spectra->SetBottomMargin(0.1);
csignificance_mee_spectra->SetLeftMargin(0.11);
csignificance_mee_spectra->SetTicks(1,1);
csignificance_mee_spectra->SetLogy();
csignificance_mee_spectra->Draw();

gsignificance_err_mee->Draw("AP");
gsignificance_err_mee->SetStats(0);
gsignificance_err_mee->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})"); 
gsignificance_err_mee->GetYaxis()->SetTitle("significance S/#sqrt{S+B}");
gsignificance_err_mee->SetTitle("");
gsignificance_err_mee->GetYaxis()->SetRangeUser(1e1, 1e5);  
gsignificance_err_mee->GetXaxis()->SetRangeUser(0,4);

latex.SetNDC();
latex.SetTextSize(0.03);
latex.SetTextAlign(13);
latex.DrawLatex(0.15, 0.88, "Pb-Pb #sqrt{s_{NN}} = 5.36 TeV");
latex.DrawLatex(0.15, 0.85, "0_90");

csignificance_mee_spectra->SaveAs(("mee_spectra_significance_" + cenLabels[ip]  + ".png").c_str());

TCanvas* cNZrebinned_pt = new TCanvas("cNZrebinned_pt", "", 800, 800);

cNZrebinned_pt->SetBottomMargin(0.1);
cNZrebinned_pt->SetLeftMargin(0.13);
cNZrebinned_pt->SetTicks(1,1);
cNZrebinned_pt->SetLogy();
cNZrebinned_pt->Draw();


gZerr_rebinned_pt->Draw("AP");
gZerr_rebinned_pt->SetStats(0);
gZerr_rebinned_pt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
gZerr_rebinned_pt->GetYaxis()->SetTitle("#frac{1}{N_{ev}} #frac{dN}{dp_{T}} (GeV/c)^{-1}");
gZerr_rebinned_pt->SetTitle("");
gZerr_rebinned_pt->GetYaxis()->SetRangeUser(3e-8, 1e-1);  
gZerr_rebinned_pt->GetXaxis()->SetRangeUser(0,10);

gNerr_rebin_pt->Draw("P SAME");
gBerr_rebin_pt->Draw("P SAME");



TLegend* legendRZ_rebinned_pt = new TLegend(0.6, 0.7, 0.79, 0.89);
legendRZ_rebinned_pt->AddEntry(gZerr_rebinned_pt, "signal", "lp");
legendRZ_rebinned_pt->AddEntry(gBerr_rebin_pt, "background", "lp");
legendRZ_rebinned_pt->AddEntry(gNerr_rebin_pt, "background + signal", "lp");
legendRZ_rebinned_pt->SetBorderSize(0);
legendRZ_rebinned_pt->SetTextSize(0.03);
legendRZ_rebinned_pt->Draw();


latex.SetNDC();
latex.SetTextSize(0.03);
latex.SetTextAlign(13);
latex.DrawLatex(0.15, 0.88, "Pb-Pb #sqrt{s_{NN}} = 5.36 TeV");
latex.DrawLatex(0.15, 0.85, "0_90");

cNZrebinned_pt->SaveAs(("pt_spectra_" + cenLabels[ip] + ".png").c_str());



TCanvas* cNZrebinned_pt_S_over_B = new TCanvas("cNZrebinned_pt_S_over_B", "", 800, 800);

cNZrebinned_pt_S_over_B->SetBottomMargin(0.1);
cNZrebinned_pt_S_over_B->SetLeftMargin(0.11);
cNZrebinned_pt_S_over_B->SetTicks(1,1);
cNZrebinned_pt_S_over_B->Draw();
cNZrebinned_pt_S_over_B->SetLogy();


gerrS_over_B_rebinned_pt->Draw("AP");
gerrS_over_B_rebinned_pt->SetStats(0);
gerrS_over_B_rebinned_pt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
gerrS_over_B_rebinned_pt->GetYaxis()->SetTitle("S/B"); 
gerrS_over_B_rebinned_pt->SetTitle("");
gerrS_over_B_rebinned_pt->GetXaxis()->SetRangeUser(0,10);     
gerrS_over_B_rebinned_pt->GetYaxis()->SetRangeUser(6e-2, 11);

latex.SetNDC();
latex.SetTextSize(0.03);
latex.SetTextAlign(13);
latex.DrawLatex(0.15, 0.88, "Pb-Pb #sqrt{s_{NN}} = 5.36 TeV");
latex.DrawLatex(0.15, 0.85, "0_90");

cNZrebinned_pt_S_over_B->SaveAs(("pt_spectra_S_over_B_" + cenLabels[ip] + ".png").c_str());



TCanvas* csignificance_pt_spectra = new TCanvas("csignificance_pt_spectra", "", 800, 800);

csignificance_pt_spectra->SetTicks(1,1);
csignificance_pt_spectra->SetLogy();
gsignificance_err_pt->Draw("AP");
gsignificance_err_pt->SetStats(0);

gsignificance_err_pt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
gsignificance_err_pt->GetYaxis()->SetTitle("Significance S/#sqrt{S+B}");
gsignificance_err_pt->SetTitle("");
gsignificance_err_pt->GetYaxis()->SetRangeUser(1, 5e3);
gsignificance_err_pt->GetXaxis()->SetRangeUser(0.0, 10.0);
latex.SetNDC();
latex.SetTextSize(0.03);
latex.SetTextAlign(13);
latex.DrawLatex(0.15, 0.88, "Pb-Pb #sqrt{s_{NN}} = 5.36 TeV");
latex.DrawLatex(0.15, 0.85, "0_90");

csignificance_pt_spectra->SaveAs(("pt_spectra_significance"+ cenLabels[ip]  + ".png").c_str());



TCanvas* cMcompare2 = new TCanvas("cMcompare", "", 800, 900);

TPad *pad1 = new TPad("pad1","pad1",0,0.30,1,1);
TPad *pad2 = new TPad("pad2","pad2",0,0.00,1,0.30);

pad1->SetBottomMargin(0.0);
pad1->SetLeftMargin(0.15);
pad1->SetLogy();
pad1->SetTicks(1,1);

pad2->SetTopMargin(0.0);
pad2->SetBottomMargin(0.30);
pad2->SetLeftMargin(0.15);
pad2->SetTicks(1,1);

pad1->Draw();
pad2->Draw();

pad1->cd();

// -----------------------------
// Baseline: 1D
// -----------------------------
hZ_1D_mee_visual->SetMarkerStyle(27);
hZ_1D_mee_visual->SetMarkerColor(kMagenta);

// IMPORTANT:
// black lines + black error bars
hZ_1D_mee_visual->SetLineColor(kBlack);
hZ_1D_mee_visual->SetMarkerSize(1.2);

hZ_1D_mee_visual->SetStats(0);

hZ_1D_mee_visual->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
hZ_1D_mee_visual->GetYaxis()->SetTitle("#frac{1}{N_{ev}} #frac{dN}{dm_{ee}} (GeV/c^{2})^{-1}");
hZ_1D_mee_visual->SetTitle("");
hZ_1D_mee_visual->GetYaxis()->SetTitleOffset(1.9);

hZ_1D_mee_visual->GetXaxis()->SetRangeUser(0,4);
hZ_1D_mee_visual->GetYaxis()->SetRangeUser(1.1e-6,1e-1);

// draw axis first
hZ_1D_mee_visual->Draw("E1");

// draw markers separately (keeps error bars black)
hZ_1D_mee_visual->Draw("P SAME");


// ----- Olaf
hZ_1D_rebinned_visual_M->SetMarkerStyle(24);
hZ_1D_rebinned_visual_M->SetMarkerColor(kBlue);
hZ_1D_rebinned_visual_M->SetLineColor(kBlack);

// ----- pT slice
hZ_1D_rebinned_M_pT->SetMarkerStyle(25);
hZ_1D_rebinned_M_pT->SetMarkerColor(kGreen);
hZ_1D_rebinned_M_pT->SetLineColor(kBlack);

// ----- DCA slice
hZ_1D_rebinned_M_DCA->SetMarkerStyle(26);
hZ_1D_rebinned_M_DCA->SetMarkerColor(kRed);
hZ_1D_rebinned_M_DCA->SetLineColor(kBlack);

// ----- 2D slice
hZ_1D_rebinned_visual_M_2DR->SetMarkerStyle(28);
hZ_1D_rebinned_visual_M_2DR->SetMarkerColor(kCyan);
hZ_1D_rebinned_visual_M_2DR->SetLineColor(kBlack);


hZ_1D_rebinned_visual_M->Draw("E1 SAME");
hZ_1D_rebinned_visual_M->Draw("P SAME");

hZ_1D_rebinned_M_pT->Draw("E1 SAME");
hZ_1D_rebinned_M_pT->Draw("P SAME");

hZ_1D_rebinned_M_DCA->Draw("E1 SAME");
hZ_1D_rebinned_M_DCA->Draw("P SAME");

hZ_1D_rebinned_visual_M_2DR->Draw("E1 SAME");
hZ_1D_rebinned_visual_M_2DR->Draw("P SAME");

TLegend* legendM2 = new TLegend(0.5, 0.7, 0.79, 0.89);

legendM2->AddEntry(hZ_1D_mee_visual, "1D", "lp");
legendM2->AddEntry(hZ_1D_rebinned_visual_M_2DR, "2D correction", "lp");
legendM2->AddEntry(hZ_1D_rebinned_visual_M, "1D with 3D Rfactor", "lp");
legendM2->AddEntry(hZ_1D_rebinned_M_pT, "1D with 3D Rfactor via pT", "lp");
legendM2->AddEntry(hZ_1D_rebinned_M_DCA, "1D with 3D Rfactor via DCA", "lp");

legendM2->SetBorderSize(0);
legendM2->SetTextSize(0.03);
legendM2->Draw();

latex.SetNDC();
latex.SetTextSize(0.026);
latex.SetTextAlign(13);

latex.DrawLatex(0.20, 0.88, "Pb-Pb, #sqrt{s_{NN}} = 5.36 TeV");
latex.DrawLatex(0.20, 0.846, "0-90%");
latex.DrawLatex(0.20, 0.82, "p_{T} > 0.2 GeV/c, |#eta| < 0.8");

pad2->cd();

TH1D* ratio_olaf =
    (TH1D*)hZ_1D_rebinned_visual_M->Clone("ratio_olaf");

TH1D* ratio_pt =
    (TH1D*)hZ_1D_rebinned_M_pT->Clone("ratio_pt");

TH1D* ratio_dca =
    (TH1D*)hZ_1D_rebinned_M_DCA->Clone("ratio_dca");

TH1D* ratio_2dr =
    (TH1D*)hZ_1D_rebinned_visual_M_2DR->Clone("ratio_2dr");


// divide by baseline (1D)
ratio_olaf->Divide(hZ_1D_mee_visual);
ratio_pt->Divide(hZ_1D_mee_visual);
ratio_dca->Divide(hZ_1D_mee_visual);
ratio_2dr->Divide(hZ_1D_mee_visual);

ratio_olaf->SetStats(0);
ratio_olaf->SetTitle("");

ratio_olaf->GetYaxis()->SetTitle("Ratio");
ratio_olaf->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");

ratio_olaf->GetYaxis()->SetNdivisions(505);

ratio_olaf->GetYaxis()->SetTitleSize(0.08);
ratio_olaf->GetYaxis()->SetLabelSize(0.07);

ratio_olaf->GetXaxis()->SetTitleSize(0.10);
ratio_olaf->GetXaxis()->SetLabelSize(0.08);

ratio_olaf->GetYaxis()->SetTitleOffset(0.7);

ratio_olaf->GetYaxis()->SetRangeUser(0.75,1.5);


// keep black errors
ratio_olaf->SetLineColor(kBlack);
ratio_pt->SetLineColor(kBlack);
ratio_dca->SetLineColor(kBlack);
ratio_2dr->SetLineColor(kBlack);

// marker colors
ratio_olaf->SetMarkerColor(kBlue);
ratio_pt->SetMarkerColor(kGreen);
ratio_dca->SetMarkerColor(kRed);
ratio_2dr->SetMarkerColor(kCyan);

// marker styles
ratio_olaf->SetMarkerStyle(24);
ratio_pt->SetMarkerStyle(25);
ratio_dca->SetMarkerStyle(26);
ratio_2dr->SetMarkerStyle(28);


// ======================================================
// Draw ratios
// ======================================================

ratio_olaf->Draw("E1");
ratio_olaf->Draw("P SAME");

ratio_pt->Draw("E1 SAME");
ratio_pt->Draw("P SAME");

ratio_dca->Draw("E1 SAME");
ratio_dca->Draw("P SAME");

ratio_2dr->Draw("E1 SAME");
ratio_2dr->Draw("P SAME");

double xminmee = ratio_olaf->GetXaxis()->GetXmin();
double xmaxmee = ratio_olaf->GetXaxis()->GetXmax();
// unity line
TLine *line_mee = new TLine(xminmee, 1.0, xmaxmee, 1.0);
line_mee->SetLineStyle(2);   // dashed
line_mee->SetLineWidth(2);
line_mee->SetLineColor(kBlack);
line_mee->Draw("SAME");


pad1->cd();
pad1->Draw();
pad1->Modified();
pad1->Update();


// ======================================================
// Save
// ======================================================

cMcompare2->SaveAs(("mee_spectra_olaf_method_ratio_" + cenLabels[ip] + ".png").c_str());

// ======================================================
// pT comparison canvas
// ======================================================

TCanvas* cPtcompare = new TCanvas("cPtcompare", "", 800, 900);

TPad *pad1_pt = new TPad("pad1_pt","pad1_pt",0,0.30,1,1);
TPad *pad2_pt = new TPad("pad2_pt","pad2_pt",0,0.00,1,0.30);

pad1_pt->SetBottomMargin(0.0);
pad1_pt->SetLeftMargin(0.15);
pad1_pt->SetLogy();
pad1_pt->SetTicks(1,1);

pad2_pt->SetTopMargin(0.0);
pad2_pt->SetBottomMargin(0.30);
pad2_pt->SetLeftMargin(0.15);
pad2_pt->SetTicks(1,1);

pad1_pt->Draw();
pad2_pt->Draw();

pad1_pt->cd();

// ======================================================
// Baseline: 1D
// ======================================================

hZ_1D_pT_visual->SetMarkerStyle(27);
hZ_1D_pT_visual->SetMarkerColor(kMagenta);

hZ_1D_pT_visual->SetLineColor(kBlack);
hZ_1D_pT_visual->SetMarkerSize(1.2);

hZ_1D_pT_visual->SetStats(0);

hZ_1D_pT_visual->GetXaxis()->SetTitle("p_{T,ee} (GeV/c)");
hZ_1D_pT_visual->GetYaxis()->SetTitle("#frac{1}{N_{ev}} #frac{dN}{dp_{T,ee}} (GeV/c)^{-1}");

hZ_1D_pT_visual->SetTitle("");
hZ_1D_pT_visual->GetYaxis()->SetTitleOffset(1.9);

hZ_1D_pT_visual->GetXaxis()->SetRangeUser(0,10);
hZ_1D_pT_visual->GetYaxis()->SetRangeUser(2.0e-7,1e-1);

// draw axis first
hZ_1D_pT_visual->Draw("E1");

// markers separately
hZ_1D_pT_visual->Draw("P SAME");

// ======================================================
// Other methods
// ======================================================

// ----- Olaf
hZ_1D_rebinned_visual_pT->SetMarkerStyle(24);
hZ_1D_rebinned_visual_pT->SetMarkerColor(kBlue);
hZ_1D_rebinned_visual_pT->SetLineColor(kBlack);

// ----- pT slice
hZ_1D_rebinned_pT_M->SetMarkerStyle(25);
hZ_1D_rebinned_pT_M->SetMarkerColor(kGreen);
hZ_1D_rebinned_pT_M->SetLineColor(kBlack);

// ----- DCA slice
hZ_1D_rebinned_pT_DCA->SetMarkerStyle(26);
hZ_1D_rebinned_pT_DCA->SetMarkerColor(kRed);
hZ_1D_rebinned_pT_DCA->SetLineColor(kBlack);

// ----- 2D slice
hZ_1D_rebinned_visual_pT_2DR->SetMarkerStyle(28);
hZ_1D_rebinned_visual_pT_2DR->SetMarkerColor(kCyan);
hZ_1D_rebinned_visual_pT_2DR->SetLineColor(kBlack);


// ======================================================
// Draw spectra
// ======================================================

hZ_1D_rebinned_visual_pT->Draw("E1 SAME");
hZ_1D_rebinned_visual_pT->Draw("P SAME");

hZ_1D_rebinned_pT_M->Draw("E1 SAME");
hZ_1D_rebinned_pT_M->Draw("P SAME");

hZ_1D_rebinned_pT_DCA->Draw("E1 SAME");
hZ_1D_rebinned_pT_DCA->Draw("P SAME");

hZ_1D_rebinned_visual_pT_2DR->Draw("E1 SAME");
hZ_1D_rebinned_visual_pT_2DR->Draw("P SAME");


// ======================================================
// Legend
// ======================================================

TLegend* legendPt = new TLegend(0.5, 0.7, 0.79, 0.89);

legendPt->AddEntry(hZ_1D_pT_visual, "1D", "lp");
legendPt->AddEntry(hZ_1D_rebinned_visual_pT_2DR, "2D correction", "lp");
legendPt->AddEntry(hZ_1D_rebinned_visual_pT, "1D with 3D Rfactor", "lp");
legendPt->AddEntry(hZ_1D_rebinned_pT_M, "1D with 3D Rfactor via M", "lp");
legendPt->AddEntry(hZ_1D_rebinned_pT_DCA, "1D with 3D Rfactor via DCA", "lp");

legendPt->SetBorderSize(0);
legendPt->SetTextSize(0.03);
legendPt->Draw();


// ======================================================
// Labels
// ======================================================

latex.SetNDC();
latex.SetTextSize(0.026);
latex.SetTextAlign(13);

latex.DrawLatex(0.20, 0.88, "Pb-Pb, #sqrt{s_{NN}} = 5.36 TeV");
latex.DrawLatex(0.20, 0.846, "0-90%");
latex.DrawLatex(0.20, 0.82, "p_{T} > 0.2 GeV/c,|#eta| < 0.8");


// ======================================================
// Ratio pad
// ======================================================

pad2_pt->cd();

TH1D* ratio_olaf_pt =
    (TH1D*)hZ_1D_rebinned_visual_pT->Clone("ratio_olaf_pt");

TH1D* ratio_pt_pt =
    (TH1D*)hZ_1D_rebinned_pT_M->Clone("ratio_pt_pt");

TH1D* ratio_dca_pt =
    (TH1D*)hZ_1D_rebinned_pT_DCA->Clone("ratio_dca_pt");

TH1D* ratio_2dr_pt =
    (TH1D*)hZ_1D_rebinned_visual_pT_2DR->Clone("ratio_2dr_pt");


// divide by baseline
ratio_olaf_pt->Divide(hZ_1D_pT_visual);
ratio_pt_pt->Divide(hZ_1D_pT_visual);
ratio_dca_pt->Divide(hZ_1D_pT_visual);
ratio_2dr_pt->Divide(hZ_1D_pT_visual);


// ======================================================
// Ratio styling
// ======================================================

ratio_olaf_pt->SetStats(0);
ratio_olaf_pt->SetTitle("");

ratio_olaf_pt->GetYaxis()->SetTitle("Ratio");
ratio_olaf_pt->GetXaxis()->SetTitle("p_{T,ee} (GeV/c)");

ratio_olaf_pt->GetYaxis()->SetNdivisions(505);

ratio_olaf_pt->GetYaxis()->SetTitleSize(0.08);
ratio_olaf_pt->GetYaxis()->SetLabelSize(0.07);

ratio_olaf_pt->GetXaxis()->SetTitleSize(0.10);
ratio_olaf_pt->GetXaxis()->SetLabelSize(0.08);

ratio_olaf_pt->GetYaxis()->SetTitleOffset(0.7);

ratio_olaf_pt->GetYaxis()->SetRangeUser(0.4,1.8);


// black error bars
ratio_olaf_pt->SetLineColor(kBlack);
ratio_pt_pt->SetLineColor(kBlack);
ratio_dca_pt->SetLineColor(kBlack);
ratio_2dr_pt->SetLineColor(kBlack);

// marker colors
ratio_olaf_pt->SetMarkerColor(kBlue);
ratio_pt_pt->SetMarkerColor(kGreen);
ratio_dca_pt->SetMarkerColor(kRed);
ratio_2dr_pt->SetMarkerColor(kCyan);

// marker styles
ratio_olaf_pt->SetMarkerStyle(24);
ratio_pt_pt->SetMarkerStyle(25);
ratio_dca_pt->SetMarkerStyle(26);
ratio_2dr_pt->SetMarkerStyle(28);


// ======================================================
// Draw ratios
// ======================================================

ratio_olaf_pt->Draw("E1");
ratio_olaf_pt->Draw("P SAME");

ratio_pt_pt->Draw("E1 SAME");
ratio_pt_pt->Draw("P SAME");

ratio_dca_pt->Draw("E1 SAME");
ratio_dca_pt->Draw("P SAME");

ratio_2dr_pt->Draw("E1 SAME");
ratio_2dr_pt->Draw("P SAME");


// unity line
double xmin = ratio_olaf_pt->GetXaxis()->GetXmin();
double xmax = ratio_olaf_pt->GetXaxis()->GetXmax();

TLine *line_pt = new TLine(xmin, 1.0, xmax, 1.0);
line_pt->SetLineStyle(2);   // dashed
line_pt->SetLineWidth(2);
line_pt->SetLineColor(kBlack);
line_pt->Draw("SAME");


// ======================================================
// Final update
// ======================================================

pad1_pt->cd();
pad1_pt->Draw();
pad1_pt->Modified();
pad1_pt->Update();


// ======================================================
// Save
// ======================================================

cPtcompare->SaveAs(("pt_spectra_olaf_method_ratio_" + cenLabels[ip] + ".png").c_str());

// ==================================== reading in qvectors =============================== //
//  reading in the Qvectors for calculating the resolution factor Rn for the flow analysis  //
// ==================================== reading in qvectors =============================== //

int nBinsX = h_proj_uls3D->GetNbinsX();
int nBinsY = h_proj_uls3D->GetNbinsY();
int nBinsZ = h_proj_uls3D->GetNbinsZ();


TH2D* hCent = (TH2D*) d_after->Get("hCentFT0CvsMultNTracksPV");
int NBinsX = hCent->GetNbinsX();
hCent->GetXaxis()->SetRangeUser(0,90);
hCent->GetYaxis()->SetRangeUser(0,6000);

int binMin = hCent->GetXaxis()->FindBin(0.0);
int binMax = hCent->GetXaxis()->FindBin(90.0);
std::vector<TH1D*> slices;
std::vector<double> means_squared(NBinsX, 0.0);
std::vector<double> sigmas(NBinsX, 0.0);

for (int i = binMin; i <= binMax; i++) {
    
    TH1D* h1 = hCent->ProjectionY(Form("slice_xbin_%d", i), i, i);
    if (h1->GetEntries() < 10) continue;   // skip empty slices

    h1->Fit("gaus","Q");
    TF1* fit = h1->GetFunction("gaus");
    if (!fit) continue;


    double mean  = fit->GetParameter(1);
    double sigma = fit->GetParameter(2);

    means_squared[i-1] = mean;
    sigmas[i-1] = sigma;

    //h1->Write();
    slices.push_back(h1);
}


//for (size_t i = 0; i < means_squared.size(); ++i) {
  //  std::cout << "mean[" << i << "] = " << means_squared[i] << std::endl;
//}




TH2D* QTPC_pos_QFT0M = (TH2D*) d_after->Get("hPrfQ2FT0MQ2BPos_CentFT0C");
TH2D* QTPC_neg_QFT0M = (TH2D*) d_after->Get("hPrfQ2FT0MQ2BNeg_CentFT0C"); 
TH2D* QTPC_pos_QTPC_neg = (TH2D*) d_after->Get("hPrfQ2BPosQ2BNeg_CentFT0C"); 

int xbinsq = QTPC_pos_QFT0M->GetNbinsX();
std::vector<double> Q2TPC_pos_Q2TPC_neg2(xbinsq,  0.0);
std::vector<double> Q2TPC_pos_Q2FT0M2(xbinsq,  0.0);
std::vector<double> Q2TPC_neg_Q2FT0M2(xbinsq,  0.0);
std::vector<double> Rn_values(xbinsq, 0.0);
std::vector<double> Rn_errors(xbinsq, 0.0);

for (int ix = 1; ix <= xbinsq; ix++){
    double Q2TPC_pos_Q2TPC_neg = QTPC_pos_QTPC_neg->GetBinContent(ix);
    double Q2TPC_pos_Q2FT0M = QTPC_pos_QFT0M->GetBinContent(ix);
    double Q2TPC_neg_Q2FT0M = QTPC_neg_QFT0M->GetBinContent(ix);
    double Q2TPC_pos_Q2TPC_neg_error = QTPC_pos_QTPC_neg->GetBinError(ix);
    double Q2TPC_pos_Q2FT0M_error = QTPC_pos_QFT0M->GetBinError(ix);
    double Q2TPC_neg_Q2FT0M_error = QTPC_neg_QFT0M->GetBinError(ix);
    if (Q2TPC_pos_Q2TPC_neg > 0) {
    double Rn = std::sqrt ((Q2TPC_pos_Q2FT0M * Q2TPC_neg_Q2FT0M) / ( Q2TPC_pos_Q2TPC_neg));
    Q2TPC_pos_Q2TPC_neg2[ix-1] = Q2TPC_pos_Q2TPC_neg;
    Q2TPC_pos_Q2FT0M2[ix-1] = Q2TPC_pos_Q2FT0M;
    Q2TPC_neg_Q2FT0M2[ix-1] = Q2TPC_neg_Q2FT0M;
    Rn_values[ix-1] = Rn;
    double relErrSq =
    (Q2TPC_pos_Q2FT0M_error * Q2TPC_pos_Q2FT0M_error) / (Q2TPC_pos_Q2FT0M * Q2TPC_pos_Q2FT0M) +
    (Q2TPC_neg_Q2FT0M_error * Q2TPC_neg_Q2FT0M_error) / (Q2TPC_neg_Q2FT0M * Q2TPC_neg_Q2FT0M) +
    (Q2TPC_pos_Q2TPC_neg_error * Q2TPC_pos_Q2TPC_neg_error) / (Q2TPC_pos_Q2TPC_neg * Q2TPC_pos_Q2TPC_neg);

    double Rn_err = 0.5 * Rn * std::sqrt(relErrSq);
    Rn_errors[ix-1] = Rn_err;
    }
}   

double sumRn = 0.0;
int nValid = 0;


for (int ix = 0; ix < xbinsq; ++ix) {
   if (Rn_values[ix] > 0) {
       sumRn += Rn_values[ix];  // weight by the mean of the multiplicity distribution for that centrality bin
       nValid++;
   }
}


double rn_average = (nValid > 0) ? sumRn / nValid : 0.0;

double weighted_sum = 0.0;
double sum_weights = 0.0;

for (int ix = 0; ix < xbinsq; ++ix) {
    if (Rn_values[ix] > 0 && means_squared[ix] > 0) {
        std::cout << "Centrality bin " << ix << ": Rn = " << Rn_values[ix] << ", mean^2 = " << means_squared[ix] << std::endl;
        weighted_sum += Rn_values[ix] * means_squared[ix];
        sum_weights += means_squared[ix];
    }
}

double Rn_average = (sum_weights > 0) ? weighted_sum / sum_weights : 0.0;
std::cout << "Average Rn (unweighted): " << rn_average << std::endl;
std::cout << "Average Rn (weighted): " << Rn_average << std::endl;
double Rn = Rn_average;
double Rval = Rn_average; 
double Rn_variance = 0.0;
for (size_t ix = 0; ix < means_squared.size(); ++ix) {
    if (Rn_values[ix] > 0 && sigmas[ix] > 0) {
        double weight = means_squared[ix];
        double deltaRn = Rn_values[ix] * (sigmas[ix] / sqrt(means_squared[ix])); // relative fluctuation
        Rn_variance += weight * weight * deltaRn * deltaRn;
    }
}
double Rn_error_sigma = (sum_weights > 0) ? sqrt(Rn_variance) / sum_weights : 0.0;
std::cout << "Estimated uncertainty on Rn: " << Rn_error_sigma << std::endl;

TH1D* hRn = new TH1D("hRn"," ;Centrality;R_{n}",xbinsq,QTPC_pos_QTPC_neg->GetXaxis()->GetXmin(),QTPC_pos_QTPC_neg->GetXaxis()->GetXmax());
for (int ix = 0; ix < xbinsq; ix++) {

    hRn->SetBinContent(ix+1, Rn_values[ix]);
    hRn->SetBinError(ix+1, Rn_errors[ix]);

}
hRn->Write();

// ================================== flow part of the analysis remove over and undeflow ==================================//
//         remove the underflow and overflow bins from x y and z axes for the flow part of the analysis                    //
// ================================== flow part of the analysis remove over and undeflow ==================================//

std::vector<BinRange> zbins_flow2 = {
    {-5.0, 5.0, 0.05}   
};
std::vector<double> zbins_flow;
for (auto &[start, end, step] : zbins_flow2) {
    if (zbins_flow.empty()) zbins_flow.push_back(start);
    for (double z = start + step; z < end - 1e-8; z += step) {
        zbins_flow.push_back(z);
    }
    zbins.push_back(end);
}
TAxis* zax2 = h_proj_uls3D->GetZaxis();

auto FindBinLow = [&](TAxis* ax, double edge) {
    return ax->FindBin(edge + delta);
};

auto FindBinHigh = [&](TAxis* ax, double edge) {
    return ax->FindBin(edge - delta);
};

// ---- z bins (keep fine binning) ----

TH3D* h_rebinned_uls_flow = new TH3D("h_uls_flow", "ULS (flow bins);m_{ee};p_{T};#Delta#phi",xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data(),zbins_flow.size()-1, zbins_flow.data());
TH3D* h_rebinned_lspp_flow = new TH3D("h_lspp_flow", "LS++ (flow bins);m_{ee};p_{T};#Delta#phi",xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data(),zbins_flow.size()-1, zbins_flow.data());
TH3D* h_rebinned_lsmm_flow = new TH3D("h_lsmm_flow", "LS-- (flow bins);m_{ee};p_{T};#Delta#phi",xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data(),zbins_flow.size()-1, zbins_flow.data());


std::vector<std::pair<int,int>> qvecSlicesFlow;

for (size_t i = 0; i < zbins_flow.size()-1; ++i) {

    int lo = FindBinLow(zax2, zbins_flow[i]);
    int hi = FindBinHigh(zax2, zbins_flow[i+1]);

    lo = std::max(lo, 1);
    hi = std::min(hi, zax2->GetNbins());

    qvecSlicesFlow.emplace_back(lo, hi);
}



for (size_t i = 0; i < meeSlices2.size(); ++i) {
    for (size_t j = 0; j < ptSlices2.size(); ++j) {
        for (size_t k = 0; k < qvecSlicesFlow.size(); ++k) {

            double sumULS  = 0.0;
            double sumLSPP = 0.0;
            double sumLSMM = 0.0;

            int xlo = meeSlices2[i].first;
            int xhi = meeSlices2[i].second;

            int ylo = ptSlices2[j].first;
            int yhi = ptSlices2[j].second;

            int zlo = qvecSlicesFlow[k].first;
            int zhi = qvecSlicesFlow[k].second;

            for (int ix = xlo; ix <= xhi; ++ix) {
                for (int iy = ylo; iy <= yhi; ++iy) {
                    for (int iz = zlo; iz <= zhi; ++iz) {

                        sumULS  += h_proj_uls3D ->GetBinContent(ix,iy,iz);
                        sumLSPP += h_proj_lspp3D->GetBinContent(ix,iy,iz);
                        sumLSMM += h_proj_lsmm3D->GetBinContent(ix,iy,iz);
                    }
                }
            }

            // ROOT bins start at 1
            h_rebinned_uls_flow ->SetBinContent(i+1, j+1, k+1, sumULS);
            h_rebinned_lspp_flow->SetBinContent(i+1, j+1, k+1, sumLSPP);
            h_rebinned_lsmm_flow->SetBinContent(i+1, j+1, k+1, sumLSMM);
        }
    }
}
// ========================= create 2D histograms for N, Z, B, R, S/B and significance in flow bins =========================//
//             here we calculate the required variables for flow in 2D and save them in TH2D for later                       //
// ========================= create 2D histograms for N, Z, B, R, S/B and significance in flow bins =========================//

TH2D* hN_2D_flow = new TH2D("hN_2D_flow", "N (flow bins);m_{ee};p_{T}",xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data());
TH2D* hZ_2D_flow = new TH2D("hZ_2D_flow", "Z (flow bins);m_{ee};p_{T}",xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data());
TH2D* hB_2D_flow = new TH2D("hB_2D_flow", "B (flow bins);m_{ee};p_{T}",xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data());
TH2D* hR_2D_flow = new TH2D("hR_2D_flow", "R (flow bins);m_{ee};p_{T}",xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data());
TH2D* hS_over_B_2D_flow = new TH2D("hS_over_B_2D_flow", "R (flow bins);m_{ee};p_{T}",xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data());
TH2D* hSignificance_2D_flow = new TH2D("hSignificance_2D_flow","S/#sqrt{S+B};m_{ee};p_{T}",xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data());

TH2D* hSigma_N_2D_flow = new TH2D("hSigma_N_2D_flow", "σN (flow bins);m_{ee};p_{T}",xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data());
TH2D* hSigma_Z_2D_flow = new TH2D("hSigma_Z_2D_flow", "σZ (flow bins);m_{ee};p_{T}",xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data());
TH2D* hSigma_B_2D_flow = new TH2D("hSigma_B_2D_flow", "σB (flow bins);m_{ee};p_{T}",xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data());
TH2D* hSigma_R_2D_flow = new TH2D("hSigma_R_2D_flow", "σR (flow bins);m_{ee};p_{T}",xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data());
TH2D* hSigma_S_over_B_2D_flow = new TH2D("hSigma_S_over_B_2D_flow", "σR (flow bins);m_{ee};p_{T}",xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data());

for (int ix = 1; ix <= hN_2D_rebinned->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= hN_2D_rebinned->GetNbinsY(); ++iy) {

 

        double x = hN_2D_rebinned->GetXaxis()->GetBinCenter(ix);
        double y = hN_2D_rebinned->GetYaxis()->GetBinCenter(iy);

        int ix_new = hN_2D_flow->GetXaxis()->FindBin(x);
        int iy_new = hN_2D_flow->GetYaxis()->FindBin(y);

        // --- N, Z, B: sum ---
        double N = hN_2D_rebinned->GetBinContent(ix, iy);
        double Z = hZ_2D_rebinned->GetBinContent(ix, iy);
        double B = hB_2D_rebinned->GetBinContent(ix, iy);
        double R = hR_2D_rebinned->GetBinContent(ix, iy); 

        double sN = hN_2D_rebinned->GetBinError(ix, iy);
        double sZ = hZ_2D_rebinned->GetBinError(ix, iy);
        double sB = hB_2D_rebinned->GetBinError(ix, iy);
        double sR = hR_2D_rebinned->GetBinError(ix, iy);

        hN_2D_flow->AddBinContent(hN_2D_flow->GetBin(ix_new, iy_new), N);
        hZ_2D_flow->AddBinContent(hZ_2D_flow->GetBin(ix_new, iy_new), Z);
        hB_2D_flow->AddBinContent(hB_2D_flow->GetBin(ix_new, iy_new), B);

        int bin_new = hN_2D_flow->GetBin(ix_new, iy_new);
        double sN2_prev = std::pow(hSigma_N_2D_flow->GetBinContent(bin_new), 2);
        double sZ2_prev = std::pow(hSigma_Z_2D_flow->GetBinContent(bin_new), 2);
        double sB2_prev = std::pow(hSigma_B_2D_flow->GetBinContent(bin_new), 2);

        hSigma_N_2D_flow->SetBinContent(bin_new, std::sqrt(sN2_prev + sN*sN));
        hSigma_Z_2D_flow->SetBinContent(bin_new, std::sqrt(sZ2_prev + sZ*sZ));
        hSigma_B_2D_flow->SetBinContent(bin_new, std::sqrt(sB2_prev + sB*sB));

        // --- R: weighted average ---
        if (N > 0) {
            double prevR  = hR_2D_flow->GetBinContent(ix_new, iy_new);
            double prevW  = hN_2D_flow->GetBinContent(ix_new, iy_new) - N;
            double newR = (prevR * prevW + R * N) / (prevW + N);
            hR_2D_flow->SetBinContent(ix_new, iy_new, newR);
        }
    }
}

TH1D* hZ_projX_flow = hZ_2D_flow->ProjectionX("hZ_projX_flow");
TH1D* hB_projX_flow = hB_2D_flow->ProjectionX("hB_projX_flow");
TH1D* hN_projX_flow = hN_2D_flow->ProjectionX("hN_projX_flow");
TH1D* hR_projX_flow = hR_2D_flow->ProjectionX("hR_projX_flow");

TH1D* hZ_projY_flow = hZ_2D_flow->ProjectionY("hZ_projY_flow");
TH1D* hB_projY_flow = hB_2D_flow->ProjectionY("hB_projY_flow");
TH1D* hR_projY_flow = hR_2D_flow->ProjectionY("hR_projY_flow");
TH1D* hN_projY_flow = hN_2D_flow->ProjectionY("hN_projY_flow");

hN_projX_flow->Scale(norm);
hZ_projX_flow->Scale(norm);
hB_projX_flow->Scale(norm);

hN_projY_flow->Scale(norm);
hZ_projY_flow->Scale(norm);
hB_projY_flow->Scale(norm);

hN_projX_flow->Write();
hZ_projX_flow->Write();
hB_projX_flow->Write();
hR_projX_flow->Write();

hN_projY_flow->Write();
hZ_projY_flow->Write();
hB_projY_flow->Write();
hR_projY_flow->Write();

for (int ix = 1; ix <= hZ_2D_flow->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= hZ_2D_flow->GetNbinsY(); ++iy) {

        double Z = hZ_2D_flow->GetBinContent(ix, iy);
        double B = hB_2D_flow->GetBinContent(ix, iy);

        double sZ = hSigma_Z_2D_flow->GetBinContent(ix, iy);
        double sB = hSigma_B_2D_flow->GetBinContent(ix, iy);

        if (B > 0) {
            double SoB = Z / B;
            hS_over_B_2D_flow->SetBinContent(ix, iy, SoB);

            // error propagation
            double sSoB = 0.0;
            if (std::abs(Z) > 0) {
                sSoB = SoB * std::sqrt(
                    (sZ*sZ)/(Z*Z) +
                    (sB*sB)/(B*B)
                );
            }
            if (Z + B > 0) {
                double signif = Z / std::sqrt(Z + B);
                hSignificance_2D_flow->SetBinContent(ix, iy, signif);
            }
            hSigma_S_over_B_2D_flow->SetBinContent(ix, iy, sSoB);
        } else {
            hS_over_B_2D_flow->SetBinContent(ix, iy, 0.0);
            hSigma_S_over_B_2D_flow->SetBinContent(ix, iy, 0.0);
        }
    }
}
// ========================================= rebin the flow variables since low statistics ===========================================//
// here we rebin the flow variables in rough bins for better visualization done after calcualting the variables to keep R factor fine //
// ========================================= rebin the flow variables since low statistics ===========================================//
std::vector<double> xbins_flow = {
    0.0,
    0.14,
    1.20,
    2.70,
    3.20,
    4.00
};
std::vector<std::pair<int,int>> meeSlices;



for (size_t i = 0; i < xbins_flow.size() - 1; ++i) {

    int lo = FindBinLow (xax, xbins_flow[i]);
    int hi = FindBinHigh(xax, xbins_flow[i+1]);

    lo = std::max(lo, 1);
    hi = std::min(hi, xax->GetNbins());

    meeSlices.emplace_back(lo, hi);
}

std::vector<double> ybins_flow_pion_mass = {
    0.0, 0.3, 0.5, 0.7,
    1.0,
    1.5,  2.0, 2.5, 3.5, 
    5.0,
    10.0
};

std::vector<double> ybins_flow_low_mass = {
    0.0, 
    1.0,
    2.0, 4.0, 
    6.0,
    10.0
};

std::vector<double> ybins_flow_j_psi_mass = {
    0.0, 
    2.0, 4.0, 
    6.0,
    10.0
};
std::vector<int> slices_to_use = {1, 2, 4}; // slices you want

for (size_t s = 0; s < slices_to_use.size(); ++s) {
    int ix = slices_to_use[s]; // mass slice
    const std::vector<double>* ybins_flow = nullptr;

    if (ix == 1) {                 // pion mass region
        ybins_flow = &ybins_flow_pion_mass;
    } 
    else if (ix == 2) {            // low mass region
        ybins_flow = &ybins_flow_low_mass;
    } 
    else if (ix == 4) {            // J/psi mass region
        ybins_flow = &ybins_flow_j_psi_mass;
    } 
    else {
        continue; // skip other slices
    }
       
    


TH3D* h_rebinned_uls_flow_visual = new TH3D("h_uls_flow", "ULS (flow bins);m_{ee};p_{T};#Delta#phi",xbins_flow.size()-1, xbins_flow.data(),ybins_flow->size()-1, ybins_flow->data(),zbins_flow.size()-1, zbins_flow.data());
TH3D* h_rebinned_lspp_flow_visual = new TH3D("h_lspp_flow", "LS++ (flow bins);m_{ee};p_{T};#Delta#phi",xbins_flow.size()-1, xbins_flow.data(),ybins_flow->size()-1, ybins_flow->data(),zbins_flow.size()-1, zbins_flow.data());
TH3D* h_rebinned_lsmm_flow_visual = new TH3D("h_lsmm_flow", "LS-- (flow bins);m_{ee};p_{T};#Delta#phi",xbins_flow.size()-1, xbins_flow.data(),ybins_flow->size()-1, ybins_flow->data(),zbins_flow.size()-1, zbins_flow.data());
TH2D* hN_2D_flow_visual = new TH2D("hN_2D_flow_visual", "N (flow bins);m_{ee};p_{T}",xbins_flow.size()-1, xbins_flow.data(),ybins_flow->size()-1, ybins_flow->data());
TH2D* hZ_2D_flow_visual = new TH2D("hZ_2D_flow_visual", "Z (flow bins);m_{ee};p_{T}",xbins_flow.size()-1, xbins_flow.data(),ybins_flow->size()-1, ybins_flow->data());
TH2D* hB_2D_flow_visual = new TH2D("hB_2D_flow_visual", "B (flow bins);m_{ee};p_{T}",xbins_flow.size()-1, xbins_flow.data(),ybins_flow->size()-1, ybins_flow->data());
TH2D* hSigma_Z_2D_flow_visual = new TH2D("hSigma_Z_2D_flow_visual", "σZ (flow bins);m_{ee};p_{T}",xbins_flow.size()-1, xbins_flow.data(),ybins_flow->size()-1, ybins_flow->data());
TH2D* hSigma_B_2D_flow_visual = new TH2D("hSigma_B_2D_flow_visual", "σB (flow bins);m_{ee};p_{T}",xbins_flow.size()-1, xbins_flow.data(),ybins_flow->size()-1, ybins_flow->data());
TH2D* hSigma_N_2D_flow_visual = new TH2D("hSigma_N_2D_flow_visual", "σN (flow bins);m_{ee};p_{T}",xbins_flow.size()-1, xbins_flow.data(),ybins_flow->size()-1, ybins_flow->data());

std::vector<std::pair<int,int>> ptSlices;

for (size_t i = 0; i < ybins_flow->size() - 1; ++i) {

    int lo = FindBinLow(yax, ybins_flow->at(i));       // use ->at(i)
    int hi = FindBinHigh(yax, ybins_flow->at(i + 1));  // use ->at(i+1)

    lo = std::max(lo, 1);
    hi = std::min(hi, yax->GetNbins());

    ptSlices.emplace_back(lo, hi);
}
for (size_t i = 0; i < meeSlices.size(); ++i) {
    for (size_t j = 0; j < ptSlices.size(); ++j) {
        for (size_t k = 0; k < qvecSlicesFlow.size(); ++k) {

            double sum_uls  = 0.0;
            double sum_lspp = 0.0;
            double sum_lsmm = 0.0;


            int xlo = meeSlices[i].first;
            int xhi = meeSlices[i].second;

            int ylo = ptSlices[j].first;
            int yhi = ptSlices[j].second;

            int zlo = qvecSlicesFlow[k].first;
            int zhi = qvecSlicesFlow[k].second;

            for (int ix = xlo; ix <= xhi; ++ix) {
                for (int iy = ylo; iy <= yhi; ++iy) {
                    for (int iz = zlo; iz <= zhi; ++iz) {

                        sum_uls  += h_rebinned_uls_flow ->GetBinContent(ix,iy,iz);
                        sum_lspp += h_rebinned_lspp_flow->GetBinContent(ix,iy,iz);
                        sum_lsmm += h_rebinned_lsmm_flow->GetBinContent(ix,iy,iz);
                    }
                }
            }
            h_rebinned_uls_flow_visual->SetBinContent(i+1, j+1, k+1, sum_uls);
            h_rebinned_lspp_flow_visual->SetBinContent(i+1, j+1, k+1, sum_lspp);
            h_rebinned_lsmm_flow_visual->SetBinContent(i+1, j+1, k+1, sum_lsmm);
        }
    }
}

for(size_t i=0;i<meeSlices.size();++i){
    for(size_t j=0;j<ptSlices.size();++j){

        double sumZ=0;
        double sumB=0;
        double sumN=0;

        double errZ2=0;
        double errB2=0;
        double errN2=0;

        int xlo = meeSlices[i].first;
        int xhi = meeSlices[i].second;

        int ylo = ptSlices[j].first;
        int yhi = ptSlices[j].second;

        for(int ix=xlo; ix<=xhi; ++ix){
            for(int iy=ylo; iy<=yhi; ++iy){

                sumZ += hZ_2D_flow->GetBinContent(ix,iy);
                sumB += hB_2D_flow->GetBinContent(ix,iy);
                sumN += hN_2D_flow->GetBinContent(ix,iy);

                double eZ = hSigma_Z_2D_flow->GetBinContent(ix,iy);
                double eB = hSigma_B_2D_flow->GetBinContent(ix,iy);
                double eN = hSigma_N_2D_flow->GetBinContent(ix,iy);

                errZ2 += eZ*eZ;
                errB2 += eB*eB;
                errN2 += eN*eN;
            }
        }
        hZ_2D_flow_visual->SetBinContent(i+1,j+1,sumZ);
        hB_2D_flow_visual->SetBinContent(i+1,j+1,sumB);
        hN_2D_flow_visual->SetBinContent(i+1,j+1,sumN);

        hSigma_Z_2D_flow_visual->SetBinContent(i+1,j+1,std::sqrt(errZ2));
        hSigma_B_2D_flow_visual->SetBinContent(i+1,j+1,std::sqrt(errB2));
        hSigma_N_2D_flow_visual->SetBinContent(i+1,j+1,std::sqrt(errN2));
    }
}





int nbz2 = h_rebinned_lsmm_flow_visual->GetNbinsZ();
int nBinsYvariable = h_rebinned_lsmm_flow_visual->GetNbinsY();
std::vector<double> yedges2(nBinsYvariable+1);
for (int i = 0; i <= nBinsYvariable; ++i)
yedges2[i] = h_rebinned_uls_flow_visual->GetYaxis()->GetBinLowEdge(i+1);
// Loop over slices (each slice was previously a single X-bin range)

for (size_t s = 0; s < slices_to_use.size(); ++s) {
    int ix = slices_to_use[s]; // <- use the value from slices_to_use

    int nY = nBinsYvariable;

    for (int iy = 1; iy <= nY; ++iy) {
        for (int iz = 1; iz <= nbz2; ++iz) {
            double N_pm = h_rebinned_uls_flow_visual->GetBinContent(ix, iy, iz);
            double N_pp = h_rebinned_lspp_flow_visual->GetBinContent(ix, iy, iz);
            double N_mm = h_rebinned_lsmm_flow_visual->GetBinContent(ix, iy, iz);

            if (N_pm > 0 || N_pp > 0 || N_mm > 0) {

            }
        }
    }
}


int nBinsX_flow = h_rebinned_uls_flow_visual->GetNbinsX();
int nBinsY_flow = h_rebinned_uls_flow_visual->GetNbinsY();

// Vectors to store content and errors for this plot_standard_large_cuts_cuts_testing
std::vector<double> Z_pt_per_slice(nBinsY_flow, 0.0);
std::vector<double> Z_pt_err_per_slice(nBinsY_flow, 0.0);

// Fill Z values and errors by summing over X bins
for (int iy = 0; iy < nBinsY_flow; ++iy) {
    double Z_sum = 0.0;
    double sigma2_sum = 0.0;

    for (int ix = 0; ix < nBinsX_flow; ++ix) {
        double binContent = hZ_2D_flow_visual->GetBinContent(ix+1, iy+1);
        double binError   = hSigma_Z_2D_flow_visual->GetBinContent(ix+1, iy+1);

        Z_sum += binContent;
        sigma2_sum += binError * binError;
    }
    double binWidth = h_rebinned_uls_flow_visual->GetYaxis()->GetBinWidth(iy+1);
    Z_pt_per_slice[iy]     = Z_sum/(multiplicity*binWidth);
    Z_pt_err_per_slice[iy] = std::sqrt(sigma2_sum)/(multiplicity*binWidth);
}

// Vectors for pT bin centers and half-widths
std::vector<double> pT_bin_centers(nBinsY_flow, 0.0);
std::vector<double> pT_bin_halfwidths(nBinsY_flow, 0.0);

for (int iy = 0; iy < nBinsY_flow; ++iy) {
    double lowEdge  = h_rebinned_lsmm_flow_visual->GetYaxis()->GetBinLowEdge(iy+1);
    double highEdge = h_rebinned_lsmm_flow_visual->GetYaxis()->GetBinUpEdge(iy+1);

    pT_bin_centers[iy]     = 0.5 * (lowEdge + highEdge);
    pT_bin_halfwidths[iy]  = 0.5 * (highEdge - lowEdge);
}

// Create TGraphErrors for this plot_standard_large_cuts_cuts_testing
TGraphErrors* g_pt_per_slice = new TGraphErrors(
    nBinsY_flow,
    pT_bin_centers.data(),
    Z_pt_per_slice.data(),
    pT_bin_halfwidths.data(),
    Z_pt_err_per_slice.data()
);

g_pt_per_slice->SetLineColor(kBlack);
g_pt_per_slice->SetMarkerColor(kBlue);
g_pt_per_slice->SetMarkerStyle(26);
g_pt_per_slice->SetTitle("");

std::vector<double> N_flow_pt(nBins_rebin_Y, 0.0);
std::vector<double> eyN_flow_pt(nBins_rebin_Y, 0.0);


// Loop over X bins
for (int ix = 0; ix < nBins_rebin_Y; ix++) {
double Z_sum = 0.0, B_sum = 0.0, N_sum = 0.0;
double sigmaZ2 = 0.0, sigmaB2 = 0.0, sigmaN2 = 0.0;

for (int iz = 0; iz < nBins_rebin_X; iz++) {
    Z_sum += hN_2D_flow_visual->GetBinContent(iz+1,ix+1);


    sigmaZ2 += hSigma_N_2D_flow_visual->GetBinContent(iz+1,ix+1)*hSigma_N_2D_flow_visual->GetBinContent(iz+1,ix+1);

}
double binWidth = h_rebinned_uls_3D->GetYaxis()->GetBinWidth(ix+1);
N_flow_pt[ix] = Z_sum/(multiplicity*binWidth);
eyN_flow_pt[ix] = std::sqrt(sigmaZ2)/(multiplicity*binWidth);
}

std::vector<double> B_flow_pt(nBins_rebin_Y, 0.0);
std::vector<double> eyB_flow_pt(nBins_rebin_Y, 0.0);


// Loop over X bins
for (int ix = 0; ix < nBins_rebin_Y; ix++) {
double Z_sum = 0.0;
double sigmaZ2 = 0.0;

for (int iz = 0; iz < nBins_rebin_X; iz++) {
    Z_sum += hB_2D_flow_visual->GetBinContent(iz+1,ix+1);


    sigmaZ2 += hSigma_B_2D_flow_visual->GetBinContent(iz+1,ix+1)*hSigma_B_2D_flow_visual->GetBinContent(iz+1,ix+1);

}    
double binWidth = h_rebinned_uls_3D->GetYaxis()->GetBinWidth(ix+1);
B_flow_pt[ix] = Z_sum/(multiplicity*binWidth);
eyB_flow_pt[ix] = std::sqrt(sigmaZ2)/(multiplicity*binWidth);
}

std::vector<double> Significance_pt_spectra_flow(nBins_rebin_Y, 0.0);
std::vector<double> eSignificance_pt_spectra_flow(nBins_rebin_Y, 0.0);

for (int ix = 0; ix < nBins_rebin_Y; ++ix) {

    double Z_sum = 0.0;
    double N_sum = 0.0;
    for (int iz = 0; iz < nBins_rebin_X; ++iz) {

        double Z = hZ_2D_flow_visual->GetBinContent(iz+1, ix+1);
        double N = hN_2D_flow_visual->GetBinContent(iz+1, ix+1);


        Z_sum += Z;
        N_sum += N;
    }

    if (N_sum > 0 && Z_sum > 0) {
        double S = Z_sum / std::sqrt(N_sum);
        double eS = std::sqrt(
            (Z_sum / N_sum) +
            (Z_sum * Z_sum) / (4.0 * N_sum * N_sum)
        );

        Significance_pt_spectra_flow[ix]  = S;
        eSignificance_pt_spectra_flow[ix] = eS;

    }
}

std::vector<double> S_over_B_flow(nBins_rebin_Y, 0.0);
std::vector<double> eyS_over_B_flow(nBins_rebin_Y, 0.0);


for (int iz = 0; iz < nBins_rebin_Y; ++iz) {

    double Z_sum = 0.0;
    double B_sum = 0.0;
    double sigmaZ2 = 0.0;
    double sigmaB2 = 0.0;

    for (int ix = 0; ix < nBins_rebin_X; ++ix) {

        double Z = hZ_2D_flow_visual->GetBinContent(ix+1, iz+1);
        double B = hB_2D_flow_visual->GetBinContent(ix+1, iz+1);

        double sigmaZ = hSigma_Z_2D_flow_visual->GetBinContent(ix+1, iz+1);
        double sigmaB = hSigma_B_2D_flow_visual->GetBinContent(ix+1, iz+1);

        Z_sum += Z;
        B_sum += B;

        sigmaZ2 += sigmaZ * sigmaZ;
        sigmaB2 += sigmaB * sigmaB;
    }

    if (B_sum > 0 && Z_sum > 0) {
        S_over_B_flow[iz] = Z_sum / B_sum;

        eyS_over_B_flow[iz] =
            (Z_sum / B_sum) *
            std::sqrt(
                (sigmaZ2 / (Z_sum * Z_sum)) +
                (sigmaB2 / (B_sum * B_sum))
            );
    }
}
TGraphErrors* gNerr_flow_pt = new TGraphErrors(
    nBins_rebin_Y,
    pT_bin_centers.data(),
    N_flow_pt.data(),   // Z values for chosen slice
    pT_bin_halfwidths.data(),
    eyN_flow_pt.data() 
);
gNerr_flow_pt->SetLineColor(kBlack);
gNerr_flow_pt->SetMarkerColor(kGreen);
gNerr_flow_pt->SetMarkerStyle(24);
gNerr_flow_pt->SetTitle("");

TGraphErrors* gBerr_flow_pt = new TGraphErrors(
    nBins_rebin_Y,
    pT_bin_centers.data(),
    B_flow_pt.data(),   // Z values for chosen slice
    pT_bin_halfwidths.data(),
    eyB_flow_pt.data() 
);
gBerr_flow_pt->SetLineColor(kBlack);
gBerr_flow_pt->SetMarkerColor(kRed);
gBerr_flow_pt->SetMarkerStyle(25);
gBerr_flow_pt->SetTitle(""); 

TGraphErrors* gS_over_B_err_flow_pt = new TGraphErrors(
    nBins_rebin_Y,
    pT_bin_centers.data(),
    S_over_B_flow.data(),   // Z values for chosen slice
    pT_bin_halfwidths.data(),
    eyS_over_B_flow.data() 
);
gS_over_B_err_flow_pt->SetLineColor(kBlack);
gS_over_B_err_flow_pt->SetMarkerColor(kBlue);
gS_over_B_err_flow_pt->SetMarkerStyle(25);
gS_over_B_err_flow_pt->SetTitle(""); 

TGraphErrors* gSignificance_err_flow_pt = new TGraphErrors(
    nBins_rebin_Y,
    pT_bin_centers.data(),
    Significance_pt_spectra_flow.data(),   // Z values for chosen slice
    pT_bin_halfwidths.data(),
    eSignificance_pt_spectra_flow.data() 
);
gSignificance_err_flow_pt->SetLineColor(kBlack);
gSignificance_err_flow_pt->SetMarkerColor(kRed);
gSignificance_err_flow_pt->SetMarkerStyle(25);
gSignificance_err_flow_pt->SetTitle("");

// Create canvas
TCanvas* c_pt_per_slice = new TCanvas("c_pt_per_slice", "", 800, 600);

c_pt_per_slice->SetBottomMargin(0.1); // small margin between pads
c_pt_per_slice->SetLeftMargin(0.11);
c_pt_per_slice->SetTicks(1,1);
c_pt_per_slice->SetLogy();
c_pt_per_slice->Draw();
c_pt_per_slice->cd();
g_pt_per_slice->Draw("AP");
g_pt_per_slice->SetStats(0);
g_pt_per_slice->GetXaxis()->SetTitle("p_{T} (GeV/c)");
g_pt_per_slice->GetYaxis()->SetTitle("#frac{1}{N_{ev}} #frac{dN}{dp_{T}} (GeV/c)^{-1}");
g_pt_per_slice->SetTitle("");
g_pt_per_slice->GetYaxis()->SetRangeUser(3e-8, 1);  
g_pt_per_slice->GetXaxis()->SetRangeUser(0, 10.0);  
gBerr_flow_pt->Draw("P same");
gNerr_flow_pt->Draw("P same");

// Legend
TLegend* legend_pt_per_slice = new TLegend(0.6, 0.7, 0.89, 0.89);
legend_pt_per_slice->AddEntry(g_pt_per_slice, "signal", "lp");
legend_pt_per_slice->AddEntry(gBerr_flow_pt, "background", "lp");
legend_pt_per_slice->AddEntry(gNerr_flow_pt, "background+signal", "lp");
legend_pt_per_slice->SetTextSize(0.03);
legend_pt_per_slice->SetBorderSize(0);
legend_pt_per_slice->Draw();

// Latex text
latex.SetNDC();
latex.SetTextSize(0.03);
latex.SetTextAlign(13);
latex.DrawLatex(0.15, 0.88, "Pb-Pb #sqrt{s_{NN}} = 5.36 TeV");
latex.DrawLatex(0.15, 0.85, cenLabels[ip].c_str());

c_pt_per_slice->SaveAs(("Pt_spectra_per_centrality"+ cenLabels[ip]  + ".png").c_str());


TCanvas* c_pt_per_slice_S_over_B = new TCanvas("c_pt_per_slice_S_over_B", "", 800, 600);
c_pt_per_slice_S_over_B->SetBottomMargin(0.1);
c_pt_per_slice_S_over_B->SetLeftMargin(0.11);
c_pt_per_slice_S_over_B->SetTicks(1,1);
c_pt_per_slice_S_over_B->SetLogy();
c_pt_per_slice_S_over_B->Draw();
c_pt_per_slice_S_over_B->cd();

gS_over_B_err_flow_pt->Draw("AP");
gS_over_B_err_flow_pt->SetStats(0);
gS_over_B_err_flow_pt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
gS_over_B_err_flow_pt->GetYaxis()->SetTitle("S/B");
gS_over_B_err_flow_pt->SetTitle("");
gS_over_B_err_flow_pt->GetXaxis()->SetRangeUser(0,10.0);
gS_over_B_err_flow_pt->GetYaxis()->SetRangeUser(5e-2, 1);

latex.SetNDC();
latex.SetTextSize(0.03);
latex.SetTextAlign(13);
latex.DrawLatex(0.15, 0.88, "Pb-Pb #sqrt{s_{NN}} = 5.36 TeV");
latex.DrawLatex(0.15, 0.85, cenLabels[ip].c_str());


// Save plot_standard_large_cuts_cuts_testing
c_pt_per_slice_S_over_B->SaveAs(("Pt_spectra_per_centrality_S_over_B"+ cenLabels[ip]  + ".png").c_str());
std::vector<std::string> massrange = {"pion_", "intermediate_", "j_psi_"};

TCanvas* csignificance_pt_flow_spectra = new TCanvas("csignificance_pt_spectra", "", 800, 800);

csignificance_pt_flow_spectra->SetTicks(1,1);
csignificance_pt_flow_spectra->SetLogy();

gSignificance_err_flow_pt->Draw("AP");
gSignificance_err_flow_pt->SetStats(0);

gSignificance_err_flow_pt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
gSignificance_err_flow_pt->GetYaxis()->SetTitle("Significance per unit #it{p}_{T} (GeV/#it{c})^{-1}");
gSignificance_err_flow_pt->SetTitle("");
gSignificance_err_flow_pt->GetYaxis()->SetRangeUser(1e-2, 1e4);
gSignificance_err_flow_pt->GetXaxis()->SetRangeUser(0.0, 10.0);
latex.SetNDC();
latex.SetTextSize(0.03);
latex.SetTextAlign(13);
latex.DrawLatex(0.15, 0.88, "Pb-Pb #sqrt{s_{NN}} = 5.36 TeV");
latex.DrawLatex(0.15, 0.85, cenLabels[ip].c_str());

// Save
csignificance_pt_flow_spectra->SaveAs(("pt_spectra_flow_significance_" + cenLabels[ip]  + ".png").c_str());
    double mass_center = h_rebinned_uls_flow_visual->GetXaxis()->GetBinCenter(ix);
    double binWidthX = h_rebinned_uls_flow_visual->GetXaxis()->GetBinWidth(ix);

    int nY = nBinsYvariable; // number of Y bins

    std::vector<double> v2Pt(nY, 0.0);
    std::vector<double> v2PtErr(nY, 0.0);
    std::vector<double> v2BPt(nY, 0.0);
    std::vector<double> v2BPtErr(nY, 0.0);
    std::vector<double> v2SPt(nY, 0.0);
    std::vector<double> v2SPtErr(nY, 0.0);

    // Loop over Y bins
    for (int iy = 1; iy <= nY; ++iy) {


        double sum_wzSB = 0.0, sum_wSB = 0.0, sum_wz2SB = 0.0;
        double sum_wzB  = 0.0, sum_wB  = 0.0, sum_wz2B  = 0.0;
        double sum_pp   = 0.0, sum_z_pp   = 0.0, sum_z2_pp   = 0.0;
        double sum_mm   = 0.0, sum_z_mm   = 0.0, sum_z2_mm   = 0.0;

        // Loop over Z bins
        for (int iz = 1; iz <= nbz2; ++iz) {
            // --- S+B ---
            double N_pm_same = h_rebinned_uls_flow_visual->GetBinContent(ix, iy, iz);

            if (N_pm_same > 0) {
                double zcenter = h_rebinned_uls_flow_visual->GetZaxis()->GetBinCenter(iz);
                sum_wSB   += N_pm_same;
                sum_wzSB  += N_pm_same * zcenter;
                sum_wz2SB += N_pm_same * zcenter * zcenter;
            }

            // --- Background B only ---
            double zcenter = h_rebinned_lspp_flow_visual->GetZaxis()->GetBinCenter(iz);

            double N_pp = h_rebinned_lspp_flow_visual->GetBinContent(ix, iy, iz);
            if (N_pp > 0) {
                sum_pp   += N_pp;
                sum_z_pp += N_pp * zcenter;
                sum_z2_pp += N_pp * zcenter * zcenter;
            }
        
            double N_mm = h_rebinned_lsmm_flow_visual->GetBinContent(ix, iy, iz);
            if (N_mm > 0) {
                sum_mm   += N_mm;
                sum_z_mm += N_mm * zcenter;
                sum_z2_mm += N_mm * zcenter * zcenter;
            }
        }

        // Compute v2(S+B)
        if (sum_wSB > 0) {
            double mean_z = sum_wzSB / sum_wSB;
            double var_z  = sum_wz2SB / sum_wSB - mean_z * mean_z;
            var_z = (var_z < 0 && var_z > -1e-12) ? 0.0 : std::max(var_z, 0.0);
            double err_mean = std::sqrt(var_z / sum_wSB);

            v2Pt[iy - 1]    = mean_z / Rval;
            v2PtErr[iy - 1] = err_mean / Rval;
        }

        // Compute v2(B)
        if (sum_pp > 0 && sum_mm > 0) {

            double mean_pp = sum_z_pp / sum_pp;
            double mean_mm = sum_z_mm / sum_mm;
        
            double var_pp = sum_z2_pp / sum_pp - mean_pp * mean_pp;
            double var_mm = sum_z2_mm / sum_mm - mean_mm * mean_mm;
        
            var_pp = std::max(var_pp, 0.0);
            var_mm = std::max(var_mm, 0.0);
        
            double err_pp = std::sqrt(var_pp / sum_pp);
            double err_mm = std::sqrt(var_mm / sum_mm);
        
            // background v_n
            v2BPt[iy - 1] = (mean_pp + mean_mm) / (2.0 * Rval);
        
            // error propagation
            v2BPtErr[iy - 1] =
                std::sqrt(err_pp * err_pp + err_mm * err_mm) / (2.0 * Rval);
        }
        

        // Compute v2(S)
        double S = hZ_2D_flow_visual->GetBinContent(ix, iy);
        double B = hB_2D_flow_visual->GetBinContent(ix, iy);
        double sigmaS = hSigma_Z_2D_flow_visual->GetBinContent(ix,iy);
        double sigmaB = hSigma_B_2D_flow_visual->GetBinContent(ix, iy);
        double T = S + B;

        if (T > 0) {
            double a = S / T;

            double da_dS = B / (T*T);
            double da_dB = -S / (T*T);
            double var_a = da_dS*da_dS*sigmaS*sigmaS + da_dB*da_dB*sigmaB*sigmaB;
            double sigma_a = (var_a > 0) ? std::sqrt(var_a) : 0.0;

            double v2S_val = (a > 0) ? (v2Pt[iy-1] - (1.0 - a)*v2BPt[iy-1])/a : 0.0;
            double term_SB = v2PtErr[iy-1] / a;
            double term_B  = ((1.0 - a) * v2BPtErr[iy-1]) / a;
            double term_a  = ((v2BPt[iy-1] - v2Pt[iy-1]) * sigma_a) / (a*a);
            double var_v2S = term_SB*term_SB + term_B*term_B + term_a*term_a;
            v2SPt[iy-1] = v2S_val;
            v2SPtErr[iy-1] = (var_v2S > 0) ? std::sqrt(var_v2S) : 0.0;
            std::cout << "Bin " << iy  << ", v2(S) = " << v2SPt[iy-1] << std::endl;
        }
    }
    std::cout << "Finished Y-bin loop for X-bin " << ix << std::endl;
    for (int iz = 0; iz < nBinsYvariable; ++iz) {

        double Z_sum = 0.0;
        double B_sum = 0.0;
        double sigmaZ2 = 0.0;
        double sigmaB2 = 0.0;
    
        
    
            double Z = v2SPt[iz];

            double B = v2BPt[iz];

            double sigmaZ = v2SPtErr[iz];

            double sigmaB = v2BPtErr[iz];
    
            Z_sum += Z;
            B_sum += B;
    
            sigmaZ2 += sigmaZ * sigmaZ;
            sigmaB2 += sigmaB * sigmaB;
        
    
        if (B_sum > 0 && Z_sum > 0) {
            S_over_B_flow[iz] = Z_sum / B_sum;
    
            eyS_over_B_flow[iz] =
                (Z_sum / B_sum) *
                std::sqrt(
                    (sigmaZ2 / (Z_sum * Z_sum)) +
                    (sigmaB2 / (B_sum * B_sum))
                );
        }
    }

        // --- Fill TH1Ds for this slice ---
        std::vector<std::pair<std::vector<double>, std::vector<double>>> vecs = {
            {v2Pt, v2PtErr},
            {v2BPt, v2BPtErr},
            {v2SPt, v2SPtErr}
        };


    
    
        std::vector<std::string> labels = {"Signal+Background", "Background_", "Signal_"};

        std::vector<std::string> massrange2 = {"0 < m_{ee} < 0.14", "0. 14 < m_{ee} < 1.2", "2.7 < m_{ee} < 3.2"};
            

        for (size_t i=0; i<vecs.size(); ++i) {


            std::string name = "plot_standard_large_cuts_v2_pt_cuts_test_" + massrange[s] + "mass_range" + labels[i];
            TH1D* h = new TH1D(name.c_str(), ";p_{T} GeV/c ;v_{2}", nBinsYvariable, yedges2.data());
            for (int iy=0; iy<nBinsYvariable; ++iy) {
                h->SetBinContent(iy+1, vecs[i].first[iy]);
                h->SetBinError(iy+1, vecs[i].second[iy]);
            }


            TCanvas* c = new TCanvas(("c"+name).c_str(), (name).c_str(), 800, 800);
            c->SetTicks(1,1);
            c->SetLeftMargin(0.16);
            c->SetBottomMargin(0.12);
            c->SetRightMargin(0.04);
            c->SetTopMargin(0.05);
            h->SetMarkerColor(kRed+1);
            h->SetMarkerSize(1.5);
            h->SetDirectory(0);              // VERY IMPORTANT (prevents deletion)
            hFlowVec.push_back(h);
            histLabels.push_back(labels[i] + ", " + massrange2[s] + ", " + cenLabels[ip]);
            histTypeIndex.push_back(i); // i = type index
            h->SetLineColor(colors[ip]);
            h->SetMarkerColor(colors[ip]);
            h->SetMarkerStyle(24 + ip);      
            h->SetMarkerSize(1.2);
            h->SetLineWidth(2);
            h->GetXaxis()->SetTitleSize(0.045);
            h->GetYaxis()->SetTitleSize(0.045);
            h->GetXaxis()->SetLabelSize(0.04);
            h->GetYaxis()->SetLabelSize(0.04);
            h->SetMaximum(0.8);
            h->SetMinimum(-0.1);
            h->SetStats(0);
            h->Draw("E1");

            TLine* l0 = new TLine(0, 0, 10, 0);
            l0->SetLineStyle(3);
            l0->SetLineWidth(2);
            l0->SetLineColor(kGray+2);
            l0->Draw("same");

            latex.SetNDC();
            latex.SetTextSize(0.03);
            latex.SetTextAlign(13);
            latex.DrawLatex(0.18, 0.84, ("mass region " + massrange2[s] + " GeV/c^{2} " + labels[i]).c_str());
            latex.DrawLatex(0.18, 0.89, (cenLabels[ip] + " Pb-Pb #sqrt{s_{NN}} = 5.36 TeV").c_str());
            c->SaveAs((name + cenLabels[ip] +".png").c_str());

            h->Write(("flow" + massrange[s] + labels[i] + cenLabels[ip]).c_str());   // histogram

            delete c;

        }


        int nBins_rebin_X_big = h_rebinned_uls_flow_visual->GetNbinsX();
        int nBins_rebin_Y_big = h_rebinned_uls_flow_visual->GetNbinsY();
        
        // Prepare vectors
        std::vector<double> Z_rebinned_pt_slice(nBins_rebin_Y_big, 0.0);
        std::vector<double> eyZ_rebinned_pt_slice(nBins_rebin_Y_big, 0.0);
        std::vector<double> y_rebinned_slice(nBins_rebin_Y_big, 0.0);
        std::vector<double> ey_rebinned_slice(nBins_rebin_Y_big, 0.0);
        std::vector<double> B_pt(nBins_rebin_Y_big, 0.0),  B_pt_err(nBins_rebin_Y_big, 0.0);
        std::vector<double> N_pt(nBins_rebin_Y_big, 0.0),  N_pt_err(nBins_rebin_Y_big, 0.0);
        TH1D* hS_over_B_slice = new TH1D("hSigma_S_over_B_1D_slice", "σR (flow bins);m_{ee};p_{T}",ybins_flow->size()-1, ybins_flow->data());


        
        for (int iy = 0; iy < nBins_rebin_Y_big; ++iy) {
            // Take content from 2D histogram for this slice only
            double Z = hZ_2D_flow_visual->GetBinContent(ix+1, iy+1);
            double N = hN_2D_flow_visual->GetBinContent(ix+1, iy+1);
            double B = hB_2D_flow_visual->GetBinContent(ix+1, iy+1);
            double sigma_Z = hSigma_Z_2D_flow_visual->GetBinContent(ix+1, iy+1);
            double sigma_N = hSigma_N_2D_flow_visual->GetBinContent(ix+1, iy+1);
            double sigma_B = hSigma_B_2D_flow_visual->GetBinContent(ix+1, iy+1);
            double binWidth = h_rebinned_uls_flow_visual->GetYaxis()->GetBinWidth(iy+1);


            
            Z_rebinned_pt_slice[iy] = Z/(multiplicity*binWidth*binWidthX);
            N_pt[iy] = N/(multiplicity*binWidth*binWidthX); 
            B_pt[iy] = B/(multiplicity*binWidth*binWidthX);
            eyZ_rebinned_pt_slice[iy] = sigma_Z/(multiplicity*binWidth*binWidthX);
            N_pt_err[iy] = sigma_N/(multiplicity*binWidth*binWidthX);
            B_pt_err[iy] = sigma_B/(multiplicity*binWidth*binWidthX);

            double S_over_B = 0.0;
            double S_over_B_err = 0.0;
            if (B_pt[iy] > 0) {
                S_over_B = Z_rebinned_pt_slice[iy] / B_pt[iy];  // signal / background
        
                // propagate errors: (Z/B)^2 * [(σ_Z/Z)^2 + (σ_B/B)^2]
                S_over_B_err = S_over_B * sqrt(
                    pow(eyZ_rebinned_pt_slice[iy]/Z_rebinned_pt_slice[iy], 2) +
                    pow(B_pt_err[iy]/B_pt[iy], 2)
                );
            }
            hS_over_B_slice->SetBinContent(iy+1, S_over_B);
            hS_over_B_slice->SetBinError(iy+1, S_over_B_err);
        
            // pT bin centers / half-widths
            double low = h_rebinned_lsmm_flow_visual->GetYaxis()->GetBinLowEdge(iy+1);
            double high = h_rebinned_lsmm_flow_visual->GetYaxis()->GetBinUpEdge(iy+1);
            y_rebinned_slice[iy] = 0.5 * (low + high);
            ey_rebinned_slice[iy] = 0.5 * (high - low);
        
        }
        hS_over_B_slice->Write();
        // Create TGraphErrors for this slice
        TGraphErrors* gZerr_slice = new TGraphErrors(
            nBins_rebin_Y_big,
            y_rebinned_slice.data(),
            Z_rebinned_pt_slice.data(),
            ey_rebinned_slice.data(),
            eyZ_rebinned_pt_slice.data()
        );
        
        gZerr_slice->SetLineColor(kBlack);
        gZerr_slice->SetMarkerColor(kBlue);
        gZerr_slice->SetMarkerStyle(26);
        gZerr_slice->SetTitle(" ");
        
        TGraphErrors* gNerr_flow_pt_slice = new TGraphErrors(
            nBins_rebin_Y_big,
            y_rebinned_slice.data(),
            N_pt.data(),   // Z values for chosen slice
            ey_rebinned_slice.data(),
            N_pt_err.data() 
        );
        gNerr_flow_pt_slice->SetLineColor(kBlack);
        gNerr_flow_pt_slice->SetMarkerColor(kGreen);
        gNerr_flow_pt_slice->SetMarkerStyle(24);
        gNerr_flow_pt_slice->SetTitle("");
        
        TGraphErrors* gBerr_flow_pt_slice = new TGraphErrors(
            nBins_rebin_Y_big,
            y_rebinned_slice.data(),
            B_pt.data(),   // Z values for chosen slice
            ey_rebinned_slice.data(),
            B_pt_err.data() 
        );
        gBerr_flow_pt_slice->SetLineColor(kBlack);
        gBerr_flow_pt_slice->SetMarkerColor(kRed);
        gBerr_flow_pt_slice->SetMarkerStyle(25);
        gBerr_flow_pt_slice->SetTitle(""); 
     
        // Draw
        TCanvas* c_slice = new TCanvas(("cPT_slice_" + std::to_string(ix)).c_str(), "", 800, 600);
        c_slice->SetLogy();
        c_slice->SetTicks(1,1);
        gZerr_slice->Draw("AP");
        gZerr_slice->SetStats(0);
        gZerr_slice->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        gZerr_slice->GetYaxis()->SetTitle("#frac{1}{N_{ev}} #frac{dN}{dp_{T}} (GeV/c)^{-1}");
        gZerr_slice->SetTitle("");
        gZerr_slice->GetYaxis()->SetRangeUser(3e-8, 1e-1); 
        gZerr_slice->GetXaxis()->SetRangeUser(0, 10); 
        gBerr_flow_pt_slice->Draw("P same");
        gNerr_flow_pt_slice->Draw("P same");
        
        latex.SetNDC();
        latex.SetTextSize(0.03);
        latex.SetTextAlign(13);
        latex.DrawLatex(0.15, 0.82, ("mass region " + massrange2[s] + " GeV/c^2").c_str());

        latex.DrawLatex(0.15, 0.88, (cenLabels[ip] + "Pb-Pb #sqrt{s_{NN}} = 5.36 TeV").c_str());
        
        TLegend* legend = new TLegend(0.7, 0.7, 0.89, 0.89);
        legend->AddEntry(gZerr_slice, "signal", "lp");
        legend->AddEntry(gBerr_flow_pt_slice, "background", "lp");
        legend->AddEntry(gNerr_flow_pt_slice, "signal + background", "lp");
        legend->SetBorderSize(0);
        legend->Draw();
        
        // Save
        c_slice->SaveAs(("plot_pt_spectra_" + cenLabels[ip] + std::to_string(ix) + ".png").c_str());

        delete c_slice;




        




    double xcenter = h_rebinned_uls_flow_visual->GetXaxis()->GetBinCenter(ix);
}
    }

    std::vector<std::string> labels = {"Signal+Background", "Background_", "Signal_"};
    std::vector<std::string> massrange = {"pion_", "intermediate_", "j_psi_"};
    std::vector<std::string> massrange2 = {"0 < m_{ee} < 0.14", "0. 14 < m_{ee} < 1.2", "2.7 < m_{ee} < 3.2"};
// Loop over mass regions
for (size_t s = 0; s < massrange2.size(); ++s) {
    // Loop over types
    for (size_t t = 0; t < labels.size(); ++t) {

        TCanvas* cAll = new TCanvas(Form("cAll_mass%d_type%d", int(s), int(t)),
                                    Form("Flow - %s - %s", labels[t].c_str(), massrange2[s].c_str()),
                                    800, 800);
        cAll->cd();
        cAll->SetTicks(1,1);        // ticks on all sides
cAll->SetLeftMargin(0.16);
cAll->SetBottomMargin(0.12);
cAll->SetRightMargin(0.05);
cAll->SetTopMargin(0.05);




        bool first = true;
        std::string excludeCent = cenLabels.back(); // last centrality, e.g., "50_100" 
        // Loop over all histograms and select only the ones matching this mass region and type
        for (size_t i = 0; i < hFlowVec.size(); ++i) {
            if (histTypeIndex[i] != t) continue;


            // Check if this histogram belongs to current mass region and type

            if (histLabels[i].find(massrange2[s]) == std::string::npos) continue;
            if (histLabels[i].find(excludeCent) != std::string::npos) continue;
            if (first) {
                hFlowVec[i]->Draw("E1");
                hFlowVec[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
                hFlowVec[i]->GetYaxis()->SetTitle("v_{2}");
                hFlowVec[i]->SetTitle("");
                first = false;
            } else {
                hFlowVec[i]->Draw("E1 SAME");
            }
        }

        // legend
        TLegend* leg = new TLegend(0.75, 0.78, 0.9, 0.92);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.03);

        for (size_t i = 0; i < hFlowVec.size(); ++i) {

            if (histTypeIndex[i] != t) continue;
            if (histLabels[i].find(massrange2[s]) == std::string::npos) continue;
            if (histLabels[i].find(excludeCent) != std::string::npos) continue;
        
            // --- extract ONLY the centrality ---
            std::string centLabel = "";
            for (const auto& cen : cenLabels) {
                if (histLabels[i].find(cen) != std::string::npos) {
                    centLabel = cen;
                    break;
                }
            }
        
            // Safety check
            if (centLabel.empty()) continue;
            std::string legendLabel = centLabel + " %";
            leg->AddEntry(hFlowVec[i], legendLabel.c_str(), "lp");
        }

        leg->Draw();

        TLine* l0 = new TLine(0, 0, 10, 0);
        l0->SetLineStyle(3);
        l0->SetLineWidth(2);
        l0->SetLineColor(kGray+2);
        l0->Draw("same");

        TLatex latex;
        latex.SetNDC();
        latex.SetTextFont(42);
        latex.SetTextSize(0.03);
        
        latex.DrawLatex(0.2, 0.89, "Pb-Pb #sqrt{s_{NN}} = 5.36 TeV");
        latex.DrawLatex(0.2, 0.85, "flow signal p_{T,e} > 0.2 GeV/c, |#eta_{e}| < 0.8");
        latex.DrawLatex(0.2, 0.81, massrange2[s].c_str());

        // Save canvas
        cAll->SaveAs(Form("Flow_%s_mass%d.png", labels[t].c_str(), int(s)));
        //cAll->Write(Form("Flow_%s_mass%d", labels[t].c_str(), int(s)));

        delete cAll; // optional cleanup
    }
}
file->Close();

    return 0;
}
