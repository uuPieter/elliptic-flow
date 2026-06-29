
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

    TFile *f = TFile::Open("/home/pieter/Downloads/Oxygen_train_669682.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return 1;
    }

    // create new Root file to save figures and th1d //
    TFile* file = TFile::Open("pieter_oxygen_train_669682_occu_base.root", "RECREATE");
    file->cd();

    // save between loop
    std::vector<int> histTypeIndex; // 0 = S+B, 1 = B, 2 = S
    std::vector<std::string> histLabels;
    std::vector<TH1D*> hFlowVec;
    std::vector<std::string> cenLabels = {"base_0", "base_10", "base_20", "base_30", "base_50"};
    std::vector<int> colors = {kRed+1, kBlue, kBlack, kBlue+1, kCyan+1, kGreen+2, kOrange+7};

    std::vector<std::vector<double>> weightedSum_v2S_2D;    // numerator: S*v2
    std::vector<std::vector<double>> sum_S_2D;              // denominator: sum(S)
    std::vector<std::vector<double>> weightedErr2_v2S_2D;   // sum((S*err)^2)
    
    bool firstCentrality = true; // flag to resize vectors

    // start the loop over the different centralities //
    for (int ip = 0; ip < 5; ++ip) {
    
    // enter the correct subfolder of root file //
    TDirectory* dir = file->mkdir(Form("centrality_%s", cenLabels[ip].c_str()));
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

     // get the multiplicity 
     if (!d_dielectron) { std::cerr << "dielectron not found!\n"; return 1; }
     TDirectory *d_Event = (TDirectory*) d_dielectron->Get("Event");
     if (!d_Event) { std::cerr << "Pair not found!\n"; return 1; }
     TDirectory *d_after = (TDirectory*) d_Event->Get("after");
     if (!d_after) { std::cerr << "Pair not found!\n"; return 1; }
     TH2D* hZvtx = (TH2D*) d_after->Get("hZvtx"); 
     if (!hZvtx) { std::cerr << "hZvtx histogram not found!\n"; return 1; }
     double multiplicity = hZvtx->GetEntries();
     std::cout << "multiplicity" << multiplicity << std::endl;

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
    double dcamin = 0.0;
    double dcamax = 40.0;
    double delta = 1e-4;
    
    int dcaLowBin  = hs_uls->GetAxis(2)->FindBin(dcamin + delta);
    int dcaHighBin = hs_uls->GetAxis(2)->FindBin(dcamax - delta);

    // Make a clone to not touch original
THnSparseD* hs_uls_slice  = (THnSparseD*) hs_uls->Clone("hs_uls_slice");
THnSparseD* hs_lspp_slice = (THnSparseD*) hs_lspp->Clone("hs_lspp_slice");
THnSparseD* hs_lsmm_slice = (THnSparseD*) hs_lsmm->Clone("hs_lsmm_slice");

// Apply DCA range
hs_uls_slice->GetAxis(2)->SetRange(dcaLowBin, dcaHighBin);
hs_lspp_slice->GetAxis(2)->SetRange(dcaLowBin, dcaHighBin);
hs_lsmm_slice->GetAxis(2)->SetRange(dcaLowBin, dcaHighBin);
TH3D* h_proj_uls3D = (TH3D*) hs_uls_slice->Projection(0 , 1 , 4);  
TH3D* h_proj_lspp3D = (TH3D*) hs_lspp_slice->Projection(0 , 1 , 4);  
TH3D* h_proj_lsmm3D = (TH3D*) hs_lsmm_slice->Projection(0 , 1 , 4); 

delete hs_uls_slice;   hs_uls_slice  = nullptr;
delete hs_lspp_slice;  hs_lspp_slice = nullptr;
delete hs_lsmm_slice;  hs_lsmm_slice = nullptr;

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
    
\
    // --- X-axis bins ---
    std::vector<double> xbins;
    for (auto &[start, end, step] : xRanges) {
        if (xbins.empty()) xbins.push_back(start); 
        for (double x = start + step; x < end - 1e-8; x += step) {
            xbins.push_back(x);
        }
        xbins.push_back(end); 
    }
    
    // --- Y-axis bins ---
    std::vector<double> ybins;
    for (auto &[start, end, step] : yRanges) {
        if (ybins.empty()) ybins.push_back(start);
        for (double y = start + step; y < end - 1e-8; y += step) {
            ybins.push_back(y);
        }
        ybins.push_back(end);
    }
    
    // --- Z-axis bins ---
    std::vector<double> zbins;
    for (auto &[start, end, step] : zRanges) {
        if (zbins.empty()) zbins.push_back(start);
        for (double z = start + step; z < end - 1e-8; z += step) {
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
    for (size_t i = 0; i < zbins.size() - 1; ++i) {
    
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

TH2D* h_rebinned_uls_mix_2D  = (TH2D*)h_rebinned_uls_mix_3D->Project3D("yx");
TH2D* h_rebinned_lspp_mix_2D = (TH2D*)h_rebinned_lspp_mix_3D->Project3D("yx");
TH2D* h_rebinned_lsmm_mix_2D = (TH2D*)h_rebinned_lsmm_mix_3D->Project3D("yx");

int nBinsXrebin = h_rebinned_uls_3D->GetNbinsX();
int nBinsYrebin = h_rebinned_uls_3D->GetNbinsY();
int nBinsZrebin = h_rebinned_uls_3D->GetNbinsZ();


TH2D *hZ_2D_rebinned = new TH2D("hZ_2D_rebinned", "Z;X;Y",xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data());
TH2D *hB_2D_rebinned = new TH2D("hB_2D_rebinned", "B;X;Y",xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data());
TH2D *hS_over_B_2D_rebinned = new TH2D("hSoverB_2D_rebinned", "S/B;X;Y",xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data());
TH2D *hR_2D_rebinned = new TH2D("hR_2D_rebinned", "R;X;Y", xbins.size()-1, xbins.data(), ybins.size()-1, ybins.data());
TH2D *hN_2D_rebinned = new TH2D("hN_2D_rebinned", "N;X;Y",xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data());
// Uncertainties
TH2D *hSigma_Z_2D_rebinned = new TH2D("hSigma_Z_2D_rebinned", "σZ;X;Y",xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data());
TH2D *hSigma_B_2D_rebinned = new TH2D("hSigma_B_2D_rebinned", "σB;X;Y", xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data());
TH2D *hSigma_R_2D_rebinned = new TH2D("hSigma_R_2D_rebinned", "σR;X;Y",xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data());
TH2D *hSigma_N_2D_rebinned = new TH2D("hSigma_N_2D_rebinned", "σN;X;Y",xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data());
TH2D *hSigma_S_over_B_2D_rebinned = new TH2D("hSigma_S_over_B_2D_rebinned", "σN;X;Y",xbins.size()-1, xbins.data(),ybins.size()-1, ybins.data());


// get the entries from each individual bin and calculate the variables we are interessted in 
for (int ix = 1; ix <= nBinsXrebin; ix++) {
    for (int iy = 1; iy <= nBinsYrebin; iy++) {
        double N_pm_same = h_rebinned_uls_2D->GetBinContent(ix, iy);
        double N_pp_same = h_rebinned_lspp_2D->GetBinContent(ix, iy);
        double N_mm_same = h_rebinned_lsmm_2D->GetBinContent(ix, iy);

        double N_pm_mix = h_rebinned_uls_mix_2D->GetBinContent(ix, iy);
        double N_pp_mix = h_rebinned_lspp_mix_2D->GetBinContent(ix, iy);
        double N_mm_mix = h_rebinned_lsmm_mix_2D->GetBinContent(ix, iy);

        if (N_pp_mix > 0 && N_mm_mix > 0) {
            double R = N_pm_mix / (2.0 * std::sqrt(N_pp_mix * N_mm_mix));
            double B = 2.0 * R * std::sqrt(N_pp_same * N_mm_same);
            double Z = N_pm_same - B;
            double N = N_pm_same;
            double S_over_B = (B > 0) ? Z / B : 0.0;
            hN_2D_rebinned     ->SetBinContent(ix, iy, N );
            hR_2D_rebinned     ->SetBinContent(ix, iy, R );
            hZ_2D_rebinned     ->SetBinContent(ix, iy, Z );
            hB_2D_rebinned     ->SetBinContent(ix, iy, B );
            hS_over_B_2D_rebinned->SetBinContent(ix, iy, S_over_B);
            double S2 = N_pp_mix * N_mm_mix;
            double S  = std::sqrt(S2);
            

            // R uncertainty
            double term1 = N_pm_mix / (4.0 * S2);
            double term2 = (R*R / 4.0) * ((N_pp_mix>0 ? 1.0/N_pp_mix : 0.0) + (N_mm_mix>0 ? 1.0/N_mm_mix : 0.0));
            double varR = term1 + term2;
            if (varR < 0 || std::isnan(varR) || std::isinf(varR)) varR = 0;
            hSigma_R_2D_rebinned->SetBinContent(ix, iy, std::sqrt(varR));

            // Background uncertainty
            double sigmaB2 = 4.0 * N_pp_same * N_mm_same * varR + (R*R) * (N_pp_same + N_mm_same);
            if (sigmaB2 < 0 || std::isnan(sigmaB2) || std::isinf(sigmaB2)) sigmaB2 = 0;
            hSigma_B_2D_rebinned->SetBinContent(ix, iy, std::sqrt(sigmaB2));

            // Z uncertainty
            double varZ = N_pm_same + sigmaB2;
            if (varZ < 0 || std::isnan(varZ) || std::isinf(varZ)) varZ = 0;
            hSigma_Z_2D_rebinned->SetBinContent(ix, iy, std::sqrt(varZ));

            // N uncertainty
            hSigma_N_2D_rebinned->SetBinContent(ix, iy, std::sqrt(N_pm_same));
                            // Signal over background uncertainty
            double sigmaSB = 0.0;

            if (Z > 0 && B > 0) {
                double sigmaZ = hSigma_Z_2D_rebinned->GetBinContent(ix, iy);
                double sigmaB = hSigma_B_2D_rebinned->GetBinContent(ix, iy);
            
                double termZ = sigmaZ / Z;
                double termB = sigmaB / B;
                sigmaSB = (Z/B) * std::sqrt(termZ*termZ + termB*termB);
            if (std::isnan(sigmaSB) || std::isinf(sigmaSB)) sigmaSB = 0.0;
            }
            hSigma_S_over_B_2D_rebinned->SetBinContent(ix, iy, sigmaSB);
        }
    }
} 
std::cout << "test" << std::endl;
// 2D mix projections no longer needed after the variable calculation loop above
delete h_rebinned_uls_2D;      h_rebinned_uls_2D     = nullptr;
delete h_rebinned_lspp_2D;     h_rebinned_lspp_2D    = nullptr;
delete h_rebinned_lsmm_2D;     h_rebinned_lsmm_2D    = nullptr;
delete h_rebinned_uls_mix_2D;  h_rebinned_uls_mix_2D  = nullptr;
delete h_rebinned_lspp_mix_2D; h_rebinned_lspp_mix_2D = nullptr;
delete h_rebinned_lsmm_mix_2D; h_rebinned_lsmm_mix_2D = nullptr;

hZ_2D_rebinned->Write();



double xMin_rebin = h_rebinned_uls_3D->GetXaxis()->GetXmin();
double xMax_rebin = h_rebinned_uls_3D->GetXaxis()->GetXmax();
double yMin_rebin = h_rebinned_uls_3D->GetYaxis()->GetXmin();
double yMax_rebin = h_rebinned_uls_3D->GetYaxis()->GetXmax();

 // ======================== calculate variables 1D Mee ========================//
 //       calculate uls, signal, background, S/B and R in 1D                    //
 // ======================== calculate variables 1D Mee ========================//

TH1D* h_rebinned_uls_1D_x  = h_rebinned_uls_3D->ProjectionX();  
TH1D* h_rebinned_lspp_1D_x = h_rebinned_lspp_3D->ProjectionX();
TH1D* h_rebinned_lsmm_1D_x = h_rebinned_lsmm_3D->ProjectionX();
TH1D* h_rebinned_uls_mix_1D_x  = h_rebinned_uls_mix_3D->ProjectionX();
TH1D* h_rebinned_lspp_mix_1D_x = h_rebinned_lspp_mix_3D->ProjectionX();
TH1D* h_rebinned_lsmm_mix_1D_x = h_rebinned_lsmm_mix_3D->ProjectionX();

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

            if (N_pp_mix > 0 && N_mm_mix > 0) {
                double R = N_pm_mix / (2.0 * std::sqrt(N_pp_mix * N_mm_mix));
                double B = 2.0 * R * std::sqrt(N_pp_same * N_mm_same);
                double Z = N_pm_same - B;
                double N = N_pm_same;
                double S_over_B = (B > 0) ? Z / B : 0.0;

                hN_1D_mee_rebin     ->SetBinContent(ix, N );
                hR_1D_mee_rebin     ->SetBinContent(ix, R );
                hZ_1D_mee_rebin     ->SetBinContent(ix, Z );
                hB_1D_mee_rebin     ->SetBinContent(ix, B );
                hSoverB_1D_mee_rebin->SetBinContent(ix, S_over_B);

                double S2 = N_pp_mix * N_mm_mix;
                double S  = std::sqrt(S2);
                
                // R uncertainty
                double term1 = N_pm_mix / (4.0 * S2);
                double term2 = (R*R / 4.0) * ((N_pp_mix>0 ? 1.0/N_pp_mix : 0.0) + (N_mm_mix>0 ? 1.0/N_mm_mix : 0.0));
                double varR = term1 + term2;
                if (varR < 0 || std::isnan(varR) || std::isinf(varR)) varR = 0;
                hR_1D_mee_rebin->SetBinError(ix, std::sqrt(varR));

                // Background uncertainty
                double sigmaB2 = 4.0 * N_pp_same * N_mm_same * varR + (R*R) * (N_pp_same + N_mm_same);
                if (sigmaB2 < 0 || std::isnan(sigmaB2) || std::isinf(sigmaB2)) sigmaB2 = 0;
                hB_1D_mee_rebin->SetBinError(ix, std::sqrt(sigmaB2));

                // Z uncertainty
                double varZ = N_pm_same + sigmaB2;
                if (varZ < 0 || std::isnan(varZ) || std::isinf(varZ)) varZ = 0;
                hZ_1D_mee_rebin->SetBinError(ix, std::sqrt(varZ));

                // N uncertainty
                hN_1D_mee_rebin->SetBinError(ix,  std::sqrt(N_pm_same));
                
                // Signal over background uncertainty
                double sigmaSB = 0.0;

                if (Z > 0 && B > 0) {
                    double sigmaZ = std::sqrt(varZ);
                    double sigmaB = std::sqrt(sigmaB2);
            
                    double termZ = sigmaZ / Z;
                    double termB = sigmaB / B;
                    sigmaSB = (Z/B) * std::sqrt(termZ*termZ + termB*termB);
                if (std::isnan(sigmaSB) || std::isinf(sigmaSB)) sigmaSB = 0.0;
                }
                hSoverB_1D_mee_rebin->SetBinError(ix, sigmaSB);
            
        }
    }
// 1D mee projections no longer needed after the loop above
delete h_rebinned_uls_1D_x;      h_rebinned_uls_1D_x     = nullptr;
delete h_rebinned_lspp_1D_x;     h_rebinned_lspp_1D_x    = nullptr;
delete h_rebinned_lsmm_1D_x;     h_rebinned_lsmm_1D_x    = nullptr;
delete h_rebinned_uls_mix_1D_x;  h_rebinned_uls_mix_1D_x  = nullptr;
delete h_rebinned_lspp_mix_1D_x; h_rebinned_lspp_mix_1D_x = nullptr;
delete h_rebinned_lsmm_mix_1D_x; h_rebinned_lsmm_mix_1D_x = nullptr;
 
// 1D mee histograms written, no longer needed
delete hZ_1D_mee_rebin;    hZ_1D_mee_rebin    = nullptr;
delete hB_1D_mee_rebin;    hB_1D_mee_rebin    = nullptr;
delete hSoverB_1D_mee_rebin; hSoverB_1D_mee_rebin = nullptr;
delete hN_1D_mee_rebin;    hN_1D_mee_rebin    = nullptr;
 
hR_1D_mee_rebin->Write();

 
 // ======================== calculate variables 1D pT ========================//
 //       calculate uls, signal, background, S/B and R in 1D                   //
 // ======================== calculate variables 1D pT ========================//

TH1D* h_rebinned_uls_1D_pt  = h_rebinned_uls_3D->ProjectionY();
TH1D* h_rebinned_lspp_1D_pt = h_rebinned_lspp_3D->ProjectionY();
TH1D* h_rebinned_lsmm_1D_pt = h_rebinned_lsmm_3D->ProjectionY();

TH1D* h_rebinned_uls_mix_1D_pt  = h_rebinned_uls_mix_3D->ProjectionY();
TH1D* h_rebinned_lspp_mix_1D_pt = h_rebinned_lspp_mix_3D->ProjectionY();
TH1D* h_rebinned_lsmm_mix_1D_pt = h_rebinned_lsmm_mix_3D->ProjectionY();


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

    if (N_pp_mix > 0 && N_mm_mix > 0) {

        double R = N_pm_mix / (2.0 * std::sqrt(N_pp_mix * N_mm_mix));
        double B = 2.0 * R * std::sqrt(N_pp_same * N_mm_same);
        double Z = N_pm_same - B;
        double N = N_pm_same;
        double S_over_B = (B > 0.0) ? Z / B : 0.0;

        hN_1D_pt_rebin_initial     ->SetBinContent(iy, N);
        hR_1D_pt_rebin_initial     ->SetBinContent(iy, R);
        hZ_1D_pt_rebin_initial     ->SetBinContent(iy, Z);
        hB_1D_pt_rebin_initial     ->SetBinContent(iy, B);
        hSoverB_1D_pt_rebin_initial->SetBinContent(iy, S_over_B);

        double S2 = N_pp_mix * N_mm_mix;
        double S  = std::sqrt(S2);
        
        // R uncertainty
        double term1 = N_pm_mix / (4.0 * S2);
        double term2 = (R*R / 4.0) * ((N_pp_mix>0 ? 1.0/N_pp_mix : 0.0) + (N_mm_mix>0 ? 1.0/N_mm_mix : 0.0));
        double varR = term1 + term2;
        if (varR < 0 || std::isnan(varR) || std::isinf(varR)) varR = 0;
        hR_1D_pt_rebin_initial->SetBinError(iy, std::sqrt(varR));

        // Background uncertainty
        double sigmaB2 = 4.0 * N_pp_same * N_mm_same * varR + (R*R) * (N_pp_same + N_mm_same);
        if (sigmaB2 < 0 || std::isnan(sigmaB2) || std::isinf(sigmaB2)) sigmaB2 = 0;
        hB_1D_pt_rebin_initial->SetBinError(iy, std::sqrt(sigmaB2));

        // Z uncertainty
        double varZ = N_pm_same + sigmaB2;
        if (varZ < 0 || std::isnan(varZ) || std::isinf(varZ)) varZ = 0;
        hZ_1D_pt_rebin_initial->SetBinError(iy, std::sqrt(varZ));

        // N uncertainty
        hN_1D_pt_rebin_initial->SetBinError(iy,  std::sqrt(N_pm_same));
        
        // Signal over background uncertainty
        double sigmaSB = 0.0;

        if (Z > 0 && B > 0) {
            double sigmaZ = std::sqrt(varZ);
            double sigmaB = std::sqrt(sigmaB2);
    
            double termZ = sigmaZ / Z;
            double termB = sigmaB / B;
            sigmaSB = (Z/B) * std::sqrt(termZ*termZ + termB*termB);
        if (std::isnan(sigmaSB) || std::isinf(sigmaSB)) sigmaSB = 0.0;
        }
        hSoverB_1D_pt_rebin_initial->SetBinError(iy, sigmaSB);
    }
}

// 1D pT projections no longer needed after the loop above
delete h_rebinned_uls_1D_pt;      h_rebinned_uls_1D_pt     = nullptr;
delete h_rebinned_lspp_1D_pt;     h_rebinned_lspp_1D_pt    = nullptr;
delete h_rebinned_lsmm_1D_pt;     h_rebinned_lsmm_1D_pt    = nullptr;
delete h_rebinned_uls_mix_1D_pt;  h_rebinned_uls_mix_1D_pt  = nullptr;
delete h_rebinned_lspp_mix_1D_pt; h_rebinned_lspp_mix_1D_pt = nullptr;
delete h_rebinned_lsmm_mix_1D_pt; h_rebinned_lsmm_mix_1D_pt = nullptr;
 
// Write and delete 1D pT histograms
hR_1D_pt_rebin_initial->Write("rfactor_pt");
hR_1D_mee_rebin->Write("rfactor_mee");
 
delete hZ_1D_pt_rebin_initial;      hZ_1D_pt_rebin_initial     = nullptr;
delete hB_1D_pt_rebin_initial;      hB_1D_pt_rebin_initial     = nullptr;
delete hSoverB_1D_pt_rebin_initial; hSoverB_1D_pt_rebin_initial = nullptr;
delete hR_1D_pt_rebin_initial;      hR_1D_pt_rebin_initial     = nullptr;
delete hN_1D_pt_rebin_initial;      hN_1D_pt_rebin_initial     = nullptr;
std::cout << "test604" << std::endl;
// ======================================== Rebinning again ==================================== //
// here we rebin the already projected 1D histograms for better visualization and then plot them //
// ======================================= Rebinning again ===================================== // 





std::vector<BinRange> xRanges_visual = {
    {0.0, 0.44, 0.04},   
    {0.44, 1.2, 0.02},   
    {1.2, 2.8, 0.1},    
    {2.8, 3.2, 0.02},    
    {3.2, 4.0, 0.1}      
};
std::vector<BinRange> yRanges_visual = {
    {0.0, 5.0, 0.5},   
    {5.0, 10.0, 1.0}  
};

// --- X-axis bins ---
std::vector<double> xbins_visual;
for (auto &[start, end, step] : xRanges_visual) {
    if (xbins_visual.empty()) xbins_visual.push_back(start); 
    for (double x = start + step; x < end - 1e-8; x += step) {
        xbins_visual.push_back(x);
    }
    xbins_visual.push_back(end); 
}

// --- Y-axis bins ---
std::vector<double> ybins_visual;
for (auto &[start, end, step] : yRanges_visual) {
    if (ybins_visual.empty()) ybins_visual.push_back(start);
    for (double y = start + step; y < end - 1e-8; y += step) {
        ybins_visual.push_back(y);
    }
    ybins_visual.push_back(end);
}

TAxis* xax_rebin = hZ_2D_rebinned->GetXaxis();
TAxis* yax_rebin = hZ_2D_rebinned->GetYaxis();

std::vector<std::pair<int,int>> meeSlices2_visual;
for (size_t i = 0; i < xbins_visual.size() - 1; ++i) {

    int lo = xax_rebin->FindBin(xbins_visual[i] + 1e-6);
    int hi = xax_rebin->FindBin(xbins_visual[i+1] - 1e-6);

    lo = std::max(lo, 1);
    hi = std::min(hi, xax->GetNbins());

    meeSlices2_visual.emplace_back(lo, hi);
}

std::vector<std::pair<int,int>> ptSlices2_visual;
for (size_t i = 0; i < ybins_visual.size() - 1; ++i) {

    int lo = yax_rebin->FindBin(ybins_visual[i] + 1e-6);
    int hi = yax_rebin->FindBin(ybins_visual[i+1] - 1e-6);

    lo = std::max(lo, 1);
    hi = std::min(hi, yax->GetNbins());

    ptSlices2_visual.emplace_back(lo, hi);
}



TH2D *hZ_2D_rebinned_vis  = new TH2D("hZ_2D_rebinned_vis","Z;X;Y",xbins_visual.size()-1,xbins_visual.data(),ybins_visual.size()-1,ybins_visual.data());
TH2D *hB_2D_rebinned_vis  = new TH2D("hB_2D_rebinned_vis","B;X;Y",xbins_visual.size()-1,xbins_visual.data(),ybins_visual.size()-1,ybins_visual.data());
TH2D *hN_2D_rebinned_vis  = new TH2D("hN_2D_rebinned_vis","N;X;Y",xbins_visual.size()-1,xbins_visual.data(),ybins_visual.size()-1,ybins_visual.data());
TH2D *hSoverB_2D_rebinned_vis  = new TH2D("hSB_2D_rebinned_vis","N;X;Y",xbins_visual.size()-1,xbins_visual.data(),ybins_visual.size()-1,ybins_visual.data());
TH2D *hSigma_Z_2D_rebinned_vis = new TH2D("hSigma_Z_2D_rebinned_vis","σZ;X;Y",xbins_visual.size()-1,xbins_visual.data(),ybins_visual.size()-1,ybins_visual.data());
TH2D *hSigma_B_2D_rebinned_vis = new TH2D("hSigma_B_2D_rebinned_vis","σB;X;Y",xbins_visual.size()-1,xbins_visual.data(),ybins_visual.size()-1,ybins_visual.data());
TH2D *hSigma_N_2D_rebinned_vis = new TH2D("hSigma_N_2D_rebinned_vis","σN;X;Y",xbins_visual.size()-1,xbins_visual.data(),ybins_visual.size()-1,ybins_visual.data());
TH2D *hSigma_SoverB_2D_rebinned_vis = new TH2D("hSigma_SB_2D_rebinned_vis","σN;X;Y",xbins_visual.size()-1,xbins_visual.data(),ybins_visual.size()-1,ybins_visual.data());

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

                double eZ = hSigma_Z_2D_rebinned->GetBinContent(ix,iy);
                double eB = hSigma_B_2D_rebinned->GetBinContent(ix,iy);
                double eN = hSigma_N_2D_rebinned->GetBinContent(ix,iy);

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
hZ_2D_rebinned_vis->Write("hZ_2D_rebinned_vis");

for(int ix=1; ix<=hZ_2D_rebinned_vis->GetNbinsX(); ++ix){
    for(int iy=1; iy<=hZ_2D_rebinned_vis->GetNbinsY(); ++iy){

        double Z = hZ_2D_rebinned_vis->GetBinContent(ix,iy);
        double B = hB_2D_rebinned_vis->GetBinContent(ix,iy);

        double eZ = hSigma_Z_2D_rebinned_vis->GetBinContent(ix,iy);
        double eB = hSigma_B_2D_rebinned_vis->GetBinContent(ix,iy);

        if(B>0){
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

// Fine 2D rebinned histograms no longer needed after visual rebinning loop above
delete hS_over_B_2D_rebinned;       hS_over_B_2D_rebinned      = nullptr;
delete hSigma_S_over_B_2D_rebinned; hSigma_S_over_B_2D_rebinned = nullptr;



 // ========================== projection from 2D to 1D ==========================//
 //                   here we project the 2D histograms to 1D                     //
 // ========================== projection from 2D to 1D ==========================//

int nBins_rebin_X = hZ_2D_rebinned_vis->GetNbinsX();
int nBins_rebin_Y = hZ_2D_rebinned_vis->GetNbinsY();
std::cout << "test773" << std::endl;
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
std::cout << "test913" << std::endl;
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
 // Visual 2D rebinned histograms no longer needed after all 1D projection loops above
delete hZ_2D_rebinned_vis;          hZ_2D_rebinned_vis          = nullptr;
delete hB_2D_rebinned_vis;          hB_2D_rebinned_vis          = nullptr;
delete hN_2D_rebinned_vis;          hN_2D_rebinned_vis          = nullptr;
delete hSoverB_2D_rebinned_vis;     hSoverB_2D_rebinned_vis     = nullptr;
delete hSigma_Z_2D_rebinned_vis;    hSigma_Z_2D_rebinned_vis    = nullptr;
delete hSigma_B_2D_rebinned_vis;    hSigma_B_2D_rebinned_vis    = nullptr;
delete hSigma_N_2D_rebinned_vis;    hSigma_N_2D_rebinned_vis    = nullptr;
delete hSigma_SoverB_2D_rebinned_vis; hSigma_SoverB_2D_rebinned_vis = nullptr;

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

delete hN;           hN          = nullptr;
delete hB;           hB          = nullptr;
delete hZ;           hZ          = nullptr;
delete hS_over_B;    hS_over_B   = nullptr;
delete hsignificance; hsignificance = nullptr;

delete hN_pt;           hN_pt          = nullptr;
delete hB_pt;           hB_pt          = nullptr;
delete hZ_pt;           hZ_pt          = nullptr;
delete hS_over_B_pt;    hS_over_B_pt   = nullptr;
delete hSignificance_pt; hSignificance_pt = nullptr;
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

TLatex latex;
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
gerrS_over_B_rebinned_pt->GetXaxis()->SetRangeUser(0,10);     
gerrS_over_B_rebinned_pt->GetYaxis()->SetRangeUser(1e-3, 1e-1);

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
gsignificance_err_pt->GetYaxis()->SetRangeUser(1, 5e3);
gsignificance_err_pt->GetXaxis()->SetRangeUser(0.0, 10.0);
latex.SetNDC();
latex.SetTextSize(0.03);
latex.SetTextAlign(13);
latex.DrawLatex(0.15, 0.88, "Pb-Pb #sqrt{s_{NN}} = 5.36 TeV");
latex.DrawLatex(0.15, 0.85, "0_90");

csignificance_pt_spectra->SaveAs(("pt_spectra_significance"+ cenLabels[ip]  + ".png").c_str());

// All mee/pt spectra TGraphErrors are now done with
delete gNerr_tot;               gNerr_tot              = nullptr;
delete gBerr_tot;               gBerr_tot              = nullptr;
delete gSoverBerr_tot;          gSoverBerr_tot         = nullptr;
delete gsignificance_err_mee;   gsignificance_err_mee  = nullptr;
delete gZerr_rebinned;          gZerr_rebinned         = nullptr;
delete gerrS_over_B_rebinned_pt; gerrS_over_B_rebinned_pt = nullptr;
delete gZerr_rebinned_pt;       gZerr_rebinned_pt      = nullptr;
delete gNerr_rebin_pt;          gNerr_rebin_pt         = nullptr;
delete gsignificance_err_pt;    gsignificance_err_pt   = nullptr;
delete gBerr_rebin_pt;          gBerr_rebin_pt         = nullptr;

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
hCent->SaveAs(Form("hCentFT0CvsMultNTracksPV%s.png", cenLabels[ip].c_str()));
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

    
    slices.push_back(h1);
}


for (size_t i = 0; i < means_squared.size(); ++i) {
    std::cout << "mean[" << i << "] = " << means_squared[i] << std::endl;
}




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
delete hRn; hRn = nullptr;
 
// Gaussian fit slices no longer needed
for (auto* h : slices) delete h;
slices.clear();
std::cout << "test1596" << std::endl;
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
std::cout << "test1621" << std::endl;
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
delete h_proj_uls3D;  h_proj_uls3D  = nullptr;
delete h_proj_lspp3D; h_proj_lspp3D = nullptr;
delete h_proj_lsmm3D; h_proj_lsmm3D = nullptr;
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

        double sN = hSigma_N_2D_rebinned->GetBinContent(ix, iy);
        double sZ = hSigma_Z_2D_rebinned->GetBinContent(ix, iy);
        double sB = hSigma_B_2D_rebinned->GetBinContent(ix, iy);
        double sR = hSigma_R_2D_rebinned->GetBinContent(ix, iy);

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
std::cout << "test1743" << std::endl;
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
std::cout << "test1776" << std::endl;
delete hN_2D_rebinned;       hN_2D_rebinned      = nullptr;
delete hZ_2D_rebinned;       hZ_2D_rebinned      = nullptr;
delete hB_2D_rebinned;       hB_2D_rebinned      = nullptr;
delete hR_2D_rebinned;       hR_2D_rebinned      = nullptr;
delete hSigma_N_2D_rebinned; hSigma_N_2D_rebinned = nullptr;
delete hSigma_Z_2D_rebinned; hSigma_Z_2D_rebinned = nullptr;
delete hSigma_B_2D_rebinned; hSigma_B_2D_rebinned = nullptr;
delete hR_2D_flow;              hR_2D_flow             = nullptr;
delete hS_over_B_2D_flow;       hS_over_B_2D_flow      = nullptr;
delete hSignificance_2D_flow;   hSignificance_2D_flow   = nullptr;
delete hSigma_R_2D_flow;        hSigma_R_2D_flow       = nullptr;
delete hSigma_S_over_B_2D_flow; hSigma_S_over_B_2D_flow = nullptr;
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


std::vector<double> ybins_flow = {
    0.0, 0.3, 0.5, 0.7,
    1.0,
    1.5,  2.0, 2.5, 3.5, 
    5.0,
    10.0
};
std::cout << "test1809" << std::endl;
TH3D* h_rebinned_uls_flow_visual = new TH3D("h_uls_flow_visual", "ULS (flow bins);m_{ee};p_{T};#Delta#phi",xbins_flow.size()-1, xbins_flow.data(),ybins_flow.size()-1, ybins_flow.data(),zbins_flow.size()-1, zbins_flow.data());
TH3D* h_rebinned_lspp_flow_visual = new TH3D("h_lspp_flow_visual", "LS++ (flow bins);m_{ee};p_{T};#Delta#phi",xbins_flow.size()-1, xbins_flow.data(),ybins_flow.size()-1, ybins_flow.data(),zbins_flow.size()-1, zbins_flow.data());
TH3D* h_rebinned_lsmm_flow_visual = new TH3D("h_lsmm_flow_visual", "LS-- (flow bins);m_{ee};p_{T};#Delta#phi",xbins_flow.size()-1, xbins_flow.data(),ybins_flow.size()-1, ybins_flow.data(),zbins_flow.size()-1, zbins_flow.data());
TH2D* hN_2D_flow_visual = new TH2D("hN_2D_flow_visual", "N (flow bins);m_{ee};p_{T}",xbins_flow.size()-1, xbins_flow.data(),ybins_flow.size()-1, ybins_flow.data());
TH2D* hZ_2D_flow_visual = new TH2D("hZ_2D_flow_visual", "Z (flow bins);m_{ee};p_{T}",xbins_flow.size()-1, xbins_flow.data(),ybins_flow.size()-1, ybins_flow.data());
TH2D* hB_2D_flow_visual = new TH2D("hB_2D_flow_visual", "B (flow bins);m_{ee};p_{T}",xbins_flow.size()-1, xbins_flow.data(),ybins_flow.size()-1, ybins_flow.data());
TH2D* hSigma_Z_2D_flow_visual = new TH2D("hSigma_Z_2D_flow_visual", "σZ (flow bins);m_{ee};p_{T}",xbins_flow.size()-1, xbins_flow.data(),ybins_flow.size()-1, ybins_flow.data());
TH2D* hSigma_B_2D_flow_visual = new TH2D("hSigma_B_2D_flow_visual", "σB (flow bins);m_{ee};p_{T}",xbins_flow.size()-1, xbins_flow.data(),ybins_flow.size()-1, ybins_flow.data());
TH2D* hSigma_N_2D_flow_visual = new TH2D("hSigma_N_2D_flow_visual", "σN (flow bins);m_{ee};p_{T}",xbins_flow.size()-1, xbins_flow.data(),ybins_flow.size()-1, ybins_flow.data());
std::vector<std::pair<int,int>> meeSlices;
std::cout << "test" << std::endl;
for (size_t i = 0; i < xbins_flow.size() - 1; ++i) {
    std::cout << "test1822" << std::endl;
    int lo = FindBinLow (xax, xbins_flow[i]);
    int hi = FindBinHigh(xax, xbins_flow[i+1]);
    std::cout << "test1825" << std::endl;
    lo = std::max(lo, 1);
    hi = std::min(hi, xax->GetNbins());

    meeSlices.emplace_back(lo, hi);
}
std::cout << "test1830" << std::endl;
std::vector<std::pair<int,int>> ptSlices;

for (size_t i = 0; i < ybins_flow.size() - 1; ++i) {

    int lo = FindBinLow (yax, ybins_flow[i]);
    int hi = FindBinHigh(yax, ybins_flow[i+1]);

    lo = std::max(lo, 1);
    hi = std::min(hi, yax->GetNbins());

    ptSlices.emplace_back(lo, hi);
}
std::cout << "test1842" << std::endl;
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
delete h_rebinned_uls_flow;  h_rebinned_uls_flow  = nullptr;
delete h_rebinned_lspp_flow; h_rebinned_lspp_flow = nullptr;
delete h_rebinned_lsmm_flow; h_rebinned_lsmm_flow = nullptr;
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
std::cout << "test1919" << std::endl;
delete hN_2D_flow;       hN_2D_flow      = nullptr;
delete hZ_2D_flow;       hZ_2D_flow      = nullptr;
delete hB_2D_flow;       hB_2D_flow      = nullptr;
delete hSigma_N_2D_flow; hSigma_N_2D_flow = nullptr;
delete hSigma_Z_2D_flow; hSigma_Z_2D_flow = nullptr;
delete hSigma_B_2D_flow; hSigma_B_2D_flow = nullptr;



int nbz2 = h_rebinned_lsmm_flow_visual->GetNbinsZ();
int nBinsYvariable = h_rebinned_lsmm_flow_visual->GetNbinsY();
std::vector<double> yedges2(nBinsYvariable+1);
for (int i = 0; i <= nBinsYvariable; ++i)
yedges2[i] = h_rebinned_uls_flow_visual->GetYaxis()->GetBinLowEdge(i+1);
// Loop over slices (each slice was previously a single X-bin range)
std::vector<int> slices_to_use = {1, 2, 4}; // slices you want
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
std::cout << "test2074" << std::endl;

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
gSignificance_err_flow_pt->GetYaxis()->SetRangeUser(1e-2, 1e4);
gSignificance_err_flow_pt->GetXaxis()->SetRangeUser(0.0, 10.0);
latex.SetNDC();
latex.SetTextSize(0.03);
latex.SetTextAlign(13);
latex.DrawLatex(0.15, 0.88, "Pb-Pb #sqrt{s_{NN}} = 5.36 TeV");
latex.DrawLatex(0.15, 0.85, cenLabels[ip].c_str());

// Save
csignificance_pt_flow_spectra->SaveAs(("pt_spectra_flow_significance_" + cenLabels[ip]  + ".png").c_str());
delete csignificance_pt_flow_spectra; csignificance_pt_flow_spectra = nullptr;
 
// All flow pT overview TGraphErrors done
delete g_pt_per_slice;          g_pt_per_slice         = nullptr;
delete gNerr_flow_pt;           gNerr_flow_pt          = nullptr;
delete gBerr_flow_pt;           gBerr_flow_pt          = nullptr;
delete gS_over_B_err_flow_pt;   gS_over_B_err_flow_pt  = nullptr;
delete gSignificance_err_flow_pt; gSignificance_err_flow_pt = nullptr;

/*for (int ix = 1; ix <= hZ_1D_pt_flow->GetNbinsX(); ++ix) {
    std::cout << ix << "  "
              << hZ_1D_pt_flow->GetBinContent(ix) << " "
              << hB_1D_pt_flow->GetBinContent(ix) << " "
              << hSigma_Z_1D_pt_flow->GetBinContent(ix) << " "
              << hSigma_B_1D_pt_flow->GetBinContent(ix) <<std::endl;
        
} */
if (firstCentrality) {
    int nSlices = slices_to_use.size();   // e.g. 3
    int nBinsPt = nBinsYvariable;         // pt bins

    weightedSum_v2S_2D.assign(nSlices, std::vector<double>(nBinsPt, 0.0));
    sum_S_2D.assign(nSlices, std::vector<double>(nBinsPt, 0.0));
    weightedErr2_v2S_2D.assign(nSlices, std::vector<double>(nBinsPt, 0.0));

    firstCentrality = false;
}

int nMass = hZ_2D_flow_visual->GetNbinsX();;  // number of mass ranges (for example)
int nBinsPt = hZ_2D_flow_visual->GetNbinsY();


for (size_t s = 0; s < slices_to_use.size(); ++s) {
    int ix = slices_to_use[s]; // <- use the value from slices_to_use
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



            /*std::cout << "Slice=" << s << " iy=" << iy
                      << "  v2(S+B)=" << v2Pt[iy-1]
                      << "  v2(B)=" << v2BPt[iy-1]
                      << "  v2(S)=" << v2SPt[iy-1] << std::endl; */
        }

        hZ_2D_flow_visual->Write();
        double S_val   = hZ_2D_flow_visual->GetBinContent(ix, iy);  // signal for this mass & pt
        double v2S_val = v2SPt[iy-1];
        double v2S_err = v2SPtErr[iy-1];

        weightedSum_v2S_2D[s][iy-1]   += S_val * v2S_val;
        sum_S_2D[s][iy-1]            += S_val;
        weightedErr2_v2S_2D[s][iy-1] += (S_val * v2S_err) * (S_val * v2S_err);
        std::cout << "Slice=" << s << " iy=" << iy
                  << "  weigted_sum =" <<     weightedSum_v2S_2D[s][iy-1] 
                  << "  weight =" <<         sum_S_2D[s][iy-1]
                  << std::endl; 

    }
    std::cout << "" << s << std::endl;

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
            h->SetMaximum(0.4);
            h->SetMinimum(-0.1);
            h->GetXaxis()->SetRangeUser(0,5);
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
            latex.DrawLatex(0.18, 0.84, ("mass region " + massrange2[s] + " GeV/c^{2} ").c_str());
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
        
            // pT bin centers / half-widths
            double low = h_rebinned_lsmm_flow_visual->GetYaxis()->GetBinLowEdge(iy+1);
            double high = h_rebinned_lsmm_flow_visual->GetYaxis()->GetBinUpEdge(iy+1);
            y_rebinned_slice[iy] = 0.5 * (low + high);
            ey_rebinned_slice[iy] = 0.5 * (high - low);
        }
        
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
        c_slice->SaveAs(("plot_standard_large_cuts_cuts_testing_pt_slice_" + cenLabels[ip] + std::to_string(ix) + ".png").c_str());
        delete c_slice;           c_slice           = nullptr;
        delete gZerr_slice;       gZerr_slice       = nullptr;
        delete gNerr_flow_pt_slice; gNerr_flow_pt_slice = nullptr;
        delete gBerr_flow_pt_slice; gBerr_flow_pt_slice = nullptr;

        /*for (int ix = 1; ix <= 11; ++ix) {
            std::cout << ix << "  "
                      << Z_rebinned_pt_slice[ix-1] << " "
                      << y_rebinned_slice [ix-1] << std::endl;
                
        }*/
        
}

for (size_t s = 0; s < slices_to_use.size(); ++s) {
    int ix = slices_to_use[s]; // the X-bin index
    double xcenter = h_rebinned_uls_flow_visual->GetXaxis()->GetBinCenter(ix);

   /* std::cout << "Slice index=" << s
              << "  X-bin=" << ix
              << "  X-center=" << xcenter
              << std::endl; */
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

        int nMass = weightedSum_v2S_2D.size();
        int nBinsPt = (nMass > 0) ? weightedSum_v2S_2D[0].size() : 0;
        std::vector<std::vector<double>> v2S_weighted_2D(nMass, std::vector<double>(nBinsPt, 0.0));
        std::vector<std::vector<double>> v2S_weighted_err_2D(nMass, std::vector<double>(nBinsPt, 0.0));
        
        for (int mass = 0; mass < nMass; ++mass) {
            for (int pt = 0; pt < nBinsPt; ++pt) {
                if (sum_S_2D[mass][pt] != 0.0) {
                    v2S_weighted_2D[mass][pt]     = weightedSum_v2S_2D[mass][pt] / sum_S_2D[mass][pt];
                    v2S_weighted_err_2D[mass][pt] = std::sqrt(weightedErr2_v2S_2D[mass][pt]) / sum_S_2D[mass][pt];
                }
            }
        }

        


                

        std::vector<double> ybins_flow = {
            0.0, 0.3, 0.5, 0.7,
            1.0,
            1.5,  2.0, 2.5, 3.5, 
            5.0,
            10.0
        };

            

        
        std::vector<std::string> massrange2 = {"0 < m_{ee} < 0.14", "0.14 < m_{ee} < 1.2", "2.7 < m_{ee} < 3.2"};

        for (int mass = 0; mass < nMass; ++mass) {
            std::string name = "plot_of_v2S_massrange_" + std::to_string(mass);
        
            TH1D* h = new TH1D(
                name.c_str(),
                ";p_{T} [GeV/c]; v_{2}",
                nBinsPt,
                ybins_flow.data()  // your pt bin edges
            );
            
            h->SetDirectory(0); // detach from ROOT memory management
        
            // Fill histogram
            for (int pt = 0; pt < nBinsPt; ++pt) {
                h->SetBinContent(pt+1, v2S_weighted_2D[mass][pt]);
                h->SetBinError(pt+1, v2S_weighted_err_2D[mass][pt]);
            }

            // Canvas
            TCanvas* c = new TCanvas(("c"+name).c_str(), name.c_str(), 800, 800);
            c->SetTicks(1,1);
            c->SetLeftMargin(0.16);
            c->SetBottomMargin(0.12);
            c->SetRightMargin(0.04);
            c->SetTopMargin(0.05);
        
            // Styling
            h->SetMarkerSize(1.5);
            h->SetLineColor(kBlack);
            h->SetMarkerColor(kBlack);
            h->SetMarkerStyle(24);      
            h->SetLineWidth(2);
            h->SetStats(0);
            h->SetMaximum(0.4);
            h->SetMinimum(-0.1);
        
            h->Draw("E1");
        
            TLine* l0 = new TLine(0,0,10,0);
            l0->SetLineStyle(3);
            l0->SetLineWidth(2);
            l0->SetLineColor(kGray+2);
            l0->Draw("same");

            latex.DrawLatex(0.18, 0.84, ("mass region 0 < m_{ee} < 0.14 GeV/c^{2} "));
            latex.DrawLatex(0.18, 0.89, ("10_30 Pb-Pb #sqrt{s_{NN}} = 5.36 TeV"));
        

            h->Write();
            c->SaveAs((name + ".png").c_str());
            delete c; c = nullptr;
            delete h; h = nullptr;

        

            
        }



    }
}
// Free hFlowVec entries (kept alive for overlay plots, now done)
for (auto* h : hFlowVec) delete h;
hFlowVec.clear();
file->Close();

    return 0;
}
