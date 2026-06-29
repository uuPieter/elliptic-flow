#include <iostream>
#include <vector>
#include <cmath>

#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include <TFile.h>

using namespace std;


////////////////////////////////////////////////////////////
// Add variation helper
////////////////////////////////////////////////////////////

void AddVariation(vector<TH1D*>& variations, TH1D* hVar)
{
    if (hVar)
        variations.push_back(hVar);
}


////////////////////////////////////////////////////////////
// Compute max deviation systematic
////////////////////////////////////////////////////////////

TH1D* ComputeALICESystematic(TH1D* hNominal, vector<TH1D*> variations, const char* name)
{
    int nBins = hNominal->GetNbinsX();

    TH1D* hSyst = (TH1D*)hNominal->Clone(name);
    hSyst->Reset();

    for (int i = 1; i <= nBins; i++)
    {
        double nominal = hNominal->GetBinContent(i);

        vector<double> deltas;

        // collect all trial deviations
        for (auto& hVar : variations)
        {
            double var      = hVar->GetBinContent(i);
        
            double sigmaNom = hNominal->GetBinError(i);
            double sigmaVar = hVar->GetBinError(i);
        
            double delta = var - nominal;
        
            // Barlow uncertainty
            double sigmaDelta2 =
                fabs(sigmaVar*sigmaVar - sigmaNom*sigmaNom);
        
            double sigmaDelta =
                (sigmaDelta2 > 0)
                ? sqrt(sigmaDelta2)
                : 0.0;
        
            double nSigma = 0.0;
        
            if (sigmaDelta > 0)
                nSigma = fabs(delta) / sigmaDelta;
        
            // Apply Barlow criterion
            if (nSigma > 1.0)
                deltas.push_back(delta);
            else
                deltas.push_back(0.0);
        }

        if (deltas.empty())
        {
            hSyst->SetBinContent(i, 0);
            continue;
        }

        // mean shift
        double mean = 0;
        for (double d : deltas) mean += d;
        mean /= deltas.size();

        // RMS
        double rms = 0;
        for (double d : deltas)
            rms += (d - mean) * (d - mean);
        rms = sqrt(rms / deltas.size());

        // ALICE combination
        double syst = sqrt(mean * mean + rms * rms);

        hSyst->SetBinContent(i, syst);
    }

    return hSyst;
}


////////////////////////////////////////////////////////////
// Combine systematic sources in quadrature
////////////////////////////////////////////////////////////

TH1D* CombineSystematics(vector<TH1D*> systSources, const char* name)
{
    int nBins = systSources[0]->GetNbinsX();

    TH1D* hTotal = (TH1D*)systSources[0]->Clone(name);
    hTotal->Reset();

    for (int i = 1; i <= nBins; i++)
    {
        double sum = 0;

        for (auto& h : systSources)
        {
            double val = h->GetBinContent(i);
            sum += val * val;
        }

        hTotal->SetBinContent(i, sqrt(sum));
    }

    return hTotal;
}

 
////////////////////////////////////////////////////////////
// Convert histogram to systematic band
////////////////////////////////////////////////////////////

TGraphErrors* MakeSystematicBand(TH1D* hNominal, TH1D* hSyst)
{
    int nBins = hNominal->GetNbinsX();

    TGraphErrors* g = new TGraphErrors(nBins);

    for (int i = 1; i <= nBins; i++)
    {
        double x  = hNominal->GetBinCenter(i);
        double y  = hNominal->GetBinContent(i);
        double ex = hNominal->GetBinWidth(i) / 8.0;
        double ey = hSyst->GetBinContent(i);

        g->SetPoint(i-1, x, y);
        g->SetPointError(i-1, ex, ey);
    }

    g->SetFillColorAlpha(kGray+1,0.35);
    g->SetLineWidth(0);

    return g;
}


////////////////////////////////////////////////////////////
// MAIN
////////////////////////////////////////////////////////////

int main()
{   
    TFile* fOut = new TFile("systematics_output.root", "RECREATE");
    vector<int> centBins =
    {0,10,20,30};
    std::vector<int> markerStyles = {
        20, // open circle
        21, // open square
        24, // open triangle up
        25  // open star
    };
    

    vector<TH1D*> centralityTotals;
    vector<double> centralityWeights;
    vector<TH1D*> centralityHists;
    vector<TGraphErrors*> centralityBands;



    ////////////////////////////////////////////////////////////
    // Open file
    ////////////////////////////////////////////////////////////

    TFile *fTPC = TFile::Open("Oxygen_Oxygen_systematics_TPC_train_661261_1.root");

    if (!fTPC || fTPC->IsZombie())
    {
        cout << "Error opening file\n";
        return 1;
    }

    TFile *fPID = TFile::Open("Oxygen_Oxygen_systematics_PID_train_663089_1.root");

    if (!fPID || fPID->IsZombie())
    {
        cout << "Error opening file\n";
        return 1;
    }

    TFile *fITS = TFile::Open("Oxygen_Oxygen_systematics_ITS_train_666987_1.root");

    if (!fITS || fITS->IsZombie())
    {
        cout << "Error opening file\n";
        return 1;
    }
    TFile *foccu = TFile::Open("pieter_train_669682_occu.root");

    if (!foccu || foccu->IsZombie())
    {
        cout << "Error opening file\n";
        return 1;
    }
    TFile *fTotal = TFile::Open("Oxygen_Oxygen_TPC_weighted_flow_train_661261_jpsi.root");

    if (!fTotal || fTotal->IsZombie())
    {
        cout << "Error opening file\n";
        return 1;
    }
    TDirectory* dir = (TDirectory*) fTotal->Get("centrality_base50");

    if (!dir) {
        std::cout << "Directory centrality_base50 not found!" << std::endl;
        return 1;
    }


    

    for (auto cent : centBins)
    {
        cout << "\nProcessing centrality "
             << cent << endl;

             string centStr = 
             to_string(cent);

             string centLabel;

if (cent == 0)
    centLabel = "0#minus10%";
else if (cent == 10)
    centLabel = "10#minus20%";
else if (cent == 20)
    centLabel = "20#minus30%";
else if (cent == 30)
    centLabel = "30#minus50%";
    ////////////////////////////////////////////////////////////
    // Load nominal histogram
    ////////////////////////////////////////////////////////////

    TH1D* hTPCBase =(TH1D*)((TDirectory*)fTPC->Get(("centralitybase" + centStr).c_str()))->Get(("flowpion_Signal_base"+centStr).c_str());
    TH1D* hTPCML =(TH1D*)((TDirectory*)fTPC->Get(("centralityTPCML_" + centStr).c_str()))->Get(("flowpion_Signal_TPCML_"+centStr).c_str());
    TH1D* hTPCMS =(TH1D*)((TDirectory*)fTPC->Get(("centralityTPCMS_" + centStr).c_str()))->Get(("flowpion_Signal_TPCMS_"+centStr).c_str());

    TH1D* hITSBase =(TH1D*)((TDirectory*)fITS->Get(("centralitybase" + centStr).c_str()))->Get(("flowpion_Signal_base"+centStr).c_str());
    TH1D* hITSML =(TH1D*)((TDirectory*)fITS->Get(("centralityITSML" + centStr).c_str()))->Get(("flowpion_Signal_ITSML"+centStr).c_str());
    TH1D* hITSMS =(TH1D*)((TDirectory*)fITS->Get(("centralityITSMS" + centStr).c_str()))->Get(("flowpion_Signal_ITSMS"+centStr).c_str());

    TH1D* hPIDBase =(TH1D*)((TDirectory*)fPID->Get(("centralitybase" + centStr).c_str()))->Get(("flowpion_Signal_base"+centStr).c_str());
    TH1D* hPIDML =(TH1D*)((TDirectory*)fPID->Get(("centralityPIDML_" + centStr).c_str()))->Get(("flowpion_Signal_PIDML_"+centStr).c_str());
    TH1D* hPIDMS =(TH1D*)((TDirectory*)fPID->Get(("centralityPIDMS_" + centStr).c_str()))->Get(("flowpion_Signal_PIDMS_"+centStr).c_str());

    TH1D* hoccuBase =(TH1D*)((TDirectory*)foccu->Get(("centralitybase_" + centStr).c_str()))->Get(("flowpion_Signal_base_"+centStr).c_str());
    TH1D* hoccuML =(TH1D*)((TDirectory*)foccu->Get(("centralityoccuML_" + centStr).c_str()))->Get(("flowpion_Signal_occuML_"+centStr).c_str());
    TH1D* hoccuMS =(TH1D*)((TDirectory*)foccu->Get(("centralityoccuMS_" + centStr).c_str()))->Get(("flowpion_Signal_occuMS_"+centStr).c_str());

    ////////////////////////////////////////////////////////////
    // TPC systematics
    ////////////////////////////////////////////////////////////


    vector<TH1D*> tpcVars;

    AddVariation(tpcVars, hTPCML);
    AddVariation(tpcVars, hTPCMS);
    
    TH1D* hTPCSyst =
    ComputeALICESystematic(
        hTPCBase,
        tpcVars,
        ("TPCSyst_" + centStr).c_str());



    ////////////////////////////////////////////////////////////
    //  ITS  systematics
    ////////////////////////////////////////////////////////////

    vector<TH1D*> itsVars;

    AddVariation(itsVars, hITSML);
    AddVariation(itsVars, hITSMS);

    TH1D* hITSSyst =
        ComputeALICESystematic(
            hITSBase,
            itsVars,
            ("ITSSyst_" + centStr).c_str());
            std::cout << "Loaded histograms for centrality " << centStr << endl;

    ////////////////////////////////////////////////////////////
    //  PID systematics
    ////////////////////////////////////////////////////////////

        vector<TH1D*> pidVars;

        AddVariation(pidVars, hPIDML);
        AddVariation(pidVars, hPIDMS);
        std::cout << "test 1 " << centStr << endl;
        TH1D* hPIDSyst =
            ComputeALICESystematic(
                hPIDBase,
                pidVars,
                ("PIDSyst_" + centStr).c_str());

                std::cout << "test 2 " << centStr << endl;


        vector<TH1D*> occuVars;

        AddVariation(occuVars, hoccuML);
        AddVariation(occuVars, hoccuMS);
        std::cout << "test 1 " << centStr << endl;
        TH1D* hoccuSyst =
            ComputeALICESystematic(
                hoccuBase,
                occuVars,
                ("occuSyst_" + centStr).c_str());                



    ////////////////////////////////////////////////////////////
    // Combine systematic sources
    ////////////////////////////////////////////////////////////

    vector<TH1D*> systs =
    {
        hTPCSyst,
        hPIDSyst,
        hITSSyst,
        hoccuSyst
    };

    TH1D* hTotal =
    CombineSystematics(
        systs,
        ("Total_" + centStr).c_str());

            centralityTotals.push_back(hTotal);

            cout << "\n========================================\n";
            cout << "Centrality " << centStr << "\n";
            
            cout << "Bin  Nominal  "
                 << "TPC(abs)  PID(abs)  ITS(abs) occu(abs) Total(abs)\n";
            
            int nBins = hTPCBase->GetNbinsX();
            
            for (int i = 1; i <= nBins; i++)
            {
                double nominal = hTPCBase->GetBinContent(i);
            
                double tpcAbs   = hTPCSyst->GetBinContent(i);
                double pidAbs   = hPIDSyst->GetBinContent(i);
                double itsAbs   = hITSSyst->GetBinContent(i);
                double occuAbs  = hoccuSyst->GetBinContent(i);
                double totalAbs = hTotal->GetBinContent(i);
            

                cout
                    << i << "   "
                    << nominal << "   "
                    << tpcAbs << "   "
                    << pidAbs << "   "
                    << itsAbs << "   "
                    << occuAbs << "   "
                    << totalAbs << "   "
                    << endl;
            }

            cout << "\n========================================\n";
cout << "Centrality " << centStr << "\n";

cout << "Bin  Nominal  "
     << "TPC(%)  PID(%)  ITS(%) occu(%) Total(%)\n";



for (int i = 1; i <= nBins; i++)
{
    double nominal = hTPCBase->GetBinContent(i);

    double tpcAbs   = hTPCSyst->GetBinContent(i);
    double pidAbs   = hPIDSyst->GetBinContent(i);
    double itsAbs   = hITSSyst->GetBinContent(i);
    double occuAbs  = hoccuSyst->GetBinContent(i);
    double totalAbs = hTotal->GetBinContent(i);

    double tpcPct   = (nominal != 0.0) ? 100.0 * tpcAbs   / fabs(nominal) : 0.0;
    double pidPct   = (nominal != 0.0) ? 100.0 * pidAbs   / fabs(nominal) : 0.0;
    double itsPct   = (nominal != 0.0) ? 100.0 * itsAbs   / fabs(nominal) : 0.0;
    double occuPct   = (nominal != 0.0) ? 100.0 * occuAbs   / fabs(nominal) : 0.0;
    double totalPct = (nominal != 0.0) ? 100.0 * totalAbs / fabs(nominal) : 0.0;

    cout
        << i << "   "
        << nominal << "   "
        << tpcPct << "   "
        << pidPct << "   "
        << itsPct << "   "
        << occuPct << "   "
        << totalPct
        << endl;
}
    ////////////////////////////////////////////////////////////
    // Plot
    ////////////////////////////////////////////////////////////

    TGraphErrors* gBand =
    MakeSystematicBand(hTPCBase, hTotal);

    TCanvas* c = new TCanvas(
        ("c_" + centStr).c_str(),
        ("Centrality " + centStr).c_str(),
        800,800);
        c->SetTicks(1,1);
        c->SetLeftMargin(0.14);
        c->SetBottomMargin(0.12);
        c->SetRightMargin(0.04);
        c->SetTopMargin(0.02);
    
    int color = kBlack;
    
    // Marker style (open)
    hTPCBase->SetMarkerStyle(24);
    hTPCBase->SetMarkerSize(1.2);
    hTPCBase->SetMarkerColor(color);
    hTPCBase->SetLineColor(color);
    hTPCBase->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hTPCBase->GetYaxis()->SetTitle("v_{2}");
    hTPCBase->GetYaxis()->SetTitleOffset(1.2);

    
    
    // Axis range
    hTPCBase->GetXaxis()->SetRangeUser(0, 5);
    hTPCBase->GetYaxis()->SetRangeUser(0, 0.3);
    
    // Systematic band style
    gBand->SetFillColorAlpha(kGray+1, 0.35);
    gBand->SetLineColor(kGray+2);
    gBand->SetLineWidth(1);
    
    // Draw order
    hTPCBase->Draw("E1");
    gBand->Draw("E2 SAME");
    hTPCBase->Draw("E1 SAME");
    
    // -----------------------------------
    // Legend
    // -----------------------------------
    TLegend* leg = new TLegend(0.58, 0.82, 0.88, 0.96);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.03);
    
    leg->AddEntry(hTPCBase, "Data", "lep");
    leg->AddEntry(gBand, "Systematic uncertainty", "f");
    
    leg->Draw();
    
    // -----------------------------------
    // Latex text
    // -----------------------------------
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    

    latex.DrawLatex(0.17, 0.92, "O-O, #sqrt{s_{NN}} = 5.36 TeV");
    latex.DrawLatex(0.17, 0.88, Form("%s, mass < 0.14 GeV/c^{2}", centLabel.c_str()));
    latex.DrawLatex(0.17, 0.84, "p_{T} > 0.2 GeV/c, |#eta| < 0.8");

    
    // Save
    c->SaveAs(("Systematics_" + centStr + ".pdf").c_str());
    fOut->cd();
c->Write(("c_" + centStr).c_str());
gBand->Write(("gBand_" + centStr).c_str());
hTPCBase->Write(("hTPCBase_" + centStr).c_str());
hTotal->Write(("hTotalSyst_" + centStr).c_str());
TCanvas* cTPCVar = new TCanvas(
    ("cTPCVar_" + centStr).c_str(),
    ("TPC variations " + centStr).c_str(),
    800, 800);

cTPCVar->SetTicks(1,1);
cTPCVar->SetLeftMargin(0.14);
cTPCVar->SetBottomMargin(0.12);
cTPCVar->SetRightMargin(0.04);
cTPCVar->SetTopMargin(0.02);

// Style nominal
hTPCBase->SetMarkerStyle(20);
hTPCBase->SetMarkerSize(1.2);
hTPCBase->SetMarkerColor(kBlack);
hTPCBase->SetLineColor(kBlack);
hTPCBase->SetTitle(";p_{T} (GeV/c);v_{2}");

// Style variations
hTPCML->SetMarkerStyle(24);
hTPCML->SetMarkerSize(1.1);
hTPCML->SetMarkerColor(kRed);
hTPCML->SetLineColor(kRed);

hTPCMS->SetMarkerStyle(25);
hTPCMS->SetMarkerSize(1.1);
hTPCMS->SetMarkerColor(kBlue);
hTPCMS->SetLineColor(kBlue);

// Axis range (match your analysis)
hTPCBase->GetXaxis()->SetRangeUser(0, 5);
hTPCBase->GetYaxis()->SetRangeUser(0, 0.3);

// Draw order
hTPCBase->Draw("E1");
hTPCML->Draw("E1 SAME");
hTPCMS->Draw("E1 SAME");

// Legend
TLegend* legTPC = new TLegend(0.58, 0.80, 0.88, 0.95);
legTPC->SetBorderSize(0);
legTPC->SetFillStyle(0);
legTPC->SetTextSize(0.03);

legTPC->AddEntry(hTPCBase, "TPC standard", "lep");
legTPC->AddEntry(hTPCML, "TPC less strict variation", "lep");
legTPC->AddEntry(hTPCMS, "TPC more strict variation", "lep");

legTPC->Draw();

// Label

latex.SetNDC();
latex.SetTextSize(0.035);
latex.DrawLatex(0.17, 0.92, "TPC systematics check");
latex.DrawLatex(0.17, 0.88, ("Centrality " + centLabel).c_str());

// Save + write
cTPCVar->SaveAs(("TPC_variations_" + centStr + ".png").c_str());

fOut->cd();
cTPCVar->Write(("cTPCVar_" + centStr).c_str());

TCanvas* cITSVar = new TCanvas(
    ("cITSVar_" + centStr).c_str(),
    ("ITS variations " + centStr).c_str(),
    800, 800);

cITSVar->SetTicks(1,1);
cITSVar->SetLeftMargin(0.14);
cITSVar->SetBottomMargin(0.12);
cITSVar->SetRightMargin(0.04);
cITSVar->SetTopMargin(0.02);

// Style nominal
hITSBase->SetMarkerStyle(20);
hITSBase->SetMarkerSize(1.2);
hITSBase->SetMarkerColor(kBlack);
hITSBase->SetLineColor(kBlack);
hITSBase->SetTitle(";p_{T} (GeV/c);v_{2}");

// ITS variations
hITSML->SetMarkerStyle(24);
hITSML->SetMarkerSize(1.1);
hITSML->SetMarkerColor(kRed);
hITSML->SetLineColor(kRed);

hITSMS->SetMarkerStyle(25);
hITSMS->SetMarkerSize(1.1);
hITSMS->SetMarkerColor(kBlue);
hITSMS->SetLineColor(kBlue);

// Axis range
hITSBase->GetXaxis()->SetRangeUser(0, 5);
hITSBase->GetYaxis()->SetRangeUser(0, 0.3);

// Draw order
hITSBase->Draw("E1");
hITSML->Draw("E1 SAME");
hITSMS->Draw("E1 SAME");

// Legend
TLegend* legITS = new TLegend(0.58, 0.8, 0.88, 0.95);
legITS->SetBorderSize(0);
legITS->SetFillStyle(0);
legITS->SetTextSize(0.03);

legITS->AddEntry(hITSBase, "ITS standard", "lep");
legITS->AddEntry(hITSML, "ITS less strict variation", "lep");
legITS->AddEntry(hITSMS, "ITS more strict variation", "lep");

legITS->Draw();

// Labels

latex.SetNDC();
latex.SetTextSize(0.035);
latex.DrawLatex(0.17, 0.92, "ITS systematics check");
latex.DrawLatex(0.17, 0.88, ("Centrality " + centLabel).c_str());

// Save to file
cITSVar->SaveAs(("ITS_variations_" + centStr + ".png").c_str());

fOut->cd();
cITSVar->Write(("cITSVar_" + centStr).c_str());

centralityHists.push_back(hTPCBase);
centralityBands.push_back(gBand);
    }
    TFile* fTPC90_base = TFile::Open("pieter_oxygen_train_661261_TPC_base_low.root");
    TFile* fTPC90_ML   = TFile::Open("pieter_oxygen_train_661261_TPC_ML_low.root");
    TFile* fTPC90_MS   = TFile::Open("pieter_oxygen_train_661261_TPC_MS_low.root");
    
    TFile* fITS90_base = TFile::Open("pieter_oxygen_train_666987_ITS_base_low.root");
    TFile* fITS90_ML   = TFile::Open("pieter_oxygen_train_666987_ITS_ML_low.root");
    TFile* fITS90_MS   = TFile::Open("pieter_oxygen_train_666987_ITS_MS_low.root");
    
    TFile* fPID90_base = TFile::Open("pieter_oxygen_train_663089_PID_base_low.root");
    TFile* fPID90_ML   = TFile::Open("pieter_oxygen_train_663089_PID_ML_low.root");
    TFile* fPID90_MS   = TFile::Open("pieter_oxygen_train_663089_PID_MS_low.root");

    TFile* foccu90_base = TFile::Open("pieter_oxygen_train_669682_occu_base_low.root");
    TFile* foccu90_ML   = TFile::Open("pieter_oxygen_train_669682_occu_ML_low.root");
    TFile* foccu90_MS   = TFile::Open("pieter_oxygen_train_669682_occu_MS_low.root");


    TH1D* hTPCBase90 = (TH1D*)((TDirectory*)fTPC90_base->Get("centrality_base50"))->Get("plot_of_v2S_massrange_1");
    TH1D* hTPCML90 = (TH1D*)((TDirectory*)fTPC90_ML->Get("centrality_TPCML_50"))->Get("plot_of_v2S_massrange_1");
    TH1D* hTPCMS90 = (TH1D*)((TDirectory*)fTPC90_MS->Get("centrality_TPCMS_50"))->Get("plot_of_v2S_massrange_1");

    TH1D* hITSBase90 = (TH1D*)((TDirectory*)fITS90_base->Get("centrality_base50"))->Get("plot_of_v2S_massrange_1");
    TH1D* hITSML90 = (TH1D*)((TDirectory*)fITS90_ML->Get("centrality_ITSML50"))->Get("plot_of_v2S_massrange_1");
    TH1D* hITSMS90 = (TH1D*)((TDirectory*)fITS90_MS->Get("centrality_ITSMS50"))->Get("plot_of_v2S_massrange_1");

    TH1D* hPIDBase90 = (TH1D*)((TDirectory*)fPID90_base->Get("centrality_base50"))->Get("plot_of_v2S_massrange_1");
    TH1D* hPIDML90 = (TH1D*)((TDirectory*)fPID90_ML->Get("centrality_PIDML_50"))->Get("plot_of_v2S_massrange_1");
    TH1D* hPIDMS90 = (TH1D*)((TDirectory*)fPID90_MS->Get("centrality_PIDMS_50"))->Get("plot_of_v2S_massrange_1");

    TH1D* hoccuBase90 = (TH1D*)((TDirectory*)foccu90_base->Get("centrality_base_50"))->Get("plot_of_v2S_massrange_1");
    TH1D* hoccuML90 = (TH1D*)((TDirectory*)foccu90_ML->Get("centrality_occuML_50"))->Get("plot_of_v2S_massrange_1");
    TH1D* hoccuMS90 = (TH1D*)((TDirectory*)foccu90_MS->Get("centrality_occuMS_50"))->Get("plot_of_v2S_massrange_1");
    double ptBins[] = {0,1,2,4,6,10};

TH1D *hTPCBase90_fixed = new TH1D(
    "hTPCBase90_fixed",
    hTPCBase90->GetTitle(),
    5,
    ptBins);

for (int i = 1; i <= hTPCBase90->GetNbinsX(); ++i) {
    hTPCBase90_fixed->SetBinContent(i, hTPCBase90->GetBinContent(i));
    hTPCBase90_fixed->SetBinError(i, hTPCBase90->GetBinError(i));
} 
    vector<TH1D*> tpcVars90;
    AddVariation(tpcVars90, hTPCML90);
    AddVariation(tpcVars90, hTPCMS90);
    
    TH1D* hTPCSyst90 =
    ComputeALICESystematic(hTPCBase90_fixed, tpcVars90, "TPCSyst_90");
    
    vector<TH1D*> itsVars90;
    AddVariation(itsVars90, hITSML90);
    AddVariation(itsVars90, hITSMS90);
    
    TH1D* hITSSyst90 =
    ComputeALICESystematic(hITSBase90, itsVars90, "ITSSyst_90");
    
    
    vector<TH1D*> pidVars90;
    AddVariation(pidVars90, hPIDML90);
    AddVariation(pidVars90, hPIDMS90);
    
    TH1D* hPIDSyst90 =
    ComputeALICESystematic(hPIDBase90, pidVars90, "PIDSyst_90");

    vector<TH1D*> occuVars90;
    AddVariation(occuVars90, hoccuML90);
    AddVariation(occuVars90, hoccuMS90);
    
    TH1D* hoccuSyst90 =
    ComputeALICESystematic(hoccuBase90, occuVars90, "occuSyst_90");
    
    
    vector<TH1D*> syst90 = {
        hTPCSyst90,
        hITSSyst90,
        hPIDSyst90,
        hoccuSyst90
    };
    
    TH1D* hTotal90 =
        CombineSystematics(syst90, "Total_90");
    
        TGraphErrors* gBand90 =
    MakeSystematicBand(hTPCBase90_fixed, hTotal90);
    
    
    TCanvas* c90 = new TCanvas("c90", "0-90 systematics", 800, 800);
    c90->SetTicks(1,1);
    c90->SetLeftMargin(0.14);
    c90->SetBottomMargin(0.12);
    c90->SetRightMargin(0.04);
    c90->SetTopMargin(0.02);
    
    // -------------------------
    // Style nominal histogram
    // -------------------------
    hTPCBase90_fixed->SetMarkerStyle(24);
    hTPCBase90_fixed->SetMarkerSize(1.2);
    hTPCBase90_fixed->SetMarkerColor(kBlack);
    hTPCBase90_fixed->SetLineColor(kBlack);
    hTPCBase90_fixed->SetTitle(";p_{T} (GeV/c);v_{2}");
    hTPCBase90_fixed->SetStats(0);
    
    hTPCBase90_fixed->GetXaxis()->SetRangeUser(0, 10);
    hTPCBase90_fixed->GetYaxis()->SetRangeUser(-0.1, 0.3);
    hTPCBase90_fixed->GetYaxis()->SetTitleOffset(1.2);
    
    // -------------------------
    // Style systematic band
    // -------------------------
    gBand90->SetFillColorAlpha(kGray+1, 0.35);
    gBand90->SetLineColor(kGray+2);
    gBand90->SetLineWidth(1);
    
    // -------------------------
    // Draw order
    // -------------------------
    hTPCBase90_fixed->Draw("E1");
    gBand90->Draw("E2 SAME");
    hTPCBase90_fixed->Draw("E1 SAME");
    
    // -------------------------
    // Legend
    // -------------------------
    TLegend* leg90 = new TLegend(0.58, 0.75, 0.88, 0.9);
    leg90->SetBorderSize(0);
    leg90->SetFillStyle(0);
    leg90->SetTextSize(0.03);
    
    leg90->AddEntry(hTPCBase90_fixed, "data", "lep");
    leg90->AddEntry(gBand90, "Systematic uncertainty", "f");
    
    leg90->Draw();
    
    // -------------------------
    // Labels
    // -------------------------
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    
    latex.DrawLatex(0.17, 0.92, "O-O, #sqrt{s_{NN}} = 5.36 TeV");
    latex.DrawLatex(0.17, 0.88, "0-90\%, 2.7 < mass < 3.2 GeV/c^{2}");
    latex.DrawLatex(0.17, 0.84, "p_{T} > 0.2 GeV/c, |#eta| < 0.8");
    
    // -------------------------
    // Save outputs
    // -------------------------
    c90->SaveAs("Systematics_0_90_low.pdf");
    c90->SaveAs("Systematics_0_90.png");
    
    // optional ROOT output
    fOut->cd();
    c90->Write("c90");
    gBand90->Write("gBand90");
    hTPCBase90_fixed->Write("hTPCBase90_fixed");

std:cout << "test" << endl; 
    cout << "Bin  Nominal  TPC(%)  PID(%)  ITS(%) occu(%) Total(%)\n";

    for (int i = 1; i <= hTPCBase90_fixed->GetNbinsX(); i++) {
    
        double nominal = hTPCBase90_fixed->GetBinContent(i);
    
        double tpcAbs  = hTPCSyst90->GetBinContent(i);
        double pidAbs  = hPIDSyst90->GetBinContent(i);
        double itsAbs  = hITSSyst90->GetBinContent(i);
        double occuAbs = hoccuSyst90->GetBinContent(i);
        double totalAbs= hTotal90->GetBinContent(i);
    
        double denom = fabs(nominal);
    
        double tpcPct  = (denom != 0.0) ? 100.0 * tpcAbs  / denom : 0.0;
        double pidPct  = (denom != 0.0) ? 100.0 * pidAbs  / denom : 0.0;
        double itsPct  = (denom != 0.0) ? 100.0 * itsAbs  / denom : 0.0;
        double occuPct = (denom != 0.0) ? 100.0 * occuAbs / denom : 0.0;
        double totPct  = (denom != 0.0) ? 100.0 * totalAbs/ denom : 0.0;
    
        cout << i << "  "
             << nominal << "  "
             << tpcPct << "  "
             << pidPct << "  "
             << itsPct << "  "
             << occuPct << "  "
             << totPct
             << endl;
    }

   
TCanvas* cAll = new TCanvas("cAll", "All centralities", 800, 800);
cAll->SetTicks(1,1);
cAll->SetLeftMargin(0.14);
cAll->SetBottomMargin(0.12);
cAll->SetRightMargin(0.04);
cAll->SetTopMargin(0.02);

// nice color palette
vector<int> colors = {
    kBlack,
    kBlue,
    kGreen+1,
    kRed
};

// Legend
vector<string> centLabels = {
    "0#minus10%",
    "10#minus20%",
    "20#minus30%",
    "30#minus50%"
};
TLegend* legAll = new TLegend(0.65, 0.75, 0.88, 0.96);
legAll->SetBorderSize(0);
legAll->SetFillStyle(0);

for (size_t i = 0; i < centralityHists.size(); i++)
{
    int color = colors[i % colors.size()];
    int marker = markerStyles[i % markerStyles.size()];


    auto h = centralityHists[i];
    auto g = centralityBands[i];

    // Style
    h->SetMarkerStyle(marker);
    h->SetMarkerColor(color);
    h->SetLineColor(color);

    g->SetFillColorAlpha(color, 0.25);
    g->SetLineColor(color);

    if (i == 0)
    {
        h->SetTitle(";p_{T} (GeV/c);v_{2}");
        h->GetXaxis()->SetRangeUser(0,5);
        h->GetYaxis()->SetRangeUser(0,0.3);
        h->Draw("E1");
    }
    else
    {
        h->Draw("E1 SAME");
    }

    g->Draw("E2 SAME");

    legAll->AddEntry(
        h,
        ((centLabels[i])).c_str(),
        "lep");
}

legAll->Draw();


latex.SetNDC();
latex.SetTextSize(0.035);
latex.DrawLatex(0.17, 0.92, "O-O #sqrt{s_{NN}} = 5.36 TeV");
latex.DrawLatex(0.17, 0.88, "mass < 0.14 GeV/c^{2}");
latex.DrawLatex(0.17, 0.84, "p_{T} > 0.2 GeV/c, |#eta| < 0.8");

cAll->SaveAs("AllCentralities.pdf");
fOut->cd();
cAll->Write("cAll");
fOut->Write();
fOut->Close();
delete fOut;


return 0;
}