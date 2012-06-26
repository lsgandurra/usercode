//#include "TROOT.h "
#include "TF1.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TBranch.h"
#include "TChain.h"
#include "TFile.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include "TBits.h"
#include "TMath.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TPad.h"
#include "TGraph.h"
#include "TText.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TProfile.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <iomanip>
#include <vector>
#include "RooAbsPdf.h"
#include "RooCBShape.h"
#include "RooLandau.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooHistPdf.h"
#include "RooHist.h"
#include "RooMsgService.h"
#include "RooNLLVar.h"
#include "RooPlot.h"
#include "RooRandom.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooFFTConvPdf.h"
#include "RooLognormal.h"
#include "RooGaussian.h"
#include "RooGamma.h"
#include "RooBifurGauss.h"
#include "RooAddPdf.h"
#include "RooGenericPdf.h"
#include "RooProdPdf.h"
#include "RooDerivative.h"
#include "TROOT.h"
#include "TRint.h"
#include "CrystalBall.C"
#include "setTDRStyle.C"
#include "CMSStyle.C"
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
#include <algorithm>
#include <vector>
#include <set>
#include <list>
#include <time.h>
#include <stdio.h>
#pragma optimize 0

using namespace RooFit;
using namespace std;

int NbLignesFichier(string fichier);
string DoubleToString(double x);
void enregistrementPlots(string nomDossier, string nomFichier, int EndCaps, int iteration, TCanvas * c1);

int NbLignesFichier(string fichier)
{
        ifstream in(fichier.c_str()); //Ouverture en mode lecture de fichier

        string ligne; //Création d'une chaine de caracterec
        int nbLignes = 0;

        while(std::getline(in, ligne)) nbLignes++;

        in.close(); //On ferme le fichier

        return nbLignes;
}


string DoubleToString(double x)
{

        std::string s;
        {
                std::ostringstream oss;
                oss << x;
                s = oss.str();
        }
        std::cout << "x = " << x << " s = " << s << std::endl;

        return s;
}

void enregistrementPlots(string nomDossier, string nomFichier, int EndCaps, int iteration, TCanvas * c1)
{

        if(EndCaps == 0) nomDossier += "EB/";
        if(EndCaps == 1) nomDossier += "EE/";
        //mkdir(nomDossier.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        system(Form("mkdir -p %s", nomDossier.c_str()));
        if(iteration != 10000) nomFichier += Form("%d",iteration);
        c1->Print(Form("%s%s.root",nomDossier.c_str(),nomFichier.c_str()));
        c1->Print(Form("%s%s.C",nomDossier.c_str(),nomFichier.c_str()));
        c1->Print(Form("%s%s.pdf",nomDossier.c_str(),nomFichier.c_str()));
        c1->Print(Form("%s%s.ps",nomDossier.c_str(),nomFichier.c_str()));
        c1->Print(Form("%s%s.png",nomDossier.c_str(),nomFichier.c_str()));
        return;
}

int SVsPtProfile(string variableX = "Photon_Et", string version = "true", int phiCracks = 1, int etaCracks = 1, string MZbinning = "05GeV")
{

	gROOT->Reset();
        TGaxis::SetMaxDigits(3);
        setTDRStyle();
        gStyle->SetPalette(1,0);

	double yminErecoOverEtrue, ymaxErecoOverEtrue;
        double xminErecoOverEtrue, xmaxErecoOverEtrue;

        yminErecoOverEtrue = -15.0;
        ymaxErecoOverEtrue = 15.0;
        xminErecoOverEtrue = 0.0;
        xmaxErecoOverEtrue = 100.0;

	int nBins = 6;
        double t_PhotonPt[7] = {0,12,15,20,25,30,100};

        double xError[6] = {0};
        double xValueTab[6] = {0};

	for(int b = 0; b < nBins; b++)
        {
                if(b == 0)
                {
                        xError[b] = (t_PhotonPt[b+1] - 10.0) / 2.0;
                        xValueTab[b] = (t_PhotonPt[b+1] + 10.0) / 2.0;

                }

                if(b != 0)
                {
                        xError[b] = (t_PhotonPt[b+1] - t_PhotonPt[b]) / 2.0;
                        xValueTab[b] = (t_PhotonPt[b+1] + t_PhotonPt[b]) / 2.0;
                }
                //cout<<endl<<"xError[b] = "<<xError[b]<<endl<<"xValueTab[b] = "<<xValueTab[b];

        }

	TString temp = "";
        TString tempVarChain = "";

	double r9infEBTab[6] = {0};
        double r9infEBErrorTab[6] = {0};
        double r9infEETab[6] = {0};
        double r9infEEErrorTab[6] = {0};
        double r9supEBTab[6] = {0};
        double r9supEBErrorTab[6] = {0};
        double r9supEETab[6] = {0};
        double r9supEEErrorTab[6] = {0};



	TChain * chain = new TChain("miniTree");

	string LowMmumuLim = "40";
	string HightMmumuLim = "80";
	string compressSetOfCorrectionsChain = "ETHZ";
	string nomDossier = "";
	string nomFichier = "";
	string SetOfCorrections = "ETHZCorrections";
	string r9supC = "";
	string EndCapsC = "";



	//chain->Add(Form("/sps/cms/sgandurr/CMSSW_4_2_8/src/UserCode/IpnTreeProducer/ZMuMuGammaMinitrees_Mars_2012/miniTree_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_PU_S6_v2_1_%s_%s_%s_partALL.root",compressSetOfCorrectionsChain.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str()));
        //chain->Add(Form("/sps/cms/sgandurr/CMSSW_4_2_8/src/UserCode/IpnTreeProducer/ZMuMuGammaMinitrees_Mars_2012/miniTree_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_PU_S6_v2_2_%s_%s_%s_partALL.root",compressSetOfCorrectionsChain.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str()));
        chain->Add(Form("/sps/cms/sgandurr/CMSSW_4_2_8/src/UserCode/IpnTreeProducer/miniTree_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_PU_S6_v2_1_%s_%s_%s_%s_partALL.root",compressSetOfCorrectionsChain.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),MZbinning.c_str()));
        chain->Add(Form("/sps/cms/sgandurr/CMSSW_4_2_8/src/UserCode/IpnTreeProducer/miniTree_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_PU_S6_v2_2_%s_%s_%s_%s_partALL.root",compressSetOfCorrectionsChain.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),MZbinning.c_str()));


	for(int EndCaps = 0; EndCaps < 2; EndCaps++)
        {
        for(int r9sup = 0; r9sup < 2; r9sup++)
        {

                TCanvas *c2 = new TCanvas("c2", "c2",0,0,600,600);
		TH2F* MC = new TH2F("MC", "MC", nBins, &t_PhotonPt[0], 10000,-100,100);

                temp.Clear();
                tempVarChain.Clear();

                temp += "isTightMMG == 1 && isMultipleCandidate == 0";

		if(EndCaps == 0 && r9sup == 1 && variableX != "Photon_SC_Eta") temp += " && (abs(Photon_SC_Eta)) < 1.4442 && Photon_r9 > 0.94";
                if(EndCaps == 0 && r9sup == 0 && variableX != "Photon_SC_Eta") temp += " && (abs(Photon_SC_Eta)) < 1.4442 && Photon_r9 < 0.94";
                if(EndCaps == 1 && r9sup == 1 && variableX != "Photon_SC_Eta") temp += " && (abs(Photon_SC_Eta)) > 1.566 && Photon_r9 > 0.95";
                if(EndCaps == 1 && r9sup == 0 && variableX != "Photon_SC_Eta") temp += " && (abs(Photon_SC_Eta)) > 1.566 && Photon_r9 < 0.95";

                if(EndCaps == 0 && r9sup == 1 && variableX == "Photon_SC_Eta") temp += " && Photon_r9 > 0.94";
                if(EndCaps == 0 && r9sup == 0 && variableX == "Photon_SC_Eta") temp += " && Photon_r9 < 0.94";
                if(EndCaps == 1 && r9sup == 1 && variableX == "Photon_SC_Eta") temp += " && Photon_r9 > 0.95";
                if(EndCaps == 1 && r9sup == 0 && variableX == "Photon_SC_Eta") temp += " && Photon_r9 < 0.95";


                if(EndCaps == 0 && etaCracks == false) temp += " && ( (abs(Photon_SC_Eta)) > 0.018 || ( (abs(Photon_SC_Eta)) < 0.423 && (abs(Photon_SC_Eta)) > 0.461 ) || ( (abs(Photon_SC_Eta)) < 0.770 && (abs(Photon_SC_Eta)) > 0.806 ) || ( (abs(Photon_SC_Eta)) < 1.127 && (abs(Photon_SC_Eta)) > 1.163 ) )";

                if(EndCaps == 0 && phiCracks == false)
                {
                        temp += Form(" && ((abs(Photon_SC_Phi) + (%f)) < (%f) ||",phiOffset,phiCrackSize);
                        temp += Form(" ( (1.0*(%f) - (%f)) < (abs(Photon_SC_Phi) + (%f)) && (abs(Photon_SC_Phi) + (%f)) < (1.0*(%f) + (%f)) ) ||",phiCrackPosition,phiCrackSize,phiOffset,phiOffset,phiCrackPosition,phiCrackSize);
                        temp += Form(" ( (2.0*(%f) - (%f)) < (abs(Photon_SC_Phi) + (%f)) && (abs(Photon_SC_Phi) + (%f)) < (2.0*(%f) + (%f)) ) ||",phiCrackPosition,phiCrackSize,phiOffset,phiOffset,phiCrackPosition,phiCrackSize);
                        temp += Form(" ( (3.0*(%f) - (%f)) < (abs(Photon_SC_Phi) + (%f)) && (abs(Photon_SC_Phi) + (%f)) < (3.0*(%f) + (%f)) ) ||",phiCrackPosition,phiCrackSize,phiOffset,phiOffset,phiCrackPosition,phiCrackSize);
                        temp += Form(" ( (4.0*(%f) - (%f)) < (abs(Photon_SC_Phi) + (%f)) && (abs(Photon_SC_Phi) + (%f)) < (4.0*(%f) + (%f)) ) ||",phiCrackPosition,phiCrackSize,phiOffset,phiOffset,phiCrackPosition,phiCrackSize);
                        temp += Form(" ( (5.0*(%f) - (%f)) < (abs(Photon_SC_Phi) + (%f)) && (abs(Photon_SC_Phi) + (%f)) < (5.0*(%f) + (%f)) ) ||",phiCrackPosition,phiCrackSize,phiOffset,phiOffset,phiCrackPosition,phiCrackSize);
                        temp += Form(" ( (6.0*(%f) - (%f)) < (abs(Photon_SC_Phi) + (%f)) && (abs(Photon_SC_Phi) + (%f)) < (6.0*(%f) + (%f)) ) ||",phiCrackPosition,phiCrackSize,phiOffset,phiOffset,phiCrackPosition,phiCrackSize);
                        temp += Form(" ( (7.0*(%f) - (%f)) < (abs(Photon_SC_Phi) + (%f)) && (abs(Photon_SC_Phi) + (%f)) < (7.0*(%f) + (%f)) ) ||",phiCrackPosition,phiCrackSize,phiOffset,phiOffset,phiCrackPosition,phiCrackSize);
                        temp += Form(" ( (8.0*(%f) - (%f)) < (abs(Photon_SC_Phi) + (%f)) && (abs(Photon_SC_Phi) + (%f)) < (8.0*(%f) + (%f)) ) ||",phiCrackPosition,phiCrackSize,phiOffset,phiOffset,phiCrackPosition,phiCrackSize);
                        temp += Form(" ( (9.0*(%f) - (%f)) < (abs(Photon_SC_Phi) + (%f)) ) )",phiCrackPosition,phiCrackSize,phiOffset);

                }


		cout<<endl<<"temp = "<<temp<<endl;
		if(version == "reco") chain->Draw("mmg_s*100:Photon_Et>>MC",temp);
		if(version == "recoSurface") chain->Draw("mmg_s_MZ_Surface*100:Photon_Et>>MC",temp);
		if(version == "true") chain->Draw("((Photon_E_o_MC_E - 1)*100):Photon_Et>>MC",temp);
	        TProfile* MCP = (TProfile*)(MC->ProfileX());
		c2->Clear();

                //MC->Draw("");

                MCP->Draw("");

                if(r9sup == 0) r9supC = "inf";
                if(r9sup == 1) r9supC = "sup";

                if(EndCaps == 0) EndCapsC = "EB";
                if(EndCaps == 1) EndCapsC = "EE";

                MCP->GetXaxis()->SetLabelFont(42);
                MCP->GetXaxis()->SetTitleFont(42);
                MCP->GetXaxis()->SetLabelSize(0.03);
                if(version == "reco") MCP->GetYaxis()->SetTitle("s_{RECO} = 1/k_{RECO} - 1 (%)");
		if(version == "recoSurface") MCP->GetYaxis()->SetTitle("s_{RECO Surface} = 1/k_{RECO Surface} - 1 (%)");
		if(version == "true") MCP->GetYaxis()->SetTitle("s_{TRUE} = E_{RECO}/E_{TRUE} - 1 (%)");

                MCP->GetXaxis()->SetTitle("P_{T}^{#gamma}");
                MCP->GetYaxis()->SetLabelFont(42);
                MCP->GetYaxis()->SetTitleOffset(1.24);
                MCP->GetYaxis()->SetTitleFont(42);
                MCP->GetYaxis()->SetLabelSize(0.03);

                MCP->SetLineColor(413);
                MCP->SetMarkerStyle(20);
                MCP->SetMarkerSize(0.5);
                MCP->SetMarkerColor(413);
                MCP->GetYaxis()->SetRangeUser(-15,15);

		TLatex *text = new TLatex();
                text = new TLatex();
                text->SetNDC();
                text->SetTextAlign(11);
                text->SetTextFont(42);
                text->SetTextSizePixels(17);
                text->SetTextSize(0.028);
                text->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
                text->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");


                string * EndCapsR9Chain(0);
                EndCapsR9Chain = new string;
                if(EndCaps == 0) *EndCapsR9Chain = "Barrel, ";
                if(EndCaps == 1) *EndCapsR9Chain = "Endcaps, ";
                if(r9sup == 0) *EndCapsR9Chain += "Low r9";
                if(r9sup == 1) *EndCapsR9Chain += "High r9";
                text->DrawLatex(0.16, 0.80, EndCapsR9Chain->c_str());

                double entries = MC->GetEntries();

                int entriesInt = (int) entries;

		text->DrawLatex(0.16, 0.75, "ETHZ Corrections");
		
		nomDossier = Form("SrecoVsPt_v2/%s/2011/MmumuSelection%s_%s/WithCracks/%s/MZbinning_%s/",variableX.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),SetOfCorrections.c_str(),MZbinning.c_str());

		nomFichier = Form("SVsPt_%s%s_%s",r9supC.c_str(),EndCapsC.c_str(),version.c_str());

        	enregistrementPlots(nomDossier, nomFichier, 2, 10000, c2);


		if(EndCaps == 0 && r9sup == 0)
                {
                        for(int a = 1; a <= 6; a++)
                        {
                                r9infEBTab[a-1] = MCP->GetBinContent(a);
                                r9infEBErrorTab[a-1] = MCP->GetBinError(a);
                                //cout<<endl<<"a = "<<a<<endl<<"r9infEBTab[a] = "<<r9infEBTab[a-1]<<endl<<"r9infEBErrorTab[a] = "<<r9infEBErrorTab[a-1];
                        }
                }

                if(EndCaps == 1 && r9sup == 0)
                {
                        for(int a = 1; a <= 6; a++)
                        {
                                r9infEETab[a-1] = MCP->GetBinContent(a);
                                r9infEEErrorTab[a-1] = MCP->GetBinError(a);
                        }
                }

                if(EndCaps == 0 && r9sup == 1)
                {
                        for(int a = 1; a <= 6; a++)
                        {
                                r9supEBTab[a-1] = MCP->GetBinContent(a);
                                r9supEBErrorTab[a-1] = MCP->GetBinError(a);
                        }
                }

                if(EndCaps == 1 && r9sup == 1)
                {
                        for(int a = 1; a <= 6; a++)
                        {
                                r9supEETab[a-1] = MCP->GetBinContent(a);
                                r9supEEErrorTab[a-1] = MCP->GetBinError(a);
                        }
                }

		delete EndCapsR9Chain;
                text->Delete();
                MC->Delete();
                MCP->Delete();
                //palette->Delete();
                delete c2;
                //}


        }
        }


	TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);

        TMultiGraph *mg = new TMultiGraph();

        TGraphAsymmErrors * ErecoOverEtruer9infEB = new TGraphAsymmErrors(nBins,xValueTab, r9infEBTab, xError, xError, r9infEBErrorTab, r9infEBErrorTab);

        TGraphAsymmErrors * ErecoOverEtruer9supEB = new TGraphAsymmErrors(nBins,xValueTab, r9supEBTab, xError, xError, r9supEBErrorTab, r9supEBErrorTab);

        TGraphAsymmErrors * ErecoOverEtruer9infEE = new TGraphAsymmErrors(nBins,xValueTab, r9infEETab, xError, xError, r9infEEErrorTab, r9infEEErrorTab);

        TGraphAsymmErrors * ErecoOverEtruer9supEE = new TGraphAsymmErrors(nBins,xValueTab, r9supEETab, xError, xError, r9supEEErrorTab, r9supEEErrorTab);


        c1->ToggleEventStatus();
        //ErecoOverEtruer9infEB->Draw("AP");
        //ErecoOverEtruer9supEB->Draw("SAMES");
        //ErecoOverEtruer9infEE->Draw("SAMES");
        //ErecoOverEtruer9supEE->Draw("SAMES");

        mg->Add(ErecoOverEtruer9infEB);
        mg->Add(ErecoOverEtruer9supEB);
        mg->Add(ErecoOverEtruer9infEE);
        mg->Add(ErecoOverEtruer9supEE);
        mg->Draw("AP");


        mg->GetXaxis()->SetTitle("P_{T}^{#gamma}");

        mg->GetXaxis()->SetLabelFont(42);
        mg->GetXaxis()->SetTitleFont(42);
        mg->GetXaxis()->SetLabelSize(0.03);

        mg->GetYaxis()->SetTitle("s_{RECO} = 1/k_{RECO} - 1 (%)");
	if(version == "reco") mg->GetYaxis()->SetTitle("s_{RECO} = 1/k_{RECO} - 1 (%)");
        if(version == "recoSurface") mg->GetYaxis()->SetTitle("s_{RECO Surface} = 1/k_{RECO Surface} - 1 (%)");
        if(version == "true") mg->GetYaxis()->SetTitle("s_{TRUE} = E_{RECO}/E_{TRUE} - 1 (%)");


        mg->GetYaxis()->SetLabelFont(42);
        mg->GetYaxis()->SetTitleOffset(1.24);
        mg->GetYaxis()->SetTitleFont(42);
        mg->GetYaxis()->SetLabelSize(0.03);

        ErecoOverEtruer9infEB->SetLineColor(629);
        ErecoOverEtruer9supEB->SetLineColor(597);
        ErecoOverEtruer9infEE->SetLineColor(613);
        ErecoOverEtruer9supEE->SetLineColor(429);
	//413
        //397
        ErecoOverEtruer9infEB->SetMarkerColor(629);
        ErecoOverEtruer9supEB->SetMarkerColor(597);
        ErecoOverEtruer9infEE->SetMarkerColor(613);
        ErecoOverEtruer9supEE->SetMarkerColor(429);
        //413
        //397   

        ErecoOverEtruer9infEB->SetMarkerStyle(20);
        ErecoOverEtruer9supEB->SetMarkerStyle(24);
        ErecoOverEtruer9infEE->SetMarkerStyle(21);
        ErecoOverEtruer9supEE->SetMarkerStyle(25);


        ErecoOverEtruer9infEB->SetName("ErecoOverEtruer9infEB");
        ErecoOverEtruer9supEB->SetName("ErecoOverEtruer9supEB");
        ErecoOverEtruer9infEE->SetName("ErecoOverEtruer9infEE");
        ErecoOverEtruer9supEE->SetName("ErecoOverEtruer9supEE");

        text = new TLatex();
        text->SetNDC();
        text->SetTextAlign(11);
        text->SetTextFont(42);
        text->SetTextSizePixels(17);
        text->SetTextSize(0.028);
        text->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        text->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        text->DrawLatex(0.16, 0.80, "TProfile");
        text->DrawLatex(0.16, 0.75, "ETHZ Corrections");

        TLegend * leg = new TLegend(0.6,0.7,0.9,0.9,"","brNDC");
        //leg->SetHeader("The Legend Title");
        //leg->AddEntry(h1,"Histogram filled with random numbers","f");
        //leg->AddEntry("f1","Function abs(#frac{sin(x)}{x})","l");
        leg->SetTextSize(0.030);
        leg->SetFillColor(kWhite);
        leg->SetLineColor(kWhite);
        leg->SetShadowColor(kWhite);
        leg->AddEntry(ErecoOverEtruer9infEB->GetName(),"Barrel R_{9 #gamma} < 0.94","lep");
        leg->AddEntry(ErecoOverEtruer9supEB->GetName(),"Barrel R_{9 #gamma} > 0.94","lep");
        leg->AddEntry(ErecoOverEtruer9infEE->GetName(),"Endcaps R_{9 #gamma} < 0.95","lep");
        leg->AddEntry(ErecoOverEtruer9supEE->GetName(),"Endcaps R_{9 #gamma} > 0.95","lep");
        leg->Draw();

	gStyle->SetPadBorderMode(0);

        mg->GetYaxis()->SetRangeUser(yminErecoOverEtrue,ymaxErecoOverEtrue);
        mg->GetXaxis()->SetLimits(xminErecoOverEtrue,xmaxErecoOverEtrue);

        c1->SetTickx(1);
        c1->SetTicky(1);
        c1->SetGridx(1);
        c1->SetGridy(1);
        c1->Modified();
        c1->cd();
        c1->SetSelected(c1);
        c1->ToggleToolBar();

        //nomDossier = Form("SrecoVsPt/%s/2011/MmumuSelection%s_%s/WithCracks/%s/",variableX.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),SetOfCorrections.c_str());

        nomFichier = Form("SrecoVsPt_4Cat_%s",version.c_str());

        enregistrementPlots(nomDossier, nomFichier, 2, 10000, c1);
	
	return 0;

}
