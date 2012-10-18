// -------------------------------------------------------
// program implemented by Louis Sgandurra (October 2012)
// -------------------------------------------------------

#ifndef __CINT__
#include "RooGlobalFunc.h"
#else
// Refer to a class implemented in libRooFit to force its loading
// via the autoloader.
class Roo2DKeysPdf;
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TMath.h"
#include "TF1.h"
#include "Math/DistFunc.h"
#ifndef __CINT__
#include "RooCFunction1Binding.h" 
#include "RooCFunction3Binding.h"
#endif
#include "RooTFnBinding.h"
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
#include "TSystem.h"
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
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <iomanip>
#include "RooAbsPdf.h"
#include "RooCBShape.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooMsgService.h"
#include "RooNLLVar.h"
#include "RooRandom.h"
#include "RooWorkspace.h"
#include "RooFFTConvPdf.h"
#include "RooLognormal.h"
#include "RooGamma.h"
#include "RooBifurGauss"
#include "RooAddPdf"
#include "RooGenericPdf"
#include "RooProdPdf"
#include "TROOT.h"
#include "TRint.h"
#include "TMultiGraph.h"
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
#include <sstream>
#pragma optimize 0

using namespace RooFit;
using namespace std;

void enregistrementPlots(string nomDossier, string nomFichier, int EndCaps, int iteration, TCanvas * c1);


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


int SameStatBins(string eta = "EB_EE", string r9 = "AllR9", string variable = "nVertices", int nBins = 100)
{
	TChain * chain = new TChain("miniTree");
	chain->Add("/sps/cms/sgandurr/CMSSW_5_2_5_patch3/src/Zmumugamma/Selection/miniTree_DoubleMu_Run2012B-Zmmg-PromptSkim-v1_13JulyJSON_xi0_RECO_5_2_4_v3_NoMuonCor_40_80_v1_partALL.root");
	TString cuts = "isLooseMMG == 1";	
	
	if(eta == "EB" && r9 == "highR9") cuts += " && (abs(Photon_SC_Eta)) < 1.4442 && Photon_r9 > 0.94";
        if(eta == "EB" && r9 == "lowR9") cuts += " && (abs(Photon_SC_Eta)) < 1.4442 && Photon_r9 < 0.94";
        if(eta == "EE" && r9 == "highR9") cuts += " && (abs(Photon_SC_Eta)) > 1.566 && Photon_r9 > 0.95";
        if(eta == "EE" && r9 == "lowR9") cuts += " && (abs(Photon_SC_Eta)) > 1.566 && Photon_r9 < 0.95";
        if(eta == "EB" && r9 == "AllR9") cuts += " && (abs(Photon_SC_Eta)) < 1.4442";
        if(eta == "EE" && r9 == "AllR9") cuts += " && (abs(Photon_SC_Eta)) > 1.566";
        if(eta == "EB_EE" && r9 == "AllR9") cuts = "isLooseMMG == 1";	


	TChain * ReducedChain = (TChain *) chain->CopyTree(cuts);

	float nVertices;
        float Photon_Et;
	float Photon_E;
	float Photon_SC_brem;
	float Photon_SC_Eta;
	float Photon_SC_rawEt;
	
	ReducedChain->SetBranchAddress("nVertices",&nVertices);
	ReducedChain->SetBranchAddress("Photon_Et",&Photon_Et);
	ReducedChain->SetBranchAddress("Photon_E",&Photon_E);
	ReducedChain->SetBranchAddress("Photon_SC_brem",&Photon_SC_brem);
	ReducedChain->SetBranchAddress("Photon_SC_Eta",&Photon_SC_Eta);
	ReducedChain->SetBranchAddress("Photon_SC_rawEt",&Photon_SC_rawEt);

	vector <double> variableTab;

	cout<<endl<<"ReducedChain->GetEntries() = "<<ReducedChain->GetEntries()<<endl;
	if(variable == "nVertices")
	{
		for (int ievt = 0 ; ievt < ReducedChain->GetEntries() ; ievt++)
		{
			ReducedChain->GetEntry(ievt);
			variableTab.push_back(nVertices);
		}
	}
	if(variable == "Photon_Et")
        {
                for (int ievt = 0 ; ievt < ReducedChain->GetEntries() ; ievt++)
                {
			ReducedChain->GetEntry(ievt);
                        variableTab.push_back(Photon_Et);
                }
        }
	if(variable == "Photon_E")
        {
                for (int ievt = 0 ; ievt < ReducedChain->GetEntries() ; ievt++)
                {
			ReducedChain->GetEntry(ievt);
                        variableTab.push_back(Photon_E);
                }
        }
	if(variable == "Photon_SC_brem")
        {
                for (int ievt = 0 ; ievt < ReducedChain->GetEntries() ; ievt++)
                {
			ReducedChain->GetEntry(ievt);
                        variableTab.push_back(Photon_SC_brem);
                }
        }
	if(variable == "Photon_SC_Eta")
        {
                for (int ievt = 0 ; ievt < ReducedChain->GetEntries() ; ievt++)
                {
			ReducedChain->GetEntry(ievt);
                        variableTab.push_back(Photon_SC_Eta);
                }
        }
	if(variable == "Photon_SC_rawEt")
        {
                for (int ievt = 0 ; ievt < ReducedChain->GetEntries() ; ievt++)
                {
			ReducedChain->GetEntry(ievt);
                        variableTab.push_back(Photon_SC_rawEt);
                }
        }	


        sort(variableTab.begin(), variableTab.end());

	int entriesPerBins = variableTab.size() / nBins;

	cout<<endl<<"entriesPerBins = "<<entriesPerBins<<endl;

	ofstream fichier(Form("Limites_%s_%s%s_%dbins.txt",variable.c_str(),r9.c_str(),eta.c_str(),nBins));

	fichier<<0.0<<endl;	

	for(int i = 1; i <= nBins; i++)
	{
		cout<<endl<<"i = "<<i;
		if(i != nBins) fichier<<variableTab[i*entriesPerBins -1]<<endl;
		if(i == nBins) fichier<<variableTab[variableTab.size() - 1];
	}

	
	fichier.close();	


	return 0;
}	













