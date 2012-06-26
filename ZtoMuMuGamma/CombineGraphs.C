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



int CombineGraphsF_Juin_2012(string scale = "mmg_ik_MZ_Surface", string nomFitMethode = "RooLandau2", string PileUpVersion = "2011", string LowMmumuLim = "40", string HightMmumuLim = "80", bool phiCracks = true, bool etaCracks = true, string SetOfCorrections = "ETHZCorrections", string variableX = "Photon_Et", int isMC = 1, int FitRange = 0, string MZbinning = "05GeV",string MuonCorrection = "Rochester", string SurfaceMethod = "Fitted_Surface", string Category = "Vgamma8", string date = "16Jan")
{
	//ElectronTunedCorrections
	gROOT->Reset();
	setTDRStyle();
	TGaxis::SetMaxDigits(3);

	double yminErecoOverEtrue, ymaxErecoOverEtrue;
	double xminErecoOverEtrue, xmaxErecoOverEtrue;
	
	yminErecoOverEtrue = -15.0;
	ymaxErecoOverEtrue = 15.0;
	xminErecoOverEtrue = 0.0;
	if(variableX == "Photon_SC_Eta") xminErecoOverEtrue = -3.0;
        xmaxErecoOverEtrue = 100.0;
        if(variableX == "Photon_SC_rawEt") xmaxErecoOverEtrue = 250.0;
        if(variableX == "Photon_Et") xmaxErecoOverEtrue = 100.0;
        if(variableX == "Photon_E") xmaxErecoOverEtrue = 350.0;
        if(variableX == "Photon_SC_Eta") xmaxErecoOverEtrue = 3.0;
        if(variableX == "Photon_SC_brem") xmaxErecoOverEtrue = 15.0;


	int n = 6;
	
	string compressFitMethodeName;
	if(nomFitMethode == "RooCrystalBall") compressFitMethodeName = "CB";
	if(nomFitMethode == "RooBifurcatedGauss") compressFitMethodeName = "BFG";
	if(nomFitMethode == "RooLandau2") compressFitMethodeName = "Landau";
        if(nomFitMethode == "RooLandauConvGaussian") compressFitMethodeName = "LandauConvGaus";

	string compressScaleName;
        if(scale == "Photon_E_o_MC_E") compressScaleName = "ErecoOverEtrue";
        if(scale == "ComparaisonFits1overKrecoDiffMmumuJan") compressScaleName = "1overKreco";
	if(scale == "mmg_ik_MZ_Surface") compressScaleName = "1overKreco";

	string ChainCracks = ""; 
        if(phiCracks == true && etaCracks == true) ChainCracks = "WithCracks";
        if(phiCracks == true && etaCracks == false) ChainCracks = "WithoutEtaCracks";
        if(phiCracks == false && etaCracks == true) ChainCracks = "WithoutPhiCracks";
        if(phiCracks == false && etaCracks == false) ChainCracks = "WithoutCracks";

	string isMCChain = ""; 
        if(isMC == 0) isMCChain = "Data";
        if(isMC == 1) isMCChain = "MC";


	string DossierEnregistrement = "";
	string FichierEnregistrement = "";

	string fichierTXT = "";
	string fichierxValueTXT = "";
	string fichierxErrorLTXT = "";
	string fichierxErrorRTXT = "";
	string fichierMeanTabTXT = "";
	string fichierMeanErrorTabTXT = "";


	double xValueTabR9infEB[6] = {0.0};
	double xErrorLTabR9infEB[6] = {0.0};
	double xErrorRTabR9infEB[6] = {0.0};
	double MeanTabR9infEB[6] = {0.0};
	double MeanErrorTabR9infEB[6] = {0.0};

	double xValueTabR9supEB[6] = {0.0};
        double xErrorLTabR9supEB[6] = {0.0};
        double xErrorRTabR9supEB[6] = {0.0};
        double MeanTabR9supEB[6] = {0.0};
        double MeanErrorTabR9supEB[6] = {0.0};

	double xValueTabR9infEE[6] = {0.0};
        double xErrorLTabR9infEE[6] = {0.0};
        double xErrorRTabR9infEE[6] = {0.0};
        double MeanTabR9infEE[6] = {0.0};
        double MeanErrorTabR9infEE[6] = {0.0};

	double xValueTabR9supEE[6] = {0.0};
        double xErrorLTabR9supEE[6] = {0.0};
        double xErrorRTabR9supEE[6] = {0.0};
        double MeanTabR9supEE[6] = {0.0};
        double MeanErrorTabR9supEE[6] = {0.0};

	int rangeOpt = 0;
        if(FitRange == 0) rangeOpt = 1;

	for(int i = 0; i < 4; i++)
	{

		if(rangeOpt == 1){


		////////// Results_Mai_2012_v3_SurfaceGen //////////
			
		/// BFG ///
		// S True //
/*		if(scale == "Photon_E_o_MC_E" && nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MuonCorrection == "Rochester")
                {
                        if(i == 0) FitRange = 68; //low EB 
                        if(i == 1) FitRange = 81; //high EB
                        if(i == 2) FitRange = 72; //low EE
                        if(i == 3) FitRange = 84; //high EE
                }	
		if(scale == "Photon_E_o_MC_E" && nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MuonCorrection == "NoMuonCorrection")
                {
                        if(i == 0) FitRange = 60; //low EB 
                        if(i == 1) FitRange = 82; //high EB
                        if(i == 2) FitRange = 73; //low EE
                        if(i == 3) FitRange = 84; //high EE
                }
*/

		// S Surface //
		if(scale == "mmg_ik_MZ_Surface" && nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MZbinning == "05GeV" && MuonCorrection == "Rochester" && SurfaceMethod == "Profile_Surface")
                {
                        if(i == 0) FitRange = 85; //low EB 
                        if(i == 1) FitRange = 72; //high EB
                        if(i == 2) FitRange = 81; //low EE
                        if(i == 3) FitRange = 85; //high EE
                }	
		if(scale == "mmg_ik_MZ_Surface" && nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MZbinning == "05GeV" && MuonCorrection == "NoMuonCorrection" && SurfaceMethod == "Profile_Surface")
                {
                        if(i == 0) FitRange = 84; //low EB 
                        if(i == 1) FitRange = 78; //high EB
                        if(i == 2) FitRange = 89; //low EE
                        if(i == 3) FitRange = 81; //high EE
                }
		if(scale == "mmg_ik_MZ_Surface" && nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MZbinning == "1GeV" && MuonCorrection == "Rochester" && SurfaceMethod == "Profile_Surface")
                {
                        if(i == 0) FitRange = 94; //low EB 
                        if(i == 1) FitRange = 87; //high EB
                        if(i == 2) FitRange = 90; //low EE
                        if(i == 3) FitRange = 81; //high EE
                }
                if(scale == "mmg_ik_MZ_Surface" && nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MZbinning == "1GeV" && MuonCorrection == "NoMuonCorrection" && SurfaceMethod == "Profile_Surface")
                {
                        if(i == 0) FitRange = 94; //low EB 
                        if(i == 1) FitRange = 83; //high EB
                        if(i == 2) FitRange = 94; //low EE
                        if(i == 3) FitRange = 91; //high EE
                }
		if(scale == "mmg_ik_MZ_Surface" && nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MZbinning == "2GeV" && MuonCorrection == "Rochester" && SurfaceMethod == "Profile_Surface")
                {
                        if(i == 0) FitRange = 94; //low EB 
                        if(i == 1) FitRange = 94; //high EB
                        if(i == 2) FitRange = 90; //low EE
                        if(i == 3) FitRange = 92; //high EE
                }
                if(scale == "mmg_ik_MZ_Surface" && nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MZbinning == "2GeV" && MuonCorrection == "NoMuonCorrection" && SurfaceMethod == "Profile_Surface")
                {
                        if(i == 0) FitRange = 88; //low EB 
                        if(i == 1) FitRange = 74; //high EB
                        if(i == 2) FitRange = 88; //low EE
                        if(i == 3) FitRange = 93; //high EE
                }
		if(scale == "mmg_ik_MZ_Surface" && nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MZbinning == "5GeV" && MuonCorrection == "Rochester" && SurfaceMethod == "Profile_Surface")
                {
                        if(i == 0) FitRange = 79; //low EB 
                        if(i == 1) FitRange = 86; //high EB
                        if(i == 2) FitRange = 88; //low EE
                        if(i == 3) FitRange = 92; //high EE
                }
                if(scale == "mmg_ik_MZ_Surface" && nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MZbinning == "5GeV" && MuonCorrection == "NoMuonCorrection" && SurfaceMethod == "Profile_Surface")
                {
                        if(i == 0) FitRange = 84; //low EB 
                        if(i == 1) FitRange = 82; //high EB
                        if(i == 2) FitRange = 91; //low EE
                        if(i == 3) FitRange = 92; //high EE
                }
		if(scale == "mmg_ik_MZ_Surface" && nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MZbinning == "2GeV" && MuonCorrection == "NoMuonCorrection" && SurfaceMethod == "Fitted_Surface")
                {
                        if(i == 0) FitRange = 88; //low EB 
                        if(i == 1) FitRange = 76; //high EB
                        if(i == 2) FitRange = 93; //low EE
                        if(i == 3) FitRange = 86; //high EE
                }	
		if(scale == "mmg_ik_MZ_Surface" && nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MZbinning == "2GeV" && MuonCorrection == "Rochester" && SurfaceMethod == "Fitted_Surface")
                {
                        if(i == 0) FitRange = 85; //low EB 
                        if(i == 1) FitRange = 83; //high EB
                        if(i == 2) FitRange = 94; //low EE
                        if(i == 3) FitRange = 95; //high EE
                }	
		if(scale == "mmg_ik_MZ_Surface" && nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MZbinning == "5GeV" && MuonCorrection == "NoMuonCorrection" && SurfaceMethod == "Fitted_Surface")
                {
                        if(i == 0) FitRange = 82; //low EB 
                        if(i == 1) FitRange = 77; //high EB
                        if(i == 2) FitRange = 86; //low EE
                        if(i == 3) FitRange = 82; //high EE
                }
                if(scale == "mmg_ik_MZ_Surface" && nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MZbinning == "5GeV" && MuonCorrection == "Rochester" && SurfaceMethod == "Fitted_Surface")
                {
                        if(i == 0) FitRange = 93; //low EB 
                        if(i == 1) FitRange = 72; //high EB
                        if(i == 2) FitRange = 93; //low EE
                        if(i == 3) FitRange = 82; //high EE
                }


                }


		//fichierTXT = Form("Results_Avril_2012_v1/%s/%s/MmumuSelection%s_%s/%s/%s/%s/%s/MZbinning_%s/%dPourcents/Fits/",variableX.c_str(),PileUpVersion.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),ChainCracks.c_str(),SetOfCorrections.c_str(),MuonCorrection.c_str(),isMCChain.c_str(),MZbinning.c_str(),FitRange);
		//fichierTXT = Form("Results_Avril_2012_v2/%s/%s/MmumuSelection%s_%s/%s/%s/%s/%s/MZbinning_%s/%s/%dPourcents/Fits/",variableX.c_str(),PileUpVersion.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),ChainCracks.c_str(),SetOfCorrections.c_str(),MuonCorrection.c_str(),isMCChain.c_str(),MZbinning.c_str(),SurfaceMethod.c_str(),FitRange);
		fichierTXT = Form("Results_Mai_2012_v3_SurfaceGen/%s/%s/MmumuSelection%s_%s/%s/%s/%s/%s/MZbinning_%s/%s/%s/%s/%dPourcents/Fits/",variableX.c_str(),PileUpVersion.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),ChainCracks.c_str(),SetOfCorrections.c_str(),MuonCorrection.c_str(),isMCChain.c_str(),MZbinning.c_str(),SurfaceMethod.c_str(),Category.c_str(),date.c_str(),FitRange);


                fichierTXT += Form("%s%s",compressScaleName.c_str(),compressFitMethodeName.c_str());



		if(i == 0) fichierTXT += Form("r9infEB/",compressFitMethodeName.c_str());
		if(i == 1) fichierTXT += Form("r9supEB/",compressFitMethodeName.c_str());
		if(i == 2) fichierTXT += Form("r9infEE/",compressFitMethodeName.c_str());
		if(i == 3) fichierTXT += Form("r9supEE/",compressFitMethodeName.c_str());


		fichierxValueTXT = fichierTXT + "xValue.txt";
		fichierxErrorLTXT = fichierTXT + "xErrorL.txt";
		fichierxErrorRTXT = fichierTXT + "xErrorR.txt";
		fichierMeanTabTXT = fichierTXT + "MeanTab.txt";
		fichierMeanErrorTabTXT = fichierTXT + "MeanErrorTab.txt";


        	cout<<"fichierTXT = "<<fichierTXT<<endl;
        	ifstream monFlux1(fichierxValueTXT.c_str());	
		ifstream monFlux2(fichierxErrorLTXT.c_str());
		ifstream monFlux3(fichierxErrorRTXT.c_str());
		ifstream monFlux4(fichierMeanTabTXT.c_str());
		ifstream monFlux5(fichierMeanErrorTabTXT.c_str());

		for(int j = 0; j < n; j++)
		{ 

			if(i == 0)
			{
				cout<<"fichierxValueTXT = "<<fichierxValueTXT<<endl;
				monFlux1 >> xValueTabR9infEB[j];
				cout<<"xValueTabR9infEB[j] = "<<xValueTabR9infEB[j]<<endl;
				monFlux2 >> xErrorLTabR9infEB[j];
				monFlux3 >> xErrorRTabR9infEB[j];
				monFlux4 >> MeanTabR9infEB[j];
				monFlux5 >> MeanErrorTabR9infEB[j];
			}
			
			if(i == 1)
                        {
                                monFlux1 >> xValueTabR9supEB[j];
                                monFlux2 >> xErrorLTabR9supEB[j];
                                monFlux3 >> xErrorRTabR9supEB[j];
                                monFlux4 >> MeanTabR9supEB[j];
                                monFlux5 >> MeanErrorTabR9supEB[j];
                        }

			if(i == 2)
                        {
                                monFlux1 >> xValueTabR9infEE[j];
                                monFlux2 >> xErrorLTabR9infEE[j];
                                monFlux3 >> xErrorRTabR9infEE[j];
                                monFlux4 >> MeanTabR9infEE[j];
                                monFlux5 >> MeanErrorTabR9infEE[j];
                        }

			if(i == 3)
                        {
                                monFlux1 >> xValueTabR9supEE[j];
                                monFlux2 >> xErrorLTabR9supEE[j];
                                monFlux3 >> xErrorRTabR9supEE[j];
                                monFlux4 >> MeanTabR9supEE[j];
                                monFlux5 >> MeanErrorTabR9supEE[j];
                        }
	
		}	

	}


	TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);

	TMultiGraph *mg = new TMultiGraph();


	TGraphAsymmErrors * ErecoOverEtruer9infEB = new TGraphAsymmErrors(n,xValueTabR9infEB, MeanTabR9infEB, xErrorLTabR9infEB, xErrorRTabR9infEB, MeanErrorTabR9infEB, MeanErrorTabR9infEB);

	TGraphAsymmErrors * ErecoOverEtruer9supEB = new TGraphAsymmErrors(n,xValueTabR9supEB, MeanTabR9supEB, xErrorLTabR9supEB, xErrorRTabR9supEB, MeanErrorTabR9supEB, MeanErrorTabR9supEB);

	TGraphAsymmErrors * ErecoOverEtruer9infEE = new TGraphAsymmErrors(n,xValueTabR9infEE, MeanTabR9infEE, xErrorLTabR9infEE, xErrorRTabR9infEE, MeanErrorTabR9infEE, MeanErrorTabR9infEE);

	TGraphAsymmErrors * ErecoOverEtruer9supEE = new TGraphAsymmErrors(n,xValueTabR9supEE, MeanTabR9supEE, xErrorLTabR9supEE, xErrorRTabR9supEE, MeanErrorTabR9supEE, MeanErrorTabR9supEE);


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

	//mg->GetXaxis()->SetTitle("Pt_{#gamma}");
        //mg->GetXaxis()->SetTitle("P_{T}^{#gamma}");
	
	if(variableX == "Photon_SC_rawEt") mg->GetXaxis()->SetTitle("P_{T RAW}^{#gamma}");
        if(variableX == "Photon_Et") mg->GetXaxis()->SetTitle("P_{T}^{#gamma}");
        if(variableX == "Photon_E") mg->GetXaxis()->SetTitle("E^{#gamma}");
        if(variableX == "Photon_SC_Eta") mg->GetXaxis()->SetTitle("#eta");
        if(variableX == "Photon_SC_brem") mg->GetXaxis()->SetTitle("#sigma_{#phi}/#sigma_{#eta}");

	mg->GetXaxis()->SetLabelFont(42);
        mg->GetXaxis()->SetTitleFont(42);
        mg->GetXaxis()->SetLabelSize(0.03);
        //mg->GetYaxis()->SetTitle("E_{RECO}/E_{TRUE} - 1 (%)");
       	if(scale == "Photon_E_o_MC_E") mg->GetYaxis()->SetTitle("s_{TRUE} = E_{RECO}/E_{TRUE} - 1 (%)");
        if(scale == "ComparaisonFits1overKrecoDiffMmumuJan") mg->GetYaxis()->SetTitle("s_{RECO} = 1/k_{RECO} - 1 (%)");
	if(scale == "mmg_ik_MZ_Surface") mg->GetYaxis()->SetTitle("s_{RECO Surface} (%)");
	mg->GetYaxis()->SetLabelFont(42);
        mg->GetYaxis()->SetTitleOffset(1.24);
        mg->GetYaxis()->SetTitleFont(42);
        mg->GetYaxis()->SetLabelSize(0.03);
	//mg->SetTitle("");
        //mg->SetMarkerColor(4);
        //mg->SetMarkerStyle(21);
        //mg->SetMarkerSize(0.6);




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


	TLatex *text = new TLatex();

        text = new TLatex();
        text->SetNDC();
        text->SetTextAlign(11);
        text->SetTextFont(42);
        text->SetTextSizePixels(17);
	text->SetTextSize(0.028);
	text->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	if(isMC == 1) text->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) text->DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4,66892 fm^{-1}");
	if(nomFitMethode == "RooCrystalBall") text->DrawLatex(0.16, 0.80, "Crystal Ball");
        if(nomFitMethode == "RooBifurcatedGauss") text->DrawLatex(0.16, 0.80, "Bifurcated Gaussian");
	if(nomFitMethode == "RooLandau2") text->DrawLatex(0.16, 0.80, "Landau");
	if(SetOfCorrections == "ETHZCorrections") text->DrawLatex(0.16, 0.75, "ETHZ Corrections");
	if(SetOfCorrections == "ElectronTunedCorrections") text->DrawLatex(0.16, 0.75, "ElectronTuned Corrections");

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

	//DossierEnregistrement = Form("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsErecoOverEtrueDiffMmumuJan/Olivier9/%s/MmumuSelection%s_%s/Fits/",PileUpVersion.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str());
	

	if(rangeOpt == 1){


                //DossierEnregistrement = Form("Results_Avril_2012_v1/%s/%s/MmumuSelection%s_%s/%s/%s/%s/%s/MZbinning_%s/CombinedGraphs/",variableX.c_str(),PileUpVersion.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),ChainCracks.c_str(),SetOfCorrections.c_str(),MuonCorrection.c_str(),isMCChain.c_str(),MZbinning.c_str());
		//DossierEnregistrement = Form("Results_Avril_2012_v2/%s/%s/MmumuSelection%s_%s/%s/%s/%s/%s/MZbinning_%s/%s/CombinedGraphs/",variableX.c_str(),PileUpVersion.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),ChainCracks.c_str(),SetOfCorrections.c_str(),MuonCorrection.c_str(),isMCChain.c_str(),MZbinning.c_str(),SurfaceMethod.c_str());
		DossierEnregistrement = Form("Results_Mai_2012_v3_SurfaceGen/%s/%s/MmumuSelection%s_%s/%s/%s/%s/%s/MZbinning_%s/%s/%s/%s/CombinedGraphs/",variableX.c_str(),PileUpVersion.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),ChainCracks.c_str(),SetOfCorrections.c_str(),MuonCorrection.c_str(),isMCChain.c_str(),MZbinning.c_str(),SurfaceMethod.c_str(),Category.c_str(),date.c_str());



        }

	//FichierEnregistrement = Form("ErecoOverEtrueVsPt%s",compressFitMethodeName.c_str());

	FichierEnregistrement = Form("%sVsPt%s_%s",compressScaleName.c_str(),compressFitMethodeName.c_str(),SetOfCorrections.c_str());

        enregistrementPlots(DossierEnregistrement, FichierEnregistrement, 2, 10000, c1);




	return 0;
}
