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
#include "TProfile3D.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <iomanip>
#include "RooProduct.h"
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


int SurfaceGenerator(string type = "MC", string MuonCorrection = "NoMuonCorrection", string date = "Fall11", string choix = "plot3D", string SurfaceBinning = "5GeV", string Fits = "0", string FitMethod = "RooCrystalBall", int Gen = 1, int density = 1)
{
	//int Gen = 0;
        //int Gen = 1;
	
	//string type = "Data"; 
	//string type = "MC";

	//string MuonCorrection = "NoMuonCor";
	//string MuonCorrection = "RochCor";

	//string date = "16Jan2012";
	//string date = "30Nov2011";
	//string date = "Vgamma";
	
	if(type == "MC") date = "Fall11";

	//string choix = "plot3D";
	//string choix = "plot2D";
	//string choix = "Mmumu.txt";
	
	//string SurfaceBinning = "05GeV";
	//string SurfaceBinning = "1GeV";
	//string SurfaceBinning = "2GeV";
	//string SurfaceBinning = "5GeV";


	int nbOfBins = 0;
        if(SurfaceBinning == "05GeV") nbOfBins = 400;
        if(SurfaceBinning == "1GeV") nbOfBins = 200;
        if(SurfaceBinning == "2GeV") nbOfBins = 100;
        if(SurfaceBinning == "5GeV") nbOfBins = 40;

	string nomDossier = "";
	string nomDossierTXT = "";
	string nomFichier = "";
	string nomDossierSurface = "";

	if(Gen == 0)
        {

		if(Fits == "1")
		{
			
			nomDossier = Form("SurfacePlots_v2/%s/%s/%s/Surface_%s/%s/",type.c_str(),MuonCorrection.c_str(),date.c_str(),SurfaceBinning.c_str(),FitMethod.c_str());
	        	nomDossierTXT = nomDossier + "TXT/";
			nomDossierSurface = nomDossier;
			nomFichier = "Surface_Fit"; 
		}
	
		if(Fits == "0")
	        {
			nomDossier = "SurfacePlots_v1b/"; 
	        	nomFichier = Form("Surface_%s_%s_%s_%s",type.c_str(),MuonCorrection.c_str(),date.c_str(),SurfaceBinning.c_str());
		}
	}

	if(Gen == 1)
	{
		if(Fits == "0")
		{
			nomDossier = Form("SurfacePlots_MMGgen_v1/%s/%s/%s/Surface_%s/Profile_Surface/",type.c_str(),MuonCorrection.c_str(),date.c_str(),SurfaceBinning.c_str());
			nomFichier = "Surface_Profile";	
		}	
		if(Fits == "1")
                {
                        nomDossier = Form("SurfacePlots_MMGgen_v1/%s/%s/%s/Surface_%s/Fitted_Surface/%s/",type.c_str(),MuonCorrection.c_str(),date.c_str(),SurfaceBinning.c_str(),FitMethod.c_str());
			nomDossierTXT = nomDossier + "TXT/";
                	nomDossierSurface = nomDossier;
                	nomFichier = "Surface_Fit";
                }


	}


	gROOT->Reset();
        TGaxis::SetMaxDigits(3);
        setTDRStyle();
	gStyle->SetPalette(1,0);
	//gStyle->SetCanvasDefW(600);


	TCanvas *c2 = new TCanvas("c2", "c2",0,0,600,600);

	TChain *chain = new TChain("miniTree");
	
	if(Gen == 0)
	{	
		if(type == "MC")
		{
			if(MuonCorrection == "NoMuonCor") chain->Add("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/Selection/miniTreeMuons_v08_Zmumu_partALL.root");
			if(MuonCorrection == "RochCor") chain->Add("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/Selection/miniTreeMuons_v09_Zmumu_rochcor_partALL.root"); 
	
		}
		if(type == "Data")
	        {
	                if(MuonCorrection == "NoMuonCor")
			{
				if(date == "16Jan2012")
				{
					chain->Add("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/Selection/miniTreeMuons_v10_Zmumu_Run2011A-16Jan2012-v1_partALL.root");
					chain->Add("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/Selection/miniTreeMuons_v10_Zmumu_Run2011B-16Jan2012-v1_partALL.root");
				}
				if(date == "30Nov2011")
	                        {
					chain->Add("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/Selection/miniTreeMuons_v10_Zmumu_Run2011A-30Nov2011-v1_partALL.root");
					chain->Add("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/Selection/miniTreeMuons_v10_Zmumu_Run2011B-30Nov2011-v1_partALL.root");
				}
				if(date == "Vgamma")
				{
					chain->Add("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/Selection/miniTreeMuons_v10_Zmumu_Run2011A-ZMu-PromptSkim-v4_partALL.root");
					chain->Add("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/Selection/miniTreeMuons_v10_Zmumu_Run2011A-ZMu-May10ReReco-v1_partALL.root");
					chain->Add("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/Selection/miniTreeMuons_v10_Zmumu_Run2011A-ZMu-05Aug2011-v1_V04_partALL.root");
					chain->Add("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/Selection/miniTreeMuons_v10_Zmumu_Run2011A-ZMu-03Oct2011-v1_partALL.root");
					chain->Add("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/Selection/miniTreeMuons_v10_Zmumu_Run2011B-ZMu-PromptSkim-v1_finalJson_partALL.root");
	
				}
	
			}
	
			if(MuonCorrection == "RochCor")
	                {
	
				if(date == "16Jan2012")
	                        {
	                                chain->Add("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/Selection/miniTreeMuons_v10_Zmumu_Run2011A-16Jan2012-v1_rochcor_partALL.root");
	                                chain->Add("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/Selection/miniTreeMuons_v10_Zmumu_Run2011B-16Jan2012-v1_rochcor_partALL.root");
	                        }
	                        if(date == "30Nov2011")
	                        {
	                                chain->Add("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/Selection/miniTreeMuons_v10_Zmumu_Run2011A-30Nov2011-v1_rochcor_partALL.root");
	                                chain->Add("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/Selection/miniTreeMuons_v10_Zmumu_Run2011B-30Nov2011-v1_rochcor_partALL.root");
	                        }
	                        if(date == "Vgamma")
	                        {
	                                chain->Add("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/Selection/miniTreeMuons_v10_Zmumu_Run2011A-ZMu-PromptSkim-v4_rochcor_partALL.root");
	                                chain->Add("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/Selection/miniTreeMuons_v10_Zmumu_Run2011A-ZMu-May10ReReco-v1_rochcor_partALL.root");
	                                chain->Add("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/Selection/miniTreeMuons_v10_Zmumu_Run2011A-ZMu-05Aug2011-v1_V04_rochcor_partALL.root");
	                                chain->Add("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/Selection/miniTreeMuons_v10_Zmumu_Run2011A-ZMu-03Oct2011-v1_rochcor_partALL.root");
	                                chain->Add("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/Selection/miniTreeMuons_v10_Zmumu_Run2011B-ZMu-PromptSkim-v1_finalJson_rochcor_partALL.root");
	
	                        }
	
	
	                }
		}
	}

	string compressSetOfCorrectionsChain = "ETHZ";
	string LowMmumuLim = "40";
	string HightMmumuLim = "80";
	if(Gen == 1)
        {
		if(type == "MC")
                {
			chain->Add(Form("/sps/cms/sgandurr/CMSSW_4_2_8/src/UserCode/IpnTreeProducer/miniTree_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_PU_S6_v2_%s_%s_%s_%s_%s_1_v3_partALL.root",compressSetOfCorrectionsChain.c_str(),MuonCorrection.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),SurfaceBinning.c_str()));
        		chain->Add(Form("/sps/cms/sgandurr/CMSSW_4_2_8/src/UserCode/IpnTreeProducer/miniTree_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_PU_S6_v2_%s_%s_%s_%s_%s_2_v3_partALL.root",compressSetOfCorrectionsChain.c_str(),MuonCorrection.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),SurfaceBinning.c_str()));
		cout<<endl<<"coucou"<<endl;
		}

	}

	if(choix == "plot3D" || choix == "plot2D")
	{

		if(SurfaceBinning == "2GeV")
		{
			TH2D *h__ = new TH2D("h__", "h__", 50, 0, 100, 50, 0, 100);
        		TH2D *h2__ = new TH2D("h2__", "h2__", 50, 0, 100, 50, 0, 100);
		}
		if(SurfaceBinning == "05GeV")
                {
        		TH2D *h__ = new TH2D("h__", "h__", 200, 0, 100, 200, 0, 100);
        		TH2D *h2__ = new TH2D("h2__", "h2__", 200, 0, 100, 200, 0, 100);      
		}
		if(SurfaceBinning == "1GeV")
                {
        		TH2D *h__ = new TH2D("h__", "h__", 100, 0, 100, 100, 0, 100);
        		TH2D *h2__ = new TH2D("h2__", "h2__", 100, 0, 100, 100, 0, 100);      
		}
		if(SurfaceBinning == "5GeV")
                {
        		TH2D *h__ = new TH2D("h__", "h__", 20, 0, 100, 20, 0, 100);
        		TH2D *h2__ = new TH2D("h2__", "h2__", 20, 0, 100, 20, 0, 100);	
		}
	}

	if(choix == "Mmumu.txt")
        {

                if(SurfaceBinning == "2GeV")
                {

			TH2D *h__ = new TH2D("h__", "h__", 100, 0, 200, 100, 0, 200);
			TH2D *h2__ = new TH2D("h2__", "h2__", 100, 0, 200, 100, 0, 200);
		}
		if(SurfaceBinning == "05GeV")
                {
			TH2D *h__ = new TH2D("h__", "h__", 400, 0, 200, 400, 0, 200);
       			TH2D *h2__ = new TH2D("h2__", "h2__", 400, 0, 200, 400, 0, 200);	
		}
		if(SurfaceBinning == "1GeV")
                {
			TH2D *h__ = new TH2D("h__", "h__", 200, 0, 200, 200, 0, 200);
        		TH2D *h2__ = new TH2D("h2__", "h2__", 200, 0, 200, 200, 0, 200);	
		}
		if(SurfaceBinning == "5GeV")
                {
			TH2D *h__ = new TH2D("h__", "h__", 40, 0, 200, 40, 0, 200);
        		TH2D *h2__ = new TH2D("h2__", "h2__", 40, 0, 200, 40, 0, 200);	
		}
	}


	if(Fits == "0")
        {

	if(Gen == 0) chain->Draw("MuonL_Pt:MuonS_Pt>>h__", "Mmumu*(MuonL_Pt < 200 && MuonS_Pt < 200 && isMM == 1 && Mmumu > 80.0 && Mmumu < 100.0)", "LEGO2");

	//if(Gen == 1) chain->Draw("sqrt(MuonL_MC_Px*MuonL_MC_Px + MuonL_MC_Py*MuonL_MC_Py):sqrt(MuonS_MC_Px*MuonS_MC_Px + MuonS_MC_Py*MuonS_MC_Py)>>h__", "Mmumu_Muons_MC*(sqrt(MuonL_MC_Px*MuonL_MC_Px + MuonL_MC_Py*MuonL_MC_Py) < 200 && sqrt(MuonS_MC_Px*MuonS_MC_Px + MuonS_MC_Py*MuonS_MC_Py) < 200)", "LEGO2");
	if(Gen == 1) chain->Draw("((sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2))>MuonF_Pt)?sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2)):MuonF_Pt):((sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2))<MuonF_Pt)?sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2)):MuonF_Pt) >>h__", "Mmumugamma_MMG_MC*(sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2)) < 200 &&  MuonF_Pt < 200 && isTightMMG)", "LEGO2");

	TH2D *_h = (TH2D*)gDirectory->Get("h__");

	if(Gen == 0) chain->Draw("MuonL_Pt:MuonS_Pt>>h2__", "(MuonL_Pt < 200 && MuonS_Pt < 200 && isMM == 1 && Mmumu > 80.0 && Mmumu < 100.0)", "LEGO2");
	
	//if(Gen == 1) chain->Draw("sqrt(MuonL_MC_Px*MuonL_MC_Px + MuonL_MC_Py*MuonL_MC_Py):sqrt(MuonS_MC_Px*MuonS_MC_Px + MuonS_MC_Py*MuonS_MC_Py)>>h2__", "(sqrt(MuonL_MC_Px*MuonL_MC_Px + MuonL_MC_Py*MuonL_MC_Py) < 200 && sqrt(MuonS_MC_Px*MuonS_MC_Px + MuonS_MC_Py*MuonS_MC_Py) < 200)", "LEGO2");
	if(Gen == 1) chain->Draw("((sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2))>MuonF_Pt)?sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2)):MuonF_Pt):((sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2))<MuonF_Pt)?sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2)):MuonF_Pt) >>h2__", "(sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2)) < 200 && MuonF_Pt < 200 && isTightMMG)", "LEGO2"); 

	TH2D *_h2 = (TH2D*)gDirectory->Get("h2__");

	_h->Divide(_h2); // Comme ça on divise le bin weighté avec Mmumu par lenombre d'entrées dans le bin : on obtient la moyenne des Mmumu pour cebin

	}

	string nomFichierDensity = nomFichier;
	if(density == 1)
	{
		if(choix == "plot2D") nomFichierDensity += "_Density_2D";
		if(choix == "plot3D") nomFichierDensity += "_Density_3D";

		if(choix == "plot2D") h2__->Draw("COLZ"); 
		if(choix == "plot3D") h2__->Draw("LEGO2");
		h2__->GetYaxis()->SetTitle("P_{T_{#mu leading}}");
		h2__->GetXaxis()->SetTitle("P_{T_{#mu trailing}}");
		if(choix == "plot3D") h2__->GetZaxis()->SetTitle("events");
		h2__->GetXaxis()->SetLabelFont(42);
		if(choix == "plot2D") h2__->GetXaxis()->SetTitleOffset(1.5);
		if(choix == "plot3D") _h->GetXaxis()->SetTitleOffset(1.8);
		h2__->GetXaxis()->SetTitleFont(42);
        	h2__->GetXaxis()->SetLabelSize(0.03);
		h2__->GetYaxis()->SetLabelFont(42);
		if(choix == "plot2D") h2__->GetYaxis()->SetTitleOffset(1.5);
		if(choix == "plot3D") _h->GetYaxis()->SetTitleOffset(1.8);
		h2__->GetYaxis()->SetTitleFont(42);
		h2__->GetYaxis()->SetLabelSize(0.03);
		h2__->GetZaxis()->SetLabelFont(42);
		h2__->GetZaxis()->SetTitleFont(42);
        	h2__->GetZaxis()->SetLabelSize(0.03);
		if(choix == "plot3D") _h->GetZaxis()->SetTitleOffset(1.5);
        	//h2__->GetZaxis()->SetRangeUser(87,97);
        	//h2__->GetZaxis()->SetRangeUser(70,110); 

		if(choix == "plot2D")
        	{
	                TPaletteAxis *densityPalette = new TPaletteAxis(96,0,100,100,h2__);
	                densityPalette->SetLabelColor(1);
	                densityPalette->SetLabelFont(42);
	                densityPalette->SetLabelOffset(0.007);
	                densityPalette->SetLabelSize(0.02);
	                densityPalette->SetTitleOffset(1);
	                densityPalette->SetTitleSize(0.02);
	                densityPalette->SetFillColor(100);
	                densityPalette->SetFillStyle(1001);
	                h2__->GetListOfFunctions()->Add(densityPalette,"br");
                }

		enregistrementPlots(nomDossier, nomFichierDensity, 2, 10000, c2);

		
	}
	
	c2->Clear();

	double content = 0;


	if(Fits == "1")
        {
		
		if(choix == "plot3D" || choix == "plot2D")
        	{

                	if(SurfaceBinning == "2GeV")
                	{
                        	TH2D *_h = new TH2D("_h", "_h", 50, 0, 100, 50, 0, 100);
                	}
                	if(SurfaceBinning == "05GeV")
                	{
                        	TH2D *_h = new TH2D("_h", "_h", 200, 0, 100, 200, 0, 100);
               	 	}
                	if(SurfaceBinning == "1GeV")
                	{
                        	TH2D *_h = new TH2D("_h", "_h", 100, 0, 100, 100, 0, 100);
                	}
                	if(SurfaceBinning == "5GeV")
                	{
                        	TH2D *_h = new TH2D("_h", "_h", 20, 0, 100, 20, 0, 100);
                	}
		}

		ifstream monFlux3(Form("%sMmumu.txt",nomDossierTXT.c_str()));
		if(Gen == 0) chain->Draw("MuonL_Pt:MuonS_Pt>>_h","(MuonL_Pt < 200 && MuonS_Pt < 200 && isMM == 1 && Mmumu > 80.0 && Mmumu < 100.0)");
		if(Gen == 1) chain->Draw("((sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2))>MuonF_Pt)?sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2)):MuonF_Pt):((sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2))<MuonF_Pt)?sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2)):MuonF_Pt) >>h__", "(sqrt(pow(MuonN_Px + Photon_Px,2) + pow(MuonN_Py + Photon_Py,2)) < 200 &&  MuonF_Pt < 200 && isTightMMG)");
		for(int i = 0; i < nbOfBins; i++)
		{
			for(int j = 0; j < nbOfBins; j++)
			{

				monFlux3>>content;
				if(_h->GetBinContent(j,i) != 0)
				{
					_h->SetBinContent(j,i,content);
				}

			}
		}		

		monFlux3.close();
	}

	if(choix == "plot3D" || choix == "Mmumu.txt") _h->Draw("LEGO2");
	if(choix == "plot2D") _h->Draw("COLZ");	
	//_h->Draw("SURF1");
	//_h->Draw("SURF2");

	_h->GetYaxis()->SetTitle("P_{T_{#mu leading}}");
	_h->GetXaxis()->SetTitle("P_{T_{#mu trailing}}");
	if(Gen == 0) _h->GetZaxis()->SetTitle("M_{#mu#mu_{RECO}}");
	if(Gen == 1) _h->GetZaxis()->SetTitle("M_{#mu#mu_{GEN}}");

	_h->GetXaxis()->SetLabelFont(42);
        if(choix == "plot2D") _h->GetXaxis()->SetTitleOffset(1.5);
	if(choix == "plot3D") _h->GetXaxis()->SetTitleOffset(1.8);
        _h->GetXaxis()->SetTitleFont(42);
        _h->GetXaxis()->SetLabelSize(0.03);

	_h->GetYaxis()->SetLabelFont(42);
        if(choix == "plot2D") _h->GetYaxis()->SetTitleOffset(1.5);
	if(choix == "plot3D") _h->GetYaxis()->SetTitleOffset(1.8);
        _h->GetYaxis()->SetTitleFont(42);
        _h->GetYaxis()->SetLabelSize(0.03);
	
	_h->GetZaxis()->SetLabelFont(42);
        if(choix == "plot3D") _h->GetZaxis()->SetTitleOffset(1.4);
        _h->GetZaxis()->SetTitleFont(42);
        _h->GetZaxis()->SetLabelSize(0.03);
	_h->GetZaxis()->SetRangeUser(87,97);
	//_h->GetZaxis()->SetRangeUser(70,110);	


	if(choix == "plot2D")
	{
		TPaletteAxis *palette = new TPaletteAxis(96,0,100,100,_h);
        	palette->SetLabelColor(1);
        	palette->SetLabelFont(42);
        	palette->SetLabelOffset(0.007);
        	palette->SetLabelSize(0.02);
        	palette->SetTitleOffset(1);
        	palette->SetTitleSize(0.02);
        	palette->SetFillColor(100);
        	palette->SetFillStyle(1001);
        	_h->GetListOfFunctions()->Add(palette,"br");

	}


	if(choix == "plot3D" || choix == "plot2D")
	{	
		if(choix == "plot3D") nomFichier += "_3D";
		if(choix == "plot2D") nomFichier += "_2D";
		enregistrementPlots(nomDossier, nomFichier, 2, 10000, c2);	
	}

	//string nameMmumu = Form("SurfacePlots_v1b/MmumuTXT/Mmumu_%s_%s_%s_%s.txt",type.c_str(),MuonCorrection.c_str(),date.c_str(),SurfaceBinning.c_str());
	//string nameMmumu = Form("Mmumu_%s_%s_%s_%s.txt",type.c_str(),MuonCorrection.c_str(),date.c_str(),SurfaceBinning.c_str());
	string nameMmumu = "Mmumu.txt";
	cout<<endl<<"nameMmumu = "<<nameMmumu<<endl;
	//ofstream monFlux(Form("%s",nameMmumu.c_str()));
	ofstream monFlux(nameMmumu.c_str());
	//ofstream monFlux("MmumuTest.txt");	

	if(choix == "Mmumu.txt")
	{
		
		for(int i = 0; i <nbOfBins; i++)
        	{

                	for(int j = 0; j <nbOfBins; j++)
                	{
			
				if(_h->GetBinContent(j,i) != 0) monFlux<<_h->GetBinContent(j,i)<<endl;
				else monFlux<<91.187<<endl;
				//cout<<endl<<j<<" "<<i<<" "<<_h->GetBinContent(j,i);
			}
			if(i%10 == 0) cout<<"i = "<<i<<endl; 
		}	


	}
	monFlux.close();

	if(Gen == 0)
	{
		if(choix == "Mmumu.txt") system("mkdir SurfacePlots_v1b/MmumuTXT/");
		if(choix == "Mmumu.txt") system(Form("mv %s SurfacePlots_v1b/MmumuTXT/%s", nameMmumu.c_str(),nameMmumu.c_str()));
	}
	if(Gen == 1)
        {
		if(choix == "Mmumu.txt") system(Form("mkdir %sTXT/",nomDossier.c_str()));
                if(choix == "Mmumu.txt") system(Form("mv %s %sTXT/%s", nameMmumu.c_str(),nomDossier.c_str(),nameMmumu.c_str()));

	}
	return 0;

}	





















