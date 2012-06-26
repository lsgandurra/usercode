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



int Tableau_Strue_Sreco_Avril_2012(string scale = "mmg_ik_MZ_Surface", string nomFitMethode = "RooBifurcatedGauss", string PileUpVersion = "2011", string LowMmumuLim = "40", string HightMmumuLim = "80", bool phiCracks = true, bool etaCracks = true, string SetOfCorrections = "ETHZCorrections", string variableX = "Photon_Et", int isMC = 1, int FitRange = 0, string MZbinning = "2GeV",string MuonCorrection = "NoMuonCorrection", string SurfaceMethod = "Fitted_Surface", string Category = "Vgamma8", string date = "16Jan") //"Profile_Surface")
{
	//Fitted_Surface
	//Profile_Surface
	//ElectronTunedCorrections
	gROOT->Reset();
	setTDRStyle();
	TGaxis::SetMaxDigits(3);
	gStyle->SetOptStat("emr");

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
        if(scale == "ComparaisonFitsErecoOverEtrueDiffMmumuJan") compressScaleName = "ErecoOverEtrue";
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


		////////// Results_Avril_2012_v2 //////////
			
		/// BFG ///
		if(scale == "mmg_ik_MZ_Surface" && nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MZbinning == "05GeV" && MuonCorrection == "Rochester" && SurfaceMethod == "Profile_Surface")
                {
                        if(i == 0) FitRange = 80; //low EB 
                        if(i == 1) FitRange = 92; //high EB
                        if(i == 2) FitRange = 81; //low EE
                        if(i == 3) FitRange = 81; //high EE
                }	
		if(scale == "mmg_ik_MZ_Surface" && nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MZbinning == "05GeV" && MuonCorrection == "NoMuonCorrection" && SurfaceMethod == "Profile_Surface")
                {
                        if(i == 0) FitRange = 80; //low EB 
                        if(i == 1) FitRange = 88; //high EB
                        if(i == 2) FitRange = 94; //low EE
                        if(i == 3) FitRange = 83; //high EE
                }
		if(scale == "mmg_ik_MZ_Surface" && nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MZbinning == "1GeV" && MuonCorrection == "Rochester" && SurfaceMethod == "Profile_Surface")
                {
                        if(i == 0) FitRange = 84; //low EB 
                        if(i == 1) FitRange = 85; //high EB
                        if(i == 2) FitRange = 82; //low EE
                        if(i == 3) FitRange = 93; //high EE
                }
                if(scale == "mmg_ik_MZ_Surface" && nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MZbinning == "1GeV" && MuonCorrection == "NoMuonCorrection" && SurfaceMethod == "Profile_Surface")
                {
                        if(i == 0) FitRange = 80; //low EB 
                        if(i == 1) FitRange = 90; //high EB
                        if(i == 2) FitRange = 93; //low EE
                        if(i == 3) FitRange = 86; //high EE
                }
		if(scale == "mmg_ik_MZ_Surface" && nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MZbinning == "2GeV" && MuonCorrection == "Rochester" && SurfaceMethod == "Profile_Surface")
                {
                        if(i == 0) FitRange = 82; //low EB 
                        if(i == 1) FitRange = 85; //high EB
                        if(i == 2) FitRange = 88; //low EE
                        if(i == 3) FitRange = 88; //high EE
                }
                if(scale == "mmg_ik_MZ_Surface" && nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MZbinning == "2GeV" && MuonCorrection == "NoMuonCorrection" && SurfaceMethod == "Profile_Surface")
                {
                        if(i == 0) FitRange = 79; //low EB 
                        if(i == 1) FitRange = 74; //high EB
                        if(i == 2) FitRange = 91; //low EE
                        if(i == 3) FitRange = 90; //high EE
                }
		if(scale == "mmg_ik_MZ_Surface" && nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MZbinning == "5GeV" && MuonCorrection == "Rochester" && SurfaceMethod == "Profile_Surface")
                {
                        if(i == 0) FitRange = 81; //low EB 
                        if(i == 1) FitRange = 92; //high EB
                        if(i == 2) FitRange = 87; //low EE
                        if(i == 3) FitRange = 89; //high EE
                }
                if(scale == "mmg_ik_MZ_Surface" && nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MZbinning == "5GeV" && MuonCorrection == "NoMuonCorrection" && SurfaceMethod == "Profile_Surface")
                {
                        if(i == 0) FitRange = 93; //low EB 
                        if(i == 1) FitRange = 93; //high EB
                        if(i == 2) FitRange = 92; //low EE
                        if(i == 3) FitRange = 93; //high EE
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
		fichierTXT = Form("Results_Mai_2012_v2/%s/%s/MmumuSelection%s_%s/%s/%s/%s/%s/MZbinning_%s/%s/%s/%s/%dPourcents/Fits/",variableX.c_str(),PileUpVersion.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),DossierCracks.c_str(),SetOfCorrections.c_str(),MuonCorrection.c_str(),isMCChain.c_str(),MZbinning.c_str(),SurfaceMethod.c_str(),Category.c_str(),date.c_str(),FitRange);


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


	double xValueTab_strue_R9infEB[6] = {0.0};
        double xErrorLTab_strue_R9infEB[6] = {0.0};
        double xErrorRTab_strue_R9infEB[6] = {0.0};
        double MeanTab_strue_R9infEB[6] = {0.0};
        double MeanErrorTab_strue_R9infEB[6] = {0.0};

        double xValueTab_strue_R9supEB[6] = {0.0};
        double xErrorLTab_strue_R9supEB[6] = {0.0};
        double xErrorRTab_strue_R9supEB[6] = {0.0};
        double MeanTab_strue_R9supEB[6] = {0.0};
        double MeanErrorTab_strue_R9supEB[6] = {0.0};

        double xValueTab_strue_R9infEE[6] = {0.0};
        double xErrorLTab_strue_R9infEE[6] = {0.0};
        double xErrorRTab_strue_R9infEE[6] = {0.0};
        double MeanTab_strue_R9infEE[6] = {0.0};
        double MeanErrorTab_strue_R9infEE[6] = {0.0};

        double xValueTab_strue_R9supEE[6] = {0.0};
        double xErrorLTab_strue_R9supEE[6] = {0.0};
        double xErrorRTab_strue_R9supEE[6] = {0.0};
        double MeanTab_strue_R9supEE[6] = {0.0};
        double MeanErrorTab_strue_R9supEE[6] = {0.0};


	string compressScaleName_2 = "ErecoOverEtrue";
	string isMCChain_2 = "MC";
	string SurfaceMethod_2 = "Profile_Surface";
	string MZbinning_2 = "5GeV";
	string date_2 = "Fall11";

	for(int i = 0; i < 4; i++)
	{

		if(rangeOpt == 1){


		////////// Results_Avril_2012_v2 //////////
			
		/// BFG ///
                if(nomFitMethode == "RooBifurcatedGauss" && MuonCorrection == "Rochester")
		{
                        if(i == 0) FitRange = 68; //low EB 
                        if(i == 1) FitRange = 81; //high EB
                        if(i == 2) FitRange = 72; //low EE
                        if(i == 3) FitRange = 84; //high EE
                }    
                if(nomFitMethode == "RooBifurcatedGauss" && MuonCorrection == "NoMuonCorrection")
		{
                        if(i == 0) FitRange = 60; //low EB 
                        if(i == 1) FitRange = 82; //high EB
                        if(i == 2) FitRange = 73; //low EE
                        if(i == 3) FitRange = 84; //high EE
                }

                }


		//fichierTXT = Form("Results_Avril_2012_v1/%s/%s/MmumuSelection%s_%s/%s/%s/%s/%s/MZbinning_2_%s/%dPourcents/Fits/",variableX.c_str(),PileUpVersion.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),ChainCracks.c_str(),SetOfCorrections.c_str(),MuonCorrection.c_str(),isMCChain.c_str(),MZbinning_2.c_str(),FitRange);
		//fichierTXT = Form("Results_Avril_2012_v2/%s/%s/MmumuSelection%s_%s/%s/%s/%s/%s/MZbinning_%s/%s/%dPourcents/Fits/",variableX.c_str(),PileUpVersion.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),ChainCracks.c_str(),SetOfCorrections.c_str(),MuonCorrection.c_str(),isMCChain_2.c_str(),MZbinning_2.c_str(),SurfaceMethod_2.c_str(),FitRange);
		fichierTXT = Form("Results_Mai_2012_v2/%s/%s/MmumuSelection%s_%s/%s/%s/%s/%s/MZbinning_%s/%s/%s/%s/%dPourcents/Fits/",variableX.c_str(),PileUpVersion.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),ChainCracks.c_str(),SetOfCorrections.c_str(),MuonCorrection.c_str(),isMCChain_2.c_str(),MZbinning_2.c_str(),SurfaceMethod_2.c_str(),Category.c_str(),date_2.c_str(),FitRange);

                fichierTXT += Form("%s%s",compressScaleName_2.c_str(),compressFitMethodeName.c_str());



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
				monFlux1 >> xValueTab_strue_R9infEB[j];
				cout<<"xValueTab_strue_R9infEB[j] = "<<xValueTab_strue_R9infEB[j]<<endl;
				monFlux2 >> xErrorLTab_strue_R9infEB[j];
				monFlux3 >> xErrorRTab_strue_R9infEB[j];
				monFlux4 >> MeanTab_strue_R9infEB[j];
				monFlux5 >> MeanErrorTab_strue_R9infEB[j];
			}
			
			if(i == 1)
                        {
                                monFlux1 >> xValueTab_strue_R9supEB[j];
                                monFlux2 >> xErrorLTab_strue_R9supEB[j];
                                monFlux3 >> xErrorRTab_strue_R9supEB[j];
                                monFlux4 >> MeanTab_strue_R9supEB[j];
                                monFlux5 >> MeanErrorTab_strue_R9supEB[j];
                        }

			if(i == 2)
                        {
                                monFlux1 >> xValueTab_strue_R9infEE[j];
                                monFlux2 >> xErrorLTab_strue_R9infEE[j];
                                monFlux3 >> xErrorRTab_strue_R9infEE[j];
                                monFlux4 >> MeanTab_strue_R9infEE[j];
                                monFlux5 >> MeanErrorTab_strue_R9infEE[j];
                        }

			if(i == 3)
                        {
                                monFlux1 >> xValueTab_strue_R9supEE[j];
                                monFlux2 >> xErrorLTab_strue_R9supEE[j];
                                monFlux3 >> xErrorRTab_strue_R9supEE[j];
                                monFlux4 >> MeanTab_strue_R9supEE[j];
                                monFlux5 >> MeanErrorTab_strue_R9supEE[j];
                        }
	
		}	

	}


	//DossierEnregistrement = Form("Results_Avril_2012_v2/%s/%s/MmumuSelection%s_%s/%s/%s/%s/%s/MZbinning_%s/%s/CombinedGraphs/",variableX.c_str(),PileUpVersion.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),ChainCracks.c_str(),SetOfCorrections.c_str(),MuonCorrection.c_str(),isMCChain.c_str(),MZbinning.c_str(),SurfaceMethod.c_str());
	DossierEnregistrement = Form("Results_Mai_2012_v2/%s/%s/MmumuSelection%s_%s/%s/%s/%s/%s/MZbinning_%s/%s/%s/%s/CombinedGraphs/",variableX.c_str(),PileUpVersion.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),DossierCracks.c_str(),SetOfCorrections.c_str(),MuonCorrection.c_str(),isMCChain.c_str(),MZbinning.c_str(),SurfaceMethod.c_str(),Category.c_str(),date.c_str());

	FichierEnregistrement = Form("Tableau_Strue_Sreco_%s",compressFitMethodeName.c_str());

	//FILE *f1 = fopen(Form("%s%s.txt",DossierEnregistrement.c_str(),FichierEnregistrement.c_str()), "w" );
        //delete f1; 

	cout<<endl<<"DossierEnregistrement = "<<DossierEnregistrement<<endl;
	cout<<endl<<"FichierEnregistrement = "<<FichierEnregistrement<<endl;

        //ofstream monFlux20(Form("%s%s.txt",DossierEnregistrement.c_str(),FichierEnregistrement.c_str()), ios::app);
	ofstream monFlux20(Form("%s%s.txt",DossierEnregistrement.c_str(),FichierEnregistrement.c_str()));


	for(int i = 0; i < n; i++)
	{	

		if(i == 0) monFlux20<<"\\cellcolor{Dandelion} 10-12 GeV";
		if(i == 1) monFlux20<<"\\cellcolor{Dandelion} 12-15 GeV";
		if(i == 2) monFlux20<<"\\cellcolor{Dandelion} 15-20 GeV";
		if(i == 3) monFlux20<<"\\cellcolor{Dandelion} 20-25 GeV";
		if(i == 4) monFlux20<<"\\cellcolor{Dandelion} 25-30 GeV";
		if(i == 5) monFlux20<<"\\cellcolor{Dandelion} 30-100 GeV";
		monFlux20<<" & \\cellcolor{Periwinkle} low r9 EB"<<" & "<<MeanTab_strue_R9infEB[i]<<" & "<<MeanErrorTab_strue_R9infEB[i]<<" & "<<MeanTabR9infEB[i]<<" & "<<MeanErrorTabR9infEB[i]<<" & "<<MeanTab_strue_R9infEB[i] - MeanTabR9infEB[i]<<" \\\\";
		monFlux20<<endl;
		monFlux20<<" & low r9 EE"<<" & "<<MeanTab_strue_R9infEE[i]<<" & "<<MeanErrorTab_strue_R9infEE[i]<<" & "<<MeanTabR9infEE[i]<<" & "<<MeanErrorTabR9infEE[i]<<" & "<<MeanTab_strue_R9infEE[i] - MeanTabR9infEE[i]<<" \\\\";
		monFlux20<<endl;
		monFlux20<<"\\cellcolor{white} & \\cellcolor{Periwinkle} high r9 EB"<<" & "<<MeanTab_strue_R9supEB[i]<<" & "<<MeanErrorTab_strue_R9supEB[i]<<" & "<<MeanTabR9supEB[i]<<" & "<<MeanErrorTabR9supEB[i]<<" & "<<MeanTab_strue_R9supEB[i] - MeanTabR9supEB[i]<<" \\\\";
                monFlux20<<endl;
                monFlux20<<" & high r9 EE"<<" & "<<MeanTab_strue_R9supEE[i]<<" & "<<MeanErrorTab_strue_R9supEE[i]<<" & "<<MeanTabR9supEE[i]<<" & "<<MeanErrorTabR9supEE[i]<<" & "<<MeanTab_strue_R9supEE[i] - MeanTabR9supEE[i]<<" \\\\";
		monFlux20<<endl;	
	}

	monFlux20.close();

	TCanvas *c2 = new TCanvas("c2", "c2",0,0,600,600);

	TH1D * strueMenusSreco = new TH1D("strueMenusSreco","strueMenusSreco",26,-13,13);

	for(int j = 0; j < n; j++)
	{
		strueMenusSreco->Fill(MeanTabR9infEB[j] - MeanTab_strue_R9infEB[j]);
		strueMenusSreco->Fill(MeanTabR9infEE[j] - MeanTab_strue_R9infEE[j]);
		strueMenusSreco->Fill(MeanTabR9supEB[j] - MeanTab_strue_R9supEB[j]);
		strueMenusSreco->Fill(MeanTabR9supEE[j] - MeanTab_strue_R9supEE[j]);

	}
		



	strueMenusSreco->Draw("");

	strueMenusSreco->GetXaxis()->SetTitle("S_{RECO Surface} - S_{TRUE} (%)");
	strueMenusSreco->GetYaxis()->SetTitle("entries / 1 (%)");

	cout<<endl<<"strueMenusSreco->GetEntries() = "<<strueMenusSreco->GetEntries()<<endl;
	

	double rangeSup = strueMenusSreco->GetMean() + 4 * strueMenusSreco->GetRMS();
	double rangeInf = strueMenusSreco->GetMean() - 4 * strueMenusSreco->GetRMS();

	cout<<endl<<"rangeSup = "<<rangeSup<<", rangeInf = "<<rangeInf<<endl;
		

	TF1 * f1 = new TF1("f1","gaus",-13,13);
	
	//TF1 * f1 = new TF1("f1","gaus",-7,7);

	//f1->SetParameters(10,strueMenusSreco->GetMean(), 1);
	f1->SetParameters(5,strueMenusSreco->GetMean(), strueMenusSreco->GetRMS());	
	
	f1->SetLineColor(2);	
	f1->SetFillColor(19);
	f1->SetFillStyle(0);
	f1->SetLineWidth(2);
	//strueMenusSreco->Fit("f1","","",rangeInf,rangeSup);
	strueMenusSreco->Fit("f1");	


	//gStyle->SetOptStat("emr");

	FichierEnregistrement = Form("Plot1D_SrecoSurface_minus_Strue_%s",compressFitMethodeName.c_str());

	enregistrementPlots(DossierEnregistrement,FichierEnregistrement,2,10000, c2);

	//c2->Print("DiffTest2.png");

/*
	delete c2;
	c2 = 0;
	delete strueMenusSreco;
	strueMenusSreco = 0;
	delete f1;
	f1 = 0;	
*/

	return 0;
}
