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



int CombineGraphsRecoTrueF_Mai_2012_v2(string scale = "mmg_ik_MZ_Surface", string nomFitMethode = "RooBifurcatedGauss", string PileUpVersion = "2011", string LowMmumuLim = "40", string HightMmumuLim = "80", bool phiCracks = true, bool etaCracks = true, string SetOfCorrections = "ETHZCorrections", string variableX = "Photon_Et", int isMC = 1, int FitRange = 0, string MZbinning = "05GeV",string MuonCorrection = "Rochester", string SurfaceMethod = "Profile_Surface", string Category = "Vgamma8", string date = "16Jan")
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


	string compressEndCapsName;
        string compressR9Name;

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


	double xValueTabRecoInfEB[6] = {0.0};
	double xErrorLTabRecoInfEB[6] = {0.0};
	double xErrorRTabRecoInfEB[6] = {0.0};
	double MeanTabRecoInfEB[6] = {0.0};
	double MeanErrorTabRecoInfEB[6] = {0.0};

	double xValueTabRecoInfEE[6] = {0.0};
        double xErrorLTabRecoInfEE[6] = {0.0};
        double xErrorRTabRecoInfEE[6] = {0.0};
        double MeanTabRecoInfEE[6] = {0.0};
        double MeanErrorTabRecoInfEE[6] = {0.0};

	double xValueTabRecoSupEB[6] = {0.0};
        double xErrorLTabRecoSupEB[6] = {0.0};
        double xErrorRTabRecoSupEB[6] = {0.0};
        double MeanTabRecoSupEB[6] = {0.0};
        double MeanErrorTabRecoSupEB[6] = {0.0};

	double xValueTabRecoSupEE[6] = {0.0};
        double xErrorLTabRecoSupEE[6] = {0.0};
        double xErrorRTabRecoSupEE[6] = {0.0};
        double MeanTabRecoSupEE[6] = {0.0};
        double MeanErrorTabRecoSupEE[6] = {0.0};



	double xValueTabTrueInfEB[6] = {0.0};
        double xErrorLTabTrueInfEB[6] = {0.0};
        double xErrorRTabTrueInfEB[6] = {0.0};
        double MeanTabTrueInfEB[6] = {0.0};
        double MeanErrorTabTrueInfEB[6] = {0.0};

        double xValueTabTrueInfEE[6] = {0.0};
        double xErrorLTabTrueInfEE[6] = {0.0};
        double xErrorRTabTrueInfEE[6] = {0.0};
        double MeanTabTrueInfEE[6] = {0.0};
        double MeanErrorTabTrueInfEE[6] = {0.0};

        double xValueTabTrueSupEB[6] = {0.0};
        double xErrorLTabTrueSupEB[6] = {0.0};
        double xErrorRTabTrueSupEB[6] = {0.0};
        double MeanTabTrueSupEB[6] = {0.0};
        double MeanErrorTabTrueSupEB[6] = {0.0};

        double xValueTabTrueSupEE[6] = {0.0};
        double xErrorLTabTrueSupEE[6] = {0.0};
        double xErrorRTabTrueSupEE[6] = {0.0};
        double MeanTabTrueSupEE[6] = {0.0};
        double MeanErrorTabTrueSupEE[6] = {0.0};




	int rangeOpt = 0;
        if(FitRange == 0) rangeOpt = 1;
	int FitRangeTrue = 0;
	int FitRangeReco = 0;
	
	

	
	for(int i = 0; i < 4; i++)
	{
		if(i == 0 || i == 1) compressEndCapsName = "EB";
        	if(i == 2 || i == 3) compressEndCapsName = "EE";

        	if(i == 0 || i == 2) compressR9Name = "r9inf";
        	if(i == 1 || i == 3) compressR9Name = "r9sup";



                if(rangeOpt == 1){


                ////////// Results_Avril_2012_v2 //////////

                /// BFG ///
                // S True //
                if(nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MuonCorrection == "Rochester")
                {
                        if(i == 0) FitRangeTrue = 68; //low EB 
                        if(i == 1) FitRangeTrue = 81; //high EB
                        if(i == 2) FitRangeTrue = 72; //low EE
                        if(i == 3) FitRangeTrue = 84; //high EE
                }
                if(nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MuonCorrection == "NoMuonCorrection")
                {
                        if(i == 0) FitRangeTrue = 60; //low EB 
                        if(i == 1) FitRangeTrue = 82; //high EB
                        if(i == 2) FitRangeTrue = 73; //low EE
                        if(i == 3) FitRangeTrue = 84; //high EE
                }
		
		 // S Surface //
                if(nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MZbinning == "05GeV" && MuonCorrection == "Rochester" && SurfaceMethod == "Profile_Surface")
                {
                        if(i == 0) FitRangeReco = 80; //low EB 
                        if(i == 1) FitRangeReco = 92; //high EB
                        if(i == 2) FitRangeReco = 81; //low EE
                        if(i == 3) FitRangeReco = 81; //high EE
                }
                if(nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MZbinning == "05GeV" && MuonCorrection == "NoMuonCorrection" && SurfaceMethod == "Profile_Surface")
                {
                        if(i == 0) FitRangeReco = 80; //low EB 
                        if(i == 1) FitRangeReco = 88; //high EB
                        if(i == 2) FitRangeReco = 94; //low EE
                        if(i == 3) FitRangeReco = 83; //high EE
                }
                if(nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MZbinning == "1GeV" && MuonCorrection == "Rochester" && SurfaceMethod == "Profile_Surface")
                {
                        if(i == 0) FitRangeReco = 84; //low EB 
                        if(i == 1) FitRangeReco = 85; //high EB
                        if(i == 2) FitRangeReco = 82; //low EE
                        if(i == 3) FitRangeReco = 93; //high EE
                }
                if(nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MZbinning == "1GeV" && MuonCorrection == "NoMuonCorrection" && SurfaceMethod == "Profile_Surface")
                {
                        if(i == 0) FitRangeReco = 80; //low EB 
                        if(i == 1) FitRangeReco = 90; //high EB
                        if(i == 2) FitRangeReco = 93; //low EE
                        if(i == 3) FitRangeReco = 86; //high EE
                }
                if(nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MZbinning == "2GeV" && MuonCorrection == "Rochester" && SurfaceMethod == "Profile_Surface")
                {
                        if(i == 0) FitRangeReco = 82; //low EB 
                        if(i == 1) FitRangeReco = 85; //high EB
                        if(i == 2) FitRangeReco = 88; //low EE
                        if(i == 3) FitRangeReco = 88; //high EE
                }
                if(nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MZbinning == "2GeV" && MuonCorrection == "NoMuonCorrection" && SurfaceMethod == "Profile_Surface")
                {
                        if(i == 0) FitRangeReco = 79; //low EB 
                        if(i == 1) FitRangeReco = 74; //high EB
                        if(i == 2) FitRangeReco = 91; //low EE
                        if(i == 3) FitRangeReco = 90; //high EE
                }
		if(nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MZbinning == "5GeV" && MuonCorrection == "Rochester" && SurfaceMethod == "Profile_Surface")
                {
                        if(i == 0) FitRangeReco = 81; //low EB 
                        if(i == 1) FitRangeReco = 92; //high EB
                        if(i == 2) FitRangeReco = 87; //low EE
                        if(i == 3) FitRangeReco = 89; //high EE
                }
                if(nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MZbinning == "5GeV" && MuonCorrection == "NoMuonCorrection" && SurfaceMethod == "Profile_Surface")
                {
                        if(i == 0) FitRangeReco = 93; //low EB 
                        if(i == 1) FitRangeReco = 93; //high EB
                        if(i == 2) FitRangeReco = 92; //low EE
                        if(i == 3) FitRangeReco = 93; //high EE
                }
                if(nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MZbinning == "2GeV" && MuonCorrection == "NoMuonCorrection" && SurfaceMethod == "Fitted_Surface")
                {
                        if(i == 0) FitRangeReco = 88; //low EB 
                        if(i == 1) FitRangeReco = 76; //high EB
                        if(i == 2) FitRangeReco = 93; //low EE
                        if(i == 3) FitRangeReco = 86; //high EE
                }
                if(nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MZbinning == "2GeV" && MuonCorrection == "Rochester" && SurfaceMethod == "Fitted_Surface")
                {
                        if(i == 0) FitRangeReco = 85; //low EB 
                        if(i == 1) FitRangeReco = 83; //high EB
                        if(i == 2) FitRangeReco = 94; //low EE
                        if(i == 3) FitRangeReco = 95; //high EE
                }
                if(nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MZbinning == "5GeV" && MuonCorrection == "NoMuonCorrection" && SurfaceMethod == "Fitted_Surface")
                {
                        if(i == 0) FitRangeReco = 82; //low EB 
                        if(i == 1) FitRangeReco = 77; //high EB
                        if(i == 2) FitRangeReco = 86; //low EE
                        if(i == 3) FitRangeReco = 82; //high EE
                }
		if(nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC" && MZbinning == "5GeV" && MuonCorrection == "Rochester" && SurfaceMethod == "Fitted_Surface")
                {
                        if(i == 0) FitRangeReco = 93; //low EB 
                        if(i == 1) FitRangeReco = 72; //high EB
                        if(i == 2) FitRangeReco = 93; //low EE
                        if(i == 3) FitRangeReco = 82; //high EE
                }


                }
	
		for(int k = 0; k < 2; k++)
		{
	
			if(k == 0) FitRange = FitRangeTrue;
			if(k == 1) FitRange = FitRangeReco;

			cout<<endl<<"FitRange = "<<FitRange<<endl;
	
			if(k == 0)  scale ="Photon_E_o_MC_E";
			if(k == 1)  scale = "mmg_ik_MZ_Surface";

			if(scale == "mmg_ik_MZ_Surface") compressScaleName = "1overKreco";
			if(scale == "Photon_E_o_MC_E") compressScaleName = "ErecoOverEtrue";


			if(scale == "mmg_ik_MZ_Surface") fichierTXT = Form("Results_Mai_2012_v2/%s/%s/MmumuSelection%s_%s/%s/%s/%s/%s/MZbinning_%s/%s/%s/%s/%dPourcents/Fits/",variableX.c_str(),PileUpVersion.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),ChainCracks.c_str(),SetOfCorrections.c_str(),MuonCorrection.c_str(),isMCChain.c_str(),MZbinning.c_str(),SurfaceMethod.c_str(),Category.c_str(),date.c_str(),FitRange);
			
			if(scale == "Photon_E_o_MC_E") fichierTXT = Form("Results_Mai_2012_v2/%s/%s/MmumuSelection%s_%s/%s/%s/%s/MC/MZbinning_5GeV/Profile_Surface/%s/Fall11/%dPourcents/Fits/",variableX.c_str(),PileUpVersion.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),ChainCracks.c_str(),SetOfCorrections.c_str(),MuonCorrection.c_str(),Category.c_str(),FitRange);


                	fichierTXT += Form("%s%s%s%s/",compressScaleName.c_str(),compressFitMethodeName.c_str(),compressR9Name.c_str(),compressEndCapsName.c_str());


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
					
					if(k == 0)
                                        {
                                                cout<<"fichierxValueTXT = "<<fichierxValueTXT<<endl;
                                                monFlux1 >> xValueTabTrueInfEB[j];
                                                cout<<"xValueTabTrueInfEB[j] = "<<xValueTabTrueInfEB[j]<<endl;
                                                monFlux2 >> xErrorLTabTrueInfEB[j];
                                                monFlux3 >> xErrorRTabTrueInfEB[j];
                                                monFlux4 >> MeanTabTrueInfEB[j];
                                                monFlux5 >> MeanErrorTabTrueInfEB[j];
                                        }

					if(k == 1)
					{
						cout<<"fichierxValueTXT = "<<fichierxValueTXT<<endl;
						monFlux1 >> xValueTabRecoInfEB[j];
						cout<<"xValueTabRecoInfEB[j] = "<<xValueTabRecoInfEB[j]<<endl;
						monFlux2 >> xErrorLTabRecoInfEB[j];
						monFlux3 >> xErrorRTabRecoInfEB[j];
						monFlux4 >> MeanTabRecoInfEB[j];
						monFlux5 >> MeanErrorTabRecoInfEB[j];
					}
				}
				
				if(i == 1)
                                {

					if(k == 0)
                                        {
                                                monFlux1 >> xValueTabTrueSupEB[j];
                                                monFlux2 >> xErrorLTabTrueSupEB[j];
                                                monFlux3 >> xErrorRTabTrueSupEB[j];
                                                monFlux4 >> MeanTabTrueSupEB[j];
                                                monFlux5 >> MeanErrorTabTrueSupEB[j];
                                        }


					if(k == 1)
                                        {
	                                        monFlux1 >> xValueTabRecoSupEB[j];
	                                        monFlux2 >> xErrorLTabRecoSupEB[j];
	                                        monFlux3 >> xErrorRTabRecoSupEB[j];
	                                        monFlux4 >> MeanTabRecoSupEB[j];
	                                        monFlux5 >> MeanErrorTabRecoSupEB[j];
                                	}
				}			

				if(i == 2)
	                        {
					if(k == 0)
                                        {
                                                monFlux1 >> xValueTabTrueInfEE[j];
                                                monFlux2 >> xErrorLTabTrueInfEE[j];
                                                monFlux3 >> xErrorRTabTrueInfEE[j];
                                                monFlux4 >> MeanTabTrueInfEE[j];
                                                monFlux5 >> MeanErrorTabTrueInfEE[j];
                                        }
	
					if(k == 1)
                                        {
		                                monFlux1 >> xValueTabRecoInfEE[j];
		                                monFlux2 >> xErrorLTabRecoInfEE[j];
		                                monFlux3 >> xErrorRTabRecoInfEE[j];
		                                monFlux4 >> MeanTabRecoInfEE[j];
		                                monFlux5 >> MeanErrorTabRecoInfEE[j];
	                        	}
				}
	
				if(i == 3)
	                        {
		
					if(k == 0)
                                        {
                                                monFlux1 >> xValueTabTrueSupEE[j];
                                                monFlux2 >> xErrorLTabTrueSupEE[j];
                                                monFlux3 >> xErrorRTabTrueSupEE[j];
                                                monFlux4 >> MeanTabTrueSupEE[j];
                                                monFlux5 >> MeanErrorTabTrueSupEE[j];

                                        }


					if(k == 1)
                                        {
		                                monFlux1 >> xValueTabRecoSupEE[j];
		                                monFlux2 >> xErrorLTabRecoSupEE[j];
		                                monFlux3 >> xErrorRTabRecoSupEE[j];
		                                monFlux4 >> MeanTabRecoSupEE[j];
		                                monFlux5 >> MeanErrorTabRecoSupEE[j];
	                        
					}
				}
			}
		}	



		////////// Plots /////////

		
		TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);

        	TMultiGraph *mg = new TMultiGraph();

		TGraphAsymmErrors * ErecoOverEtrueRecoInfEB = new TGraphAsymmErrors(n,xValueTabRecoInfEB, MeanTabRecoInfEB, xErrorLTabRecoInfEB, xErrorRTabRecoInfEB, MeanErrorTabRecoInfEB, MeanErrorTabRecoInfEB);

        	TGraphAsymmErrors * ErecoOverEtrueRecoInfEE = new TGraphAsymmErrors(n,xValueTabRecoInfEE, MeanTabRecoInfEE, xErrorLTabRecoInfEE, xErrorRTabRecoInfEE, MeanErrorTabRecoInfEE, MeanErrorTabRecoInfEE);

        	TGraphAsymmErrors * ErecoOverEtrueRecoSupEB = new TGraphAsymmErrors(n,xValueTabRecoSupEB, MeanTabRecoSupEB, xErrorLTabRecoSupEB, xErrorRTabRecoSupEB, MeanErrorTabRecoSupEB, MeanErrorTabRecoSupEB);

        	TGraphAsymmErrors * ErecoOverEtrueRecoSupEE = new TGraphAsymmErrors(n,xValueTabRecoSupEE, MeanTabRecoSupEE, xErrorLTabRecoSupEE, xErrorRTabRecoSupEE, MeanErrorTabRecoSupEE, MeanErrorTabRecoSupEE);
	
		TGraphAsymmErrors * ErecoOverEtrueTrueInfEB = new TGraphAsymmErrors(n,xValueTabTrueInfEB, MeanTabTrueInfEB, xErrorLTabTrueInfEB, xErrorRTabTrueInfEB, MeanErrorTabTrueInfEB, MeanErrorTabTrueInfEB);

                TGraphAsymmErrors * ErecoOverEtrueTrueInfEE = new TGraphAsymmErrors(n,xValueTabTrueInfEE, MeanTabTrueInfEE, xErrorLTabTrueInfEE, xErrorRTabTrueInfEE, MeanErrorTabTrueInfEE, MeanErrorTabTrueInfEE);

                TGraphAsymmErrors * ErecoOverEtrueTrueSupEB = new TGraphAsymmErrors(n,xValueTabTrueSupEB, MeanTabTrueSupEB, xErrorLTabTrueSupEB, xErrorRTabTrueSupEB, MeanErrorTabTrueSupEB, MeanErrorTabTrueSupEB);

                TGraphAsymmErrors * ErecoOverEtrueTrueSupEE = new TGraphAsymmErrors(n,xValueTabTrueSupEE, MeanTabTrueSupEE, xErrorLTabTrueSupEE, xErrorRTabTrueSupEE, MeanErrorTabTrueSupEE, MeanErrorTabTrueSupEE);

		c1->ToggleEventStatus();
		
		if(i == 0)
		{
			mg->Add(ErecoOverEtrueRecoInfEB);
			mg->Add(ErecoOverEtrueTrueInfEB);
		
		}
		if(i == 1)
                {
                        mg->Add(ErecoOverEtrueRecoSupEB);
                        mg->Add(ErecoOverEtrueTrueSupEB);

                }
		if(i == 2)
                {
                        mg->Add(ErecoOverEtrueRecoInfEE);
                        mg->Add(ErecoOverEtrueTrueInfEE);

                }
		if(i == 3)
                {
                        mg->Add(ErecoOverEtrueRecoSupEE);
                        mg->Add(ErecoOverEtrueTrueSupEE);

                }

		mg->Draw("AP");

		mg->GetXaxis()->SetTitle("P_{T}^{#gamma}");
        	mg->GetXaxis()->SetLabelFont(42);
        	mg->GetXaxis()->SetTitleFont(42);
        	mg->GetXaxis()->SetLabelSize(0.03);
        	mg->GetYaxis()->SetTitle("s (%)");
		
		mg->GetYaxis()->SetLabelFont(42);
        	mg->GetYaxis()->SetTitleOffset(1.24);
        	mg->GetYaxis()->SetTitleFont(42);
        	mg->GetYaxis()->SetLabelSize(0.03);
		

		ErecoOverEtrueTrueInfEB->SetLineColor(597);
		ErecoOverEtrueRecoInfEB->SetLineColor(613);
		ErecoOverEtrueTrueInfEB->SetMarkerColor(597);
		ErecoOverEtrueRecoInfEB->SetMarkerColor(613);
		ErecoOverEtrueTrueInfEB->SetMarkerStyle(20);
		ErecoOverEtrueRecoInfEB->SetMarkerStyle(20);
		ErecoOverEtrueTrueInfEB->SetName("ErecoOverEtrueTrueInfEB");
		ErecoOverEtrueRecoInfEB->SetName("ErecoOverEtrueRecoInfEB");

		ErecoOverEtrueTrueSupEB->SetLineColor(597);
                ErecoOverEtrueRecoSupEB->SetLineColor(613);
                ErecoOverEtrueTrueSupEB->SetMarkerColor(597);
                ErecoOverEtrueRecoSupEB->SetMarkerColor(613);
                ErecoOverEtrueTrueSupEB->SetMarkerStyle(20);
                ErecoOverEtrueRecoSupEB->SetMarkerStyle(20);
                ErecoOverEtrueTrueSupEB->SetName("ErecoOverEtrueTrueSupEB");
                ErecoOverEtrueRecoSupEB->SetName("ErecoOverEtrueRecoSupEB");

		ErecoOverEtrueTrueInfEE->SetLineColor(597);
                ErecoOverEtrueRecoInfEE->SetLineColor(613);
                ErecoOverEtrueTrueInfEE->SetMarkerColor(597);
                ErecoOverEtrueRecoInfEE->SetMarkerColor(613);
                ErecoOverEtrueTrueInfEE->SetMarkerStyle(20);
                ErecoOverEtrueRecoInfEE->SetMarkerStyle(20);
                ErecoOverEtrueTrueInfEE->SetName("ErecoOverEtrueTrueInfEE");
                ErecoOverEtrueRecoInfEE->SetName("ErecoOverEtrueRecoInfEE");

		ErecoOverEtrueTrueSupEE->SetLineColor(597);
                ErecoOverEtrueRecoSupEE->SetLineColor(613);
                ErecoOverEtrueTrueSupEE->SetMarkerColor(597);
                ErecoOverEtrueRecoSupEE->SetMarkerColor(613);
                ErecoOverEtrueTrueSupEE->SetMarkerStyle(20);
                ErecoOverEtrueRecoSupEE->SetMarkerStyle(20);
                ErecoOverEtrueTrueSupEE->SetName("ErecoOverEtrueTrueSupEE");
                ErecoOverEtrueRecoSupEE->SetName("ErecoOverEtrueRecoSupEE");

		TLatex *text = new TLatex();

        	text->SetNDC();
        	text->SetTextAlign(11);
        	text->SetTextFont(42);
        	text->SetTextSizePixels(17);
        	text->SetTextSize(0.028);
        	text->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");

		
		if(i == 0) text->DrawLatex(0.16, 0.85, "Barrel, Low R9");
		if(i == 1) text->DrawLatex(0.16, 0.85, "Barrel, High R9");
		if(i == 2) text->DrawLatex(0.16, 0.85, "Endcaps, Low R9");
		if(i == 3) text->DrawLatex(0.16, 0.85, "Endcaps, High R9");	

		TLegend * leg = new TLegend(0.6,0.7,0.9,0.9,"","brNDC");
		leg->SetTextSize(0.030);
        	leg->SetFillColor(kWhite);
        	leg->SetLineColor(kWhite);
        	leg->SetShadowColor(kWhite);
		
		if(i == 0)
		{
			leg->AddEntry(ErecoOverEtrueTrueInfEB->GetName(),"s_{TRUE}","lep");
			leg->AddEntry(ErecoOverEtrueRecoInfEB->GetName(),"s_{RECO Surface}","lep");
			
		}
		if(i == 1)
                {
                        leg->AddEntry(ErecoOverEtrueTrueSupEB->GetName(),"s_{TRUE}","lep");
                        leg->AddEntry(ErecoOverEtrueRecoSupEB->GetName(),"s_{RECO Surface}","lep");

                }
		if(i == 2)
                {
                        leg->AddEntry(ErecoOverEtrueTrueInfEE->GetName(),"s_{TRUE}","lep");
                        leg->AddEntry(ErecoOverEtrueRecoInfEE->GetName(),"s_{RECO Surface}","lep");

                }
		if(i == 3)
                {
                        leg->AddEntry(ErecoOverEtrueTrueSupEE->GetName(),"s_{TRUE}","lep");
                        leg->AddEntry(ErecoOverEtrueRecoSupEE->GetName(),"s_{RECO Surface}","lep");

                }
		
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

		//DossierEnregistrement = Form("Results_Avril_2012_v2/%s/%s/MmumuSelection%s_%s/%s/%s/%s/%s/MZbinning_%s/%s/CombinedGraphs/",variableX.c_str(),PileUpVersion.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),ChainCracks.c_str(),SetOfCorrections.c_str(),MuonCorrection.c_str(),isMCChain.c_str(),MZbinning.c_str(),SurfaceMethod.c_str());	
		DossierEnregistrement = Form("Results_Mai_2012_v2/%s/%s/MmumuSelection%s_%s/%s/%s/%s/%s/MZbinning_%s/%s/%s/%s/CombinedGraphs/",variableX.c_str(),PileUpVersion.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),DossierCracks.c_str(),SetOfCorrections.c_str(),MuonCorrection.c_str(),isMCChain.c_str(),MZbinning.c_str(),SurfaceMethod.c_str(),Category.c_str(),date.c_str());	
	
		FichierEnregistrement = Form("CombineSRecoStrueVsPt_%s_%s_%s",compressFitMethodeName.c_str(),compressR9Name.c_str(),compressEndCapsName.c_str());

		enregistrementPlots(DossierEnregistrement, FichierEnregistrement, 2, 10000, c1);

		c1->Clear();
		/*delete mg;
		mg = 0;
		delete ErecoOverEtrueRecoInfEB;
		ErecoOverEtrueRecoInfEB = 0;
		delete ErecoOverEtrueRecoSupEB;
                ErecoOverEtrueRecoSupEB = 0;
		delete ErecoOverEtrueRecoInfEE;
                ErecoOverEtrueRecoInfEE = 0;
		delete ErecoOverEtrueRecoSupEE;
                ErecoOverEtrueRecoSupEE = 0;

		delete ErecoOverEtrueTrueInfEB;
                ErecoOverEtrueTrueInfEB = 0;
                delete ErecoOverEtrueTrueSupEB;
                ErecoOverEtrueTrueSupEB = 0;
                delete ErecoOverEtrueTrueInfEE;
                ErecoOverEtrueTrueInfEE = 0;
                delete ErecoOverEtrueTrueSupEE;
                ErecoOverEtrueTrueSupEE = 0;

		delete text;
		text = 0;
		delete leg;
		leg = 0;

		delete c1;
		c1 = 0;
*/
	}







	return 0;
}
