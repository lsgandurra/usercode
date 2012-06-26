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



int CombinePValuesRangeOptF_Mars_2012(int EndCaps = 0, int r9sup = 1, string scale = "ComparaisonFitsErecoOverEtrueDiffMmumuJan", string nomFitMethode = "RooBifurcatedGauss", string PileUpVersion = "2011", string LowMmumuLim = "40", string HightMmumuLim = "80", bool phiCracks = true, bool etaCracks = true, int FitRange = 83, string SetOfCorrections = "ETHZCorrections", string variableX = "Photon_Et", int isMC = 1, string MZbinning = "05GeV")
{
	//ElectronTunedCorrections
	gROOT->Reset();
	setTDRStyle();
	TGaxis::SetMaxDigits(3);

	double yminPValues, ymaxPValues;
	double xminPValues, xmaxPValues;
	
	yminPValues = 0.0;
	ymaxPValues = 1.0;
	xminPValues = 0.0;
        xmaxPValues = 100.0;

        if(variableX == "Photon_SC_Eta") xminPValues = -3.0;
        xmaxPValues = 100.0;
        if(variableX == "Photon_SC_rawEt") xmaxPValues = 250.0;
        if(variableX == "Photon_Et") xmaxPValues = 100.0;
        if(variableX == "Photon_E") xmaxPValues = 350.0;
        if(variableX == "Photon_SC_Eta") xmaxPValues = 3.0;
        if(variableX == "Photon_SC_brem") xmaxPValues = 15.0;


	int n = 6;

	string compressEndCapsName = ""; 
        if(EndCaps == 0) compressEndCapsName = "EB";
        if(EndCaps == 1) compressEndCapsName = "EE";

        string compressR9Name = ""; 
        if(r9sup == 0) compressR9Name = "r9inf";
        if(r9sup == 1) compressR9Name = "r9sup";
	
	string compressFitMethodeName;
	if(nomFitMethode == "RooCrystalBall") compressFitMethodeName = "CB";
	if(nomFitMethode == "RooBifurcatedGauss") compressFitMethodeName = "BFG";
	if(nomFitMethode == "RooLandauConvGaussian") compressFitMethodeName = "LandauConvGaus";
        if(nomFitMethode == "RooLandau2") compressFitMethodeName = "Landau";

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
	string fichierPValueTabTXT = "";
	string fichierPValueErrorTabTXT = "";


	double xValueTabOptRangeM2[6] = {0.0};
	double xErrorLTabOptRangeM2[6] = {0.0};
	double xErrorRTabOptRangeM2[6] = {0.0};
	double PValueTabOptRangeM2[6] = {0.0};
	double PValueErrorTabOptRangeM2[6] = {0.0};

	double xValueTabOptRangeM1[6] = {0.0};
        double xErrorLTabOptRangeM1[6] = {0.0};
        double xErrorRTabOptRangeM1[6] = {0.0};
        double PValueTabOptRangeM1[6] = {0.0};
        double PValueErrorTabOptRangeM1[6] = {0.0};

	double xValueTabOptRange[6] = {0.0};
        double xErrorLTabOptRange[6] = {0.0};
        double xErrorRTabOptRange[6] = {0.0};
        double PValueTabOptRange[6] = {0.0};
        double PValueErrorTabOptRange[6] = {0.0};

	double xValueTabOptRangeP1[6] = {0.0};
        double xErrorLTabOptRangeP1[6] = {0.0};
        double xErrorRTabOptRangeP1[6] = {0.0};
        double PValueTabOptRangeP1[6] = {0.0};
        double PValueErrorTabOptRangeP1[6] = {0.0};

	double xValueTabOptRangeP2[6] = {0.0};
        double xErrorLTabOptRangeP2[6] = {0.0};
        double xErrorRTabOptRangeP2[6] = {0.0};
        double PValueTabOptRangeP2[6] = {0.0};
        double PValueErrorTabOptRangeP2[6] = {0.0};


	for(int i = FitRange - 2; i <= FitRange + 2; i++)
	{

		fichierTXT = Form("Results_Mars_2012_v2/%s/%s/MmumuSelection%s_%s/%s/%s/%s/MZbinning_%s/%dPourcents/Fits/",variableX.c_str(),PileUpVersion.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),ChainCracks.c_str(),SetOfCorrections.c_str(),isMCChain.c_str(),MZbinning.c_str(),i);


                fichierTXT += Form("%s%s%s%s/",compressScaleName.c_str(),compressFitMethodeName.c_str(),compressR9Name.c_str(),compressEndCapsName.c_str());


		fichierxValueTXT = fichierTXT + "xValue.txt";
                fichierxErrorLTXT = fichierTXT + "xErrorL.txt";
                fichierxErrorRTXT = fichierTXT + "xErrorR.txt";
                fichierPValueTabTXT = fichierTXT + "PValueTab.txt";


        	cout<<"fichierTXT = "<<fichierTXT<<endl;
        	ifstream monFlux1(fichierxValueTXT.c_str());	
		ifstream monFlux2(fichierxErrorLTXT.c_str());
		ifstream monFlux3(fichierxErrorRTXT.c_str());
		ifstream monFlux4(fichierPValueTabTXT.c_str());

		for(int j = 0; j < n; j++)
		{ 

			if(i == FitRange - 2)
			{
				cout<<"fichierxValueTXT = "<<fichierxValueTXT<<endl;
				monFlux1 >> xValueTabOptRangeM2[j];
				cout<<"xValueTabOptRangeM2[j] = "<<xValueTabOptRangeM2[j]<<endl;
				monFlux2 >> xErrorLTabOptRangeM2[j];
				monFlux3 >> xErrorRTabOptRangeM2[j];
				monFlux4 >> PValueTabOptRangeM2[j];
			}
			
			if(i == FitRange - 1)
                        {
                                monFlux1 >> xValueTabOptRangeM1[j];
                                monFlux2 >> xErrorLTabOptRangeM1[j];
                                monFlux3 >> xErrorRTabOptRangeM1[j];
                                monFlux4 >> PValueTabOptRangeM1[j];
                        }

			if(i == FitRange)
                        {
                                monFlux1 >> xValueTabOptRange[j];
                                monFlux2 >> xErrorLTabOptRange[j];
                                monFlux3 >> xErrorRTabOptRange[j];
                                monFlux4 >> PValueTabOptRange[j];
                        }

			if(i == FitRange + 1)
                        {
                                monFlux1 >> xValueTabOptRangeP1[j];
                                monFlux2 >> xErrorLTabOptRangeP1[j];
                                monFlux3 >> xErrorRTabOptRangeP1[j];
                                monFlux4 >> PValueTabOptRangeP1[j];
                        }

			if(i == FitRange + 2)
                        {
                                monFlux1 >> xValueTabOptRangeP2[j];
				monFlux2 >> xErrorLTabOptRangeP2[j];
                                monFlux3 >> xErrorRTabOptRangeP2[j];
                                monFlux4 >> PValueTabOptRangeP2[j];
                        }
	
		}	

	}


	TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);

	TMultiGraph *mg = new TMultiGraph();


	TGraphAsymmErrors * PValuesOptRangeM2 = new TGraphAsymmErrors(n,xValueTabOptRangeM2, PValueTabOptRangeM2, xErrorLTabOptRangeM2, xErrorRTabOptRangeM2, PValueErrorTabOptRangeM2, PValueErrorTabOptRangeM2);

	TGraphAsymmErrors * PValuesOptRangeM1 = new TGraphAsymmErrors(n,xValueTabOptRangeM1, PValueTabOptRangeM1, xErrorLTabOptRangeM1, xErrorRTabOptRangeM1, PValueErrorTabOptRangeM1, PValueErrorTabOptRangeM1);

	TGraphAsymmErrors * PValuesOptRange = new TGraphAsymmErrors(n,xValueTabOptRange, PValueTabOptRange, xErrorLTabOptRange, xErrorRTabOptRange, PValueErrorTabOptRange, PValueErrorTabOptRange);

	TGraphAsymmErrors * PValuesOptRangeP1 = new TGraphAsymmErrors(n,xValueTabOptRangeP1, PValueTabOptRangeP1, xErrorLTabOptRangeP1, xErrorRTabOptRangeP1, PValueErrorTabOptRangeP1, PValueErrorTabOptRangeP1);

	TGraphAsymmErrors * PValuesOptRangeP2 = new TGraphAsymmErrors(n,xValueTabOptRangeP2, PValueTabOptRangeP2, xErrorLTabOptRangeP2, xErrorRTabOptRangeP2, PValueErrorTabOptRangeP2, PValueErrorTabOptRangeP2);


        c1->ToggleEventStatus();
        //PValuesOptRangeM2->Draw("AP");
	//PValuesOptRangeM1->Draw("SAMES");
	//PValuesOptRange->Draw("SAMES");
	//PValuesOptRangeP1->Draw("SAMES");

	mg->Add(PValuesOptRangeM2);
	mg->Add(PValuesOptRangeM1);
	mg->Add(PValuesOptRange);
	mg->Add(PValuesOptRangeP1);
	mg->Add(PValuesOptRangeP2);	
	mg->Draw("AP");

	//mg->GetXaxis()->SetTitle("Pt_{#gamma}");
	
	if(variableX == "Photon_SC_rawEt") mg->GetXaxis()->SetTitle("P_{T RAW}^{#gamma}");
        if(variableX == "Photon_Et") mg->GetXaxis()->SetTitle("P_{T}^{#gamma}");
        if(variableX == "Photon_E") mg->GetXaxis()->SetTitle("E^{#gamma}");
        if(variableX == "Photon_SC_Eta") mg->GetXaxis()->SetTitle("#eta");
        if(variableX == "Photon_SC_brem") mg->GetXaxis()->SetTitle("#sigma_{#phi}/#sigma_{#eta}");

	mg->GetXaxis()->SetLabelFont(42);
        mg->GetXaxis()->SetTitleFont(42);
        mg->GetXaxis()->SetLabelSize(0.03);
        //mg->GetYaxis()->SetTitle("E_{RECO}/E_{TRUE} - 1 (%)");
		
	mg->GetYaxis()->SetTitle("p-value");

	mg->GetYaxis()->SetLabelFont(42);
        mg->GetYaxis()->SetTitleOffset(1.24);
        mg->GetYaxis()->SetTitleFont(42);
        mg->GetYaxis()->SetLabelSize(0.03);
	//mg->SetTitle("");
        //mg->SetMarkerColor(4);
        //mg->SetMarkerStyle(21);
        //mg->SetMarkerSize(0.6);

	PValuesOptRangeM2->SetLineColor(629);
        PValuesOptRangeM1->SetLineColor(613);//613
        PValuesOptRange->SetLineColor(1);
        PValuesOptRangeP1->SetLineColor(429);//429
        PValuesOptRangeP2->SetLineColor(597);//597
     
        //413
        //397
        PValuesOptRangeM2->SetMarkerColor(629);
        PValuesOptRangeM1->SetMarkerColor(613);
        PValuesOptRange->SetMarkerColor(1);
        PValuesOptRangeP1->SetMarkerColor(429);
        PValuesOptRangeP2->SetMarkerColor(597);

	PValuesOptRangeM2->SetMarkerStyle(20);
        PValuesOptRangeM1->SetMarkerStyle(21);
        PValuesOptRange->SetMarkerStyle(34);
        PValuesOptRangeP1->SetMarkerStyle(25);
        PValuesOptRangeP2->SetMarkerStyle(24);	
 
        //413
        //397	

	PValuesOptRangeM2->SetName("PValuesOptRangeM2");
	PValuesOptRangeM1->SetName("PValuesOptRangeM1");
	PValuesOptRange->SetName("PValuesOptRange");
	PValuesOptRangeP1->SetName("PValuesOptRangeP1");
	PValuesOptRangeP2->SetName("PValuesOptRangeP2");

	TLatex *text = new TLatex();

        text = new TLatex();
        text->SetNDC();
        text->SetTextAlign(11);
        text->SetTextFont(42);
        text->SetTextSizePixels(17);
	text->SetTextSize(0.028);
	text->DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
	//if(nomFitMethode == "RooCrystalBall") text->DrawLatex(0.16, 0.85, "Crystal Ball");
        //if(nomFitMethode == "RooBifurcatedGauss") text->DrawLatex(0.16, 0.85, "Bifurcated Gaussian");
	if(EndCaps == 0  &&  r9sup == 1 && nomFitMethode == "RooCrystalBall") text->DrawLatex(0.16, 0.85, "Barrel, High R9, Crystal Ball");
        if(EndCaps == 0  &&  r9sup == 0 && nomFitMethode == "RooCrystalBall") text->DrawLatex(0.16, 0.85, "Barrel, Low R9, Crystal Ball");
        if(EndCaps == 0  &&  r9sup == 1 && nomFitMethode == "RooBifurcatedGauss") text->DrawLatex(0.16, 0.85, "Barrel, High R9, Bifurcated Gaussian");
        if(EndCaps == 0  &&  r9sup == 0 && nomFitMethode == "RooBifurcatedGauss") text->DrawLatex(0.16, 0.85, "Barrel, Low R9, Bifurcated Gaussian");
        if(EndCaps == 1  &&  r9sup == 1 && nomFitMethode == "RooCrystalBall") text->DrawLatex(0.16, 0.85, "Endcaps, High R9, Crystal Ball");
        if(EndCaps == 1  &&  r9sup == 0 && nomFitMethode == "RooCrystalBall") text->DrawLatex(0.16, 0.85, "Endcaps, Low R9, Crystal Ball");
        if(EndCaps == 1  &&  r9sup == 1 && nomFitMethode == "RooBifurcatedGauss") text->DrawLatex(0.16, 0.85, "Endcaps, High R9, Bifurcated Gaussian");
        if(EndCaps == 1  &&  r9sup == 0 && nomFitMethode == "RooBifurcatedGauss") text->DrawLatex(0.16, 0.85, "Endcaps, Low R9, Bifurcated Gaussian");

	TLegend * leg = new TLegend(0.6,0.7,0.9,0.9,"","brNDC");
   	//leg->SetHeader("The Legend Title");
   	//leg->AddEntry(h1,"Histogram filled with random numbers","f");
   	//leg->AddEntry("f1","Function abs(#frac{sin(x)}{x})","l");
   	leg->SetTextSize(0.030);
	leg->SetFillColor(kWhite);
	leg->SetLineColor(kWhite);
	leg->SetShadowColor(kWhite);
	string temp = "";
	temp = Form("%d",FitRange-2); temp += "%";
	leg->AddEntry(PValuesOptRangeM2->GetName(),temp.c_str(),"lep");
	temp = Form("%d",FitRange-1); temp += "%";
	leg->AddEntry(PValuesOptRangeM1->GetName(),temp.c_str(),"lep");
	temp = Form("%d",FitRange); temp += "%";
	leg->AddEntry(PValuesOptRange->GetName(),temp.c_str(),"lep");
	temp = Form("%d",FitRange+1); temp += "%";
	leg->AddEntry(PValuesOptRangeP1->GetName(),temp.c_str(),"lep");
	temp = Form("%d",FitRange+2); temp += "%";
	leg->AddEntry(PValuesOptRangeP2->GetName(),temp.c_str(),"lep");
   	leg->Draw();


        
	gStyle->SetPadBorderMode(0);

        mg->GetYaxis()->SetRangeUser(yminPValues,ymaxPValues);
        mg->GetXaxis()->SetLimits(xminPValues,xmaxPValues);

        c1->SetTickx(1);
        c1->SetTicky(1);
        c1->SetGridx(1);
        c1->SetGridy(1);
        c1->Modified();
        c1->cd();
        c1->SetSelected(c1);
        c1->ToggleToolBar();

	//DossierEnregistrement = Form("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/ComparaisonFitsPValuesDiffMmumuJan/Olivier9/%s/MmumuSelection%s_%s/Fits/",PileUpVersion.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str());
	
	//if(StudyFolder == "Ranges") DossierEnregistrement = Form("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/%s/%s/%s/MmumuSelection%s_%s/%dPourcents/%s/CombinedGraphs/",scale.c_str(),DossierStudyFolder.c_str(),PileUpVersion.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),FitRange,ChainCracks.c_str());

	//if(StudyFolder == "FinalResults") DossierEnregistrement = Form("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/%s/%s/%s/MmumuSelection%s_%s/%s/%s/CombinedGraphs/",scale.c_str(),DossierStudyFolder.c_str(),PileUpVersion.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),ChainCracks.c_str(),SetOfCorrections.c_str());

	DossierEnregistrement = Form("Results_Mars_2012_v2/%s/%s/MmumuSelection%s_%s/%s/%s/%s/MZbinning_%s/%dPourcents/CombinedGraphs/",variableX.c_str(),PileUpVersion.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),ChainCracks.c_str(),SetOfCorrections.c_str(),isMCChain.c_str(),MZbinning.c_str(),FitRange);

	//FichierEnregistrement = Form("PValuesVsPt%s",compressFitMethodeName.c_str());



	FichierEnregistrement = Form("PValuesVsPt%s%s%s_%s_%s",compressFitMethodeName.c_str(),compressR9Name.c_str(),compressEndCapsName.c_str(),ChainCracks.c_str(),SetOfCorrections.c_str());

        enregistrementPlots(DossierEnregistrement, FichierEnregistrement, 2, 10000, c1);




	return 0;
}
