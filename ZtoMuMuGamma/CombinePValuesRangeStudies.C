//Script by Louis Sgandurra (January 2012)
//combine P-Value vs Ptgamma graphs, for different ranges of fit, in one graph.

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

int NbLignesFichier(string fichier);
void enregistrementPlots(string nomDossier, string nomFichier, int EndCaps, int iteration, TCanvas * c1);

int NbLignesFichier(string fichier)
{
        ifstream in(fichier.c_str()); //Ouverture en mode lecture de fichier

        string ligne; //Création d'une chaine de caracterece
        int nbLignes = 0;  

        while(std::getline(in, ligne)) nbLignes++;

        in.close(); //On ferme le fichier

        return nbLignes;
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


int CombinePValuesRangeStudiesF_Juin_2012(int EndCaps = 0, int r9sup = 1, string scale = "mmg_ik_MZ_Surface", string SetOfCorrections = "ETHZCorrections", string nomFitMethode = "RooCrystalBall", string PileUpVersion = "2011", string LowMmumuLim = "40", string HightMmumuLim = "80", bool phiCracks = true, bool etaCracks = true, string variableX = "Photon_Et", int isMC = 1, int LowPourcent = 80, int HighPourcent = 90, string MZbinning = "05GeV", string MuonCorrection = "NoMuonCorrection", string SurfaceMethod = "Fitted_Surface", string Category = "Vgamma8", string date = "16Jan")
{

	gROOT->Reset();
	setTDRStyle();
	TGaxis::SetMaxDigits(3);

	double xValueTab[6];
	double xErrorLTab[6];
	double xErrorRTab[6];

	double PValueTab1[6];
	double PValueTab2[6];
	double PValueTab3[6];
	double PValueTab4[6];
	double PValueTab5[6];
	double PValueTab6[6];
        double PValueTab7[6];
        double PValueTab8[6];
        double PValueTab9[6];
        double PValueTab10[6];
	double PValueTab11[6];

	double mean1 = 0;
	double mean2 = 0;
	double mean3 = 0;
	double mean4 = 0;
	double mean5 = 0;
	double mean6 = 0;
        double mean7 = 0;
        double mean8 = 0;
        double mean9 = 0;
        double mean10 = 0;
	double mean11 = 0;	


	
	double PValueErrorTab[6] = {0.0};
	
	double yminPValue, ymaxPValue;
	double xminPValue, xmaxPValue;
	
	yminPValue = 0.0;
	ymaxPValue = 1.0;
	xminPValue = 0.0;
        xmaxPValue = 60.0;

	if(variableX == "Photon_SC_Eta") xminPValue = -3.0;
        xmaxPValue = 100.0;
        if(variableX == "Photon_SC_rawEt") xmaxPValue = 100.0;
        if(variableX == "Photon_Et") xmaxPValue = 100.0;
        if(variableX == "Photon_E") xmaxPValue = 200.0;
        if(variableX == "Photon_SC_Eta") xmaxPValue = 3.0;
        if(variableX == "Photon_SC_brem") xmaxPValue = 15.0;



	string compressScaleName = "";
        if(scale == "Photon_E_o_MC_E") compressScaleName = "ErecoOverEtrue";
        if(scale == "ComparaisonFits1overKrecoDiffMmumuJan") compressScaleName = "1overKreco";
	if(scale == "mmg_ik_MZ_Surface") compressScaleName = "1overKreco";
	if(scale == "mmg_ik") compressScaleName = "1overKrecoClassical";
	if(scale == "mmg_s") compressScaleName = "1overKrecoClassical";

        string compressEndCapsName = "";
        if(EndCaps == 0) compressEndCapsName = "EB";
        if(EndCaps == 1) compressEndCapsName = "EE";

        string compressR9Name = "";
        if(r9sup == 0) compressR9Name = "r9inf";
        if(r9sup == 1) compressR9Name = "r9sup";
	if(r9sup == 2) compressR9Name = "r9All";
	
	string compressFitMethodeName = "";
	if(nomFitMethode == "RooCrystalBall") compressFitMethodeName = "CB";
	if(nomFitMethode == "RooBifurcatedGauss") compressFitMethodeName = "BFG";
	if(nomFitMethode == "RooLandau2") compressFitMethodeName = "Landau";
	if(nomFitMethode == "RooLandauConvGaussian") compressFitMethodeName = "LandauConvGaus";
	if(nomFitMethode == "RooVoigtian2") compressFitMethodeName = "Voigtian";

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
	string fichierPValueTXT = "";
	string fichierXValueTXT = "";
	string fichierXErrorLTXT = "";
	string fichierXErrorRTXT = "";

	int n = 0;


	for(int i = LowPourcent; i < HighPourcent; i++)
	{
		//fichierTXT = Form("Results_Avril_2012_v2/%s/%s/MmumuSelection%s_%s/%s/%s/%s/%s/MZbinning_%s/%s/%dPourcents/Fits/",variableX.c_str(),PileUpVersion.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),ChainCracks.c_str(),SetOfCorrections.c_str(),MuonCorrection.c_str(),isMCChain.c_str(),MZbinning.c_str(),SurfaceMethod.c_str(),i); 	
		fichierTXT = Form("Results_Juin_2012_v6_Approval/%s/%s/MmumuSelection%s_%s/%s/%s/%s/%s/MZbinning_%s/%s/%s/%s/%dPourcents/Fits/",variableX.c_str(),PileUpVersion.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),ChainCracks.c_str(),SetOfCorrections.c_str(),MuonCorrection.c_str(),isMCChain.c_str(),MZbinning.c_str(),SurfaceMethod.c_str(),Category.c_str(),date.c_str(),i);	


		fichierTXT += Form("%s%s%s%s/",compressScaleName.c_str(),compressFitMethodeName.c_str(),compressR9Name.c_str(),compressEndCapsName.c_str());

		fichierPValueTXT = fichierTXT + "PValueTab.txt";

		fichierXValueTXT = fichierTXT + "xValue.txt"; 
        	fichierXErrorLTXT = fichierTXT + "xErrorL.txt"; 
        	fichierXErrorRTXT = fichierTXT + "xErrorR.txt";


		cout<<endl<<"fichierXValueTXT = "<<fichierXValueTXT<<endl;			

        	//cout<<"fichierXErrorLTXT = "<<fichierXErrorLTXT<<endl;
		n = NbLignesFichier(fichierPValueTXT.c_str()); 
                cout<<"n = "<<n<<endl;
		ifstream monFlux1(fichierPValueTXT.c_str());
/*
		if(monFlux1)    //On teste si tout est OK.
		{
			cout << "OK : Ouverture du fichier" << endl;
		}
		else
		{
    			cout << "ERREUR: Impossible d'ouvrir le fichier." << endl;
		}
*/
		for(int j = 0; j < n; j++)
		{ 
			if(i == LowPourcent) monFlux1 >> PValueTab1[j];
			if(i == (LowPourcent+1)) monFlux1 >> PValueTab2[j];
			if(i == (LowPourcent+2)) monFlux1 >> PValueTab3[j];
			if(i == (LowPourcent+3)) monFlux1 >> PValueTab4[j];
			if(i == (LowPourcent+4)) monFlux1 >> PValueTab5[j];
			if(i == (LowPourcent+5)) monFlux1 >> PValueTab6[j];
			if(i == (LowPourcent+6)) monFlux1 >> PValueTab7[j];
			if(i == (LowPourcent+7)) monFlux1 >> PValueTab8[j];
			if(i == (LowPourcent+8)) monFlux1 >> PValueTab9[j];
			if(i == (LowPourcent+9)) monFlux1 >> PValueTab10[j];
			if(HighPourcent == 101 && i == (LowPourcent+10)) monFlux1 >> PValueTab11[j];
		}

		if(i == LowPourcent)
		{
			ifstream monFlux2(fichierXValueTXT.c_str());
			for(int k = 0; k < n; k++)
                	{
				monFlux2 >> xValueTab[k];
				//cout<<endl<<"xValueTab[k] = "<<xValueTab[k]<<endl;			
			}

			ifstream monFlux3(fichierXErrorLTXT.c_str());
                        for(int l = 0; l < n; l++)
                        {
                                monFlux3 >> xErrorLTab[l];
				//cout<<endl<<"xErrorLTab[l] = "<<xErrorLTab[l]<<endl;
                        }

			ifstream monFlux4(fichierXErrorRTXT.c_str());
                        for(int m = 0; m < n; m++)
                        {
                                monFlux4 >> xErrorRTab[m];
				//cout<<endl<<"xErrorRTab[m] = "<<xErrorRTab[m]<<endl;
                        }

		}

	}

	for(int q = 0; q < n; q++)
	{
		mean1 += PValueTab1[q];
		mean2 += PValueTab2[q];
                mean3 += PValueTab3[q];
                mean4 += PValueTab4[q];
                mean5 += PValueTab5[q];
                mean6 += PValueTab6[q];
		mean7 += PValueTab7[q];
                mean8 += PValueTab8[q];
                mean9 += PValueTab9[q];
                mean10 += PValueTab10[q];
		if(HighPourcent == 101) mean11 += PValueTab11[q];

	}

	mean1 *= 1.0 / n;
	mean2 *= 1.0 / n;
	mean3 *= 1.0 / n;
	mean4 *= 1.0 / n;
	mean5 *= 1.0 / n;
	mean6 *= 1.0 / n;
        mean7 *= 1.0 / n;
        mean8 *= 1.0 / n;
        mean9 *= 1.0 / n;
        mean10 *= 1.0 / n;
	mean11 *= 1.0 / n;


	TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);

	TMultiGraph *mg = new TMultiGraph();

/*
	TGraphAsymmErrors * PValueVsPt1 = new TGraphAsymmErrors(n,xValueTab, PValueTab1, xErrorLTab, xErrorRTab, PValueErrorTab, PValueErrorTab);
	TGraphAsymmErrors * PValueVsPt3 = new TGraphAsymmErrors(n,xValueTab, PValueTab3, xErrorLTab, xErrorRTab, PValueErrorTab, PValueErrorTab);
	TGraphAsymmErrors * PValueVsPt5 = new TGraphAsymmErrors(n,xValueTab, PValueTab5, xErrorLTab, xErrorRTab, PValueErrorTab, PValueErrorTab);
	TGraphAsymmErrors * PValueVsPt7 = new TGraphAsymmErrors(n,xValueTab, PValueTab7, xErrorLTab, xErrorRTab, PValueErrorTab, PValueErrorTab);
	TGraphAsymmErrors * PValueVsPt10 = new TGraphAsymmErrors(n,xValueTab, PValueTab10, xErrorLTab, xErrorRTab, PValueErrorTab, PValueErrorTab);
	TGraphAsymmErrors * PValueVsPtm1 = new TGraphAsymmErrors(n,xValueTab, PValueTabm1, xErrorLTab, xErrorRTab, PValueErrorTab, PValueErrorTab);
        TGraphAsymmErrors * PValueVsPtm3 = new TGraphAsymmErrors(n,xValueTab, PValueTabm3, xErrorLTab, xErrorRTab, PValueErrorTab, PValueErrorTab);
        TGraphAsymmErrors * PValueVsPtm5 = new TGraphAsymmErrors(n,xValueTab, PValueTabm5, xErrorLTab, xErrorRTab, PValueErrorTab, PValueErrorTab);
        TGraphAsymmErrors * PValueVsPtm7 = new TGraphAsymmErrors(n,xValueTab, PValueTabm7, xErrorLTab, xErrorRTab, PValueErrorTab, PValueErrorTab);
        TGraphAsymmErrors * PValueVsPtm10 = new TGraphAsymmErrors(n,xValueTab, PValueTabm10, xErrorLTab, xErrorRTab, PValueErrorTab, PValueErrorTab);
	TGraphAsymmErrors * PValueVsPtm13 = new TGraphAsymmErrors(n,xValueTab, PValueTabm13, xErrorLTab, xErrorRTab, PValueErrorTab, PValueErrorTab);
        TGraphAsymmErrors * PValueVsPtm15 = new TGraphAsymmErrors(n,xValueTab, PValueTabm15, xErrorLTab, xErrorRTab, PValueErrorTab, PValueErrorTab);
        TGraphAsymmErrors * PValueVsPtm18 = new TGraphAsymmErrors(n,xValueTab, PValueTabm18, xErrorLTab, xErrorRTab, PValueErrorTab, PValueErrorTab);
        TGraphAsymmErrors * PValueVsPtm20 = new TGraphAsymmErrors(n,xValueTab, PValueTabm20, xErrorLTab, xErrorRTab, PValueErrorTab, PValueErrorTab);
        TGraphAsymmErrors * PValueVsPt0 = new TGraphAsymmErrors(n,xValueTab, PValueTab0, xErrorLTab, xErrorRTab, PValueErrorTab, PValueErrorTab);

*/

	TGraph * PValueVsPt1 = new TGraph(n,xValueTab, PValueTab1);
	TGraph * PValueVsPt2 = new TGraph(n,xValueTab, PValueTab2);
	TGraph * PValueVsPt3 = new TGraph(n,xValueTab, PValueTab3);
	TGraph * PValueVsPt4 = new TGraph(n,xValueTab, PValueTab4);
	TGraph * PValueVsPt5 = new TGraph(n,xValueTab, PValueTab5);
	TGraph * PValueVsPt6 = new TGraph(n,xValueTab, PValueTab6);
	TGraph * PValueVsPt7 = new TGraph(n,xValueTab, PValueTab7);
	TGraph * PValueVsPt8 = new TGraph(n,xValueTab, PValueTab8);
	TGraph * PValueVsPt9 = new TGraph(n,xValueTab, PValueTab9);
	TGraph * PValueVsPt10 = new TGraph(n,xValueTab, PValueTab10);
	if(HighPourcent == 101) TGraph * PValueVsPt11 = new TGraph(n,xValueTab, PValueTab11);

        c1->ToggleEventStatus();

	mg->Add(PValueVsPt1);
	mg->Add(PValueVsPt2);
	mg->Add(PValueVsPt3);
	mg->Add(PValueVsPt4);
	mg->Add(PValueVsPt5);
	mg->Add(PValueVsPt6);
	mg->Add(PValueVsPt7);
	mg->Add(PValueVsPt8);
	mg->Add(PValueVsPt9);
	mg->Add(PValueVsPt10);
	if(HighPourcent == 101) mg->Add(PValueVsPt11);
	
	mg->Draw("AP");

        if(variableX == "Photon_SC_rawEt") mg->GetXaxis()->SetTitle("P_{T RAW}^{#gamma}");
        if(variableX == "Photon_Et") mg->GetXaxis()->SetTitle("P_{T}^{#gamma}");
        if(variableX == "Photon_E") mg->GetXaxis()->SetTitle("E^{#gamma}");
        if(variableX == "Photon_SC_Eta") mg->GetXaxis()->SetTitle("#eta");
        if(variableX == "Photon_SC_brem") mg->GetXaxis()->SetTitle("#sigma_{#phi}/#sigma_{#eta}");

	mg->GetXaxis()->SetLabelFont(42);
        mg->GetXaxis()->SetTitleFont(42);
        mg->GetXaxis()->SetLabelSize(0.03);
        mg->GetYaxis()->SetTitle("P-Value");
        mg->GetYaxis()->SetLabelFont(42);
        mg->GetYaxis()->SetTitleOffset(1.24);
        mg->GetYaxis()->SetTitleFont(42);
        mg->GetYaxis()->SetLabelSize(0.03);
	//mg->SetTitle("");
        //mg->SetMarkerColor(4);
        //mg->SetMarkerStyle(21);
        //mg->SetMarkerSize(0.6);


	PValueVsPt1->SetLineColor(629);
	PValueVsPt1->SetMarkerColor(629);
	PValueVsPt1->SetMarkerStyle(20);
	PValueVsPt1->SetMarkerSize(1.15);
	PValueVsPt1->SetName("PValueVsPt1");
        PValueVsPt2->SetLineColor(629);
        PValueVsPt2->SetMarkerColor(629);
	PValueVsPt2->SetMarkerStyle(24);
	PValueVsPt2->SetMarkerSize(1.15);
	PValueVsPt2->SetName("PValueVsPt2");
	PValueVsPt3->SetLineColor(397);
	PValueVsPt3->SetMarkerColor(397);
	PValueVsPt3->SetMarkerStyle(21);
	PValueVsPt3->SetMarkerSize(1.15);
	PValueVsPt3->SetName("PValueVsPt3");
        PValueVsPt4->SetLineColor(397);
	PValueVsPt4->SetMarkerColor(397);
	PValueVsPt4->SetMarkerStyle(25);
	PValueVsPt4->SetMarkerSize(1.15);
	PValueVsPt4->SetName("PValueVsPt4");
        PValueVsPt5->SetLineColor(413);
	PValueVsPt5->SetMarkerColor(413);
	PValueVsPt5->SetMarkerStyle(22);
	PValueVsPt5->SetMarkerSize(1.15);
	PValueVsPt5->SetName("PValueVsPt5");
	PValueVsPt6->SetLineColor(413);
        PValueVsPt6->SetMarkerColor(413);
        PValueVsPt6->SetMarkerStyle(26);
        PValueVsPt6->SetMarkerSize(1.15);
        PValueVsPt6->SetName("PValueVsPt6");
	        PValueVsPt7->SetLineColor(597);
	        PValueVsPt7->SetMarkerColor(597);
	        PValueVsPt7->SetMarkerStyle(23);
	        PValueVsPt7->SetMarkerSize(1.15);
	        PValueVsPt7->SetName("PValueVsPt7");
	        PValueVsPt8->SetLineColor(597);
	        PValueVsPt8->SetMarkerColor(597);
	        PValueVsPt8->SetMarkerStyle(32);
	        PValueVsPt8->SetMarkerSize(1.15);
	        PValueVsPt8->SetName("PValueVsPt8");
	        PValueVsPt9->SetLineColor(613);
	        PValueVsPt9->SetMarkerColor(613);
	        PValueVsPt9->SetMarkerStyle(34);
	        PValueVsPt9->SetMarkerSize(1.15);
	        PValueVsPt9->SetName("PValueVsPt9");
	        PValueVsPt10->SetLineColor(613);
	        PValueVsPt10->SetMarkerColor(613);
	        PValueVsPt10->SetMarkerStyle(28);
	        PValueVsPt10->SetMarkerSize(1.15);
	        PValueVsPt10->SetName("PValueVsPt10");
		if(HighPourcent == 101)
                {
                        PValueVsPt11->SetLineColor(1);
                        PValueVsPt11->SetMarkerColor(1);
                        PValueVsPt11->SetMarkerStyle(29);
                        PValueVsPt11->SetMarkerSize(1.15);
                        PValueVsPt11->SetName("PValueVsPt11");
                }


	TLegend * leg = new TLegend(0.65,0.3,0.9,0.9,"","brNDC");
        //leg->SetHeader("The Legend Title");
        //leg->AddEntry(h1,"Histogram filled with random numbers","f");
        //leg->AddEntry("f1","Function abs(#frac{sin(x)}{x})","l");
        leg->SetTextSize(0.023);
        leg->SetFillColor(kWhite);
        leg->SetLineColor(kWhite);
        leg->SetShadowColor(kWhite);
	string temp = "";
	temp = Form("%d",LowPourcent); temp += "%, m = "; temp += Form("%f",mean1);
        leg->AddEntry(PValueVsPt1->GetName(),temp.c_str(),"p");
	temp = Form("%d",LowPourcent+1); temp += "%, m = "; temp += Form("%f",mean2);
        leg->AddEntry(PValueVsPt2->GetName(),temp.c_str(),"p");
	temp = Form("%d",LowPourcent+2); temp += "%, m = "; temp += Form("%f",mean3);
        leg->AddEntry(PValueVsPt3->GetName(),temp.c_str(),"p");
	temp = Form("%d",LowPourcent+3); temp += "%, m = "; temp += Form("%f",mean4);
        leg->AddEntry(PValueVsPt4->GetName(),temp.c_str(),"p");
	temp = Form("%d",LowPourcent+4); temp += "%, m = "; temp += Form("%f",mean5);
        leg->AddEntry(PValueVsPt5->GetName(),temp.c_str(),"p");
	temp = Form("%d",LowPourcent+5); temp += "%, m = "; temp += Form("%f",mean6);
        leg->AddEntry(PValueVsPt6->GetName(),temp.c_str(),"p");	
		temp = Form("%d",LowPourcent+6); temp += "%, m = "; temp += Form("%f",mean7);
        	leg->AddEntry(PValueVsPt7->GetName(),temp.c_str(),"p");
		temp = Form("%d",LowPourcent+7); temp += "%, m = "; temp += Form("%f",mean8);
        	leg->AddEntry(PValueVsPt8->GetName(),temp.c_str(),"p");
		temp = Form("%d",LowPourcent+8); temp += "%, m = "; temp += Form("%f",mean9);
        	leg->AddEntry(PValueVsPt9->GetName(),temp.c_str(),"p");
		temp = Form("%d",LowPourcent+9); temp += "%, m = "; temp += Form("%f",mean10);
        	leg->AddEntry(PValueVsPt10->GetName(),temp.c_str(),"p");
		if(HighPourcent == 101)
                {
                        temp = Form("%d",LowPourcent+10); temp += "%, m = "; temp += Form("%f",mean11);
                        leg->AddEntry(PValueVsPt11->GetName(),temp.c_str(),"p");
                }
	leg->Draw();







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
      
	if(nomFitMethode == "RooCrystalBall")
	{ 
		if(EndCaps == 0 && r9sup == 0) text->DrawLatex(0.16, 0.85, "Crystal Ball, Barrel, Low r9");
       		if(EndCaps == 0 && r9sup == 1) text->DrawLatex(0.16, 0.85, "Crystal Ball, Barrel, High r9");
       		if(EndCaps == 1 && r9sup == 0) text->DrawLatex(0.16, 0.85, "Crystal Ball, Endcaps, Low r9");
        	if(EndCaps == 1 && r9sup == 1) text->DrawLatex(0.16, 0.85, "Crystal Ball, Endcaps, High r9");
		if(EndCaps == 0 && r9sup == 2) text->DrawLatex(0.16, 0.85, "Crystal Ball, Barrel, All r9");
		if(EndCaps == 1 && r9sup == 2) text->DrawLatex(0.16, 0.85, "Crystal Ball, Endcaps, All r9");
	}

	if(nomFitMethode == "RooBifurcatedGauss")
        {
                if(EndCaps == 0 && r9sup == 0) text->DrawLatex(0.16, 0.85, "Bifurcated Gaussian, Barrel, Low r9");
                if(EndCaps == 0 && r9sup == 1) text->DrawLatex(0.16, 0.85, "Bifurcated Gaussian, Barrel, High r9");
                if(EndCaps == 1 && r9sup == 0) text->DrawLatex(0.16, 0.85, "Bifurcated Gaussian, Endcaps, Low r9");
                if(EndCaps == 1 && r9sup == 1) text->DrawLatex(0.16, 0.85, "Bifurcated Gaussian, Endcaps, High r9");
        	if(EndCaps == 0 && r9sup == 2) text->DrawLatex(0.16, 0.85, "Bifurcated Gaussian, Barrel, All r9");
                if(EndCaps == 1 && r9sup == 2) text->DrawLatex(0.16, 0.85, "Bifurcated Gaussian, Endcaps, All r9");
	}	
 

	gStyle->SetPadBorderMode(0);

        mg->GetYaxis()->SetRangeUser(yminPValue,ymaxPValue);
        mg->GetXaxis()->SetLimits(xminPValue,xmaxPValue);

        c1->SetTickx(1);
        c1->SetTicky(1);
        c1->SetGridx(1);
        c1->SetGridy(1);
        c1->Modified();
        c1->cd();
        c1->SetSelected(c1);
        c1->ToggleToolBar();

	//DossierEnregistrement = Form("Results_Avril_2012_v2/%s/%s/MmumuSelection%s_%s/%s/%s/%s/%s/MZbinning_%s/%s/CombinedGraphs/",variableX.c_str(),PileUpVersion.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),ChainCracks.c_str(),SetOfCorrections.c_str(),MuonCorrection.c_str(),isMCChain.c_str(),MZbinning.c_str(),SurfaceMethod.c_str());
	DossierEnregistrement = Form("Results_Juin_2012_v6_Approval/%s/%s/MmumuSelection%s_%s/%s/%s/%s/%s/MZbinning_%s/%s/%s/%s/CombinedGraphs/",variableX.c_str(),PileUpVersion.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),ChainCracks.c_str(),SetOfCorrections.c_str(),MuonCorrection.c_str(),isMCChain.c_str(),MZbinning.c_str(),SurfaceMethod.c_str(),Category.c_str(),date.c_str());
	FichierEnregistrement = Form("%s_PValueVsPt%s%s%s_%d_%d",compressScaleName.c_str(),compressFitMethodeName.c_str(),compressR9Name.c_str(),compressEndCapsName.c_str(),LowPourcent,HighPourcent-1);


        enregistrementPlots(DossierEnregistrement, FichierEnregistrement, 2, 10000, c1);




	return 0;
}
