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


int PullsPloterF2_Mars_2012(string scale = "mmg_ik_MZ_Surface", string nomFitMethode = "RooLandauConvGaussian", string PileUpVersion = "2011", string LowMmumuLim = "40", string HightMmumuLim = "80", bool phiCracks = true, bool etaCracks = true, int FitRange = 0, string SetOfCorrections = "ETHZCorrections", string variableX = "Photon_Et", int isMC = 1, int FitRangeVal = 0, string MZbinning = "05GeV")
{
	gROOT->Reset();
	setTDRStyle();
	TGaxis::SetMaxDigits(3);
	gStyle->SetOptFit(0);

	double yminErecoOverEtrue, ymaxErecoOverEtrue;
	double xminErecoOverEtrue, xmaxErecoOverEtrue;
	
	yminErecoOverEtrue = -15.0;
	ymaxErecoOverEtrue = 15.0;
	xminErecoOverEtrue = 0.0;
        xmaxErecoOverEtrue = 100.0;


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
	string fichierPullsErrorXTXT = "";
	string fichierPullsErrorYTXT = "";
	string fichierPullsXTXT = "";
	string fichierPullsYTXT = "";
	
	double PullsErrorXTab[100] = {0.0};
	double PullsErrorYTab[100] = {0.0};
	double PullsXTab[100] = {0.0};
	double PullsYTab[100] = {0.0};
	vector<double> myvectorX;
	vector<double> myvectorY;

	double BinX[10];
	double BinY[10];
	double SuperBinX[60];
	double SuperBinY[60];


	int nb_lignes1 = 0;
	int nb_lignes2 = 0;
	string line;
	char buffer [1000];
	int MAXSIZE = 1000;
	int PosBin = 0;
        double NegaProb = 0.0;
	int iter = 0;
	int Inttemp = 0; 
	TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);

	string nomFichier = "";

	int rangeOpt = 0;
	if(FitRange == 0) rangeOpt = 1;

	for(int i = 0; i < 4; i++)
	{

		if(rangeOpt == 1){

		if(scale == "ComparaisonFitsErecoOverEtrueDiffMmumuJan" && nomFitMethode == "RooBifurcatedGauss") 
		{
			if(i == 0) FitRange = 65; //infEB
			if(i == 1) FitRange = 82; //supEB
			if(i == 2) FitRange = 72; //infEE
			if(i == 3) FitRange = 84; //supEE

		}
        	if(scale == "ComparaisonFits1overKrecoDiffMmumuJan" && nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC") 
		{
			if(i == 0) FitRange = 97;
                        if(i == 1) FitRange = 87;
                        if(i == 2) FitRange = 93;
                        if(i == 3) FitRange = 96;
		}
		
		if(scale == "ComparaisonFits1overKrecoDiffMmumuJan" && nomFitMethode == "RooBifurcatedGauss" && isMCChain == "Data")
                {
                        if(i == 0) FitRange = 86;
                        if(i == 1) FitRange = 86;
                        if(i == 2) FitRange = 93;
                        if(i == 3) FitRange = 96;
                }

		

		if(scale == "ComparaisonFitsErecoOverEtrueDiffMmumuJan" && nomFitMethode == "RooCrystalBall") 
                {
                        if(i == 0) FitRange = 84; 
                        if(i == 1) FitRange = 93; 
                        if(i == 2) FitRange = 72; 
                        if(i == 3) FitRange = 84; 

                }
                if(scale == "ComparaisonFits1overKrecoDiffMmumuJan" && nomFitMethode == "RooCrystalBall" && isMCChain == "MC") 
                {
                        if(i == 0) FitRange = 89; 
                        if(i == 1) FitRange = 82; 
                        if(i == 2) FitRange = 92; 
                        if(i == 3) FitRange = 89; 
                }

		if(scale == "ComparaisonFits1overKrecoDiffMmumuJan" && nomFitMethode == "RooCrystalBall" && isMCChain == "Data") 
                {
                        if(i == 0) FitRange = 86;
                        if(i == 1) FitRange = 88;
                        if(i == 2) FitRange = 92;
                        if(i == 3) FitRange = 87;
                }
	
		if(scale == "mmg_ik_MZ_Surface" && nomFitMethode == "RooBifurcatedGauss" && isMCChain == "MC") 
                {
                        if(i == 0) FitRange = 81; 
                        if(i == 1) FitRange = 89; 
                        if(i == 2) FitRange = 76; 
                        if(i == 3) FitRange = 80; 
                }    
     
                if(scale == "mmg_ik_MZ_Surface" && nomFitMethode == "RooCrystalBall" && isMCChain == "MC")
                {
                        if(i == 0) FitRange = 76; 
                        if(i == 1) FitRange = 90; 
                        if(i == 2) FitRange = 85; 
                        if(i == 3) FitRange = 87; 
                }

                if(scale == "mmg_ik_MZ_Surface" && (nomFitMethode == "RooLandau2" || nomFitMethode == "RooLandauConvGaussian") && isMCChain == "MC") 
                {
                        if(i == 0) FitRange = 81; 
                        if(i == 1) FitRange = 67; 
                        if(i == 2) FitRange = 80; 
                        if(i == 3) FitRange = 67; 
                }




		}


		FitRange += FitRangeVal;

		fichierTXT = Form("Results_Mars_2012_v2/%s/%s/MmumuSelection%s_%s/%s/%s/%s/MZbinning_%s/%dPourcents/Fits/",variableX.c_str(),PileUpVersion.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),ChainCracks.c_str(),SetOfCorrections.c_str(),isMCChain.c_str(),MZbinning.c_str(),FitRange);
		fichierTXT += Form("%s%s",compressScaleName.c_str(),compressFitMethodeName.c_str());


                if(i == 0) fichierTXT += Form("r9infEB/",compressFitMethodeName.c_str());
                if(i == 1) fichierTXT += Form("r9supEB/",compressFitMethodeName.c_str());
                if(i == 2) fichierTXT += Form("r9infEE/",compressFitMethodeName.c_str());
                if(i == 3) fichierTXT += Form("r9supEE/",compressFitMethodeName.c_str());


		
		for(int j = 0; j < n; j++)
		{

			for(int r = 0; r < 100; r++)
			{

				if(r<10) BinY[r] = 0.0;
				PullsErrorXTab[r] = 0.0;
                        	PullsErrorYTab[r] = 0.0;
                        	PullsXTab[r] = 0.0;
                        	PullsYTab[r] = 0.0;		

			}

			nb_lignes1 = 0;
			nb_lignes2 = 0;

			fichierPullsXTXT = fichierTXT + Form("PullsX%d.txt",j);
			fichierPullsYTXT = fichierTXT + Form("PullsY%d.txt",j);
        		
			//cout<<"fichierPullsXTXT = "<<fichierPullsXTXT<<endl;
        		ifstream monFlux1(fichierPullsXTXT.c_str());	
			ifstream monFlux2(fichierPullsYTXT.c_str());		


			while(!monFlux1.eof())
        		{
            			//fgets(buffer, MAXSIZE, monFlux1);
				getline (monFlux1,line);	
				nb_lignes1++;
			
        		}

			//cout<<endl<<"nb_lignes1 = "<<nb_lignes1<<endl;

			while(!monFlux2.eof())
                        {
                                //fgets(buffer, MAXSIZE, monFlux1);
                                getline (monFlux2,line);    
                                nb_lignes2++;
     
                        }	

			monFlux1.close();
			monFlux2.close();	

			ifstream monFlux3(fichierPullsXTXT.c_str());
			for(int k = 0; k < nb_lignes1-1; k++)
			{
				monFlux3 >> PullsXTab[k];
				myvectorX.push_back(PullsXTab[k]);
			} 

			monFlux3.close();

			ifstream monFlux4(fichierPullsYTXT.c_str());
                        for(int l = 0; l < nb_lignes2-1; l++)
                        {
                                monFlux4 >> PullsYTab[l];
				myvectorY.push_back(PullsYTab[k]);
                        }

                        monFlux4.close();

		
			
			TH1D *h1=new TH1D("h1","h1",10,-6.0,6.0);
		/*	
			for(int t = 0; t < nb_lignes1; t++)
        		{
                		NegaProb = PullsXTab[t] + 6.0;
                		for(int numBin1 = 0; numBin1 < 10; numBin1++)
                		{
                        		if(fabs(((1.2 * numBin1) - NegaProb)) < fabs(((1.2 * (numBin1 + 1)) - NegaProb)))
                        		{
                                		PosBin = numBin1;
                                		break;
                        		}
                		}
                		//cout<<"PosBin = "<<PosBin<<endl;
                		BinY[PosBin] += PullsYTab[t];
        		}

		*/	for(int t = 0; t < nb_lignes1-1; t++)
                        {
				iter = 0;
				//NegaProb = (PullsXTab[t]) + 6.0;
				//cout<<endl<<"PullsXTab[t] = "<<PullsXTab[t];
                                for(double  numBin1 = -6.0; numBin1 < 6.0; numBin1 += 1.2)
                                {
					//cout<<endl<<"numBin1 = "<<numBin1;
					//cout<<endl<<"(numBin1 - (1.2 /2.0)) = "<<(numBin1 - (1.2 /2.0));
					//cout<<endl<<"(numBin1 + (1.2 /2.0)) = "<<(numBin1 + (1.2 /2.0));
					//cout<<endl<<"iter1 = "<<iter;
					//if(fabs(NegaProb - (1.2*(numBin1+1))) < fabs(NegaProb - (1.2*(numBin1+2))))
					if(PullsXTab[t] > numBin1 && PullsXTab[t] < (numBin1 + 1.2))
					{
						//if((i == 3 && j == 0) || (i == 2 && j == 5)) cout<<endl<<"PullsXTab = "<<PullsXTab[t];	
						//cout<<endl<<"PullsXTab2[t] = "<<PullsXTab[t];
						//cout<<endl<<"NegaProb = "<<NegaProb;
						//cout<<endl<<"numBin1 = "<<numBin1;
						//cout<<endl<<"Diff = "<<fabs(NegaProb - (1.2*numBin1));
						//cout<<endl<<"Diff = "<<fabs(NegaProb - (1.2*(numBin1+1)))<<endl;
						//PosBin = numBin1;
						BinY[iter+1] += 1;
						break;
					}
					iter++;

                                }
				//BinY[iter+1] += 1;
                        }			

	
			for(int o = 0; o < 10; o++)
                        {
				h1->SetBinContent(o,BinY[o]);
				//cout<<endl<<"BinY[o] = "<<BinY[o];

			}

			

			h1->Draw("E");
			h1->GetXaxis()->SetTitle("#chi^{2} Pulls");
        		h1->GetXaxis()->SetLabelFont(42);
        		h1->GetXaxis()->SetTitleFont(42);
        		h1->GetXaxis()->SetLabelSize(0.03);
        		if(scale == "ComparaisonFitsErecoOverEtrueDiffMmumuJan") h1->GetYaxis()->SetTitle("Bins of s_{TRUE} = E_{RECO}/E_{TRUE} - 1 (%)");
        		if(scale == "ComparaisonFits1overKrecoDiffMmumuJan") h1->GetYaxis()->SetTitle("Bins of s_{RECO} = 1/k_{RECO} - 1 (%)");
        		h1->GetYaxis()->SetLabelFont(42);
        		h1->GetYaxis()->SetTitleOffset(1.24);
        		h1->GetYaxis()->SetTitleFont(42);
        		h1->GetYaxis()->SetLabelSize(0.03);
			h1->SetMarkerStyle(20);

			TF1 * f1 = new TF1("f1","gaus",-6,6);
			f1->FixParameter(1,0.0);
			f1->FixParameter(2,1.0);
			h1->Fit(f1,"b");
			f1->SetLineColor(kBlue);
			f1->SetLineWidth(2);
			f1->Draw("SAMES");	

			c1->Clear();
			h1->Draw("E");
			f1->Draw("SAMES");
			
			nomFichier = "Pulls1D";

			enregistrementPlots(fichierTXT, nomFichier, 2, j, c1);			
	
			//system("mkdir SuperTestDeLaMort");
			//c1->Print(Form("SuperTestDeLaMort/test%d_%d.png",i,j));
			c1->Clear();
			h1->Delete();
			f1->Delete();
		}

	}


	TH1D *h2=new TH1D("h2","h2",60,-6.0,6.0);

	for(int q = 0; q < myvectorX.size() ; q++)
        {
		iter = 0;
		//NegaProb = myvectorX[q] + 6.0;
		for(double numBin = -6.0; numBin < 6.0; numBin += 0.2)
		{
			//if(fabs(((0.2 * (numBin+1)) - NegaProb)) < fabs(((0.2 * (numBin + 2)) - NegaProb))) 
			if(myvectorX[q] > numBin && myvectorX[q] < (numBin + 0.2))
			{
				//PosBin = numBin;		
				//cout<<"fabs((0.2 * numBin - NegaProb)) = "<<fabs((0.2 * numBin - NegaProb))<<endl;
				//cout<<"fabs((0.2 * (numBin + 1)) - NegaProb)) = "<<fabs((0.2 * (numBin + 1) - NegaProb))<<endl;
				SuperBinY[iter+1] += 1;
				break;
			}
			iter++;
		}
		//cout<<"PosBin = "<<PosBin<<endl;
        	//SuperBinY[PosBin] += myvectorY[q];
        	//SuperBinY[iter+1] += 1;
	}


	for(int p = 0; p < 60; p++)
        {
		//cout<<SuperBinY[p]<<endl;
		Inttemp = p-1;
        	h2->SetBinContent(p,SuperBinY[p]);


        }




	h2->Draw("E");
        h2->GetXaxis()->SetTitle("#chi^{2} Pulls");
        h2->GetXaxis()->SetLabelFont(42);
        h2->GetXaxis()->SetTitleFont(42);
        h2->GetXaxis()->SetLabelSize(0.03);
	if(scale == "ComparaisonFitsErecoOverEtrueDiffMmumuJan") h2->GetYaxis()->SetTitle("Bins of s_{TRUE} = E_{RECO}/E_{TRUE} - 1 (%)");
        if(scale == "ComparaisonFits1overKrecoDiffMmumuJan") h2->GetYaxis()->SetTitle("Bins of s_{RECO} = 1/k_{RECO} - 1 (%)");
        h2->GetYaxis()->SetLabelFont(42);
        h2->GetYaxis()->SetTitleOffset(1.24);
        h2->GetYaxis()->SetTitleFont(42);
        h2->GetYaxis()->SetLabelSize(0.03);
	h2->SetMarkerStyle(20);

	TF1 * f2 = new TF1("f2","gaus",-6,6);
        f2->FixParameter(1,0.0);
        f2->FixParameter(2,1.0);
        h2->Fit(f2,"b");
        f2->SetLineColor(kBlue);
	f2->SetLineWidth(2);
        f2->Draw("SAMES");    
        c1->Clear();
        h2->Draw("E");
        f2->Draw("SAMES");                              

	string DossierEnregistrement = "";
	string FichierEnregistrement = "";

	//DossierEnregistrement = Form("/sps/cms/sgandurr/CMSSW_4_2_3_patch2/src/UserCode/IpnTreeProducer/%s/%s/%s/MmumuSelection%s_%s/CombinedGraphs/",scale.c_str(),DossierStudyFolder.c_str(),PileUpVersion.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str());
        


	if(rangeOpt == 1){ 


		DossierEnregistrement = Form("Results_Mars_2012_v2/%s/%s/MmumuSelection%s_%s/%s/%s/%s/MZbinning_%s/CombinedGraphs/",variableX.c_str(),PileUpVersion.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),ChainCracks.c_str(),SetOfCorrections.c_str(),isMCChain.c_str(),MZbinning.c_str());

        }



        FichierEnregistrement = Form("Pulls%s",compressFitMethodeName.c_str());
        if(rangeOpt == 1) FichierEnregistrement = Form("PullsRangeOpt%s_Val%d",compressFitMethodeName.c_str(),FitRangeVal);



        enregistrementPlots(DossierEnregistrement, FichierEnregistrement, 2, 10000, c1);	

 
        //c1->Print("SuperTestDeLaMort/testALL.png");
        c1->Clear();





	return 0;
}
