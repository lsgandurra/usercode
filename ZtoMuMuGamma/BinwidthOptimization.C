// -------------------------------------------------------
// program implemented by Louis Sgandurra (September 2012)
// -------------------------------------------------------

// for more informations about the algorithm see : http://toyoizumilab.brain.riken.jp/hideaki/res/histogram.html

#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
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
#include "TROOT.h"
#include "TRint.h"
#include "TMultiGraph.h"
#include "setTDRStyle.C"
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



int BinwidthOptimization(string sample = "", string variable = "", TString cuts = "", double rangeMin = -0.5, double rangeMax = 0.5, int nBinsMax = 1000)
{

	sample = "/sps/cms/sgandurr/CMSSW_4_2_8_patch7/src/Zmumugamma/RegressionMiniTrees_June_2012/MC_Fall11_Rochester/Surface_5GeV/miniTree_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1_40_80_1_v11_partALL.root";
	cuts = "isLooseMMG == 1 && Photon_isEB == 1";
	variable = "mmg_s";
	
	cout<<endl<<"// --------------------" <<endl
		  <<"// BinwidthOptimization" <<endl
		  <<"// --------------------" <<endl;

	cout<<endl<<"sample = "<<sample;
	cout<<endl<<"variable = "<<variable;
	cout<<endl<<"cuts = "<<cuts;
	cout<<endl<<"rangeMin = "<<rangeMin<<", rangeMax = "<<rangeMax;
	cout<<endl<<"nBinsMax = "<<nBinsMax<<endl;

	if(rangeMax <= rangeMin)
	{
		cout<<endl<<"error : rangeMax <= rangeMin"<<endl;
		return 0;
	}
	

        TChain * chain = new TChain("miniTree");
	
	chain->Add(Form("%s",sample.c_str()));
	
	TH1D * histo = 0;	
	
	double * kiTab = 0;
	double k = 0;
	double v = 0;
	double C = 0;
	double delta = 0;

	int nBinBest = 1;
	double CBest = 9999999;


	for(int nBins = 1; nBins <= nBinsMax; nBins++)
	{

		//cout<<endl<<"nBins = "<<nBins;
		histo = new TH1D("histo", "histo", nBins, rangeMin, rangeMax);	
		chain->Draw(Form("%s>>histo",variable.c_str()),cuts);

		kiTab = new double[nBins];
		
		k = 0;
		v = 0;

		for(int i = 0; i < nBins; i++)
		{
			kiTab[i] = histo->GetBinContent(i+1);	
			k += kiTab[i];	
		}
		
		k *= 1.0 / nBins;		
	
		for(int i = 0; i < nBins; i++)
		{
			v += (kiTab[i] - k) * (kiTab[i] - k);
		}

		v *= 1.0 / nBins;

		
		delta = (rangeMax - rangeMin) / (nBins * 1.0);

		C = (2 * k - v) / (delta * delta);

		if(C < CBest)
		{
			CBest = C;
			nBinBest = nBins;
		} 


		delete kiTab;
		kiTab = 0;
		histo->Delete();
        	histo = 0;

	}
	
	cout<<endl<<"CBest = "<<CBest;
	cout<<endl<<"nBinBest = "<<nBinBest<<endl;


	return 0;

}






 
