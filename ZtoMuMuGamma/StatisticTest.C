// -------------------------------------------------------
// program implemented by Louis Sgandurra (September 2012)
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
#include <sstream>
#pragma optimize 0

using namespace RooFit;
using namespace std;

void StatisticTest(string statTest, double & testEstimator, double & testPvalue, TString cut, TChain * chain, TString variable, TF1 * f, double RangeMin, double RangeMax) //, TH1D * hh, int nbOfBins)
{
	//for Kolmogorov-Smirnov test you have put statTest = "KS"
	//for Anderson-Darling statTest = "AD"


	// ------------------------------------------ //	
	// create a table of our variable (with cuts) //
	// ------------------------------------------ //

	TChain * ReducedChain = (TChain *) chain->CopyTree(cut);
	float value;
        ReducedChain->SetBranchAddress(variable,&value);

	int nEvents = ReducedChain->GetEntries();

	Double_t * sample = new Double_t[nEvents];

	//Double_t * sample2 = new Double_t[nbOfBins];

	for (int ievt = 0 ; ievt < nEvents; ievt++)
        {
                ReducedChain->GetEntry(ievt);
                //sample.push_back(value);
        	sample[ievt] = value;
	}

/*
	for (int ievt = 0 ; ievt < nbOfBins; ievt++)
	{
		sample2[ievt] = hh->GetBinContent(ievt + 1);
		//cout<<endl<<"sample2[ievt] = "<<sample2[ievt]<<endl;
	}
*/


	// -------------------------------------- //        
        // Transform the TF1 and create a GoFTest //
        // -------------------------------------- //	

	ROOT::Math::WrappedTF1 wf1(*f);

	ROOT::Math::GoFTest* goftest_1 = new ROOT::Math::GoFTest(nEvents, sample, wf1, ROOT::Math::GoFTest::kPDF, RangeMin, RangeMax);
	//ROOT::Math::GoFTest* goftest_1 = new ROOT::Math::GoFTest(nbOfBins, sample2, wf1, ROOT::Math::GoFTest::kCDF, RangeMin, RangeMax);

	// ---------------------------------------- //        
        // Return the test statistic and the pvalue //
        // ---------------------------------------- //


	if(statTest == "KS")
	{
		testEstimator = goftest_1-> KolmogorovSmirnovTest("t");
		testPvalue = goftest_1-> KolmogorovSmirnovTest();
	}

	if(statTest == "AD")
        {
                testEstimator = goftest_1-> AndersonDarlingTest("t");
                testPvalue = goftest_1-> AndersonDarlingTest();
        }


	delete sample;
	delete goftest_1;
	sample = 0;
	goftest_1 = 0;

}






 
