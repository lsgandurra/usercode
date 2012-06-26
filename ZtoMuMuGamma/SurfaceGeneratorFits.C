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

string DoubleToString(double x);

Double_t effSigma(TH1 * hist);

double SigmaR(TF1* crystalBall, double Xmin, double Xmax);

double SigmaL(TF1* crystalBall, double Xmin, double Xmax);

Double_t chiSquare(RooPlot* plot_, char* pdfname, char* histname, int nFitParam, double* JanChi2, double* DegreesOfFreedom, double* pValue, int* fewBins);

void RooCrystalBall(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, double MinVar1, double MaxVar1, double MinVar2, double MaxVar2, string variableX, int isMC);


void RooBifurcatedGauss(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, double MinVar1, double MaxVar1, double MinVar2, double MaxVar2, string variableX, int isMC);


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

string DoubleToString(double x)
{

        std::string s;
        {
                std::ostringstream oss;
                oss << x;
                s = oss.str();
        }
        std::cout << "x = " << x << " s = " << s << std::endl;

        return s;
}

Double_t effSigma(TH1 * hist)
{

  TAxis *xaxis = hist->GetXaxis();
  Int_t nb = xaxis->GetNbins();
  if(nb < 10) {
    cout << "effsigma: Not a valid histo. nbins = " << nb << endl;
    return 0.;
  }

  Double_t bwid = xaxis->GetBinWidth(1);
  if(bwid == 0) {
    cout << "effsigma: Not a valid histo. bwid = " << bwid << endl;
    return 0.;
  }
  Double_t xmax = xaxis->GetXmax();
  Double_t xmin = xaxis->GetXmin();
  Double_t ave = hist->GetMean();
  Double_t rms = hist->GetRMS();

  Double_t total=0.;
  for(Int_t i=0; i<nb+2; i++) {
    total+=hist->GetBinContent(i);
  }
  if(total < 100.) {
    cout << "effsigma: Too few entries " << total << endl;
    return 0.;
  }
  Int_t ierr=0;
  Int_t ismin=999;

  Double_t rlim=0.683*total;
  Int_t nrms=rms/(bwid);    // Set scan size to +/- rms
  if(nrms > nb/10) nrms=nb/10; // Could be tuned...

  Double_t widmin=9999999.;
  for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
    Int_t ibm=(ave-xmin)/bwid+1+iscan;
    Double_t x=(ibm-0.5)*bwid+xmin;
    Double_t xj=x;
    Double_t xk=x;
    Int_t jbm=ibm;
    Int_t kbm=ibm;
    Double_t bin=hist->GetBinContent(ibm);
    total=bin;
    for(Int_t j=1;j<nb;j++){
      if(jbm < nb) {
        jbm++;
        xj+=bwid;
        bin=hist->GetBinContent(jbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
      if(kbm > 0) {
        kbm--;
        xk-=bwid;
        bin=hist->GetBinContent(kbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
    }
    Double_t dxf=(total-rlim)*bwid/bin;
    Double_t wid=(xj-xk+bwid-dxf)*0.5;
    if(wid < widmin) {
      widmin=wid;
      ismin=iscan;
    }
  }
  if(ismin == nrms || ismin == -nrms) ierr=3;
  if(ierr != 0) cout << "effsigma: Error of type " << ierr << endl;

  return widmin;

}

double SigmaR(TF1* crystalBall, double Xmin, double Xmax)
{

        double y = crystalBall->GetMaximum(Xmin, Xmax) * 1.0/exp(1.0);
        double MaxX = crystalBall->GetMaximumX(Xmin, Xmax);
        double sigmaR = crystalBall->GetX(y, MaxX, Xmax) - MaxX;

        return sigmaR;
}

double SigmaL(TF1* crystalBall, double Xmin, double Xmax)
{

        double y = crystalBall->GetMaximum(Xmin, Xmax) * 1.0/exp(1.0);
        double MaxX = crystalBall->GetMaximumX(Xmin, Xmax);
        double sigmaL = MaxX - crystalBall->GetX(y, Xmin, MaxX);

        return sigmaL;
}

Double_t chiSquare(RooPlot* plot_, char* pdfname, char* histname, int nFitParam, double* JanChi2, double* DegreesOfFreedom, double* pValue, int* fewBins)
{
  // Calculate the chi^2/NDOF of this curve with respect to the histogram
  // 'hist' accounting nFitParam floating parameters in case the curve
  // was the result of a fit

  // Find curve object
  RooCurve* curve = (RooCurve*) plot_->findObject(pdfname, RooCurve::Class());
  //RooCurve* curve = plot_->getCurve(pdfname);  
  //curve->Print();

  if (!curve) {
    cout<<endl << "cit::RooChi2Calculator(plotname=" << plot_->GetName()
         << ")::chiSquare(..) cannot find curve" << endl ;
    return 0 ;
  }

  // Find histogram object
  RooHist* hist = (RooHist*) plot_->findObject(histname, RooHist::Class()) ;
  //RooHist* hist = plot_->getHist(histname);

  if (!hist) {
    cout<<endl << "cit::RooChi2Calculator(plotname=" << plot_->GetName()
         << ")::chiSquare(..) cannot find histogram" << endl ;
    return 0 ;
  }


  Int_t i,np = hist->GetN() ;
  Double_t x,y,/*eyl,eyh,*/ xl,xh ;

  // Find starting and ending bin of histogram based on range of RooCurve
  Double_t xstart,xstop ;

#if ROOT_VERSION_CODE >= ROOT_VERSION(4,0,1)
  curve->GetPoint(0,xstart,y) ;
  curve->GetPoint(curve->GetN()-1,xstop,y) ;
#else
  const_cast<RooCurve*>(curve)->GetPoint(0,xstart,y) ;
  const_cast<RooCurve*>(curve)->GetPoint(curve->GetN() - 1,xstop,y) ;
#endif

  Int_t nbin(0) ;

  Double_t chisq(0) ;
  for (i=0 ; i<np ; i++) {

    // Retrieve histogram contents
    hist->GetPoint(i,x,y) ;
    xl = x - hist->GetEXlow()[i] ;
    xh = x + hist->GetEXhigh()[i] ;
    // eyl = hist->GetEYlow()[i] ;
    // eyh = hist->GetEYhigh()[i] ;

    // Check if the whole bin is in range of curve
    if (xl < xstart || xstop < xh) continue ;

    if(y != 0 && y < 35.0)
    {
        cout<<endl<<"Trop peu d'entree : "<<y<<" dans le bin : "<<i<<"  >>>Need to reduce the binning for the p-value calculation!"<<endl;
        *fewBins = 1;
        break;

    }
    else *fewBins = 0;

    nbin++ ;

    // Integrate function over this bin.
    // Start a hack to work around a bug in RooCurve::interpolate
    // that sometimes gives a wrong result.
    Double_t avg = curve->average(xl, xh);
    Double_t avg2 = 0.5 * (curve->average(xl, x) + curve->average(x, xh));
    if (avg + avg2 > 0 &&
        (avg2 - avg) / (avg2 + avg) > 0.1) {
      avg = curve->interpolate(x);
    }
    // End of hack around the bug in RooCurve::interpolate

    // JV: Adjust observed and expected number of events for bin width to represent
    // number of events.
    Double_t norm = (xh - xl) / plot_->getFitRangeBinW();
    y *= norm;
    avg *= norm;

    if (avg < 5.) {
      cout << "cit::RooChi2Calculator(plotname=" << plot_->GetName()
                            << ")::chiSquare(..) expectation in bin "
                            << i << " is " << avg << " < 5!" << endl ;
    }

    // JV: Use the expected number of events for the y uncertainty,
    // See (33.34) of http://pdg.lbl.gov/2011/reviews/rpp2011-rev-statistics.pdf

    // Add pull^2 to chisq
    if (avg != 0) {
      Double_t resid = y - avg;
      chisq += (resid * resid / avg) ;
    }
  }

  // Return chisq/nDOF 
  *JanChi2 = chisq / (nbin - nFitParam);
  *DegreesOfFreedom = (nbin - nFitParam);
  *pValue =  TMath::Prob(chisq, nbin - nFitParam);

  return chisq / (nbin - nFitParam) ;
}

void RooCrystalBall(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, double MinVar1, double MaxVar1, double MinVar2, double MaxVar2, string variableX, int isMC)
{

        RooRealVar * MuonL_Pt = new RooRealVar("MuonL_Pt", "MuonL_Pt", 0.0, 200.0, "");
	RooRealVar * MuonS_Pt = new RooRealVar("MuonS_Pt", "MuonS_Pt", 0.0, 200.0, "");
	RooRealVar * isMM = new RooRealVar("isMM", "isMM", 0.0, 2.0, ""); 
        RooRealVar * Mmumu = new RooRealVar("Mmumu", "Mmumu", 80.0, 100.0, ""); 

        RooRealVar * weight_pileUp = new RooRealVar("weight_pileUp", "weight_pileUp", 0.0, 10000000);

        RooArgSet *ntplVars = new RooArgSet(*MuonL_Pt,*MuonS_Pt,*isMM,*Mmumu);

	RooDataSet *Data = 0;
        if(isMC == 1) ntplVars->add(*weight_pileUp);
	if(isMC == 1) Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars, "", "weight_pileUp");
	if(isMC == 0) Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars, "");
	RooDataSet* Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);


	RooRealVar * Eraw_CB_m0 = new RooRealVar("Eraw_CB_m0", "CB #Delta m_{0}", 91.187,91.187-rms/2.0,91.187+rms/2.0, "GeV");
        RooRealVar * Eraw_CB_sigma = new RooRealVar("Eraw_CB_sigma", "CB ", rms, 0.0, 2 * rms, "GeV");
        RooRealVar * Eraw_CB_alpha = new RooRealVar("Eraw_CB_alpha", "CB #alpha", 1.0, 0.0, 10.0);
        RooRealVar * Eraw_CB_n = new RooRealVar("Eraw_CB_n", "CB n", 5.0, 1.0, 20.0);

        RooCBShape * Eraw_CrystalBall = new RooCBShape("Eraw_CrystalBall","mmg_ik_MZ_Surface_CrystalBall", *Mmumu, *Eraw_CB_m0, *Eraw_CB_sigma, *Eraw_CB_alpha, *Eraw_CB_n);



        int fewBins = 1;
        Double_t Chi2J;
	Erawframe->Clear();
	//Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);
        //Eraw_CrystalBall = new RooBifurGauss("Eraw_CrystalBall","mmg_ik_MZ_Surface_BifurGauss",*mmg_ik_MZ_Surface,*BifurGauss_mean,*BifurGauss_sigmaR,*BifurGauss_sigmaL);

        Erawframe = Mmumu->frame();
        //Data_subset->plotOn(Erawframe,Name("myhist"),Binning(RightBinning));
        Data_subset->plotOn(Erawframe,Name("myhist"));
	Eraw_CrystalBall->fitTo(*Data_subset, Range(RangeMin, RangeMax));
        Eraw_CrystalBall->plotOn(Erawframe,Name("mycurve"));
        Erawframe->Draw();
	Chi2J = chiSquare(Erawframe,(char *)"mycurve",(char *)"myhist",4, JanChi2, DegreesOfFreedom, pValue, &fewBins);

	f = Eraw_CrystalBall->asTF( RooArgList(*Mmumu) );
        *mean_value  = Eraw_CB_m0->getVal();
        *mean_error  = Eraw_CB_m0->getError();
        *sigmaEff_value = effSigma(hh);
        *sigma_value = Eraw_CB_sigma->getVal();
        *sigmaR_value = SigmaR(f, 0.0, 2.0);
        *sigmaL_value = SigmaL(f, 0.0, 2.0);
        *ChiSquare = Erawframe->chiSquare();


	double entries = hh->GetEntries();

        int entriesInt = (int) entries;

        double fxmax = f->GetMaximumX(0.8,1.2,1.E-10,100,false);
        //cout<<"////////// ---- FXMAX BG = "<<f->GetMaximumX(0.8,1.2,1.E-10,100,false)<<" ---- //////////"<<endl;

        TLatex latexLabel;
        latexLabel.SetTextSize(0.026);
        latexLabel.SetNDC();
        latexLabel.DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        if(isMC == 1) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4,66892 fb^{-1}");
        //if(EndCaps == 0) latexLabel.DrawLatex(0.16, 0.85, "Barrel");
        //if(EndCaps == 1) latexLabel.DrawLatex(0.16, 0.85, "EndCaps");
        string * tempLegChain(0);
        tempLegChain = new string;
        *tempLegChain = DoubleToString(MinVar1) + " < MuonL_Pt  < "  + DoubleToString(MaxVar1);
	string * tempLegChain2(0);
	tempLegChain2 = new string;
	*tempLegChain2 = DoubleToString(MinVar2) + " < MuonS_Pt  < "  + DoubleToString(MaxVar2);
        latexLabel.DrawLatex(0.16, 0.80, tempLegChain->c_str());
	latexLabel.DrawLatex(0.16, 0.75, tempLegChain2->c_str());

        latexLabel.DrawLatex(0.16, 0.70, Form("Entries = %d",entriesInt));
        cout<<endl<<"entries = "<<entries<<endl;
        cout<<endl<<"Mean = "<<hh->GetMean()<<endl;
        latexLabel.DrawLatex(0.16, 0.65, Form("Mean = %f +/- %f",hh->GetMean(), hh->GetMeanError()));
        latexLabel.DrawLatex(0.16, 0.60, Form("RMS = %f +/- %f",hh->GetRMS(), hh->GetRMSError()));
	latexLabel.DrawLatex(0.61, 0.90, Form("#color[4]{m_{0} = %f +/- %f}",Eraw_CB_m0->getVal(), Eraw_CB_m0->getError()));
        latexLabel.DrawLatex(0.61, 0.85, Form("#color[4]{#sigma = %f +/- %f}",Eraw_CB_sigma->getVal(), Eraw_CB_sigma->getError()));
        latexLabel.DrawLatex(0.61, 0.80, Form("#color[4]{#alpha = %f +/- %f}",Eraw_CB_alpha->getVal(), Eraw_CB_alpha->getError()));
        latexLabel.DrawLatex(0.61, 0.75, Form("#color[4]{n = %f +/- %f}",Eraw_CB_n->getVal(), Eraw_CB_n->getError()));
        latexLabel.DrawLatex(0.61, 0.70, Form("#color[4]{Default #chi^{2}/ndf = %f}",Erawframe->chiSquare()));
        latexLabel.DrawLatex(0.61, 0.65, Form("#color[4]{Jan #chi^{2}/ndf = %f}",Chi2J));
        latexLabel.DrawLatex(0.61, 0.60, Form("#color[4]{p-value = %0.6e}",*pValue));


        enregistrementPlots(nomDossier, nomFichier, 2, j, c2);

        c2->Clear();
	

	
	delete tempLegChain; tempLegChain = 0;
	delete tempLegChain2; tempLegChain2 = 0;
	delete MuonL_Pt; MuonL_Pt = 0; 
	delete MuonS_Pt; MuonS_Pt = 0;
	delete Mmumu; Mmumu = 0;
	delete weight_pileUp; weight_pileUp = 0;
	delete ntplVars; ntplVars = 0;
	delete Data; Data = 0;
	delete Data_subset; Data_subset = 0;
	delete Eraw_CB_m0; Eraw_CB_m0 = 0;
	delete Eraw_CB_sigma; Eraw_CB_sigma = 0;
	delete Eraw_CB_alpha; Eraw_CB_alpha = 0;
	delete Eraw_CB_n; Eraw_CB_n = 0;
	delete Eraw_CrystalBall; Eraw_CrystalBall = 0;

	


	return ;

}



void RooBifurcatedGauss(double * mean_value, double * mean_error, Double_t * sigmaEff_value, double * sigma_value, double * sigma_value_error, double * sigmaR_value, double * sigmaL_value, double * ChiSquare, TH1D * hh, int j, TString temp, TTree* Tree_Data, double mean, double rms, TCanvas * c2, RooPlot* Erawframe, TF1 * f, string nomDossier, string nomFichier, double RangeMin, double RangeMax, double * JanChi2, double * DegreesOfFreedom, double * pValue, double MinVar1, double MaxVar1, double MinVar2, double MaxVar2, string variableX, int isMC)
{

        RooRealVar * MuonL_Pt = new RooRealVar("MuonL_Pt", "MuonL_Pt", 0.0, 200.0, "");
	RooRealVar * MuonS_Pt = new RooRealVar("MuonS_Pt", "MuonS_Pt", 0.0, 200.0, "");
	RooRealVar * isMM = new RooRealVar("isMM", "isMM", 0.0, 2.0, "");
	RooRealVar * Mmumu = new RooRealVar("Mmumu", "Mmumu", 80.0, 100.0, "");

        RooRealVar * weight_pileUp = new RooRealVar("weight_pileUp", "weight_pileUp", 0.0, 10000000);

        RooArgSet *ntplVars = new RooArgSet(*MuonL_Pt,*MuonS_Pt,*isMM,*Mmumu);

	RooDataSet *Data = 0;
        if(isMC == 1) ntplVars->add(*weight_pileUp);
	if(isMC == 1) Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars, "", "weight_pileUp");
	if(isMC == 0) Data = new RooDataSet("Data", "Data", Tree_Data, *ntplVars, "");
	RooDataSet* Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);
	RooRealVar * BifurGauss_mean = new RooRealVar("BifurGauss_mean","BifurGauss_mean",91.187,91.187-rms,91.187+rms);
        RooRealVar * BifurGauss_sigmaR = new RooRealVar("BifurGauss_sigmaR","BifurGauss_sigmaR",rms,0.0,2*rms);
        RooRealVar * BifurGauss_sigmaL = new RooRealVar("BifurGauss_sigmaL","BifurGauss_sigmaL",rms,0.0,2*rms);


        // --- Build Bifurcated Gaussian PDF ---
        RooBifurGauss * Eraw_BifurGauss = new RooBifurGauss("Eraw_BifurGauss","mmg_ik_MZ_Surface_BifurGauss",*Mmumu,*BifurGauss_mean,*BifurGauss_sigmaR,*BifurGauss_sigmaL);


        int fewBins = 1;
        Double_t Chi2J;
	Erawframe->Clear();
	//Data_subset = (RooDataSet*)Data->reduce(*ntplVars, temp);
        //Eraw_BifurGauss = new RooBifurGauss("Eraw_BifurGauss","mmg_ik_MZ_Surface_BifurGauss",*mmg_ik_MZ_Surface,*BifurGauss_mean,*BifurGauss_sigmaR,*BifurGauss_sigmaL);

        Erawframe = Mmumu->frame();
        //Data_subset->plotOn(Erawframe,Name("myhist"),Binning(RightBinning));
        Data_subset->plotOn(Erawframe,Name("myhist"));
	Eraw_BifurGauss->fitTo(*Data_subset, Range(RangeMin, RangeMax));
        Eraw_BifurGauss->plotOn(Erawframe,Name("mycurve"));
        Erawframe->Draw();
	Chi2J = chiSquare(Erawframe,(char *)"mycurve",(char *)"myhist",3, JanChi2, DegreesOfFreedom, pValue, &fewBins);

	f = Eraw_BifurGauss->asTF( RooArgList(*Mmumu) );
        *mean_value  = BifurGauss_mean->getVal();
        *mean_error  = BifurGauss_mean->getError();
        *sigmaEff_value = effSigma(hh);
        //*sigma_value = sigwidth.getVal();
        *sigmaR_value = SigmaR(f, 0.0, 2.0);
        *sigmaL_value = SigmaL(f, 0.0, 2.0);
        *ChiSquare = Erawframe->chiSquare();


	double entries = hh->GetEntries();

        int entriesInt = (int) entries;

        double fxmax = f->GetMaximumX(0.8,1.2,1.E-10,100,false);
        //cout<<"////////// ---- FXMAX BG = "<<f->GetMaximumX(0.8,1.2,1.E-10,100,false)<<" ---- //////////"<<endl;

        TLatex latexLabel;
        latexLabel.SetTextSize(0.026);
        latexLabel.SetNDC();
        latexLabel.DrawLatex(0.16, 0.90, "CMS Preliminary 2011");
        if(isMC == 1) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV");
        if(isMC == 0) latexLabel.DrawLatex(0.16, 0.85, "#sqrt{s} = 7 TeV, #int L = 4,66892 fm^{-1}");
        //if(EndCaps == 0) latexLabel.DrawLatex(0.16, 0.85, "Barrel");
        //if(EndCaps == 1) latexLabel.DrawLatex(0.16, 0.85, "EndCaps");
        string * tempLegChain(0);
        tempLegChain = new string;
        *tempLegChain = DoubleToString(MinVar1) + " < MuonL_Pt  < "  + DoubleToString(MaxVar1);
	string * tempLegChain2(0);
	tempLegChain2 = new string;
	*tempLegChain2 = DoubleToString(MinVar2) + " < MuonS_Pt  < "  + DoubleToString(MaxVar2);
        latexLabel.DrawLatex(0.16, 0.80, tempLegChain->c_str());
	latexLabel.DrawLatex(0.16, 0.75, tempLegChain2->c_str());

        latexLabel.DrawLatex(0.16, 0.70, Form("Entries = %d",entriesInt));
        cout<<endl<<"entries = "<<entries<<endl;
        cout<<endl<<"Mean = "<<hh->GetMean()<<endl;
        latexLabel.DrawLatex(0.16, 0.65, Form("Mean = %f +/- %f",hh->GetMean(), hh->GetMeanError()));
        latexLabel.DrawLatex(0.16, 0.60, Form("RMS = %f +/- %f",hh->GetRMS(), hh->GetRMSError()));
        latexLabel.DrawLatex(0.61, 0.90, Form("#color[4]{mean = %f +/- %f}",BifurGauss_mean->getVal(), BifurGauss_mean->getError()));
        latexLabel.DrawLatex(0.61, 0.85, Form("#color[4]{#sigma_{R} = %f +/- %f}",BifurGauss_sigmaR->getVal(), BifurGauss_sigmaR->getError()));
        latexLabel.DrawLatex(0.61, 0.80, Form("#color[4]{#sigma_{L} = %f +/- %f}",BifurGauss_sigmaL->getVal(), BifurGauss_sigmaL->getError()));
        //latexLabel.DrawLatex(0.61, 0.75, Form("#color[4]{#chi^{2} = %f}",Erawframe->chiSquare()));
        latexLabel.DrawLatex(0.61, 0.75, Form("#color[4]{Default #chi^{2}/ndf = %f}",Erawframe->chiSquare()));
        latexLabel.DrawLatex(0.61, 0.70, Form("#color[4]{Jan #chi^{2}/ndf = %f}",Chi2J));
        latexLabel.DrawLatex(0.61, 0.65, Form("#color[4]{p-value = %0.6e}",*pValue));

        enregistrementPlots(nomDossier, nomFichier, 2, j, c2);

        c2->Clear();
	

	
	delete tempLegChain; tempLegChain = 0;
	delete tempLegChain2; tempLegChain2 = 0;
	delete MuonL_Pt; MuonL_Pt = 0; 
	delete MuonS_Pt; MuonS_Pt = 0;
	delete Mmumu; Mmumu = 0;
	delete weight_pileUp; weight_pileUp = 0;
	delete ntplVars; ntplVars = 0;
	delete Data; Data = 0;
	delete Data_subset; Data_subset = 0;
	delete BifurGauss_mean; BifurGauss_mean = 0;
	delete BifurGauss_sigmaR; BifurGauss_sigmaR = 0;
	delete BifurGauss_sigmaL; BifurGauss_sigmaL = 0;
	delete Eraw_BifurGauss; Eraw_BifurGauss = 0;

	


	return ;

}

//./SurfaceGeneratorFit.exe MC NoMuonCor Fall11 Mmumu.txt 5GeV RooCrystalBall           
//./SurfaceGeneratorFit.exe MC NoMuonCor Fall11 Mmumu.txt 5GeV RooBifurcatedGauss

//int SurfaceGeneratorFit(string type = "Data", string MuonCorrection = "NoMuonCor", string date = "16Jan2012", string choix = "Mmumu.txt", string SurfaceBinning = "5GeV")
int main(int argc, char *argv[])
{

	// RAJOUTER isMM == 1 && Mmumu > 80.0 && Mmumu < 100.0"

        cout << "argc= " << argc << endl;
        for(int iarg = 0 ; iarg < argc; iarg++)
        {
                cout << "argv[" << iarg << "]= " << argv[iarg] << endl;
        }

        if( argc == 1 )
        {
                cerr << "arguments :  type, MuonCorrection, date, choix, SurfaceBinning, FitMethod, ntotjob, ijob" <<endl;
                return 1;

        }
	string type = "MC";
	string MuonCorrection = "NoMuonCor";
	string date = "16Jan2012";
	string choix = "Mmumu.txt";
	string SurfaceBinning = "5GeV";
	string FitMethod = "RooCrystalBall"; //RooBifurcatedGauss
	int ntotjob = 18;
	int ijob = 0;
	int Gen = 0;

	if( argc > 1 )
        {
                type = argv[1];
        }
	if( argc > 2 )
        {
                MuonCorrection = argv[2];
        }
	if( argc > 3 )
        {
                date = argv[3];
        }
	if( argc > 4 )
        {
                choix = argv[4];
        }
	if( argc > 5 )
        {
                SurfaceBinning = argv[5];
        }
	if( argc > 6 )
        {
                FitMethod = argv[6];
        }
        if( argc > 7 )
        {
                std::stringstream ss ( argv[7] );
                ss >> ntotjob;
        }
        if( argc > 8 )
        {
                std::stringstream ss ( argv[8] );
                ss >> ijob;
        }
	if( argc > 9 )
        {
                std::stringstream ss ( argv[9] );
                ss >> Gen;
        }	


	
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


	gROOT->Reset();
        TGaxis::SetMaxDigits(3);
        setTDRStyle();
	gStyle->SetPalette(1,0);
	//gStyle->SetCanvasDefW(600);


	TCanvas *c2 = new TCanvas("c2", "c2",0,0,600,600);

	TChain *chain = new TChain("miniTree");

	// /sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/Selection/


	if(Gen == 0)
        {
	
		if(type == "MC")
		{
			if(MuonCorrection == "NoMuonCor") chain->Add("miniTreeMuons_v08_Zmumu_partALL.root");
			//if(MuonCorrection == "NoMuonCor") chain->Add("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/Selection/miniTreeMuons_v08_Zmumu_partALL.root");
			if(MuonCorrection == "RochCor") chain->Add("miniTreeMuons_v09_Zmumu_rochcor_partALL.root"); 
	
		}
		if(type == "Data")
	        {
	                if(MuonCorrection == "NoMuonCor")
			{
				if(date == "16Jan2012")
				{
					chain->Add("miniTreeMuons_v10_Zmumu_Run2011A-16Jan2012-v1_partALL.root");
					chain->Add("miniTreeMuons_v10_Zmumu_Run2011B-16Jan2012-v1_partALL.root");
				}
				if(date == "30Nov2011")
	                        {
					chain->Add("miniTreeMuons_v10_Zmumu_Run2011A-30Nov2011-v1_partALL.root");
					chain->Add("miniTreeMuons_v10_Zmumu_Run2011B-30Nov2011-v1_partALL.root");
				}
				if(date == "Vgamma")
				{
					chain->Add("miniTreeMuons_v10_Zmumu_Run2011A-ZMu-PromptSkim-v4_partALL.root");
					chain->Add("miniTreeMuons_v10_Zmumu_Run2011A-ZMu-May10ReReco-v1_partALL.root");
					chain->Add("miniTreeMuons_v10_Zmumu_Run2011A-ZMu-05Aug2011-v1_V04_partALL.root");
					chain->Add("miniTreeMuons_v10_Zmumu_Run2011A-ZMu-03Oct2011-v1_partALL.root");
					chain->Add("miniTreeMuons_v10_Zmumu_Run2011B-ZMu-PromptSkim-v1_finalJson_partALL.root");
	
	
				}
	
			}
	
			if(MuonCorrection == "RochCor")
	                {
	
				if(date == "16Jan2012")
	                        {
	                                chain->Add("miniTreeMuons_v10_Zmumu_Run2011A-16Jan2012-v1_rochcor_partALL.root");
	                                chain->Add("miniTreeMuons_v10_Zmumu_Run2011B-16Jan2012-v1_rochcor_partALL.root");
	                        }
	                        if(date == "30Nov2011")
	                        {
	                                chain->Add("miniTreeMuons_v10_Zmumu_Run2011A-30Nov2011-v1_rochcor_partALL.root");
	                                chain->Add("miniTreeMuons_v10_Zmumu_Run2011B-30Nov2011-v1_rochcor_partALL.root");
	                        }
	                        if(date == "Vgamma")
	                        {
	                                chain->Add("miniTreeMuons_v10_Zmumu_Run2011A-ZMu-PromptSkim-v4_rochcor_partALL.root");
	                                chain->Add("miniTreeMuons_v10_Zmumu_Run2011A-ZMu-May10ReReco-v1_rochcor_partALL.root");
	                                chain->Add("miniTreeMuons_v10_Zmumu_Run2011A-ZMu-05Aug2011-v1_V04_rochcor_partALL.root");
	                                chain->Add("miniTreeMuons_v10_Zmumu_Run2011A-ZMu-03Oct2011-v1_rochcor_partALL.root");
	                                chain->Add("miniTreeMuons_v10_Zmumu_Run2011B-ZMu-PromptSkim-v1_finalJson_rochcor_partALL.root");
	
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
                //cout<<endl<<"coucou"<<endl;
                }

        }


	cout<<endl<<"chain->GetEntries() = "<<chain->GetEntries()<<endl;
/*
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

*/
	////////// Allow jobs to do the same number of fits //////////

	int nbOfBins = 0;
        if(SurfaceBinning == "05GeV") nbOfBins = 400;
        if(SurfaceBinning == "1GeV") nbOfBins = 200;
        if(SurfaceBinning == "2GeV") nbOfBins = 100;
        if(SurfaceBinning == "5GeV") nbOfBins = 40;

	double BinSize = 200.0 / nbOfBins;
	
	cout<<endl<<"BinSize = "<<BinSize<<endl;

	double totalNbOfFits = 0;
	int iterMax = nbOfBins - 10 / BinSize; 
        for(int i = 1; i <= iterMax; i++) 
        {
                totalNbOfFits += i;

        }
        totalNbOfFits -= ( (10 / BinSize) * (10 / BinSize) / 2 + (10 / BinSize) / 2);	

	double nbOfFitsPerJob = totalNbOfFits / ntotjob;

	cout<<endl<<"totalNbOfFits = "<<totalNbOfFits<<", nbOfFitsPerJob = "<<nbOfFitsPerJob<<endl;

	double nbOfPerformedFits = nbOfFitsPerJob * ijob; 
	double nbOfPerformedFitsAfterJob = nbOfFitsPerJob * (ijob + 1);
	
	cout<<endl<<"nbOfPerformedFits = "<<nbOfPerformedFits<<", nbOfPerformedFitsAfterJob = "<<nbOfPerformedFitsAfterJob<<endl;

	double N0 = 20.0 / BinSize * nbOfBins + 10.0 / BinSize + 1.0;
	double alpha = 0;
	double beta = 0;
	double somme1 = 0;
	double somme2 = 0;
	int iniAlpha = 0;
	if(SurfaceBinning == "2GeV") iniAlpha = 6;
	if(SurfaceBinning == "5GeV") iniAlpha = 3;
	for(int k = iniAlpha; k <= totalNbOfFits; k++)
	{
		somme1 += k;
		if(nbOfPerformedFits <= somme1) 
		{
			alpha = k - iniAlpha;
			break;			
		}

	}

	/*
	somme2 = alpha * 2;
	for(int k = 0; k <= alpha; k++)
	{
		somme2 += k;
	}	
	beta = nbOfPerformedFits - somme2;
*/

	somme2 = 0; 
        int newIter = iniAlpha - 1;     
        for(int k = 0; k <= alpha; k++) 
        {
                newIter++;
                somme2 += newIter;
                //cout<<endl<<"k = "<<k<<", newIter = "<<newIter<<", somme2 = "<<somme2<<endl;
        }
        beta = newIter - (somme2 - nbOfPerformedFits);



	double iXj = N0 + alpha * nbOfBins + beta;

	cout<<endl<<"alpha = "<<alpha<<", beta = "<<beta<<endl;
	cout<<endl<<"iXj = "<<iXj<<endl;
	
	int initialI = 1;
	int initialJ = 1;
		
	double somme3 = 0;
	int iXjInt = (int) iXj;
	if(ijob != 0)
	{	
        	if((iXjInt % nbOfBins) == 0) 
        	{
                	initialI = iXjInt / nbOfBins;
                	initialJ = nbOfBins;
        	}
       		else
        	{
      
                	initialI = iXjInt / nbOfBins + 1; 
                	initialJ = iXjInt % nbOfBins;

        	}

		/*
		for(int k = 1; k < nbOfBins; k++)
		{
			//somme3 = k * nbOfBins;
			if(somme3 >= iXj)
			{
				initialI = k - 1; 
				initialJ = iXj - (initialI-1) * nbOfBins;
				break;
			}
			somme3 = k * nbOfBins;

		}
		*/
	}

	alpha = 0;
	beta = 0;
	somme1 = 0;
	for(int k = iniAlpha; k <= totalNbOfFits; k++)
        {
        	//cout<<endl<<"k = "<<k<<endl;
	        somme1 += k;
                if(nbOfPerformedFitsAfterJob <= somme1)
                {
                        alpha = k - iniAlpha;
                        break;
                }

        }

	somme2 = 0;
	newIter = iniAlpha - 1;	
	for(int k = 0; k <= alpha; k++)
	{
		newIter++;
		somme2 += newIter;
		//cout<<endl<<"k = "<<k<<", newIter = "<<newIter<<", somme2 = "<<somme2<<endl;
	}
	beta = newIter - (somme2 - nbOfPerformedFitsAfterJob) -1;

	/*
        somme2 = alpha * 2;
        for(int k = 0; k <= alpha; k++)
        {
                somme2 += k;
        }
        beta = nbOfPerformedFitsAfterJob - somme2 - 1;
	*/
	double iXjLast = N0 + alpha * nbOfBins + beta;

	cout<<endl<<"alpha = "<<alpha<<", beta = "<<beta<<endl;

	cout<<endl<<"iXjLast = "<<iXjLast<<endl;
	
	int lastI = 0;
	int lastJ = 0;

	somme3 = 0;
	int iXjLastInt = (int) iXjLast;
	if((iXjLastInt % nbOfBins) == 0) 
	{
		lastI = iXjLastInt / nbOfBins;
		lastJ = nbOfBins;
	}
	else
	{
		
		lastI = iXjLastInt / nbOfBins + 1;
		lastJ = iXjLastInt % nbOfBins;

	}

	if(ntotjob - ijob == 1) 
	{
		lastI = nbOfBins;
		lastJ = nbOfBins;
	}

/*
	for(int k = 1; k < nbOfBins; k++)
        {
        	//somme3 = k * nbOfBins;
                if(somme3 >= iXjLast)
                {
                	lastI = k - 1;
                        lastJ = iXjLast - (lastI-1) * nbOfBins;
                        break;
                }
		somme3 = k * nbOfBins;
        }
*/

		


	cout<<endl<<"initialI = "<<initialI<<", initialJ = "<<initialJ<<endl<<"lastI = "<<lastI<<", lastJ = "<<lastJ<<endl;

	int initialJTemp = 0;
	int lastJTemp = 0;


	////////// //////////
	string nomDossier = "";

	if(Gen == 0) nomDossier = Form("SurfacePlots_v2/%s/%s/%s/Surface_%s/%s/",type.c_str(),MuonCorrection.c_str(),date.c_str(),SurfaceBinning.c_str(),FitMethod.c_str());
	if(Gen == 1) nomDossier = Form("SurfacePlots_MMGgen_v1/%s/%s/%s/Surface_%s/Fitted_Surface/%s/",type.c_str(),MuonCorrection.c_str(),date.c_str(),SurfaceBinning.c_str(),FitMethod.c_str()); 
	string nomDossierFits = nomDossier + "Fits/";
	string nomDossierTXT = nomDossier + "TXT/";

	system(Form("mkdir -p %s", nomDossierFits.c_str()));
	system(Form("mkdir -p %s", nomDossierTXT.c_str()));

	TString temp = "";

	TH1D *MmumuBin = new TH1D("MmumuBin", "MmumuBin", 200, 80, 100);

	double infLead, infTrail, supLead, supTrail;
	double mean_value = 0;
	double mean_error = 0;
	double sigmaEff_value = 0;
	double sigma_value_error =0;
	double sigma_value = 0;
	double sigmaR_value = 0;
	double sigmaL_value = 0;
	double ChiSquare = 0;
	double mean = 0;
	double rms = 0;
        TF1 * f = new TF1();
        RooPlot* Erawframe = new RooPlot(80.0,100.0);
        string nomFichier = "Mmumu";	
	double MinRange = 86.0; //change
	double MaxRange = 95.0; //change
	double JanChi2 = 0;
	double DegreesOfFreedom = 0;
	double pValue = 0;
	double MinVar = 0;
	double MaxVar = 0;
	//double MeanTab[2000];
	vector<double> MeanTab;
	string variableX = "Mmumu";
	int isMC = 0;
	int iter = 0; 
	if(type == "Data") isMC = 0;
	if(type == "MC") isMC = 1;
	if(FitMethod == "RooBifurcatedGauss") nomFichier += "BFG";
	if(FitMethod == "RooCrystalBall") nomFichier += "CB";	

	for(int i = initialI; i <= lastI; i++)
	{
		if( i == initialI )
		{
			initialJTemp = initialJ;
		}
		else 
		{
			initialJTemp = 1;
		}
	
		if( i == lastI )
		{
			lastJTemp = lastJ;
		}
		else	
		{
			lastJTemp = nbOfBins;
		}
 
		for(int j = initialJTemp; j <= lastJTemp; j++)
		{
			iter = (i-1) * nbOfBins + j;
			if(Gen == 0) MmumuBin = new TH1D("MmumuBin", "MmumuBin", 200, 80, 100);
			if(Gen == 1) MmumuBin = new TH1D("MmumuBin", "MmumuBin", 200, 70, 110);
			infLead = (i-1)*200.0/nbOfBins;
			supLead = i*200.0/nbOfBins;
			infTrail = (j-1)*200.0/nbOfBins;
			supTrail = j*200.0/nbOfBins;
			if(Gen == 0) temp = Form("isMM == 1 && Mmumu > 80.0 && Mmumu < 100.0 && MuonL_Pt > %f && MuonL_Pt <= %f && MuonS_Pt > %f && MuonS_Pt <= %f",infLead,supLead,infTrail,supTrail);
			if(Gen == 1) temp = Form("",infLead,supLead,infTrail,supTrail);	//???
			cout<<endl<<"temp = "<<temp;
			if(Gen == 0) chain->Draw("Mmumu>>MmumuBin",temp);
			if(Gen == 1) chain->Draw("Mmumugamma_MMG_MC>>MmumuBin",temp);	
			mean = MmumuBin->GetMean();
                        rms = MmumuBin->GetRMS();
			cout<<endl<<"MmumuBin->GetEntries() = "<<MmumuBin->GetEntries()<<endl;	
			if(MmumuBin->GetEntries() != 0)
			{
				cout<<endl<<"coucou1"<<endl;
				if(FitMethod == "RooBifurcatedGauss") 
				{
					RooBifurcatedGauss(&mean_value, &mean_error, &sigmaEff_value, &sigma_value_error, &sigma_value, &sigmaR_value, &sigmaL_value, &ChiSquare, MmumuBin, iter, temp, chain, mean, rms, c2, Erawframe, f, nomDossierFits, nomFichier, MinRange, MaxRange, &JanChi2, &DegreesOfFreedom, &pValue, infLead, supLead, infTrail, supTrail, variableX, isMC);
				}
				if(FitMethod == "RooCrystalBall") 
				{
					cout<<endl<<"coucou2"<<endl;
					RooCrystalBall(&mean_value, &mean_error, &sigmaEff_value, &sigma_value_error, &sigma_value, &sigmaR_value, &sigmaL_value, &ChiSquare, MmumuBin, iter, temp, chain, mean, rms, c2, Erawframe, f, nomDossierFits, nomFichier, MinRange, MaxRange, &JanChi2, &DegreesOfFreedom, &pValue, infLead, supLead, infTrail, supTrail, variableX, isMC);
				} 
				cout<<endl<<"je fitte";
			}
			else mean_value = 91.187;
			//MeanTab[i*j] = mean_value;
			MeanTab.push_back(mean_value);
			cout<<endl<<"mean_value = "<<mean_value;
			//c2->Print("FitSurfaceTest.png");
			temp.Clear();	
			c2->Clear();
			delete MmumuBin;
			MmumuBin = 0;	
		}
	}

	int Size = MeanTab.size();
	cout<<endl<<"Size = "<<Size<<endl;
	
	if(FitMethod == "RooBifurcatedGauss")
	{
		ofstream monFlux(Form("%sMmumuFit_%d.txt",nomDossierTXT.c_str(),ijob), ios::trunc);
		cout<<endl<<"monFlux openned"<<endl;

       		for(int i = 0; i <Size; i++)
        	{
			/*
        		for(int j = 0; j <nbOfBins; j++)
                	{
                        
                		monFlux<<MeanTab[i*j]<<endl;
                	}
                	if(i%10 == 0) cout<<"i = "<<i<<endl; 
			*/
			 monFlux<<MeanTab[i]<<endl;	
        	}       

        	monFlux.close();

	}
	string TestTXT = Form("%sMmumuFit_%d.txt",nomDossierTXT.c_str(),ijob);
	cout<<endl<<"TestTXT = "<<TestTXT<<endl;
	if(FitMethod == "RooCrystalBall")
	{
		ofstream monFlux(Form("%sMmumuFit_%d.txt",nomDossierTXT.c_str(),ijob), ios::trunc);
		cout<<endl<<"monFlux openned"<<endl;
	
		for(int i = 0; i <Size; i++)
                {
			/*
                        for(int j = 0; j <nbOfBins; j++)
                        {
                     
                                monFlux<<MeanTab[i*j]<<endl;
                        }
                        if(i%10 == 0) cout<<"i = "<<i<<endl; 
               		*/
			monFlux<<MeanTab[i]<<endl;
		}    

                monFlux.close();	



	}
	cout<<endl<<"////////// END OF PROGRAM //////////"<<endl;	

/*

	if(Gen == 0) chain->Draw("MuonL_Pt:MuonS_Pt>>h__", "Mmumu*(MuonL_Pt < 200 && MuonS_Pt < 200)", "LEGO2");

	if(Gen == 1) chain->Draw("sqrt(MuonL_MC_Px*MuonL_MC_Px + MuonL_MC_Py*MuonL_MC_Py):sqrt(MuonS_MC_Px*MuonS_MC_Px + MuonS_MC_Py*MuonS_MC_Py)>>h__", "Mmumu_Muons_MC*(sqrt(MuonL_MC_Px*MuonL_MC_Px + MuonL_MC_Py*MuonL_MC_Py) < 200 && sqrt(MuonS_MC_Px*MuonS_MC_Px + MuonS_MC_Py*MuonS_MC_Py) < 200)", "LEGO2");

	TH2D *_h = (TH2D*)gDirectory->Get("h__");

	if(Gen == 0) chain->Draw("MuonL_Pt:MuonS_Pt>>h2__", "(MuonL_Pt < 200 && MuonS_Pt < 200)", "LEGO2");

	if(Gen == 1) chain->Draw("sqrt(MuonL_MC_Px*MuonL_MC_Px + MuonL_MC_Py*MuonL_MC_Py):sqrt(MuonS_MC_Px*MuonS_MC_Px + MuonS_MC_Py*MuonS_MC_Py)>>h2__", "(sqrt(MuonL_MC_Px*MuonL_MC_Px + MuonL_MC_Py*MuonL_MC_Py) < 200 && sqrt(MuonS_MC_Px*MuonS_MC_Px + MuonS_MC_Py*MuonS_MC_Py) < 200)", "LEGO2");

	TH2D *_h2 = (TH2D*)gDirectory->Get("h2__");

	_h->Divide(_h2); // Comme ça on divise le bin weighté avec Mmumu par lenombre d'entrées dans le bin : on obtient la moyenne des Mmumu pour cebin

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
	_h->GetZaxis()->SetRangeUser(80,110);
	

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

	string nomDossier = "SurfacePlots_v1/"; 
	string nomFichier = Form("Surface_%s_%s_%s_%s",type.c_str(),MuonCorrection.c_str(),date.c_str(),SurfaceBinning.c_str());

	if(choix == "plot3D" || choix == "plot2D")
	{	
		if(choix == "plot3D") nomFichier += "_3D";
		if(choix == "plot2D") nomFichier += "_2D";
		enregistrementPlots(nomDossier, nomFichier, 2, 10000, c2);	
	}

	//string nameMmumu = Form("SurfacePlots_v1/MmumuTXT/Mmumu_%s_%s_%s_%s.txt",type.c_str(),MuonCorrection.c_str(),date.c_str(),SurfaceBinning.c_str());
	string nameMmumu = Form("Mmumu_%s_%s_%s_%s.txt",type.c_str(),MuonCorrection.c_str(),date.c_str(),SurfaceBinning.c_str());
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

	if(choix == "Mmumu.txt") system("mkdir SurfacePlots_v1/MmumuTXT/");
	if(choix == "Mmumu.txt") system(Form("mv %s SurfacePlots_v1/MmumuTXT/%s", nameMmumu.c_str(),nameMmumu.c_str()));
*/	
	return 0;

}	





















