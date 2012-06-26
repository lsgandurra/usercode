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

bool is_readable( const std::string & file );

bool is_readable( const std::string & file ) 
{ 
    std::ifstream fichier( file.c_str() ); 
    return !fichier.fail(); 
} 

//./SurfaceGeneratorFit.exe MC NoMuonCor Fall11 Mmumu.txt 5GeV RooCrystalBall           
//./SurfaceGeneratorFit.exe MC NoMuonCor Fall11 Mmumu.txt 5GeV RooBifurcatedGauss

//int SurfaceGeneratorFit(string type = "Data", string RochesterCorrections = "NoMuonCor", string date = "16Jan2012", string choix = "Mmumu.txt", string SurfaceBinning = "5GeV")
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
                cerr << "arguments :  type, RochesterCorrections, date, choix, SurfaceBinning, FitMethod, ntotjob, ijob" <<endl;
                return 1;

        }
	string type = "MC";
	string RochesterCorrections = "NoMuonCor";
	string date = "Vgamma";
	string choix = "Mmumu.txt";
	string SurfaceBinning = "5GeV";
	string FitMethod = "RooCrystalBall"; //RooBifurcatedGauss
	int ntotjob = 18;
	int ijob = 17;

	if( argc > 1 )
        {
                type = argv[1];
        }
	if( argc > 2 )
        {
                RochesterCorrections = argv[2];
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


	int Gen = 0;
        //int Gen = 1;
	
	//string type = "Data"; 
	//string type = "MC";

	//string RochesterCorrections = "NoMuonCor";
	//string RochesterCorrections = "RochCor";

	//string date = "16Jan2012";
	//string date = "30Nov2011";
	//string date = "Vgamma";
	
	if(type == "MC") date = "Fall11";

	int ntotjobM1 = 0;
	if(SurfaceBinning == "5GeV") ntotjobM1 = 17;
	if(SurfaceBinning == "2GeV") ntotjobM1 = 44;

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



	string nomDossier = Form("SurfacePlots_v2/%s/%s/%s/Surface_%s/%s/",type.c_str(),RochesterCorrections.c_str(),date.c_str(),SurfaceBinning.c_str(),FitMethod.c_str());
	string nomDossierTXT = nomDossier + "TXT/";

	cout<<endl<<"nomDossierTXT = "<<nomDossierTXT<<endl; 

	string temp = "";	

	if ( is_readable( Form("%sliste.txt",nomDossierTXT.c_str()) ) ) 
        { 
                cout << "Fichier existant et lisible.\n"; 
                system(Form("rm %sliste.txt",nomDossierTXT.c_str()));
        } 
        else 
        { 
                cout << "Fichier inexistant ou non lisible.\n"; 
        }	

	system(Form("source Sort.sh %s %d",nomDossierTXT.c_str(),ntotjobM1));

	if ( is_readable( Form("%sMmumu.txt",nomDossierTXT.c_str()) ) )
        {
                cout << "Fichier existant et lisible.\n";
                system(Form("rm %sMmumu.txt",nomDossierTXT.c_str()));
        }
        else 
        {
                cout << "Fichier inexistant ou non lisible.\n"; 
        }    
	


	ifstream monFlux(Form("%sliste.txt",nomDossierTXT.c_str()));
	for(int i = 0; i <= ntotjobM1; i++)
	{
		temp.clear();
		monFlux>>temp;
		cout<<"temp = "<<temp<<endl;
		system(Form("cat %s%s >> %sMmumu.txt", nomDossierTXT.c_str(),temp.c_str(),nomDossierTXT.c_str()));

	}	
	
	monFlux.close();


	return 0;

}	





















