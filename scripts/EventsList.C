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
#include "TGaxis"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <iomanip>
#pragma optimize 0
using namespace std;

int EventsList(string inputMiniTree = "miniTree_partALL")
{
	TChain* plaf = new TChain("miniTree");
	plaf->Add(Form("%s.root",inputMiniTree.c_str()));
	plaf->SetScanField(0);
	plaf->Scan("iRunID:iLumiID:iEventID", "isLooseMMG == 1", "colsize=50"); 

	return 0;
}

