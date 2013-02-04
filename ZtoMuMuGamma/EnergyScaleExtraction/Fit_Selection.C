//Script by Louis Sgandurra (December 2012)
//Select the best fit for all categories and different fit ranges.

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <iomanip>
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


int Fit_Selection(int EndCaps = 0, int r9sup = 2, string scale = "mmg_s", string SetOfCorrections = "Regression", string nomFitMethode = "RooVoigtian2", string PileUpVersion = "2011", string LowMmumuLim = "40", string HightMmumuLim = "80", bool phiCracks = true, bool etaCracks = true, string variableX = "Photon_Et", int isMC = 1, string MZbinning = "5GeV", string MuonCorrection = "Rochester", string SurfaceMethod = "Profile_Surface", string Category = "OneBin", int lowRange = 70, int highRange = 100)
{

	gROOT->Reset();
	TGaxis::SetMaxDigits(3);

	string compressScaleName = "";
        if(scale == "Photon_E_o_MC_E") compressScaleName = "ErecoOverEtrue";
        if(scale == "ComparaisonFits1overKrecoDiffMmumuJan") compressScaleName = "1overKreco";
	if(scale == "mmg_ik_MZ_Surface") compressScaleName = "1overKreco";
	if(scale == "mmg_ik") compressScaleName = "mmg_s";
	if(scale == "mmg_s") compressScaleName = "mmg_s";
	if(scale == "mmg_s_true") compressScaleName = "mmg_s_true";

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
	if(nomFitMethode == "RooLogNormal") compressFitMethodeName = "LN";
	if(nomFitMethode == "Cruijff") compressFitMethodeName = "Cruijff";	
	if(nomFitMethode == "roofit_bwcb_fft2") compressFitMethodeName = "BWxCB";

	string ChainCracks = "";
	if(phiCracks == true && etaCracks == true) ChainCracks = "WithCracks";
        if(phiCracks == true && etaCracks == false) ChainCracks = "WithoutEtaCracks";
        if(phiCracks == false && etaCracks == true) ChainCracks = "WithoutPhiCracks";
        if(phiCracks == false && etaCracks == false) ChainCracks = "WithoutCracks";

	string isMCChain = "";
        if(isMC == 0) isMCChain = "Data";
        if(isMC == 1) isMCChain = "MC";
	
	string fichierTXT = "";
	string fichierPValueTXT = "";
	string fichierMeanTXT = "";
	string fichierMeanErrorTXT = "";

	int n = 0;

	double average_pvalue = 0;
	double average_mean_error = 0;
	double temp_var = 0;
	int best_fit = 0;
	double min_error = 100;
	
	int fit_number = highRange - lowRange + 1;	

	double * meanTab = new double[fit_number];
	double * mean_errorTab = new double[fit_number];
	double * pvalue_tab = new double[fit_number];

	string directory_name = Form("Results_November_2012_v1_Rochv4_NewSelection_v3/%s/MmumuSelection%s_%s/%s/%s/%s/%s/MZbinning_%s/%s/%s/",variableX.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),ChainCracks.c_str(),SetOfCorrections.c_str(),MuonCorrection.c_str(),isMCChain.c_str(),MZbinning.c_str(),SurfaceMethod.c_str(),Category.c_str());

	cout<<endl<<"directory_name = "<<directory_name<<endl;

	for(int i = lowRange; i <= highRange; i++)
	{
		fichierTXT = directory_name + Form("%dPourcents/Fits/%s%s%s%s/",i,compressScaleName.c_str(),compressFitMethodeName.c_str(),compressR9Name.c_str(),compressEndCapsName.c_str());

		//fichierTXT = Form("Results_November_2012_v1_Rochv4_NewSelection_v3/%s/MmumuSelection%s_%s/%s/%s/%s/%s/MZbinning_%s/%s/%s/%dPourcents/Fits/",variableX.c_str(),LowMmumuLim.c_str(),HightMmumuLim.c_str(),ChainCracks.c_str(),SetOfCorrections.c_str(),MuonCorrection.c_str(),isMCChain.c_str(),MZbinning.c_str(),SurfaceMethod.c_str(),Category.c_str(),i);	


		//fichierTXT += Form("%s%s%s%s/",compressScaleName.c_str(),compressFitMethodeName.c_str(),compressR9Name.c_str(),compressEndCapsName.c_str());

		fichierPValueTXT = fichierTXT + "PValueTab.txt";
		
		fichierMeanTXT = fichierTXT + "MeanTab.txt";

		fichierMeanErrorTXT = fichierTXT + "MeanErrorTab.txt"; 

		ifstream monFlux1(fichierPValueTXT.c_str());
		
		if(monFlux1)    
		{
			//cout << "OK : Ouverture du fichier" << endl;
		}
		else
		{
    			cout << "ERREUR: Impossible d'ouvrir le fichier." << endl;
		}
	
		monFlux1 >> pvalue_tab[i-lowRange];
		average_pvalue += pvalue_tab[i-lowRange];
			
		monFlux1.close();

		ifstream monFlux2(fichierMeanTXT.c_str());

                if(monFlux2)    
                {
                        //cout << "OK : Ouverture du fichier" << endl;
                }
                else
                {
                        cout << "ERREUR: Impossible d'ouvrir le fichier." << endl;
                }
	

		monFlux2 >> meanTab[i-lowRange];

                monFlux2.close();

	
		ifstream monFlux3(fichierMeanErrorTXT.c_str());

                if(monFlux3)    
                {
                        //cout << "OK : Ouverture du fichier" << endl;
                }
                else
                {
                        cout << "ERREUR: Impossible d'ouvrir le fichier." << endl;
                }
     

                monFlux3 >> mean_errorTab[i-lowRange];
                average_mean_error += mean_errorTab[i-lowRange];

                monFlux3.close();	

	}

	average_pvalue = average_pvalue / fit_number;
	average_mean_error = average_mean_error / fit_number;
		
	min_error = average_mean_error;

	///// Fit range selection /////
/*
	for(int i = lowRange; i <= highRange; i++)
	{
		if(pvalue_tab[i-lowRange] > average_pvalue && average_pvalue < 0.5)
		{
			if(mean_errorTab[i-lowRange] < min_error && mean_errorTab[i-lowRange] < average_mean_error)	
			{
				min_error = mean_errorTab[i-lowRange];
				best_fit = i;
			}
		}	
		else if(pvalue_tab[i-lowRange] < average_pvalue &&  pvalue_tab[i-lowRange] > 0.1 && average_pvalue > 0.5)
		{
			if(mean_errorTab[i-lowRange] < min_error && mean_errorTab[i-lowRange] < average_mean_error) 
                        {
                                min_error = mean_errorTab[i-lowRange];
                                best_fit = i;
                        }	
		}

		if(best_fit == 0)
		{
			if(pvalue_tab[i-lowRange] > 0.1 && mean_errorTab[i-lowRange] < min_error)
			{
                                min_error = mean_errorTab[i-lowRange];
                                best_fit = i;
                        }

		}

	}
*/
	for(int i = highRange; i >= lowRange; i--)
	{
		//if(i == 100 || i == 89 || i = 79 || i = 69) min_error = average_mean_error;

		if(pvalue_tab[i-lowRange] > 0.01)
		{
			if(mean_errorTab[i-lowRange] < min_error && i >= 90)
			{
				min_error = mean_errorTab[i-lowRange];
                                best_fit = i;
			}
			if(i >= 80 && i < 90 && best_fit == 0)
			{
                                min_error = mean_errorTab[i-lowRange];
                                best_fit = i;
                        }
			if(i >= 70 && i < 80 && best_fit == 0)
                        {
                                min_error = mean_errorTab[i-lowRange];
                                best_fit = i;
                        }
			if(i >= 60 && i < 70 && best_fit == 0)
                        {
                                min_error = mean_errorTab[i-lowRange];
                                best_fit = i;
                        }
		}
	}
	if(best_fit == 0) 
	{
		for(int i = highRange; i >= lowRange; i--)
        	{
			if(pvalue_tab[i-lowRange] > 0.001)
	                {
	                        if(mean_errorTab[i-lowRange] < min_error && i >= 90) 
	                        {
	                                min_error = mean_errorTab[i-lowRange];
	                                best_fit = i;
	                        }
	                        if(i >= 80 && i < 90 && best_fit == 0)
	                        {
	                                min_error = mean_errorTab[i-lowRange];
	                                best_fit = i;
	                        }
	                        if(i >= 70 && i < 80 && best_fit == 0)
	                        {
	                                min_error = mean_errorTab[i-lowRange];
	                                best_fit = i;
	                        }
	                        if(i >= 60 && i < 70 && best_fit == 0)
	                        {
	                                min_error = mean_errorTab[i-lowRange];
	                                best_fit = i;
	                        }
	                }
	
		}
	}

	if(best_fit == 0) 
	{
	
		for(int i = highRange; i >= lowRange; i--)
                {
                        if(pvalue_tab[i-lowRange] > average_pvalue)
                        {
                                if(mean_errorTab[i-lowRange] < min_error && i >= 90) 
                                {
                                        min_error = mean_errorTab[i-lowRange];
                                        best_fit = i;
                                }
                                if(i >= 80 && i < 90 && best_fit == 0)
                                {
                                        min_error = mean_errorTab[i-lowRange];
                                        best_fit = i;
                                }
                                if(i >= 70 && i < 80 && best_fit == 0)
                                {
                                        min_error = mean_errorTab[i-lowRange];
                                        best_fit = i;
                                }
                                if(i >= 60 && i < 70 && best_fit == 0)
                                {
                                        min_error = mean_errorTab[i-lowRange];
                                        best_fit = i;
                                }
                        }
     
                }	
	}

	if(best_fit == 0)
        {
		cout<<endl<<"Unable to find a good fit"<<endl;
                best_fit = 200;
	}

	cout<<endl<<"average_pvalue = "<<average_pvalue<<", average_mean_error = "<<average_mean_error<<endl;

	cout<<endl<<"Best fit range = "<<best_fit<<" %"<<endl;
	cout<<endl<<"mean_errorTab[best_fit] = "<<mean_errorTab[best_fit-lowRange]<<" %, pvalue_tab[best_fit] = "<<pvalue_tab[best_fit-lowRange]<<endl;

	///// Systematics /////

	double systematics = 0;
	double temp_error = 0;

	for(int i = lowRange; i <= highRange; i++)
	{
		//if(mean_errorTab[i-lowRange] > (10 * mean_errorTab[best_fit-lowRange]) || pvalue_tab[i-lowRange] < 0.01)
		if(mean_errorTab[i-lowRange] > fabs(10 * mean_errorTab[best_fit-lowRange]) || pvalue_tab[i-lowRange] < 0.000001)
		{
			continue;
		}	
		
		if(meanTab[i-lowRange] < meanTab[best_fit-lowRange])
		{
			//temp_error = fabs( (meanTab[best_fit-lowRange] - mean_errorTab[best_fit-lowRange]) - (meanTab[i-lowRange] - mean_errorTab[i-lowRange]) );
			temp_error = fabs(meanTab[best_fit-lowRange] - meanTab[i-lowRange]);
		}
		else
		{
			//temp_error = fabs( (meanTab[i-lowRange] + mean_errorTab[i-lowRange]) - (meanTab[best_fit-lowRange] + mean_errorTab[best_fit-lowRange]));	
			temp_error = fabs(meanTab[i-lowRange] - meanTab[best_fit-lowRange]);	
		}
		if(temp_error > systematics) systematics = temp_error;
		//cout<<endl<<"temp_error = "<<temp_error;
	}

	cout<<endl<<"mean = "<<meanTab[best_fit-lowRange]<<" +- "<<mean_errorTab[best_fit-lowRange]<<" (stat) +- "<<systematics<<" (range) %"<<endl;

	system(Form("source create_directory.sh %sSelected_Fits/",directory_name.c_str()));

	system(Form("cp -r %s%dPourcents/Fits/%s%s%s%s/ %sSelected_Fits/%s%s%s%s/",directory_name.c_str(),best_fit,compressScaleName.c_str(),compressFitMethodeName.c_str(),compressR9Name.c_str(),compressEndCapsName.c_str(),directory_name.c_str(),compressScaleName.c_str(),compressFitMethodeName.c_str(),compressR9Name.c_str(),compressEndCapsName.c_str()));	


	ofstream monFlux4(Form("%sSelected_Fits/Summary_%s.txt",directory_name.c_str(),scale.c_str()), ios::app);

        monFlux4 << scale << " " << isMCChain << " >> " << compressR9Name << " " << compressEndCapsName << ", " << compressFitMethodeName <<" : mu = " << meanTab[best_fit-lowRange] << " +- " << mean_errorTab[best_fit-lowRange] << " (stat) +- " << systematics << " (range) %" << ", fit range = "<< best_fit <<"%, p-value = " << pvalue_tab[best_fit-lowRange] <<endl;

        monFlux4.close();
	
	delete [] meanTab;
	meanTab = 0;
	delete [] mean_errorTab;
	mean_errorTab = 0;
	delete [] pvalue_tab;
	pvalue_tab = 0;

	return 0;
}
