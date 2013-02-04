#include "fitFunctions.h"


void voigtian(RooDataSet *dataset, RooDataSet *dataset2, RooRealVar &variable, RooPlot *fitFrame, RooBinning b, double rangeMin, double rangeMax, vector <double> &fitParameters)
{
        RooRealVar mean("mean","mean",0.0,-0.1,0.1);
        RooRealVar sigma("sigma","sigma",0.5,0.0,1.0);
        RooRealVar width("width","width",0.5,0.0,1.0);

        RooVoigtian * pdf = new RooVoigtian("pdf","Voigtian",variable,mean,sigma,width);

        dataset->plotOn(fitFrame,Name("myhist"),Binning(b),DataError(RooAbsData::SumW2));

        RooFitResult * res = pdf->fitTo(*dataset, Range(rangeMin, rangeMax),Save(),SumW2Error(kTRUE));
	res->Print();

	//minNll = res->minNll();
        
	pdf->plotOn(fitFrame,Name("mycurve"));


	fitParameters.push_back(3); // nb of fit parameters

	fitParameters.push_back(mean.getVal());
	fitParameters.push_back(mean.getError());

	fitParameters.push_back(sigma.getVal());
        fitParameters.push_back(sigma.getError());
	
	fitParameters.push_back(width.getVal());
        fitParameters.push_back(width.getError());		
	

}

