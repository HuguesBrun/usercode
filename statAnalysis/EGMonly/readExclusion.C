

using namespace RooFit;
using namespace RooStats;


readExclusion(){
	
	TFile * file = new TFile("theExclusion_M115_EGMLOOSE006.root"); 
 HypoTestInverterResult * r = dynamic_cast<HypoTestInverterResult*>( file->Get("result_r") );
	double upperLimit = r->UpperLimit();
	double ulError = r->UpperLimitEstimatedError();
	
	
	std::cout << "The computed upper limit is: " << upperLimit << " +/- " << ulError << std::endl;
	
	const int nEntries = r->ArraySize();
	
	
	TCanvas *c0 = new TCanvas("c0","coucou",600,600);
	c0->SetFillColor(0);
	TString plotTitle = "CL Scan for workspace";
	HypoTestInverterPlot *plot = new HypoTestInverterPlot("HTI_Result_Plot",plotTitle,r);
	plot->Draw("CLb 2CL");  // plot all and Clb
	c0->Print("limit_plot.gif");
	
	const int nEntries = r->ArraySize();
	
	cout << "N entries = " << nEntries << endl;

	TCanvas * c2 = new TCanvas();
	c2->Divide( 2, TMath::Ceil(nEntries/2));
	for (int i=0; i<nEntries; i++) {
		c2->cd(i+1);
		SamplingDistPlot * pl = plot->MakeTestStatPlot(i);
		pl->SetLogYaxis(true);
		pl->Draw();
	}
	
}