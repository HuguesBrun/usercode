

using namespace RooFit;
using namespace RooStats;

TString IntToString(int number){
	ostringstream oss;
	oss << number;
	return oss.str();
}


int number = 4;
int massPoint = 140; 
TString nomPlot = "CICnew";


readExclusion_M140(){
	
    TString massPointString = IntToString(massPoint);
    
	TFile * file = new TFile("theExclusion_M"+massPointString+"_"+nomPlot+".root"); 
 HypoTestInverterResult * r = dynamic_cast<HypoTestInverterResult*>( file->Get("result_r") );
	double upperLimit = r->UpperLimit();
	double ulError = r->UpperLimitEstimatedError();
	
	
	std::cout << "The computed upper limit is: " << upperLimit << " +/- " << ulError << std::endl;
	
	const int nEntries = r->ArraySize();
	
	
	TCanvas *c0 = new TCanvas("c0","coucou",600,600);
	c0->SetFillColor(0);
	TString plotTitle = "CL Scan for workspace";
	HypoTestInverterPlot *plot = new HypoTestInverterPlot("HTI_Result_Plot",plotTitle,r);
	plot->Draw("CLb 2CL");//("CLb 2CL");  // plot all and Clb
	c0->Print("limit_plot_M"+massPointString+".gif");
	
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
    
    std::cout << "The computed upper limit is: " << upperLimit << " +/- " << ulError << std::endl;
    std::cout << "The expected limi is : " << r->GetExpectedUpperLimit(0) << endl;
    
    float expected = r->GetExpectedUpperLimit(0);
    float observed = r->UpperLimit();
    float observedError = r->UpperLimitEstimatedError();
    float expected1sP = r->GetExpectedUpperLimit(1);
    float expected2sP = r->GetExpectedUpperLimit(2);
    float expected1sM = r->GetExpectedUpperLimit(-1);
    float expected2sM = r->GetExpectedUpperLimit(-2);

    
    cout << "ZZZZ graph->SetPoint(" <<number << "," << massPoint<< "," << expected   << ");" << endl;
    cout << "ZZZZ grae->SetPoint(" <<number << "," << massPoint<< "," << expected   << ");" << endl;
    cout << "ZZZZ grae->SetPointError(" << number << ",0,0," << (expected-expected1sM) << "," << (expected1sP-expected) << ");" << endl;
    cout << "ZZZZ grae2->SetPoint(" <<number << "," << massPoint<< "," << expected   << ");" << endl;
    cout << "ZZZZ grae2->SetPointError(" << number << ",0,0," << (expected-expected2sM) << "," << (expected2sP-expected) << ");" << endl;
    cout << "ZZZZ gre->SetPoint(" <<number << "," << massPoint<< "," << observed   << ");" << endl;
    cout << "ZZZZ gre->SetPointError(" <<number << ",0," << observedError <<");" << endl;	
}