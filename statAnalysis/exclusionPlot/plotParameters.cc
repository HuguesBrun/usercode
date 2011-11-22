void plotParameters(RooArgSet *r2_cat0_param, TCanvas *c, int canvasDivision, RooPlot* frame0, bool isSignal)
{
	c->cd(canvasDivision);
	TLatex latexLabel;
	latexLabel.SetNDC();
  latexLabel.SetTextSize(0.03);
  latexLabel.DrawLatex(0.13, 0.91, "CMS Preliminary 2011");
  latexLabel.DrawLatex(0.45, 0.91, "#sqrt{s} = 7 TeV");
  latexLabel.DrawLatex(0.67, 0.91, isSignal ? "Simulation" : "DATA");
  latexLabel.SetTextSize(0.02);
  TIterator *it = (TIterator*) r2_cat0_param->createIterator();
  RooRealVar* obj = new RooRealVar();
  double position = 0.92;
  position -= 0.04;
  latexLabel.DrawLatex(0.55, position, isSignal ? "Fit: Sum of Gaussians" : "Fit: Second Order Bernstein Polynomial");
  while(( (RooRealVar*)obj = it->Next()) != 0)
  {
   if( ! strcmp(((char*)obj->GetName()), "mgg") ) continue;
//    cout << "obj->getTitle()= " << obj->getTitle() << endl; // char*
//    cout << "obj->GetName()= " << obj->GetName() << endl; // char*
//    cout << "obj->getVal()= " << obj->getVal() << endl; // Double_t
//    cout << "obj->getError()= " << obj->getError() << endl; // Double_t
//    cout << "obj->getUnit()= " << obj->getUnit() << endl; // Text_t
//    cout << endl;
    position -= 0.04;
    std::ostringstream valueStream;
    if( (double)obj->getError() != 0.0 )
    {
      valueStream << setprecision (2) << fixed << (double)obj->getVal() << " +- " << (double)obj->getError();
    } else {
       valueStream << setprecision (2) << fixed << (double)obj->getVal();
    }
    string valueString = valueStream.str();
    latexLabel.DrawLatex(0.60, position, Form("%s = %s %s", obj->GetTitle(), valueString.c_str(), (char*)obj->getUnit()));
  }
//  cout << "it->Next()->GetName()= " << it->Next()->GetName() << "\tit->Next()->getVal()= " << it->Next()->getVal() << endl;

  position -= 0.04;
  std::ostringstream valueStream;
  valueStream << setprecision (3) << fixed << (double)(frame0->chiSquare("model", isSignal ? "mc" : "data", isSignal ? 7 : 3));
  string valueString = valueStream.str();
  latexLabel.DrawLatex(0.60, position, Form("#chi^{2} / ndf = %s", valueString.c_str()));

	return;
}
