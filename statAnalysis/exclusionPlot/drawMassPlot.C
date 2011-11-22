drawMassPlot(){

TCanvas *c0 = new TCanvas("c0", "coucou",700,500);
c0->SetFillColor(0);
c0->SetBorderMode(0);



TGraphAsymmErrors *grae = new TGraphAsymmErrors(5);
grae->SetFillColor(kTeal+9); 




TGraphAsymmErrors *grae2 = new TGraphAsymmErrors(5);
grae2->SetTitle("");
grae2->GetXaxis()->SetTitle("M_{H} GeV/c^{2}");
grae2->SetFillColor(kYellow);


TGraph *graph = new TGraph(5);
graph->SetFillColor(1);
graph->SetLineStyle(2);
graph->SetLineWidth(2);
graph->SetLineColor(kRed);

TGraphErrors *gre = new TGraphErrors(5);
gre->SetFillColor(1);
gre->SetLineWidth(2);
gre->SetMarkerColor(kBlack);
gre->SetMarkerStyle(20);

 graph->SetPoint(0,110,4.91073);
 grae->SetPoint(0,110,4.91073);
 grae->SetPointError(0,0,0,1.2935,1.08999);
 grae2->SetPoint(0,110,4.91073);
 grae2->SetPointError(0,0,0,2.12693,2.69447);
 gre->SetPoint(0,110,3.76726);
 gre->SetPointError(0,0,0.0739304);
 graph->SetPoint(1,115,3.26991);
 grae->SetPoint(1,115,3.26991);
 grae->SetPointError(1,0,0,0.767079,1.31215);
 grae2->SetPoint(1,115,3.26991);
 grae2->SetPointError(1,0,0,1.30523,2.2696);
 gre->SetPoint(1,115,4.53826);
 gre->SetPointError(1,0,0.0316776);
 graph->SetPoint(2,120,2.85985);
 grae->SetPoint(2,120,2.85985);
 grae->SetPointError(2,0,0,0.880408,0.972118);
 grae2->SetPoint(2,120,2.85985);
 grae2->SetPointError(2,0,0,1.25489,2.17879);
 gre->SetPoint(2,120,2.94032);
 gre->SetPointError(2,0,0.0295871);
 graph->SetPoint(3,130,2.9124);
 grae->SetPoint(3,130,2.9124);
 grae->SetPointError(3,0,0,0.840572,1.00115);
 grae2->SetPoint(3,130,2.9124);
 grae2->SetPointError(3,0,0,1.16341,2.16229);
 gre->SetPoint(3,130,2.83868);
 gre->SetPointError(3,0,0.0323039);
 graph->SetPoint(4,140,3.35217);
 grae->SetPoint(4,140,3.35217);
 grae->SetPointError(4,0,0,0.963504,1.25723);
 grae2->SetPoint(4,140,3.35217);
 grae2->SetPointError(4,0,0,1.48371,2.21231);
 gre->SetPoint(4,140,3.41541);
 gre->SetPointError(4,0,0.0343649);

grae2->Draw("a3::same");
grae2->GetXaxis()->SetTitle("M_{H} GeV/c^{2}");
grae2->GetYaxis()->SetTitle("#sigma(H #rightarrow #gamma #gamma)_{95% CL}/#sigma(H #rightarrow #gamma #gamma)_{SM}");
grae2->SetMinimum(0);
grae2->SetMaximum(10);
TLine *line = new TLine(110,1,140,1);
line->SetLineWidth(3);
line->SetLineColor(kBlue);
line->Draw("same");


grae->Draw("3::same");
graph->Draw("L::same");
gre->Draw("pl::same");


TLatex latexLabel;
latexLabel.SetTextSize(0.04);
latexLabel.SetNDC();
latexLabel.DrawLatex(0.16, 0.85, "CMS Preliminary 2011");
latexLabel.DrawLatex(0.16, 0.8, "#sqrt{s} = 7 TeV, L = 1.143 fb^{-1}");

TLegend *t = new TLegend(0.52,0.66,0.84,0.86);
t->SetLineColor(0);
t->SetFillColor(0);
t->AddEntry(gre, "observed CLs limit", "Pl");
t->AddEntry(graph, "expected CLs limit","l");
t->AddEntry(grae,"#pm #sigma expected limit","F");
t->AddEntry(grae2,"#pm 2#sigma expected limit","F");
t->Draw();

c0->Print("gif/theExclusionMassPlot.gif");

}

