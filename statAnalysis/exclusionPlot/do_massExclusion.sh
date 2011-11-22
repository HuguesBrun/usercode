#!/bin/bach

rm drawMassPlot.C

cat >> drawMassPlot.C << EOF
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

EOF

root -b -q readExclusion_M110.C | grep ZZZZ | tr -d ZZZZ >> drawMassPlot.C
root -b -q readExclusion_M115.C | grep ZZZZ | tr -d ZZZZ >> drawMassPlot.C
root -b -q readExclusion_M120.C | grep ZZZZ | tr -d ZZZZ >> drawMassPlot.C
root -b -q readExclusion_M130.C | grep ZZZZ | tr -d ZZZZ >> drawMassPlot.C
root -b -q readExclusion_M140.C | grep ZZZZ | tr -d ZZZZ >> drawMassPlot.C


cat >> drawMassPlot.C << EOF

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

EOF
