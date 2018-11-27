#include <iostream>
#include <fstream>

#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TStyle.h"

//#include "external/utils.h" // for the timer


using namespace std;

extern "C"
{
  void gmuini_();
  void gmucha_();
  void gmubeg_();
  void gmugna_();
  void gmufil_();
  extern struct
  {
    double sigma, err_sigma;
    double pl[7][4];
  } output_;
  extern struct
  {
    int ilepton, proc, ngen;
  } settings_;
}

int main()
{
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.12);

  //Timer tmr;
  TLorentzVector l1, l2, ip1, ip2, op1, op2;
  TTree *t;
  const Int_t max_particles = 9;
  Int_t npart;
  Double_t xsect, err_xsect;
  Double_t px[max_particles], py[max_particles], pz[max_particles], E[max_particles];
  Double_t weight;
  Double_t x1, x2, t1, t2;
  Double_t yy;
  Float_t gen_time;

  t = new TTree("events", "Events collection");
  t->Branch("xsect", &xsect, "xsect/D");
  t->Branch("errxsect", &err_xsect, "errxsect/D");
  t->Branch("npart", &npart, "npart/I");
  t->Branch("gen_time", &gen_time, "gen_time/F");
  t->Branch("px", px, "px[npart]/D");
  t->Branch("py", py, "py[npart]/D");
  t->Branch("pz", pz, "pz[npart]/D");
  t->Branch("E", E, "E[npart]/D");

  typedef enum
  {
    invm, ptpair, ptpairzoom,
    pt1, pt2, y1, y2,
    mx, my,
    logx1, logx2,
    num_dists
  } dists;

  const Int_t num_distributions = num_dists;
  const Int_t num_2d_distributions = 7; //FIXME

  TCanvas *c[num_distributions], *c2[2];
  TH1D *h[num_distributions];
  TString unit[num_distributions];
  TH2D *h2[num_2d_distributions];

  gmuini_(); // initializes parameters
  gmucha_(); // loads the input card
  gmubeg_(); // computes cross-section and prepares the events generation's grid

  cout << "Computed cross-section = " << output_.sigma*1.e3 << " +/- " << output_.err_sigma*1.e3 << " pb" << endl;

  h[invm] = new TH1D("invm", "m(pair)", 200, 0., 500.);       unit[invm] = "GeV";
  h[ptpair] = new TH1D("ptpair", "p_{T}(pair)", 100, 0., 50.);  unit[ptpair] = "GeV";
  h[ptpairzoom] = new TH1D("ptpair_zoom", "p_{T}(pair)", 100, 0., 5.);  unit[ptpairzoom] = "GeV";
  h[y1] = new TH1D("y1", "y(lepton 1)", 100, -4., 4.);        unit[y1] = "";
  h[y2] = new TH1D("y2", "y(lepton 2)", 100, -4., 4.);        unit[y2] = "";
  h[pt1] = new TH1D("pt1", "p_{T}(lepton 1)", 150, 0., 150.); unit[pt1] = "GeV";
  h[pt2] = new TH1D("pt2", "p_{T}(lepton 2)", 150, 0., 150.); unit[pt2] = "GeV";
  h[mx] = new TH1D("mx", "M_{X}", 300, 0., 1000.);             unit[mx] = "GeV";
  h[my] = new TH1D("my", "M_{Y}", 300, 0., 1000.);             unit[my] = "GeV";
  h[logx1] = new TH1D("logx1", "-log_{10}(x_{1})", 100, 0., 5.); unit[logx1] = "";
  h[logx2] = new TH1D("logx2", "-log_{10}(x_{2})", 100, 0., 5.); unit[logx2] = "";
  
  h2[0] = new TH2D("pt1_vs_pt2", "", 50, 0., 150., 50, 0., 150.);
  h2[1] = new TH2D("y1_vs_y2", "", 50, -2.5, 2.5, 50, -2.5, 2.5);
  h2[2] = new TH2D("x1_vs_x2", "", 100, 0., 0.1, 100, 0., 0.1);
//  h2[3] = new TH2D("logx1_vs_logx2", "", 100, 0., 5., 100, 0., 5.);
  h2[3] = new TH2D("MX_vs_MY", "", 150, 0., 1000., 150, 0., 1000.);
  h2[4] = new TH2D("t1_vs_t2", "", 100, 0., 0.5, 100, 0., 0.5);
//  h2[5] = new TH2D("logt1_vs_logt2", "", 100, 0., 10., 100, 0., 10.);
  h2[5] = new TH2D("logt1_vs_logt2", "", 100, -3., 3., 100, -3., 3.);
  h2[6] = new TH2D("ptpair_vs_dphi", "", 25, 0., 50., 36, 0., 180.);

  weight = 1./settings_.ngen*output_.sigma;

  for (Int_t i=0; i<settings_.ngen; i++) {
    //tmr.reset();

    gmugna_(); // generates the event
    gmufil_(); // fills common block with the event's kinematics

    //gen_time = tmr.elapsed();
    
    if (i%10000==0) cout << "Generating event #" << i << endl;

    // first lepton
    l1.SetXYZT(output_.pl[5][1], output_.pl[5][2], output_.pl[5][3], output_.pl[5][0]);
    // second lepton
    l2.SetXYZT(output_.pl[6][1], output_.pl[6][2], output_.pl[6][3], output_.pl[6][0]);
    // first incoming proton
    ip1.SetXYZT(output_.pl[0][1], output_.pl[0][2], output_.pl[0][3], output_.pl[0][0]);
    // second incoming proton
    ip2.SetXYZT(output_.pl[1][1], output_.pl[1][2], output_.pl[1][3], output_.pl[1][0]);
    // first outgoing proton (or remnant)
    op1.SetXYZT(output_.pl[2][1], output_.pl[2][2], output_.pl[2][3], output_.pl[2][0]);
    // second outgoing proton (or remnant)
    op2.SetXYZT(output_.pl[4][1], output_.pl[4][2], output_.pl[4][3], output_.pl[4][0]);

    h[invm]->Fill((l1+l2).M(), weight);
    h[ptpair]->Fill((l1+l2).Pt(), weight);
    h[ptpairzoom]->Fill((l1+l2).Pt(), weight);
    h[y1]->Fill(l1.Rapidity(), weight);
    h[y2]->Fill(l2.Rapidity(), weight);
    h[pt1]->Fill(l1.Pt(), weight);
    h[pt2]->Fill(l2.Pt(), weight);
    h[mx]->Fill(op1.M(), weight);
    h[my]->Fill(op2.M(), weight);

    h2[0]->Fill(l1.Pt(), l2.Pt(), weight);
    h2[1]->Fill(l1.Rapidity(), l2.Rapidity(), weight);

    x1 = 1.-(op1.E()+fabs(op1.Pz()))/(ip1.E()+fabs(ip1.Pz()));
    x2 = 1.-(op2.E()+fabs(op2.Pz()))/(ip2.E()+fabs(ip2.Pz()));
    t1 = -(op1-ip1).M2();
    t2 = -(op2-ip2).M2();

    h[logx1]->Fill(-log10(x1), weight);
    h[logx2]->Fill(-log10(x2), weight);

    //cout << "x1=" << x1 << ", EX1=" << op1.E() << ", PXz=" << op1.Pz() << ", E1=" << ip1.E() << ", P1z=" << ip1.Pz() << endl;
    //cout << "x2=" << x2 << ", EX2=" << op2.E() << ", PXz=" << op2.Pz() << ", E2=" << ip2.E() << ", P2z=" << ip2.Pz() << endl;

    h2[2]->Fill(x1, x2, weight);
//    h2[3]->Fill(-log10(x1), -log10(x2), weight);
    h2[3]->Fill(op1.M(), op2.M(), weight);
    h2[4]->Fill(t1, t2, weight);
//    h2[5]->Fill(-log10(t1), -log10(t2), weight);
    h2[5]->Fill(log10(t1), log10(t2), weight);
    h2[6]->Fill((l1+l2).Pt(), fabs(l1.DeltaPhi(l2))/TMath::Pi()*180., weight);

    // We fill the tree
    xsect = output_.sigma*1.e3;
    err_xsect = output_.err_sigma*1.e3;
    npart = 0;
    for (Int_t i=0; i<7; i++) {
      px[npart] = output_.pl[i][1];
      py[npart] = output_.pl[i][2];
      pz[npart] = output_.pl[i][3];
      E[npart] = output_.pl[i][0];
      npart++;
    }
    t->Fill();
  }

  t->SaveAs("events.root");

  TString xtitle, ytitle, tmp;

//  ofstream outfile("mx_contents.dat");
//  h[mx]->Print("all"); 

    Int_t nbins = h[mx]->GetNbinsX();
		          
    for (Int_t i=1; i<=nbins; i++) {
    	printf("%g %g\n",h[mx]->GetBinLowEdge(i)+h[mx]->GetBinWidth(i)/2,h[mx]->GetBinContent(i));
    }
  

	ofstream output("output.dat");
	
		          
    for (Int_t i=1; i<=nbins; i++) {
	    yy = (h[mx]->GetBinContent(i)*h[mx]->GetBinWidth(i));
//    	output << h[mx]->GetBinLowEdge(i)+h[mx]->GetBinWidth(i)/2 << "\t" << h[mx]->GetBinContent(i) << endl;
    	output << h[mx]->GetBinLowEdge(i)+h[mx]->GetBinWidth(i)/2 << "\t" << yy << endl;
    }
  


  for (Int_t i=0; i<num_distributions; i++) {
    c[i] = new TCanvas;
    h[i]->Draw();
    xtitle = h[i]->GetTitle();
    tmp.Form("%3.1f", (h[i]->GetXaxis()->GetXmax()-h[i]->GetXaxis()->GetXmin())/h[i]->GetXaxis()->GetNbins());
    ytitle = TString("d#sigma/d")+h[i]->GetTitle()+TString(" (nb / ")+tmp;
    if (unit[i]!="") {
      xtitle += TString(" (")+unit[i]+TString(")");
      ytitle += TString(" ")+unit[i];
    }
    ytitle += TString(")");
    h[i]->GetXaxis()->SetTitle(xtitle);
    h[i]->GetYaxis()->SetTitle(ytitle);
    h[i]->SetTitle("");
    c[i]->SaveAs("output/"+TString(h[i]->GetName())+".png");
    c[i]->SetLogy();
    //c[i]->SetLogx();
    c[i]->SaveAs("output/"+TString(h[i]->GetName())+"_logscale.png");
    
  }

  // 2-dimensional distribution(s)
  TString x2title[num_2d_distributions], y2title[num_2d_distributions];
  x2title[0] = "p_{T}(lepton 1) (GeV)"; y2title[0] = "p_{T}(lepton 2) (GeV)";
  x2title[1] = "y(lepton 1)"; y2title[1] = "y(lepton 2)";
  x2title[2] = "x_{1}"; y2title[2] = "x_{2}";
//  x2title[3] = "-log_{10}(x_{1})"; y2title[3] = "-log_{10}(x_{2})";
  x2title[3] = "MX"; y2title[3] = "MY";
  x2title[4] = "-t_{1}"; y2title[4] = "-t_{2}";
//  x2title[5] = "-log_{10}(-t_{1})"; y2title[5] = "-log_{10}(-t_{2})";
  x2title[5] = "log_{10}(-t_{1})"; y2title[5] = "log_{10}(-t_{2})";

  for (Int_t i=0; i<num_2d_distributions; i++) {
    c2[i] = new TCanvas(h2[i]->GetName(), "", 600, 600);
    h2[i]->Draw("colz");
    h2[i]->GetXaxis()->SetTitle(x2title[i]);
    h2[i]->GetYaxis()->SetTitle(y2title[i]);
    c2[i]->SaveAs("output/"+TString(h2[i]->GetName())+"_colz.png");
    c2[i]->SetLogz();
    c2[i]->SaveAs("output/"+TString(h2[i]->GetName())+"_colz_logscale.png");
    c2[i]->SetLogz(0);
    h2[i]->Draw("lego");
    c2[i]->SetPhi(-120);
    c2[i]->SaveAs("output/"+TString(h2[i]->GetName())+"_lego.png");
    c2[i]->SetLogz();
    c2[i]->SaveAs("output/"+TString(h2[i]->GetName())+"_lego_logscale.png");
    h2[i]->Draw("colz");
    /*c2[i]->SetLogx();
    c2[i]->SetLogy();
    c2[i]->SaveAs(TString(h2[i]->GetName())+"_colz_logscalexy.png");
    c2[i]->SetLogz();
    c2[i]->SaveAs(TString(h2[i]->GetName())+"_colz_logscalexyz.png");
    c2[i]->SetLogz(0);
    h2[i]->Draw("lego");
    c2[i]->SaveAs(TString(h2[i]->GetName())+"_lego_logscalexy.png");
    c2[i]->SetLogz();
    c2[i]->SaveAs(TString(h2[i]->GetName())+"_lego_logscalexyz.png");*/
  }

  return 0;
}
