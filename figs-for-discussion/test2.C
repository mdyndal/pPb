#include <fstream>
#define max_events 2500000
#define invm 0
#define ptpair 1
#define pt1 2
#define pt2 3
#define ptsqr 4
#define dpt 5
#define dphi 6
#define ptsingle 7
#define y1 8
#define y2 9
#define logx1 10
#define logx2 11
#define t1 12
#define t2 13
#define logt1 14
#define logt2 15
#define msqbal 16
#define mllsq 17
//#define mggsq 18
#define mx 18
#define s1 19
#define th1 20
#define u1 21
#define a1 22
#define b1 23
#define c1 24
#define d1 25

#define max_dists 26


void
test2() // this should match your file name (e.g. plot() if file is called plot.C... thanks a lot, ROOT...)
{
  TFile *f_pptoll[16];
  TTree *t_pptoll[16];
  TH1D *h[16][max_dists];
  TLorentzVector l1, l2, op1, op2, ip1, ip2;
  Double_t weight, x1, x2;
  TCanvas *c[max_dists];
  TString unit[max_dists];
  TLine *line;
  TString xtitle, ytitle, tmp;
  Double_t max;
  Double_t t1val, s1val, tval, u1val, mval, m1val, t2val, amw2, c1val, dval1, dval2, yy;


  gStyle->SetOptStat(0);

  f_pptoll[0] = new TFile("events.root"); // Replace me !!!
//  f_pptoll[1] = new TFile("pptoll_output2.root"); // Replace me !!!

  for (Int_t i=0; i<2; i++) {
    tmp.Form("invm_%i", i); h[i][invm] = new TH1D(tmp, "m(#mu#mu)", 200, 0., 100.);   unit[invm] = "GeV";
    tmp.Form("ptpair_%i", i); h[i][ptpair] = new TH1D(tmp, "p_{T}(#mu#mu)", 100, 0., 50.);      unit[ptpair] = "GeV"; 
    tmp.Form("y1_%i", i); h[i][y1] = new TH1D(tmp, "y_{1}(#mu#mu)",  48, -2.4, 2.4);      unit[y1] = "";
    tmp.Form("y2_%i", i); h[i][y2] = new TH1D(tmp, "y_{2}(#mu#mu)",  48, -2.4, 2.4);      unit[y2] = "";
    tmp.Form("pt1_%i", i); h[i][pt1] = new TH1D(tmp, "pt_{1}(#mu#mu)", 200, 0., 50.);  unit[pt1] = "GeV";
    tmp.Form("pt2_%i", i); h[i][pt2] = new TH1D(tmp, "pt_{2}(#mu#mu)", 200, 0., 50.);  unit[pt2] = "GeV";
    tmp.Form("mx_%i", i); h[i][mx] = new TH1D(tmp, "MX", 400, 0., 100.);       unit[mx] = "GeV";
  
    tmp.Form("t1_%i", i); h[i][t1] = new TH1D(tmp, "t_{1}", 100, 0., 10000.);    unit[t1] = "";

    tmp.Form("s1_%i", i); h[i][s1] = new TH1D(tmp, "s_{1}", 100, 0., 200000.);    unit[s1] = "";
    tmp.Form("th1_%i", i); h[i][th1] = new TH1D(tmp, "t_{h1}", 100, -1000000., 0.);    unit[th1] = "";
    tmp.Form("u1_%i", i); h[i][u1] = new TH1D(tmp, "u_{1}", 100, -1000000., 0.);    unit[u1] = "";
    tmp.Form("a1_%i", i); h[i][a1] = new TH1D(tmp, "a_{1}", 20, -2., 2.);    unit[a1] = "";
    tmp.Form("b1_%i", i); h[i][b1] = new TH1D(tmp, "b_{1}", 50, -100000., 100000.);    unit[b1] = "";

    tmp.Form("c1_%i", i); h[i][c1] = new TH1D(tmp, "c_{1}", 100, -1., 1.);    unit[c1] = "";

    tmp.Form("d1_%i", i); h[i][d1] = new TH1D(tmp, "d_{1}", 40, -4., 4.);    unit[d1] = "";


    tmp.Form("dphi_%i", i); h[i][dphi] = new TH1D(tmp, "|#Delta#phi(#mu,#mu)|", 30, 2., 3.2);    unit[dphi] = "";
  }

  Double_t px[9], py[9], pz[9], E[9];
  Double_t xsect;

  /////////////////////////////////////////////////////////////////////////////
  // PPTOLL (1)
  /////////////////////////////////////////////////////////////////////////////
  
  cout << "Processing PPTOLL (sample 1)" << endl;

  t_pptoll[0] = (TTree*)(f_pptoll[0]->Get("events"));
  t_pptoll[0]->SetBranchAddress("xsect", &xsect);
  t_pptoll[0]->SetBranchAddress("px", px);
  t_pptoll[0]->SetBranchAddress("py", py);
  t_pptoll[0]->SetBranchAddress("pz", pz);
  t_pptoll[0]->SetBranchAddress("E", E);

  for (Int_t i=0; i<max_events; i++) {
    t_pptoll[0]->GetEntry(i);
    weight = xsect/max_events;
    l1.SetXYZT(px[5], py[5], pz[5], E[5]);
    l2.SetXYZT(px[6], py[6], pz[6], E[6]);
    ip1.SetXYZT(px[0], py[0], pz[0], E[0]);
    ip2.SetXYZT(px[1], py[1], pz[1], E[1]);
    op1.SetXYZT(px[2], py[2], pz[2], E[2]);
    op2.SetXYZT(px[4], py[4], pz[4], E[4]);
    //if ((1.-fabs(l1.DeltaPhi(l2))/TMath::Pi()>0.1)) continue;
    //if (fabs(l1.Pt()-l2.Pt())>1.) continue;
    h[0][invm]->Fill((l1+l2).M(), weight);
    h[0][mx]->Fill(op1.M(), weight);
    h[0][ptpair]->Fill((l1+l2).Pt(), weight);
    h[0][dphi]->Fill(fabs(l1.DeltaPhi(l2)), weight);
    h[0][y1]->Fill(l1.Rapidity(), weight);
    h[0][y2]->Fill(l2.Rapidity(), weight);
    h[0][pt1]->Fill(l1.Pt(), weight);
    h[0][pt2]->Fill(l2.Pt(), weight);
    t1val= -(op1-ip1).M2();
    h[0][t1]->Fill(t1val, weight);
    s1val= (l1+l2).M2();
    h[0][s1]->Fill(s1val, weight);
    tval = (l1-(ip1-op1)).M2();
    h[0][th1]->Fill(tval, weight);
    u1val = (l1-(ip2-op2)).M2();
    h[0][u1]->Fill(u1val, weight);
//    shatval = (l1+l2).M2();
//    thatval = (l1-(ip1-op1)).M2();
//    uhatval = (l1-(ip2-op2)).M2();
//    sigval = s1val+tval+u1val
    mval = (l1.M2()+ l1.M2())/(l1.M2()+l1.M2());
    m1val = l1.M2();
//    sigval = s1val+tval+u1val- mval
//    allval = sigval - mval
    h[0][a1]->Fill(s1val+tval+u1val-mval, weight);
    h[0][b1]->Fill(s1val+tval+u1val, weight);
    
    amw2 = l1.M2();

//    t1val= -(op1-ip1).M2();
     t2val= -(op2-ip2).M2();


//    t1val = - q12
//    c1 = (that-uhat)/(shat-2.d0*am_W**2-q12-q22)

    c1val=(tval-u1val)/(s1val-2.*amw2-t1val-t2val);

    h[0][c1]->Fill(c1val, weight);

     dval1 = l1.Rapidity();
     dval2 = l2.Rapidity();

     h[0][d1]->Fill(dval1-dval2, weight);

//    h[0][uhatval]->Fill(uhatval, weight);
 }

 	Int_t nbins = h[0][mx]->GetNbinsX();
		          
    	ofstream output("output_MX.dat");
        
		          
    for (Int_t i=1; i<=nbins; i++) {
	    yy = (h[0][mx]->GetBinContent(i)/h[0][mx]->GetBinWidth(i));
// in pb 
    	output << h[0][mx]->GetBinLowEdge(i)+h[0][mx]->GetBinWidth(i)/2 << "\t" << yy << endl;
    }

    output.close();
    output.clear();

/////////////////////////////////////////////////////////////////////////////
//...................file dsig/dm_l+l-.....................................
/////////////////////////////////////////////////////////////////////////////
	nbins = h[0][invm]->GetNbinsX();
		          
    	output.open("output_Mll.dat");
        	
		          
    for (Int_t i=1; i<=nbins; i++) {
	    yy = (h[0][invm]->GetBinContent(i)/h[0][invm]->GetBinWidth(i));
// in pb 
//    	output << h[mx]->GetBinLowEdge(i)+h[mx]->GetBinWidth(i)/2 << "\t" << h[mx]->GetBinContent(i) << endl;
    	output << h[0][invm]->GetBinLowEdge(i)+h[0][invm]->GetBinWidth(i)/2 << "\t" << yy << endl;
    }

    output.close();
    output.clear();

//==========================================================================
//...................file dsig/dptpair.....................................
//==========================================================================
    nbins = h[0][ptpair]->GetNbinsX();
		          
    	output.open("output_ptpair.dat");
	
		          
    for (Int_t i=1; i<=nbins; i++) {
	    yy = (h[0][ptpair]->GetBinContent(i)/h[0][ptpair]->GetBinWidth(i));
// in pb 
//    	output << h[mx]->GetBinLowEdge(i)+h[mx]->GetBinWidth(i)/2 << "\t" << h[mx]->GetBinContent(i) << endl;
    	output << h[0][ptpair]->GetBinLowEdge(i)+h[0][ptpair]->GetBinWidth(i)/2 << "\t" << yy << endl;
    }
    output.close();
    output.clear();

//==========================================================================
//...................file dsig/dphi.....................................
//==========================================================================
	nbins = h[0][dphi]->GetNbinsX();
		          
    	output.open("output_dphi.dat");
	
		          
    for (Int_t i=1; i<=nbins; i++) {
	    yy = (h[0][dphi]->GetBinContent(i)/h[0][dphi]->GetBinWidth(i));
// in pb 
//    	output << h[mx]->GetBinLowEdge(i)+h[mx]->GetBinWidth(i)/2 << "\t" << h[mx]->GetBinContent(i) << endl;
    	output << h[0][dphi]->GetBinLowEdge(i)+h[0][dphi]->GetBinWidth(i)/2 << "\t" << yy << endl;
    }
    output.close();
    output.clear();
///==========================================================================
//...................file dsig/dt1.....................................
//==========================================================================
	nbins = h[0][t1]->GetNbinsX();
		          
    	output.open("output_t1.dat");
	
		          
    for (Int_t i=1; i<=nbins; i++) {
	    yy = (h[0][t1]->GetBinContent(i)/h[0][t1]->GetBinWidth(i));
// in pb 
//    	output << h[mx]->GetBinLowEdge(i)+h[mx]->GetBinWidth(i)/2 << "\t" << h[mx]->GetBinContent(i) << endl;
    	output << h[0][t1]->GetBinLowEdge(i)+h[0][t1]->GetBinWidth(i)/2 << "\t" << yy << endl;
    }	

    output.close();
    output.clear();
///==========================================================================
//...................file dsig/dt1.....................................
//==========================================================================
	nbins = h[0][a1]->GetNbinsX();
		          
    	output.open("output_a1.dat");
	
		          
    for (Int_t i=1; i<=nbins; i++) {
	    yy = (h[0][a1]->GetBinContent(i)/h[0][a1]->GetBinWidth(i));
// in pb 
//    	output << h[mx]->GetBinLowEdge(i)+h[mx]->GetBinWidth(i)/2 << "\t" << h[mx]->GetBinContent(i) << endl;
    	output << h[0][a1]->GetBinLowEdge(i)+h[0][a1]->GetBinWidth(i)/2 << "\t" << yy << endl;
    }	

    output.close();
    output.clear();


///==========================================================================
//...................file dsig/dshat.....................................
//==========================================================================
	nbins = h[0][s1]->GetNbinsX();
		          
    	output.open("output_s1.dat");
	
		          
    for (Int_t i=1; i<=nbins; i++) {
	    yy = (h[0][s1]->GetBinContent(i)/h[0][s1]->GetBinWidth(i));
// in pb 
//    	output << h[mx]->GetBinLowEdge(i)+h[mx]->GetBinWidth(i)/2 << "\t" << h[mx]->GetBinContent(i) << endl;
    	output << h[0][s1]->GetBinLowEdge(i)+h[0][s1]->GetBinWidth(i)/2 << "\t" << yy << endl;
    }	

    output.close();
    output.clear();
///==========================================================================
//...................file dsig/dtval.....................................
//==========================================================================
	nbins = h[0][th1]->GetNbinsX();
		          
    	output.open("output_th1.dat");
	
		          
    for (Int_t i=1; i<=nbins; i++) {
	    yy = (h[0][th1]->GetBinContent(i)/h[0][th1]->GetBinWidth(i));
// in pb 
//    	output << h[mx]->GetBinLowEdge(i)+h[mx]->GetBinWidth(i)/2 << "\t" << h[mx]->GetBinContent(i) << endl;
    	output << h[0][th1]->GetBinLowEdge(i)+h[0][th1]->GetBinWidth(i)/2 << "\t" << yy << endl;
    }	

    output.close();
    output.clear();

///==========================================================================
//...................file dsig/dshat.....................................
//==========================================================================
	nbins = h[0][u1]->GetNbinsX();
		          
    	output.open("output_u1.dat");
	
		          
    for (Int_t i=1; i<=nbins; i++) {
	    yy = (h[0][u1]->GetBinContent(i)/h[0][u1]->GetBinWidth(i));
// in pb 
//    	output << h[mx]->GetBinLowEdge(i)+h[mx]->GetBinWidth(i)/2 << "\t" << h[mx]->GetBinContent(i) << endl;
    	output << h[0][u1]->GetBinLowEdge(i)+h[0][u1]->GetBinWidth(i)/2 << "\t" << yy << endl;
    }

    output.close();
    output.clear();

///==========================================================================
//...................file dsig/dcostheta.....................................
//==========================================================================
	nbins = h[0][c1]->GetNbinsX();
		          
    	output.open("output_c1.dat");
	
		          
    for (Int_t i=1; i<=nbins; i++) {
	    yy = (h[0][c1]->GetBinContent(i)/h[0][c1]->GetBinWidth(i));
// in pb 
//    	output << h[mx]->GetBinLowEdge(i)+h[mx]->GetBinWidth(i)/2 << "\t" << h[mx]->GetBinContent(i) << endl;
    	output << h[0][c1]->GetBinLowEdge(i)+h[0][c1]->GetBinWidth(i)/2 << "\t" << yy << endl;
    }

    output.close();
    output.clear();

///==========================================================================
//...................file dsig/dydiff.....................................
//==========================================================================
	nbins = h[0][d1]->GetNbinsX();
		          
    	output.open("output_d1.dat");
	
		          
    for (Int_t i=1; i<=nbins; i++) {
	    yy = (h[0][d1]->GetBinContent(i)/h[0][d1]->GetBinWidth(i));
// in pb 
//    	output << h[mx]->GetBinLowEdge(i)+h[mx]->GetBinWidth(i)/2 << "\t" << h[mx]->GetBinContent(i) << endl;
    	output << h[0][d1]->GetBinLowEdge(i)+h[0][d1]->GetBinWidth(i)/2 << "\t" << yy << endl;
    }

    output.close();
    output.clear();

///==========================================================================
//...................file dsig/dshat.....................................
//==========================================================================
	nbins = h[0][b1]->GetNbinsX();
		          
    	output.open("output_b1.dat");
	
		          
    for (Int_t i=1; i<=nbins; i++) {
	    yy = (h[0][b1]->GetBinContent(i)/h[0][b1]->GetBinWidth(i));
// in pb 
//    	output << h[mx]->GetBinLowEdge(i)+h[mx]->GetBinWidth(i)/2 << "\t" << h[mx]->GetBinContent(i) << endl;
    	output << h[0][b1]->GetBinLowEdge(i)+h[0][b1]->GetBinWidth(i)/2 << "\t" << yy << endl;
    }

    output.close();
    output.clear();
	

//==========================================================================
//...................file dsig/dpt1.....................................
//==========================================================================
	nbins = h[0][pt1]->GetNbinsX();
		          
    	output.open("output_pt1.dat");
	
		          
    for (Int_t i=1; i<=nbins; i++) {
	    yy = (h[0][pt1]->GetBinContent(i)/h[0][pt1]->GetBinWidth(i));
// in pb 
//    	output << h[mx]->GetBinLowEdge(i)+h[mx]->GetBinWidth(i)/2 << "\t" << h[mx]->GetBinContent(i) << endl;
    	output << h[0][pt1]->GetBinLowEdge(i)+h[0][pt1]->GetBinWidth(i)/2 << "\t" << yy << endl;
    }

    output.close();
    output.clear();
//==========================================================================
//...................file dsig/dpt2.....................................
//==========================================================================
	nbins = h[0][pt2]->GetNbinsX();
		          
    	output.open("output_pt2.dat");
	
		          
    for (Int_t i=1; i<=nbins; i++) {
	    yy = (h[0][pt2]->GetBinContent(i)/h[0][pt2]->GetBinWidth(i));
// in pb 
//    	output << h[mx]->GetBinLowEdge(i)+h[mx]->GetBinWidth(i)/2 << "\t" << h[mx]->GetBinContent(i) << endl;
    	output << h[0][pt2]->GetBinLowEdge(i)+h[0][pt2]->GetBinWidth(i)/2 << "\t" << yy << endl;
    }

    output.close();
    output.clear();
//==========================================================================
//...................file dsig/dy1.....................................
//==========================================================================
	nbins = h[0][y1]->GetNbinsX();
		          
    	output.open("output_y1.dat");
	
		          
    for (Int_t i=1; i<=nbins; i++) {
	    yy = (h[0][y1]->GetBinContent(i)/h[0][y1]->GetBinWidth(i));
// in pb 
//    	output << h[mx]->GetBinLowEdge(i)+h[mx]->GetBinWidth(i)/2 << "\t" << h[mx]->GetBinContent(i) << endl;
    	output << h[0][y1]->GetBinLowEdge(i)+h[0][y1]->GetBinWidth(i)/2 << "\t" << yy << endl;
    }
    output.close();
    output.clear();
//==========================================================================
//...................file dsig/dy2.....................................
//==========================================================================
	nbins = h[0][y2]->GetNbinsX();
		          
    	output.open("output_y2.dat");
	
		          
    for (Int_t i=1; i<=nbins; i++) {
	    yy = (h[0][y2]->GetBinContent(i)/h[0][y2]->GetBinWidth(i));
// in pb 
//    	output << h[mx]->GetBinLowEdge(i)+h[mx]->GetBinWidth(i)/2 << "\t" << h[mx]->GetBinContent(i) << endl;
    	output << h[0][y2]->GetBinLowEdge(i)+h[0][y2]->GetBinWidth(i)/2 << "\t" << yy << endl;
    }

    output.close();
    output.clear();

  /////////////////////////////////////////////////////////////////////////////
  // PPTOLL (2)
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  // PLOTTING
  /////////////////////////////////////////////////////////////////////////////

  // We plot everything together...



  for (Int_t i=0; i<max_dists; i++) {
    c[i] = new TCanvas;

    h[0][i]->Draw();
    xtitle = h[0][i]->GetTitle();
    tmp.Form("%3.2f", (h[0][i]->GetXaxis()->GetXmax()-h[0][i]->GetXaxis()->GetXmin())/h[0][i]->GetXaxis()->GetNbins());
    ytitle = TString("d#sigma/d(")+h[0][i]->GetTitle()+TString(") (pb / ")+tmp;
    if (unit[i]!="") {
      xtitle += TString(" (")+unit[i]+TString(")");
      ytitle += TString(" ")+unit[i];
    }
    ytitle += TString(")");
    h[0][i]->GetXaxis()->SetTitle(xtitle);
    h[0][i]->GetYaxis()->SetTitle(ytitle);
    //h[0][i]->SetLineStyle(2);
    h[0][i]->SetTitle("");
    h[0][i]->SetLineColor(kBlack);
    max = -1.;
    for (Int_t j=0; j<2; j++) {
      //h[j][i]->SetLineWidth(2);
      max = TMath::Max(max, h[j][i]->GetMaximum());
    }
    h[0][i]->GetYaxis()->SetRangeUser(0.00005, max*1.1);

//    c[i]->SaveAs("output/"+TString(h[0][i]->GetName())".png");
//    c[i]->SaveAs("output/"+TString(h[0][i]->GetName())".pdf");

    c[i]->SaveAs(TString(h[0][i]->GetName())+".png");
    c[i]->SaveAs(TString(h[0][i]->GetName())+".eps");

  }

}
