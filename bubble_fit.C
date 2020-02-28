#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TF2.h"
#include "TFile.h"
#include "cstring"
#include "iostream"
#include "fstream"
#include "TMath.h"
#include "TROOT.h"

using namespace std;


//Float_t livetime = 788.018037037; //without 2019
Float_t livetime = 100.767767037; //Oct to Feb 2020
//Float_t livetime = 63.439809069;
//Float_t livetime = 124.022611437;
//Float_t livetime  = 921.322731437; //livetime until Feb, 2020
//Float_t livetime = 730.;
Float_t M = 20;
/*Float_t r = 1.8;
Float_t b =  pow((M/(3.14*r*r*0.878)),1./3.);
Float_t a = r*b ; // 2.08474;
//Float_t b = 1.61810;
//Float_t z0 = 0.0333;//0.085;
*/

Float_t a = 0;
Float_t z0 = 0;
Float_t bn = 0;
Float_t bp = 0;
Float_t b = 0;

Float_t a1 = 0;
Float_t z01 = 0;
Float_t bn1 = 0;
Float_t bp1 = 0;
Float_t b1 = 0;

Float_t a2 = 0;
Float_t z02 = 0;
Float_t bn2 = 0;
Float_t bp2 = 0;
Float_t b2 = 0;

//Float_t leak =0;
Float_t leak = 1.9*livetime*(10./50.)*(25./50.)*0.878/100.;//2.72051;
//	Float_t leak = 2.72051;
Float_t eff = 0.9712;//0.9699;
//	Float_t eff = 0.9699;
Float_t leak_err =0;//0.179538
///	Float_t leak_err =0.179538;

Float_t MT = 3.14*(10./50.)*(25./50.)*(livetime/100.)*0.878;
Double_t fitf(Double_t *x, Double_t *par){
	Double_t f;

	if ((x[0]/a2**2 + (x[1] - z02)**2/b2**2)  >  pow((M*3./(4.*3.14*a2*a2*b2*0.878)), 2./3.)){
		TF2::RejectPoint();
		return 0;
	}

	f=  MT*((par[0] + leak)/eff)*(1 + x[0]/par[1]**2 + (x[1]-par[2])**2/par[3]**2);
	return f;

}
Double_t fitf40(Double_t *x, Double_t *par){
	Double_t f;

	if ((x[0]/a1**2 + (x[1] - z01)**2/b1**2)  >  pow((40.*3./(4.*3.14*a1*a1*b1*0.878)), 2./3.)){
		TF2::RejectPoint();
		return 0;
	}

	f=  MT*((par[0] + leak)/eff)*(1 + x[0]/par[1]**2 + (x[1]-par[2])**2/par[3]**2);
	return f;

}
Double_t fitf70(Double_t *x, Double_t *par){
	Double_t f;

	if ((x[0]/a**2 + (x[1] - z0)**2/b**2)  >  pow((70.*3./(4.*3.14*a*a*b*0.878)), 2./3.)){
		TF2::RejectPoint();
		return 0;
	}

	f=  MT*((par[0] + leak)/eff)*(1 + x[0]/par[1]**2 + (x[1]-par[2])**2/par[3]**2);
	return f;

}
Double_t fitf_full(Double_t *x, Double_t *par){
	Double_t ff;
	if ((x[0]/2.**2 + (x[1])**2/1.**2)  >  pow((100.*3./(4.*3.14*2.*2.*1.*0.878)), 2./3.)){
		TF2::RejectPoint();
		return 0;
	}


	ff = MT*((par[0] + leak)/eff)*(1 + x[0]/par[1]**2 + (x[1]-par[2])**2/par[3]**2);
	return ff;

}
Double_t parabola(Double_t *x, Double_t *par){
	Double_t p;
	
	p = (par[0] + par[1]*(x[0]**2));
	return p;
}


void bubble_fit(){
//		TFile *f = new TFile("BiPo_aligned_201619_full_8m.root", "READ");
//		TFile *f = new TFile("aligned_data_2011-202001.root", "READ");
		TFile *f = new TFile("Po_201601-20200902_aligned_1902.root", "READ");
//		TFile *f = new TFile("Dtoymc.root", "READ");



		TH2D *rhoz =  new TH2D("rhoz", "rhoz", 50, 0, 25, 50, -5, 5);
		TH1D *plat = new TH1D("plat", "plat", 400, -4, 4);
//		rhoz->Sumw2();
		TTree *t;
		f->GetObject("t",t);
//		f->GetObject("bub",t);
//		t->Draw("z:(x*x+y*y)>>rhoz", "MLPv8 < 0.3 && Charge_Geo < 270 && Charge_Geo > 150 ");
		t->Draw("(z -z0):(x*x+y*y)>>rhoz", "MLPv8 < 0.3 && Charge_Geo < 270 && Charge_Geo > 150 && DstNumber >= 387");
//		t->Draw("(z -z0)>>plat", "MLPv8 < 0.3 && Charge_Geo < 270 && Charge_Geo > 150 && DstNumber >= 247 && (x**2 + y**2 < 0.96)");
//		t->Draw("z:(x*x+y*y)>>rhoz");
//		t->Draw("z>>plat");
		TF2 *func_full = new TF2("func_full",fitf_full, 0, 25, -5,5 ,4);	//1.209468
		func_full->SetParameters(10,1,0,1);
		//func_full->SetParameters(10,0.5,0,0.5,0.5);
		rhoz->Fit(func_full," L");

		a = func_full->GetParameter(1);
		b = func_full->GetParameter(3);
		z0 = func_full->GetParameter(2);


		TF2 *func_70 = new TF2("func_70",fitf70, 0, 25, -5,5 ,4 );	//1.209468
		func_70->SetParameters(10,1,0,1);
//		func_70->SetNpx(1000);
//		func_70->SetNpy(1000);
		rhoz->Fit(func_70,"L");

		a1 = func_70->GetParameter(1);
		b1 = func_70->GetParameter(3);
		z01 = func_70->GetParameter(2);

 		cout << func_70->GetChisquare()/func_70->GetNDF() << " " << func_70->GetProb()<< endl;
		TF2 *func_40 = new TF2("func_40",fitf40, 0, 25, -5,5 , 4);	//1.209468
		func_40->SetParameters(10, 1, 0,1);
		rhoz->Fit(func_40,"L");

		a2 = func_40->GetParameter(1);
		b2 = func_40->GetParameter(3);
		z02 = func_40->GetParameter(2);
 		cout << func_40->GetChisquare()/func_40->GetNDF() << " " << func_40->GetProb()<< endl;
		TF2 *func2 = new TF2("func2",fitf, 0,25, -5, 5, 4);
//		TF2 *func2 = new TF2("func2",fitf, 0,10, -1.5, 1.5, 4);

		func2->SetNpx(1000);
		func2->SetNpy(1000);
		func2->SetParameters(10,1,0,1);
		func2->SetParName(0, "R_{min}");
		func2->SetParName(1, "a");
		func2->SetParName(2, "z_{0}");
	//	func2->SetParName(2, "b+");
		func2->SetParName(3, "b");
//while(M < 35){
		TFitResultPtr fitres =  rhoz->Fit(func2, " L && S");
	Float_t chi2 =  func2->GetChisquare();
	Float_t ndf = func2->GetNDF();
	Float_t pvalue = func2->GetProb();
	Float_t fcn = fitres->MinFcnValue();
//	Float_t Rmin = (func2->GetParameter(0) - leak)/eff;
//	Float_t Rmin_err = sqrt(func2->GetParError(0)**2 + leak_err**2)/eff;
	cout << "The Rmin obtained is " << func2->GetParameter(0) << " +/- " << func2->GetParError(0) << " with a chi2/ndf of " << chi2/ndf << " and a p-value of " << pvalue << " and -2lnL of " << 2*fcn << endl;
//		ofstream file("fit_param_mass_1602.txt", ofstream::app);
//		file <<M << " " << func2->GetParameter(0) << " " << func2->GetParError(0) << " " << func2->GetParameter(1) << " " << func2->GetParError(1)<< " " << func2->GetParameter(2) << " " << func2->GetParError(2) << " " << func2->GetParameter(3) << " " << func2->GetParError(3) << " " << chi2/ndf << " " << pvalue << endl;
//M = M + 1;
//	}
		rhoz->Scale(100./(livetime*0.878));
		rhoz->Draw("COLZ");

/*		TH1D *plat = rhoz->ProjectionY("",2, 2);
		TF1 *fit = func2->PrjectionY("",2,2);
	plat->Scale(MT);
		TF1 *para = new TF1("para", parabola, -0.96, 0.96, 2);
		para->SetParameters(15,0.5);
		plat->Fit(para, "R && L");
	 cout << para->GetChisquare()/para->GetNDF() << " " << para->GetProb() << endl;	
	plat->Draw();
*/

// say something
int x123 = 5;		
		
}
