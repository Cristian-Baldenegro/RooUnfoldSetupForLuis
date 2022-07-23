#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TInterpreter.h>
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TString.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TDatabasePDG.h>
#include <TGraph.h>
#include <TTree.h>
#include <TMath.h>
#endif







void DrawLatex(Float_t x, Float_t y, Int_t color, const char* text, Float_t textSize = 0.06)
{
  TLatex* latex = new TLatex(x, y, text);
  latex->SetNDC();
  latex->SetTextSize(textSize);
  latex->SetTextColor(color);
  latex->SetTextFont(42);
  latex->Draw();
}



void plotBasic_LundPlane(int flagdummy,std:: string file, std::string fileother, int bintruexj1, int bintruexj2, int bintruerg1, int bintruerg2, int iter_ref){
 
   TCanvas *canv1;
   TCanvas *canv2;
   TCanvas *canv3;
   TCanvas *canv4;
   TCanvas *canv5;
   TCanvas *canv6;
   TCanvas *canv7;
   TCanvas *canv8;
    TCanvas *canv9;
   TLegend *leg;
   TLegend *lego;

   TPad *pad;
   TPad *pad2;
   
   TH2D* raw;
   TH2D* smeared;
   TH1D* raw1;
   TH1D* ptraw1;
 
   

   TH1D* ptdet;
   TH2D* det;
   TH2D* htrue;
   TH1D* htrue1;
   TH2D* htrue_ineff;
   TH2D* pthtrue;
   TH1D* pthtrue1;
   TH1D* smeared1;
   
   TH2D* unfold2d[20];
   TH2D* fold2d;
   TH1D* unfold[20];
    TH1D* unfold_extreme1[20];
      TH1D* unfold_extreme2[20];
   TH1D* fold[20];

   TH1D* ptunfold[20];
   TH1D* ptfold[20];
   TH1D* denom;
   TH1D* denom1;
   TH1D* denom2;
   TH2D* data;
   TH1D* data1;
 
   TH1D* denompt;

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
   
    TFile* f = new TFile(file.c_str());
    TFile* fprior = new TFile(fileother.c_str());


    std:stringstream ssn;
    ssn << "correff"<< bintruexj1<<"-"<<bintruexj2;           
 
    if(flagdummy>=1) raw=(TH2D*)f->Get("smeared");
    if(flagdummy==0) raw=(TH2D*)f->Get("raw");
    smeared=(TH2D*)f->Get("smeared");
    data=(TH2D*)f->Get("raw"); 


	       raw1=(TH1D*)raw->ProjectionX("raw1",1,-1);
	       smeared1=(TH1D*)smeared->ProjectionX("smeared1",1,1);
               ptraw1=(TH1D*)raw->ProjectionY("ptraw1",1,-1);
	       data1=(TH1D*)data->ProjectionX("data1",1,1);
              	      
	       htrue=(TH2D*)f->Get("truef");
	       if(flagdummy==2) htrue=(TH2D*)f->Get("true");
	       htrue1=(TH1D*)htrue->ProjectionX("true1",bintruexj1,bintruexj2);
               htrue_ineff=(TH2D*)f->Get("true");
            
               TH1D* eff=htrue_ineff->ProjectionX("eff",bintruexj1,bintruexj2);
                TH1D* effdenom=htrue->ProjectionX("effdenom",bintruexj1,bintruexj2);
		eff->Divide(effdenom);


                TH1D* effxj=htrue_ineff->ProjectionY("effxj",bintruerg1,bintruerg2);
                TH1D* effdenomxj=htrue->ProjectionY("effdenomxj",bintruerg1,bintruerg2);
		effxj->Divide(effdenomxj);

	       pthtrue=(TH2D*)f->Get("truef");
               if(flagdummy==2) pthtrue=(TH2D*)f->Get("true");
               pthtrue1=(TH1D*)pthtrue->ProjectionY("true1pt",bintruerg1,bintruerg2);
	      
	       htrue1->Scale(1./htrue1->Integral(1,-1));                       
               htrue1->Scale(1,"width");
	      
	       pthtrue1->Scale(1./pthtrue1->Integral(1,-1));                       
               pthtrue1->Scale(1,"width");
	     	        

               //loop over iterations
	       for(Int_t j=0;j<15;j++){
		 unfold2d[j]=(TH2D*)f->Get(Form("Bayesian_Unfoldediter%d.root",j+1));
		 fold2d=(TH2D*)f->Get(Form("Bayesian_Foldediter%d.root",j+1));
                 unfold[j]=(TH1D*)unfold2d[j]->ProjectionX(Form("unfold%d",j),bintruexj1,bintruexj2,"");
                 unfold_extreme1[j]=(TH1D*)unfold2d[j]->ProjectionX(Form("unfold_extreme1%d",j),2,2,"");
                  unfold_extreme2[j]=(TH1D*)unfold2d[j]->ProjectionX(Form("unfold_extreme2%d",j),5,5,"");
		 fold[j]=(TH1D*)fold2d->ProjectionX(Form("fold%d",j),1,-1,"");
		 fold[j]->Divide(raw1);

                 ptunfold[j]=(TH1D*)unfold2d[j]->ProjectionY(Form("ptunfold%d",j),bintruerg1,bintruerg2);
		 ptfold[j]=(TH1D*)fold2d->ProjectionY(Form("ptfold%d",j),1,-1,"");
                 ptfold[j]->Divide(ptraw1); 
		 
		 
		 if(flagdummy!=2)unfold[j]->Divide(eff);
		 unfold[j]->SetLineColor(j+1);
		 fold[j]->SetLineColor(j+1);
		 unfold[j]->Scale(1./unfold[j]->Integral(1,-1));
		 unfold[j]->Scale(1,"width");
                 if(flagdummy!=2)ptunfold[j]->Divide(effxj); 
		 ptunfold[j]->Scale(1./ptunfold[j]->Integral(1,-1));
		 ptunfold[j]->SetLineColor(j+1);
		 ptfold[j]->SetLineColor(j+1);
		 ptunfold[j]->Scale(1,"width");
		 
		
	       }
	       


                TH1D* eff1=htrue_ineff->ProjectionX("eff1",2,2);
                TH1D* effdenom1=htrue->ProjectionX("effdenom1",2,2);
		eff1->Divide(effdenom1);

                TH1D* eff2=htrue_ineff->ProjectionX("eff2",5,5);
                TH1D* effdenom2=htrue->ProjectionX("effdenom2",5,5);
		eff2->Divide(effdenom1);
	      
  	       denom=(TH1D*)unfold[iter_ref]->Clone("denom");   
               denom1=(TH1D*)unfold_extreme1[iter_ref]->Clone("denom1"); 
               denom2=(TH1D*)unfold_extreme2[iter_ref]->Clone("denom2"); 

	       denom1->Divide(eff1);
               denom2->Divide(eff2);
                 denom1->Scale(1./denom1->Integral(1,-1));
		 denom1->Scale(1,"width");


      denom2->Scale(1./denom2->Integral(1,-1));
		 denom2->Scale(1,"width");
	      
	 ////////////////crossing of reg and statistical uncerts//////////

        TH1D *def1, *def2, *def3, *def4;
        TH2D *itera, *iterad, *iterau, *iterp;
        Double_t errprior, errreg1, errreg2, errreg, errstat, errtot;
        TH1D *histotot(0);
        TH1D *historeg(0);
	TH1D *histoprior(0);
	TH1D *histostat(0);
        histotot=new TH1D("histotot","histot",12,0,15);
        historeg=new TH1D("historeg","historeg",12,0,15);
	histoprior=new TH1D("histoprior","histoprior",12,0,15);
        histostat=new TH1D("histostat","histostat",12,0,15);


	  for(Int_t k=2;k<13;k++){
	    
	   
	  itera=(TH2D*)f->Get(Form("Bayesian_Unfoldediter%d.root",k));
	  iterad=(TH2D*)f->Get(Form("Bayesian_Unfoldediter%d.root",k-1));
          iterau=(TH2D*)f->Get(Form("Bayesian_Unfoldediter%d.root",k+2));
	  iterp=(TH2D*)fprior->Get(Form("Bayesian_Unfoldediter%d.root",k));
	  
            def1=(TH1D*)itera->ProjectionX(Form("def1_%i",k),bintruexj1,bintruexj2);
	    def2=(TH1D*)iterad->ProjectionX(Form("def2_%i",k),bintruexj1,bintruexj2);
	    def3=(TH1D*)iterau->ProjectionX(Form("def3_%i",k),bintruexj1,bintruexj2);
            def4=(TH1D*)iterp->ProjectionX(Form("def4_%i",k),bintruexj1,bintruexj2);

	    errprior=0;
	    errreg1=0;
	    errreg2=0;
	    errreg=0;
	    errstat=0;
	    errtot=0;
	    
	    for(Int_t i=1;i<=def1->GetNbinsX();i++){
	      errprior=errprior+TMath::Abs(def4->GetBinContent(i)-def1->GetBinContent(i))/def1->GetBinContent(i);
              errreg1=TMath::Abs(def2->GetBinContent(i)-def1->GetBinContent(i));
	      errreg2=TMath::Abs(def3->GetBinContent(i)-def1->GetBinContent(i));
	      errreg=errreg+TMath::Max(errreg1,errreg2)/def1->GetBinContent(i);
	      errstat=errstat+def1->GetBinError(i)/def1->GetBinContent(i);
	      errtot=TMath::Sqrt(errprior*errprior+errreg*errreg+errstat*errstat);
}
              
           
 	     histotot->SetBinContent(k,errtot);
	     historeg->SetBinContent(k,errreg);
             histoprior->SetBinContent(k,errprior);
	     histostat->SetBinContent(k,errstat);
	  }





          canv1= new TCanvas(Form("canvas1"),Form("canvas1") ,1100,1100);
          canv1->SetTicks();
  	  canv1->cd();
  	  pad=new TPad("pad0","this is pad",0,0,1,1);
  	  pad->SetFillColor(0);
	 
  	  pad->SetMargin(0.15,0.12,0.25,0.9);
  	  pad->Draw();
  	  pad->SetTicks(1,1);
  	  pad->cd();
  	   
    
           histotot->GetYaxis()->SetTitleOffset(0.9);
           histotot->GetXaxis()->SetTitleOffset(0.9);
  
           histotot->GetXaxis()->SetLabelFont(42);
           histotot->GetYaxis()->SetLabelFont(42);
           histotot->GetXaxis()->SetLabelSize(0.04);
           histotot->GetYaxis()->SetLabelSize(0.04); 
    
           histotot->GetXaxis()->SetTitleFont(42);
           histotot->GetYaxis()->SetTitleFont(42);
 
           histotot->GetXaxis()->SetTitleSize(0.065);
           histotot->GetYaxis()->SetTitleSize(0.065);
 
  	   histotot->GetXaxis()->SetTitle("Iterations");
  	   histotot->GetYaxis()->SetTitle("Summed errors");
	
      
      histotot->SetMarkerSize(1.3);
      histotot->SetMarkerStyle(21);
      histotot->SetMarkerColor(kBlack);
      histotot->SetLineColor(1);
      histotot->Draw("");
      historeg->SetLineColor(4);
      historeg->Draw("same");
      histostat->SetLineColor(3);
      histostat->Draw("same");
      histoprior->SetLineColor(2);
      histoprior->Draw("same");
      
     
   lego = new TLegend(0.6, 0.45, 0.75, 0.67);
   lego->SetBorderSize(0);
   lego->SetTextSize(0.025);
   lego->SetTextFont(42);
   lego->AddEntry(histotot,"total", "L");
   lego->AddEntry(historeg,"Regularization", "L");
   lego->AddEntry(histostat,"Statistical", "L");
   lego->AddEntry(histoprior,"Prior", "L");
   
   lego->Draw();
   lego->SetFillColor(0);  
   DrawLatex(0.2, 0.85, 1,Form("Bayes unfolded,x_{#gammaJ}^{true} bin [%i,%i]",bintruexj1,bintruexj2),0.03);
   DrawLatex(0.2, 0.8, 1,"PbPb #sqrt{#it{s}} = 5.02 TeV, E_{#gamma}>100 GeV",0.03); 
   DrawLatex(0.2, 0.75, 1,Form("Anti-#it{k}_{T} #it{R} = 0.2"),0.03); 
   canv1->SaveAs("IterChoice_rg_GJ.pdf");





  	  canv2= new TCanvas(Form("canvas2"),Form("canvas2") ,1100,1100);
          canv2->SetTicks();
  	  canv2->cd();
  	  pad=new TPad("pad0","this is pad",0,0,1,1);
  	  pad->SetFillColor(0);
	 
  	  pad->SetMargin(0.15,0.12,0.25,0.9);
  	  pad->Draw();
  	  pad->SetTicks(1,1);
  	  pad->cd();
  	
	     
     unfold[0]->GetYaxis()->SetTitleOffset(0.9);
     unfold[0]->GetXaxis()->SetTitleOffset(0.9);
  
   unfold[0]->GetXaxis()->SetLabelFont(42);
   unfold[0]->GetYaxis()->SetLabelFont(42);
   unfold[0]->GetXaxis()->SetLabelSize(0.04);
   unfold[0]->GetYaxis()->SetLabelSize(0.04); 
    
  unfold[0]->GetXaxis()->SetTitleFont(42);

 
  unfold[0]->GetXaxis()->SetTitleSize(0.065);
  unfold[0]->GetYaxis()->SetTitleSize(0.065);
  unfold[0]->GetXaxis()->SetRangeUser(-0.05,0.35);
   unfold[0]->GetYaxis()->SetRangeUser(0.,2);
   
  	   unfold[0]->GetXaxis()->SetTitle("R_{g}");
	   if(flagdummy>=1)  unfold[0]->GetYaxis()->SetTitle(Form("Unfolded/True"));
	   if(flagdummy==0)  unfold[0]->GetYaxis()->SetTitle(Form("Unfolded/Iter %i",iter_ref));
	   if(flagdummy>=1) unfold[0]->Divide(htrue1);
            if(flagdummy==0) unfold[0]->Divide(denom);
  	  unfold[0]->Draw("");
  	 
  	   leg = new TLegend(0.6, 0.28, 0.7, 0.53);
     	  leg->SetBorderSize(0);
          leg->SetTextSize(0.02);
          leg->SetTextFont(42);
	  
     
         
           leg->AddEntry(unfold[0]," iter 1", "PEL");
	   
    	  
     	
	  
  	  for(Int_t k=1;k<15;k++){
	 
  	    if(k==9) unfold[k]->SetLineColor(kMagenta);
	    if(flagdummy>=1) unfold[k]->Divide(htrue1);
             if(flagdummy==0) unfold[k]->Divide(denom);
  	    leg->AddEntry(unfold[k],Form(" iter %d",k+1), "PEL");
      unfold[k]->Draw("same");}

      	        leg->Draw();
      		leg->SetFillColor(0);
	  
 
   	  
    DrawLatex(0.2, 0.85, 1,Form("Bayes unfolded,  X_{#gamma,jet}^{true} bin [ %i,%i] ",bintruexj1,bintruexj2),0.03);
    DrawLatex(0.2, 0.8, 1,"Pb-Pb 5.02 TeV, E_{#gamma}>1009 GeV",0.03);
    DrawLatex(0.2, 0.75, 1,"Anti-k_{T}, R=0.2",0.03);
    canv2->SaveAs("UnfvsIterRgGJ.pdf");


  canv3= new TCanvas(Form("canvas3"),Form("canvas3") ,1100,1100);
          canv3->SetTicks();
  	  canv3->cd();
  	  pad2=new TPad("pad0","this is pad",0,0,1,1);
  	  pad2->SetFillColor(0);
	 
  	  pad2->SetMargin(0.15,0.12,0.25,0.9);
  	  pad2->Draw();
  	  pad2->SetTicks(1,1);
  	  pad2->cd();
  	     gPad->SetLogz();
     fold[0]->GetYaxis()->SetTitleOffset(0.9);
     fold[0]->GetXaxis()->SetTitleOffset(0.9);
  
   fold[0]->GetXaxis()->SetLabelFont(42);
   fold[0]->GetYaxis()->SetLabelFont(42);
   fold[0]->GetXaxis()->SetLabelSize(0.04);
   fold[0]->GetYaxis()->SetLabelSize(0.04); 
    
  fold[0]->GetXaxis()->SetTitleFont(42);
  fold[0]->GetYaxis()->SetTitleFont(42);
  fold[0]->GetYaxis()->SetRangeUser(0.,2);
  fold[0]->GetXaxis()->SetTitleSize(0.065);
  fold[0]->GetYaxis()->SetTitleSize(0.065);
  
  	   fold[0]->GetXaxis()->SetTitle("R_{g}");
  	  fold[0]->GetYaxis()->SetTitle("Folded/Raw");
  	  fold[0]->Draw("");

           lego = new TLegend(0.2, 0.3, 0.3, 0.55);
     	  lego->SetBorderSize(0);
          lego->SetTextSize(0.02);
          lego->SetTextFont(42);
	  
     
         
           lego->AddEntry(fold[0]," iter 1", "PEL");

	  
           for(Int_t k=1;k<15;k++){
  	      if(k==9) fold[k]->SetLineColor(kMagenta);
  	      lego->AddEntry(fold[k],Form(" iter %d",k+1), "PEL");
  	    fold[k]->Draw("same");} 
	  
  	   lego->Draw();
     	   lego->SetFillColor(0);
	   DrawLatex(0.4, 0.4, 1,Form("Refolding test x_{#gamma jet}^{det} bin [%i,%i]",1,-1),0.03);
   DrawLatex(0.4, 0.35, 1,"pp 5.02 TeV",0.03);
   DrawLatex(0.4, 0.3, 1,"Anti-k_{T}, R=0.4, SD zcut 0.4",0.03);
  	  canv3->SaveAs("RefoldvsIterRgGJ.pdf");


  
  	  denompt=(TH1D*)ptunfold[iter_ref]->Clone("denompt");
          canv4= new TCanvas(Form("canvas4"),Form("canvas4") ,1100,1100);
          canv4->SetTicks();
  	  canv4->cd();
  	  pad=new TPad("pad0","this is pad",0,0,1,1);
  	  pad->SetFillColor(0);
	 
  	  pad->SetMargin(0.15,0.12,0.25,0.9);
  	  pad->Draw();
  	  pad->SetTicks(1,1);
  	  pad->cd();
  
	     
     ptunfold[0]->GetYaxis()->SetTitleOffset(1.1);
     ptunfold[0]->GetXaxis()->SetTitleOffset(0.9);
  
   ptunfold[0]->GetXaxis()->SetLabelFont(42);
   ptunfold[0]->GetYaxis()->SetLabelFont(42);
   ptunfold[0]->GetXaxis()->SetLabelSize(0.04);
   ptunfold[0]->GetYaxis()->SetLabelSize(0.04); 
    
  ptunfold[0]->GetXaxis()->SetTitleFont(42);

 
  ptunfold[0]->GetXaxis()->SetTitleSize(0.065);
  ptunfold[0]->GetYaxis()->SetTitleSize(0.065);
 
  ptunfold[0]->GetXaxis()->SetRangeUser(40,120);
  ptunfold[0]->GetXaxis()->SetRangeUser(40,120);ptunfold[0]->GetYaxis()->SetRangeUser(0,2);
  	   ptunfold[0]->GetXaxis()->SetTitle("x_{#gammajet}");
	   if(flagdummy>=1) ptunfold[0]->GetYaxis()->SetTitle(Form("Unfolded/True"));
	   if(flagdummy==0) ptunfold[0]->GetYaxis()->SetTitle(Form("Unfolded/Iter%i",iter_ref));
	   if(flagdummy>=1)ptunfold[0]->Divide(pthtrue1);
           if(flagdummy==0)ptunfold[0]->Divide(denompt);
  	  ptunfold[0]->Draw("");

  	  // htrue1->Draw("same");
  	   leg = new TLegend(0.25, 0.3, 0.3, 0.55);
     	  leg->SetBorderSize(0);
          leg->SetTextSize(0.025);
          leg->SetTextFont(42);
	  
     
         
           leg->AddEntry(ptunfold[0]," iter 1", "PEL");
	   
    	  
     	
	  
  	  for(Int_t k=1;k<15;k++){
	   
  	    if(k==9) ptunfold[k]->SetLineColor(kMagenta);
  	    ptunfold[k]->Divide(denompt);
  	    leg->AddEntry(ptunfold[k],Form(" iter %d",k+1), "PEL");
  	   ptunfold[k]->Draw("same");}

     	        leg->Draw();
     		leg->SetFillColor(0);
  	   
	  	  
            DrawLatex(0.2, 0.85, 1,Form("Bayes unfolded, R_{g}^{true} bin  [%i,%i]",bintruerg1,bintruerg2),0.03);
    DrawLatex(0.2, 0.8, 1,"Pb-Pb 5.02 TeV, E_{#gamma}>1009 GeV",0.03);
    DrawLatex(0.2, 0.75, 1,"Anti-k_{T}, R=0.2",0.03);
	
  	  canv4->SaveAs("UnfoldvsIterxj.pdf");


          canv5= new TCanvas(Form("canvas5"),Form("canvas5") ,1100,1100);
          canv5->SetTicks();
  	  canv5->cd();
  	  pad2=new TPad("pad0","this is pad",0,0,1,1);
  	  pad2->SetFillColor(0);
	 
  	  pad2->SetMargin(0.15,0.12,0.25,0.9);
  	  pad2->Draw();
  	  pad2->SetTicks(1,1);
  	  pad2->cd();
  	     gPad->SetLogz();
     ptfold[0]->GetYaxis()->SetTitleOffset(0.9);
     ptfold[0]->GetXaxis()->SetTitleOffset(0.9);
  
   ptfold[0]->GetXaxis()->SetLabelFont(42);
   ptfold[0]->GetYaxis()->SetLabelFont(42);
   ptfold[0]->GetXaxis()->SetLabelSize(0.04);
   ptfold[0]->GetYaxis()->SetLabelSize(0.04); 
    
  ptfold[0]->GetXaxis()->SetTitleFont(42);
  ptfold[0]->GetYaxis()->SetTitleFont(42);
  ptfold[0]->GetYaxis()->SetRangeUser(0.6,1.6);
  ptfold[0]->GetXaxis()->SetTitleSize(0.065);
  ptfold[0]->GetYaxis()->SetTitleSize(0.065);
  
  	   ptfold[0]->GetXaxis()->SetTitle("x_{#gammaJ}");
  	  ptfold[0]->GetYaxis()->SetTitle("Folded/Raw");
  	  ptfold[0]->Draw("");

           lego = new TLegend(0.6, 0.55, 0.8, 0.85);
     	  lego->SetBorderSize(0);
          lego->SetTextSize(0.03);
          lego->SetTextFont(42);
	  
     
         
           lego->AddEntry(ptfold[0]," iter 1", "PEL");

	  
           for(Int_t k=1;k<15;k++){
  	      if(k==9) ptfold[k]->SetLineColor(kMagenta);
  	      lego->AddEntry(ptfold[k],Form(" iter %d",k+1), "PEL");
  	    ptfold[k]->Draw("same");} 
	  
  	   lego->Draw();
     	   lego->SetFillColor(0);
         
 DrawLatex(0.2, 0.85, 1,Form("Refolding test, R_{g}^{det} bin  [%i,%i]",1,-1),0.03);
             DrawLatex(0.2, 0.8, 1,"Pb-Pb 5.02 TeV, E_{#gamma}>1009 GeV",0.03);
             DrawLatex(0.2, 0.75, 1,"Anti-k_{T}, R=0.2",0.03);
          
  	  canv5->SaveAs("RefoldvsIterXJ.pdf");


      
          canv6= new TCanvas(Form("canvas6"),Form("canvas6") ,1100,1100);
          canv6->SetTicks();
  	  canv6->cd();
  	  pad2=new TPad("pad0","this is pad",0,0,1,1);
  	  pad2->SetFillColor(0);
	 
  	  pad2->SetMargin(0.15,0.12,0.25,0.9);
  	  pad2->Draw();
  	  pad2->SetTicks(1,1);
  	  pad2->cd();
	   
     raw->GetYaxis()->SetTitleOffset(0.9);
     raw->GetXaxis()->SetTitleOffset(0.9);
  
   raw->GetXaxis()->SetLabelFont(42);
   raw->GetYaxis()->SetLabelFont(42);
   raw->GetXaxis()->SetLabelSize(0.04);
   raw->GetYaxis()->SetLabelSize(0.04); 
    
  raw->GetXaxis()->SetTitleFont(42);
  raw->GetYaxis()->SetTitleFont(42);
  raw->GetYaxis()->SetRangeUser(0.,1.6);
  raw->GetXaxis()->SetTitleSize(0.065);
  raw->GetYaxis()->SetTitleSize(0.065);
  
  	  raw->GetYaxis()->SetTitle("x_{#gammaJ}");
  	  raw->GetXaxis()->SetTitle("R_{g}");
  	  raw->Draw("text");

          
  
  	  canv6->SaveAs("RawStats_GJ.pdf");




  	  canv7= new TCanvas(Form("canvas7"),Form("canvas7") ,1100,1100);
          canv7->SetTicks();
  	  canv7->cd();
  	  pad2=new TPad("pad0","this is pad",0,0,1,1);
  	  pad2->SetFillColor(0);
	 
  	  pad2->SetMargin(0.15,0.12,0.25,0.9);
  	  pad2->Draw();
  	  pad2->SetTicks(1,1);
  	  pad2->cd();
	   
     eff->GetYaxis()->SetTitleOffset(0.9);
     eff->GetXaxis()->SetTitleOffset(0.9);
  
   eff->GetXaxis()->SetLabelFont(42);
   eff->GetYaxis()->SetLabelFont(42);
   eff->GetXaxis()->SetLabelSize(0.04);
   eff->GetYaxis()->SetLabelSize(0.04); 
    
  eff->GetXaxis()->SetTitleFont(42);
  eff->GetYaxis()->SetTitleFont(42);
  eff->GetYaxis()->SetRangeUser(0.,2);
  eff->GetXaxis()->SetRangeUser(-0.05,0.35);
  eff->GetXaxis()->SetTitleSize(0.065);
  eff->GetYaxis()->SetTitleSize(0.065);
  
  	   eff->GetYaxis()->SetTitle("Kine Efficiency");
  	  eff->GetXaxis()->SetTitle("R_{g}");
  	  eff->Draw("");
    DrawLatex(0.2, 0.8, 1,Form("Kinematic efficiency for the true bin x_{#gammaJ} [%i,%i]",bintruexj1,bintruexj2),0.025);
	  canv7->SaveAs("Kineff_rg_GJ.pdf");
          

       
 canv8= new TCanvas(Form("canvas8"),Form("canvas8") ,1100,1100);
          canv8->SetTicks();
  	  canv8->cd();
  	  pad2=new TPad("pad0","this is pad",0,0,1,1);
  	  pad2->SetFillColor(0);
	 
  	  pad2->SetMargin(0.15,0.12,0.25,0.9);
  	  pad2->Draw();
  	  pad2->SetTicks(1,1);
  	  pad2->cd();
	   
     effxj->GetYaxis()->SetTitleOffset(0.9);
     effxj->GetXaxis()->SetTitleOffset(0.9);
  
   effxj->GetXaxis()->SetLabelFont(42);
   effxj->GetYaxis()->SetLabelFont(42);
   effxj->GetXaxis()->SetLabelSize(0.04);
   effxj->GetYaxis()->SetLabelSize(0.04); 
    
  effxj->GetXaxis()->SetTitleFont(42);
  effxj->GetYaxis()->SetTitleFont(42);
  effxj->GetYaxis()->SetRangeUser(0.,2);

  effxj->GetXaxis()->SetTitleSize(0.065);
  effxj->GetYaxis()->SetTitleSize(0.065);
  
  	   effxj->GetYaxis()->SetTitle("Kine Efficiency");
  	  effxj->GetXaxis()->SetTitle("x_{#gamma J}");
  	  effxj->Draw("");
    DrawLatex(0.2, 0.8, 1,Form("Kinematic efficiency for the true bin R_{g} [%i,%i]",bintruerg1,bintruerg2),0.025);
	  canv8->SaveAs("Kineff_xj_GJ.pdf");
          

     
canv9= new TCanvas(Form("canvas9"),Form("canvas9") ,1100,1100);
          canv9->SetTicks();
  	  canv9->cd();
  	  pad2=new TPad("pad0","this is pad",0,0,1,1);
  	  pad2->SetFillColor(0);
	 
  	  pad2->SetMargin(0.15,0.12,0.25,0.9);
  	  pad2->Draw();
  	  pad2->SetTicks(1,1);
  	  pad2->cd();
	   
     denom->GetYaxis()->SetTitleOffset(0.9);
     denom->GetXaxis()->SetTitleOffset(0.9);
  
   denom->GetXaxis()->SetLabelFont(42);
   denom->GetYaxis()->SetLabelFont(42);
   denom->GetXaxis()->SetLabelSize(0.04);
   denom->GetYaxis()->SetLabelSize(0.04); 
    
  denom->GetXaxis()->SetTitleFont(42);
  denom->GetYaxis()->SetTitleFont(42);
 denom->GetYaxis()->SetRangeUser(0.,15);

  denom->GetXaxis()->SetTitleSize(0.065);
  denom->GetYaxis()->SetTitleSize(0.065);
  
  if(flagdummy==0)  	   denom->GetYaxis()->SetTitle("Unfolded");
  if(flagdummy>=1)         denom->GetYaxis()->SetTitle("True MC");
  	  denom->GetXaxis()->SetTitle("R_{g}");
  	  denom->Draw("");
          denom->SetLineColor(1); 
          denom1->SetLineColor(1);
	  denom1->Draw("same");
          denom2->SetLineColor(2);
          denom2->Draw("same");
          smeared1->Scale(1./smeared1->Integral(1,-1)); 
			  data1->Scale(1./data1->Integral(1,-1));
	  //    DrawLatex(0.2, 0.8, 1,Form("UnfoldedSolution in x_{#gammaJ} bin [%i,%i]",bintruexj1,bintruexj2),0.025);
          smeared1->SetLineStyle(2);
	  smeared1->Draw("same");
	  data1->SetLineStyle(2);
	  data1->SetLineColor(2);
         data1->Draw("same");

    lego = new TLegend(0.6, 0.45, 0.75, 0.67);
    lego->SetBorderSize(0);
    lego->SetTextSize(0.025);
    lego->SetTextFont(42);
    lego->AddEntry(denom,"0.4<x_{#gamma j}<0.6", "L");
    lego->AddEntry(denom2,"1<x_{#gamma j} <1.2", "L");
    

    lego->Draw();
    lego->SetFillColor(0);
	  canv9->SaveAs("Rg_unfolded.pdf");
          

	  
	  
	  
}

