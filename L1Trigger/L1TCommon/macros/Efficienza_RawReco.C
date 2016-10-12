#include <set>
#include <algorithm>
#include <string>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <memory>
#include <TFile.h>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMultiGraph.h"

#include <iostream>

#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1ExtraDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoMuon2DataFormat.h"


void Efficienza_RawReco()
{
    //gROOT->Reset();
    //gROOT->SetBatch();

    TFile* f2 = new TFile("./L1_ZMu_RawReco_v3.root","READ");

    cout<<"ciao"<<endl;

    TTree* t_Reco = (TTree*) f2->Get("l1MuonRecoTree/Muon2RecoTree");
    L1Analysis::L1AnalysisRecoMuon2DataFormat    *Muon_ = new L1Analysis::L1AnalysisRecoMuon2DataFormat();
    t_Reco->SetBranchAddress("Muon", &Muon_);

    TTree* t_Trigger = (TTree*) f2->Get("l1ExtraTree/L1ExtraTree");
    L1Analysis::L1AnalysisL1ExtraDataFormat    *legacy_ = new L1Analysis::L1AnalysisL1ExtraDataFormat();
    t_Trigger->SetBranchAddress("L1Extra", &legacy_);



    //TH1F *eta_trigger = new TH1F("#eta_trigger","#eta_trigger", 80, -1, 1);;
    //eta_trigger->GetXaxis()->SetTitle("#eta");
    
    TH1F *eta_reco_match = new TH1F("#eta reco match","#eta reco match", 80, -1, 1);;
    eta_reco_match->GetXaxis()->SetTitle("#eta");
   
    TH1F *eta_reco = new TH1F("#eta reco","#eta reco", 80, -1, 1);;
    eta_reco->GetXaxis()->SetTitle("#eta");
   
    TH1F *efficiency =  new TH1F("efficiency","efficiency", 80, -1, 1);;
    efficiency->GetXaxis()->SetTitle("#eta");
    
    
    //TH1F *pt_trigger = new TH1F("p_{T}_trigger","p_{T}_trigger", 24, 0, 120);
    //pt_trigger->GetXaxis()->SetTitle("p_{T}");
    
    TH1F *pt_reco_match = new TH1F("p_{T} reco match","p_{T} reco match", 24, 0, 120);
    pt_reco_match->GetXaxis()->SetTitle("p_{T}");
   
    TH1F *pt_reco = new TH1F("p_{T} reco","p_{T} reco", 24, 0, 120);
    pt_reco->GetXaxis()->SetTitle("p_{T}");
   
    TH1F *efficiency_pt =  new TH1F("p_{T} efficiency","p_{T} efficiency", 24, 0, 120);
    efficiency_pt->GetXaxis()->SetTitle("p_{T}");
    
    
    TH2F *scatter_R = new TH2F("#DeltaR vs p_{T}","#DeltaR vs p_{T}", 100, 0, 100, 30, 0, 1.5);
    scatter_R->GetXaxis()->SetTitle("p_{T}");
    scatter_R->GetYaxis()->SetTitle("#DeltaR");
    scatter_R->SetMarkerSize(1);
    
    TH2F *scatter_eta_pt = new TH2F("p_{T} vs #eta","p_{T} vs #eta", 20, -1, 1, 40, 0, 200);
    scatter_eta_pt->GetYaxis()->SetTitle("p_{T}");
    scatter_eta_pt->GetXaxis()->SetTitle("#eta");
    scatter_eta_pt->SetMarkerSize(1);
    scatter_eta_pt->GetYaxis()->SetRangeUser(0, 200);
    
    
    TH1F *DeltaR_cut = new TH1F("#DeltaR_cut", "#DeltaR_cut", 30, 0, 1.5);
    DeltaR_cut->GetXaxis()->SetTitle("#DeltaR");
    DeltaR_cut->GetYaxis()->SetTitle("Cut efficiency");



    int j;
    int i;
    int h;
    float deltaR;
    float min;
    float reco_pt;
    float reco_eta;
    float reco_phi;

    TVector3 trigger_muon;
    TVector3 reco_muon;

    Long64_t nentries = t_Trigger->GetEntriesFast();
	cout<<nentries<<endl;
	
	for (Long64_t jentry=0; jentry<nentries;jentry++){
    	
    	if(jentry%100000 ==0)
    		cout<<jentry<<endl;
        
        t_Reco->GetEntry(jentry);
        t_Trigger->GetEntry(jentry);
        
        for(h = 0; h < Muon_->nMuons; h++){
        	if(Muon_->pt.at(h) < 12) continue;
        	eta_reco->Fill(Muon_->eta.at(h));
        	min = 100;
            reco_muon.SetPtThetaPhi(Muon_->pt.at(h), Muon_->eta.at(h), Muon_->phi.at(h));
        	for(j = 0; j < legacy_->muonEta.size(); j++){
        		trigger_muon.SetPtThetaPhi(legacy_->muonEt.at(j), legacy_->muonEta.at(j), legacy_->muonPhi.at(j));
        		deltaR = trigger_muon.DeltaR(reco_muon);
                if(deltaR < min)
                    min = deltaR;
            }
        	if(min < 1){
	        	eta_reco_match->Fill(Muon_->eta.at(h));
            	}
         }
    }

    for (Long64_t jentry=0; jentry<nentries;jentry++){  
     
         	if(jentry%100000 ==0)
    			cout<<jentry<<endl;
        
        t_Reco->GetEntry(jentry);
        t_Trigger->GetEntry(jentry);
          
         for(h = 0; h < Muon_->nMuons; h++){
        	if(abs(Muon_->eta.at(h)) > 0.8) continue;
        	pt_reco->Fill(Muon_->pt.at(h));
        	scatter_eta_pt->Fill(Muon_->eta.at(h), Muon_->pt.at(h));
        	min = 100;
            reco_muon.SetPtThetaPhi(Muon_->pt.at(h), Muon_->eta.at(h), Muon_->phi.at(h));
        	for(j = 0; j < legacy_->muonEta.size(); j++){
        		trigger_muon.SetPtThetaPhi(legacy_->muonEt.at(j), legacy_->muonEta.at(j), legacy_->muonPhi.at(j));
        		deltaR = trigger_muon.DeltaR(reco_muon);
                if(deltaR < min)
                    min = deltaR;
            }
        	if(min < 0.45){
	            pt_reco_match->Fill(Muon_->pt.at(h));
            	}
         }
    } //ciclo eventi
    

         

        
        /*
         for(j = 0; j < legacy_->muonEta.size(); j++){
            if(abs(legacy_->muonEta.at(j)) > 0.8) continue;
            //cout<<legacy_->muonEta.at(j)<<endl;
            min = 100;
            eta_trigger->Fill(legacy_->muonEta.at(j));
            trigger_muon.SetPtThetaPhi(legacy_->muonEt.at(j), legacy_->muonEta.at(j), legacy_->muonPhi.at(j));
            for(h = 0; h < Muon_->eta.size(); h++){
                reco_muon.SetPtThetaPhi(Muon_->pt.at(h), Muon_->eta.at(h), Muon_->phi.at(h));
                deltaR = trigger_muon.DeltaR(reco_muon);
                if(deltaR < min){
                    min = deltaR;
                    reco_eta = Muon_->eta.at(h);
                    reco_pt = Muon_->pt.at(h);
                    }
            }
            scatter_R->Fill(reco_pt, min);
            if(min < 0.45){
	            eta_reco->Fill(reco_eta);
            	eta_reco_match->Fill(legacy_->muonEta.at(j));
            	}
         }
		*/



  

        
        
 /*       for(i = 0; i < scatter_R->GetNbinsY(); i++){
        	float numerator = scatter_R->Integral(0, scatter_R->GetNbinsX(), 0, i);
        	float denominator = scatter_R->Integral(0, scatter_R->GetNbinsX(), 0, scatter_R->GetNbinsY());
        	float entry = numerator / denominator;
        	DeltaR_cut->Fill(i*0.1, entry);
        	}


    TCanvas* c_eta_trigger_prima = new TCanvas("#eta_trigger","#eta_trigger",140,30,700,500);
    eta_trigger->Draw();
    
    TCanvas* c_eta_reco_match = new TCanvas("#eta_reco_match","#eta_reco_match",140,30,700,500);
    eta_reco_match->Draw();

    TCanvas* c_eta_reco = new TCanvas("#eta_reco","#eta_reco",140,30,700,500);
    eta_reco->Draw();
    
    TCanvas* c_scatter_R = new TCanvas("scatter","scatter",140,30,700,500);
    scatter_R->Draw();
    
    TCanvas* c_DeltaR_cut = new TCanvas("#DeltaR_cut","#DeltaR_cut",140,30,700,500);
    DeltaR_cut->Draw();
    
    TCanvas* c_eta_sovra = new TCanvas("#eta_sovra","#eta_sovra",140,30,700,500);
    eta_reco_match->SetStats(0);
    eta_reco_match->Draw();
    eta_reco->SetStats(0);
    eta_reco->Draw("same");
    eta_reco_match->SetLineColor(kRed);
    TLegend* leg1 = new TLegend(0.7,0.9, 0.95, 1);
    leg1->AddEntry(eta_reco_match, "#eta trigger muon matched","lp");
    leg1->AddEntry(eta_reco,"#eta reco muon matched","lp");
    leg1->Draw();
*/
    
    
    eta_reco->Sumw2(kTRUE);
    eta_reco_match->Sumw2(kTRUE);
    efficiency->Divide(eta_reco_match, eta_reco,1,1,"B");
    
    TCanvas* c_eff = new TCanvas("efficiency","efficiency",140,30,700,500);
    efficiency->SetStats(0);
    efficiency->Draw();
    TLatex n3;
  	n3.SetNDC();
  	n3.SetTextFont(52);
  	n3.SetTextSize(0.04);
    n3.DrawLatex(0.55, 0.45, "Cut: #DeltaR < 0.45");
    n3.DrawLatex(0.55, 0.4, "p_{T} > 12");
    c_eff->Print("eta_efficiency_RawReco.pdf");
    
    pt_reco->Sumw2(kTRUE);
    pt_reco_match->Sumw2(kTRUE);
    efficiency_pt->Divide(pt_reco_match, pt_reco,1,1,"B");
    
    TCanvas* c_eff_pt = new TCanvas("efficiency p_{T}","efficiency p_{T}",140,30,700,500);
    efficiency_pt->SetStats(0);
    efficiency_pt->Draw();
    TLatex n4;
  	n4.SetNDC();
  	n4.SetTextFont(52);
  	n4.SetTextSize(0.04);
    n4.DrawLatex(0.65, 0.55, "Cut: #DeltaR < 0.45");
    n4.DrawLatex(0.65, 0.5, "|#eta| < 0.8");
    c_eff_pt->Print("pt_efficiency_RawReco.pdf");
    
    //TCanvas* c_eta_reco_match = new TCanvas("#eta reco match","#eta reco match",140,30,700,500);
    //eta_reco_match->Draw();

    //TCanvas* c_eta_reco = new TCanvas("#eta reco","#eta reco",140,30,700,500);
    //eta_reco->Draw();
    
    TCanvas* c_eta_sovra = new TCanvas("#eta sovra","#eta sovra",140,30,700,500);
    eta_reco_match->GetYaxis()->SetRangeUser(5e3, 35e3);
    eta_reco->GetYaxis()->SetRangeUser(3e3, 20e3);
    eta_reco_match->SetStats(0);
    eta_reco->Draw();
    eta_reco->SetStats(0);
    eta_reco_match->Draw("same");
    eta_reco_match->SetLineColor(kRed);
    TLegend* leg1 = new TLegend(0.7,0.9, 0.95, 1);
    leg1->AddEntry(eta_reco_match, "#eta reco muon matched","lp");
    leg1->AddEntry(eta_reco,"#eta reco muon","lp");
    leg1->Draw();
    c_eta_sovra->Print("eta_sovra_RawReco.pdf");
    
    
    //TCanvas* c_pt_reco_match = new TCanvas("p_{T} reco match","#eta reco match",140,30,700,500);
    //pt_reco_match->Draw();

    //TCanvas* c_pt_reco = new TCanvas("#p_{T} reco","#eta reco",140,30,700,500);
    //pt_reco->Draw();
    
    TCanvas* c_pt_sovra = new TCanvas("p_{T} sovra","p_{T} sovra",140,30,700,500);
    pt_reco_match->SetStats(0);
    pt_reco->Draw();
    pt_reco->SetStats(0);
    pt_reco_match->Draw("same");
    pt_reco_match->SetLineColor(kRed);
    TLegend* leg2 = new TLegend(0.7,0.9, 0.95, 1);
    leg2->AddEntry(pt_reco_match, "p_{T} reco muon matched","lp");
    leg2->AddEntry(pt_reco,"p_{T} reco muon","lp");
    leg2->Draw();
    c_pt_sovra->Print("pt_sovra_RawReco.pdf");
    
    TCanvas* c_scatter_eta_pt = new TCanvas("scatter","scatter",140,30,700,500);
    scatter_eta_pt->Draw();
    
    

    
 
 /*   TH1D* h2_ETA_GEM;
    TString NOME;
    float MIN, MAX;
    std::vector<float> MEAN_ETA, RMS_ETA;
    std::vector<float> SIGMA_RES_ETA_GEM;
    std::vector<float> SIGMA_RES_ETA_ERR_GEM;
    std::vector<float> RMS_RES_ETA_GEM;
    std::vector<float> RMS_RES_ETA_ERR_GEM;
    std::vector<float> ETA_GEM;
    std::vector<float> ETA_ERR_GEM;
    
    for(int i=1; i<=scatter_R->GetNbinsX(); i++){
      cout<<" ---------------------------------------------------------------------------- "<<i<<endl;
      NOME = Form("ETA_GEM_%d", i);
      TCanvas *ETA_CANVAS_GEM = new TCanvas(NOME, NOME, 200,10,700,500);
      h2_ETA_GEM = scatter_R->ProjectionX(NOME,i,i);
      MEAN_ETA.push_back(h2_ETA_GEM->GetMean());
      RMS_ETA.push_back(h2_ETA_GEM->GetRMS());
      MIN = MEAN_ETA.at(i-1) - RMS_ETA.at(i-1);
      MAX = MEAN_ETA.at(i-1) + RMS_ETA.at(i-1);
      TF1 *f1 = new TF1("f1","gaus",MIN,MAX);
      h2_ETA_GEM->Fit("f1","R");
      RMS_RES_ETA_GEM.push_back(h2_ETA_GEM->GetRMS());
      RMS_RES_ETA_ERR_GEM.push_back(h2_ETA_GEM->GetRMSError());
      SIGMA_RES_ETA_GEM.push_back(f1->GetParameter(2));
      SIGMA_RES_ETA_ERR_GEM.push_back(f1->GetParError(2));
      ETA_GEM.push_back(scatter_R->GetXaxis()->GetBinCenter(i));
      ETA_ERR_GEM.push_back(0.05);

      if(i==1)
          ETA_CANVAS_GEM->Print("./ETA_NoAging_GEM_RawReco.pdf[");

          ETA_CANVAS_GEM->Print("./ETA_NoAging_GEM_RawReco.pdf");

      if(i == scatter_R->GetNbinsX())
          ETA_CANVAS_GEM->Print("./ETA_NoAging_GEM_RawReco.pdf]");

  }
    TGraphErrors *gr1_GEM =  new TGraphErrors(scatter_R->GetNbinsX(), &(ETA_GEM.front()), &(SIGMA_RES_ETA_GEM.front()),
                                         &(ETA_ERR_GEM.front()), &(SIGMA_RES_ETA_ERR_GEM.front()));
    TCanvas *c1_GEM = new TCanvas("PtRes vs |Eta|: Gem - Sigma","PtRes vs |Eta|: Gem - Sigma",200,10,700,500);
    gr1_GEM->Draw("A*");
    gr1_GEM->GetXaxis()->SetTitle("|Eta|");
    gr1_GEM->GetYaxis()->SetTitle("(p_{T} Global - p_{T} Gen) / p_{T} Gen");
    gr1_GEM->SetTitle("PtGlobal Res vs |Eta|: Gem - Sigma");
    
    
    TGraphErrors *gr11_GEM =  new TGraphErrors(scatter_R->GetNbinsX(), &(ETA_GEM.front()), &(RMS_RES_ETA_GEM.front()),
                                         &(ETA_ERR_GEM.front()), &(RMS_RES_ETA_ERR_GEM.front()));
    TCanvas *c11_GEM = new TCanvas("PtRes vs |Eta|: Gem - Rms","PtRes vs |Eta|: Gem - Rms",200,10,700,500);
    gr11_GEM->Draw("A*");
    gr11_GEM->GetXaxis()->SetTitle("|Eta|");
    gr11_GEM->GetYaxis()->SetTitle("(p_{T} Global - p_{T} Gen) / p_{T} Gen");
    gr11_GEM->SetTitle("PtGlobal Res vs |Eta|: Gem - Rms");

*/







} //funzione principale








