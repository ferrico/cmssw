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
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoMuon2DataFormat.h"


void Efficienza_EmuReco()
{
    //gROOT->Reset();
    //gROOT->SetBatch();

    TFile* f2 = new TFile("/cmshome/filippoe/trigger/CMSSW_8_0_0_pre6/src/L1Ntuple_273450_SingleMuon.root","READ");
    TFile* f = new TFile("/cmshome/filippoe/trigger/CMSSW_8_0_0_pre6/src//L1_ZMu_RawReco_v3.root","READ");
    
    TFile* Efficiency = new TFile("./Efficiency_273450.root", "recreate");
    Efficiency->cd();

	
    TTree* t_Reco = (TTree*) f2->Get("l1MuonRecoTree/Muon2RecoTree");
    L1Analysis::L1AnalysisRecoMuon2DataFormat    *Muon_ = new L1Analysis::L1AnalysisRecoMuon2DataFormat();
    t_Reco->SetBranchAddress("Muon", &Muon_);

    TTree* t_Trigger = (TTree*) f2->Get("l1UpgradeTree/L1UpgradeTree");
    L1Analysis::L1AnalysisL1UpgradeDataFormat    *upgrade_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
    t_Trigger->SetBranchAddress("L1Upgrade", &upgrade_);
    
    TTree* t_Extra = (TTree*) f->Get("l1ExtraTree/L1ExtraTree");
    L1Analysis::L1AnalysisL1ExtraDataFormat    *legacy_ = new L1Analysis::L1AnalysisL1ExtraDataFormat();
    t_Extra->SetBranchAddress("L1Extra", &legacy_);
    
    TTree* t_Reco_Raw = (TTree*) f->Get("l1MuonRecoTree/Muon2RecoTree");
    L1Analysis::L1AnalysisRecoMuon2DataFormat    *Muon2_ = new L1Analysis::L1AnalysisRecoMuon2DataFormat();
    t_Reco_Raw->SetBranchAddress("Muon", &Muon2_);


    Double_t PT_BINS[] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85};
    Int_t  binnum = sizeof(PT_BINS)/sizeof(Double_t)-1;
    
    TH1F *eta_reco_match = new TH1F("#eta reco match","#eta reco match", 200,-5, 5);
    eta_reco_match->GetXaxis()->SetTitle("#eta");
    
    TH1F *eta_reco = new TH1F("#eta reco","#eta reco", 200,-5, 5);
    eta_reco->GetXaxis()->SetTitle("#eta");
    
    TH1F *eta_reco_match_Raw = new TH1F("#eta reco match Raw","#eta reco match Raw", 200,-5, 5);
    eta_reco_match_Raw->GetXaxis()->SetTitle("#eta");
   
	TH1F *eta_reco_Raw = new TH1F("#eta reco Raw","#eta reco Raw", 200,-5, 5);
    eta_reco_Raw->GetXaxis()->SetTitle("#eta");
   
    TH1F *efficiency_eta =  new TH1F("#eta efficiency","#eta efficiency", 200, -5, 5);
    efficiency_eta->GetXaxis()->SetTitle("#eta");
    
    TH1F *efficiency_eta_Raw =  new TH1F("#eta efficiency Raw","#eta efficiency Raw", 200, -5, 5);
    efficiency_eta_Raw->GetXaxis()->SetTitle("#eta");
    
    
    
    TH1F *pt_reco_match= new TH1F("p_{T} reco match","p_{T} reco match", binnum, PT_BINS);
    pt_reco_match->GetXaxis()->SetTitle("p_{T}");
    
    TH1F *pt_reco = new TH1F("p_{T} reco","p_{T} reco", binnum, PT_BINS);
    pt_reco->GetXaxis()->SetTitle("p_{T}");
    
    TH1F *pt_reco_match_Raw = new TH1F("p_{T} reco match Raw","p_{T} reco match Raw", binnum, PT_BINS);
    pt_reco_match_Raw->GetXaxis()->SetTitle("p_{T}");

    TH1F *pt_reco_Raw = new TH1F("p_{T} reco Raw","p_{T} reco Raw", binnum, PT_BINS);
    pt_reco_Raw->GetXaxis()->SetTitle("p_{T}");   

    TH1F *efficiency_pt =  new TH1F("p_{T} efficiency","p_{T} efficiency", binnum, PT_BINS);
    efficiency_pt->GetXaxis()->SetTitle("p_{T}");
    
    TH1F *efficiency_pt_Raw =  new TH1F("p_{T} efficiency Raw","p_{T} efficiency Raw", binnum, PT_BINS);
    efficiency_pt_Raw->GetXaxis()->SetTitle("p_{T}");
    
    
  
    TH1F *eta_Emu = new TH1F("#eta Emu ","#eta Emu ", 200,-5, 5);
    eta_Emu->GetXaxis()->SetTitle("#eta");
  
    TH1F *eta_Raw= new TH1F("#eta Raw ","#eta Raw ", 200,-5, 5);
    eta_Raw->GetXaxis()->SetTitle("#eta");
    
    TH1F *eta_Emu_match = new TH1F("#eta Emu match","#eta Emu match", 200,-5, 5);
    eta_Emu_match->GetXaxis()->SetTitle("#eta");
  
    TH1F *eta_Raw_match = new TH1F("#eta Raw match","#eta Raw match", 200,-5, 5);
    eta_Raw_match->GetXaxis()->SetTitle("#eta");
    
    TH1F *pt_Emu = new TH1F("p_{T} Emu ","p_{T} Emu ", binnum, PT_BINS);
    pt_Emu->GetXaxis()->SetTitle("p_{T}");
  
    TH1F *pt_Raw= new TH1F("p_{T} Raw ","p_{T} Raw ", binnum, PT_BINS);
    pt_Raw->GetXaxis()->SetTitle("p_{T}");
    
    TH1F *pt_Emu_match = new TH1F("p_{T} Emu match","p_{T} Emu match", binnum, PT_BINS);
    pt_Emu_match->GetXaxis()->SetTitle("p_{T}");
  
    TH1F *pt_Raw_match = new TH1F("p_{T} Raw match","p_{T} Raw match", binnum, PT_BINS);
    eta_Raw_match->GetXaxis()->SetTitle("p_{T}");
    
    
    
    
    TH2F *scatter_R = new TH2F("#DeltaR vs p_{T}","#DeltaR vs p_{T}", 100, 0, 100, 30, 0, 1.5);
    scatter_R->GetXaxis()->SetTitle("p_{T}");
    scatter_R->GetYaxis()->SetTitle("#DeltaR");
    scatter_R->SetMarkerSize(1);
    
    TH2F *scatter_eta_pt = new TH2F("p_{T} vs #eta","p_{T} vs #eta", 20, -1, 1, 40, 0, 200);
    scatter_eta_pt->GetYaxis()->SetTitle("p_{T}");
    scatter_eta_pt->GetXaxis()->SetTitle("#eta");
    scatter_eta_pt->SetMarkerSize(1);
    scatter_eta_pt->GetYaxis()->SetRangeUser(0, 200);
    
    TH2F *Raw_scatter_eta_pt = new TH2F("Raw p_{T} vs #eta","Raw  p_{T} vs #eta", 20, -1, 1, 40, 0, 200);
    Raw_scatter_eta_pt->GetYaxis()->SetTitle("p_{T}");
    Raw_scatter_eta_pt->GetXaxis()->SetTitle("#eta");
    Raw_scatter_eta_pt->SetMarkerSize(1);
    Raw_scatter_eta_pt->GetYaxis()->SetRangeUser(0, 200);
    
    
    TH1F *DeltaR_cut = new TH1F("#DeltaR_cut", "#DeltaR_cut", 30, 0, 1.5);
    DeltaR_cut->GetXaxis()->SetTitle("#DeltaR");
    DeltaR_cut->GetYaxis()->SetTitle("Cut efficiency");

  
    int j;
    int i;
    int h;
    float deltaR;
    float min;
    float emu_eta;
    float raw_eta;
    float emu_pt;
    float raw_pt;

    TVector3 trigger_muon;
    TVector3 reco_muon;
    TVector3 Raw_muon;
    TVector3 reco_muon_Raw;


    // Raw - Reco   

/*	Long64_t nentries_Raw = t_Extra->GetEntriesFast();
	cout<<nentries_Raw<<endl;
	
	for (Long64_t jentry = 0; jentry<nentries_Raw; jentry++){
    	
    	if(jentry%100000 == 0) cout<<jentry<<endl;
        
        t_Reco_Raw->GetEntry(jentry);
        t_Extra->GetEntry(jentry);
        
        
        for(j = 0; j < legacy_->muonEt.size(); j++){
        	if(abs(legacy_->muonEta.at(j)) > 0.8) continue;
        	pt_Raw->Fill(legacy_->muonEt.at(j));
        	min = 100;
            Raw_muon.SetPtEtaPhi(legacy_->muonEt.at(j), legacy_->muonEta.at(j), legacy_->muonPhi.at(j));
        	for(h = 0; h < Muon2_->pt.size(); h++){
        		//pt_reco_Raw->Fill(Muon2_->pt.at(h));
        		reco_muon_Raw.SetPtEtaPhi(Muon2_->pt.at(h), Muon2_->eta.at(h), Muon2_->phi.at(h));
        		deltaR = Raw_muon.DeltaR(reco_muon_Raw);
                if(deltaR < min){
                    min = deltaR;
                    raw_pt = Muon2_->pt.at(h);
                }
            }
        	if(min < 0.3){
	        	//eta_reco_match_Raw->Fill(raw_eta);
	        	pt_Raw_match->Fill(legacy_->muonEt.at(j));
            }
         }    
    }	
    
	for (Long64_t jentry = 0; jentry<nentries_Raw; jentry++){
    	
    	if(jentry%100000 == 0) cout<<jentry<<endl;
        
        t_Reco_Raw->GetEntry(jentry);
        t_Extra->GetEntry(jentry);
        
        
        for(j = 0; j < legacy_->muonEta.size(); j++){
        	if(legacy_->muonEt.at(j) < 12) continue;
        	eta_Raw->Fill(legacy_->muonEta.at(j));
        	min = 100;
            Raw_muon.SetPtEtaPhi(legacy_->muonEt.at(j), legacy_->muonEta.at(j), legacy_->muonPhi.at(j));
        	for(h = 0; h < Muon2_->eta.size(); h++){
        		eta_reco_Raw->Fill(Muon2_->eta.at(h));
        		reco_muon_Raw.SetPtEtaPhi(Muon2_->pt.at(h), Muon2_->eta.at(h), Muon2_->phi.at(h));
        		deltaR = Raw_muon.DeltaR(reco_muon_Raw);
                if(deltaR < min){
                    min = deltaR;
                    raw_eta = Muon2_->eta.at(h);
                }
            }
        	if(min < 0.3){
	        	eta_reco_match_Raw->Fill(raw_eta);
	        	eta_Raw_match->Fill(legacy_->muonEta.at(j));
            }
         }
         		   		
         
    }

	for (Long64_t jentry = 0; jentry<nentries_Raw; jentry++){
	    if(jentry%100000 == 0) cout<<jentry<<endl;
        t_Reco_Raw->GetEntry(jentry);
        t_Extra->GetEntry(jentry);
        for(h = 0; h < legacy_->nMuons; h++)
        	Raw_scatter_eta_pt->Fill(legacy_->muonEta.at(h), legacy_->muonEt.at(h));
	} 
*/		
    // Emu - Reco

    Long64_t nentries = t_Trigger->GetEntriesFast();
	cout<<nentries<<endl;

    for (Long64_t jentry = 0; jentry<nentries;jentry++){
    	
    	if(jentry%100000 == 0)
    		cout<<jentry<<endl;
        
        t_Reco->GetEntry(jentry);
        t_Trigger->GetEntry(jentry);
        
        for(j = 0; j < upgrade_->muonEt.size(); j++){
        	if(abs(upgrade_->muonEta.at(j)) > 0.8) continue;
        	pt_Emu->Fill(upgrade_->muonEt.at(j));
        	min = 100;
	        trigger_muon.SetPtEtaPhi(upgrade_->muonEt.at(j), upgrade_->muonEta.at(j), upgrade_->muonPhi.at(j));
	        for(h = 0; h < Muon_->pt.size(); h++){
	        	//eta_reco->Fill(Muon_->eta.at(h));
	            reco_muon.SetPtEtaPhi(Muon_->pt.at(h), Muon_->eta.at(h), Muon_->phi.at(h));
        		deltaR = trigger_muon.DeltaR(reco_muon);
                if(deltaR < min){
                    min = deltaR;
                    emu_pt = Muon_->pt.at(h);
                }
            }
        	if(min < 0.3){
	        	//eta_reco_match->Fill(emu_eta);
	        	pt_Emu_match->Fill(upgrade_->muonEt.at(j));
            }
         }
    }
 
    for (Long64_t jentry = 0; jentry<nentries;jentry++){
    	
    	if(jentry%100000 == 0)
    		cout<<jentry<<endl;
        
        t_Reco->GetEntry(jentry);
        t_Trigger->GetEntry(jentry);
        
        for(j = 0; j < upgrade_->muonEta.size(); j++){
        	if(upgrade_->muonEt.at(j) < 12) continue;
        	eta_Emu->Fill(upgrade_->muonEta.at(j));
        	min = 100;
	        trigger_muon.SetPtEtaPhi(upgrade_->muonEt.at(j), upgrade_->muonEta.at(j), upgrade_->muonPhi.at(j));
	        for(h = 0; h < Muon_->eta.size(); h++){
	        	eta_reco->Fill(Muon_->eta.at(h));
	            reco_muon.SetPtEtaPhi(Muon_->pt.at(h), Muon_->eta.at(h), Muon_->phi.at(h));
        		deltaR = trigger_muon.DeltaR(reco_muon);
                if(deltaR < min){
                    min = deltaR;
                    emu_eta = Muon_->eta.at(h);
                }
            }
        	if(min < 0.3){
	        	eta_reco_match->Fill(emu_eta);
	        	eta_Emu_match->Fill(upgrade_->muonEta.at(j));
            }
         }
    }

    for (Long64_t jentry = 0; jentry<nentries;jentry++){
    	if(jentry%100000 == 0) cout<<jentry<<endl;      
        t_Reco->GetEntry(jentry);
        t_Trigger->GetEntry(jentry);  
        for(h = 0; h < upgrade_->nMuons; h++)
        	scatter_eta_pt->Fill(upgrade_->muonEta.at(h), upgrade_->muonEt.at(h)); 
    }   
 
/*TCanvas* c_eta_emu = new TCanvas("#eta emu  ","#eta emu  ",140,30,700,500);
eta_Emu_match->SetStats(0);
eta_Emu->SetStats(0);
eta_reco->SetStats(0);
eta_Emu_match->SetLineColor(kRed);
eta_Emu->SetLineColor(kGreen);
eta_Emu_match->Draw();
eta_Emu->Draw("same");
//eta_reco->Draw();
TLegend* leg2222 = new TLegend(0.7,0.9, 0.95, 1);
leg2222->AddEntry(eta_reco, "#eta Reco","lp");
leg2222->AddEntry(eta_Emu,"#eta Emu","lp");
leg2222->AddEntry(eta_Emu_match,"#eta Emu match","lp");
leg2222->Draw();
c_eta_emu->Print("eta_emu.pdf");

TCanvas* c_eta_raw = new TCanvas("#eta raw  ","#eta raw  ",140,30,700,500);
eta_Raw_match->SetStats(0);
eta_Raw->SetStats(0);
eta_reco_Raw->SetStats(0); 
eta_Raw_match->SetLineColor(kRed);
eta_Raw->SetLineColor(kGreen);
eta_Raw_match->Draw();
eta_Raw->Draw("same");
//eta_reco_Raw->Draw("same"); 
TLegend* leg22222 = new TLegend(0.7,0.9, 0.95, 1);
leg22222->AddEntry(eta_reco_Raw, "#eta Reco","lp");
leg22222->AddEntry(eta_Raw,"#eta Raw","lp");
leg22222->AddEntry(eta_Raw_match,"#eta Raw match","lp");
leg22222->Draw();
c_eta_raw->Print("eta_raw.pdf");  

  
    

    eta_reco->Sumw2(kTRUE);
    eta_reco_Raw->Sumw2(kTRUE);
    eta_reco_match->Sumw2(kTRUE);
    eta_reco_match_Raw->Sumw2(kTRUE);
    efficiency_eta->Divide(eta_reco_match, eta_reco,1,1,"B");
    efficiency_eta_Raw->Divide(eta_reco_match_Raw, eta_reco_Raw,1,1,"B");
    */
    eta_Raw->Sumw2(kTRUE);
    eta_Raw_match->Sumw2(kTRUE);
    eta_Emu->Sumw2(kTRUE);
    eta_Emu_match->Sumw2(kTRUE);
    efficiency_eta->Divide(eta_Emu_match, eta_Emu,1,1,"B");
    efficiency_eta_Raw->Divide(eta_Raw_match, eta_Raw,1,1,"B");
          
    TCanvas* c_eff = new TCanvas("#eta efficiency Emu Reco","#eta efficiency Emu Reco",140,30,700,500);
    efficiency_eta->SetStats(0);
    efficiency_eta->GetXaxis()->SetRangeUser(-1., 1.);
    efficiency_eta->Draw();
    TLatex n3;
  	n3.SetNDC();
  	n3.SetTextFont(52);
  	n3.SetTextSize(0.04);
    n3.DrawLatex(0.55, 0.45, "Cut: #DeltaR < 0.45");
    n3.DrawLatex(0.55, 0.4, "p_{T} > 12");
    c_eff->Print("eta_efficiency_EmuReco.pdf");

	TCanvas* c_eff_Raw = new TCanvas("#eta efficiency Raw Reco","#eta efficiency Raw Reco",140,30,700,500);
    efficiency_eta_Raw->SetStats(0);
    efficiency_eta_Raw->GetXaxis()->SetRangeUser(-1., 1.);
    efficiency_eta_Raw->Draw();
    TLatex n33;
  	n33.SetNDC();
  	n33.SetTextFont(52);
  	n33.SetTextSize(0.04);
    n33.DrawLatex(0.55, 0.45, "Cut: #DeltaR < 0.45");
    n33.DrawLatex(0.55, 0.4, "p_{T} > 12");
    c_eff_Raw->Print("eta_efficiency_RawReco.pdf");
    /*
    pt_reco->Sumw2(kTRUE);
    pt_reco_Raw->Sumw2(kTRUE);
    pt_reco_match->Sumw2(kTRUE);
    pt_reco_match_Raw->Sumw2(kTRUE);
    efficiency_pt->Divide(pt_reco_match, pt_reco,1,1,"B");
    efficiency_pt_Raw->Divide(pt_reco_match_Raw, pt_reco_Raw,1,1,"B");
    */
    
    pt_Raw->Sumw2(kTRUE);
    pt_Raw_match->Sumw2(kTRUE);
    pt_Emu->Sumw2(kTRUE);
    pt_Emu_match->Sumw2(kTRUE);
    efficiency_pt->Divide(pt_Emu_match, pt_Emu,1,1,"B");
    efficiency_pt_Raw->Divide(pt_Raw_match, pt_Raw,1,1,"B");
    
    TCanvas* c_eff_pt = new TCanvas("p_{T} efficiency Emu Reco","p_{T} efficiency Emu Reco",140,30,700,500);
    efficiency_pt->SetStats(0);
    efficiency_pt->Draw();
    TLatex n4;
  	n4.SetNDC();
  	n4.SetTextFont(52);
  	n4.SetTextSize(0.04);
    n4.DrawLatex(0.65, 0.55, "Cut: #DeltaR < 0.45");
    n4.DrawLatex(0.65, 0.5, "|#eta| < 0.8");
    c_eff_pt->Print("pt_efficiency_EmuReco.pdf");
    
    TCanvas* c_eff_pt_Raw = new TCanvas("p_{T} efficiency Raw Reco","p_{T} efficiency Raw Reco",140,30,700,500);
    efficiency_pt_Raw->SetStats(0);
    efficiency_pt_Raw->Draw();
    TLatex n44;
  	n44.SetNDC();
  	n44.SetTextFont(52);
  	n44.SetTextSize(0.04);
    n44.DrawLatex(0.65, 0.55, "Cut: #DeltaR < 0.45");
    n44.DrawLatex(0.65, 0.5, "|#eta| < 0.8");
    c_eff_pt_Raw->Print("pt_efficiency_RawReco.pdf");
    
    //TCanvas* c_eta_reco_match = new TCanvas("#eta reco match","#eta reco match",140,30,700,500);
    //eta_reco_match->Draw();

    //TCanvas* c_eta_reco = new TCanvas("#eta reco","#eta reco",140,30,700,500);
    //eta_reco->Draw();
    
    /*TCanvas* c_eta_sovra = new TCanvas("#eta sovra","#eta sovra",140,30,700,500);
    eta_reco_match->GetYaxis()->SetRangeUser(60e3, 135e3);
    eta_reco->GetYaxis()->SetRangeUser(60e3, 135e3);
    eta_reco_match->SetStats(0);
    eta_reco->Draw();
    eta_reco->SetStats(0);
    eta_reco_match->Draw("same");
    eta_reco_match->SetLineColor(kRed);
    TLegend* leg1 = new TLegend(0.7,0.9, 0.95, 1);
    leg1->AddEntry(eta_reco_match, "#eta reco muon matched","lp");
    leg1->AddEntry(eta_reco,"#eta reco muon","lp");
    leg1->Draw();
    c_eta_sovra->Print("eta_sovra_EmuReco.pdf");*/
     
    //TCanvas* c_pt_reco_match = new TCanvas("p_{T} reco match","#eta reco match",140,30,700,500);
    //pt_reco_match->Draw();

    //TCanvas* c_pt_reco = new TCanvas("#p_{T} reco","#eta reco",140,30,700,500);
    //pt_reco->Draw();
    
    /*TCanvas* c_pt_sovra = new TCanvas("p_{T} sovra","p_{T} sovra",140,30,700,500);
    pt_reco_match->SetStats(0);
    pt_reco->Draw();
    pt_reco->SetStats(0);
    pt_reco_match->Draw("same");
    pt_reco_match->SetLineColor(kRed);
    TLegend* leg2 = new TLegend(0.7,0.9, 0.95, 1);
    leg2->AddEntry(pt_reco_match, "p_{T} reco muon matched","lp");
    leg2->AddEntry(pt_reco,"p_{T} reco muon","lp");
    leg2->Draw();
    c_pt_sovra->Print("pt_sovra_EmuReco.pdf");*/
    
    //TCanvas* c_scatter_eta_pt = new TCanvas("scatter","scatter",140,30,700,500);
    //scatter_eta_pt->Draw();
/*    
    TCanvas* c_pt_sovra = new TCanvas("p_{T} sovra","p_{T} sovra",140,30,700,500);
    efficiency_pt->SetStats(0);
    efficiency_pt->Draw();
    efficiency_pt_Raw->SetStats(0);
    efficiency_pt_Raw->Draw("same");
    efficiency_pt_Raw->SetLineColor(kRed);
    TLatex n444;
  	n444.SetNDC();
  	n444.SetTextFont(52);
  	n444.SetTextSize(0.04);
    n444.DrawLatex(0.65, 0.55, "Cut: #DeltaR < 0.45");
    n444.DrawLatex(0.65, 0.5, "|#eta| < 0.8");
    TLegend* leg22 = new TLegend(0.7,0.9, 0.95, 1);
    leg22->AddEntry(efficiency_pt, "p_{T} efficiency  Emu Reco","lp");
    leg22->AddEntry(efficiency_pt_Raw,"p_{T} efficiency Raw Reco","lp");
    leg22->Draw();
    c_pt_sovra->Print("pt_sovra.pdf");
    
    TCanvas* c_eta_sovra = new TCanvas("#eta sovra","#eta sovra",140,30,700,500);
    efficiency_eta->SetStats(0);
    efficiency_eta->Draw();
    efficiency_eta_Raw->SetStats(0);
    efficiency_eta_Raw->Draw("same");
    efficiency_eta_Raw->SetLineColor(kRed);
    TLatex n333;
  	n333.SetNDC();
  	n333.SetTextFont(52);
  	n333.SetTextSize(0.04);
    n333.DrawLatex(0.55, 0.45, "Cut: #DeltaR < 0.45");
    n333.DrawLatex(0.55, 0.4, "p_{T} > 12");
    efficiency_eta->GetYaxis()->SetRangeUser(0.8, 1.1);
    efficiency_eta_Raw->GetYaxis()->SetRangeUser(0.8, 1.1);
    TLegend* leg222 = new TLegend(0.7,0.9, 0.95, 1);
    leg222->AddEntry(efficiency_eta, "#eta efficiency Emu Reco","lp");
    leg222->AddEntry(efficiency_eta_Raw,"#eta efficiency Raw Reco","lp");
    leg222->Draw();
    c_eta_sovra->Print("eta_sovra.pdf");

    
TCanvas* c_pt_reco_match = new TCanvas("#pt reco match","#pt reco match",140,30,700,500);
pt_Emu_match->SetStats(1);
pt_Emu_match->Draw();
c_pt_reco_match->Print("pt_emu_match.pdf");
TCanvas* c_pt_reco = new TCanvas("#pt reco","#pt reco",140,30,700,500);
pt_Emu->SetStats(1);
pt_Emu->Draw();
c_pt_reco->Print("pt_emu.pdf");

TCanvas* c_pt_reco_match_Raw = new TCanvas("#pt reco match Raw","#pt reco match Raw",140,30,700,500);
pt_Raw_match->SetStats(1);
pt_Raw_match->Draw();
c_pt_reco_match_Raw->Print("pt_raw_match.pdf");
TCanvas* c_pt_reco_Raw = new TCanvas("#pt reco  Raw","#pt reco  Raw",140,30,700,500);
pt_Raw->SetStats(1);
pt_Raw->Draw();
c_pt_reco_Raw->Print("pt_raw.pdf");

TCanvas* c_eta_reco_match = new TCanvas("#eta reco match","#eta reco match",140,30,700,500);
eta_Emu_match->SetStats(1);
eta_Emu_match->Draw();
c_eta_reco_match->Print("eta_emu_match.pdf");
TCanvas* c_eta_reco = new TCanvas("#eta reco","#eta reco",140,30,700,500);
eta_Emu->SetStats(1);
eta_Emu->Draw();
c_eta_reco->Print("eta_emu.pdf");

TCanvas* c_eta_reco_match_Raw = new TCanvas("#eta reco match Raw","#eta reco match Raw",140,30,700,500);
eta_Raw_match->SetStats(1);
eta_Raw_match->Draw();
c_eta_reco_match_Raw->Print("eta_raw_match.pdf");
TCanvas* c_eta_reco_Raw = new TCanvas("#eta reco  Raw","#eta reco  Raw",140,30,700,500);
eta_Raw->SetStats(1);
eta_Raw->Draw();
c_eta_reco_Raw->Print("eta_raw.pdf");
*/

    TCanvas* c_scatter_eta_pt = new TCanvas("Emu scatter","Emu scatter",140,30,700,500);
    scatter_eta_pt->Draw("colz");

//    TCanvas* c_Raw_scatter_eta_pt = new TCanvas("Raw scatter","Raw scatter",140,30,700,500);
//    Raw_scatter_eta_pt->Draw("colz");

//TCanvas* c_eta_raw_diff = new TCanvas("#eta raw diff ","#eta raw diff",140,30,700,500);
//TCanvas* c_eta_emu_diff = new TCanvas("#eta emu diff ","#eta emu diff",140,30,700,500);







/*c_eta_raw->Write();
c_eta_emu->Write();
c_eta_sovra->Write();
c_pt_sovra->Write();
c_eff_pt_Raw->Write();
c_eff_pt->Write();
c_eff_Raw->Write();
c_eff->Write();
*/



} //funzione principale








