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


void Efficienza_EmuReco_L()
{
    //gROOT->Reset();
    //gROOT->SetBatch();

    TFile* f2 = new TFile("/cmshome/filippoe/trigger/CMSSW_8_0_0_pre6/src/L1Ntuple_273450_SingleMuon.root","READ");
    TFile* f = new TFile("/cmshome/filippoe/trigger/CMSSW_8_0_0_pre6/src/L1_ZMu_RawReco_v3.root","READ");
    
    TFile* Efficiency = new TFile("./Efficiency.root", "recreate");
    Efficiency->cd();

	
    TTree* t_Reco_Emu = (TTree*) f2->Get("l1MuonRecoTree/Muon2RecoTree");
    L1Analysis::L1AnalysisRecoMuon2DataFormat    *RecoMuonEmu_ = new L1Analysis::L1AnalysisRecoMuon2DataFormat();
    t_Reco_Emu->SetBranchAddress("Muon", &RecoMuonEmu_);

    TTree* t_Upgrade = (TTree*) f2->Get("l1UpgradeTree/L1UpgradeTree");
    L1Analysis::L1AnalysisL1UpgradeDataFormat    *upgrade_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
    t_Upgrade->SetBranchAddress("L1Upgrade", &upgrade_);
    
/*    TTree* t_Extra = (TTree*) f->Get("l1ExtraTree/L1ExtraTree");
    L1Analysis::L1AnalysisL1ExtraDataFormat    *extra_ = new L1Analysis::L1AnalysisL1ExtraDataFormat();
    t_Extra->SetBranchAddress("L1Extra", &extra_);
    
    TTree* t_Reco_Raw = (TTree*) f->Get("l1MuonRecoTree/Muon2RecoTree");
    L1Analysis::L1AnalysisRecoMuon2DataFormat    *RecoMuonRaw_ = new L1Analysis::L1AnalysisRecoMuon2DataFormat();
    t_Reco_Raw->SetBranchAddress("Muon", &RecoMuonRaw_);
  */  
    
    Double_t PT_BINS[] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85};
    Int_t  binnum = sizeof(PT_BINS)/sizeof(Double_t)-1;
    
     Double_t PT_BINS_1[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    Int_t  binnum_1 = sizeof(PT_BINS_1)/sizeof(Double_t)-1;

    
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
    
    
    TH1F *ratio_efficiency_pt =  new TH1F("ratio p_{T} efficiency","ratio p_{T} efficiency", binnum, PT_BINS);
    ratio_efficiency_pt->GetXaxis()->SetTitle("p_{T}");
  
  
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
    
    
    
    // TH1 for muon in the dump: pt < 10 
    TH1F *eta_Raw_A = new TH1F("#eta Reco (Raw)","#eta Reco (Raw)", 40,-1, 1);
    eta_Raw_A->GetXaxis()->SetTitle("#eta");
    
    TH1F *phi_Raw_A = new TH1F("#phi Reco (Raw)","#phi Reco (Raw)", 128,-3.2, 3.2);
    phi_Raw_A->GetXaxis()->SetTitle("#phi");
    
    TH1F *pt_Raw_A = new TH1F("p_{T} Reco (Raw)","p_{T} Reco (Raw)", binnum_1, PT_BINS_1);
    pt_Raw_A->GetXaxis()->SetTitle("p_{T}");
    
    
    
    
    TH1F *eta_Emu_A = new TH1F("#eta Reco (Emu)","#eta Reco (Emu)", 40,-1, 1);
    eta_Emu_A->GetXaxis()->SetTitle("#eta");
    
    TH1F *pt_Emu_A = new TH1F("p_{T} Reco (Emu)","p_{T} Reco (Emu)", binnum_1, PT_BINS_1);
    pt_Emu_A->GetXaxis()->SetTitle("p_{T}");
    
    TH1F *phi_Emu_A = new TH1F("#phi Reco (Emu)","#phi Reco (Emu)", 128,-3.2, 3.2);
    phi_Emu_A->GetXaxis()->SetTitle("#phi");
    
    
  	TH2F *scatter_pt_Raw = new TH2F("Raw vs Reco","Raw vs Reco", 10, 0, 10, 30, 10, 40);
    scatter_pt_Raw->GetXaxis()->SetTitle("Reco p_{T}");
    scatter_pt_Raw->GetYaxis()->SetTitle("Raw p_{T}");
    scatter_pt_Raw->SetMarkerSize(1);

  	TH2F *scatter_pt_Emu = new TH2F("Emu","Emu vs Reco", 10, 0, 10, 30, 10, 40);
    scatter_pt_Emu->GetXaxis()->SetTitle("Reco p_{T}");
    scatter_pt_Emu->GetYaxis()->SetTitle("Emu p_{T}");
    scatter_pt_Emu->SetMarkerSize(1);

    TH1F *phi_Raw_eta_cut= new TH1F("#phi Reco (Raw) [|#eta| > 0.25]","#phi Reco (Raw) [|#eta| > 0.25]", 128,-3.2, 3.2);
    phi_Raw_A->GetXaxis()->SetTitle("#phi");
    
    TH1F *phi_Emu_eta_cut = new TH1F("#phi Reco (Emu) [|#eta| > 0.25]","#phi Reco (Emu) [|#eta| > 0.25]", 128,-3.2, 3.2);
    phi_Emu_A->GetXaxis()->SetTitle("#phi");
  
    int j;
    int i;
    int h;
    int z;
    float deltaR;
    float min;
    float emu_eta;
    float raw_eta;
    float emu_pt;
    float raw_pt;
    float deltaR_2;
    float min_2;

    TVector3 trigger_muon;
    TVector3 reco_muon;
    TVector3 Raw_muon;
    TVector3 reco_muon_Raw;
    TVector3 reco;


    // Raw - Reco  
    float extra_pt; 
/*    
	Long64_t nentries_Raw = t_Extra->GetEntriesFast();
	cout<<nentries_Raw<<endl;
	
	for (Long64_t jentry = 0; jentry<nentries_Raw; jentry++){
    	
    	if(jentry%100000 == 0) cout<<jentry<<endl;
        
        t_Reco_Raw->GetEntry(jentry);
        t_Extra->GetEntry(jentry);
        
        for(h = 0; h < RecoMuonRaw_->pt.size(); h++){
        	if(abs(RecoMuonRaw_->eta.at(h)) > 0.8) continue;
        	if(!(RecoMuonRaw_->isMediumMuon.at(h))) continue;
        	if(RecoMuonRaw_->iso.at(h) < 0.3) continue;

			min = 100;
        	min_2 = 100;
        	reco_muon_Raw.SetPtEtaPhi(RecoMuonRaw_->pt.at(h), RecoMuonRaw_->eta.at(h), RecoMuonRaw_->phi.at(h));
        	for(z = h+1; z < RecoMuonRaw_->pt.size(); z++){
        		reco.SetPtEtaPhi(RecoMuonRaw_->pt.at(z), RecoMuonRaw_->eta.at(z), RecoMuonRaw_->phi.at(z));
        		deltaR_2 = reco.DeltaR(reco_muon_Raw);
                if(deltaR_2 < min_2){
                    min_2 = deltaR_2;
                }
        	}
        	if(min_2 < 0.45) continue;
        	pt_reco_Raw->Fill(RecoMuonRaw_->pt.at(h));
        	for(j = 0; j < extra_->muonEt.size(); j++){
        		if(extra_->muonEt.at(j) < 15) continue;        		
        		if((extra_->muonQuality.at(j) != 6) && (extra_->muonQuality.at(j) != 7)) continue;
				Raw_muon.SetPtEtaPhi(extra_->muonEt.at(j), extra_->muonEta.at(j), extra_->muonPhi.at(j));
        		deltaR = Raw_muon.DeltaR(reco_muon_Raw);
                if(deltaR < min){
                    min = deltaR;
                    extra_pt = extra_->muonEt.at(j);
                }
            }
        	if(min < 0.45){
	        	pt_reco_match_Raw->Fill(RecoMuonRaw_->pt.at(h));
	        	if((extra_pt > 12) && (RecoMuonRaw_->pt.at(h) < 8))
					scatter_pt_Raw->Fill(RecoMuonRaw_->pt.at(h), extra_pt);
	        	if(RecoMuonRaw_->pt.at(h) < 12){
    	    		pt_Raw_A->Fill(RecoMuonRaw_->pt.at(h));
        			eta_Raw_A->Fill(RecoMuonRaw_->eta.at(h));
        			phi_Raw_A->Fill(RecoMuonRaw_->phi.at(h));
        			if(abs(RecoMuonRaw_->eta.at(h)) > 0.25)
        				phi_Raw_eta_cut->Fill(RecoMuonRaw_->phi.at(h));
	        	}
            }
         }    
    }	
    
	for (Long64_t jentry = 0; jentry<nentries_Raw; jentry++){
    	
    	if(jentry%100000 == 0) cout<<jentry<<endl;
        
        t_Reco_Raw->GetEntry(jentry);
        t_Extra->GetEntry(jentry);
        
        
        for(h = 0; h < RecoMuonRaw_->eta.size(); h++){
        	if(RecoMuonRaw_->pt.at(h) < 12) continue;
        	if(!(RecoMuonRaw_->isMediumMuon.at(h))) continue;
        	if(RecoMuonRaw_->iso.at(h) < 0.3) continue;
        	min = 100;
        	min_2 = 100;
        	reco_muon_Raw.SetPtEtaPhi(RecoMuonRaw_->pt.at(h), RecoMuonRaw_->eta.at(h), RecoMuonRaw_->phi.at(h));
        	for(z = h+1; z < RecoMuonRaw_->pt.size(); z++){
        		reco.SetPtEtaPhi(RecoMuonRaw_->pt.at(z), RecoMuonRaw_->eta.at(z), RecoMuonRaw_->phi.at(z));
        		deltaR_2 = reco.DeltaR(reco_muon_Raw);
                if(deltaR_2 < min_2){
                    min_2 = deltaR_2;
                }
        	}
        	if(min_2 < 0.45) continue;
        	eta_reco_Raw->Fill(RecoMuonRaw_->eta.at(h));
        	for(j = 0; j < extra_->muonEta.size(); j++){
        		if((extra_->muonQuality.at(j) != 6) && (extra_->muonQuality.at(j) != 7)) continue;
	            Raw_muon.SetPtEtaPhi(extra_->muonEt.at(j), extra_->muonEta.at(j), extra_->muonPhi.at(j));
        		deltaR = Raw_muon.DeltaR(reco_muon_Raw);
                if(deltaR < min){
                    min = deltaR;
                }
            }
        	if(min < 0.45){
	        	eta_reco_match_Raw->Fill(RecoMuonRaw_->eta.at(h));
            }
         }
         		   		
         
    }
*/		
    // Emu - Reco
	float upgrade_pt;
	
    Long64_t nentries = t_Upgrade->GetEntriesFast();
	cout<<nentries<<endl;

    for (Long64_t jentry = 0; jentry<nentries;jentry++){
    	
    	if(jentry%100000 == 0)
    		cout<<jentry<<endl;
        
        t_Reco_Emu->GetEntry(jentry);
        t_Upgrade->GetEntry(jentry);
        
        for(h = 0; h < RecoMuonEmu_->pt.size(); h++){
        	if(abs(RecoMuonEmu_->eta.at(h)) > 0.8) continue;
        	if(!(RecoMuonEmu_->isMediumMuon.at(h))) continue;
        	if(RecoMuonEmu_->iso.at(h) < 0.3) continue;
        	min = 100;
        	min_2 = 100;
        	reco_muon.SetPtEtaPhi(RecoMuonEmu_->pt.at(h), RecoMuonEmu_->eta.at(h), RecoMuonEmu_->phi.at(h));
        	for(z = h+1; z < RecoMuonEmu_->pt.size(); z++){
        		reco.SetPtEtaPhi(RecoMuonEmu_->pt.at(z), RecoMuonEmu_->eta.at(z), RecoMuonEmu_->phi.at(z));
        		deltaR_2 = reco.DeltaR(reco_muon);
                if(deltaR_2 < min_2){
                    min_2 = deltaR_2;
                }
        	}
        	if(min_2 < 0.45) continue;        
        	pt_reco->Fill(RecoMuonEmu_->pt.at(h));
	        for(j = 0; j < upgrade_->muonEt.size(); j++){
	        	//if(upgrade_->muonEt.at(j) < 15) continue;
        		//if((upgrade_->muonQual.at(j) != 6) && (upgrade_->muonQual.at(j) != 7)) continue;
	            trigger_muon.SetPtEtaPhi(upgrade_->muonEt.at(j), upgrade_->muonEta.at(j), upgrade_->muonPhi.at(j));
        		deltaR = trigger_muon.DeltaR(reco_muon);
                if(deltaR < min){
                    min = deltaR;
                    upgrade_pt = upgrade_->muonEt.at(j);
                }
            }
           
        	if(min < 0.45){
	        	pt_reco_match->Fill(RecoMuonEmu_->pt.at(h));
	        	if((upgrade_pt > 12) && (RecoMuonEmu_->pt.at(h) < 8))
		       		scatter_pt_Emu->Fill(RecoMuonEmu_->pt.at(h), upgrade_pt);
        		if(RecoMuonEmu_->pt.at(h) < 12){
        			pt_Emu_A->Fill(RecoMuonEmu_->pt.at(h));
	        		eta_Emu_A->Fill(RecoMuonEmu_->eta.at(h));
    	    		phi_Emu_A->Fill(RecoMuonEmu_->phi.at(h));
        			if(abs(RecoMuonEmu_->eta.at(h)) > 0.25)
        				phi_Emu_eta_cut->Fill(RecoMuonEmu_->phi.at(h));       	
        		}
            }
         }
    }
 
    for (Long64_t jentry = 0; jentry<nentries;jentry++){
    	
    	if(jentry%100000 == 0)
    		cout<<jentry<<endl;
        
        t_Reco_Emu->GetEntry(jentry);
        t_Upgrade->GetEntry(jentry);
        
        for(h = 0; h < RecoMuonEmu_->eta.size(); h++){
        	if(RecoMuonEmu_->pt.at(h) < 12) continue;
        	if(!(RecoMuonEmu_->isMediumMuon.at(h))) continue;
        	if(RecoMuonEmu_->iso.at(h) < 0.3) continue;
        	min = 100;
        	min_2 = 100;
        	reco_muon.SetPtEtaPhi(RecoMuonEmu_->pt.at(h), RecoMuonEmu_->eta.at(h), RecoMuonEmu_->phi.at(h));
        	for(z = h+1; z < RecoMuonEmu_->pt.size(); z++){
        		reco.SetPtEtaPhi(RecoMuonEmu_->pt.at(z), RecoMuonEmu_->eta.at(z), RecoMuonEmu_->phi.at(z));
        		deltaR_2 = reco.DeltaR(reco_muon);
                if(deltaR_2 < min_2){
                    min_2 = deltaR_2;
                }
        	}
        	if(min_2 < 0.45) continue; 
        	eta_reco->Fill(RecoMuonEmu_->eta.at(h));
	        for(j = 0; j < upgrade_->muonEta.size(); j++){
		        trigger_muon.SetPtEtaPhi(upgrade_->muonEt.at(j), upgrade_->muonEta.at(j), upgrade_->muonPhi.at(j));
        		deltaR = trigger_muon.DeltaR(reco_muon);
                if(deltaR < min){
                    min = deltaR;
                }
            }
        	if(min < 0.45){
	        	eta_reco_match->Fill(RecoMuonEmu_->eta.at(h));
            }
         }
    }




 


  
    

 
 
/*  
    TCanvas* c_eta_sovra_raw = new TCanvas("#eta sovra raw","#eta sovra raw",140,30,700,500);
    eta_reco_match_Raw->GetXaxis()->SetRangeUser(-3, 3);
    eta_reco_Raw->GetXaxis()->SetRangeUser(-3, 3);
    eta_reco_match_Raw->SetStats(0);
    eta_reco_Raw->Draw();
    eta_reco_Raw->SetStats(0);
    eta_reco_match_Raw->Draw("same");
    eta_reco_match_Raw->SetLineColor(kRed);
    TLegend* leg11 = new TLegend(0.65,0.6, 0.95, 0.8);
    leg11->AddEntry(eta_reco_match_Raw, "#eta reco muon matched - Stage1","lp");
    leg11->AddEntry(eta_reco_Raw,"#eta reco muon - Stage1","lp");
    leg11->Draw();
    c_eta_sovra_raw->Print("eta_sovra_RawReco.pdf");
     

    
    TCanvas* c_pt_sovra_raw = new TCanvas("p_{T} sovra raw ","p_{T} sovra raw",140,30,700,500);
    pt_reco_match_Raw->SetStats(0);
    pt_reco_match_Raw->SetTitle(0);
    pt_reco_Raw->Draw();
    pt_reco_Raw->SetStats(0);
    pt_reco_Raw->SetTitle(0);
    pt_reco_match_Raw->Draw("same");
    pt_reco_match_Raw->SetLineColor(kRed);
    TLegend* leg22 = new TLegend(0.65,0.6, 0.95, 0.8);
    leg22->AddEntry(pt_reco_match_Raw, "p_{T} reco muon matched - Stage1","lp");
    leg22->AddEntry(pt_reco_Raw,"p_{T} reco muon - Stage1","lp");
    leg22->Draw();
    c_pt_sovra_raw->Print("pt_sovra_RawReco.pdf");
    

    
    TCanvas* c_eta_sovra = new TCanvas("#eta sovra","#eta sovra",140,30,700,500);
    eta_reco_match->GetXaxis()->SetRangeUser(-3, 3);
    eta_reco->GetXaxis()->SetRangeUser(-3, 3);
    eta_reco_match->SetStats(0);
    eta_reco->Draw();
    eta_reco->SetStats(0);
    eta_reco_match->Draw("same");
    eta_reco_match->SetLineColor(kRed);
    TLegend* leg1 = new TLegend(0.65,0.6, 0.95, 0.8);
    leg1->AddEntry(eta_reco_match, "#eta reco muon matched - Stage2","lp");
    leg1->AddEntry(eta_reco,"#eta reco muon - Stage2","lp");
    leg1->Draw();
    c_eta_sovra->Print("eta_sovra_EmuReco.pdf");
     

    
    TCanvas* c_pt_sovra = new TCanvas("p_{T} sovra","p_{T} sovra",140,30,700,500);
    pt_reco_match->SetStats(0);
    pt_reco->Draw();
    pt_reco->SetStats(0);
    pt_reco_match->Draw("same");
    pt_reco_match->SetLineColor(kRed);
    TLegend* leg2 = new TLegend(0.65,0.6, 0.95, 0.8);
    leg2->AddEntry(pt_reco_match, "p_{T} reco muon matched - Stage2","lp");
    leg2->AddEntry(pt_reco,"p_{T} reco muon - Stage2","lp");
    leg2->Draw();
    c_pt_sovra->Print("pt_sovra_EmuReco.pdf");
 */   
    eta_reco->Sumw2(kTRUE);
    eta_reco_Raw->Sumw2(kTRUE);
    eta_reco_match->Sumw2(kTRUE);
    eta_reco_match_Raw->Sumw2(kTRUE);
    efficiency_eta->Divide(eta_reco_match, eta_reco,1,1,"B");
    efficiency_eta_Raw->Divide(eta_reco_match_Raw, eta_reco_Raw,1,1,"B");
    
    TCanvas* c_eff_eta = new TCanvas("#eta efficiency Emu Reco","#eta efficiency Emu Reco",140,30,700,500);
    efficiency_eta->SetStats(0);
    efficiency_eta->GetYaxis()->SetRangeUser(0.7, 1.05);
    efficiency_eta->GetXaxis()->SetRangeUser(-1, 1);
    efficiency_eta->Draw();
    TLatex n4;
  	n4.SetNDC();
  	n4.SetTextFont(52);
  	n4.SetTextSize(0.04);
    n4.DrawLatex(0.65, 0.55, "Cut: #DeltaR < 0.45");
    n4.DrawLatex(0.65, 0.5, "p_{T} > 12");
    c_eff_eta->Print("eta_efficiency_EmuReco.pdf");
    
    TCanvas* c_eff_eta_Raw = new TCanvas("#eta efficiency Raw Reco","eta efficiency Raw Reco",140,30,700,500);
    efficiency_eta_Raw->SetStats(0);
    efficiency_eta_Raw->Draw();
    TLatex n44;
  	n44.SetNDC();
  	n44.SetTextFont(52);
  	n44.SetTextSize(0.04);
    n44.DrawLatex(0.65, 0.55, "Cut: #DeltaR < 0.45");
    n44.DrawLatex(0.65, 0.5, "p_{T} > 12");
    //c_eff_eta_Raw->Print("eta_efficiency_RawReco.pdf");


    TCanvas* c_eff_eta_sovra = new TCanvas("#eta efficiency","#eta efficiency",140,30,700,500);
    efficiency_eta_Raw->SetStats(0);
    efficiency_eta->SetStats(0);
    efficiency_eta->Draw();
    efficiency_eta_Raw->Draw("same");    
    efficiency_eta_Raw->SetLineColor(kRed);
    efficiency_eta->GetXaxis()->SetRangeUser(-1, 1);
    efficiency_eta_Raw->GetXaxis()->SetRangeUser(-1, 1);
    efficiency_eta->GetYaxis()->SetRangeUser(0.7, 1.1);
    efficiency_eta_Raw->GetYaxis()->SetRangeUser(0.7, 1.1);
    TLegend* leg222 = new TLegend(0.7,0.9, 0.95, 1);
    leg222->AddEntry(efficiency_eta, "#eta efficiency - Stage2","lp");
    leg222->AddEntry(efficiency_eta_Raw,"#eta efficiency - Stage1","lp");
    leg222->Draw();
    n44.DrawLatex(0.65, 0.35, "Cut: #DeltaR < 0.45");
    n44.DrawLatex(0.65, 0.3, "p_{T} > 12");
    //c_eff_eta_sovra->Print("eta_eff_sovra.pdf");
    
    

    pt_reco->Sumw2(kTRUE);
    pt_reco_Raw->Sumw2(kTRUE);
    pt_reco_match->Sumw2(kTRUE);
    pt_reco_match_Raw->Sumw2(kTRUE);
    efficiency_pt->Divide(pt_reco_match, pt_reco,1,1,"B");
    efficiency_pt_Raw->Divide(pt_reco_match_Raw, pt_reco_Raw,1,1,"B");
 
    TCanvas* c_eff_pt = new TCanvas("p_{T} efficiency Emu Reco","p_{T} efficiency Emu Reco",140,30,700,500);
    efficiency_pt->SetStats(0);
    efficiency_pt->Draw();
    TLatex n14;
  	n14.SetNDC();
  	n14.SetTextFont(52);
  	n14.SetTextSize(0.04);
    n14.DrawLatex(0.65, 0.55, "Cut: #DeltaR < 0.45");
    n14.DrawLatex(0.65, 0.5, "|#eta| < 0.8");
    c_eff_pt->Print("pt_efficiency_EmuReco.pdf");
    
    TCanvas* c_eff_pt_Raw = new TCanvas("p_{T} efficiency Raw Reco","p_{T} efficiency Raw Reco",140,30,700,500);
    efficiency_pt_Raw->SetStats(0);
    efficiency_pt_Raw->Draw();
    TLatex n144;
  	n144.SetNDC();
  	n144.SetTextFont(52);
  	n144.SetTextSize(0.04);
    n144.DrawLatex(0.65, 0.55, "Cut: #DeltaR < 0.45");
    n144.DrawLatex(0.65, 0.5, "|#eta| < 0.8");
    //c_eff_pt_Raw->Print("pt_efficiency_RawReco.pdf");


    TCanvas* c_eff_pt_sovra = new TCanvas("p_{T} efficiency","p_{T} efficiency",140,30,700,500);
    efficiency_pt_Raw->SetStats(0);
    efficiency_pt->SetStats(0);
    efficiency_pt->Draw();
    efficiency_pt_Raw->Draw("same");    
    efficiency_pt_Raw->SetLineColor(kRed);
    TLegend* leg1222 = new TLegend(0.65, 0.9, 0.95, 1.);
    leg1222->AddEntry(efficiency_pt, "p_{T} efficiency - Stage2","lp");
    leg1222->AddEntry(efficiency_pt_Raw,"p_{T} efficiency - Stage1","lp");
    leg1222->Draw();
    n14.DrawLatex(0.65, 0.55, "Cut: #DeltaR < 0.45");
    n14.DrawLatex(0.65, 0.5, "|#eta| < 0.8");
    //c_eff_pt_sovra->Print("pt_eff_sovra.pdf");
    
 /*   
    ratio_efficiency_pt->Divide(efficiency_pt, efficiency_pt_Raw,1,1,"B");
    TCanvas* c_eff_pt_ratio = new TCanvas("ratio p_{T} efficiency","ratio p_{T} efficiency",140,30,700,500);
    ratio_efficiency_pt->GetXaxis()->SetRangeUser(12, 85);
    ratio_efficiency_pt->SetStats(0);
    ratio_efficiency_pt->Draw();
    TLatex n1144;
  	n1144.SetNDC();
  	n1144.SetTextFont(52);
  	n1144.SetTextSize(0.04);
    n1144.DrawLatex(0.55, 0.35, "Efficiency ratio: 2016 Stage2 / 2015 Stage1");
    TLine *line = new TLine(12.,1,85,1);
    line->SetLineColor(kRed);
    line->SetLineWidth(3);
    line->Draw();
    c_eff_pt_ratio->Print("pt_efficiency_ratio.pdf");



    TCanvas* c_eta_raw_A = new TCanvas("#eta raw A","#eta raw A",140,30,700,500);
    eta_Raw_A->Draw();
    c_eta_raw_A->Print("eta_RawReco_A.pdf");
    
    TCanvas* c_pt_raw_A = new TCanvas("#pt raw A","#pt raw A",140,30,700,500);
    pt_Raw_A->Draw();
    c_pt_raw_A->Print("pt_RawReco_A.pdf");
    
    TCanvas* c_phi_raw_A = new TCanvas("#phi raw A","#phi raw A",140,30,700,500);
    phi_Raw_A->Draw();
    c_phi_raw_A->Print("phi_RawReco_A.pdf");
    
    
    
    TCanvas* c_eta_Emu_A = new TCanvas("#eta Emu A","#eta Emu A",140,30,700,500);
    eta_Emu_A->Draw();
    c_eta_Emu_A->Print("eta_EmuReco_A.pdf");
    
    TCanvas* c_pt_Emu_A = new TCanvas("#pt Emu A","#pt Emu A",140,30,700,500);
    pt_Emu_A->Draw();
    c_pt_Emu_A->Print("pt_EmuReco_A.pdf");
    
    TCanvas* c_phi_Emu_A = new TCanvas("#phi Emu A","#phi Emu A",140,30,700,500);
    phi_Emu_A->Draw();
    c_phi_Emu_A->Print("phi_EmuReco_A.pdf");
    
    TCanvas* c_phi_raw_eta_cut = new TCanvas("#phi raw eta_cut","#phi raw eta_cut",140,30,700,500);
    phi_Raw_eta_cut->Draw();
    c_phi_raw_eta_cut->Print("phi_RawReco_eta_cut.pdf");

    TCanvas* c_phi_Emu_eta_cut= new TCanvas("#phi Emu eta_cut","#phi Emu eta_cut",140,30,700,500);
    phi_Emu_eta_cut->Draw();
    c_phi_Emu_eta_cut->Print("phi_EmuReco_eta_cut.pdf");

	gStyle->SetOptStat(0);    
    TCanvas* c_pt_scatter_Raw = new TCanvas("Scatter Raw vs Reco","Scatter Raw vs Reco",140,30,700,500);
	scatter_pt_Raw->Draw();
	c_pt_scatter_Raw->Print("pt_scatter_Raw.pdf");

	TCanvas* c_pt_scatter_Emu = new TCanvas("Scatter Emu vs Reco", "Scatter Emu vs Reco",140,30,700,500);
	scatter_pt_Emu->Draw();
	c_pt_scatter_Emu->Print("pt_scatter_Emu.pdf");
	
	

    
    //eta_reco_match_Raw->GetXaxis()->SetRangeUser(-3, 3);
    //eta_reco_Raw->GetXaxis()->SetRangeUser(-3, 3);
    //eta_reco_match_Raw->SetStats(0);
    //eta_reco_Raw->Draw();
    //eta_reco_Raw->SetStats(0);
    //eta_reco_match_Raw->Draw("same");
    //eta_reco_match_Raw->SetLineColor(kRed);
    //TLegend* leg11 = new TLegend(0.65,0.6, 0.95, 0.8);
    //leg11->AddEntry(eta_reco_match_Raw, "#eta reco muon matched - Stage1","lp");
    //leg11->AddEntry(eta_reco_Raw,"#eta reco muon - Stage1","lp");
    //leg11->Draw();



c_eta_raw->Write();
c_eta_emu->Write();
c_eta_sovra->Write();
c_pt_sovra->Write();
c_eff_pt_Raw->Write();
c_eff_pt->Write();
c_eff_Raw->Write();
c_eff->Write();
*/



} //funzione principale








