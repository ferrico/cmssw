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
#include "TChain.h"

#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeDataFormat.h"


void Prova(){

	//TFile * file = new TFile("/lustre/cms/store/user/ferrico/ZeroBias1/reEmul_ZeroBias1_Run2015D_v1_RAW/160228_163327/0000/L1Ntuple_1.root");
        //TTree * treeL1Up  = (TTree*) file->Get("l1UpgradeEmuTree/L1UpgradeTree");
	//L1Analysis::L1AnalysisL1UpgradeDataFormat    *upgrade_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
        //treeL1Up->SetBranchAddress("L1Upgrade", &upgrade_);
	TChain* chain = new TChain("l1UpgradeEmuTree/L1UpgradeTree"); 	
	TH1F *Global_eta_GEM = new TH1F("Global_Eta_GEM","Global_Eta_GEM", 70, -3.5, 3.5);
	TString nome_file;
        TFile *file;
	
	for(int i = 1; i <10; i++){
	    nome_file = Form("/lustre/cms/store/user/ferrico/ZeroBias1/reEmul_ZeroBias1_Run2015D_v1_RAW/160228_163327/0000/L1Ntuple_%d.root", i);
	    file= new TFile(nome_file);
            //TTree * treeL1Up  = (TTree*) file->Get("l1UpgradeEmuTree/L1UpgradeTree"); 
            //L1Analysis::L1AnalysisL1UpgradeDataFormat    *upgrade_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();    
            //treeL1Up->SetBranchAddress("L1Upgrade", &upgrade_);	    
            chain->Add(nome_file);
	    //chain->SetBranchAddress("L1Upgrade", &upgrade_);
	    cout<<i<<endl;
 	    
	}



	cout<<"prova"<<endl;
	Int_t nevent = chain->GetEntries();
	cout<<nevent<<endl;

	for(int i = 0; i < nevent; i++){
	  chain->GetEntry(i);
            TTree * treeL1Up  = (TTree*) file->Get("l1UpgradeEmuTree/L1UpgradeTree"); 
            L1Analysis::L1AnalysisL1UpgradeDataFormat    *upgrade_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();    
            treeL1Up->SetBranchAddress("L1Upgrade", &upgrade_);           
            chain->SetBranchAddress("L1Upgrade", &upgrade_);
            cout<<" ---- "<<upgrade_->muonEta.size()<<endl; 
	  for(int j = 0; j < upgrade_->muonEta.size(); j++){
		Global_eta_GEM->Fill(upgrade_->muonEta.at(j));
	  }
	}


	    TCanvas* cGlobal_Eta_GEM_prima = new TCanvas("cGlobal_Eta_GEM_prima","cGlobal_Eta_GEM_prima",140,30,700,500);
	        Global_eta_GEM->Draw();
		//chain->Draw("muonEta");



}
