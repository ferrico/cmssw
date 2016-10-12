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

void rates_2(){
  nevents=1000;

  // make trees
  TFile * file = new TFile("/lustre/cms/store/user/ferrico/ZeroBias3/reEmul_ZeroBias3_Run2015D_v1_RAW/160228_163611/ZeroBias3.root");

  TTree * treeL1Up  = (TTree*) file->Get("l1UpgradeEmuTree/L1UpgradeTree");
  treeL1Up->Print();
  L1Analysis::L1AnalysisL1UpgradeDataFormat    *upgrade_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  treeL1Up->SetBranchAddress("L1Upgrade", &upgrade_);

  TTree * treeL1Up_Raw  = (TTree*) file->Get("l1UpgradeTree/L1UpgradeTree");
  treeL1Up_Raw->Print();
  L1Analysis::L1AnalysisL1UpgradeDataFormat    *upgrade_raw_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  treeL1Up_Raw->SetBranchAddress("L1Upgrade", &upgrade_raw_);




    TFile *file_out = new TFile("./Rates_ZeroBias3/ZeroBia.root","recreate");
    file_out->cd();



  // mu bins
  int nMuBins = 50;
  float muLo = 10.;
  float muHi = 50.;
  float muBinWidth = (muHi-muLo)/nMuBins;

  // eg bins
  int nEgBins = 50;
  float egLo = 10.;
  float egHi = 260.;
  float egBinWidth = (egHi-egLo)/nEgBins;

  // tau bins
  int nTauBins = 50;
  float tauLo = 10.;
  float tauHi = 260.;
  float tauBinWidth = (tauHi-tauLo)/nTauBins;
 
  // jet bins
  int nJetBins = 50;
  float jetLo = 10.;
  float jetHi = 260.;
  float jetBinWidth = (jetHi-jetLo)/nJetBins;

  // etSum bins
  int nEtSumBins = 100;
  float etSumLo = 10.;
  float etSumHi = 510.;
  float etSumBinWidth = (etSumHi-etSumLo)/nEtSumBins;

  // htSum bins
  int nHtSumBins = 100;
  float htSumLo = 10.;
  float htSumHi = 510.;
  float htSumBinWidth = (htSumHi-htSumLo)/nHtSumBins;

  // metSum bins
  int nMetSumBins = 50;
  float metSumLo = 10.;
  float metSumHi = 260.;
  float metSumBinWidth = (metSumHi-metSumLo)/nMetSumBins;

  // mhtSum bins
  int nMhtSumBins = 50;
  float mhtSumLo = 10.;
  float mhtSumHi = 260.;
  float mhtSumBinWidth = (mhtSumHi-mhtSumLo)/nMhtSumBins;

  //make histos
  TH1D* muRates = new TH1D("muRates", "", nMuBins, muLo-2.5, muHi-2.5);
  TH1D* egRates = new TH1D("egRates", "", nEgBins, egLo-2.5, egHi-2.5);
  TH1D* tauRates = new TH1D("tauRates", "", nTauBins, tauLo-2.5, tauHi-2.5);
  TH1D* jetRates = new TH1D("jetRates", "", nJetBins, jetLo-2.5, jetHi-2.5);
  TH1D* etSumRates = new TH1D("etSumRates","", nEtSumBins, etSumLo-2.5, etSumHi-2.5);
  TH1D* htSumRates = new TH1D("htSumRates","", nHtSumBins, htSumLo-2.5, htSumHi-2.5);
  TH1D* metSumRates = new TH1D("metSumRates","", nMetSumBins, metSumLo-2.5, metSumHi-2.5);
  TH1D* mhtSumRates = new TH1D("mhtSumRates","", nMhtSumBins, mhtSumLo-2.5, mhtSumHi-2.5);

  TH1D* hMuEt = new TH1D("muEt", "", 400, 0.,400.);
  TH1D* hEgEt = new TH1D("egEt", "", 400, 0.,400.);
  TH1D* hTauEt = new TH1D("tauEt","",400, 0.,400.);
  TH1D* hJetEt = new TH1D("jetEt","",400, 0.,400.);
  TH1D* hEtSum = new TH1D("etSum","",800, 0.,800.);
  TH1D* hHtSum = new TH1D("htSum","",800, 0.,800.);

  // get entries
  Long64_t nentries = treeL1Up->GetEntriesFast();


  for (Long64_t jentry=0; jentry<nentries; jentry++){
    if((jentry%100000)==0) std::cout << "Done " << jentry  << " events..." << std::endl;
    treeL1Up->GetEntry(jentry);
    
    //cout << upgrade_->nJets << "\n";

    // get Mu rates
    double muEt(0);
    for(uint it=0; it<upgrade_->nMuons; ++it){
      // work around a muon bug:
      int offset = upgrade_->muonQual.size() - upgrade_->nMuons;
      //cout << "INFO:  " << upgrade_->nMuons << "\n";
      //cout << "INFO:  " << upgrade_->muonEt.size() << "\n";
      //cout << "INFO:  " << upgrade_->muonQual.size() << "\n";
      if (upgrade_->muonQual[it+offset]==0) continue;
      if(abs(upgrade_->muonEta[it]) > 0.8) continue;
      hMuEt->Fill(upgrade_->muonEt[it]);
      //std::cout<<muEt<<" ----- "<<upgrade_->muonEt[it]<<endl;
      muEt = upgrade_->muonEt[it] > muEt ?  upgrade_->muonEt[it]  : muEt;
    }
    for(int bin=0; bin<nMuBins; bin++)
      if(muEt >= muLo + (bin*muBinWidth) ) muRates->Fill(muLo+(bin*muBinWidth)); //GeV

    // get Eg rates
    int egEt(0);
    for(uint it=0; it<upgrade_->nEGs; ++it){
      hEgEt->Fill(0.5*upgrade_->egIEt[it]);
      egEt = upgrade_->egIEt[it] > egEt ?  upgrade_->egIEt[it]  : egEt;
    }
    for(int bin=0; bin<nEgBins; bin++)
      if(egEt*0.5 >= egLo + (bin*egBinWidth) ) egRates->Fill(egLo+(bin*egBinWidth)); //GeV
    
    // get Tau rates
    int tauEt(0);
    for(uint it=0; it<upgrade_->nTaus; ++it){
      hTauEt->Fill(0.5*upgrade_->tauIEt[it]);
      tauEt = upgrade_->tauIEt[it] > tauEt ? upgrade_->tauIEt[it] : tauEt;
    }
    for(int bin=0; bin<nTauBins; bin++)
      if( (tauEt*0.5) >= tauLo + (bin*tauBinWidth) ) tauRates->Fill(tauLo+(bin*tauBinWidth)); //GeV
        
    // get Jet rates
    int jetEt(0);
    for(uint it=0; it<upgrade_->nJets; ++it){
      hJetEt->Fill(0.5*upgrade_->jetIEt[it]);
      jetEt =  upgrade_->jetIEt[it] > jetEt ? upgrade_->jetIEt[it] : jetEt;
    }
    for(int bin=0; bin<nJetBins; bin++)
      if( (jetEt*0.5) >= jetLo + (bin*jetBinWidth) ) jetRates->Fill(jetLo+(bin*jetBinWidth));  //GeV

    double etSum  = -1.0;
    double htSum  = -1.0;
    double metSum = -1.0;
    double mhtSum = -1.0;
    for(uint it=0; it<upgrade_->nSums; ++it){
      double et = upgrade_->sumEt[it];
      if (upgrade_->sumType[it] == L1Analysis::kTotalEt)   etSum  = et;
      if (upgrade_->sumType[it] == L1Analysis::kTotalHt)   htSum  = et;
      if (upgrade_->sumType[it] == L1Analysis::kMissingEt) metSum = et;
      if (upgrade_->sumType[it] == L1Analysis::kMissingHt) mhtSum = et;
    }
    //std::cout << "mht:  " << mhtSum << "\n";
    //std::cout << "ht sum:  " << htSum << "\n";

    hEtSum->Fill(0.5*etSum);
    //std::cout << "et sum = " << etSum << std::endl;
    for(int bin=0; bin<nEtSumBins; bin++)
      if( (etSum*0.5) >= etSumLo+(bin*etSumBinWidth) ) etSumRates->Fill(etSumLo+(bin*etSumBinWidth)); //GeV
    
    hHtSum->Fill(0.5*htSum);
    //std::cout << "ht sum = " << htSum << std::endl;
    for(int bin=0; bin<nHtSumBins; bin++){
      //std::cout << "Ht? " << upgrade_->sumEt[1]->getType() << std::endl;
      if( (htSum*0.5) >= htSumLo+(bin*htSumBinWidth) ) htSumRates->Fill(htSumLo+(bin*htSumBinWidth)); //GeV
    }

    //hMetSum->Fill(0.5*metSum);
    //std::cout << "met sum = " << metSum << std::endl;
    for(int bin=0; bin<nMetSumBins; bin++)
      if( (metSum*0.5) >= metSumLo+(bin*metSumBinWidth) ) metSumRates->Fill(metSumLo+(bin*metSumBinWidth)); //GeV
        
    //hMhtSum->Fill(0.5*mhtSum]);
    //std::cout << "mht sum = " << mhtSum << std::endl;
    for(int bin=0; bin<nMhtSumBins; bin++){
      //std::cout << "Mht? " << upgrade_->sumEt[1]->getType() << std::endl;
      if( (mhtSum*0.5) >= mhtSumLo+(bin*mhtSumBinWidth) ) mhtSumRates->Fill(mhtSumLo+(bin*mhtSumBinWidth)); //GeV
    }

  }


  //normalisation factor
  double avrgInstLumi = 4.5e33; 
  double sigmaPP = 6.9e-26;
  //double norm = (avrgInstLumi*sigmaPP)/(nevents*1000); //kHz
  double norm = (11.*2244.)/nevents; // zb rate = n_colliding * 11 kHz 
  std::cout << "norm = " << norm << std::endl;

  //make TGraphs
  TGraph* gMuRate = new TGraph(nMuBins);
  TGraph* gEgRate = new TGraph(nEgBins);
  TGraph* gTauRate = new TGraph(nTauBins);
  TGraph* gJetRate = new TGraph(nJetBins);
  TGraph* gEtSumRate = new TGraph(nEtSumBins);
  TGraph* gHtSumRate = new TGraph(nHtSumBins);
  TGraph* gMetSumRate = new TGraph(nMetSumBins);
  TGraph* gMhtSumRate = new TGraph(nMhtSumBins);

  //norm=1;
  for(int bin=0;bin<nMuBins;bin++) gMuRate->SetPoint(bin,muLo+muBinWidth*bin,muRates->GetBinContent(bin+1)*norm);
  for(int bin=0;bin<nEgBins;bin++) gEgRate->SetPoint(bin,egLo+egBinWidth*bin,egRates->GetBinContent(bin+1)*norm);
  for(int bin=0;bin<nTauBins;bin++) gTauRate->SetPoint(bin,tauLo+tauBinWidth*bin,tauRates->GetBinContent(bin+1)*norm);
  for(int bin=0;bin<nJetBins;bin++) gJetRate->SetPoint(bin,jetLo+jetBinWidth*bin,jetRates->GetBinContent(bin+1)*norm);
  for(int bin=0;bin<nEtSumBins;bin++) gEtSumRate->SetPoint(bin,etSumLo+etSumBinWidth*bin,etSumRates->GetBinContent(bin+1)*norm);
  for(int bin=0;bin<nHtSumBins;bin++) gHtSumRate->SetPoint(bin,htSumLo+htSumBinWidth*bin,htSumRates->GetBinContent(bin+1)*norm);
  for(int bin=0;bin<nMetSumBins;bin++) gMetSumRate->SetPoint(bin,metSumLo+metSumBinWidth*bin,metSumRates->GetBinContent(bin+1)*norm);
  for(int bin=0;bin<nMhtSumBins;bin++) gMhtSumRate->SetPoint(bin,mhtSumLo+mhtSumBinWidth*bin,mhtSumRates->GetBinContent(bin+1)*norm);





  ///////////////// RAW DATA ////////////////////
  /// \brief c1
  ///
  TH1D* muRates_RAW = new TH1D("muRates_RAW", "", nMuBins, muLo-2.5, muHi-2.5);
  TH1D* egRates_RAW = new TH1D("egRates_RAW", "", nEgBins, egLo-2.5, egHi-2.5);
  TH1D* tauRates_RAW = new TH1D("tauRates_RAW", "", nTauBins, tauLo-2.5, tauHi-2.5);
  TH1D* jetRates_RAW = new TH1D("jetRates_RAW", "", nJetBins, jetLo-2.5, jetHi-2.5);
  TH1D* etSumRates_RAW = new TH1D("etSumRates_RAW","", nEtSumBins, etSumLo-2.5, etSumHi-2.5);
  TH1D* htSumRates_RAW = new TH1D("htSumRates_RAW","", nHtSumBins, htSumLo-2.5, htSumHi-2.5);
  TH1D* metSumRates_RAW = new TH1D("metSumRates_RAW","", nMetSumBins, metSumLo-2.5, metSumHi-2.5);
  TH1D* mhtSumRates_RAW = new TH1D("mhtSumRates_RAW","", nMhtSumBins, mhtSumLo-2.5, mhtSumHi-2.5);

  TH1D* hMuEt_RAW = new TH1D("muEt_RAW", "", 400, 0.,400.);
  TH1D* hEgEt_RAW = new TH1D("egEt_RAW", "", 400, 0.,400.);
  TH1D* hTauEt_RAW = new TH1D("tauEt_RAW","",400, 0.,400.);
  TH1D* hJetEt_RAW = new TH1D("jetEt_RAW","",400, 0.,400.);
  TH1D* hEtSum_RAW = new TH1D("etSum_RAW","",800, 0.,800.);
  TH1D* hHtSum_RAW = new TH1D("htSum_RAW","",800, 0.,800.);

  // get entries
  Long64_t nentries = treeL1Up_RAW->GetEntriesFast();


  for (Long64_t jentry=0; jentry<nentries; jentry++){
    if((jentry%100000)==0) std::cout << "Done " << jentry  << " events..." << std::endl;
    treeL1Up->GetEntry(jentry);

    //cout << upgrade_->nJets << "\n";

    // get Mu Rates_RAW
    double muEt_RAW(0);
    for(uint it=0; it<upgrade_RAW_->nMuons; ++it){
      // work around a muon bug:
      int offset = upgrade_RAW_->muonQual.size() - upgrade_RAW_->nMuons;
      //cout << "INFO:  " << upgrade_RAW_->nMuons << "\n";
      //cout << "INFO:  " << upgrade_RAW_->muonEt.size() << "\n";
      //cout << "INFO:  " << upgrade_RAW_->muonQual.size() << "\n";
      if (upgrade_RAW_->muonQual[it+offset]==0) continue;
      if(abs(upgrade_RAW_->muonEta[it]) > 0.8) continue;
      hMuEt_RAW->Fill(upgrade_RAW_->muonEt[it]);
      //std::cout<<muEt<<" ----- "<<upgrade_RAW_->muonEt[it]<<endl;
      muEt_RAW = upgrade_RAW_->muonEt[it] > muEt_RAW ?  upgrade_RAW_->muonEt[it]  : muEt_RAW;
    }
    for(int bin=0; bin<nMuBins; bin++)
      if(muEt_RAW >= muLo + (bin*muBinWidth) ) muRates_RAW->Fill(muLo+(bin*muBinWidth)); //GeV

    // get Eg Rates_RAW
    int egEt_RAW(0);
    for(uint it=0; it<upgrade_RAW_->nEGs; ++it){
      hEgEt_RAW->Fill(0.5*upgrade_RAW_->egIEt[it]);
      egEt_RAW = upgrade_RAW_->egIEt[it] > egEt_RAW ?  upgrade_RAW_->egIEt[it]  : egEt_RAW;
    }
    for(int bin=0; bin<nEgBins; bin++)
      if(egEt_RAW*0.5 >= egLo + (bin*egBinWidth) ) egRates_RAW->Fill(egLo+(bin*egBinWidth)); //GeV

    // get Tau Rates_RAW
    int tauEt_RAW(0);
    for(uint it=0; it<upgrade_RAW_->nTaus; ++it){
      hTauEt_RAW->Fill(0.5*upgrade_RAW_->tauIEt[it]);
      tauEt_RAW = upgrade_RAW_->tauIEt[it] > tauEt_RAW ? upgrade_RAW_->tauIEt[it] : tauEt_RAW;
    }
    for(int bin=0; bin<nTauBins; bin++)
      if( (tauEt_RAW*0.5) >= tauLo + (bin*tauBinWidth) ) tauRates_RAW->Fill(tauLo+(bin*tauBinWidth)); //GeV

    // get Jet Rates_RAW
    int jetEt_RAW(0);
    for(uint it=0; it<upgrade_RAW_->nJets; ++it){
      hJetEt_RAW->Fill(0.5*upgrade_RAW_->jetIEt[it]);
      jetEt_RAW =  upgrade_RAW_->jetIEt[it] > jetEt_RAW ? upgrade_RAW_->jetIEt[it] : jetEt_RAW;
    }
    for(int bin=0; bin<nJetBins; bin++)
      if( (jetEt_RAW*0.5) >= jetLo + (bin*jetBinWidth) ) jetRates_RAW->Fill(jetLo+(bin*jetBinWidth));  //GeV

    double etSum_RAW  = -1.0;
    double htSum_RAW  = -1.0;
    double metSum_RAW = -1.0;
    double mhtSum_RAW = -1.0;
    for(uint it=0; it<upgrade_RAW_->nSums; ++it){
      double et = upgrade_RAW_->sumEt[it];
      if (upgrade_RAW_->sumType[it] == L1Analysis::kTotalEt)   etSum_RAW  = et;
      if (upgrade_RAW_->sumType[it] == L1Analysis::kTotalHt)   htSum_RAW  = et;
      if (upgrade_RAW_->sumType[it] == L1Analysis::kMissingEt) metSum_RAW = et;
      if (upgrade_RAW_->sumType[it] == L1Analysis::kMissingHt) mhtSum_RAW = et;
    }
    //std::cout << "mht:  " << mhtSum << "\n";
    //std::cout << "ht sum:  " << htSum << "\n";

    hEtSum_RAW->Fill(0.5*etSum_RAW);
    //std::cout << "et sum = " << etSum << std::endl;
    for(int bin=0; bin<nEtSumBins; bin++)
      if( (etSum_RAW*0.5) >= etSumLo+(bin*etSumBinWidth) ) etSumRates_RAW->Fill(etSumLo+(bin*etSumBinWidth)); //GeV

    hHtSum->Fill(0.5*htSum_RAW);
    //std::cout << "ht sum = " << htSum << std::endl;
    for(int bin=0; bin<nHtSumBins; bin++){
      //std::cout << "Ht? " << upgrade_RAW_->sumEt[1]->getType() << std::endl;
      if( (htSum_RAW*0.5) >= htSumLo+(bin*htSumBinWidth) ) htSumRates_RAW->Fill(htSumLo+(bin*htSumBinWidth)); //GeV
    }

    //hMetSum->Fill(0.5*metSum);
    //std::cout << "met sum = " << metSum << std::endl;
    for(int bin=0; bin<nMetSumBins; bin++)
      if( (metSum_RAW*0.5) >= metSumLo+(bin*metSumBinWidth) ) metSumRates_RAW->Fill(metSumLo+(bin*metSumBinWidth)); //GeV

    //hMhtSum->Fill(0.5*mhtSum]);
    //std::cout << "mht sum = " << mhtSum << std::endl;
    for(int bin=0; bin<nMhtSumBins; bin++){
      //std::cout << "Mht? " << upgrade_RAW_->sumEt[1]->getType() << std::endl;
      if( (mhtSum_RAW*0.5) >= mhtSumLo+(bin*mhtSumBinWidth) ) mhtSumRates_RAW->Fill(mhtSumLo+(bin*mhtSumBinWidth)); //GeV
    }

  }


  //normalisation factor
  double avrgInstLumi = 4.5e33;
  double sigmaPP = 6.9e-26;
  //double norm = (avrgInstLumi*sigmaPP)/(nevents*1000); //kHz
  double norm = (11.*2244.)/nevents; // zb rate = n_colliding * 11 kHz
  std::cout << "norm = " << norm << std::endl;

  //make TGraphs
  TGraph* gMuRate_RAW = new TGraph(nMuBins);
  TGraph* gEgRate_RAW = new TGraph(nEgBins);
  TGraph* gTauRate_RAW = new TGraph(nTauBins);
  TGraph* gJetRate_RAW = new TGraph(nJetBins);
  TGraph* gEtSumRate_RAW = new TGraph(nEtSumBins);
  TGraph* gHtSumRate_RAW = new TGraph(nHtSumBins);
  TGraph* gMetSumRate_RAW = new TGraph(nMetSumBins);
  TGraph* gMhtSumRate_RAW = new TGraph(nMhtSumBins);

  //norm=1;
  for(int bin=0;bin<nMuBins;bin++) gMuRate_RAW->SetPoint(bin,muLo+muBinWidth*bin,muRates_RAW->GetBinContent(bin+1)*norm);
  for(int bin=0;bin<nEgBins;bin++) gEgRate_RAW->SetPoint(bin,egLo+egBinWidth*bin,egRates_RAW->GetBinContent(bin+1)*norm);
  for(int bin=0;bin<nTauBins;bin++) gTauRate_RAW->SetPoint(bin,tauLo+tauBinWidth*bin,tauRates_RAW->GetBinContent(bin+1)*norm);
  for(int bin=0;bin<nJetBins;bin++) gJetRate_RAW->SetPoint(bin,jetLo+jetBinWidth*bin,jetRates_RAW->GetBinContent(bin+1)*norm);
  for(int bin=0;bin<nEtSumBins;bin++) gEtSumRate_RAW->SetPoint(bin,etSumLo+etSumBinWidth*bin,etSumRates_RAW->GetBinContent(bin+1)*norm);
  for(int bin=0;bin<nHtSumBins;bin++) gHtSumRate_RAW->SetPoint(bin,htSumLo+htSumBinWidth*bin,htSumRates_RAW->GetBinContent(bin+1)*norm);
  for(int bin=0;bin<nMetSumBins;bin++) gMetSumRate_RAW->SetPoint(bin,metSumLo+metSumBinWidth*bin,metSumRates_RAW->GetBinContent(bin+1)*norm);
  for(int bin=0;bin<nMhtSumBins;bin++) gMhtSumRate_RAW->SetPoint(bin,mhtSumLo+mhtSumBinWidth*bin,mhtSumRates_RAW->GetBinContent(bin+1)*norm);




  ///////////////////////////////////////////////////////////////////////////




  TMultiGraph *mg_raw = new TMultiGraph();
  gMuRate->SetLineColor(kOrange);
  gMuRate_RAW->SetLineColor(kBlue);
  mg->Add(gMuRate);
  mg->Add(gMuRate_RAW);

  mg->Draw("APL");


















  TCanvas* c1 = new TCanvas;
  c1->SetLogy();

  TLatex n1;
  n1.SetNDC();
  n1.SetTextFont(42);
  n1.SetTextSize(0.04);
   
  TLatex n2;
  n2.SetNDC();
  n2.SetLineWidth(2);
  n2.SetTextFont(61);
  n2.SetTextSize(0.05);
   
  TLatex n3;
  n3.SetNDC();
  n3.SetTextFont(52);
  n3.SetTextSize(0.04);

  TLatex n4;
  n4.SetNDC();
  n4.SetTextFont(52);
  n4.SetTextSize(0.04);

  gEgRate->SetLineWidth(2);
  gEgRate->SetLineColor(kRed);
  gEgRate->GetXaxis()->SetTitle("Threshold [GeV]");
  gEgRate->GetYaxis()->SetTitle("Rate");
  gEgRate->SetMarkerStyle(23);
  gEgRate->SetMarkerColor(kRed);
  gEgRate->GetYaxis()->SetRangeUser(1, 1e7);
  //gEgRate->Draw("APL");
  
  gTauRate->SetLineWidth(2);
  gTauRate->SetLineColor(kBlue);
  gTauRate->GetXaxis()->SetTitle("Threshold [GeV]");
  gTauRate->GetYaxis()->SetTitle("Rate");
  gTauRate->SetMarkerStyle(23);
  gTauRate->SetMarkerColor(kBlue);
  gTauRate->GetYaxis()->SetRangeUser(1, 1e7);
  //gTauRate->Draw("sameAPL");

  gJetRate->SetLineWidth(2);
  gJetRate->SetLineColor(kGreen);
  gJetRate->GetXaxis()->SetTitle("Threshold [GeV]");
  gJetRate->GetYaxis()->SetTitle("Rate");
  gJetRate->SetMarkerStyle(23);
  gJetRate->SetMarkerColor(kGreen);
  gJetRate->GetYaxis()->SetRangeUser(1, 1e7);
  gJetRate->SetTitle("");
  //gJetRate->Draw("sameAPL");

  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gEgRate);
  mg->Add(gTauRate);
  mg->Add(gJetRate);
  mg->SetMinimum(0.1);
  mg->SetMaximum(3E3);
  mg->Draw("APL");
  mg->GetXaxis()->SetTitle("Threshold [GeV]");
  mg->GetYaxis()->SetTitle("Rate [kHz]");
  mg->GetYaxis()->SetRangeUser(1, 1e10);
  gPad->Modified();

  TLegend* leg1 = new TLegend(0.5,0.73,0.7,0.88);
  leg1->SetFillColor(0);
  leg1->AddEntry(gEgRate,"EGamma","lp");
  leg1->AddEntry(gTauRate,"Tau","lp");
  leg1->AddEntry(gJetRate,"Jets","lp");
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->Draw();


  n3.DrawLatex(0.5, 0.6, "Run 260627 #sqrt{s} = 13 TeV");
  n4.DrawLatex(0.5, 0.55, "Zero Bias");
 
  c1->SaveAs("ratesJetEgTau.pdf");


  TCanvas* c2 = new TCanvas;
  c2->SetLogy();

  gMuRate->SetTitle("");
  gMuRate->SetLineWidth(2);
  gMuRate->SetLineColor(kOrange);
  gMuRate->GetXaxis()->SetTitle("Threshold [GeV]");
  gMuRate->GetYaxis()->SetTitle("Rate");
  gMuRate->SetMarkerStyle(23);
  gMuRate->SetMarkerColor(kOrange);
  gMuRate->GetYaxis()->SetRangeUser(1, 1e7);

  gMuRate->Draw("APL");
  gMuRate->GetXaxis()->SetTitle("Threshold [GeV]");
  gMuRate->GetYaxis()->SetTitle("Rate [kHz]");
  gPad->Modified();

    //leg1->Draw();
  TLegend* leg7 = new TLegend(0.5,0.73,0.7,0.88);
  leg7->SetFillColor(0);
  leg7->AddEntry(gMuRate,"Mu","lp");
  leg7->SetBorderSize(0);
  leg7->SetFillStyle(0);
  leg7->Draw();

  n3.DrawLatex(0.5, 0.6, "Run 260627 #sqrt{s} = 13 TeV");
  n4.DrawLatex(0.5, 0.55, "Zero Bias");

  c2->SaveAs("ratesMuon.pdf");

  TCanvas* c3 = new TCanvas;
  c3->SetLogy();

  gEtSumRate->SetLineWidth(2);
  gEtSumRate->SetLineColor(kMagenta);
  gEtSumRate->GetXaxis()->SetTitle("Threshold [GeV]");
  gEtSumRate->GetYaxis()->SetTitle("Rate");
  gEtSumRate->SetMarkerStyle(23);
  gEtSumRate->SetMarkerColor(kMagenta);
  gEtSumRate->GetYaxis()->SetRangeUser(1, 1e7);
  //gEtSumRate->Draw("APL");

  gHtSumRate->SetLineWidth(2);
  gHtSumRate->SetLineColor(kTeal);
  gHtSumRate->GetXaxis()->SetTitle("Threshold [GeV]");
  gHtSumRate->GetYaxis()->SetTitle("Rate");
  gHtSumRate->SetMarkerStyle(23);
  gHtSumRate->SetMarkerColor(kTeal);
  gHtSumRate->GetYaxis()->SetRangeUser(1, 1e7);
  gHtSumRate->SetTitle("");
  //gHtSumRate->Draw("sameAPL");

  TMultiGraph *mgSums = new TMultiGraph();
  mgSums->Add(gEtSumRate);
  mgSums->Add(gHtSumRate);
  mgSums->Draw("APL");
  mgSums->GetXaxis()->SetTitle("Threshold [GeV]");
  mgSums->GetYaxis()->SetTitle("Rate [kHz]");
  gPad->Modified();
  
  TLegend* leg2 = new TLegend(0.7,0.78,0.9,0.88);
  leg2->SetFillColor(0);
  leg2->AddEntry(gEtSumRate,"E_{T}^{total}","lp");
  leg2->AddEntry(gHtSumRate,"H_{T}","lp");
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->Draw("same");

  n3.DrawLatex(0.6, 0.4, "Run 260627 #sqrt{s} = 13 TeV");
  n4.DrawLatex(0.6, 0.25, "Zero Bias");
  
  c3->SaveAs("ratesSums.pdf");





  TCanvas* c4 = new TCanvas;
  c4->SetLogy();

  gMetSumRate->SetLineWidth(2);
  gMetSumRate->SetLineColor(kViolet);
  gMetSumRate->GetXaxis()->SetTitle("Threshold [GeV]");
  gMetSumRate->GetYaxis()->SetTitle("Rate");
  gMetSumRate->SetMarkerStyle(23);
  gMetSumRate->SetMarkerColor(kViolet);
  gMetSumRate->GetYaxis()->SetRangeUser(1, 1e7);
  //gEtSumRate->Draw("APL");

  gMhtSumRate->SetLineWidth(2);
  gMhtSumRate->SetLineColor(kOrange);
  gMhtSumRate->GetXaxis()->SetTitle("Threshold [GeV]");
  gMhtSumRate->GetYaxis()->SetTitle("Rate");
  gMhtSumRate->SetMarkerStyle(23);
  gMhtSumRate->SetMarkerColor(kOrange);
  gMhtSumRate->GetYaxis()->SetRangeUser(1, 1e7);
  gMhtSumRate->SetTitle("");
  //gMhtSumRate->Draw("sameAPL");

  TMultiGraph *mgMsums = new TMultiGraph();
  mgMsums->Add(gMetSumRate);
  mgMsums->Add(gMhtSumRate);
  mgMsums->Draw("APL");
  mgMsums->GetXaxis()->SetTitle("Threshold [GeV]");
  mgMsums->GetYaxis()->SetTitle("Rate [kHz]");
  gPad->Modified();
  
  TLegend* leg3 = new TLegend(0.7,0.78,0.9,0.88);
  leg3->SetFillColor(0);
  leg3->AddEntry(gMetSumRate,"E_{T}^{miss}","lp");
  leg3->AddEntry(gMhtSumRate,"H_{T}^{miss}","lp");
  leg3->SetBorderSize(0);
  leg3->SetFillStyle(0);
  leg3->Draw("same");

  n3.DrawLatex(0.3, 0.4, "Run 260627 #sqrt{s} = 13 TeV");
  n4.DrawLatex(0.3, 0.25, "Zero Bias");
  
 
  
  c4->SaveAs("ratesMETMHT.pdf");






 muRates->Write();
 egRates->Write();
 tauRates->Write();
 jetRates->Write();
 etSumRates->Write();
 htSumRates->Write();
 metSumRates->Write();
 mhtSumRates->Write();

 hMuEt->Write();
 hEgEt->Write();
 hTauEt->Write();
 hJetEt->Write();
 hEtSum->Write();
 hHtSum ->Write();

file_out->Close();
 
}
