#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>

analysisClass::analysisClass(string * inputList, string * cutFile, string * treeName, string * outputFileName, string * cutEfficFile)
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile)
{
  std::cout << "analysisClass::analysisClass(): begins " << std::endl;
  std::cout << "analysisClass::analysisClass(): ends " << std::endl;
}

analysisClass::~analysisClass()
{
  std::cout << "analysisClass::~analysisClass(): begins " << std::endl;
  std::cout << "analysisClass::~analysisClass(): ends " << std::endl;
}

void analysisClass::Loop()
{
   std::cout << "analysisClass::Loop() begins" <<std::endl;   
    
   if (fChain == 0) return;
   
   //////////book histos here
   
   // No additional cuts
   TH1F *h_HPDHits = new TH1F("h_HPDHits",";HPD Hits",18,0.,18.);
   TH1F *h_HPDNoOtherHits = new TH1F("h_HPDNoOtherHits",";HPDNoOtherHits",18,0.,18.);
   TH1F *h_ADCMaxZeros = new TH1F("h_ADCMaxZeros",";ADCMaxZeros",100,0.,100.);
   TH2F *h_Charge45_R45 = new TH2F("h_Charge45_R45",";Charge45 (fC); R45",100,0,1000,100,-1.5,1.5);
   TH2F *h_NEF_Spike = new TH2F("h_NEF_Spike",";Rejected by SpikeFilter;Rejected by NEF",2,-0.5,1.5,2,-0.5,1.5); 
   TH2F *h_iEta_iPhi_Depth1 = new TH2F("h_iEta_iPhi_Depth1",";iEta;iPhi",61,-30.5,30.5,73,-0.5,72.5);
   TH2F *h_iEta_iPhi_Depth2 = new TH2F("h_iEta_iPhi_Depth2",";iEta;iPhi",61,-30.5,30.5,73,-0.5,72.5);
   TH1F *h_NTuple_HBSumET = new TH1F("h_NTuple_HBSumET",";HBSumET (GeV)",100,0,1000.);
   TH1F *h_NTuple_HESumET = new TH1F("h_NTuple_HESumET",";HESumET (GeV)",100,0,1000.);
   TH1F *h_Calc_HBHESumET = new TH1F("h_Calc_HBHESumET",";HBHESumET (GeV)",100,0,1000.);
   TH1F *h_MuonCSC2DRecHitsSize = new TH1F("h_MuonCSC2DRecHitsSize",";MuonCSC2DRecHitsSize",100,0,150);
   TH1F *h_MuonDT1DRecHitsSize = new TH1F("h_MuonDT1DRecHitsSize",";MuonDT1DRecHitsSize",100,0,150);
   TH1F *h_MuonRPCRecHitsSize = new TH1F("h_MuonRPCRecHitsSize",";MuonRPCRecHitsSize",100,0,150);
   TH1F *h_NominalMET = new TH1F("h_NominalMET",";Nominal MET (GeV)",100,0,1000.);
   TH1F *h_MuonVetoed_NominalMET = new TH1F("h_MuonVetoed_NominalMET",";Nominal MET (GeV)",100,0,40.);
   TH1F *h_MuonThreshold_HBHEmet = new TH1F("h_MuonThreshold_HBHEmet",";HBHE MET (GeV)",100,0,40.);
   TH1F *h_HBHEmet = new TH1F("h_HBHEmet",";HBHE MET (GeV)",100,0,1000.);
   TH1F *h_MuonVetoed_HBHEmet = new TH1F("h_MuonVetoed_HBHEmet",";HBHE MET (GeV)",100,0,40.);

   // Loose cuts
   TH1F *h_Loose_HPDHits = new TH1F("h_Loose_HPDHits",";HPD Hits",19,-0.5,18.5);
   TH1F *h_Loose_HPDNoOtherHits = new TH1F("h_Loose_HPDNoOtherHits",";HPDNoOtherHits",19,-0.5,18.5);
   TH1F *h_Loose_ADCMaxZeros = new TH1F("h_Loose_ADCMaxZeros",";ADCMaxZeros",100,0.,100.);
   TH2F *h_Loose_Charge45_R45 = new TH2F("h_Loose_Charge45_R45",";Charge45 (fC); R45",100,0,1000,100,-1.5,1.5);
   TH2F *h_Loose_NEF_Spike = new TH2F("h_Loose_NEF_Spike",";Rejected by SpikeFilter;Rejected by NEF",2,-0.5,1.5,2,-0.5,1.5); 
   TH2F *h_LoosePlus_NEF_Spike = new TH2F("h_LoosePlus_NEF_Spike",";Rejected by SpikeFilter;Rejected by NEF",2,-0.5,1.5,2,-0.5,1.5); 
   TH2F *h_Loose_iEta_iPhi_Depth1 = new TH2F("h_Loose_iEta_iPhi_Depth1",";iEta;iPhi",61,-30.5,29.5,73,-0.5,71.5);
   TH2F *h_Loose_iEta_iPhi_Depth2 = new TH2F("h_Loose_iEta_iPhi_Depth2",";iEta;iPhi",61,-30.5,29.5,73,-0.5,71.5);
   TH1F *h_Loose_Calc_HBHESumET = new TH1F("h_Loose_Calc_HBHESumET",";HBHESumET (GeV)",100,0,1000.);
   TH1F *h_LoosePlus_Calc_HBHESumET = new TH1F("h_LoosePlus_Calc_HBHESumET",";HBHESumET (GeV)",100,0,1000.);
   TH1F *h_Loose_NominalMET = new TH1F("h_Loose_NominalMET",";Nominal MET (GeV)",100,0,1000.);
   TH1F *h_Loose_MuonVetoed_NominalMET = new TH1F("h_Loose_MuonVetoed_NominalMET",";Nominal MET (GeV)",100,0,40.);
   TH1F *h_Loose_MuonThreshold_HBHEmet = new TH1F("h_Loose_MuonThreshold_HBHEmet",";HBHE MET (GeV)",100,0,40.);
   TH1F *h_Loose_HBHEmet = new TH1F("h_Loose_HBHEmet",";HBHE MET (GeV)",100,0,1000.);
   TH1F *h_Loose_MuonVetoed_HBHEmet = new TH1F("h_Loose_MuonVetoed_HBHEmet",";HBHE MET (GeV)",100,0,40.);
   

   /////////initialize variables

   Long64_t nentries = fChain->GetEntries();
   std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   

   ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
   ////// If the root version is updated and rootNtupleClass regenerated,     /////
   ////// these lines may need to be updated.                                 /////    
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if(jentry < 10 || jentry%1000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
     // if (Cut(ientry) < 0) continue;

     ////////////////////// User's code starts here ///////////////////////

     ///Stuff to be done every event

     // No additional cuts
     h_HPDHits->Fill(HPDHits->at(0));
     h_HPDNoOtherHits->Fill(HPDNoOtherHits->at(0));
     h_ADCMaxZeros->Fill(MaxZeros->at(0));
     
     Double_t calc_HBHESumET = 0.;
     for (unsigned int i = 0; i < HBHERecHitIEta->size(); i++)
     {
     vector<double> FC = HBHERecHitAuxFC->at(i);
     if (FC[4]+FC[5] > 10) h_Charge45_R45->Fill(FC[4]+FC[5], (FC[4]-FC[5])/(FC[4]+FC[5]));

     Int_t FlagWord = HBHERecHitFlags->at(i);
     Bool_t bit13 = FlagWordDecoder(FlagWord,13);
     Bool_t bit27 = FlagWordDecoder(FlagWord,27);
     if (abs(HBHERecHitIEta->at(i)) < 28 && (FC[4]+FC[5]+FC[6])>20)
        {
         if (!bit13 && !bit27) h_NEF_Spike->Fill(0.,0.);
         if (bit13 && !bit27) h_NEF_Spike->Fill(1.,0.);
         if (!bit13 && bit27) h_NEF_Spike->Fill(0.,1.);
         if (bit13 && bit27) h_NEF_Spike->Fill(1.,1.);
        }
     if (HBHERecHitDepth->at(i) == 1) h_iEta_iPhi_Depth1->Fill(HBHERecHitIEta->at(i),HBHERecHitIPhi->at(i));
     if (HBHERecHitDepth->at(i) == 2) h_iEta_iPhi_Depth2->Fill(HBHERecHitIEta->at(i),HBHERecHitIPhi->at(i));
     calc_HBHESumET += HBHERecHitEnergy->at(i)/cosh(HBHERecHitIEta->at(i)); 
     }
     
     h_NTuple_HBSumET->Fill(HBSumET->at(0));
     h_NTuple_HESumET->Fill(HESumET->at(0));
     if (HBSumET->at(0) > 0) h_Calc_HBHESumET->Fill(calc_HBHESumET);
     Double_t met = sqrt(pow(NominalMET->at(0),2)+pow(NominalMET->at(1),2));
     h_NominalMET->Fill(met);
     h_MuonCSC2DRecHitsSize->Fill(MuonCSC2DRecHitsSize->at(0));
     h_MuonDT1DRecHitsSize->Fill(MuonDT1DRecHitsSize->at(0));
     h_MuonRPCRecHitsSize->Fill(MuonRPCRecHitsSize->at(0));

     TVector3 METvector(HBET->at(0)+HEET->at(0), HBET->at(1)+HEET->at(1), 0);
     h_HBHEmet->Fill(METvector.Mag());

     if (MuonCSC2DRecHitsSize->at(0)==0 && MuonDT1DRecHitsSize->at(0)==0 && MuonRPCRecHitsSize->at(0)==0)
     {
         h_MuonVetoed_HBHEmet->Fill(METvector.Mag()); 
         h_MuonVetoed_NominalMET->Fill(met);
     }
     if (MuonCSC2DRecHitsSize->at(0)<=4 && MuonDT1DRecHitsSize->at(0)<=12 && MuonRPCRecHitsSize->at(0)<=4)
         h_MuonThreshold_HBHEmet->Fill(METvector.Mag());
     // With OfficialDecisionRun2L cuts
     if (OfficialDecisionRun2L->at(0) == 1)
     {
     h_Loose_HPDHits->Fill(HPDHits->at(0));
     h_Loose_HPDNoOtherHits->Fill(HPDNoOtherHits->at(0));
     h_Loose_ADCMaxZeros->Fill(MaxZeros->at(0));
     
     calc_HBHESumET = 0.;
     Double_t calcPlus_HBHESumET = 0.;
     for (unsigned int i = 0; i < HBHERecHitIEta->size(); i++)
     {
     vector<double> FC = HBHERecHitAuxFC->at(i);
     if (FC[4]+FC[5] > 10) h_Loose_Charge45_R45->Fill(FC[4]+FC[5], (FC[4]-FC[5])/(FC[4]+FC[5]));

     Int_t FlagWord = HBHERecHitFlags->at(i);
     Bool_t bit11 = FlagWordDecoder(FlagWord,11); //HBHEIsolatedNoise
     Bool_t bit13 = FlagWordDecoder(FlagWord,13); //HBHESpikeNoise
     Bool_t bit27 = FlagWordDecoder(FlagWord,27); //HBHENegativeNoise

     if (abs(HBHERecHitIEta->at(i)) < 28 && (FC[4]+FC[5]+FC[6])>20)
        {
         if (!bit13 && !bit27) h_Loose_NEF_Spike->Fill(0.,0.);
         if (bit13 && !bit27) h_Loose_NEF_Spike->Fill(1.,0.);
         if (!bit13 && bit27) h_Loose_NEF_Spike->Fill(0.,1.);
         if (bit13 && bit27) h_Loose_NEF_Spike->Fill(1.,1.);
        }
     if (abs(HBHERecHitIEta->at(i)) < 28 && (FC[4]+FC[5]+FC[6])>20 && !bit11)
        {
         if (!bit13 && !bit27) h_LoosePlus_NEF_Spike->Fill(0.,0.);
         if (bit13 && !bit27) h_LoosePlus_NEF_Spike->Fill(1.,0.);
         if (!bit13 && bit27) h_LoosePlus_NEF_Spike->Fill(0.,1.);
         if (bit13 && bit27) h_LoosePlus_NEF_Spike->Fill(1.,1.);
        }
     if (HBHERecHitDepth->at(i) == 1) h_Loose_iEta_iPhi_Depth1->Fill(HBHERecHitIEta->at(i),HBHERecHitIPhi->at(i));
     if (HBHERecHitDepth->at(i) == 2) h_Loose_iEta_iPhi_Depth2->Fill(HBHERecHitIEta->at(i),HBHERecHitIPhi->at(i));
     calc_HBHESumET += HBHERecHitEnergy->at(i)/cosh(HBHERecHitIEta->at(i)); 
     if (!(bit11 || bit27)) calcPlus_HBHESumET += HBHERecHitEnergy->at(i)/cosh(HBHERecHitIEta->at(i));
     }
     
     if (HBSumET->at(0) > 0) 
     {
         h_Loose_Calc_HBHESumET->Fill(calc_HBHESumET);
         h_LoosePlus_Calc_HBHESumET->Fill(calcPlus_HBHESumET);
     }
     h_Loose_NominalMET->Fill(met);
     
     h_Loose_HBHEmet->Fill(METvector.Mag());

     if (MuonCSC2DRecHitsSize->at(0)==0 && MuonDT1DRecHitsSize->at(0)==0 && MuonRPCRecHitsSize->at(0)==0) 
     {
         h_Loose_MuonVetoed_HBHEmet->Fill(METvector.Mag()); 
         h_Loose_MuonVetoed_NominalMET->Fill(met);
     }
     if (MuonCSC2DRecHitsSize->at(0)<=4 && MuonDT1DRecHitsSize->at(0)<=12 && MuonRPCRecHitsSize->at(0)<=4)
         h_Loose_MuonThreshold_HBHEmet->Fill(METvector.Mag());
    
     }
     
     ////////////////////// User's code ends here ///////////////////////

   } // End loop over events

   //////////write histos 
   h_HPDHits->Write();
   h_HPDNoOtherHits->Write();
   h_ADCMaxZeros->Write();
   h_Charge45_R45->Write();
   h_NEF_Spike->Write();
   h_iEta_iPhi_Depth1->Write();
   h_iEta_iPhi_Depth2->Write();
   h_NTuple_HBSumET->Write();
   h_NTuple_HESumET->Write();
   h_Calc_HBHESumET->Write();
   h_MuonCSC2DRecHitsSize->Write();
   h_MuonDT1DRecHitsSize->Write();
   h_MuonRPCRecHitsSize->Write();
   h_NominalMET->Write();
   h_MuonVetoed_NominalMET->Write();
   h_HBHEmet->Write();
   h_MuonVetoed_HBHEmet->Write();
   h_MuonThreshold_HBHEmet->Write();

   h_Loose_HPDHits->Write();
   h_Loose_HPDNoOtherHits->Write();
   h_Loose_ADCMaxZeros->Write();
   h_Loose_Charge45_R45->Write();
   h_Loose_NEF_Spike->Write();
   h_LoosePlus_NEF_Spike->Write();
   h_Loose_iEta_iPhi_Depth1->Write();
   h_Loose_iEta_iPhi_Depth2->Write();
   h_Loose_Calc_HBHESumET->Write();
   h_LoosePlus_Calc_HBHESumET->Write();
   h_Loose_NominalMET->Write();
   h_Loose_MuonVetoed_NominalMET->Write();
   h_Loose_HBHEmet->Write();
   h_Loose_MuonVetoed_HBHEmet->Write();
   h_Loose_MuonThreshold_HBHEmet->Write();

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
