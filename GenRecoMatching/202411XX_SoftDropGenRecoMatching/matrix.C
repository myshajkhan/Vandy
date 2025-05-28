//  example.c
//  Created by Mysha khan
#include <string>
#include "TH1F.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TCut.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TObject.h"
#include "TBranch.h"
#include "TCanvas.h"
#include <cmath>
 
#define maxJets 24000
void matrix()
{
    
    
    
    TVector3 GenJet;
    TVector3 RecoJet;

    
    

    
    
    
    
    
    // Load your declustering library
 
    TH1F *hist = new TH1F("hist", "Jet Energy Distribution", 100, 0, 50); //define histogrsm
    TH1F *eff = new TH1F("eff", "Effiency", 100, 0, 1); //define histogrsm
    TH1F *fake = new TH1F("fake", "Fake rate", 100, 0, .10); //define histogrsm
    
    TH1F *bestdistance= new TH1F("bestdistance", "Best Distance Distribution ", 100, 0, 0.5); //define
    TH1F *dijet2hist = new TH1F("dijet2hist", "DiJet 2 Energy Distribution", 100, 0, 100); //define
    std::cout << " defined hist"<< std::endl;
    TH1F *gen = new TH1F("gen", "gen", 100, 0, 100); //define
    std::cout << " defined hist"<< std::endl;
    
    TH2* h2 = new TH2F("h2", "h2 title",100, 0, 60, 100, 0, 60);
    
    //Reading input files
    TFile *input1 = TFile::Open("./LEP1MC1994_recons_aftercut-001.root", "READ");
   
    std::cout << "read files "<< std::endl;
    //Making sure the input files work
    if (!input1 || input1-> IsZombie())
    {std::cout<< "ERROR: Could not open data file" << std::endl;
        return; // Properly exit if the file can't be opened
    }
    std::cout << " reading complete"<< std::endl;
    TTree *Tree1 = (TTree *)input1->Get("akR4ESchemeJetTree;1");
    TTree *Tree2 = (TTree *)input1->Get("akR4ESchemeGenJetTree;1");

    std::cout << " all trees accessed"<< std::endl;
    if (!Tree1 || !Tree2) {
        std::cout << "ERROR: Tree not found" << std::endl;
        return; // Exit if tree is not found
    }
    
    
    // Declare variables to hold the data from the tree branches
    
    Float_t e[maxJets], jtpt[maxJets], jtphi[maxJets], jteta[maxJets], jtm[maxJets],
    jtpt2[maxJets], jtphi2[maxJets], jteta2[maxJets], jtm2[maxJets], zg_Beta0p00ZCut0p10[maxJets], zgJtPt_Beta0p00ZCut0p10[maxJets],zgJtEta_Beta0p00ZCut0p10[maxJets], zgJtPhi_Beta0p00ZCut0p10[maxJets],rg_Beta0p00ZCut0p10[maxJets], pt[maxJets], phi[maxJets], eta[maxJets], mass[maxJets];
    Int_t nref, nref2;
    std::cout << "variables defined " << std::endl;
    
    
  
    
    //Acessig necessary branches
    Tree1->SetBranchAddress( "jtpt", &jtpt );
    Tree2->SetBranchAddress( "jtpt", &jtpt2 );
    Tree1->SetBranchAddress( "jtphi", &jtphi );
    Tree2->SetBranchAddress( "jtphi", &jtphi2 );
    Tree1->SetBranchAddress( "jtm", &jtm );
    Tree2->SetBranchAddress( "jtm", &jtm2 );
    Tree1->SetBranchAddress( "jteta", &jteta );
    Tree2->SetBranchAddress( "jteta", &jteta2 );
    Tree2->SetBranchAddress( "nref", &nref2 );
    Tree1->SetBranchAddress( "nref", &nref );
    
    
    
    
    
    //Get the total number of entries in the tree. Long64_t: This is a typedef (type definition) for a signed 64-bit integer
    Long64_t entries = Tree1->GetEntries();
    Float_t pz, p, P, pz1, p1, P1,E, E1, E2 ;
    

    
    
    
   
    std::cout << "loop starting " << std::endl;
   
   // Long64_t entries = Tree1->GetEntries();
    for ( Long64_t i=0; i< entries; ++i ){
        Tree1->GetEntry(i); // get tree entries
        Tree2->GetEntry(i);
        //creating arrays for calculting efficiency
           int total_gen_jets = 0;
           int correctly_matched_gen_jets = 0;
           int total_reco_jets = 0;
           int incorrectly_matched_reco_jets = 0;

        
            total_gen_jets +=nref2;
            total_reco_jets +=nref;
            
        
        for ( Long64_t j=0; j<nref2 ; ++j ){
            
            double BestDistance = -1;
            int BestIndex = -1;
            
            GenJet.SetPtEtaPhi(jtpt2[j], jteta2[j], jtphi2[j]);
        
            
            for ( Long64_t k=0; k< nref; ++k ){
                
                RecoJet.SetPtEtaPhi(jtpt[k], jteta[k], jtphi[k]);
                
                double Distance = GenJet.Angle(RecoJet);
                if (BestIndex < 0 || Distance < BestDistance) {
                    BestIndex = k;
                    BestDistance = Distance;
                    
                 //   std::cout << "hello 2" << std::endl;
                }
              //  std::cout << "hello 3" << std::endl;
            }
            
            
         //   std::cout << "hello 4" << std::endl;
            if (BestIndex >= 0 ) {
                if (BestDistance<0.3) {
                           E1 = std::sqrt((jtpt[j] * std::cosh(jteta[j])) * (jtpt[j] * std::cosh(jteta[j])) + (jtm[j] * jtm[j]));
                           E2 = std::sqrt((jtpt2[BestIndex] * std::cosh(jteta2[BestIndex])) * (jtpt2[BestIndex] * std::cosh(jteta2[BestIndex])) + (jtm2[BestIndex] * jtm2[BestIndex]));

                           // Fill the 2D histogram
                            h2->Fill(E1, E2);
                           // bestdistance->Fill( BestDistance );
                    //fill in the array
                    bestdistance->Fill( BestDistance );
                        correctly_matched_gen_jets ++ ;
                    
                    // Calculate matching efficiency rate
                   
                    
                }
                
                else{
                    
                    
                    incorrectly_matched_reco_jets ++;
                    // Calculate fake matching rate

                  
                }
                
                
                       }
            
           
            
                   }
        
        double matching_efficiency = static_cast<double>(correctly_matched_gen_jets) / total_gen_jets;
   
     
        eff->Fill( matching_efficiency );
        
        
        
        double fake_rate = static_cast<double>(incorrectly_matched_reco_jets) / total_reco_jets;
        fake->Fill(fake_rate);
               }

    // Calculate matching efficiency and fake rate
     //  double matching_efficiency = static_cast<double>(correctly_matched_gen_jets) / total_gen_jets;
      // double fake_rate = static_cast<double>(incorrectly_matched_reco_jets) / total_reco_jets;

    
    // Fill efficiency and fake rate histograms
     
     
    std::cout << "ending nested loop " << std::endl;
    
//Normalize the histogram
/*
  if (bestdistance->Integral() != 0) {
      bestdistance->Scale(1.0 / bestdistance->Integral());
    
}
    std::cout << " normalized hist ";
 */
    
    
    
    
    
    
    
    // Create a new file to save histograms
    TFile *outputFile = new TFile("output_histograms.root", "RECREATE");
    if (!outputFile || outputFile->IsZombie())
    {
        std::cout << "ERROR: Could not create output file" << std::endl;
        return;
    }
    //hist->Write();
   // histgen->Write();
    h2->Write();
    bestdistance->Write();
    
    outputFile->Close();
    std::cout << "Histograms saved to output file" << std::endl;
    
    
    
    
    
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
    c1->Divide(2, 2);
    

   
    
    
    c1->cd(1);
   // dijet1hist->GetXaxis()->SetRangeUser(20,80);
    bestdistance->Draw();
    
    
    
    c1->cd(2);
    c1->cd(2)->SetLogz();
   // dijet1hist->GetXaxis()->SetRangeUser(20,80);
    h2->Draw();
    
    
    c1->cd(3);
  
    eff->Draw();
    
    c1->cd(4);
  
    fake->Draw();
    

    
    c1->SaveAs("example.pdf(","pdf");
    
    
    
}
