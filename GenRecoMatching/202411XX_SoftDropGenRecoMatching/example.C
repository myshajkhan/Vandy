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
void example()
{
// Load your declustering library
    gSystem->Load("libMathMore");
TH1F *hist = new TH1F("hist", "Jet Energy Distribution", 100, 0, 50); //define histogrsm
    
TH1F *dijet1hist = new TH1F("dijet1hist", "DiJet 1 Energy Distribution", 100, 0, 100); //define
TH1F *dijet2hist = new TH1F("dijet2hist", "DiJet 2 Energy Distribution", 100, 0, 100); //define
    std::cout << " defined hist"<< std::endl;
    
//Reading input files
TFile *input1 = TFile::Open("./LEP1Data1994P1_recons_aftercut-MERGED.root", "READ");
TFile *input2 = TFile::Open("./PureALLPythia91.root", "READ");
    std::cout << "read files "<< std::endl;
//Making sure the input files work
if (!input1 || input1-> IsZombie()|| !input2 || input2-> IsZombie())
{std::cout<< "ERROR: Could not open data file" << std::endl;
    return; // Properly exit if the file can't be opened
}
    std::cout << " reading complete"<< std::endl;
    TTree *Tree1 = (TTree *)input1->Get("akR4ESchemeJetTree;1");
    TTree *Tree2 = (TTree *)input2->Get("tgen;2");
    TTree *Tree3 = (TTree *)input1->Get("akR8ESchemeJetTree");
    TTree *Tree4 = (TTree *)input1->Get("ktN3ESchemeJetTree;1");
    
    
    std::cout << " all trees accessed"<< std::endl;
    if (!Tree4 || !Tree3 ) {
           std::cout << "ERROR: Tree not found" << std::endl;
           return; // Exit if tree is not found
       }
    
    
    // Declare variables to hold the data from the tree branches
    
    Float_t e[maxJets], jtpt[maxJets], jtphi[maxJets], jteta[maxJets], jtm[maxJets], zg_Beta0p00ZCut0p10[maxJets], zgJtPt_Beta0p00ZCut0p10[maxJets],zgJtEta_Beta0p00ZCut0p10[maxJets], zgJtPhi_Beta0p00ZCut0p10[maxJets],rg_Beta0p00ZCut0p10[maxJets];
    std::cout << "variables defined " << std::endl;
    
    
    
   // pz = jtpt.sinh(jteta[maxJets])
    //p= jtpt.cosh(jteta[maxJets])
  //  P= pz+p
    //E= sqrt(P^2+jtm[maxJets]^2)
    
    
    
   // Tree->SetBranchAddress( "jtpt", &jtpt );
    //Acessig necessary branches
    Tree1->SetBranchAddress( "jtpt", &jtpt );
    Tree4->SetBranchAddress( "jtpt", &jtpt );
    Tree2->SetBranchAddress( "e", &e);
    Tree1->SetBranchAddress( "jtphi", &jtphi );
    Tree1->SetBranchAddress( "jteta", &jteta );
    Tree1->SetBranchAddress( "zg_Beta0p00ZCut0p10", &zg_Beta0p00ZCut0p10 );
    Tree1->SetBranchAddress( "zgJtEta_Beta0p00ZCut0p10", &zgJtEta_Beta0p00ZCut0p10 );
    Tree1->SetBranchAddress( "zgJtPhi_Beta0p00ZCut0p10", &zgJtPhi_Beta0p00ZCut0p10 );
    Tree1->SetBranchAddress( "zgJtPt_Beta0p00ZCut0p10", &zgJtPt_Beta0p00ZCut0p10 );
    Tree1->SetBranchAddress( "rg_Beta0p00ZCut0p10", &rg_Beta0p00ZCut0p10 );
    
    
    
//Get the total number of entries in the tree. Long64_t: This is a typedef (type definition) for a signed 64-bit integer
   Long64_t entries = Tree1->GetEntries();
    Float_t pz, p, P, pz1, p1, P1,E, E1, E2 ;
 /*
    std::cout << entries << std::endl ;
    for ( Long64_t i=0; i< 100000; ++i ){
       Tree1->GetEntry(i); // get tree entries
        //now get branch entries
      //  std::cout << " recieved Tree2 entries" ;
        for ( Long64_t j=0; j< 10000; ++j ){
        if (jtpt[j]>10 && jtpt[j]<100) {
        pz = jtpt[j]*std::sinh(jteta[j]);
        p= jtpt[j]*std::cosh(jteta[j]);
       // P= pz+p;
        E= std::sqrt((p*p)+(jtm[j]*jtm[j]));
     //  for ( Long64_t j=0; j< 10000; ++j ){
          
           //std::cout << " recieved branch entries";
        //    if (jtpt[j]>0 && jtpt[j]<40) {
                hist->Fill(E);
        }
         //  std::cout << " filled hist";
        }
    }
    
*/
  //loop 2 for dijet  E 1
    Long64_t entries3 = Tree3->GetEntries();
     
    for ( Long64_t i=0; i< 1000; ++i ){
       Tree3->GetEntry(i); // get tree entries
        //now get branch entries
      //  std::cout << " recieved Tree2 entries" ;
        for ( Long64_t j=0; j< 100; ++j ){
        if (jtpt[j]>4 && jtpt[j]<100) {
        pz1 = jtpt[j]*std::sinh(jteta[j]);
        p1= jtpt[j]*std::cosh(jteta[j]);
       // P= pz+p;
        E1= std::sqrt((p*p)+(jtm[j]*jtm[j]));
     
            dijet1hist->Fill(E1);
        }
        
        }
    }
    
        
        
        //loop 3 for dijet  E 2
          Long64_t entries4 = Tree3->GetEntries();
           
          for ( Long64_t i=0; i< 10000; ++i ){
             Tree3->GetEntry(i); // get tree entries
              //now get branch entries
           
              for ( Long64_t j=0; j< 10000; ++j ){
              if (jtpt[j]>4 && jtpt[j]<60) {
                  std::cout <<"jet pt entries " << jtpt[j] << std::endl;
              pz = jtpt[j]*std::sinh(jteta[j]);
                  
              p= jtpt[j]*std::cosh(jteta[j]);
                  std::cout <<"jet p entries " << p << std::endl;
             // P= pz+p;
              E= std::sqrt((p*p)+(jtm[j]*jtm[j]));
           
                  dijet2hist->Fill(E);
              }
              
              }
          }
          
    
     
    
 
    std::cout << "ending nested loop " << std::endl;
    
     //Normalize the histogram
    if (hist->Integral() != 0) {
      hist->Scale(1.0 / hist->Integral());
    }
    
    if (dijet1hist->Integral() != 0) {
        dijet1hist->Scale(1.0 / dijet1hist->Integral());
  }
        
        
    if (dijet2hist->Integral() != 0) {
            dijet2hist->Scale(1.0 / dijet2hist->Integral());
      }
    std::cout << " normalized hist ";
    
    
    
    
    
    
    
    
    // Create a new file to save histograms
       TFile *outputFile = new TFile("output_histograms.root", "RECREATE");
       if (!outputFile || outputFile->IsZombie())
       {
           std::cout << "ERROR: Could not create output file" << std::endl;
           return;
       }
       hist->Write();
       dijet1hist->Write();
       dijet2hist->Write();
       outputFile->Close();
       std::cout << "Histograms saved to output file" << std::endl;
    
    
    
    
    
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
    c1->Divide(2, 2);
    
    c1->cd(1);
    std::cout << "drawing hist ";
    hist->GetXaxis()->SetRangeUser(20,50);
    hist->Draw();
    std::cout << " hist drawn ";
    
    
    
    
    
    c1->cd(2);
    dijet1hist->GetXaxis()->SetRangeUser(0,80);
    dijet1hist->Draw();
        
        
    c1->cd(3);
        dijet2hist->GetXaxis()->SetRangeUser(20,60);
        dijet2hist->Draw();
        
    
   c1->SaveAs("example.pdf(","pdf");
    
    


