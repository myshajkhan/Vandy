#include <vector>
using namespace std;


#include "CATree.h"
#include "TFile.h"
#include <string>
#include "TH1F.h"
#include "TChain.h"
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
#include <iostream>
#include <vector>


#define maxJets 24000


int main()
{
    // Define histograms
    vector<TH1F*> hists;
    hists.push_back(new TH1F("hist1", "30 < Jet E < 35", 100, 0, 0.5));
    hists.push_back(new TH1F("hist2", "35 < Jet E < 40", 100, 0, 0.5));
    hists.push_back(new TH1F("hist3", "40 < Jet E", 100, 0, 0.5));
    hists.push_back(new TH1F("hist4", "10 < Jet E < 15", 100, 0, 0.5));
    hists.push_back(new TH1F("hist5", "15 < Jet E < 20", 100, 0, 0.5));
    hists.push_back(new TH1F("hist6", "20 < Jet E < 25", 100, 0, 0.5));
    hists.push_back(new TH1F("hist7", "25 < Jet E < 30", 100, 0, 0.5));
    
    // Soft drop parameters
    double zcut = 0.1;
    double beta = 0.0;
    double R0 = 1.0;

    // Open the input file
    TFile *input1 = TFile::Open("./LEP1MC1994_recons_aftercut-001.root", "READ");
    
    
    if (!input1 || input1->IsZombie()) {
        cerr << "Error: Could not open input file." << endl;
        return 1;
    }

    //accessing trees
    TTree *Tree1 = (TTree *)input1->Get("akR4ESchemeJetTree;1");
    TTree *Tree2 = (TTree *)input1->Get("akR4ESchemeGenJetTree;1");
    TTree *Tree3 = (TTree *)input1->Get("t");
    
    
    if (!Tree1 ||!Tree2 ) {
        cerr << "Error: TTree not found in input file." << endl;
        return 1;
    }


    //declare variables
    Float_t jtpt[maxJets], jtphi[maxJets], jteta[maxJets], jtm[maxJets], jtpt2[maxJets], jtphi2[maxJets], jteta2[maxJets], jtm2[maxJets], Parpx[maxJets], Parpy[maxJets],Parpz[maxJets], Parmass[maxJets] ;
    
    
    //
    Int_t nref, nref2, numP;
    
    
    
    
    
    // Set branch addresses
    Tree1->SetBranchAddress("jtpt", jtpt);
    Tree1->SetBranchAddress("jtphi", jtphi);
    Tree1->SetBranchAddress("jteta", jteta);
    Tree1->SetBranchAddress("jtm", jtm);
    
    Tree2->SetBranchAddress("jtpt", jtpt2);
    Tree2->SetBranchAddress("jtphi", jtphi2);
    Tree2->SetBranchAddress("jteta", jteta2);
    Tree2->SetBranchAddress("jtm", jtm2);
    
    Tree2->SetBranchAddress( "nref", &nref2 );
    Tree1->SetBranchAddress( "nref", &nref );
    
    Tree3->SetBranchAddress( "nParticle", &numP );
    Tree3->SetBranchAddress( "px", &Parpx );
    Tree3->SetBranchAddress( "py", &Parpy );
    Tree3->SetBranchAddress( "pz", &Parpz );
    Tree3->SetBranchAddress( "mass", &Parmass );
    
    
    


    Long64_t entries = Tree1->GetEntries();
    std::cout << "Number of entries: " << entries << endl;


    // Loop over entries and read the data
    for (Long64_t i = 0; i < entries; ++i) {
        Tree1->GetEntry(i);
        Tree2->GetEntry(i);
        Tree3->GetEntry(i);

        
       
        
        
  
        
        
        
        
        
        //start of the loop
        for (int j = 0; j < nref; ++j) {
            if (jtpt[j] <= 0) continue; // Skip invalid entries
            
            
            
            
            double px1 = jtpt[j] * cos(jtphi[j]);
            double py1 = jtpt[j] * sin(jtphi[j]);
            double pz1 = jtpt[j] * sinh(jteta[j]);
            double p1 = jtpt[j] * cosh(jteta[j]);
            double E1 = sqrt((p1 * p1) + (jtm[j] * jtm[j]));
            
            //declare a jet
            
            //why am i not declaring like vector<FourVector> jetmom;
            FourVector jetmom(E1, px1, py1, pz1);
            
            //declare a particle
            vector<FourVector> Particles1;
           
            
            // Simplified for massless particles
            
            //particle loop of each jet
            for (int k=0 ; k< numP; ++k  ){
                
                double Emom = sqrt(( Parpx[k]* Parpx[k]) +( Parpy[k]* Parpy[k])+  ( Parpz[k]* Parpz[k])+ (Parmass[k] * Parmass[k]));
                
                FourVector partmom(Emom, Parpx[k], Parpy[k], Parpz[k]);
                
                double Distance = GetAngle(jetmom, partmom);
                
                if (Distance <= 0.4) {
                    Particles1.push_back(partmom);
                    
                }
                
            }
            
            // Convert particles to nodes
            vector<Node *> Nodes;
            for (size_t iP = 0; iP < Particles1.size(); ++iP) {
                Node *NewNode = new Node(Particles1[iP]);
                Nodes.push_back(NewNode);
                
                
                
            }
            
            BuildCATree(Nodes);
            
            
            // Example of soft drop declustering with z_cut 0.1, beta 0.0
            if (!Nodes.empty()) {
                Node *SDNode = FindSDNode(Nodes[0], 0.1, 0.0);  

                

                if (SDNode && SDNode->Child1 && SDNode->Child2) {
                    
                    
                    //get PT and angle
                    
                    
                    //PT
                    double E1 = SDNode->Child1->P[0];
                    double E2 = SDNode->Child2->P[0];;
                    double ZG = E2 / (E1 + E2);
                  //  double E = Nodes[0]->P[0]; // Total energy of the jet

                    //Angle between the nodes
                    double DeltaR12 = GetAngle(SDNode->Child1->P, SDNode->Child2->P);
                    
                    
                    //Now apply soft drop conditions
                   if (ZG> zcut*pow(DeltaR12/R0 , beta )){
                    double E = Nodes[0]->P[0]; // Total energy of the jet
                        
                   
                    

                    // Fill histograms based on jet energy
                    if (E >= 30 && E < 35) hists[0]->Fill(ZG);
                    else if (E >= 35 && E < 40) hists[1]->Fill(ZG);
                    else if (E >= 40) hists[2]->Fill(ZG);
                    else if (E >= 10 && E < 15) hists[3]->Fill(ZG);
                    else if (E >= 15 && E < 20) hists[4]->Fill(ZG);
                    else if (E >= 20 && E < 25) hists[5]->Fill(ZG);
                    else if (E >= 25 && E < 30) hists[6]->Fill(ZG);
                  }

                   // cout << "Soft drop parameters: PT1 = " << PT1 << ", PT2 = " << PT2 << ", ZG = " << ZG << endl;
                } else {
                   // cerr << "Error: SDNode or its children are null." << endl;
                    if (!SDNode) {
                       // cerr << "SDNode is null." << endl;
                    } else {
                       // if (!SDNode->Child1) cerr << "SDNode->Child1 is null." << endl;
                       // if (!SDNode->Child2) cerr << "SDNode->Child2 is null." << endl;
                    }
                }
            } else {
               // cerr << "Error: Nodes vector is empty." << endl;
            }

            // Cleanup
            for (Node *node : Nodes) {
                delete node;
                
                
                
            }
             
        }


        /*
       
        
       

        // Build the tree from the single-particle nodes
                BuildCATree(Nodes1);
                BuildCATree(Nodes2);
        
        
        cout << "We're here 2! " << endl;


        // Example of soft drop declustering with z_cut 0.1, beta 0.0
        if (!Nodes1.empty()) {
            Node *SDNode = FindSDNode(Nodes1[0], 0.1, 0.0);

            

            if (SDNode && SDNode->Child1 && SDNode->Child2) {
                double PT1 = SDNode->Child1->P.GetPT();
                double PT2 = SDNode->Child2->P.GetPT();
                double ZG = PT2 / (PT1 + PT2);
                double E = Nodes[0]->P[0]; // Total energy of the jet


                // Fill histograms based on jet energy
                if (E1 >= 30 && E < 35) hists[0]->Fill(ZG);
                else if (E1 >= 35 && E < 40) hists[1]->Fill(ZG);
                else if (E1 >= 40) hists[2]->Fill(ZG);
                else if (E1 >= 10 && E < 15) hists[3]->Fill(ZG);
                else if (E1 >= 15 && E < 20) hists[4]->Fill(ZG);
                else if (E1 >= 20 && E < 25) hists[5]->Fill(ZG);
                else if (E1 >= 25 && E < 30) hists[6]->Fill(ZG);


                cout << "Soft drop parameters: PT1 = " << PT1 << ", PT2 = " << PT2 << ", ZG = " << ZG << endl;
            } else {
                cerr << "Error: SDNode or its children are null." << endl;
                if (!SDNode) {
                    cerr << "SDNode is null." << endl;
                } else {
                    if (!SDNode->Child1) cerr << "SDNode->Child1 is null." << endl;
                    if (!SDNode->Child2) cerr << "SDNode->Child2 is null." << endl;
                }
            }
        } else {
            cerr << "Error: Nodes vector is empty." << endl;
        }


        // Cleanup
        for (Node *node : Nodes) {
            delete node;
        }
         
        */
    }


    // Create canvas and plot histograms
    TCanvas *c1 = new TCanvas("c1", "Jet z_G Distribution", 1600, 1200);
    c1->Divide(3, 3);


    for (int i = 0; i < 7; ++i) {
        c1->cd(i + 1);
        hists[i]->Draw();
    }


    c1->SaveAs("JetZG_Distribution.pdf");


    // Clean up
    for (TH1F* hist : hists) {
        delete hist;
    }


    input1->Close();
    delete input1;


    return 0;
}





