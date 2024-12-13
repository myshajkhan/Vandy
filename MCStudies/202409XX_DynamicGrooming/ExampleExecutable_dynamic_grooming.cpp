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
    
    
    vector<TH1F*> hard;
    hard.push_back(new TH1F("hardnessHist_b01", "Hardness Distribution b=0.1", 100, 0, 1.0));
    hard.push_back(new TH1F("hardnessHist_b1", "Hardness Distribution b=1", 100, 0, 1.0));
    hard.push_back(new TH1F("hardnessHist_b2", "Hardness Distribution b=2", 100, 0, 1.0));


    
    // Open the input file
    TFile *input1 = TFile::Open("./LEP1MC1994_recons_aftercut-002.root", "READ");
    /// change it to the
    
    
    if (!input1 || input1->IsZombie()) {
        cerr << "Error: Could not open input file." << endl;
        return 1;
    }


    TTree *Tree1 = (TTree *)input1->Get("akR4ESchemeJetTree;1");
    TTree *Tree2 = (TTree *)input1->Get("t");
    if (!Tree1|| !Tree2) {
        cerr << "Error: TTree not found in input file." << endl;
        return 1;
    }


    // Set branch addresses
    Float_t jtpt[maxJets], jtphi[maxJets], jteta[maxJets], jtm[maxJets], pt[maxJets], phi[maxJets], eta[maxJets], mass[maxJets];
    
    
    Int_t nParticle, nref;
    Tree1->SetBranchAddress("jtpt", jtpt);
    Tree1->SetBranchAddress("jtphi", jtphi);
    Tree1->SetBranchAddress("jteta", jteta);
    Tree1->SetBranchAddress("jtm", jtm);
    Tree1->SetBranchAddress("nref", &nref);
    
    Tree2->SetBranchAddress("pt", pt);
    Tree2->SetBranchAddress("phi", phi);
    Tree2->SetBranchAddress("eta", eta);
    Tree2->SetBranchAddress("mass", mass);
    Tree2->SetBranchAddress("nParticle",&nParticle);
    


    Long64_t entries = Tree1->GetEntries();
    std::cout << "Number of entries: " << entries << endl;


    // Loop over entries and read the data
    for (Long64_t i = 0; i < 1000; ++i) {
        Tree1->GetEntry(i);
        Tree2->GetEntry(i);
        for (int k =0 ; k < nref; ++k) {
            vector<FourVector> Particles;
            for (int j = 0; j < nParticle; ++j) {
                if (jtpt[k] <= 0) continue; // Skip invalid entries
                
                
                double px = jtpt[k] * cos(jtphi[k]);
                double py = jtpt[k] * sin(jtphi[k]);
                double pz = jtpt[k] * sinh(jteta[k]);
                double p = jtpt[k] * cosh(jteta[k]);
                double E = sqrt((p * p) + (jtm[k] * jtm[k])); // Simplified for massless particles
                
                double parpx = pt[j] * cos(phi[j]);
                double parpy = pt[j] * sin(phi[j]);
                double parpz = pt[j] * sinh(eta[j]);
                double parp = pt[j] * cosh(eta[j]);
                double parE = sqrt((parp * parp) + (mass[j] * mass[j]));
               // Particles.push_back(FourVector(E, px, py, pz));
                
                //why am i not declaring like vector<FourVector> jetmom;
                FourVector jetmom(E, px, py, pz);
                FourVector par(parE, parpx, parpy, parpz);
                
                double Distance = GetAngle(jetmom, par);
                
                if (Distance <= 0.4) {
                    Particles.push_back(par);
                    
                }
                
            }
            
            
            // Convert particles to nodes
            vector<Node *> Nodes;
            for (size_t iP = 0; iP < Particles.size(); ++iP) {
                Node *NewNode = new Node(Particles[iP]);
                Nodes.push_back(NewNode);
            }
            
            
            // Build the tree from the single-particle nodes
            BuildCATree(Nodes);
            cout << "We're here 2! " << endl;
         
            
            // Example of soft drop declustering with z_cut 0.1, beta 0.0
            //start of lambda funtion- we use dynamic grooming on this part
            auto FindSDNodeLambda = [](Node *HeadNode, double ZCut, double alpha, double R0) -> Node* {
                if(HeadNode == NULL)
                    return NULL;
                bool Done = false;
                Node *Current = HeadNode;
                Node *bestNode = HeadNode;
                double besthard = 0;
                double jetPT= HeadNode->P[0];
             
                
                while(Done == false)
                {
                    if(Current->N == 1)
                        Done = true;
                    else if(Current->N == 2)
                    {
                        // WTF!
                        std::cerr << "Error!  N = " << Current->N << "!" << std::endl;
                    }
                    else if(Current->Child1 == NULL || Current->Child2 == NULL)
                    {
                        // WTF!
                        std::cerr << "Error!  Child NULL while N = " << Current->N << "!" << std::endl;
                    }
                    else
                    {
                        double P1 = Current->Child1->P[0];
                        double P2 = Current->Child2->P[0];
                      
                        // I think z is
                        
                        double PRatio = std::min(P1, P2) / (P1 + P2);
                    
                        // we have the angle, this is theta(i) for dynamic grooming
                        double Angle = GetAngle(Current->Child1->P, Current->Child2->P);
                        
                        double hard= PRatio*(1 -PRatio) *(P1+P2)*pow(Angle/R0 , alpha );
                      // double Threshold = ZCut * std::pow(Angle / R0, Beta);
                        
                       if(hard > besthard)
                       {
                           besthard = hard;
                           bestNode = Current;
                           
                       }
                      //  else
                       // {
                            if(P1 > P2)
                                Current = Current->Child1;
                            else
                                Current = Current->Child2;
                        //}
                    }
                }
               
                
                return bestNode;
            };
            
            
            //end of lambda funtion
            
            
   
            
            //replace this with your function
            if (!Nodes.empty()) {
                Node *SDNode = FindSDNodeLambda(Nodes[0], 0.1, 1.0,0.4);
                
                
                if (SDNode && SDNode->Child1 && SDNode->Child2) {
                    double PT1 = SDNode->Child1->P[0];
                    double PT2 = SDNode->Child2->P[0];
                    double ZG = PT2 / (PT1 + PT2);
                    double E = Nodes[0]->P[0]; // Total energy of the jet
                    std::cout<< "Energy" << E << endl;
                    std::cout<< "ZG" << ZG << endl;
                    //for kappa calculation
    
                
                    // we have the angle, this is theta(i) for dynamic grooming
                    double Angle = GetAngle(SDNode->Child1->P, SDNode->Child2->P);
                    
                    double hardness= ZG*(1-ZG) *(PT1+PT2)*pow(Angle/0.4 ,1 );
                  // double Threshold = ZCut * std::pow(Angle / R0, Beta);
                    
                    double kappa = (1/E)*hardness;
                    
                    // Fill histograms based on jet energy
                    if (E >= 30 && E < 35) hists[0]->Fill(ZG);
                    else if (E >= 35 && E < 40) hists[1]->Fill(ZG);
                    else if (E >= 40) hists[2]->Fill(ZG);
                    else if (E >= 10 && E < 15) hists[3]->Fill(ZG);
                    else if (E >= 15 && E < 20) hists[4]->Fill(ZG);
                    else if (E >= 20 && E < 25) hists[5]->Fill(ZG);
                    else if (E >= 25 && E < 30) hists[6]->Fill(ZG);
                    
                    
                    //fill K
                    
          

                    // Fill histograms based on the jet energy range or hardness distribution
                    if (E >= 30 && E  < 35) {
                        hard[0]->Fill(kappa);
                    }
                    else if (E  >= 35 && E  < 40) {
                        hard[1]->Fill(kappa);
                    }
                    else if (E >= 40) {
                        hard[2]->Fill(kappa);
                    }

                    
                    
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
        }
        
    }
    // Create canvas and plot histograms
    TCanvas *c1 = new TCanvas("c1", "Jet z_G Distribution", 1600, 1200);
    c1->Divide(3, 3);


    for (int i = 0; i < 7; ++i) {
        c1->cd(i + 1);
        hists[i]->Draw();
    }


    TCanvas *c2 = new TCanvas("c2", "Jet Hardness Distribution", 1600, 1200);
    c2->Divide(2, 2);  // Divide based on how many histograms you have

    for (int i = 0; i < 3; ++i) {
        c2->cd(i + 1)->SetLogx();
        hard[i]->Draw();
    }

    c2->SaveAs("JetHardness_Distribution.pdf");

    
    
    c1->SaveAs("JetZG_Distribution.pdf");


    // Clean up
    for (TH1F* hist : hists) {
        delete hist;
    }


    input1->Close();
    delete input1;


    return 0;
}





