#include <vector>
#include <cmath>
#include <iostream>
#include <vector>

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
#include "SoftDropGrooming.h"

#define maxJets 24000

using namespace std;

int main()
{
    //Define histograms for Jet Energy
    vector<TH1F*> hists;
    hists.push_back(new TH1F("Zg distribution", "10 < Jet E < 15", 100, 0, 0.5));
    hists.push_back(new TH1F("Zg distribution", "15 < Jet E < 20", 100, 0, 0.5));
    hists.push_back(new TH1F("Zg distribution", "20 < Jet E < 25", 100, 0, 0.5));
    hists.push_back(new TH1F("Zg distribution", "25 < Jet E < 30", 100, 0, 0.5));
    hists.push_back(new TH1F("Zg distribution", "30 < Jet E < 35", 100, 0, 0.5));
    hists.push_back(new TH1F("Zg distribution", "35 < Jet E < 40", 100, 0, 0.5));
    hists.push_back(new TH1F("Zg distribution", "40 < Jet E", 100, 0, 0.5));
    
    // Open the input file
    TFile *input1 = TFile::Open("./LEP1MC1994_recons_aftercut-002.root", "READ");
    if (!input1 || input1->IsZombie()) {
        cerr << "Error: Could not open input file." << endl;
        return 1;
    }
    
    //Acessing Trees from our root file- the jets and particle trees
    TTree *Tree1 = (TTree *)input1->Get("akR4ESchemeJetTree;1");
    TTree *Tree2 = (TTree *)input1->Get("t");
    if (!Tree1|| !Tree2) {
        cerr << "Error: TTree not found in input file." << endl;
        return 1;
    }
    
    //Arrays for jet-level quantities (from Tree1)
    Float_t jtpt[maxJets], jtphi[maxJets], jteta[maxJets], jtm[maxJets] ;
    
    //Arrays for particle-level quantities (from Tree2)
    Float_t pt[maxJets], phi[maxJets], eta[maxJets], mass[maxJets];
    
    //Defining the single varibles- total particle number and jet number
    Int_t nParticle, nref;
    
    //Acessing branches from Tree1(Jets)
    Tree1->SetBranchAddress("jtpt", jtpt);
    Tree1->SetBranchAddress("jtphi", jtphi);
    Tree1->SetBranchAddress("jteta", jteta);
    Tree1->SetBranchAddress("jtm", jtm);
    Tree1->SetBranchAddress("nref", &nref);
    
    //Acessing branches from Tree2(Particles)
    Tree2->SetBranchAddress("pt", pt);
    Tree2->SetBranchAddress("phi", phi);
    Tree2->SetBranchAddress("eta", eta);
    Tree2->SetBranchAddress("mass", mass);
    Tree2->SetBranchAddress("nParticle",&nParticle);
    
    //Defining our variable for how many times the loop will run
    Long64_t entries = Tree1->GetEntries();
    std::cout << "Number of entries: " << entries << endl;
    
    //Begining of the loop which is over events to build jets from particles and extract substructure observables
    for (Long64_t i = 0; i < 1000; ++i) {
        
        //Get the total number of jets and total number of particle for the next nested loops
        Tree1->GetEntry(i);
        Tree2->GetEntry(i);
        
        //This loop will go through each jet so that we can look inside each jet
        for (int k =0 ; k < nref; ++k) {
            
            if (jtpt[k] <= 0) continue; // Skip invalid entries
            
            //Getting px,py,pz and E for creating each jet
            double px = jtpt[k] * cos(jtphi[k]);
            double py = jtpt[k] * sin(jtphi[k]);
            double pz = jtpt[k] * sinh(jteta[k]);
            double p = jtpt[k] * cosh(jteta[k]);
            double E = sqrt((p * p) + (jtm[k] * jtm[k])); // Simplified for massless particles
            
            //Creating individual Jet indentities. Jet is our object of type "four vector"
            FourVector jetmom(E, px, py, pz);
            
            vector<FourVector> Particles; // defining four vectors ( Four dimensional vector space) creating this to hold a list of particles for each jet
            
            //Starting the particle loop
            for (int j = 0; j < nParticle; ++j) {
                
                //getting px,py,pz and E for creating each particle
                double parpx = pt[j] * cos(phi[j]);
                double parpy = pt[j] * sin(phi[j]);
                double parpz = pt[j] * sinh(eta[j]);
                double parp = pt[j] * cosh(eta[j]);
                double parE = sqrt((parp * parp) + (mass[j] * mass[j]));
                
                //Creating identifiable individual particle indentities. Particle is our object of type "four vector". We will see if they belong to our jet and then put then in a list one by one
                FourVector par(parE, parpx, parpy, parpz);
                
                //Defining the distance between partcle and jet to see if the particle is close enough to belong to a certain jet
                double Distance = GetAngle(jetmom, par);
                
                //Our Jets have .4 rad radius. If particle is not within .4 radians from jet axis than its not a part of jet.
                if (Distance <= 0.4) {
                    Particles.push_back(par); //putting our particle object and making a list of particle for each jet.
                }
            }
            
            //Defining Nodes array. Node * is a type of object thats like an arrow ( pointer) that connects one particle to other.
            vector<Node *> Nodes; // Making a array of Nodes.
            for (size_t iP = 0; iP < Particles.size(); ++iP) { // we are going through each member of our particle array- the set of particles for each jet.
                Node *NewNode = new Node(Particles[iP]); // Create a Node (box) that wraps this particle and can later be linked in the jet tree
                Nodes.push_back(NewNode); //Putting the created Node in array
            }
            
            //Using BuildCATree to connect nodes to each other and make a tree
            BuildCATree(Nodes);
            
            //Soft drop Grooming Funtion from it's header file
            if (!Nodes.empty()) {
                Node *SDNode = SoftDrop(Nodes[0], 0.1, 0.0,0.4);
                
                if (SDNode && SDNode->Child1 && SDNode->Child2) {
                    double PT1 = SDNode->Child1->P[0]; //Momentum of child 1
                    double PT2 = SDNode->Child2->P[0]; //Momentum of child 2
                    double ZG = std::min(PT1,PT2) / (PT1 + PT2); //Fractional momentum
                    double E = Nodes[0]->P[0]; // Total energy of the jet
                    
                    // Fill histograms based on jet energy
                    if (E >= 10 && E < 15) hists[0]->Fill(ZG);
                    else if (E >= 15 && E < 20) hists[1]->Fill(ZG);
                    else if (E >= 20 && E < 25) hists[2]->Fill(ZG);
                    else if (E >= 25 && E < 30) hists[3]->Fill(ZG);
                    else if (E >= 30 && E < 35) hists[4]->Fill(ZG);
                    else if (E >= 35 && E < 40) hists[5]->Fill(ZG);
                    else if (E >= 40) hists[6]->Fill(ZG);
                    
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
    
    // Create canvas and plot ZG histograms
    TCanvas *c1 = new TCanvas("c1", "Jet z_G Distribution", 1600, 1200);
    c1->Divide(3, 3);
    
    //Draw Histograms for ZG
    for (int i = 0; i < 7; ++i) {
        c1->cd(i + 1);
        hists[i]->Draw();
    }
    
    //Saving the canvases
    c1->SaveAs("SoftDrop_JetZG_Distribution.pdf");
    
    // Clean up
    for (TH1F* hist : hists) {
        delete hist;
    }
    
    //Closing the inputs
    input1->Close();
    delete input1;
    return 0;
}





