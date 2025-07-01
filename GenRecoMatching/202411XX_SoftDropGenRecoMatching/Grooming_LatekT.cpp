#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <vector>

#include "TFile.h"
#include "TLorentzVector.h"
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
#include "LateKTGrooming.h"

#define maxJets 24000

using namespace std;

int main()
{
    //Defining histograms for Jet energy
    vector<TH1F*> hists;
    hists.push_back(new TH1F("hist1", "10 < Jet E < 15", 100, 0, 5));
    hists.push_back(new TH1F("hist2", "15 < Jet E < 20", 100, 0, 5));
    hists.push_back(new TH1F("hist3", "20 < Jet E < 25", 100, 0, 5));
    hists.push_back(new TH1F("hist4", "25 < Jet E < 30", 100, 0, 5));
    hists.push_back(new TH1F("hist5", "30 < Jet E < 35", 100, 0, 5));
    hists.push_back(new TH1F("hist6", "35 < Jet E < 40", 100, 0, 5));
    hists.push_back(new TH1F("hist7", "40 < Jet E", 100, 0, 5));
    
    //Define Histotgams for angle
    vector<TH1F*> angle;
    angle.push_back(new TH1F("angle1", "10 < Jet E < 15", 100, 0, .4));
    angle.push_back(new TH1F("angle2", "15 < Jet E < 20", 100, 0, .4));
    angle.push_back(new TH1F("angle3", "20 < Jet E < 25", 100, 0, .4));
    angle.push_back(new TH1F("angle4", "25 < Jet E < 30", 100, 0, .4));
    angle.push_back(new TH1F("angle5", "30 < Jet E < 35", 100, 0, .4));
    angle.push_back(new TH1F("angle6", "35 < Jet E < 40", 100, 0, .4));
    angle.push_back(new TH1F("angle7", "40 < Jet E", 100, 0, .4));
    
    //Define Histotgams for hardness
    vector<TH1F*> hard;
    hard.push_back(new TH1F("hardnessHist_b01", "Hardness Distribution b=0.1", 100, 0, 1.0));
    hard.push_back(new TH1F("hardnessHist_b1", "Hardness Distribution b=1", 100, 0, 1.0));
    hard.push_back(new TH1F("hardnessHist_b2", "Hardness Distribution b=2", 100, 0, 1.0));
    
    // Open the input file
    TFile *input1 = TFile::Open("./LEP1MC1994_recons_aftercut-002.root", "READ");
    if (!input1 || input1->IsZombie()) {
        cerr << "Error: Could not open input file." << endl; //Checking if the input file opens or not
        return 1;
    }
    
    //Accessing Trees inside our input file
    TTree *Tree1 = (TTree *)input1->Get("akR4ESchemeJetTree"); //Jet Tree
    TTree *Tree2 = (TTree *)input1->Get("t"); //Particle Tree
    if (!Tree1|| !Tree2) {
        cerr << "Error: TTree not found in input file." << endl; //Checking if we can access the trees or not
        return 1;
    }
    
    //Arrays for jet-level quantities (from Tree1)
    Float_t jtpt[maxJets], jtphi[maxJets], jteta[maxJets], jtm[maxJets];
    
    //Arrays for Particle-level quantities (from Tree2)
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
            
            vector<FourVector> Particles; //Defining four vectors ( Four dimensional vector space) creating this to hold a list of particles for each jet
            
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
                    Particles.push_back(par);
                }
            }
            
            //Defining Nodes array. Node * is a type of object thats like an arrow ( pointer) that connects one particle to other.
            vector<Node *> Nodes;// Making a array of Nodes.
            for (size_t iP = 0; iP < Particles.size(); ++iP) { //We are going through each member of our particle array- the set of particles for each jet.
                Node *NewNode = new Node(Particles[iP]); //Create a Node (box) that wraps this particle and can later be linked in the jet tree
                Nodes.push_back(NewNode);
            }
            
            //Using BuildCATree to connect nodes to each other and make a tree
            BuildCATree(Nodes);
            
            
            //LateKT Grooming funtion from it's header file
            if (!Nodes.empty()) {
                Node *LastKTNode = LateKTGrooming(Nodes[0], 0.1, 1.0,0.4);
                
                if (LastKTNode && LastKTNode->Child1 && LastKTNode->Child2) {
                    double PT1 = LastKTNode->Child1->P[0]; //Transverse momentum of child 1
                    double PT2 = LastKTNode->Child2->P[0]; //Transverse momentum of child 2
                    double ZG = PT2 / (PT1 + PT2); //Fracional transverse momentum
                    double E = Nodes[0]->P[0]; // Total energy of the jet
                    double Angle = GetAngle(LastKTNode->Child1->P, LastKTNode->Child2->P); // we have the angle between child 1 and 2
                    
                    // we need to define Tlorentz vector for v1
                    TLorentzVector v1 ,v2 ;
                    
                    v1.SetPxPyPzE(LastKTNode->Child1->P[1],LastKTNode->Child1->P[2],LastKTNode->Child1->P[3],LastKTNode->Child1->P[0]); // Child 1's momentum(x,y,z) and energy(PT)
                    v2.SetPxPyPzE(LastKTNode->Child2->P[1],LastKTNode->Child2->P[2],LastKTNode->Child2->P[3],LastKTNode->Child2->P[0]); // Child 2's momentum(x,y,z) and energy(PT)
                    
                    double azimuth1  = v1.Phi();  //Child 1's azimuthal angle
                    double azimuth2  = v2.Phi();  //Child 2's azimuthal angle
                    double azimuth_diff= (azimuth1-azimuth2)*(azimuth1-azimuth2);  //square of their angular difference
                    
                    // now we need to find rapidity of first child 1 and child 2
                    
                    double y1 = (LastKTNode->Child1->P[0]+LastKTNode->Child1->P[3])/(LastKTNode->Child1->P[0]-LastKTNode->Child1->P[3]); //Child 1's rapidity
                    double y2 = (LastKTNode->Child2->P[0]+LastKTNode->Child2->P[3])/(LastKTNode->Child2->P[0]-LastKTNode->Child2->P[3]); //Child 2's rapidity
                    double y_diff = (y1-y2)*(y1-y2); //Square of their rapidity difference
                    
                    //Adding rapidity difference sq and azimuthal angle difference sq
                    double delta_sq = (y_diff ) + (azimuth_diff);
                    
                    // find kt = min(pt_child1, pt_child2)∆ab,
                    double kt = min(LastKTNode->Child1->P[0],LastKTNode->Child2->P[0])* sqrt(delta_sq );
                    
                    // Fill histograms based on jet energy
                    if (E >= 10 && E < 15) hists[0]->Fill(kt);
                    else if (E >= 15 && E < 20) hists[1]->Fill(kt);
                    else if (E >= 20 && E < 25) hists[2]->Fill(kt);
                    else if (E >= 25 && E < 30) hists[3]->Fill(kt);
                    else if (E >= 30 && E < 35) hists[4]->Fill(kt);
                    else if (E >= 35 && E < 40) hists[5]->Fill(kt);
                    else if (E >= 40) hists[6]->Fill(kt);
                    
                    if ( kt> 2) // why is this 2
                    {
                        
                        // Fill histograms based on jet energy
                        if (E >= 10 && E < 15) angle[0]->Fill(Angle);
                        else if (E >= 15 && E < 20) angle[1]->Fill(Angle);
                        else if (E >= 20 && E < 25) angle[2]->Fill(Angle);
                        else if (E >= 25 && E < 30) angle[3]->Fill(Angle);
                        else if (E >= 30 && E < 35) angle[4]->Fill(Angle);
                        else if (E >= 35 && E < 40) angle[5]->Fill(Angle);
                        else if (E >= 40) angle[6]->Fill(Angle);
                    }
                    cout << "PT1 = " << PT1 << ", PT2 = " << PT2 << ", ZG = " << ZG << endl;
                } else {
                    cerr << "Error: LastKTNode or its children are null." << endl;
                    if (!LastKTNode) {
                        cerr << "LastKTNode is null." << endl;
                    } else {
                        if (!LastKTNode->Child1) cerr << "LastKTNode->Child1 is null." << endl;
                        if (!LastKTNode->Child2) cerr << "LastKTNode->Child2 is null." << endl;
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
    
    // Create canvas and plot histograms for Jet kt distribution
    TCanvas *c1 = new TCanvas("c1", "Jet Kt Distribution", 1600, 1200);
    c1->Divide(3, 3);
    for (int i = 0; i < 7; ++i) {
        c1->cd(i + 1);
        hists[i]->Draw();
    }
    
    TCanvas *c2 = new TCanvas("c2", "Jet Angle", 1600, 1200);
    c2->Divide(3, 3);  // Divide based on how many histograms you have
    for (int i = 0; i < 7; ++i) {
        c2->cd(i + 1);
        angle[i]->Draw();
    }
    
    c1->SaveAs("JetZG_Distribution.pdf");
    c2->SaveAs("Jetangle_Distribution.pdf");
    
    // Clean up
    for (TH1F* hist : hists) {
        delete hist;
    }
    
    input1->Close();
    delete input1;
    return 0;
}





