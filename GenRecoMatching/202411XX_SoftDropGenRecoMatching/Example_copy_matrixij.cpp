#include <iostream>
#include <vector>
#include <map>
#include <cmath>
using namespace std;
#include "TFile.h"
#include "TTree.h"
#include "Matching.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "JetCorrector.h"
#include "TGraph.h"
#include "TLorentzVector.h"
#include "grooming.h"
#include "TStyle.h"  // Required for gStyle

#define maxJets 2400


// We are creating construction fuction "Jet".
struct Jet
{
    
//We want to make it public so it can be accesed from anywhere- outside the code.
public:
    //Declaring the member varibales for Jet
    double PT;
    double Eta;
    double Phi;
    double E;
    double ZG;
    
    Jet(double pt, double eta, double phi, double e, double zg)
    : PT(pt), Eta(eta), Phi(phi), E(e), ZG(zg)
    {
        
    }
    
};


//Metric fuction is finding the angle between gen and reco jet
double Metric(Jet Gen, Jet Reco)
{
    double gen_px =  Gen.PT * cos(Gen.Phi);
    double gen_py =  Gen.PT * sin(Gen.Phi);
    double gen_pz =  Gen.PT * sinh(Gen.Eta);
    double gen_p = Gen.PT * cosh(Gen.Eta);
    
    
    
    
    
    double reco_px =  Reco.PT * cos(Reco.Phi);
    double reco_py =  Reco.PT * sin(Reco.Phi);
    double reco_pz =  Reco.PT * sinh(Reco.Eta);
    double reco_p = Reco.PT * cosh(Reco.Eta);
    double gen_reco_dot = (gen_px * reco_px) + (gen_py * reco_py) +  (gen_pz * reco_pz) ;
    
    
    
    double cos_angle = gen_reco_dot / (gen_p *  reco_p) ;
    
    
    // avoing cos value larger than 1, so we dont get arc cos = nan
    if ( cos_angle > 0.9999999 ) cos_angle = 0.9999 ;
    if ( cos_angle < -0.9999999 ) cos_angle = -0.9999 ;
    return acos(cos_angle);
}


int main()

{
    // Declaring vectors( Gen, reco) that will hold all the generated objects that carry the structure (Jet , TH1F, TH2F)
    vector<Jet> Gen;
    vector<Jet> Reco;
    vector<TH1F*> Angle;
    vector <TH1F*> ZG;
    vector < TH2F*> correlation;
    vector <TH2F*> EZG;
    vector<FourVector> Particles;
    
    
    // Open the input file
    TFile *input1 = TFile::Open("./LEP1MC1994_recons_aftercut-002.root", "READ");
    
    // Reading the text files for enegy correction
    JetCorrector JEC({"JEC_EEAK4_MC_20210513.txt", "JEC_EEAK4_DataL2_20210514.txt", "JEC_EEAK4_DataL3_20210514.txt"});
    
    //declaration of angle plot
    Angle.push_back(new TH1F("Angle_genreco", "0 < angle < math.pi", 100, 0, 3.5));
    
    
    
    //declaration of energy and zg plot
    correlation.push_back( new TH2F("Energy_correlation", " Energy correlation" , 100, 0 , 60 , 100, 0 , 60));
    correlation.push_back( new TH2F("ZG_correlation", " ZG correlation" , 100, 0 , 1 , 100, 0 , 1));
    
    
    if (!input1 || input1->IsZombie()) {
        cerr << "Error: Could not open input file." << endl;
        return 1;
    }
    
    
    
    //declaration of energy and zg plot

    // Define energy bin edges (0 to 60 in steps of 5, for example)
    std::vector<double> energyBins = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55};

    // Create a 2D histogram for every pair of energy bins
    std::vector<std::vector<TH2F*>> EZGMatrix; // Store 2D histograms in a matrix-like structure
    for (size_t i = 0; i < energyBins.size() - 1; ++i) {
        std::vector<TH2F*> row;
        for (size_t j = 0; j < energyBins.size() - 1; ++j) {
            TString name = Form("EZG_%lu_%lu", i, j);
            TString title = Form("%g < Gen E < %g vs %g < Reco E < %g",
                                  energyBins[i], energyBins[i + 1],
                                  energyBins[j], energyBins[j + 1]);
            row.push_back(new TH2F(name, title, 100, 0.05, .55, 100, 0.05, .55)); // ZG on both axes
        }
        EZGMatrix.push_back(row);
    }

    
    
    
    
    
    //Acessing the trees in the root file and renaming them
    TTree *Tree1 = (TTree *)input1->Get("akR4ESchemeJetTree");
    TTree *Tree2 = (TTree *)input1->Get("akR4ESchemeGenJetTree");
    TTree *Tree3 = (TTree *)input1->Get("t");
    TTree *Tree4 = (TTree *)input1->Get("tgen");
    
    if (!Tree1|| !Tree2) {
        cerr << "Error: TTree not found in input file." << endl;
        return 1;
    }
    
    
    // Set branch addresses
    Float_t jtptReco[maxJets], jtphiReco[maxJets], jtetaReco[maxJets], jtptGen[maxJets], jtphiGen[maxJets], jtetaGen[maxJets],  jtmGen[maxJets],  jtmReco[maxJets], ptGen[maxJets], phiGen[maxJets], etaGen[maxJets], massGen[maxJets], ptReco[maxJets], phiReco[maxJets],etaReco[maxJets], massReco[maxJets];
    
    
    // in qoutes is the actual branch name in the data file
    Int_t nParticle, nrefReco , nrefGen ,nParticleGen, nParticleReco  ;
    Tree1->SetBranchAddress("jtpt", jtptReco);
    Tree1->SetBranchAddress("jtphi", jtphiReco);
    Tree1->SetBranchAddress("jteta", jtetaReco);
    Tree1->SetBranchAddress("nref", &nrefReco);
    Tree1->SetBranchAddress("jtm", jtmReco);
    
    Tree2->SetBranchAddress("jtpt", jtptGen);
    Tree2->SetBranchAddress("jtphi", jtphiGen);
    Tree2->SetBranchAddress("jteta", jtetaGen);
    Tree2->SetBranchAddress("nref", &nrefGen);
    Tree2->SetBranchAddress("jtm", jtmGen);
    
    
    Tree4->SetBranchAddress("pt", ptGen);
    Tree4->SetBranchAddress("phi", phiGen);
    Tree4->SetBranchAddress("eta", etaGen);
    Tree4->SetBranchAddress("mass", massGen);
    Tree4->SetBranchAddress("nParticle",&nParticleGen);
    
    
    Tree3->SetBranchAddress("pt", ptReco);
    Tree3->SetBranchAddress("phi", phiReco);
    Tree3->SetBranchAddress("eta", etaReco);
    Tree3->SetBranchAddress("mass", massReco);
    Tree3->SetBranchAddress("nParticle",&nParticleReco);
    
    
    
    
    
    
    
    // Loop over events and match generated and reconstructed jets
    Long64_t nEntries = min(Tree1->GetEntries(), Tree2->GetEntries());
    for (Long64_t i = 0; i < nEntries; ++i) {
   // for (Long64_t i = 0; i < 10000; ++i) {
        Gen.clear() ;
        Reco.clear();
        Particles.clear();
        Tree1->GetEntry(i);
        Tree2->GetEntry(i);
        Tree3->GetEntry(i);
        Tree4->GetEntry(i);
        // cout << " print i " << i << endl;
        
        
        
        
        // creating gen jet vector with their pT eta phi
        for (int j = 0; j < nrefGen; ++j) {
            
            Particles.clear();
            vector<Node *> Nodesgen;
            //E= p^2  + m^2, we need this for energy correction
            
            double pGen= jtptGen[j] * cosh(jtetaGen[j]);
            double E_gen = sqrt(pGen*pGen + jtmGen[j] * jtmGen[j]) ;
            
            
            
            // loop for particle in the gen jet
            for ( int k= 0; k< nParticleGen ; ++k)
            {
                
                if (jtptGen[j] <= 0) continue; // Skip invalid entries
                
                
                double px = jtptGen[j] * cos(jtphiGen[j]);
                double py = jtptGen[j] * sin(jtphiGen[j]);
                double pz = jtptGen[j] * sinh(jtetaGen[j]);
                double p = jtptGen[j] * cosh(jtetaGen[j]);
                double E = sqrt((p * p) + (jtmGen[j] * jtmGen[j])); // Simplified for massless particles
                
                double parpx = ptGen[k] * cos(phiGen[k]);
                double parpy = ptGen[k] * sin(phiGen[k]);
                double parpz = ptGen[k] * sinh(etaGen[k]);
                double parp = ptGen[k] * cosh(etaGen[k]);
                double parE = sqrt((parp * parp) + (massGen[k] * massGen[k]));
                
                // Particles.push_back(FourVector(E, px, py, pz));
                
                //declaring four vector
                FourVector jetmom(E, px, py, pz);
                FourVector par(parE, parpx, parpy, parpz);
                
                double Distance = GetAngle(jetmom, par);
                
                if (Distance <= 0.4) {
                    Particles.push_back(par);
                    
                }
                
            }
            // storing the result in gen vector
            
            
            
            
            
            // Convert particles to nodes
            
            for (size_t iP = 0; iP < Particles.size(); ++iP) {
                Node *NewNode = new Node(Particles[iP]);
                Nodesgen.push_back(NewNode);
            }
            
            
            // Build the tree from the single-particle nodes
            BuildCATree(Nodesgen);
            //    cout << Nodesgen.size() << "size of nodesgen" << endl;
            
            
            // declare soft drop node
            
            Node *SDNodegen = nullptr ; // starting fro null. Particle in the jet must be greater than zero
            if (Nodesgen.size() >0) {
                // grooming
                SDNodegen =   FindSDNodeE(Nodesgen[0],.1, 0, 0.4);
                
            }
            
            //build the tree with all the collected partilcles in the loop
            // zg calc
            //    cout << SDNodegen << " SDnodegen" << endl;
            
            if (SDNodegen && SDNodegen->Child1 && SDNodegen->Child2) {
                double PT1 = SDNodegen->Child1->P[0];
                double PT2 = SDNodegen->Child2->P[0];
                
               
                    
                    double ZGgen = std::min(PT1,PT2) / (PT1 + PT2);
                    double Egen= Nodesgen[0]->P[0]; // Total energy of the jet
              
                cout<< "Egen gen " << Egen << endl;
                
                    
               
                
                
                Gen.emplace_back(Jet(jtptGen[j], jtetaGen[j],  jtphiGen[j] , E_gen , ZGgen));
                
            }
            
            else
                
                Gen.emplace_back(Jet(jtptGen[j], jtetaGen[j],  jtphiGen[j] , E_gen , -1));
            
        }
        
        
        
        
        
        // Creating reco jet vector with their pT eta phi
        for (int j = 0; j < nrefReco; ++j) {
            Particles.clear();
            //node vector decleration for particle tree creation
            vector<Node *> Nodesreco;
            
            //E= p^2  + m^2, we need this for energy correction
            
            double pReco= jtptReco[j] * cosh(jtetaReco[j]);
            double E_reco = sqrt(pReco*pReco + jtmReco[j] * jtmReco[j]) ;
            
            // eta = -ln tan(theta /2)
            
            double theta = 2*atan( exp (-jtetaReco[j]));
            
            // MJ: There is this energy correction
            // energy correction
            JEC.SetJetE(E_reco);
            JEC.SetJetTheta(theta);
            JEC.SetJetPhi(jtphiReco[j]);
            
            //getting the corrected energy
            double CorrectedEnergy = JEC.GetCorrectedE();
            
            double f_factor = CorrectedEnergy / E_reco ;
            
            
            // loop for particle in the gen jet
            for ( int k= 0; k< nParticleReco ; ++k)
                
            {
                if (jtptReco[j] <= 0) continue; // Skip invalid entries
                
                
                double px = jtptReco[j] * cos(jtphiReco[j]);
                double py = jtptReco[j] * sin(jtphiReco[j]);
                double pz = jtptReco[j] * sinh(jtetaReco[j]);
                double p = jtptReco[j] * cosh(jtetaReco[j]);
                double E = sqrt((p * p) + (jtmReco[j] * jtmReco[j])); // Simplified for massless particles
                
                double parpx = ptReco[k] * cos(phiReco[k]);
                double parpy = ptReco[k] * sin(phiReco[k]);
                double parpz = ptReco[k] * sinh(etaReco[k]);
                double parp = ptReco[k] * cosh(etaReco[k]);
                double parE = sqrt((parp * parp) + (massReco[k] * massReco[k]));
                
                // Particles.push_back(FourVector(E, px, py, pz));
                
                //why am i not declaring like vector<FourVector> jetmom;
                FourVector jetmom(E, px, py, pz);
                FourVector par(parE, parpx, parpy, parpz);
                
                double Distance = GetAngle(jetmom, par);
                
                if (Distance <= 0.4) {
                    Particles.push_back(par);
             
                }
                
            }
            
            
            //Creating nodes for particle
            for (size_t iP = 0; iP < Particles.size(); ++iP) {
                
                Node *NewNode = new Node(Particles[iP]);
                // pushback means adding new things to our Nodesreco array
                Nodesreco.push_back(NewNode);
            }
            
            
            // Build the tree from the single-particle nodes
            BuildCATree(Nodesreco);
            Node *SDNodereco = nullptr ; // starting fro null. Particle in the jet must be greater than zero
            if (Nodesreco.size() >0) {
                // grooming
                SDNodereco =   FindSDNodeE(Nodesreco[0], .1 , 0, 0.4 );
                
            }
            
            
            if (SDNodereco && SDNodereco->Child1 && SDNodereco->Child2) {
                double PT1 = SDNodereco->Child1->P[0];
                double PT2 = SDNodereco->Child2->P[0];
                
                
               
                
              
                    double ZGreco =std::min(PT1,PT2) / (PT1 + PT2);
                    double Ereco= Nodesreco[0]->P[0]; // Total energy of the jet
                // Fill histograms based on jet energy

                
                
                    
                    cout<< "Ereco " << Ereco << endl;
                    
                    
                    Reco.emplace_back(Jet(jtptReco[j]* f_factor, jtetaReco[j],  jtphiReco[j], E_reco, ZGreco ));
                    
                    
                    
                    
                
            }
            //fill hist
            else
                Reco.emplace_back(Jet(jtptReco[j]* f_factor, jtetaReco[j],  jtphiReco[j], E_reco, -1 ));
            
            
        }           // mappping gen and reco jets- collection of pairs
        map<int, int> Match = MatchJetsHungarian(&Metric, Gen, Reco);
        
        
        
        
        for(auto m : Match) {
            // angle between the already matched gen and reco jet, calling back metric function, thats why its in the loop
            if (m.first<0 || m.second<0 )
                continue ;
            double angle =   Metric ( Gen[m.first] , Reco[m.second] );
            
            
            Angle[0]->Fill(angle);
            
            if (angle < 0.4 )
            {   correlation[0]-> Fill(Gen[m.first].E, Reco[m.second].E ) ;
                
                correlation[1]-> Fill(Gen[m.first].ZG, Reco[m.second].ZG) ;
                
                ZG[0]-> Fill(Reco[m.second].ZG);
                ZG[1]-> Fill( )
                
                
                
                // Loop over energy bins
                double genE = Gen[m.first].E;
                 double recoE = Reco[m.second].E;
                 double zgGen = Gen[m.first].ZG;
                 double zgReco = Reco[m.second].ZG;

                 // Loop over all energy bin pairs to find the correct histogram
                 for (size_t i = 0; i < energyBins.size() - 1; ++i) {
                     for (size_t j = 0; j < energyBins.size() - 1; ++j) {
                         if (genE >= energyBins[j] && genE < energyBins[j + 1] &&
                             recoE >= energyBins[i] && recoE < energyBins[i + 1]) {
                             EZGMatrix[i][j]->Fill(zgReco, zgGen); // Fill the appropriate histogram
                             break;
                         }
                        }
                    }
                
            }
            
        }
        
    }
    
    TCanvas *canvas1 = new TCanvas("canvas1", "Gen Reco angle", 1600, 1200);
    
    canvas1->Divide(2, 2);
    
    canvas1->cd(1);
    Angle[0]->Draw();
    canvas1->cd(1)->SetLogy();
    
    
    canvas1->cd(2);
    correlation[1]->Draw();
    canvas1->cd(2)->SetLogz();
    
    
    
    canvas1->cd(3);
    correlation[0]->Draw();
    canvas1->cd(3)->SetLogz();
    
    // Optional TGraph usage
    TGraph G;
    G.SetPoint(0, 0, 0);
    G.SetPoint(1, 100, 100);
    
    
    canvas1->SaveAs("Gen_Reco_angle_Distribution.pdf");
    
    
    // Create a canvas for the matrix plot
    TCanvas *canvas_matrix = new TCanvas("canvas_matrix", "Smearing Matrix", 2000, 2000);
    size_t nRows = energyBins.size() - 1; // Number of energy bins
    size_t nCols = energyBins.size() - 1;
    canvas_matrix->Divide(nCols, nRows); // Create a grid layout (12x12, for example)
   
    canvas_matrix->SetFixedAspectRatio(); // Ensures square proportions

    double padSize = 0.90 / nCols; // Each pad's relative size in the canvas
    std::vector<std::vector<TPad*>> Pads; // Store pads for manual layout

    // Create tightly packed pads
    for (size_t i = 0; i < nRows; ++i) {
        std::vector<TPad*> row;
        for (size_t j = 0; j < nCols; ++j) {
            TString padName = Form("pad_%lu_%lu", i, j);
            double xLow = j * padSize;
            double xHigh = (j + 1) * padSize;
            double yLow = i * padSize; // Flip y-axis
            double yHigh =(i + 1) * padSize;

            TPad *pad = new TPad(padName, "", xLow, yLow, xHigh, yHigh);
            pad->SetMargin(0.01, 0.01, 0.01, 0.01); // Remove space
            pad->Draw();
            row.push_back(pad);
        }
        Pads.push_back(row);
    }
    // Enable color palette
   // gStyle->SetPalette(55);   // Alternative: 55 (kRainBow), 57 (kViridis)
    gStyle->SetOptStat(0);    // Disable stats boxes

    // Draw histograms inside each pad
    for (size_t i = 0; i < nRows; ++i) {
        for (size_t j = 0; j < nCols; ++j) {
            Pads[i][j]->cd();
            EZGMatrix[i][j]->SetTitle(""); // Remove title
            EZGMatrix[i][j]->Draw("COLZ");
            gPad->SetLogz(); // Optional log scale

            // X-axis labels only for bottom row
            if (i == nRows - 1)
                EZGMatrix[i][j]->GetXaxis()->SetTitle("Reco Energy ZG");
            else
                EZGMatrix[i][j]->GetXaxis()->SetLabelSize(0);

            // Y-axis labels only for leftmost column
            if (j == 0)
                EZGMatrix[i][j]->GetYaxis()->SetTitle("Gen Energy ZG");
            else
                EZGMatrix[i][j]->GetYaxis()->SetLabelSize(0);
            
            EZGMatrix[i][j]->SetStats(0); // Disable stats box

        }
    }

    // Save the properly formatted matrix
    canvas_matrix->SaveAs("SmearingMatrix_Fixed.pdf");

  


    // Save the canvas as a PDF
    canvas_matrix->SaveAs("EnergyCorrelationMatrix.pdf");

    
    
    input1->Close();
    delete input1;
    
    return 0;
}


