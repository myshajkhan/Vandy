using namespace std;

#include <iostream>
#include <vector>
#include <map>
#include <cmath>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLorentzVector.h"
#include "TStyle.h"

#include "JetCorrector.h"
#include "Matching.h"
#include "SoftDropGrooming.h"
#define maxJets 2400

// We are creating construction fuction "Jet".
struct Jet
{
    
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
    double gen_px =  Gen.PT * cos(Gen.Phi); // Calculate the x-component of the gen jet's momentum using transverse momentum and azimuthal angle
    double gen_py =  Gen.PT * sin(Gen.Phi);  // Calculate the y-component of the gen jet's momentum
    double gen_pz =  Gen.PT * sinh(Gen.Eta); // Calculate the z-component of the gen jet's momentum using pseudorapidity
    double gen_p = Gen.PT * cosh(Gen.Eta); // Calculate the magnitude of the gen jet's momentum vector
    double reco_px =  Reco.PT * cos(Reco.Phi); // Calculate the x-component of the reco jet's momentum
    double reco_py =  Reco.PT * sin(Reco.Phi); // Calculate the y-component of the reco jet's momentum
    double reco_pz =  Reco.PT * sinh(Reco.Eta);  // Calculate the z-component of the reco jet's momentum
    double reco_p = Reco.PT * cosh(Reco.Eta); // Calculate the magnitude of the reco jet's momentum vector
    double gen_reco_dot = (gen_px * reco_px) + (gen_py * reco_py) +  (gen_pz * reco_pz) ;  // Compute the dot product between the gen and reco momentum vectors
    double cos_angle = gen_reco_dot / (gen_p *  reco_p) ;  // Compute the cosine of the angle between the two vectors using the dot product formula
    
    // avoing cos value larger than 1, so we dont get arc cos = nan
    if ( cos_angle > 0.9999999999 ) cos_angle = 0.9999999999 ;
    if ( cos_angle < -0.9999999999 ) cos_angle = -0.9999999999 ;
    return acos(cos_angle);
}

int main()
{
    // Declaring vectors( Gen, reco) that will hold all the generated objects that carry the structure (Jet , TH1F, TH2F)
    vector<Jet> Gen; // Vector to store all generated (truth-level) jets for the current event
    vector<Jet> Reco; // Vector to store all reconstructed (detector-level) jets for the current event
    vector<Jet> RecoData; // Vector to store reconstructed jets from actual experimental data (not simulation)
    vector <Jet> FakeRecoJet; // Vector to store reco jets that did not match any gen jet
    vector <Jet> unmatchedGen; // Vector to store gen jets that didn’t find a matching reco jet — potentially lost or misreconstructed
    
    vector<TH1F*> Angle; // Vector of 1D histograms to record the angular distance between matched gen and reco jets
    vector <TH1F*> ZGreco; // Vector of 1D histograms for storing ZG distributions from reconstructed jets (simulated)
    vector <TH1F*> ZGgen; // Vector of 1D histograms for ZG distributions from generated jets (truth-level)
    vector <TH1F*> ZGDatareco; // Vector of 1D histograms for ZG distributions from reconstructed data jets (real data)
    vector <TH2F*> correlation; // Vector of 2D histograms to record correlations — e.g., energy or ZG between Gen and Reco jets
    vector <TH2F*> EZG; // Vector of 2D histograms used for ZG vs ZG plots, likely across energy bins (indexed as a matrix)
    
    vector<vector<TH2F*>> EZGMatrix; // Store 2D histograms in a matrix-like structure
    vector<FourVector> Particles; // Vector to temporarily store four-momentum vectors (e.g., of jet constituents) for grooming
    vector<double> energyBins = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55}; // Define energy bin edges for histogramming or bin-based analysis (12 bins from 0 to 55 in steps of 5)
    
    int numBinsPerRange = 10; // Assuming 100 bins per ZG axis
    
    //Opening our input files
    TFile *input1 = TFile::Open("./LEP1MC1994_recons_aftercut-002.root", "READ"); //MC generated data input data
    TFile *input2 = TFile::Open("./LEP1Data1994P3_recons_aftercut-MERGED.root", "READ"); // Real data input file
    if (!input1 || input1->IsZombie() || !input2 || input2->IsZombie() ) {
        cerr << "Error: Could not open input file." << endl;
        return 1;
    }
    
    // Reading the text files for enegy correction
    JetCorrector JEC({"JEC_EEAK4_MC_20210513.txt", "JEC_EEAK4_DataL2_20210514.txt", "JEC_EEAK4_DataL3_20210514.txt"});
    
    //declaration of angle plot
    Angle.push_back(new TH1F("Angle_genreco", "0 < angle < math.pi", 100, 0, 3.5));
    
    //declaration of energy and zg plot
    correlation.push_back( new TH2F("Energy_correlation", " Energy correlation" , numBinsPerRange, 0 , 60 , 100, 0 , 60));
    correlation.push_back( new TH2F("ZG_correlation", " ZG correlation" , numBinsPerRange, 0 , 1 , 100, 0 , 1));
    
    //Define one ZG histogram per reconstructed energy bin range (12 total), to analyze ZG vs MC Reco E
    for (size_t j = 0; j < energyBins.size() - 1; ++j) {
        TString name = Form("ZGreco_%lu_%lu", j, j+1);
        TString title = Form("%g < Reco E < %g ", energyBins[j], energyBins[j + 1]);
        ZGreco.push_back(new TH1F(name, title, numBinsPerRange, 0.05, .55)); // ZG on both axes
        ZGreco[j]->Sumw2(); //Turning on unceratainty bar
    }
   
    
    //Define one ZG histogram per reconstructed energy bin range (12 total), to analyze ZG vs Data Reco E
    for (size_t j = 0; j < energyBins.size() - 1; ++j) {
        TString name = Form("ZGDatareco_%lu_%lu", j, j+1);
        TString title = Form("%g < DataReco E < %g ", energyBins[j], energyBins[j + 1]);
        ZGDatareco.push_back(new TH1F(name, title, numBinsPerRange, 0.05, .55)); // ZG on both axes
        ZGDatareco[j]->Sumw2(); //Turning on unceratainty bar
    }
    
    
    //Define one ZG histogram per reconstructed energy bin range (12 total), to analyze ZG vs MC Gen
    for (size_t j = 0; j < energyBins.size() - 1; ++j) {
        TString name = Form("ZGgen_%lu_%lu", j, j+1 );
        TString title = Form("%g < Reco E < %g ",energyBins[j], energyBins[j + 1]);
        ZGgen.push_back(new TH1F(name, title, numBinsPerRange, 0.05, .55)); // ZG on both axes
        ZGgen[j]->Sumw2(); //Turning on unceratainty bar
    }
   
    
    // Create a 2D histogram for every pair of energy bins
    for (size_t i = 0; i < energyBins.size() - 1; ++i) {
        std::vector<TH2F*> row;
        for (size_t j = 0; j < energyBins.size() - 1; ++j) {
            TString name = Form("EZG_%lu_%lu", i, j);
            TString title = Form("%g < Gen E < %g vs %g < Reco E < %g",
                                 energyBins[i], energyBins[i + 1],
                                 energyBins[j], energyBins[j + 1]);
            row.push_back(new TH2F(name, title, numBinsPerRange, 0.05, .55, 10, 0.05, .55)); //ZG on both axes
        }
        EZGMatrix.push_back(row);
    }
   
    //Turning on Uncertainty bar for our matrix plot
    for (size_t i = 0; i < EZGMatrix.size(); ++i) {
        for (size_t j = 0; j < EZGMatrix[i].size(); ++j) {
            if (EZGMatrix[i][j]) EZGMatrix[i][j]->Sumw2();
        }
    }

    int numEnergyRanges = EZGMatrix.size(); // Total energy bins
    
    // Define large smearing matrix with correct binning
    TH1F *SmearingZGrecoMatrix = new TH1F("SmearZG_reco", "Reco ZG Matrix", numBinsPerRange * numEnergyRanges, 0, numBinsPerRange * numEnergyRanges);
    SmearingZGrecoMatrix ->Sumw2(); //Turning on unceratainty bar
    
    TH1F *SmearingZGgenMatrix = new TH1F("SmearZG_gen", "Gen ZG Matrix", numBinsPerRange * numEnergyRanges, 0, numBinsPerRange * numEnergyRanges);
    SmearingZGgenMatrix ->Sumw2(); //Turning on unceratainty bar
    
    TH1F *FakeRecoJetsE = new TH1F("FakeRecoJetsenergy", "Fake Reco Jets Energy", numBinsPerRange, 0, 10);
    FakeRecoJetsE ->Sumw2(); //Turning on unceratainty bar
    
    TH1F *RealRecoJetsE = new TH1F("RealRecoJetsenergy", "Real Reco Jets Energy", numBinsPerRange, 0, 10);
    RealRecoJetsE->Sumw2(); //Turning on unceratainty bar
    
    TH1F *UnmatchedGen = new TH1F("unmatchGenenergy", "Unmatched Gen Jet Energy", numBinsPerRange, 0, 10);
    UnmatchedGen->Sumw2(); //Turning on unceratainty bar
    
    TH2F *SmearingMatrix = new TH2F("SmearingMatrix", "Combined Smearing Matrix",
                                    numBinsPerRange * numEnergyRanges, 0, numBinsPerRange * numEnergyRanges, // X-axis
                                    numBinsPerRange * numEnergyRanges, 0, numBinsPerRange * numEnergyRanges); // Y-axis
    SmearingMatrix->Sumw2(); //Turning on unceratainty bar
    
    //Acessing the trees in the root file and renaming them
    TTree *Tree1 = (TTree *)input1->Get("akR4ESchemeJetTree"); // Reconstructed jet tree from MC simulation (detector-level jets)
    TTree *Tree2 = (TTree *)input1->Get("akR4ESchemeGenJetTree"); //Generated jet tree from MC simulation (truth-level jets)
    TTree *Tree3 = (TTree *)input1->Get("t"); //Particle-level information for reconstructed MC jets
    TTree *Tree4 = (TTree *)input1->Get("tgen"); //Particle-level information for generated MC jets (truth constituents)
    TTree *Tree5 = (TTree *)input2->Get("akR4ESchemeJetTree"); // Reconstructed jet tree from real experimental data
    TTree *Tree6= (TTree *)input2->Get("t"); // Particle-level information for reconstructed data jets (constituents from real data)
    
    if (!Tree1|| !Tree2 || !Tree3 || !Tree4 || !Tree5 || !Tree6   ) {
        cerr << "Error: TTree not found in input file." << endl;
        return 1;
    }
    
    // Arrays for jet-level quantities for jet and particles from input 1
    Float_t jtptReco[maxJets], jtphiReco[maxJets], jtetaReco[maxJets],  jtmReco[maxJets]; //MC Jet level arrays for reco jets
    Float_t jtptGen[maxJets], jtphiGen[maxJets], jtetaGen[maxJets], jtmGen[maxJets]; //MC Jet level arrays for gen jets
    Float_t ptReco[maxJets], phiReco[maxJets],etaReco[maxJets], massReco[maxJets]; //MC particle level arrays for reco jets
    Float_t ptGen[maxJets], phiGen[maxJets], etaGen[maxJets], massGen[maxJets]; //MC particle level arrays for gen jets
    
    // Arrays for jet-level quantities for jet and particles from input 2
    Float_t jtptDataReco[maxJets], jtphiDataReco[maxJets], jtetaDataReco[maxJets],  jtmDataReco[maxJets]; //Data Jet level arrays for reco jets
    Float_t ptDataReco[maxJets], phiDataReco[maxJets],etaDataReco[maxJets], massDataReco[maxJets]; //Data particle level arrays for reco jets
    
    //Defining the single varibles- total particle number and jet number
    Int_t  nrefReco,  nrefGen,  nrefDataReco; //Jet level
    Int_t nParticle, nParticleGen, nParticleReco,  nParticleData, nParticleDataReco ; //particle level
    
    //Setting branch address for reconstructed jet tree from MC simulation (detector-level jets)
    Tree1->SetBranchAddress("jtpt", jtptReco);
    Tree1->SetBranchAddress("jtphi", jtphiReco);
    Tree1->SetBranchAddress("jteta", jtetaReco);
    Tree1->SetBranchAddress("nref", &nrefReco);
    Tree1->SetBranchAddress("jtm", jtmReco);
    
    //Setting branch address generated jet tree from MC simulation (truth-level jets)
    Tree2->SetBranchAddress("jtpt", jtptGen);
    Tree2->SetBranchAddress("jtphi", jtphiGen);
    Tree2->SetBranchAddress("jteta", jtetaGen);
    Tree2->SetBranchAddress("nref", &nrefGen);
    Tree2->SetBranchAddress("jtm", jtmGen);
    
    //Setting branch address for particle-level information for reconstructed MC jets
    Tree3->SetBranchAddress("pt", ptReco);
    Tree3->SetBranchAddress("phi", phiReco);
    Tree3->SetBranchAddress("eta", etaReco);
    Tree3->SetBranchAddress("mass", massReco);
    Tree3->SetBranchAddress("nParticle",&nParticleReco);

    //Setting branch address for particle-level information for generated MC jets (truth constituents)
    Tree4->SetBranchAddress("pt", ptGen);
    Tree4->SetBranchAddress("phi", phiGen);
    Tree4->SetBranchAddress("eta", etaGen);
    Tree4->SetBranchAddress("mass", massGen);
    Tree4->SetBranchAddress("nParticle",&nParticleGen);
    
    //Setting branch address for reconstructed jet tree from real experimental data
    Tree5->SetBranchAddress("jtpt", jtptDataReco);
    Tree5->SetBranchAddress("jtphi", jtphiDataReco);
    Tree5->SetBranchAddress("jteta", jtetaDataReco);
    Tree5->SetBranchAddress("nref", &nrefDataReco);
    Tree5->SetBranchAddress("jtm", jtmDataReco);
    
    //Setting branch address particle-level information for reconstructed data jets (constituents from real data)
    Tree6->SetBranchAddress("pt", ptDataReco);
    Tree6->SetBranchAddress("phi", phiDataReco);
    Tree6->SetBranchAddress("eta", etaDataReco);
    Tree6->SetBranchAddress("mass", massDataReco);
    Tree6->SetBranchAddress("nParticle",&nParticleDataReco);
    
    // Loop over events and match generated and reconstructed jets
    Long64_t nEntries = min(Tree1->GetEntries(), Tree2->GetEntries() );
    for (Long64_t i = 0; i < nEntries; ++i) {
        
        Tree1->GetEntry(i); //Reconstructed jet tree from MC simulation
        Tree2->GetEntry(i); //Generated jet tree from MC simulation (truth-level jets)
        Tree3->GetEntry(i); //Particle-level information for reconstructed MC jets
        Tree4->GetEntry(i); //Particle-level information for generated MC jets (truth constituents)
        Tree5->GetEntry(i); //Reconstructed jet tree from real experimental data
        Tree6->GetEntry(i); //Reconstructed data jets
        Gen.clear() ;
        Reco.clear();
        FakeRecoJet.clear();
        unmatchedGen.clear();
        Particles.clear();
        
        // Creating gen jet vector with their pT eta phi
        for (int j = 0; j < nrefGen; ++j) {
            
            if (jtptGen[j] <= 0) continue; // Skip invalid entries
            
            //Getting px,py,pz and E for creating each jet
            double px = jtptGen[j] * cos(jtphiGen[j]);
            double py = jtptGen[j] * sin(jtphiGen[j]);
            double pz = jtptGen[j] * sinh(jtetaGen[j]);
            double pGen= jtptGen[j] * cosh(jtetaGen[j]);
            double E_gen = sqrt(pGen*pGen + jtmGen[j] * jtmGen[j]) ;
            
            //Creating individual Jet indentities. Jet is our object of type "four vector"
            FourVector jetmom(E_gen, px, py, pz);
            
            Particles.clear(); // Cleaning particle vector to have new sets of particles for each jets
            
            // Loop for particle in the gen jet
            for ( int k= 0; k< nParticleGen ; ++k)
            {
                //Getting px,py,pz and E for creating each particle
                double parpx = ptGen[k] * cos(phiGen[k]);
                double parpy = ptGen[k] * sin(phiGen[k]);
                double parpz = ptGen[k] * sinh(etaGen[k]);
                double parp = ptGen[k] * cosh(etaGen[k]);
                double parE = sqrt((parp * parp) + (massGen[k] * massGen[k]));
                
                //Declaring Particle object (four vector)
                FourVector par(parE, parpx, parpy, parpz);
                
                //Distance betwwen jet and each particle to see if the particle is close enough to belong in the jet
                double Distance = GetAngle(jetmom, par);
                
                
                if (Distance <= 0.4) { // if the particle fall under jet radius(.4) we push them in our particle vector
                    Particles.push_back(par); //Filling Particle
                    
                }
            }
            
            // Defining Nodes( for gen jet) array. Node * is a type of object thats like an arrow ( pointer) that connects one particle to other.
            vector<Node *> Nodesgen;
            
            for (size_t iP = 0; iP < Particles.size(); ++iP) { // We are going through each member of our particle array- the set of particles for each jet.
                Node *NewNode = new Node(Particles[iP]); //Each particle is now considered a node. and can later be linked in the jet tree
                Nodesgen.push_back(NewNode); //Putting the created Node in array
            }
                        
            // Build the tree from the single-particle nodes
            BuildCATree(Nodesgen);
            
            Node *SDNodegen = nullptr ; // Starting from null. Particle in the jet must be greater than zero
            if (Nodesgen.size() >0)
            {
               //Using soft drop grooming. if the jet is not empty, this is the delatration of our fuction (is already created in .h file)
                SDNodegen =  SoftDrop(Nodesgen[0],.1, 0, 0.4);
            }

            if (SDNodegen && SDNodegen->Child1 && SDNodegen->Child2) {
                double PT1 = SDNodegen->Child1->P[0]; //Transverse momentum of child 1
                double PT2 = SDNodegen->Child2->P[0]; //Transverse momentum of child 1
                double ZGgen = std::min(PT1,PT2) / (PT1 + PT2); //Fractional momentum of the children. Bigger the difference smaller the Z
                              
                Gen.emplace_back(Jet(jtptGen[j], jtetaGen[j],  jtphiGen[j] , E_gen , ZGgen));// Creating the Gen jet objects for vector entries
            }
            
            else
                
                Gen.emplace_back(Jet(jtptGen[j], jtetaGen[j],  jtphiGen[j] , E_gen , -1)); // If both child- two foward node is not present
            
        }
        
        // Creating reco jet vector with their pT eta phi
        for (int j = 0; j < nrefReco; ++j) {
            
            Particles.clear();
            //node vector decleration for particle tree creation
            vector<Node *> Nodesreco;
            
            double pReco= jtptReco[j] * cosh(jtetaReco[j]);
            double E_reco = sqrt(pReco*pReco + jtmReco[j] * jtmReco[j]) ;
            double theta = 2*atan( exp (-jtetaReco[j]));
            
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
                
                //Declaring jetmom object that holds four information
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
                
                Reco.emplace_back(Jet(jtptReco[j]* f_factor, jtetaReco[j],  jtphiReco[j], E_reco, ZGreco ));
            }
            
            //fill hist
            else
                Reco.emplace_back(Jet(jtptReco[j]* f_factor, jtetaReco[j],  jtphiReco[j], E_reco, -1 ));
            
        }
        
        // mappping gen and reco jets- collection of pairs
        map<int, int> Match = MatchJetsHungarian(&Metric, Gen, Reco);
        
        vector <bool> Exist(Reco.size(), false);
        
        for(auto m : Match) {
            // angle between the already matched gen and reco jet, calling back metric function, thats why its in the loop
            if (m.first >= 0 && m.second > 0 ) {
                
                Exist[m.second]= true; //this is to check if m.second is matched
                
                double angle =   Metric ( Gen[m.first] , Reco[m.second] );
                 
                Angle[0]->Fill(angle);
                
                if (angle < 0.4 )
                {   correlation[0]-> Fill(Gen[m.first].E, Reco[m.second].E ) ;
                    
                    correlation[1]-> Fill(Gen[m.first].ZG, Reco[m.second].ZG) ;
                    
                    // Loop over energy bins
                    double genE = Gen[m.first].E;
                    double recoE = Reco[m.second].E;
                    double zgGen = Gen[m.first].ZG;
                    double zgReco = Reco[m.second].ZG;
                    
                    RealRecoJetsE->Fill (Reco[m.second].E);
                   
                    // Loop over all energy bin pairs to find the correct histogram
                    for (size_t i = 0; i < energyBins.size() - 1; ++i) {
                        for (size_t j = 0; j < energyBins.size() - 1; ++j) {
                            if (genE >= energyBins[j] && genE < energyBins[j + 1] &&
                                recoE >= energyBins[i] && recoE < energyBins[i + 1]) {
                                EZGMatrix[i][j]->Fill(zgReco, zgGen); // Fill the appropriate histogram
                                
                                cout << zgGen << " zg gen"<< endl;
                                cout << zgReco << " zg reco " << endl;
                                cout << EZGMatrix.size() << "matrix size" << endl;
                            }
                        }
                    }
                    
                    for (size_t j = 0; j < energyBins.size() - 1; ++j) {
                        if (recoE >= energyBins[j] && recoE < energyBins[j + 1]) {
                            ZGreco[j]->Fill( zgReco); // Fill the appropriate histogram
                        
                        }
                        
                    }
                    
                    for (size_t j = 0; j < energyBins.size() - 1; ++j) {
                        if (genE >= energyBins[j] && genE < energyBins[j + 1]) {
                            ZGgen[j]->Fill( zgGen); // Fill the appropriate histogram
                            
                        }
                    }
                }
            }// if statement
            
            //gatering unmatced
            else {
                
                if (m.second <= 0  )
                    
                   //Filling up unmatched gen jets 
                    for ( int i = 0; i < Gen.size(); ++i) {
                        if( Exist[i]== false)
                            UnmatchedGen-> Fill( Gen[m.first].E);
                    }
            }
        }
        for ( int i = 0; i < Reco.size(); ++i) {
            if( Exist[i]== false)
                
                FakeRecoJetsE-> Fill (Reco[i].E);
        }
        
    }
    //tranferring hist data
    for (size_t i = 0; i < numEnergyRanges; ++i) { // Loop over energy bins
        for (size_t j = 0; j < numBinsPerRange; ++j) { // Loop over ZG bins
            for (size_t k = 0; k < numBinsPerRange; ++k) { // Loop over ZG bins in Y
                for (size_t l=0; l< numEnergyRanges; ++l){
                    
                    // Map bin j in range i to global bin k
                    int globalX = j + i * numBinsPerRange;
                    int globalY = k + l * numBinsPerRange;
                    
                    //Get bin content and Uncertainty from the small histogram
                    double binContent = EZGMatrix[i][l]->GetBinContent(j,k);
                    double binUncertainty = EZGMatrix[i][l]->GetBinError(j,k);
                    cout << binContent << " GetBinContent(k)" << endl ;
                    
                    // Fill the big histogram
                    SmearingMatrix->SetBinContent(globalX, globalY, binContent);
                    SmearingMatrix->SetBinError(globalX, globalY, binUncertainty);
                }
            }
        }
    }
    
    for (size_t i = 0; i < numEnergyRanges; ++i) { // Loop over energy bins
        for (size_t j = 0; j < numBinsPerRange; ++j) { // Loop over ZG bins
            
            // Map bin j in range i to global bin k
            int globalX = j + i * numBinsPerRange;
            
            
            // Get bin content from the small histogram
            double binContent = ZGreco[i]->GetBinContent(j);
            double binUncertainty =  ZGreco[i]->GetBinError(j);
            cout << binContent << " GetBinContent(k)" << endl ;
            
            // Fill the big histogram
            SmearingZGrecoMatrix->SetBinContent(globalX, binContent);
            SmearingZGrecoMatrix->SetBinError(globalX, binUncertainty);
        }
    }
        
    for (size_t i = 0; i < numEnergyRanges; ++i) { // Loop over energy bins
        for (size_t j = 0; j < numBinsPerRange; ++j) { // Loop over ZG bins
            
            // Map bin j in range i to global bin k
            int globalX = j + i * numBinsPerRange;
                        
            // Get bin content and uncertainty from the small histogram
            double binContent = ZGgen[i]->GetBinContent(j);
            double binUncertainty = ZGgen[i]->GetBinError(j);
            cout << binContent << " GetBinContent(k)" << endl ;
            
            // Fill the big histogram
            SmearingZGgenMatrix->SetBinContent(globalX, binContent);
            SmearingZGgenMatrix->SetBinError(globalX, binUncertainty);
        }
    }

    input1->Close();
    
    std::string outFileName = "Matched_jets.root";  // Ensure the filename is valid
    std::cout << "Attempting to open ROOT file: " << outFileName << std::endl;
    
    TFile *outHistFile = TFile::Open(outFileName.c_str(), "RECREATE");
    
    // Check if file opened successfully
    if (!outHistFile || outHistFile->IsZombie()) {
        std::cerr << "Error: Could not create output file " << outFileName << std::endl;
        return 1;  // Avoid segfault
    }
    
    std::cout << "File opened successfully!" << std::endl;
    
    // Check before writing to avoid segfault
    if (SmearingZGgenMatrix) {
        SmearingZGgenMatrix->Write();
    } else {
        std::cerr << "Error: SmearingZGgenMatrix is nullptr!" << std::endl;
    }
    
    if (SmearingZGrecoMatrix) {
        SmearingZGrecoMatrix->Write();
    } else {
        std::cerr << "Error: SmearingZGrecoMatrix is nullptr!" << std::endl;
    }
    
    if (SmearingMatrix) {
        SmearingMatrix->Write();
    } else {
        std::cerr << "Error: SmearingMatrix is nullptr!" << std::endl;
    }
    
    outHistFile->Close();
    
    TCanvas *canvas1 = new TCanvas("canvas1", "Gen Reco angle", 1600, 1200);
    
    canvas1->Divide(2, 2);
    
    canvas1->cd(1);
    SmearingZGgenMatrix->Draw();
    canvas1->cd(1)->SetLogz();
    
    canvas1->cd(2);
    SmearingZGrecoMatrix->Draw();
    canvas1->cd(2)->SetLogz();
       
    canvas1->cd(3);
    correlation[0]->Draw();
    canvas1->cd(3)->SetLogz();
    
    canvas1->SaveAs("Gen_Reco_angle_Distribution.pdf");
        
    TCanvas *canvas_matrix = new TCanvas("canvas_matrix", "Smearing Matrix", 2000, 2000);
    
    gStyle->SetPalette(55); // Use a nice color palette
  //  gStyle->SetOptStat(0);
    
    SmearingMatrix->Draw("COLZ");
    gPad->SetLogz(); // Optional log scale
    
    TCanvas *canvas3 = new TCanvas("canvas3", "Gen Reco Energy", 1600, 1200);
    
    canvas3->Divide(2, 2);
    
    canvas3->cd(1);
    FakeRecoJetsE->Draw();
    
    canvas3->cd(1)->SetLogz();
    canvas3->cd(2);
    UnmatchedGen->Draw();
    canvas3->cd(2)->SetLogz();
    canvas3->cd(3);
    RealRecoJetsE->Draw();
    canvas3->cd(3)->SetLogz();
    
    canvas3->SaveAs("Gen_Reco_energy.pdf");
    
    // Save the canvas as a PDF
    canvas_matrix->SaveAs("EnergyCorrelationMatrix.pdf");
        
    input1->Close();
    delete input1;
    return 0;
}


