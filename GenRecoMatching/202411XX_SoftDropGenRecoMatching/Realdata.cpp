#include <iostream>
#include <vector>
#include <map>
#include <cmath>
using namespace std;
#include "ProgressBar.h"
#include "TROOT.h"
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
    vector<Jet> RecoData;
    
    vector<TH1F*> Angle;
    vector < TH2F*> correlation;
    vector <TH2F*> EZG;
    vector <TH1F*> ZGreco;
    vector<FourVector> Particles;
    vector <TH1F*> ZGgen;
    vector <TH1F*> ZGDatareco;
    
    
    
    
    
    int numBinsPerRange = 100; // Assuming 100 bins per ZG axis
    // Open the input file
    
    std::vector<double> energyBins = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55};
    
    // Reading the text files for enegy correction
    JetCorrector JEC({"JEC_EEAK4_MC_20210513.txt", "JEC_EEAK4_DataL2_20210514.txt", "JEC_EEAK4_DataL3_20210514.txt"});
    
    
    //declaration of energy and zg plot
    correlation.push_back( new TH2F("Energy_correlation", " Energy correlation" , 100, 0 , 60 , 100, 0 , 60));
    correlation.push_back( new TH2F("ZG_correlation", " ZG correlation" , 100, 0 , 1 , 100, 0 , 1));
    
    
    for (size_t j = 0; j < energyBins.size() - 1; ++j) {
        TString name = Form("ZGreco_%lu_%lu", j);
        TString title = Form("%g < Reco E < %g ");
        ZGreco.push_back(new TH1F(name, title, 100, 0.05, .55)); // ZG on both axes
    }
    
    
    for (size_t j = 0; j < energyBins.size() - 1; ++j) {
        TString name = Form("ZGDatareco_%lu_%lu", j);
        TString title = Form("%g < DataReco E < %g ");
        ZGDatareco.push_back(new TH1F(name, title, 100, 0.05, .55)); // ZG on both axes
    }
    
    
    
    
    //declaration of energy and zg plot
    
    
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
    
    
    //for mapping Feb 4 3:34
    
    int numEnergyRanges = EZGMatrix.size(); // Total energy bins
    
    
    TH1F *SmearingZGrecoMatrix = new TH1F("SmearZG_reco", "RecoZG", numBinsPerRange * numEnergyRanges, 0, numBinsPerRange * numEnergyRanges);
    
    TH1F *SmearingZGDatarecoMatrix = new TH1F("SmearZG_Datareco", "DataRecoZG", numBinsPerRange * numEnergyRanges, 0, numBinsPerRange * numEnergyRanges);
    
    
    TH1F *RealRecoJetsE = new TH1F("RealRecoJetsenergy", "Real Reco Jets Energy", 100, 0, 10);
    
    
    
    TFile *input1 = TFile::Open("./LEP1MC1994_recons_aftercut-002.root", "READ");
    TFile *input2 = TFile::Open("./LEP1Data1994P3_recons_aftercut-MERGED.root", "READ");
    if (!input1 || input1->IsZombie()) {
        cerr << "Error: Could not open input file." << endl;
        return 1;
    }
    
    
    
    
    
    
    TTree *Tree5 = (TTree *)input2->Get("akR4ESchemeJetTree");
    
    TTree *Tree6= (TTree *)input2->Get("t");
    
    
    
    // Set branch addresses
   // Float_t jtptReco[maxJets], jtphiReco[maxJets], jtetaReco[maxJets], jtptGen[maxJets], jtphiGen[maxJets], jtetaGen[maxJets],  jtmGen[maxJets],  jtmReco[maxJets], ptGen[maxJets], phiGen[maxJets], etaGen[maxJets], massGen[maxJets], ptReco[maxJets], phiReco[maxJets],etaReco[maxJets], massReco[maxJets];
    
    
    //for input 2
    Float_t jtptDataReco[maxJets], jtphiDataReco[maxJets], jtetaDataReco[maxJets],   jtmDataReco[maxJets], ptDataReco[maxJets], phiDataReco[maxJets],etaDataReco[maxJets], massDataReco[maxJets];
    
    // in qoutes is the actual branch name in the data file
    Int_t nrefDataReco, nParticleDataReco ;
    
    
    Tree5->SetBranchAddress("jtpt", jtptDataReco);
    Tree5->SetBranchAddress("jtphi", jtphiDataReco);
    Tree5->SetBranchAddress("jteta", jtetaDataReco);
    Tree5->SetBranchAddress("nref", &nrefDataReco);
    Tree5->SetBranchAddress("jtm", jtmDataReco);
    
    
    Tree6->SetBranchAddress("pt", ptDataReco);
    Tree6->SetBranchAddress("phi", phiDataReco);
    Tree6->SetBranchAddress("eta", etaDataReco);
    Tree6->SetBranchAddress("mass", massDataReco);
    Tree6->SetBranchAddress("nParticle",&nParticleDataReco);
    
    
    double fraction = 1;
    // Loop over events and match generated and reconstructed jets
    //simulation loop
    Long64_t nEntries = (Tree5->GetEntries() )*fraction;
    
    
    //progress bar stuff
    // Set it print to cout, with max progress is MaxNumber
    ProgressBar Bar(cout, nEntries);
    // Set style
    Bar.SetStyle(-1);
    
    for (Long64_t i = 0; i < nEntries; ++i) {
        // for (Long64_t i = 0; i < 1000; ++i) {
       
        Reco.clear();
        
        Particles.clear();
        
        Tree5->GetEntry(i);
        Tree6->GetEntry(i);
        // cout << " print i " << i << endl;
        
        
        //progress bar stuff
        
        // update progress to iE and print to screen
        Bar.Update(i);
        Bar.Print();
        
        
        
        
        
        
        
        
        
        
        
        // Creating reco jet vector with their pT eta phi
        for (int j = 0; j < nrefDataReco; ++j) {
            Particles.clear();
            //node vector decleration for particle tree creation
            vector<Node *> Nodesreco;
            
            //E= p^2  + m^2, we need this for energy correction
            
            double pReco= jtptDataReco[j] * cosh(jtetaDataReco[j]);
            double E_reco = sqrt(pReco*pReco + jtmDataReco[j] * jtmDataReco[j]) ;
            
            // eta = -ln tan(theta /2)
            
            double theta = 2*atan( exp (-jtetaDataReco[j]));
            
            // MJ: There is this energy correction
            // energy correction
            JEC.SetJetE(E_reco);
            JEC.SetJetTheta(theta);
            JEC.SetJetPhi(jtphiDataReco[j]);
            
            //getting the corrected energy
            double CorrectedEnergy = JEC.GetCorrectedE();
            
            double f_factor = CorrectedEnergy / E_reco ;
            
            
            // loop for particle in the gen jet
            for ( int k= 0; k< nParticleDataReco ; ++k)
                
            {
                if (jtptDataReco[j] <= 0) continue; // Skip invalid entries
                
                
                double px = jtptDataReco[j] * cos(jtphiDataReco[j]);
                double py = jtptDataReco[j] * sin(jtphiDataReco[j]);
                double pz = jtptDataReco[j] * sinh(jtetaDataReco[j]);
                double p = jtptDataReco[j] * cosh(jtetaDataReco[j]);
                double E = sqrt((p * p) + (jtmDataReco[j] * jtmDataReco[j])); // Simplified for massless particles
                // particle properties inside the jer
                double parpx = ptDataReco[k] * cos(phiDataReco[k]);
                double parpy = ptDataReco[k] * sin(phiDataReco[k]);
                double parpz = ptDataReco[k] * sinh(etaDataReco[k]);
                double parp = ptDataReco[k] * cosh(etaDataReco[k]);
                double parE = sqrt((parp * parp) + (massDataReco[k] * massDataReco[k]));
                
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
                
                
                Reco.emplace_back(Jet(jtptDataReco[j]* f_factor, jtetaDataReco[j],  jtphiDataReco[j], CorrectedEnergy , ZGreco ));
                
                
                
                
                
            }
            //fill hist
            else
                Reco.emplace_back(Jet(jtptDataReco[j]* f_factor, jtetaDataReco[j],  jtphiDataReco[j], CorrectedEnergy , -1 ));
            
            
        }
        
        //Data reconstrution
        for (int j = 0; j < nrefDataReco; ++j)  {
            Particles.clear();
            //node vector decleration for particle tree creation
            vector<Node *> Nodesreco;
            
            //E= p^2  + m^2, we need this for energy correction
            
            double pDataReco= jtptDataReco[j] * cosh(jtetaDataReco[j]);
            double E_Datareco = sqrt(pDataReco*pDataReco + jtmDataReco[j] * jtmDataReco[j]) ;
            
            // eta = -ln tan(theta /2)
            
            double theta = 2*atan( exp (-jtetaDataReco[j]));
            
            // MJ: There is this energy correction
            // energy correction
            JEC.SetJetE(E_Datareco);
            JEC.SetJetTheta(theta);
            JEC.SetJetPhi(jtphiDataReco[j]);
            
            //getting the corrected energy
            double CorrectedEnergy = JEC.GetCorrectedE();
            
            double f_factor = CorrectedEnergy / E_Datareco ;
            
            
            // loop for particle in the gen jet
            for ( int k= 0; k< nParticleDataReco ; ++k)
                
            {
                if (jtptDataReco[j] <= 0) continue; // Skip invalid entries
                
                
                double px = jtptDataReco[j] * cos(jtphiDataReco[j]);
                double py = jtptDataReco[j] * sin(jtphiDataReco[j]);
                double pz = jtptDataReco[j] * sinh(jtetaDataReco[j]);
                double p = jtptDataReco[j] * cosh(jtetaDataReco[j]);
                double E = sqrt((p * p) + (jtmDataReco[j] * jtmDataReco[j])); // Simplified for massless particles
                
                double parpx = ptDataReco[k] * cos(phiDataReco[k]);
                double parpy = ptDataReco[k] * sin(phiDataReco[k]);
                double parpz = ptDataReco[k] * sinh(etaDataReco[k]);
                double parp = ptDataReco[k] * cosh(etaDataReco[k]);
                double parE = sqrt((parp * parp) + (massDataReco[k] * massDataReco[k]));
                
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
                
                
                
                
                
                double ZGDatareco =std::min(PT1,PT2) / (PT1 + PT2);
                double EDatareco= Nodesreco[0]->P[0]; // Total energy of the jet
                // Fill histograms based on jet energy
                
                
                
                
                
                
                
                RecoData.emplace_back(Jet(jtptDataReco[j]* f_factor, jtetaDataReco[j],  jtphiDataReco[j], E_Datareco, ZGDatareco ));
                
                
                
                
                
            }
            //fill hist
            else
                RecoData.emplace_back(Jet(jtptDataReco[j]* f_factor, jtetaDataReco[j],  jtphiDataReco[j], E_Datareco, -1 ));
            
            
        }
        
        
        
        
        
        // Loop over energy bins
        
        
        
        
        
        
        // Loop over energy bins
        for (size_t j = 0; j < energyBins.size() - 1; ++j) {
            
            // Loop over all reconstructed jets from data
            for (size_t i = 0; i < RecoData.size(); ++i) {
                
                // Get the energy and ZG value of the current jet
                double EDatareco = RecoData[i].E;
                double zgReco = RecoData[i].ZG;
                
                // Skip jets where grooming failed (ZG = -1)
                if (zgReco < 0) continue;
                
                // Check if the jet energy falls into the current energy bin
                if (EDatareco >= energyBins[j] && EDatareco < energyBins[j + 1]) {
                    
                    // Fill the ZG histogram corresponding to this energy bin
                    ZGDatareco[j]->Fill(zgReco);
                }
            }
        }
        
        
        
        
    }
    
    
    // finally set it to 100% so it doesn't get stuck at 99.99%
    Bar.Update(nEntries);
    Bar.Print();
    Bar.PrintLine();
    
    
    
    
    
    
    
    
    
    
    
    
    for (size_t i = 0; i < numEnergyRanges; ++i) { // Loop over energy bins
        for (size_t j = 0; j < numBinsPerRange; ++j) { // Loop over ZG bins
            
            // Map bin j in range i to global bin k
            int globalX = j + i * numBinsPerRange;
            
            
            // Get bin content from the small histogram
            double binContent = ZGDatareco[i]->GetBinContent(j);
            cout << binContent << " GetBinContent(k)" << endl ;
            
            // Fill the big histogram
            SmearingZGDatarecoMatrix->SetBinContent(globalX, binContent);
            
            
        }
    }
    
    
    
    
    
    
    std::string outFileName = "output.root";  // Ensure the filename is valid
    
    std::cout << "Attempting to open ROOT file: " << outFileName << std::endl;
    
    TFile *outHistFile = TFile::Open(outFileName.c_str(), "RECREATE");
    
    // Check if file opened successfully
    if (!outHistFile || outHistFile->IsZombie()) {
        std::cerr << "Error: Could not create output file " << outFileName << std::endl;
        return 1;  // Avoid segfault
    }
    
    std::cout << "File opened successfully!" << std::endl;
    
    
    
    outHistFile->Close();
    
    
    
    
    
    
    TCanvas *canvas3 = new TCanvas("canvas3", "Gen Reco Energy", 1600, 1200);
    
    canvas3->Divide(2, 2);
    
    //Draw real reco data
    canvas3->cd(1);
    SmearingZGDatarecoMatrix->Draw();
    canvas3->cd(1)->SetLogz();
    
    
    // Save the canvas as a PDF
    canvas3->SaveAs("Real_Reco_Jets.pdf");
    
    
    
    
    
    
    
    input1->Close();
    delete input1;
    
    return 0;
}


