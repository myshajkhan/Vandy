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
#include "TStyle.h"  // Required for gStyle

#include "JetCorrector.h"
#include "Matching.h"
#include "SoftDropGrooming.h"

#define maxJets 2400

#include "ProgressBar.h"

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
    vector<Jet> RecoData;// Vector to store reconstructed jets from actual experimental data (not simulation)
    vector<FourVector> Particles;// Vector to temporarily store four-momentum vectors (e.g., of jet constituents) for grooming
    vector<TH1F*> Angle;// Vector of 1D histograms to record the angular distance between matched gen and reco jets
    vector <TH1F*> ZGDatareco; // Vector of 1D histograms for ZG distributions from reconstructed data jets (real data)
    vector < TH2F*> correlation;  // Vector of 2D histograms to record correlations — e.g., energy or ZG between Gen and Reco jets
    vector<double> energyBins = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55}; // Define energy bin edges for histogramming or bin-based analysis (12 bins from 0 to 55 in steps of 5)
    int numBinsPerRange = 10; // Assuming 100 bins per ZG axis
    
    // Open the input file
    TFile *input1 = TFile::Open("./LEP1MC1994_recons_aftercut-002.root", "READ");
    TFile *input2 = TFile::Open("./LEP1Data1994P3_recons_aftercut-MERGED.root", "READ");
    TFile *input3 = TFile::Open("./Matched_jets.root", "READ");
    
    if (!input1 || input1->IsZombie()||!input2 || input2->IsZombie()||!input3 || input3->IsZombie()) {
        cerr << "Error: Could not open input file." << endl;
        return 1;
    }
    
    // Reading the text files for enegy correction
    JetCorrector JEC({"JEC_EEAK4_MC_20210513.txt", "JEC_EEAK4_DataL2_20210514.txt", "JEC_EEAK4_DataL3_20210514.txt"});
    
    //Define one ZG histogram per reconstructed energy bin range (12 total), to analyze ZG vs Data Reco E
    for (size_t j = 0; j < energyBins.size() - 1; ++j) {
        TString name = Form("ZGDatareco_%lu_%lu", j);
        TString title = Form("%g < DataReco E < %g ");
        ZGDatareco.push_back(new TH1F(name, title, numBinsPerRange, 0.05, .55)); // ZG on both axes
        ZGDatareco[j]->Sumw2(); //Turning in uncertainty measurment in hists
    }
    
    //Total energy bins
    int numEnergyRanges = energyBins.size() - 1;
    
    // Flattened 1D histogram to store ZG distributions across all energy bins
    TH1F *SmearingZGDatarecoMatrix = new TH1F("SmearZG_Datareco", "All Data Reco ZG Engergy Matrix ", numBinsPerRange * numEnergyRanges, 0, numBinsPerRange * numEnergyRanges);
    SmearingZGDatarecoMatrix->Sumw2();// Turning on uncertainty bar
    TH1F *SmearFaketoallRecoJetRatio = (TH1F*)input3->Get("Fake_Reco_Percentage_Matrix"); //Getting histogram that has fake reconstructed jet percentage
    SmearFaketoallRecoJetRatio->Sumw2(); //Getting error bars
    //Checkig if the histogram exists in the inut file
    if (!SmearFaketoallRecoJetRatio) {
        std::cerr << "ERROR: Could not find Fake_Reco_Percentage_Matrix in input3!" << std::endl;
        return 1;
    }
    
    SmearFaketoallRecoJetRatio->Sumw2(); //Turning on unceratainty bar
    TH1F *SelectedRealZGDatarecoMatrix =new TH1F("Selected_Real_ZG_Datareco", "Real Reco Data Jet ZG Engergy Matrix", numBinsPerRange * numEnergyRanges, 0, numBinsPerRange * numEnergyRanges); //Defining Histogram of reco jets after subtracting MC-estimated fake contribution
    SelectedRealZGDatarecoMatrix->Sumw2(); //Turning on unceratainty bar
    
    //Acessing the trees in the root file and renaming them
    TTree *Tree5 = (TTree *)input2->Get("akR4ESchemeJetTree"); // Reconstructed jet tree from real experimental data
    TTree *Tree6= (TTree *)input2->Get("t");  // Particle-level information for reconstructed data jets (constituents from real data)
    
    // Arrays for jet-level quantities for jet and particles from input 1
    Float_t jtptDataReco[maxJets], jtphiDataReco[maxJets], jtetaDataReco[maxJets],   jtmDataReco[maxJets]; //Data Jet level arrays for reco jets
    Float_t ptDataReco[maxJets], phiDataReco[maxJets],etaDataReco[maxJets], massDataReco[maxJets]; //Data particle level arrays for reco jets
    
    // Defining the single varibles- total particle number and jet number
    Int_t nrefDataReco, nParticleDataReco ;
    
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
    
    double fraction = 1;
    
    // Loop over events and match generated and reconstructed jets
    Long64_t nEntries = (Tree5->GetEntries() )*fraction;
    
    //Progress bar
    ProgressBar Bar(cout, nEntries); // Set it print to cout, with max progress is MaxNumber
    Bar.SetStyle(6);  // Set style
    
    // Loop over energy bins and fill corresponding ZG histograms for all reconstructed jets
    for (Long64_t i = 0; i < nEntries; ++i) {
        //  for (Long64_t i = 0; i < 5000; ++i) {
        Tree5->GetEntry(i); //Accessing reconstructed jet tree from real experimental data
        Tree6->GetEntry(i); //Accessing reconstructed particle data from real experimental data
        
        RecoData.clear(); //Cleaning RecoData vector for a new loop run
        Particles.clear();//Cleaning Particles vector for a new loop run
        
        //Update progress to iE and print to screen
        Bar.Update(i);
        if (i%1000 == 0) { //show progress
            Bar.Print();
        }
        
        // Creating Reco jet vector with their pT eta phi
        for (int j = 0; j < nrefDataReco; ++j) {
            
            //node vector decleration for particle tree creation
            vector<Node *> Nodesreco;
            
            //Clearning Partice vector at the start of each jet loop
            Particles.clear();
            
            double pReco= jtptDataReco[j] * cosh(jtetaDataReco[j]); // Compute reconstructed jet momentum magnitude using pT and eta
            double E_reco = sqrt(pReco*pReco + jtmDataReco[j] * jtmDataReco[j]) ;// Compute total energy of the reconstructed jet using momentum and mass
            double theta = 2*atan( exp (-jtetaDataReco[j])); // Calculate jet polar angle theta from pseudorapidity (eta)
            
            //Apply energy correction by setting jet energy, angle, and phi for JetCorrector
            JEC.SetJetE(E_reco); // Set jet energy for correction
            JEC.SetJetTheta(theta); // Set jet polar angle (theta) for correction
            JEC.SetJetPhi(jtphiDataReco[j]); // Set jet azimuthal angle (phi) for correction
            
            //Getting the corrected energy
            double CorrectedEnergy = JEC.GetCorrectedE();
            double f_factor = CorrectedEnergy / E_reco ;
            
            double px = jtptDataReco[j] * cos(jtphiDataReco[j]);
            double py = jtptDataReco[j] * sin(jtphiDataReco[j]);
            double pz = jtptDataReco[j] * sinh(jtetaDataReco[j]);
            double p = jtptDataReco[j] * cosh(jtetaDataReco[j]);
            double E = sqrt((p * p) + (jtmDataReco[j] * jtmDataReco[j])); // Simplified for massless particles
            // Loop over the particles in the Reco Data jets
            for ( int k= 0; k< nParticleDataReco ; ++k)
            {
                
                if (jtptDataReco[j] <= 0) continue; // Skip invalid entries
                
                // particle properties inside the jer
                double parpx = ptDataReco[k] * cos(phiDataReco[k]);
                double parpy = ptDataReco[k] * sin(phiDataReco[k]);
                double parpz = ptDataReco[k] * sinh(etaDataReco[k]);
                double parp = ptDataReco[k] * cosh(etaDataReco[k]);
                double parE = sqrt((parp * parp) + (massDataReco[k] * massDataReco[k]));
                
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
                // Grooming
                SDNodereco =   FindSDNodeE(Nodesreco[0], .1 , 0, 0.4 );
            }
            
            if (SDNodereco && SDNodereco->Child1 && SDNodereco->Child2) {
                double PT1 = SDNodereco->Child1->P[0];
                double PT2 = SDNodereco->Child2->P[0];
                double ZGreco =std::min(PT1,PT2) / (PT1 + PT2);
                double Ereco= Nodesreco[0]->P[0]; // Total energy of the jet
                
                RecoData.emplace_back(Jet(jtptDataReco[j]* f_factor, jtetaDataReco[j],  jtphiDataReco[j], CorrectedEnergy , ZGreco ));
                
            }
            //fill hist
            else
                RecoData.emplace_back(Jet(jtptDataReco[j]* f_factor, jtetaDataReco[j],  jtphiDataReco[j], CorrectedEnergy , -1 ));
        }
        
        for (size_t n = 0; n < energyBins.size() - 1; ++n) {
            
            // Loop over all reconstructed jets from data
            for (size_t m = 0; m < RecoData.size(); ++m) {
                
                // Get the energy and ZG value of the current jet
                double EDatareco = RecoData[m].E;
                double zgReco = RecoData[m].ZG;
                
                // Skip jets where grooming failed (ZG = -1)
                if (zgReco < 0) continue;
                
                // Check if the jet energy falls into the current energy bin
                if (EDatareco >= energyBins[n] && EDatareco < energyBins[n + 1]) {
                    
                    // Fill the ZG histogram corresponding to this energy bin
                    ZGDatareco[n]->Fill(zgReco);
                }
            }
        }
    }
    
    // finally set it to 100% so it doesn't get stuck at 99.99%
    Bar.Update(nEntries);
    Bar.Print();
    Bar.PrintLine();
    
    //Filling up 1D smearing matrix ( ZG content with energy bin) for all real data
    for (size_t i = 0; i < numEnergyRanges; ++i) { // Loop over energy bins
        for (size_t j = 0; j < numBinsPerRange; ++j) { // Loop over ZG bins
            
            // Map bin j in range i to global bin k
            int globalX =1+ j + i * numBinsPerRange;
            // Get bin content from the small histogram
            double binContent = ZGDatareco[i]->GetBinContent(j+1);
            double binUncertainty = ZGDatareco[i]->GetBinError(j+1);
            // Fill the big histogram
            SmearingZGDatarecoMatrix->SetBinContent(globalX, binContent);
            //Copy the ucertainty of the small hist
            SmearingZGDatarecoMatrix->SetBinError(globalX, binUncertainty);
        }
    }
    
    //Filling up 1D smearing matrix ( ZG content with energy bin) for all real data
    for (size_t i = 0; i < numEnergyRanges; ++i) { // Loop over energy bins
        for (size_t j = 0; j < numBinsPerRange; ++j) { // Loop over ZG bins
            
            // Map bin j in range i to global bin k
            int globalX =1+ j + i * numBinsPerRange;
            
            // Get bin content of fakepercentage and reweight
            double binContent = SmearFaketoallRecoJetRatio->GetBinContent(globalX);
            double all = SmearingZGDatarecoMatrix->GetBinContent(globalX);
            double reweightedfake = binContent * all; //mutiplying with real data bin content
            double RealPortionData =  all - reweightedfake;
            
            //Error propagation
            double FakepercentErr    = SmearFaketoallRecoJetRatio->GetBinError(globalX);    // error on fraction
            double AlldataErr   = SmearingZGDatarecoMatrix->GetBinError(globalX);      // error on data counts
            
            // Error propagation: σ² = ( (1-f)² * σD² ) + ( D² * σf² )
            double RealPortionErr = sqrt( ((1.0 - binContent) * (1.0 - binContent) * (AlldataErr * AlldataErr))
                                                + (all * all * FakepercentErr * FakepercentErr) );
                   
            
            // Fill the big histogram
            SelectedRealZGDatarecoMatrix->SetBinContent(globalX, RealPortionData);
            SelectedRealZGDatarecoMatrix->SetBinError(globalX,RealPortionErr);
        }
    }
    std::string outFileName = "Real_data.root";  // Ensure the filename is valid
    cout << "Attempting to open ROOT file: " << outFileName << std::endl;
    
    TFile *outHistFile = TFile::Open(outFileName.c_str(), "RECREATE");
    
    // Check if file opened successfully
    if (!outHistFile || outHistFile->IsZombie()) {
        std::cerr << "Error: Could not create output file " << outFileName << std::endl;
        return 1;  // Avoid segfault
    }
    
    // Check before writing to avoid segfault
    if (SmearingZGDatarecoMatrix) {
        SmearingZGDatarecoMatrix->Write();
    } else {
        std::cerr << "Error: SmearingZGDatarecoMatrix is nullptr!" << std::endl;
    }
    std::cout << "File opened successfully!" << std::endl;
    
    outHistFile->Close();
    
    TCanvas *canvas3 = new TCanvas("canvas3", "Gen Reco Energy", 1600, 1200);
    
    canvas3->Divide(2, 2);
    
    //Draw real reco data
    canvas3->cd(1);
    SmearingZGDatarecoMatrix->Draw();
    canvas3->cd(1)->SetLogz();
    
    canvas3->cd(2);
    SelectedRealZGDatarecoMatrix->Draw();
    
    canvas3->cd(3);
    SmearFaketoallRecoJetRatio->Draw();
    // Save the canvas as a PDF
    canvas3->SaveAs("Real_Reco_Jets.pdf");
    
    input1->Close();
    delete input1;
    
    return 0;
}


