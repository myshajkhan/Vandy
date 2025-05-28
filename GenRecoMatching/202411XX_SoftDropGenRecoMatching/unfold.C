#include <iostream>
#include <vector>
using namespace std;

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TError.h"
#include "TSpline.h"
#include "TGraph.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "RooUnfold.h"
#include "RooUnfoldInvert.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
// #include "TUnfold.h"
#include "TUnfoldDensity.h"

#include "RootUtilities.h"
#include "CommandLine.h"
#include "CustomAssert.h"


int main()


{
    
    
    //get input funtions
    
    // declare unfolds
    
    bool DoBayes            = CL.GetBool("DoBayes",       true); //basianunfold
    bool DoSVD              = CL.GetBool("DoSVD",         true); //single value decomposition
    
    // declare historgams
    
    
    //remove out of range stuff
    RemoveOutOfRange(HMeasured);
    RemoveOutOfRange(HTruth);
    RemoveOutOfRange(HResponse);
    RemoveOutOfRange(HRawResponse);
    RemoveOutOfRange(HInput);

    
    //bullian loop
    if(DoBayes == true)
      {
         // vector<int> Iterations{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 80, 90, 100, 125, 150, 200, 250};
         vector<int> Iterations{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 35, 40, 45, 50, 55, 60, 65, 70, 80, 90, 100, 125, 150, 200, 250};
         for(int I : Iterations)
         {
            RooUnfoldBayes BayesUnfold(Response, HInput, I);
            BayesUnfold.SetVerbose(-1);
            HUnfolded.push_back((TH1 *)(BayesUnfold.Hreco(ErrorChoice)->Clone(Form("HUnfoldedBayes%d", I))));
            Covariance.insert(pair<string, TMatrixD>(Form("MUnfoldedBayes%d", I), BayesUnfold.Ereco()));
            TH1D *HFold = ForwardFold(HUnfolded[HUnfolded.size()-1], HResponse);
            HFold->SetName(Form("HRefoldedBayes%d", I));
            HRefolded.push_back(HFold);
         }
      }
    
    
    
    //sdv loop
    if(DoSVD == true)
       {
          vector<double> SVDRegularization{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 80, 90, 100, 125, 150};
          // vector<double> SVDRegularization{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 22, 24, 26, 28, 30, 35, 40, 45, 50, 55, 60, 65, 70, 80, 90, 100, 125, 150, 200, 250};
          for(double D : SVDRegularization)
          {
             if(D >= HGen->GetNbinsX())
                continue;

             RooUnfoldSvd SVDUnfold(Response, HInput, D);
             SVDUnfold.SetVerbose(-1);
             HUnfolded.push_back((TH1 *)(SVDUnfold.Hreco(ErrorChoice)->Clone(Form("HUnfoldedSVD%.1f", D))));
             Covariance.insert(pair<string, TMatrixD>(Form("MUnfoldedSVD%.1f", D), SVDUnfold.Ereco()));
             TH1D *HFold = ForwardFold(HUnfolded[HUnfolded.size()-1], HResponse);
             HFold->SetName(Form("HRefoldedSVD%.1f", D));
             HRefolded.push_back(HFold);
          }
       }
    
    
   //plot
