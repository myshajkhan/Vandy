
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

//start of lambda funtion- we use dynamic grooming on this part
auto Find_late_kt_Lambda = [](Node *HeadNode, double ZCut, double alpha, double R0) -> Node* {
    if(HeadNode == NULL)
        return NULL;
    bool Done = false;
    Node *Current = HeadNode;
    Node *bestNode = HeadNode;
    double  bestkT= 0;
    double jetPT= HeadNode->P[0];
 
    
    while(Done == false)   {

      
 
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
            
           
            // Compute azimuthal angle difference (phi1 - phi2)
                        TLorentzVector v1, v2;
                        v1.SetPxPyPzE(Current->Child1->P[1], Current->Child1->P[2], Current->Child1->P[3], Current->Child1->P[0]);
                        v2.SetPxPyPzE(Current->Child2->P[1], Current->Child2->P[2], Current->Child2->P[3], Current->Child2->P[0]);

                        double azimuth1 = v1.Phi();
                        double azimuth2 = v2.Phi();
                        double azimuth_diff = (azimuth1 - azimuth2) * (azimuth1 - azimuth2);

                        // Compute rapidity (y) for each subjet
                        double y1 = 0.5 * log((Current->Child1->P[0] + Current->Child1->P[3]) / (Current->Child1->P[0] - Current->Child1->P[3]));
                        double y2 = 0.5 * log((Current->Child2->P[0] + Current->Child2->P[3]) / (Current->Child2->P[0] - Current->Child2->P[3]));
                        double y_diff = (y1 - y2) * (y1 - y2);

                        // Calculate delta squared (∆R^2 = ∆y^2 + ∆phi^2)
                        double delta_sq = y_diff + azimuth_diff;

                        // Calculate kT = min(PT1, PT2) * sqrt(∆R^2)
                        double kT = std::min(P1, P2) * sqrt(delta_sq);

                        // Check if this kT is the largest we've found so far
            if (kT > bestkT) {
                bestkT = kT;  // Update the best kT
                bestNode = Current;  // Update the best node
            }
            

            
            if(P1 > P2)
                Current = Current->Child1;
            else
                Current = Current->Child2;
           }
          //  else
           // {
        

 
            //}
     
}
   
    
    return bestNode;
};


//end of lambda funtion

//start of lambda funtion- we use dynamic grooming on this part
auto Find_dynamic_Lambda = [](Node *HeadNode, double ZCut, double alpha, double R0) -> Node* {
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
           
            std::cerr << "Error!  N = " << Current->N << "!" << std::endl;
        }
        else if(Current->Child1 == NULL || Current->Child2 == NULL)
        {
           
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



