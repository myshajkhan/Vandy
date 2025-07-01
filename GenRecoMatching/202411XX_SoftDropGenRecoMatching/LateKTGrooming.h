#ifndef LATEKT_GROOMING_H
#define LATEKT_GROOMING_H

#include "CATree.h"

//Late KT starts from the first genration and goes to the next genration while checking for the strongest split.
inline Node* LateKTGrooming(Node *HeadNode, double ZCut, double alpha, double R0) {
    if(HeadNode == NULL)
        return NULL;
    
    //Function will be searching for the hardest split — it'll only stop when it reaches the last particle or finds a good enough split
    bool Done = false;
    
    //Start at the starting of the jet (root node) and assume it's the hardest split/best node for now
    Node *Current = HeadNode;
    Node *bestNode = HeadNode;
    
    double  bestkT= 0;
    double jetPT= HeadNode->P[0]; //We know P[0] = pT. Total jet momentum is the first cluster momentum
    
    //Starting of the while loop for Late KT grooming method structure
    while(Done == false)   {
        
        if(Current->N == 1) //N is the number of particles. We are at the end of the tree and its no longer splitting and we should stop the function
            Done = true;
        else if(Current->N == 2) //Only two particles were clustered, and no other particles were close enough (in ΔR) to join them.
        {
            std::cerr << "Error!  N = " << Current->N << "!" << std::endl;
        }
        else if(Current->Child1 == NULL || Current->Child2 == NULL) //This node was supposed to split, but one or both children are missing — something’s broken, so we skip it
        {
            std::cerr << "Error!  Child NULL while N = " << Current->N << "!" << std::endl;
        }
        else
        {
            
            // If N == 3, that means this mom was made by merging two kids: one with 2 particles and one with 1. Child1 and Child2 are still just two branches — they each might be made of more stuff underneath.
            double P1 = Current->Child1->P[0];  //We know P[0] = pT. Child 1's transverse momentum
            double P2 = Current->Child2->P[0]; //Child 2's transverse momentum
            double Angle = GetAngle(Current->Child1->P, Current->Child2->P); //Angle between child 1 and child 2
            
            // To Compute azimuthal angle difference of the two child cluters we must define the children first with all the information of momentum. Azimuthal is the angle a particle makes around the beam axis (z-axis), measured in the transverse plane (x–y plane).
            TLorentzVector v1, v2;
            
            v1.SetPxPyPzE(Current->Child1->P[1], Current->Child1->P[2], Current->Child1->P[3], Current->Child1->P[0]); // Child 1's momentum(x,y,z) and energy(PT)
            v2.SetPxPyPzE(Current->Child2->P[1], Current->Child2->P[2], Current->Child2->P[3], Current->Child2->P[0]); // Child 2's momentum(x,y,z) and energy(PT)
            
            double azimuth1 = v1.Phi(); //Child 1's azimuthal angle
            double azimuth2 = v2.Phi(); //Child 2's azimuthal angle
            double azimuth_diff = (azimuth1 - azimuth2) * (azimuth1 - azimuth2); //square of their angular difference
            
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
            
            //Take the path of the child with the highest momentum
            if(P1 > P2)
                Current = Current->Child1;
            else
                Current = Current->Child2;
        }
    }
    return bestNode;
};
#endif
