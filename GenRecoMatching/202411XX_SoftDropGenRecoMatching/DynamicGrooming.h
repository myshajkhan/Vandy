#ifndef DYNAMIC_GROOMING_H
#define DYNAMIC_GROOMING_H

#include "CATree.h"

//Dynamic grooming funtion starts from the first genration and goes to the next genration while checking for the strongest split.
inline Node* DynamicGrooming(Node *HeadNode, double ZCut, double alpha , double R0 ) {
    if(HeadNode == NULL) // If there’s no full jet to start from, the function won’t do anything
        return NULL;
    
    //Function will be searching for the hardest split — it'll only stop when it reaches the last particle or finds a good enough split
    bool Done = false;
    
    //Start at the starting of the jet (root node) and assume it's the hardest split/best node for now
    Node *Current = HeadNode;
    Node *bestNode = HeadNode;
    
    double besthard = 0;
    double jetPT= HeadNode->P[0]; //Total jet momentum is the first cluster momentum
    
    while(Done == false)
    {
        if(Current->N == 1) //N is the number of particles. We are at the end of the tree and its no longer splitting and we should stop the function
            Done = true;
        else if(Current->N == 2) //Only two particles were clustered, and no other particles were close enough (in ΔR) to join them.
        {
            std::cerr << "Error!  N = " << Current->N << "!" << std::endl;
        }
        else if(Current->Child1 == NULL || Current->Child2 == NULL) // This node was supposed to split, but one or both children are missing — something’s broken, so we skip it
        {
            std::cerr << "Error!  Child NULL while N = " << Current->N << "!" << std::endl;
        }
        
        // If N == 3, that means this mom was made by merging two kids: one with 2 particles and one with 1. Child1 and Child2 are still just two branches — they each might be made of more stuff underneath.
        else
        {
            double P1 = Current->Child1->P[0]; //We know P[0] = pT. Child 1's transverse momentum
            double P2 = Current->Child2->P[0]; //Child 2's transverse momentum
            double PRatio = std::min(P1, P2) / (P1 + P2); //This gives us momentum sharing symmetry. We don’t care which branch is bigger, just how uneven the split is
            double Angle = GetAngle(Current->Child1->P, Current->Child2->P);  //Angle between two children, this is theta(i) for dynamic grooming
            double hard= PRatio*(1 -PRatio) *(P1+P2)*pow(Angle/R0 , alpha ); // Hardness formula — A numbring system captures how balanced, energetic, and wide the split is
            
            //If this split is the hardest one we’ve seen so far, save it.
            if(hard > besthard)
            {
                besthard = hard;
                bestNode = Current;//We got the children from the current node. It will pick any of the chilren looping back. If this hardness is bigger, we save this exact node as bestNode. Keep the node that gave us the largest hardess number
            }
            
            //If momentum is not balanced, take the bigger momentum child's node as our current node.
            if(P1 > P2)
                Current = Current->Child1;
            else
                Current = Current->Child2;
        }
    }
    return bestNode;
};

#endif
