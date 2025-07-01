#ifndef SOFT_DROP_GROOMING_H
#define SOFT_DROP_GROOMING_H

#include "CATree.h"

//start of lambda funtion
inline Node* SoftDrop(Node *HeadNode, double ZCut, double Beta, double R0){
    if(HeadNode == NULL)
        return NULL;
    bool Done = false;
    Node *Current = HeadNode;
    
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
            double PRatio = std::min(P1, P2) / (P1 + P2);
            
            double Angle = GetAngle(Current->Child1->P, Current->Child2->P);
            
            double Threshold = ZCut * std::pow(Angle / R0, Beta);
            
            if(PRatio > Threshold)
                Done = true;
            else
            {
                if(P1 > P2)
                    Current = Current->Child1;
                else
                    Current = Current->Child2;
            }
        }
    }
    
    return Current;
};

#endif
