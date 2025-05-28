#include <iostream>
#include <vector>
#include <map>
#include <deque>
#include <cmath>

#define HungarianMAX 500

template <class O, class o>
void BruteForceMatchJets(double (*Metric)(O,o), std::vector<O> &A, std::vector<o> &B,
   std::vector<int> Mapping, std::vector<int> &Best, int l, int r, double &MinDistance);
template <class O, class o>
std::map<int, int> MatchJetsGreedy(double (*Metric)(O,o), std::vector<O> GenJets, std::vector<o> RecoJets);
template <class O, class o>
std::map<int, int> MatchJetsBruteForce(double (*Metric)(O,o), std::vector<O> GenJets, std::vector<o> RecoJets);
template <class O, class o>
std::map<int, int> MatchJetsHungarian(double (*Metric)(O,o), std::vector<O> GenJets, std::vector<o> RecoJets);
bool DoHungarianAssignment(int N, double Cost[HungarianMAX][HungarianMAX], int Assignment[HungarianMAX]);
bool DoHungarianSubtraction(int N, double Cost[HungarianMAX][HungarianMAX], int Assignment[HungarianMAX]);
template <class O>
void PrintVector(int N, O V[HungarianMAX]);
template <class O>
void PrintMatrix(int N, O Cost[HungarianMAX][HungarianMAX]);

template <class O, class o>
void BruteForceMatchJets(double (*Metric)(O,o), std::vector<O> &A, std::vector<o> &B,
   std::vector<int> Mapping, std::vector<int> &Best, int l, int r, double &MinDistance)
{
   if(l == r)
   {
      double SumDistance = 0;
      for(int i = 0; i < (int)A.size(); i++)
         if(Mapping[i] < (int)B.size() && Mapping[i] >= 0)
            SumDistance += (*Metric)(A[i], B[Mapping[i]]);

      if(SumDistance < MinDistance)
      {
         MinDistance = SumDistance;
         Best = Mapping;
         // cout << "Min: " << MinDistance << ", Current: " << SumDistance << endl;
      }
   }
   else
   {
      double FixedDistance = 0;
      for(int i = 0; i < (int)A.size(); i++)
         if(Mapping[i] < B.size() && Mapping[i] >= 0)
            if(i < l)
               FixedDistance += (*Metric)(A[i], B[Mapping[i]]);

      for(int i = l; i <= r; i++)
      {
         swap(Mapping[l], Mapping[i]);
         if(FixedDistance < MinDistance)
            BruteForceMatchJets<O,o>(Metric, A, B, Mapping, Best, l + 1, r, MinDistance);
         swap(Mapping[l], Mapping[i]);
      }
   }
}

template <class O, class o>
std::map<int, int> MatchJetsGreedy(double (*Metric)(O,o), std::vector<O> GenJets, std::vector<o> RecoJets)
{
   std::map<int, int> GenReco;

   for(int iG = 0; iG < (int)GenJets.size(); iG++)
   {
      int BestIndex = -1;
      double BestDistance = -1;

      for(int iR = 0; iR < (int)RecoJets.size(); iR++)
      {
         double Distance = (*Metric)(GenJets[iG], RecoJets[iR]);
         if(BestDistance < 0 || BestDistance > Distance)
         {
            BestIndex = iR;
            BestDistance = Distance;
         }
      }

      GenReco[iG] = BestIndex;
   }

   return GenReco;
}

template <class O, class o>
std::map<int, int> MatchJetsBruteForce(double (*Metric)(O,o), std::vector<O> GenJets, std::vector<o> RecoJets)
{
   std::map<int, int> GenReco;

   double SumDistance = 0, MinDistance = 1e10;
   std::vector<int> Mapping;
   std::vector<int> Best(RecoJets.size(), -1);

   for(int i = 0; i < (int)GenJets.size(); i++)
      Mapping.push_back(i);
   BruteForceMatchJets<O,o>(Metric, GenJets, RecoJets, Mapping, Best, 0, GenJets.size() - 1, MinDistance);

   for(int i = 0; i < GenJets.size(); i++)
   {
      int Target = Best[i];
      if(Best[i] > RecoJets.size())
         Target = -1;
      GenReco.insert(std::pair<int, int>(i, Target));
   }

   return GenReco;
}

template <class O, class o>
std::map<int, int> MatchJetsHungarian(double (*Metric)(O,o), std::vector<O> JetsA, std::vector<o> JetsB)
{
   // Step 0 - construct initial cost matrix
   int NA = JetsA.size();
   int NB = JetsB.size();
   int N  = std::max(NA, NB);

   double Cost[HungarianMAX][HungarianMAX];
   for(int iA = 0; iA < N; iA++)
   {
      for(int iB = 0; iB < N; iB++)
      {
         if(iA >= NA || iB >= NB)
            Cost[iA][iB] = 99999;
         else
            Cost[iA][iB] = (*Metric)(JetsA[iA], JetsB[iB]);
      }
   }

   // Step 1 - subtract rows & columns
   for(int iA = 0; iA < N; iA++)
   {
      double MinimumCost = Cost[iA][0];
      for(int iB = 1; iB < N; iB++)
         if(MinimumCost > Cost[iA][iB])
            MinimumCost = Cost[iA][iB];
      for(int iB = 0; iB < N; iB++)
         Cost[iA][iB] = Cost[iA][iB] - MinimumCost;
   }
   for(int iB = 0; iB < N; iB++)
   {
      double MinimumCost = Cost[0][iB];
      for(int iA = 1; iA < N; iA++)
         if(MinimumCost > Cost[iA][iB])
            MinimumCost = Cost[iA][iB];
      for(int iA = 0; iA < N; iA++)
         Cost[iA][iB] = Cost[iA][iB] - MinimumCost;
   }

   int Assignment[HungarianMAX] = {-1};
 
   bool Done = false;
   bool SubtractionSuccess = true;
   while(Done == false && SubtractionSuccess == true)
   {
      // Step 2 - attempt assignment
      Done = DoHungarianAssignment(N, Cost, Assignment);
      if(Done == true)
         break;

      // Step 3 - do the ticking and subtracting on the matrix
      SubtractionSuccess = DoHungarianSubtraction(N, Cost, Assignment);
   }

   std::map<int, int> GenReco;

   if(SubtractionSuccess == false)
      return GenReco;

   for(int iA = 0; iA < NA; iA++)
   {
      int iB = Assignment[iA];
      if(iB >= NB)
         iB = -1;

      GenReco[iA] = iB;
   }

   return GenReco;
}

bool DoHungarianAssignment(int N, double Cost[HungarianMAX][HungarianMAX], int AssignmentA[HungarianMAX])
{
   double Tolerance = 1e-10;
   
   int AssignmentB[HungarianMAX] = {-1};   // the inverse of AssignmentA.  This is just for convenience
   for(int i = 0; i < N; i++)
   {
      AssignmentA[i] = -1;
      AssignmentB[i] = -1;
   }

   int AZeroCount[HungarianMAX] = {0};
   int BZeroCount[HungarianMAX] = {0};

   bool Done = false;
   while(Done == false)
   {
      for(int i = 0; i < N; i++)
      {
         AZeroCount[i] = 0;
         BZeroCount[i] = 0;
      }

      for(int iA = 0; iA < N; iA++)
      {
         if(AssignmentA[iA] >= 0)
            continue;

         for(int iB = 0; iB < N; iB++)
         {
            if(AssignmentB[iB] >= 0)
               continue;

            if(Cost[iA][iB] < Tolerance)
            {
               AZeroCount[iA] = AZeroCount[iA] + 1;
               BZeroCount[iB] = BZeroCount[iB] + 1;
            }
         }
      }

      Done = true;

      int BestA = -1, BestB = -1;
      for(int iA = 0; iA < N; iA++)
      {
         if(AZeroCount[iA] == 0)
            continue;
         if(BestA < 0 || AZeroCount[BestA] > AZeroCount[iA])
            BestA = iA;
      }
      for(int iB = 0; iB < N; iB++)
      {
         if(BZeroCount[iB] == 0)
            continue;
         if(BestB < 0 || BZeroCount[BestB] > BZeroCount[iB])
            BestB = iB;
      }

      if(BestA < 0 && BestB < 0)
         break;

      // Condition for assigning on the A direction:
      //    - if nothing on B
      //    - if both A and B have something, we check if A is better
      if(BestB < 0 || AZeroCount[BestA] <= BZeroCount[BestB])
      {
         for(int iB = 0; iB < N; iB++)
         {
            if(AssignmentB[iB] >= 0)
               continue;
            if(Cost[BestA][iB] < Tolerance)
            {
               AssignmentA[BestA] = iB;
               AssignmentB[iB] = BestA;
               break;
            }
         }

         Done = false;
      }
      else if(BestA < 0 || AZeroCount[BestA] >= BZeroCount[BestB])
      {
         for(int iA = 0; iA < N; iA++)
         {
            if(AssignmentA[iA] >= 0)
               continue;
            if(Cost[iA][BestB] < Tolerance)
            {
               AssignmentA[iA] = BestB;
               AssignmentB[BestB] = iA;
               break;
            }
         }

         Done = false;
      }
      else
      {
         cerr << "This should not happen!" << endl;
      }
   }

   bool FullyAssigned = true;
   for(int iA = 0; iA < N; iA++)
      if(AssignmentA[iA] < 0)
         FullyAssigned = false;
   return FullyAssigned;
}

bool DoHungarianSubtraction(int N, double Cost[HungarianMAX][HungarianMAX], int AssignmentA[HungarianMAX])
{
   bool TickedA[HungarianMAX] = {false};
   bool TickedB[HungarianMAX] = {false};
   double Tolerance = 1e-10;

   int AssignmentB[HungarianMAX] = {-1};
   for(int iA = 0; iA < N; iA++)
      if(AssignmentA[iA] >= 0)
         AssignmentB[AssignmentA[iA]] = iA;

   std::deque<int> ProcessA;

   // Step 1: tick unassigned row
   for(int i = 0; i < N; i++)
   {
      if(AssignmentA[i] < 0)
      {
         TickedA[i] = true;
         ProcessA.push_back(i);
      }
   }

   // Step 2 & 3
   // For all (new) ticked rows, tick corresponding (new) columns with a zero entry
   // Check if the column needs a ticked row as well
   while(ProcessA.size() > 0)
   {
      int iA = ProcessA[0];
      ProcessA.pop_front();
      for(int iB = 0; iB < N; iB++)
      {
         if(TickedB[iB] == true)
            continue;
         if(Cost[iA][iB] < Tolerance)
         {
            TickedB[iB] = true;

            if(AssignmentB[iB] >= 0 && TickedA[AssignmentB[iB]] == false)
            {
               TickedA[AssignmentB[iB]] = true;
               ProcessA.push_back(AssignmentB[iB]);
            }
         }
      }
   }

   // Now find the number to subtract
   double Theta = -1;
   for(int iA = 0; iA < N; iA++)
   {
      if(TickedA[iA] == false)
         continue;
      for(int iB = 0; iB < N; iB++)
      {
         if(TickedB[iB] == true)
            continue;
         if(Theta < 0 || Theta > Cost[iA][iB])
            Theta = Cost[iA][iB];
      }
   }

   if(Theta < 0)
   {
      // std::cerr << "No need for subtraction - we can do full assignment" << std::endl;
      return false;
   }

   // Do the subtraction
   for(int iA = 0; iA < N; iA++)
   {
      for(int iB = 0; iB < N; iB++)
      {
         if(TickedB[iB] == false && TickedA[iA] == true)   // no line pass through
            Cost[iA][iB] = Cost[iA][iB] - Theta;
         if(TickedB[iB] == true && TickedA[iA] == false)   // two lines pass through
            Cost[iA][iB] = Cost[iA][iB] + Theta;
      }
   }

   return true;
}

template <class O>
void PrintVector(int N, O V[HungarianMAX])
{
   std::cout << "(";
   for(int i = 0; i < N; i++)
   {
      if(i != 0)
         std::cout << " ";
      std::cout << V[i];
   }
   std::cout << ")" << endl;
}

template <class O>
void PrintMatrix(int N, O Cost[HungarianMAX][HungarianMAX])
{
   for(int iA = 0; iA < N; iA++)
   {
      for(int iB = 0; iB < N; iB++)
         std::cout << " " << Cost[iA][iB];
      std::cout << endl;
   }
   std::cout << endl;
}

