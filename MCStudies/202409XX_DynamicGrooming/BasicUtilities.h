#include <map>
#include <vector>
#include <iostream>

#include "TauHelperFunctions3.h"

#define PI 3.14159265358979323846264338327950288479716939937510

class XY;
class Matrix22;
double GetRho(std::vector<double> *EtaMax, std::vector<double> *Rho, double Eta, bool SuppressMessage = false);
double GetCentrality(int Bin);
bool IsDataFromTag(std::string Tag);
bool IsPPFromTag(std::string Tag);
bool IsPPHiRecoFromTag(std::string Tag);
bool IsPAFromTag(std::string Tag);
template<class T, class U>
   std::vector<int> DRMatching(int BaseN, T BaseEta, T BasePhi,
      int TargetN, U TargetEta, U TargetPhi);
int GetIndex(double X, double *Begin, double *End);
double GetDR2(double Eta1, double Phi1, double Eta2, double Phi2);
double GetDR(double Eta1, double Phi1, double Eta2, double Phi2);

class XY
{
public: 
   double X;
   double Y;
public:
   XY();
   XY(double x, double y, bool DoPrimary = false);
   XY(FourVector J, bool DoPrimary = false);
   XY(Matrix22 M, XY &C, FourVector J);
   XY operator +(XY other);
   XY operator -(XY other);
   XY operator -();
   XY operator *(double c);
   XY operator /(double c);
   double Norm();
   double Norm2();
   void SetAsAverage(FourVector J1, FourVector J2);
   void PrimaryY();
   ~XY();
};

class Matrix22
{
public:
   double T[2][2];
public:
   Matrix22();
   Matrix22(double T11, double T12, double T21, double T22);
   Matrix22(FourVector J1, FourVector J2);
   Matrix22 &operator =(Matrix22 &other);
   double *operator [](int Index);
   XY operator *(XY A);
   Matrix22 Inverse();
   void Initialize(FourVector J1, FourVector J2);
   ~Matrix22();
};



