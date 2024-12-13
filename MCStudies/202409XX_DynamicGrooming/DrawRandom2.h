//----------------------------------------------------------------------------
#ifndef DrawRandom2_H_AJGRKJVCGTZRHISCGGARHWGCA
#define DrawRandom2_H_AJGRKJVCGTZRHISCGGARHWGCA
//----------------------------------------------------------------------------
#include <cmath>
#include <cstdlib>
//----------------------------------------------------------------------------
#define PI 3.14159265358979323846264338327950288479716939937510
//----------------------------------------------------------------------------
class RandomBase;
class RandomMT;
//----------------------------------------------------------------------------
class RandomBase
{
public:
   virtual double DrawRandom();
   virtual double DrawRandom(double max);
   virtual double DrawRandom(double min, double max);
   double DrawSine(double min, double max);
   double DrawLorentzian(double center, double gamma);
   double DrawGaussian(double center, double sigma);
   double DrawGaussian(double sigma);
   double DrawTruncatedGaussian(double center, double sigma, double min, double max);
   double DrawTruncatedGaussian(double sigma, double min, double max);
   double DrawTruncatedGaussian(double min, double max);
   double DrawGaussianBoxMuller();
   double DrawCruijff(double center, double sigmal, double sigmar, double alphal, double alphar);
   double DrawCruijff(double center, double sigmal, double sigmar, double alphal, double alphar, double left, double right);
   double DrawExponential(double exponent, double left, double right);
   double DrawExponential(double exponent, double side);
   double DrawPoisson(double mean);
   double DrawPoissonDouble(double mean);
   int DrawPoissonInt(double Lambda);
   double DrawDoubleSidedCBShape(double Mean, double Sigma, double AlphaL, double AlphaR, double NL, double NR, double NormalizationL = -1, double NormalizationM = -1, double NormalizationR = -1);
   double DrawDoubleSidedCBShape(double AlphaL, double AlphaR, double NL, double NR);
   double DrawDoubleSidedCBShapeWithNormalization(double AlphaL, double AlphaR, double NL, double NR, double NormalizationL = -1, double NormalizationM = -1, double NormalizationR = -1);
   double DrawInverse(double min, double max);
   double CachedExp(double X);
   double CachedErf(double X);
};
//----------------------------------------------------------------------------
class RandomMT : public RandomBase
{
private:
   long long MT[624];
   int Index;
public:
   RandomMT(int Seed = 31426);
private:
   void GenerateNumbers();
public:
   long long DrawRandomInteger();
   double DrawRandom();
   double DrawRandom(double max);
   double DrawRandom(double min, double max);
};
//----------------------------------------------------------------------------
#endif

