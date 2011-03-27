#include <TLorentzVector.h>


float findTheMini(float , float );
float findTheMaxi(float , float );
double DeltaR(double , double , double , double );
int* findLeadAndTrail(TLorentzVector*, int);
bool EventPassDiphotonFilter(TLorentzVector*,  int, double , double , int* , int* );
double CosThetaStar(TLorentzVector, TLorentzVector);
