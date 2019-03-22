#include <cmath>
#include <map>
#include <iostream>
#include "TLorentzVector.h"


using namespace std;

float deltaPHI(float phi1, float phi2) {
    float result = phi1 - phi2;
    while (result > float(M_PI)) result -= float(2*M_PI);
    while (result <= -float(M_PI)) result += float(2*M_PI);
    return result;
}

float dphi_2(float pt1, float eta1, float phi1, float pt2, float eta2, float phi2, int ind) {
    TLorentzVector p41, p42;
    p41.SetPtEtaPhiM(pt1,eta1,phi1,0.);
    p42.SetPtEtaPhiM(pt2,eta2,phi2,0.);
    float dphi = -1.;
    dphi = deltaPHI( (p41+p42).Phi() , ind == 1 ? phi1 : phi2 );
//    std::cout<<dphi<<std::endl;
    return dphi;
}

