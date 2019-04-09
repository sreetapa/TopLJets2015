#ifndef LeptonSummary_h
#define LeptonSummary_h

#include "TLorentzVector.h"

class LeptonSummary {

 public:
 LeptonSummary(int _id,TLorentzVector &_p4) : 
  id(_id), charge(0), p4(_p4),
    chiso(0.), nhiso(0.), phoiso(0.), rho(0.), chrho(0.), nhrho(0.), phorho(0.),isofull(0.), d0(0.), d0err(0.), dz(0.),origIdx(-1) { }
  LeptonSummary(const LeptonSummary &l) {
    id=l.id;
    p4=l.p4;
    charge=l.charge;
    chiso=l.chiso;
    nhiso=l.nhiso;
    phoiso=l.phoiso;
    rho=l.rho;
    chrho=l.chrho;
    nhrho=l.nhrho;
    phorho=l.phorho;
    isofull=l.isofull;
    d0=l.d0;
    d0err=l.d0err;
    dz = l.dz;
    origIdx=l.origIdx;
  }
  int id;
  int charge;
  TLorentzVector p4;
  float chiso, nhiso, phoiso;
  float rho, chrho, nhrho, phorho;
  float isofull, d0, d0err, dz;
  int origIdx;
};

#endif
