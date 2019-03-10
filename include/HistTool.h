#ifndef _hist_tool_h_
#define _hist_tool_h_

#include "TGraph.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"

#include <map>
#include <vector>

class HistTool {

 public:
  HistTool() {}
  ~HistTool() {}
  
  inline void addHist(TString title, TH1* hist) {
    if(hist->InheritsFrom("TH2")) {
      all2dPlots_[title]=(TH2 *)hist;
    }
    else {
      allPlots_[title] = hist;
    }
  }

  inline void fill(TString title, double value, double weight,std::vector<TString> cats){
    for(auto &c : cats)
      fill(title,value,weight,c);
  }

  inline void fill(TString title, double value, double weight,TString cat="") {

    if (not allPlots_.count(title)) {
      //std::cout << "Histogram " << title << " not registered, not filling." << std::endl;
      return;
    }

    if(allPlots_[title]->InheritsFrom("TH2")) return;
    
    //category specific plot, init if needed
    if(cat!=""){
      TString newTitle=cat+"_"+title;
      if(not allPlots_.count(newTitle)) {
        //std::cout << "Histogram " << title << " for cat=" << cat << " not yet started, adding now." << std::endl;
        allPlots_[newTitle]=(TH1 *)allPlots_[title]->Clone(newTitle);
        allPlots_[newTitle]->SetDirectory(0);
        allPlots_[newTitle]->Reset("ICE");
      }
      title=newTitle;
    }
    
    allPlots_[title]->Fill(value, weight);
  }

  inline void fill2D(TString title, double valueX, double valueY, double weight,std::vector<TString> cats){
    for(auto &c : cats)
      fill2D(title,valueX,valueY,weight,c);
  }

  inline void fill2D(TString title, double valueX, double valueY, double weight,TString cat="") {
    if (not all2dPlots_.count(title)) {
      std::cout << "Histogram " << title << " not registered, not filling." << std::endl;
      return;
    }
    
    //category specific plot, init if needed
    if(cat!=""){
      TString newTitle=cat+"_"+title;
      if(not all2dPlots_.count(newTitle)) {
        //std::cout << "Histogram " << title << " for cat=" << cat << " not yet started, adding now." << std::endl;
        all2dPlots_[newTitle]=(TH2 *)all2dPlots_[title]->Clone(newTitle);
        all2dPlots_[newTitle]->SetDirectory(0);
        all2dPlots_[newTitle]->Reset("ICE");
      }
      title=newTitle;
    }
    
    all2dPlots_[title]->Fill(valueX,valueY, weight);
  }
  
  std::map<TString, TH1 *> &getPlots()   { return allPlots_; }
  std::map<TString, TH2 *> &get2dPlots() { return all2dPlots_; }
  
 private:
  std::map<TString, TH1 *> allPlots_;
  std::map<TString, TH2 *> all2dPlots_;

};

#endif
