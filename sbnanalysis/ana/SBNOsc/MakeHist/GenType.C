#include "TFile.h"
#include "TTree.h"
#include <vector>
#include <string>
#include <algorithm>

void GenType() {
  TFile *myFile = TFile::Open("output_SBNOsc_NueSelection.root");
  if (myFile==0){
    printf("File not correctly opened!\n");
    return;
  }

  TH1D* gen_nue = (TH1D*)myFile->Get("generated_nue_hist");
  TH1D* gen_numu = (TH1D*)myFile->Get("generated_numu");
  TH1D* gen_barnue = (TH1D*)myFile->Get("gen_bar_nue");
  TH1D* gen_other = (TH1D*)myFile->Get("gen_other");

  gen_nue->SetFillColor(kBlue);
  gen_numu->SetFillColor(kYellow);
  gen_barnue->SetFillColor(kGreen);
  gen_other->SetFillColor(kRed);

  THStack *nustack = new THStack("nustack","");
  nustack->Add(gen_nue);
  nustack->Add(gen_numu);
  nustack->Add(gen_barnue);
  nustack->Add(gen_other);

  nustack->SetTitle("Generated neutrino types; Neutrino energy (GeV);# Events");

  nustack->Draw();
  gPad->BuildLegend();
}
