#include "TFile.h"
#include "TTree.h"
#include <vector>
#include <string>
#include <algorithm>

void ShowerTypeStack() {
  TFile *myFile = TFile::Open("output_SBNOsc_NueSelection.root");
  if (myFile==0){
    printf("File not correctly opened!\n");
    return;
  }

  TH1D* e_shower_nu = (TH1D*)myFile->Get("e_shower");
  TH1D* gamma_shower_nu = (TH1D*)myFile->Get("gamma_shower");
  TH1D* other_shower_nu = (TH1D*)myFile->Get("other_shower");

  e_shower_nu->SetFillColor(kBlue);
  gamma_shower_nu->SetFillColor(kYellow);
  other_shower_nu->SetFillColor(kViolet);

  THStack *nustack = new THStack("nustack","");
  nustack->Add(e_shower_nu);
  nustack->Add(gamma_shower_nu);
  nustack->Add(other_shower_nu);

  nustack->SetTitle("True shower type of selected #nu_e candidates; Neutrino energy (GeV);# Events");

  nustack->Draw();
  gPad->BuildLegend();


}
