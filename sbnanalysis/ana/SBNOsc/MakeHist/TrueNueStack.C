#include "TFile.h"
#include "TTree.h"
#include <vector>
#include <string>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <algorithm>

void TrueNueStack() {
  TFile *myFile = TFile::Open("output_SBNOsc_NueSelection_0821.root");
  if (myFile==0){
    printf("File not correctly opened!\n");
    return;
  }

  TH1D* selectednue = (TH1D*)myFile->Get("true_nue_selected");
  TH1D* shower_fid_track_cg = (TH1D*)myFile->Get("final_selected_nu");

  double ntruenue = selectednue->GetEntries();
  double nselected = shower_fid_track_cg->GetEntries();

  double true_rate = ntruenue/nselected;
  std::cout<<"# true nue selected/# candidates selected = "<<true_rate<<std::endl;

  selectednue->SetFillColor(kAzure);
  shower_fid_track_cg->SetFillColor(kYellow);

  THStack *nustack = new THStack("nustack","");

  nustack->Add(selectednue);
  nustack->Add(shower_fid_track_cg);

  nustack->SetTitle("True nue in selected nue candidates; Neutrino energy (GeV);# Events ");

  nustack->Draw();
  auto legend = new TLegend();
  legend->AddEntry(selectednue, "Selected #nu_e");
  legend->AddEntry(shower_fid_track_cg, "true #nu_e in the selected candidates");
  legend->Draw();

}
