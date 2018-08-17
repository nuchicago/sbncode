#include "TFile.h"
#include "TTree.h"
#include <vector>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <algorithm>

void AnalyzeNumuCut (){
  TFile *myFile = TFile::Open("output_SBNOsc_NueSelection.root");
  if (myFile==0){
    printf("File not correctly opened!\n");
    return;
  } //if the file is not correctly opened, spit out an error message
  TH1D *gen_nuhist = (TH1D*)myFile->Get("generated_nue_hist");
  TH1D *gen_nuhist_fidvol = (TH1D*) myFile->Get("generated_nue_in_fiducial_volume");
  TH1D *shower_fid = (TH1D*)myFile->Get("shower_cut_fid_only_nu_energy");
  TH1D *shower_fid_track = (TH1D*)myFile->Get("selected_nu_hist");
  TH1D *shower_fid_track_cg = (TH1D*)myFile->Get("final_selected_nu");
  TH1D *shower_fid_track_cg_reco = (TH1D*)myFile->Get("final_selected_nu_w_reco_efficiency");


  gen_nuhist->SetFillColor(kAzure);
  gen_nuhist_fidvol->SetFillColor(kViolet);
  shower_fid->SetFillColor(kGreen);
  shower_fid_track->SetFillColor(kBlue);
  shower_fid_track_cg->SetFillColor(kYellow);
  shower_fid_track_cg_reco->SetFillColor(kRed);
  //TCanvas *c = new TCanvas ("c", "Generated and recontructed hists",10,10,1000,800);

  THStack *nustack = new THStack("nustack","Generated and reconstructed #nu_e after cuts");
  nustack->Add(nuhist);
  nustack->Add(goodnuhist);


  //nustack->Draw("nostack");

  nustack->SetTitle("Generated and recontructed #nu_e after cuts; Neutrino energy (GeV);# Events ");
  nustack->Draw("nostack");
  auto legend = new TLegend();
  legend->AddEntry(gen_nuhist, "Generated #nu_e");
  legend->AddEntry(gen_nuhist_fidvol, "Generated #nu_e after fiducial volume selection");
  legend->AddEntry(shower_fid,"+ >200MeV shower energy selection");
  legend->AddEntry(shower_fid_track, "+ < 1m track length selction");
  legend->AddEntry(shower_fid_track_cg, "+ conversion gap cut");
  legend->AddEntry(shower_fid_track_cg_reco, "+ 0.8 reco efficiency");
  legend->Draw();

  //nustack->GetXaxis()->SetTitle("Neutrino energy (GeV)");
  //nustack->GetYaxis()->SetTitle("# Events");

  //c-> Modified();






}
