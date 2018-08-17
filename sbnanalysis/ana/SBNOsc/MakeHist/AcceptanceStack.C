#include "TFile.h"
#include "TTree.h"
#include <vector>
#include <string>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <algorithm>

void AcceptanceStack (){
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

  double ngen =  gen_nuhist->GetEntries();
  double ngenfid = gen_nuhist_fidvol->GetEntries();
  double nshowerfid = shower_fid->GetEntries();
  double nshowerfidtrack = shower_fid_track->GetEntries();
  double nshowerfidtrackcg = shower_fid_track_cg->GetEntries();
  double nshowerfidtrackcgreco = shower_fid_track_cg_reco->GetEntries();

  double ngenfid_rate = ngenfid/ngen;
  double nshowerfid_rate = nshowerfid/ngenfid;
  double nshowerfidtrack_rate = nshowerfidtrack/ngenfid;
  double nshowerfidtrackcg_rate = nshowerfidtrackcg/ngenfid;
  double nshowerfidtrackcgreco_rate = nshowerfidtrackcgreco/ngenfid;

  std::cout<<"Generated nue in fid vol/generated nue = "<<ngenfid_rate<<std::endl;
  std::cout<<"Selected nue in fid vol passing shower energy cut/generated nue in fid vol = "<<nshowerfid_rate<<std::endl;
  std::cout<<"Selected nue in fid vol passing shower energy cut and track length cut/generated nue in fid vol = "<<nshowerfidtrack_rate<<std::endl;
  std::cout<<"Selected nue in fid vol passing shower energy cut, track length cut, and conversion gap cut/generated nue in fid vol = "<<nshowerfidtrackcg_rate<<std::endl;
  std::cout<<"Selected nue in fid vol passing shower energy cut, track length cut, conversion gap cut with reco efficiency applied/generated nue in fid vol = "<<nshowerfidtrackcgreco_rate<<std::endl;



  gen_nuhist->SetFillColor(kAzure);
  gen_nuhist_fidvol->SetFillColor(kViolet);
  shower_fid->SetFillColor(kGreen);
  shower_fid_track->SetFillColor(kBlue);
  shower_fid_track_cg->SetFillColor(kYellow);
  shower_fid_track_cg_reco->SetFillColor(kRed);
  //TCanvas *c = new TCanvas ("c", "Generated and recontructed hists",10,10,1000,800);

  THStack *nustack = new THStack("nustack","Generated and reconstructed #nu_e after cuts");
  nustack->Add(gen_nuhist);
  nustack->Add(gen_nuhist_fidvol);
  nustack->Add(shower_fid);
  nustack->Add(shower_fid_track);
  nustack->Add(shower_fid_track_cg);
  nustack->Add(shower_fid_track_cg_reco);


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
