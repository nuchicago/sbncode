#include "TFile.h"
#include "TTree.h"
#include <vector>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <algorithm>

void RatioPlotNue(){
  TFile *myFile = TFile::Open("output_SBNOsc_NueSelection_0821.root");
  if (myFile==0){
    printf("File not correctly opened!\n");
    return;
  } //if the file is not correctly opened, spit out an error message
  TH1D *gen_nuhist_fidvol = (TH1D*)myFile->Get("generated_nue_in_fiducial_volume");
  TH1D *shower_fid_track_cg = (TH1D*)myFile->Get("final_selected_nu");
  TH1D *h = (TH1D*)shower_fid_track_cg->Clone("h");
  h->SetLineColor(kBlue);
  h->SetMinimum(0.);
  h->SetMaximum(1.);
  //h->Sumw2();
  h->SetStats(0);
  h->Divide(gen_nuhist_fidvol);
  h->SetMarkerStyle(21);
  h->Draw();
 /*
  int nbins = nuhist->GetNbinsX();
  Double_t nuhistcontent [nbins];
  Double_t goodnuhistcontent [nbins];
  Double_t ratios [nbins];
  Double_t xcoord [nbins];
  for (int i=0; i<nbins; i++) {
    xcoord[i] = nuhist->GetXaxis()->GetBinCenter(i+1);
    nuhistcontent[i]=nuhist->GetBinContent(i+1);
    goodnuhistcontent[i]=goodnuhist->GetBinContent(i+1);
    if (nuhistcontent[i]!=0.) ratios[i]=goodnuhistcontent[i]/nuhistcontent[i];
    else ratios[i]=1.;
  }

  TGraph *ratioplot = new TGraph (nbins, xcoord, ratios);
  ratioplot->SetTitle("Ratio of numu in the events at SBND (truth level); Neutrino energy (GeV); # Events");
  ratioplot->Draw();
  */







}
