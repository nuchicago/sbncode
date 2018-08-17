#include "TFile.h"
#include "TTree.h"
#include <vector>
#include <string>
#include <algorithm>

void InteractionType() {
  TFile *myFile = TFile::Open("output_SBNOsc_NueSelection.root");
  if (myFile==0){
    printf("File not correctly opened!\n");
    return;
  }

}
