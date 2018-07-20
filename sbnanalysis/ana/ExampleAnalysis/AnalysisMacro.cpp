/**
 * An example (compiled) analysis script.
 *
 * To run it, compile with:
 *
 *   g++ -o AnalysisMacro \
 *     AnalysisMacro.cpp $SBN_LIB_DIR/libsbnanalysis_Event.so \
 *     -I../../core -I$CANVAS_INC -I$CETLIB_EXCEPT_INC -I$CETLIB_INC \
 *     `root-config --libs --cflags`
 *
 * then run:
 *
 *   ./AnalysisMacro input_file.root
 *
 * A. Mastbaum <mastbaum@uchicago.edu>, 2018/07/19
 */

#include <cassert>
#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include <Event.hh>

int main(int argc, char* argv[]) {
  // Check for command line arguments
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " input_file.root" << std::endl;
    return 1;
  }

  // Load the input file
  TFile* f = TFile::Open(argv[1]);
  assert(f && f->IsOpen());

  TTree* t = (TTree*) f->Get("sbnana");
  assert(t && t->GetEntries());

  // Set up branches
  Event* event = new Event;
  std::vector<int>* parentPDG = new std::vector<int>;
  t->SetBranchAddress("events", &event);
  t->SetBranchAddress("parentPDG", &parentPDG);

  for (long i=0; i<t->GetEntries(); i++) {
    t->GetEntry(i);
    assert(event->interactions.size() == parentPDG->size());

    std::cout << "Event " << i << ": ";

    for (size_t j=0; j<event->interactions.size(); j++) {
      double enu = event->interactions[j].neutrino.energy;
      double ptype = parentPDG->at(j);

      std::cout << "Neutrino " << j << ": "
                << "E=" << enu << " GeV, "
                << "Parent PDG=" << ptype;

      if (j < event->interactions.size() - 1) {
        std::cout << "; ";
      }
    }

    std::cout << std::endl;
  }
}

