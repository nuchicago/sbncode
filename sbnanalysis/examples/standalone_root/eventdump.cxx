#include <cassert>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <Event.hh>

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cout << "Usage: " << argv[0] << " "
              << "sbnanalysis_output_file.root" << std::endl;
    return 0;
  }

  TFile f(argv[1]);
  assert(f.IsOpen());

  std::cout << "FILE " << argv[1] << std::endl;

  TTree* t = dynamic_cast<TTree*>(f.Get("sbnana"));
  assert(t);

  std::cout << "ENTRIES " << t->GetEntries() << std::endl;

  Event* event = new Event;
  t->SetBranchAddress("events", &event);

  for (long i=0; i<t->GetEntries(); i++) {
    t->GetEntry(i);

    std::cout << "== EVENT " << i << " ==" << std::endl;

    std::cout << "Metadata:" << std::endl
              << "  run = " << event->metadata.run << std::endl
              << "  subrun = " << event->metadata.subrun << std::endl
              << "  eventID = " << event->metadata.eventID << std::endl;

    for (size_t j=0; j<event->truth.size(); j++) {
      std::cout << "Truth neutrino " << j << ": " << std::endl
                << "  iscc = " << event->truth[j].neutrino.iscc << std::endl
                << "  pdg = " << event->truth[j].neutrino.pdg << std::endl
                << "  targetPDG = " << event->truth[j].neutrino.targetPDG << std::endl
                << "  genie_intcode = " << event->truth[j].neutrino.genie_intcode << std::endl
                << "  energy = " << event->truth[j].neutrino.energy << std::endl;
    }
  }

  return 0;
}

