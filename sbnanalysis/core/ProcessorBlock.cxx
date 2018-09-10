#include <string>
#include <vector>
#include <json/json.h>
#include "ProcessorBase.hh"
#include "ProcessorBlock.hh"

namespace core {

ProcessorBlock::ProcessorBlock() {}


ProcessorBlock::~ProcessorBlock() {}


void ProcessorBlock::AddProcessor(ProcessorBase* processor,
                                  Json::Value* config) {
  fProcessors.push_back({processor, config});
}


void ProcessorBlock::ProcessFiles(std::vector<std::string> filenames) {
  // Setup
  for (auto it : fProcessors) {
    it.first->Setup(it.second);
    it.first->Initialize(it.second);
  }

  // Event loop
  for (gallery::Event ev(filenames); !ev.atEnd(); ev.next()) {
    for (unsigned processor_ind = 0; processor_ind < fProcessors.size(); processor_ind++) {
      auto &it = fProcessors[processor_ind];
      it.first->BuildEventTree(ev);

      // catch an art exceptions
      bool accept = false;
      try {
        accept = it.first->ProcessEvent(ev, it.first->fEvent->truth, *it.first->fReco);
      }
      catch (cet::exception e) {
        // in the event of an exception, gracefully fill the processor
        std::cerr << "Caught art exception. Closing and removing processor." << std::endl;
        std::cerr << e.explain_self() << std::endl;

        it.first->EventCleanup();
        it.first->Finalize();
        it.first->Teardown();
        fProcessors.erase(fProcessors.begin() + processor_ind);
        // need to bump back the index to account for erased processor
        //
        // unsigned underflow/overflow arithmetic is well-defined, so this is ok even if processor_ind == 0
        processor_ind --;
        // to the next processor
        continue;
      }

      if (accept) {
        it.first->FillTree();
      }

      it.first->EventCleanup();
    }
    // finish if no processors are left
    if (fProcessors.size() == 0) {
      std::cerr << "No remaining processors. Exiting." << std::endl;
      break;
    }
  }

  // Finalize
  for (auto it : fProcessors) {
    it.first->Finalize();
    it.first->Teardown();
  }
}


void ProcessorBlock::DeleteProcessors() {
  for (auto it : fProcessors) {
    delete it.first;
  }
}

}  // namespace core

