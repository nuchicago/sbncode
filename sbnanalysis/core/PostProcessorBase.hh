#ifndef __sbnanalysis_core_PostProcessorBase__
#define __sbnanalysis_core_PostProcessorBase__

/**
 * \file PostProcessorBase.hh
 *
 * A generic processor that writes an sbnanalysis tree.
 *
 * Author: G. Putnam <gputnam@uchicago.edu>, 2018/10/08
 */

#include <string>
#include <vector>
#include "gallery/Event.h"
#include "Loader.hh"
#include "Event.hh"

class TBranch;
class TFile;
class TTree;
class Event;

namespace Json {
  class Value;
}

/** Core framework functionality. */
namespace core {

/**
 * \class core::PostProcessorBase
 * \brief A generic tree-writing event-by-event processor.
 */
class PostProcessorBase {
public:
  /** Constructor */
  PostProcessorBase();

  /** Destructor */
  virtual ~PostProcessorBase();

  /** run */
  void Run(std::vector<std::string> filelist);

protected:

  /**
   * Process one event.
   */
  virtual void ProcessEvent(const Event *event) = 0;

  /** setup anything needed per file */
  virtual void FileSetup(TTree *eventTree) {} 
  /** any cleanup needed per file */
  virtual void FileCleanup(TTree *eventTree) {}
  /**
   * Perform user-level initialization.
   *
   * \param config A configuration, as a JSON filename.
   */
  virtual void Initialize(char* config=NULL);

  /**
   * Perform user-level initialization.
   *
   * \param config A configuration, as a JSON object.
   */
  virtual void Initialize(Json::Value* config=NULL) = 0;

  /** Perform user-level finalization. */
  virtual void Finalize() {}

  Event *fEvent;
  TTree *fEventTree;
};

}  // namespace core


/** Macro to create plugin library for user-defined PostProcessors. */
#define DECLARE_SBN_POSTPROCESSOR(classname) extern "C" { \
core::PostProcessorBase* CreatePostProcessorObject() { return new classname; } \
void DestroyPostProcessorObject(core::PostProcessorBase* o) { delete o; } \
struct core::export_table_postprocess exports = { CreatePostProcessorObject, DestroyPostProcessorObject };}

#endif  // __sbnanalysis_core_PostProcessorBase__

