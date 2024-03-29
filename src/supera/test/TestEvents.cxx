#include "supera/test/TestEvents.h"

#include "supera/data/Particle.h"

namespace supera
{
  namespace test
  {

    bool VerifyEventMeta(const ImageMeta3D &computedMeta, const ImageMeta3D &expectedMeta)
    {
      return computedMeta == expectedMeta;
    }

    bool VerifyEventLabels(const EventOutput &computedLabels, const EventOutput &expectedLabels)
    {
      return computedLabels == expectedLabels;
    }

    std::map<std::string, TestEvent> TestEvents()
    {
      return {
//          {"NumuCCIncEvt", NumuCCIncEvt()}
          {"PCMuonEvt", PCMuonEvt()}
      };
    };

  } // namespace test
} // namespace supera
