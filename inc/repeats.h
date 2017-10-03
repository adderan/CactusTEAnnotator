#ifndef _HALREPEATS_H
#define _HALREPEATS_H

#include "hal.h"

namespace hal {

class HalRepeats {
public:

   HalRepeats();
   virtual ~HalRepeats();

   void analyzeBranch(AlignmentConstPtr alignment,
                      hal_size_t gapThreshold,
                      double nThreshold);

protected:

   AlignmentConstPtr _alignment;

};

}
#endif
