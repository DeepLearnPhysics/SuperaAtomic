#ifndef SUPERA_IO_EVENTWRITER_H
#define SUPERA_IO_EVENTWRITER_H

#include "supera/data/EventRecord.h"

namespace supera {
/// A synchronous sink. write() must finish using the borrowed event before return.
/// Call close() explicitly to report finalization errors. Destructors do not throw.
class EventWriter {
public:
  virtual ~EventWriter() = default;
  virtual void write(const EventRecord &event) = 0;
  virtual void close() = 0;
};
} // namespace supera
#endif
