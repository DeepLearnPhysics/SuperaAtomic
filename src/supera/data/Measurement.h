#ifndef SUPERA_DATA_MEASUREMENT_H
#define SUPERA_DATA_MEASUREMENT_H

namespace supera {
enum class TimeUnit { Unspecified, Nanosecond, Microsecond, Tick, Second };
enum class TimeReference { Unspecified, SourceDefined, ClockEpoch, Trigger };
enum class SignalUnit { Unspecified, Photoelectron, ADC, MeV };
} // namespace supera
#endif
