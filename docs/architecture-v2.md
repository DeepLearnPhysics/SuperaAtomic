# Supera v2 architecture: first migration slice

This is a v2 development checkpoint, not a released or stable v2 API.
Repository renaming and release publication are separate decisions. The optical
checkpoint is validated; broader event and labeling migrations remain staged.

## Boundaries

```text
native event products
    -> input adapter
    -> supera::EventRecord and shared value types
    -> shared processing where required
    -> supera::EventWriter
    -> selected output implementation
```

The first slice is optical only: extraction and serialization, with no labeling
step required. The existing particle `EventInput`, `EventOutput`, and labeling
algorithms are retained unchanged. Their integration into a broader envelope
will be a later reviewed migration, not an implicit conversion between unrelated
input and result types.

| Target | Responsibility | External dependencies |
| --- | --- | --- |
| `supera::data` | Owned records, conventions, event identity, sink interface | C++17 standard library |
| `supera::labeling` | Existing Atomic algorithms | Data model, yaml-cpp |
| `supera::larcv2` | Existing optical LArCV2 compatibility projection | Data model, LArCV2, ROOT |
| `larsoft2supera::optical` (other repository) | Native OpFlash to shared records | Data model, lardataobj |
| art integration module (other repository) | Event access and chosen sink lifecycle | art, extraction target, selected sink |

`EventWriter` borrows an event synchronously; it cannot retain event references
past `write()`. The current art analyzer owns a serial writer and the per-event
record. art owns scheduling. A standalone application can own a different loop
and use the same records and writer. No second event loop or dynamic plugin ABI
has been introduced.

The root CMake entry point exports independently consumable targets. The default
build is data-only; opt into `SUPERA_BUILD_LABELING` and/or `SUPERA_BUILD_LARCV2`.
The old `src/` entry point and Python packaging are unchanged during migration.
The new entry point does not yet provide Python/ROOT bindings for the v2 records.

## Optical record contract

`OpticalFlash` preserves every stored OpFlash attribute required by legacy
SuperaOptical. IDs are collection-local; vectors retain their ordering and own
their values. There is no recalibration, filtering, time shift, or sentinel
replacement. Stored beam flags and fast fractions are not asserted to be measured.

`OpticalFlashCollection` separates the logical collection name from the configured
native source identifier. `OpticalConventions` records detector, region, coordinate
frame, PE indexing basis, and absolute-time unit/reference. Unspecified values are
explicitly unknown. This context describes data; shared algorithms must not parse
source tags or switch on experiment names to recover missing semantics.

Trigger-relative times and widths are in microseconds, spatial centers and widths
in centimeters. Absolute time remains a source value with a separate convention.
Wire vectors retain the source view convention. The initial LArSoft adapter does
not infer PE indexing or absolute-time semantics from field names. Detector,
region, and coordinate frame can be supplied in collection configuration.

This is a typed extension point, not an arbitrary string-to-value property bag.
Future detector differences should motivate explicit additions, with units,
availability, provenance, and tests. In particular, per-field measured/defaulted
status and explicit non-contiguous sensor-ID mappings are not implemented yet.
Do not infer those capabilities from the presence of a conventions record.

## LArCV2 boundary

Only `src/supera/io/larcv2/Writer.cxx` includes LArCV headers in the new production
code. Its public header hides the concrete implementation. Input conversion and
shared records do not use LArCV classes or constants.

The writer intentionally projects the existing fields, collection names, and event
identity into the old format. That format does not persist the new source and
conventions records. This is an explicit compatibility limitation, not an assertion
that metadata survives serialization. Associations, CRT matching, and optical
truth remain deferred; no new physics output is introduced here.

The writer validates collection inventory, order, names, and representable IDs
before changing output state. Each event replaces collection contents, including
empty vectors. `close()` reports finalization errors; destructor cleanup is best
effort and cannot substitute for an explicit successful close.

## Build and validation

```sh
cmake -S . -B build/data -DCMAKE_INSTALL_PREFIX="$PWD/build/install"
cmake --build build/data
ctest --test-dir build/data --output-on-failure
cmake --install build/data
```

Consumers use `find_package(Supera 2 CONFIG REQUIRED COMPONENTS data)` and link
`supera::data`. Request `labeling` or `larcv2` only when needed. Loading the data
package does not search for either optional component's dependencies.

The data test checks ownership and explicit unknown conventions. The optional
writer test checks validation, lifecycle, independent field expectations, event
identity, and clearing between events. A labeling smoke test verifies linkage,
voxel boundaries, and deposit aggregation without changing algorithm behavior.
It is not full labeling parity validation.

The coordinated larsoft2supera integration suite builds and installs both
repositories, builds a separate consumer, verifies fixture checksums, exercises
art product lookup and empty collections, and compares both flash products
against both frozen legacy references. It also rejects deliberately corrupted
comparison fixtures. Its large input and UPS products remain externally provisioned.
