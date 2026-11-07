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

The current slices cover optical flashes, auxiliary space points, and generator truth: extraction
and serialization, with no labeling step required. The existing particle `EventInput`, `EventOutput`, and labeling
algorithms are retained unchanged. Their integration into a broader envelope
will be a later reviewed migration, not an implicit conversion between unrelated
input and result types.

| Target | Responsibility | External dependencies |
| --- | --- | --- |
| `supera::data` | Owned records, conventions, event identity, sink interface | C++17 standard library |
| `supera::labeling` | Existing Atomic algorithms | Data model, yaml-cpp |
| `supera::larcv2` | Optical and supported auxiliary LArCV2 projections | Data model, LArCV2, ROOT |
| `larsoft2supera::optical` (other repository) | Native OpFlash to shared records | Data model, lardataobj |
| `larsoft2supera::auxiliary` (other repository) | Native SBN CRT to shared auxiliary records | Data model, sbnobj |
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

## Auxiliary record contract

`AuxDetSpacePoint` represents a localized reconstructed measurement, not detector
hardware or every possible auxiliary instrument. Technology and subsystem roles
are collection metadata: scintillator technology does not imply cosmic tagging,
and a cosmic-tagger measurement does not classify its particle as cosmic.

Position and uncertainty (cm), signal with explicit units, named clock readings,
source-scoped readout/channel IDs, plane/region annotations, and reconstruction
completeness are separate typed fields. Missing estimates stay absent. The initial
reader preserves numeric source sentinels without declaring them physically valid.
Unsigned integer clocks remain exact above 2^53; a clock reading also declares its
unit and reference. An uncertainty, when present, has that reading's unit.

The `SbnCrtHit` source contract preserves all five raw timestamp fields without
assuming common epoch or correction semantics. It is not an output format or a
cosmic-particle flag. Readout groups preserve empty groups and channel ordering;
readout IDs are not AuxDet geometry IDs. Geometry links, raw triggers, tracks,
TOF relations, and other instrument records remain future typed extensions.

The current LArCV2 projection supports this source contract only. It rejects
missing required values, unsupported clock conventions, ADC/energy in place of
PE, lossy float conversions, out-of-range IDs, duplicate readout groups, and
nonrepresentable timing uncertainties/completeness. A future output can retain
these quantities natively. The existing source-preserving ICARUS mapping does
not invent values to satisfy the old output schema.

## Generator truth contract

`GeneratorCollection` preserves ordered source records, availability, original
MCTruth indices, and original particle indices. Track IDs remain signed. A source
slot is stable even when a preceding product is absent; a present-empty product
is distinct from an absent one. The separate incoming-neutrino copy has no
invented particle-list index.

`GeneratorStep` names cm/ns/GeV units explicitly. Missing initial kinematics or
neutrino theta remains absent until selection. `select_generator` applies
origin and status selection independently of native extraction and serialization;
it then rejects retained records without required kinematics. An omitted origin
filter accepts every origin; explicitly supplied zero selects unknown origin.

The optional LArCV2 projection is intentionally named
`project_legacy_generator`. It preserves the pinned SuperaMCTruth output's
mixed particle units, unsigned track-ID conversion, present-producer counter,
positron shape rule, self-ancestry placeholders and object-ID sentinel behavior.
These are compatibility decisions, not shared physics or identity conventions.
The two products are prepared before output mutation and cleared on empty events.

This model does not replace the transported `Particle`/`ParticleInput` labeling
model. Generator metadata, ancestry links, and particle trajectory tails are not
newly supported here; the incoming neutrino's full trajectory is retained because
the existing output contains it. Linkage to transported truth is a later stage.

## LArCV2 boundary

Only the private implementation under `src/supera/io/larcv2` includes LArCV
headers in the new production code. Its public header hides the concrete implementation. Input conversion and
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
art product lookup and empty collections, and compares optical-only and combined optical/CRT products
against both frozen legacy references. It also rejects deliberately corrupted
comparison fixtures. Its large input and UPS products remain externally provisioned.
