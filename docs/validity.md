# Validity contract for v2

Status: implemented for the v2 Optical, AuxDet and Generator records and their
LArSoft readers. The legacy labeling types and old `src/` entry point remain
unmigrated; they must be bridged explicitly before use with these records.

## Ownership

The input adapter interprets source validity exactly once, using the native field,
producer contract and supported schema version. Shared records contain canonical
values; labeling never recognizes LArSoft or LArCV magic numbers. The output
adapter translates canonical values into its destination's field-specific encoding.

A format-wide list of suspicious numbers is forbidden. For example, 65535 can be a
valid collection index, and -1 can be a valid generator track ID. Compare native
sentinels before casting, scaling units, computing derived quantities or narrowing.
Unknown producer semantics must remain an explicit unresolved mapping, not a guess
based on the sample's value distribution.

## Canonical representations

| Kind | Representation | Rule |
| --- | --- | --- |
| Collection index or reference | `int64_t`; missing `-1` | Present values are nonnegative and within the referenced collection. Other negative values are errors. |
| Required record identity | Nonnegative index | Never missing. Scope includes event, collection and, where relevant, source/record. |
| Floating measurement | Finite value; missing quiet `NaN` | No finite magic sentinel and no infinity in shared records. Test missingness with `isnan`, never equality. |
| Signed native identifier | Optional signed integer | Preserve valid negative IDs. Do not reserve -1 without a source contract. |
| Exact integer measurement | Optional integer of sufficient width | Clocks retain exact integer precision; never convert them to double merely to obtain NaN. |
| Boolean | `optional<bool>` when unavailable is possible | False is a measured state, not missing. |
| Native categorical code | Optional native code | Unknown/unclassified codes remain meaningful if the source defines them. |
| PDG code | `0` for unspecified | Preserve signed nonzero codes; do not reject unfamiliar codes against an incomplete particle table. |
| Optional record | `optional<Record>` | Absence of a record differs from a present record with missing fields. |
| Collection | Vector plus explicit availability when needed | Present-empty differs from an unavailable input product. |

Use these rules by semantic kind, not by C++ primitive type alone. A plane number
is not necessarily an index into an event collection. Source track IDs are not
particle array offsets. Hardware identifiers are not channel array offsets.
Checked conversions guard representability; unsigned wraparound is never an
internal validity mechanism. The index range is a v2 contract, not an
implicit conversion of legacy unsigned identifiers.

Floating vectors use component-wise NaNs. Missing x does not erase measured y/z.
An algorithm requiring a complete position must check all components. An empty
trajectory means no stored steps; a present step with unavailable components is a
different case. Required geometry and binning parameters must be finite and valid;
they cannot use the missing-measurement allowance.

## Missing data versus malformed data

Known source sentinel and explicit source absence both yield the canonical missing
representation. Keep their different reasons in conversion diagnostics when needed
for audit; do not add multiple numerical sentinels for those reasons.

Nonfinite source measurements normalize to NaN with a diagnostic identifying NaN
or infinity. This is normalization, not permission to use the value: downstream
operations declare which measurements they require and reject unusable inputs
before mutating an event. Negative time, position or calibrated signal is not
intrinsically invalid. Field-specific physical consistency checks are separate
from sentinel recognition.

Out-of-range collection references, impossible narrowing, inconsistent declared
availability, and malformed associations fail conversion by default. A valid
foreign ID whose target is absent remains a valid source ID with no resolved local
reference; report the unresolved relationship. A deliberate no-match result is
also a missing local reference, but is not corrupt input. Do not synthesize a
parent or silently substitute a self-reference.

Diagnostics identify event, source product, record index, field and reason. They
are separate from the physics values and must not impose a per-hit allocation for
ordinary valid data. Reporting may aggregate repeated reasons with counts and
bounded examples. Unknown mappings block that field's migration until reviewed.

## Output and compatibility

Each writer owns a field mapping and performs checked destination conversions.
For LArCV2, the mapping names the exact field, destination type and missing value;
there is no universal LArCV invalid constant. Legacy constructor defaults are not
an implicit policy.

The existing reference comparison remains a compatibility test, not an authority
on internal semantics. A legacy projection can retain documented historical output
quirks while the shared identity is correct. If normalization loses information
needed for byte parity (for example, several native invalid encodings or a NaN
payload), record an explicit expected difference. Do not hide raw sentinel values
in the physics model to manufacture parity. Exact archival reproduction, if later
required, would need separately scoped provenance and is not this contract.

Validate all output projections before appending an event. Reject an unsupported
mapping rather than silently emitting a destination default. Future output formats
can encode these semantics directly; their schema and storage masks remain a
separate design decision.

## Evidence and rollout

The local SPINE checkout inspected at `e00f194030549bd6cd57cdb29c34d2909c75aeb8`
uses -1 for missing particle indices and NaN for floating fields in
`src/spine/data/larcv/particle.py`. Its `constants/sentinels.py` and
`io/parse/larcv/utils/particle.py` separately handle legacy LArCV encodings. This
contract adopts its useful internal conventions, with an explicit exception for
signed source IDs: SPINE's track-ID default cannot be copied indiscriminately to
generator-native identifiers.

`Validity.h` defines `RecordIndex`, `missing_index`, `missing_value`,
`unspecified_pdg`, checked index conversion and validation helpers. Source IDs and
categorical codes use `std::optional`; integer clock absence uses `std::monostate`
in the clock's tagged variant. This retains the distinction between an unavailable
integer reading and a floating NaN without converting clocks to doubles.

The three existing readers normalize native absence and nonfinite values. Their
optional report records complete reason counts and up to 16 field examples per
source/event. The art module logs these with event and product identity. Reports
are audit information, not part of the physics model or a guarantee of NaN payload
preservation.

The legacy labeling API still has its own `kINVALID_*` constants. Do not use those
in new records. Introduce an explicit checked bridge when porting the first
labeling operation; no unreviewed mass replacement of legacy constants is included.

For each step test source defaults, every documented sentinel, valid values near
sentinels, negative source IDs, partial vectors, nonfinite values, integer range
limits, missing versus empty inputs, and malformed references. Check canonical
records independently of the legacy reference, then check the output projection.
