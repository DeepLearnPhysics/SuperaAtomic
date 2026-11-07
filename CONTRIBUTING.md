# Contributing

Follow the [readability conventions](docs/readability.md), shared with
larsoft2supera. Use short purpose comments above logical steps and a blank line
between steps; group related statements so the progression is easy to follow.

Apply the convention to new and touched code. Keep legacy formatting sweeps
separate from behavior changes. The current migration review and dependency
boundaries are described in [the v2 architecture](docs/architecture-v2.md).

Use the pinned clang-format 14 configuration. The v2 workflow checks formatting,
builds the shared data and labeling targets, and runs their tests. Changes to the
optional output component also require its writer tests in a configured LArCV2
environment. Coordinated input changes require the larsoft2supera integration
suite and the frozen-reference comparisons appropriate to the affected products.
