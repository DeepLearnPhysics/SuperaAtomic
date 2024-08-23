#ifdef WITH_PYBIND
#include "process_pybind.h"

// needed so the std::map<> in the Driver::ConfigureBBoxAlgorithm() (etc.) signature gets handled right
#include "pybind11/stl.h"

#include "supera/pybind_mkdoc.h"

#include "Driver.h"

void init_process(pybind11::module& m)
{
  pybind11::class_<supera::Driver>(m, "Driver", DOC(supera, Driver))
      .def(pybind11::init(), DOC(supera, Driver, Driver))
      .def("ConfigureBBoxAlgorithm", &supera::Driver::ConfigureBBoxAlgorithm, DOC(supera, Driver, ConfigureBBoxAlgorithm))
      .def("ConfigureLabelAlgorithm", &supera::Driver::ConfigureLabelAlgorithm, DOC(supera, Driver, ConfigureLabelAlgorithm))
      .def("Reset", &supera::Driver::Reset, DOC(supera, Driver, Reset))
      .def("GenerateImageMeta", &supera::Driver::GenerateImageMeta, DOC(supera, Driver, GenerateImageMeta))
      .def("GenerateLabel", &supera::Driver::GenerateLabel, DOC(supera, Driver, GenerateLabel))
      .def("Label", &supera::Driver::Label, DOC(supera, Driver, Label))
      .def("Meta", &supera::Driver::Meta, DOC(supera, Driver, Meta));

}
#endif
