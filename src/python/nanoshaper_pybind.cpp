#include <main_functions.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(NanoShaper, m) {
  m.doc() = "";

  m.def(
      "init_config",
      &init,
      "init streams, check configuration file for errors and read variables",
      py::arg("confFile")
    );
}
