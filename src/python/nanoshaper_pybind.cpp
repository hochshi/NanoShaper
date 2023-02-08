#include <main_functions.h>
#include <pybind11/pybind11.h>
#include <spdlog/spdlog.h>

namespace py = pybind11;

PYBIND11_MODULE(NanoShaper, m) {
  m.doc() = "";

  m.def(
      "load_config",
      &load,
      "init streams, check configuration file for errors and read variables",
      py::arg("confFile")
    );
}
