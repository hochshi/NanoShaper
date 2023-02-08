#include "globals.h"
#include <main_functions.h>
#include <pybind11/pybind11.h>
#include <spdlog/spdlog.h>

namespace py = pybind11;

PYBIND11_MODULE(NanoShaper, m) {
  m.doc() = "";

  m.def("load_config", &load,
        "init streams, check configuration file for errors and read variables",
        py::arg("confFile"));
  m.def("parse_config", &parse, py::arg("cf"));
  py::class_<Configuration>(m, "Configuration")
      .def(py::init())
      .def_readwrite("cavVol", &Configuration::cavVol)
      .def_readwrite("numMol", &Configuration::numMol)

      // grid (DelPhi) params
      .def_readwrite("scale", &Configuration::scale)
      .def_readwrite("perfill", &Configuration::perfill)
      // mol file name
      .def_readwrite("molFile", &Configuration::molFile)
      // sys name
      .def_readwrite("sysName", &Configuration::sysName)

      .def_readwrite("multi_diel", &Configuration::multi_diel)

      // actions
      .def_readwrite("fillCavities", &Configuration::fillCavities)
      .def_readwrite("buildEpsmaps", &Configuration::buildEpsmaps)
      .def_readwrite("buildStatus", &Configuration::buildStatus)
      .def_readwrite("tri", &Configuration::tri)
      .def_readwrite("accTri", &Configuration::accTri)
      .def_readwrite("smoothing", &Configuration::smoothing)
      .def_readwrite("tri2balls", &Configuration::tri2balls)
      .def_readwrite("projBGP", &Configuration::projBGP)

      // save data
      .def_readwrite("saveEpsmaps", &Configuration::saveEpsmaps)
      .def_readwrite("saveIdebmap", &Configuration::saveIdebmap)
      .def_readwrite("saveBgps", &Configuration::saveBgps)
      .def_readwrite("saveStatusMap", &Configuration::saveStatusMap)
      .def_readwrite("saveCavities", &Configuration::saveCavities)

      // global parameters
      .def_readwrite("operativeMode", &Configuration::operativeMode)
      .def_readwrite("numthd", &Configuration::numthd)
      .def_readwrite("printAvailSurf", &Configuration::printAvailSurf)
      .def_readwrite("currentSeed", &Configuration::currentSeed)

      // pocket detection
      .def_readwrite("cavAndPockets", &Configuration::cavAndPockets)
      .def_readwrite("linkPockets", &Configuration::linkPockets)
      .def_readwrite("pocketRadiusBig", &Configuration::pocketRadiusBig)
      .def_readwrite("pocketRadiusSmall", &Configuration::pocketRadiusSmall)
      .def_readwrite("pocketRadiusLink", &Configuration::pocketRadiusLink)
      .def_readwrite("debug", &Configuration::debug)
      .def_readwrite("debugStatus", &Configuration::debugStatus)
      .def("stopDebug", &Configuration::stopDebug)
      .def("restartDebug", &Configuration::restartDebug);
}
