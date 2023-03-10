#include <BlobbySurface.h>
#include <ConfigFile.h>
#include <Configuration.h>
#include <ConnollySurface.h>
#include <DelphiShared.h>
#include <ExampleSurface.h>
#include <MeshSurface.h>
#include <SkinSurface.h>
#include <Surface.h>
#include <SurfaceFactory.h>
#include <globals.h>
#include <logging.h>
#include <main_functions.h>
#include <pybind11/detail/common.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tools.h>
#include <memory>

namespace py = pybind11;

using namespace nanoshaper;

namespace pybind_wrappers {
struct InputDataPy : InputData {
  py::list get_ai();
  void set_ai(const py::list& ai_py);
  py::list aiPy;
};

py::list InputDataPy::get_ai() {
  py::list ai_py;
  for (const auto& x : this->ai) {
    ai_py.append(x);
  }

  return ai_py;
}

void InputDataPy::set_ai(const py::list& ai_py) {
  this->ai.clear();
  for (int i = 0; i < py::len(ai_py); ++i) {
    this->ai.emplace_back(ai_py[i].cast<AtomInfo>());
  }
}

}  // namespace pybind_wrappers

PYBIND11_MODULE(NanoShaper, m) {
  m.doc() = "";
  py::bind_map<std::map<std::string, std::string>>(m, "MapStringString");

  m.def("load_config", &load,
        "init streams, check configuration file for errors and read variables",
        py::arg("filename"), py::arg("delimiter") = "=",
        py::arg("comment") = "#", py::arg("sentry") = "EndConfigFile",
        py::arg("format") = "plain");
  m.def("parse_config", &parse, py::arg("cf"));
  m.def("normalMode", &normalMode, py::arg("surf"), py::arg("dg"),
        py::arg("conf"));
  m.def("pocketMode", &pocketMode, py::arg("hasAtomInfo"), py::arg("cf"),
        py::arg("conf"));
  m.def(
      "createSurface",
      [](ConfigurationOP conf, DelPhiSharedOP ds) {
        SurfaceOP base = SurfaceFactory::getInstance().create(conf, ds);
        return base;
      },
      py::arg("conf"), py::arg("ds"));
  m.def("printAvailableSurfaces",
        []() { SurfaceFactory::getInstance().print(); });

  py::class_<Configuration, std::shared_ptr<Configuration>>(m, "Configuration")
      .def(py::init())
      .def_readwrite("cavVol", &Configuration::cavVol)
      .def_readwrite("numMol", &Configuration::numMol)

      // grid (DelPhi) params
      .def_readwrite("scale", &Configuration::scale)
      .def_readwrite("perfill", &Configuration::perfill)
      .def_readwrite("isAvailableAtomInfo", &Configuration::isAvailableAtomInfo)
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

      .def_readwrite("blobby_B", &Configuration::blobby_B)
      .def_readwrite("maxSESDim2D", &Configuration::maxSESDim2D)
      .def_readwrite("maxSESPatches2D", &Configuration::maxSESPatches2D)
      .def_readwrite("maxSESDim", &Configuration::maxSESDim)
      .def_readwrite("maxSESPatches", &Configuration::maxSESPatches)
      .def_readwrite("mp", &Configuration::mp)
      .def_readwrite("si_perfil", &Configuration::si_perfil)
      .def_readwrite("radius", &Configuration::radius)

      .def_readwrite("sfname", &Configuration::sfname)
      .def_readwrite("maxMeshDim", &Configuration::maxMeshDim)
      .def_readwrite("maxMeshPatches", &Configuration::maxMeshPatches)
      .def_readwrite("maxMeshDim2D", &Configuration::maxMeshDim2D)
      .def_readwrite("maxMeshPatches2D", &Configuration::maxMeshPatches2D)
      .def_readwrite("NumMSMSfiles", &Configuration::NumMSMSfiles)

      .def_readwrite("skin_s", &Configuration::skin_s)
      .def_readwrite("maxSkinDim", &Configuration::maxSkinDim)
      .def_readwrite("maxSkinPatches", &Configuration::maxSkinPatches)
      .def_readwrite("maxSkinDim2D", &Configuration::maxSkinDim2D)
      .def_readwrite("maxSkinPatches2D", &Configuration::maxSkinPatches2D)
      .def_readwrite("useFastProjection", &Configuration::useFastProjection)
      .def_readwrite("savePovRay", &Configuration::savePovRay)

      .def_readwrite("checkDuplicatedVertices",
                     &Configuration::checkDuplicatedVertices)
      .def_readwrite("wellShaped", &Configuration::wellShaped)
      .def_readwrite("probeRadius", &Configuration::probeRadius)
      .def_readwrite("lb", &Configuration::lb)
      .def_readwrite("vaFlag", &Configuration::vaFlag)
      .def_readwrite("computeNormals", &Configuration::computeNormals)
      .def_readwrite("saveMSMS", &Configuration::saveMSMS)
      .def_readwrite("sternLayer", &Configuration::sternLayer)
      .def_readwrite("Max_Atoms_Multi_Grid",
                     &Configuration::Max_Atoms_Multi_Grid)
      .def_readwrite("surfName", &Configuration::surfName)

      .def("stopDebug", &Configuration::stopDebug)
      .def("restartDebug", &Configuration::restartDebug);

  py::class_<ConfigFile, std::shared_ptr<ConfigFile>>(m, "ConfigFile")
      .def(py::init<const std::string&, const std::string&, const std::string&,
                    const std::string&, std::string&>(),
           py::arg("filename"), py::arg("delimiter") = "=",
           py::arg("comment") = "#", py::arg("sentry") = "EndConfigFile",
           py::arg("format") = "plain")
      .def_property("contents", &ConfigFile::getContents,
                    &ConfigFile::setContents)
      .def("getString", &ConfigFile::readString, py::arg("key"))
      .def("getDouble", &ConfigFile::readFloat, py::arg("key"));

  py::class_<DelPhiShared, std::shared_ptr<DelPhiShared>>(m, "DelphiShared")
      .def(py::init<const double&, const double&, const std::string&,
                    const bool&, const bool&, const bool&, const bool&>(),
           py::arg("scale"), py::arg("perfill"), py::arg("fn"),
           py::arg("eps_flag"), py::arg("stat_flag"), py::arg("multi"),
           py::arg("atinfo") = false)
      .def(py::init<const bool&, const bool&, const bool&, const bool&>(),
           py::arg("eps_flag"), py::arg("stat_flag"), py::arg("multi"),
           py::arg("atinfo") = false)
      .def(py::init<const Configuration&>(), py::arg("conf"))
      .def("init",
           py::overload_cast<double, double, std::string>(&DelPhiShared::init),
           py::arg("scale"), py::arg("perfill"), py::arg("fn"))
      .def("init",
           py::overload_cast<double, double, const InputData&>(
               &DelPhiShared::init),
           py::arg("scale"), py::arg("perfill"), py::arg("in"));

  py::class_<Surface, std::shared_ptr<Surface>>(m, "Surface")
      .def_property_readonly("triList", &Surface::gettriList)
      .def_property_readonly("vertList", &Surface::getvertList)
      .def_property_readonly("normalsList", &Surface::getnormalsList)
      .def_property_readonly("vertexAtomsMap", &Surface::getvertexAromsMap)
      .def("getNumTriangles", &Surface::getNumTriangles)
      .def("getNumVertices", &Surface::getNumVertices);

  py::class_<AtomInfo, std::shared_ptr<AtomInfo>>(m, "AtomInfo")
      .def(py::init<>())
      .def(py::init<const std::string&, const int&, const std::string&,
                    const std::string&>(),
           py::arg("name"), py::arg("resNum"), py::arg("resName"),
           py::arg("chain"))
      .def_property("name", &AtomInfo::getName, &AtomInfo::setName)
      .def_property("resNum", &AtomInfo::getResName, &AtomInfo::setResNum)
      .def_property("resName", &AtomInfo::getResName, &AtomInfo::setResName)
      .def_property("chain", &AtomInfo::getChain, &AtomInfo::setChain);

  py::class_<pybind_wrappers::InputDataPy,
             std::shared_ptr<pybind_wrappers::InputDataPy>>(m, "InputData")
      .def(py::init<>())
      .def_readwrite("na", &pybind_wrappers::InputDataPy::na)
      .def_readwrite("x", &pybind_wrappers::InputDataPy::x)
      .def_readwrite("y", &pybind_wrappers::InputDataPy::y)
      .def_readwrite("z", &pybind_wrappers::InputDataPy::z)
      .def_readwrite("r", &pybind_wrappers::InputDataPy::r)
      .def_readwrite("q", &pybind_wrappers::InputDataPy::q)
      .def_readwrite("d", &pybind_wrappers::InputDataPy::d)
      .def_property("ai", &pybind_wrappers::InputDataPy::get_ai,
                    &pybind_wrappers::InputDataPy::set_ai);
}
