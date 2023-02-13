include_directories (BEFORE include)
include_directories(.)
include_directories(./ConfigFile)
include_directories(./jama)
include_directories(./tnt)
include_directories(./sturm)
include_directories(./octree)

file(GLOB SOURCE_FILES # RELATIVE ${CMAKE_SOURCE_DIR}
  "*.h"
  "*.cpp")

################################################################################
# EXE
################################################################################
option(ENABLE_EXE "EXE" OFF)

if(ENABLE_EXE)
  add_subdirectory(exe)
endif()

################################################################################
# DELPHI_LIB
################################################################################
option(ENABLE_LIB_DELPHI "LIB_DELPHI" OFF)

if(ENABLE_LIB_DELPHI)
  add_subdirectory(libraries)
endif()

################################################################################
# SWIG_LIB
################################################################################
option(ENABLE_SWIG_LIB "SWIG_LIB" OFF)

if(ENABLE_SWIG_LIB)
  add_subdirectory(python)
endif()

################################################################################
# PYBIND11_LIB
################################################################################
option(ENABLE_PYBIND11_LIB "PYBIND11_LIB" OFF)

if(ENABLE_PYBIND11_LIB)
  add_subdirectory(python)
endif()
