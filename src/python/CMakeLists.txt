################################################################################
# PYBIND11_LIB
################################################################################
if(ENABLE_PYBIND11_LIB)
  SET(Python_FIND_STRATEGY LOCATION)
  find_package(Python 3.9 COMPONENTS Interpreter Development.Module REQUIRED)
  find_package(pybind11 CONFIG REQUIRED)
endif()

if(ENABLE_PYBIND11_LIB AND Python_FOUND AND pybind11_FOUND)
  add_definitions(-DPYTHON_ENABLED)
  message(STATUS "Using python ${Python_EXECUTABLE} ${Python_INCLUDE_DIRS} ${Python_LIBRARY_DIRS} ${Python_LIBRARIES} ${pybind11_INCLUDE_DIRS} ${pybind11_LIBRARIES}")
  
  link_directories(${Python_LIBRARY_DIRS})
  link_libraries(${Python_LIBRARIES})
  include_directories(${pybind11_INCLUDE_DIRS})
  link_libraries(${pybind11_LIBRARIES})
  
  add_library(NanoShaper MODULE
    ./nanoshaper_pybind.cpp
    ${SOURCE_FILES}
    ../main_functions.h
    ../sturm/sturm.cpp
    ../sturm/util.cpp
    ../sturm/solve.h
    ../ply/ply.cpp
    ../ply/ply.h
    ../ConfigFile/ConfigFile.cpp
    ../ConfigFile/ConfigFile.h 
    ../octree/octree.h 
    )
    
  target_link_options(NanoShaper PUBLIC -undefined dynamic_lookup)
  target_link_libraries(NanoShaper
    PRIVATE pybind11::headers
  )
  pybind11_extension(NanoShaper)

  set_target_properties(NanoShaper
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
  )
  
endif()

################################################################################
# SWIG_LIB
################################################################################
if(ENABLE_SWIG_LIB)
  find_package(SWIG REQUIRED)
  find_package(PythonLibs REQUIRED)
endif()

if(ENABLE_SWIG_LIB AND SWIG_FOUND AND PYTHONLIBS_FOUND)
  # this is added due to a bug of Cmake with Swig on linux
  execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory "src")


  SET_SOURCE_FILES_PROPERTIES(./nanoshaper.i PROPERTIES CPLUSPLUS ON)

  SWIG_ADD_MODULE(NanoShaper python 
    ./nanoshaper.i 
    ${SOURCE_FILES}
    ../ply/ply.cpp 
    ../ply/ply.h
    ../ConfigFile/ConfigFile.cpp
    ../ConfigFile/ConfigFile.h
    ../sturm/sturm.cpp
    ../sturm/util.cpp
    ../sturm/solve.h
    ../octree/octree.h
  )
  SWIG_LINK_LIBRARIES(NanoShaper ${PYTHON_LIBRARIES})
endif()
