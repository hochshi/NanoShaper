################################################################################
# DELPHI_LIB
################################################################################
option(ENABLE_LIB_DELPHI "LIB_DELPHI" OFF)

if(ENABLE_LIB_DELPHI)
  set(CMAKE_CXX_STANDARD 14)
  set(CXX_STANDARD_REQUIRED ON)
  
  add_definitions(-DDELPHI_BIND)
  
  add_library(DelphiSurface SHARED 
    ./delphi_lib.cpp
    ${SOURCE_FILES}
    ../sturm/sturm.cpp
    ../sturm/util.cpp
    ../sturm/solve.h
    ../ConfigFile/ConfigFile.cpp
    ../ConfigFile/ConfigFile.h 
    ../octree/octree.h 
    )    
  
endif()
