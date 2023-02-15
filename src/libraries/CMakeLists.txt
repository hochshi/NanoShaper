################################################################################
# DELPHI_LIB
################################################################################
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
  
  set_target_properties(DelphiSurface
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
  )
endif()


################################################################################
# LIB_BUILD
################################################################################
if(ENABLE_LIB)
  set(CMAKE_CXX_STANDARD 14)
  set(CXX_STANDARD_REQUIRED ON)
  
  add_library(NanoShaper SHARED 
    ${SOURCE_FILES}
    ../sturm/sturm.cpp
    ../sturm/util.cpp
    ../sturm/solve.h
    ../ConfigFile/ConfigFile.cpp
    ../ConfigFile/ConfigFile.h 
    ../octree/octree.h 
    ../ply/ply.cpp
    ../ply/ply.h								
    )    
  
    set_target_properties(NanoShaper
      PROPERTIES
      LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
    )
endif()
