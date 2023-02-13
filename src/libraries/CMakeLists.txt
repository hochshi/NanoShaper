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
  
  if(ENABLE_SPDLOG)
    target_link_libraries(DelphiSurface
      PRIVATE spdlog::spdlog_header_only
    )
  endif()
  if(ENABLE_NLOHMANN_JSON)
    target_link_libraries(DelphiSurface
      PRIVATE nlohmann_json::nlohmann_json
    )
  endif()
endif()

